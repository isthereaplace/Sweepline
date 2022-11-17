#include "TPolFigure.h"
#include <fstream>

int CompCount = 0;


unsignedint TimeInMilSecond()
/*Текущее время в миллисекундах*/
{
	return 1000 * (clock() / float(CLOCKS_PER_SEC));
}


void TPolFigure::ClearAll()
{
	TConnected* Com = NULL; /*уничтожение компонент*/
	while (!Components->isEmpty())
	{
		Com = (TConnected*)Components->first();
		delete Com;
	}
}


TPolFigure::~TPolFigure()
{
	ClearAll();
	delete Components;
	Components = NULL;
	// todo check:  inherited::Destroy;
}


TPolFigure::TPolFigure(BitRaster* bitRaster, double Amin)
/*Построение компонент для матрицы с отбрасыванием малых контуров*/
{
	unsigned int TT, TT1;
	TT1 = TimeInMilSecond();
	Components = new LinkedListTail < TConnected >;		/*Порождение списка компонент*/

	ContourTracer* BinIm = new ContourTracer(bitRaster, Amin);
	BinIm->traceContours();

	CreateContours(BinIm);

	/*Упорядочение контуров АВЛ деревом*/
	if (!ElementsExist)
	{
		ProduceElements(this);
		ElementsExist = true;
	}
	TT = TimeInMilSecond();
	SpaningTree(this);
	MakeComponents();

	RTab.TimeTree = TimeInMilSecond() - TT;
	RTab.ConnectComp = CompCount;
	delete BinIm;
	RTab.TimeTrace = TimeInMilSecond() - TT1 - RTab.TimeTree;
	RTab.Total = RTab.TimeTrace + RTab.TimeTree + RTab.TimeSkelet + RTab.TimePrun;
}


TPolFigure::TPolFigure(QString filepath) {
	/*Чтение компонент из текстового файла*/
	unsigned int TT, TT1;
	TT1 = TimeInMilSecond();
	Components = new LinkedListTail < TConnected >;		/*Порождение списка компонент*/

	ifstream in_file(filepath.toStdString());
	double x, y;
	int nComp;
	in_file >> nComp;
	for (int iComp = 0; iComp < nComp; iComp++) {
		int nHoles = 0;
		for (int iCont = 0; iCont <= nHoles; iCont++) {

			TContour* S = AddContour();
			S->Internal = iCont > 0;

			int nPoint;
			in_file >> nPoint;
			for (int iPoint = 0; iPoint < nPoint; iPoint++) {
				in_file >> x;
				in_file >> y;
				S->AddPoint(x, y);
			}

			if (nPoint < 3)
				delete S;
			else
			{
				S->ShiftHead();
				if (S->Internal == S->ConterClockWise())
					S->Invert();
				S->ClosedPath = true;
				RTab.Polygons++;
				RTab.Points = RTab.Points + S->ListPoints->cardinal();
			}
			
			if (iCont == 0) {
				in_file >> nHoles;
			}
		}
	}
	in_file.close();

	/*Упорядочение контуров АВЛ деревом*/
	if (!ElementsExist)
	{
		ProduceElements(this);
		ElementsExist = true;
	}
	TT = TimeInMilSecond();
	SpaningTree(this);
	MakeComponents();

	RTab.TimeTree = TimeInMilSecond() - TT;
	RTab.ConnectComp = CompCount;
	RTab.TimeTrace = TimeInMilSecond() - TT1 - RTab.TimeTree;
	RTab.Total = RTab.TimeTrace + RTab.TimeTree + RTab.TimeSkelet + RTab.TimePrun;
}


void TPolFigure::CreateContours(ContourTracer* BinIm)
/*Построение контуров для бинарного образа BinIm*/
{
	ClosedPath* CP; /* контур из растрового образа (модуль Kontur) */
	RasterPoint* r;  /* точка контура (модуль Kontur) */
	TContour* S;       /* контур для скелетизации (модуль Structure) */
	int i = 0;
	RTab.Points = 0;
	RTab.Polygons = 0;
	CP = BinIm->initialContour();
	while (CP != NULL)
	{
		S = AddContour();
		S->Internal = CP->Internal;
		r = CP->initialPoint();
		i = 0;
		while (r != NULL)
		{
			S->AddPoint(r->x, r->y);
			i++;
			r = r->getNext();
		}
		if (i < 3)
			delete S;
		else
		{
			S->ShiftHead();
			if (S->Internal == S->ConterClockWise())
				S->Invert();
			S->ClosedPath = true;
			RTab.Polygons++;
			RTab.Points = RTab.Points + S->ListPoints->cardinal();
		}
		CP = CP->getNext();
	}
}


void TPolFigure::MakeComponents()
/*Формирование компонент из контуров*/
{
	TContour* S = NULL, * S1 = NULL, * S2 = NULL;
	TConnected* Com = NULL;
	/*Временное размещение внутренних контуров во внешние*/
	S = (TContour*)Boundary->first();
	while (S != NULL)
	{
		S1 = S->getNext();
		if (S->Internal)
		{
			S2 = S->Container;
			if (S2->MySons == NULL)
				S2->MySons = new LinkedListTail < TContour >;
			S->moveIntoTail(S2->MySons);
		}
		S = S1;
	}
	/*Из каждого внешнего создается компонента*/
	S = (TContour*)Boundary->first();
	while (S != NULL)
	{
		if (!S->Internal)
		{
			Com = new TConnected;
			Com->Border = S;
			Com->moveIntoTail(Components);
			while (S->MySons != NULL)
			{
				while (!S->MySons->isEmpty())
				{
					S1 = (TContour*)S->MySons->first();
					Com->HoleList.push_back(S1);
					S1->moveIntoTail(Boundary);
				}
				delete S->MySons;
				S->MySons = NULL;
			}
		}
		S = S->getNext();
	}
}


void TPolFigure::Invert()/*Инверсия фигуры*/
{
	TContour* S = NULL;
	Point* p = NULL;
	double xmin = 0.0, Xmax = 0.0, ymin = 0.0, Ymax = 0.0;
	int W = 0;
	W = 100;
	ClearAll();
	xmin = 10000;
	Xmax = -10000;
	ymin = 10000;
	Ymax = -10000;
	S = (TContour*)Boundary->first();
	while (S != NULL)
	{
		S->Internal = !S->Internal;
		p = (Point*)S->ListPoints->first();
		while (p != NULL)
		{
			if (p->X < xmin)
				xmin = p->X;
			if (p->X > Xmax)
				Xmax = p->X;
			if (p->Y < ymin)
				ymin = p->Y;
			if (p->Y > Ymax)
				Ymax = p->Y;
			p = p->getNext();
		}
		S->Invert();
		S->Container = NULL;
		S->ClosestSite = NULL;
		if (S->MySons != NULL)
			delete S->MySons;
		S->MySons = NULL;
		S = S->getNext();
	}
	S = new TContour;
	S->Internal = false;
	p = new Point(xmin - W, ymin - W);
	p->moveIntoTail(S->ListPoints);
	p = new Point(Xmax + W, ymin - W);
	p->moveIntoTail(S->ListPoints);
	p = new Point(Xmax + W, Ymax + W);
	p->moveIntoTail(S->ListPoints);
	p = new Point(xmin - W, Ymax + W);
	p->moveIntoTail(S->ListPoints);
	S->moveIntoTail(Boundary);
	S->Fiction = true;
	S->ShiftHead();
	if (S->Internal == S->ConterClockWise())
		S->Invert();
	ElementsExist = false;
	ProduceElements(this);
	ElementsExist = true;
}


void TPolFigure::RestoreInversion()
/*Восстановление фигуры после инверсии*/
{
	TContour* S = NULL, * s0 = NULL;
	/*Удаление частей скелета, инцидентных фиктивному контуру*/
	ClearAll();
	//SkelExist = false;
	/*  Com:=Components.first AS TConnected;
	  WHILE Com<>NIL DO
	  BEGIN
	  Bone:=Com.Bones.first AS TBone;
	  WHILE Bone<>NIL DO
	  BEGIN
	  Node:=Bone.Org;
	  Node1:=Bone.Dest;
	  Bone1:=Bone.suc AS TBone;
	  IF Bone.Fiction THEN (*ребро инцидентное фиктивному контуру*)
	  BEGIN
	  Bone.DestroyWithDetach;
	  IF Node.Kind=Isolated THEN Node.Destroy;
	  IF Node1.Kind=Isolated THEN Node1.Destroy;
	  END;
	  Bone:=Bone1;
	  END;
	  Com:=Com.suc AS TConnected;
	  END; */
	  /*Инвертирование направления контуров*/
	S = (TContour*)Boundary->first();
	while (S != NULL)
	{
		S->Internal = !S->Internal;
		if (S->Fiction)
			s0 = S;
		else
		{
			S->Invert();
			S->Container = NULL;
			S->ClosestSite = NULL;
		}
		S = S->getNext();
	}
	delete s0;
}