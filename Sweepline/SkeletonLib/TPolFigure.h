#include "math.h"
#include "LinkedList.h"
#include "ContourTracer.h"
#include "StructureTD.h"
#include "TriDel.h"
#include "SpanTree.h"
#include <ctime>
#include <QDebug>
#include <QPainter>
#include <set>

class TNode;
class TBone;

typedef unsigned int unsignedint;

extern int CompCount;

/*Таблица для формирования статистики вычислений*/
struct ResultTable {
	int BMPm, BMPn;				/* размеры образа строли - столбцы */
	int Polygons;				/* количество полигональных контуров */
	int ConnectComp;			/* количество связных копонент */
	int Points;					/* количество вершин в полигонах */
	int Vertex;					/* количество вершин в скелете */
	int Edges;					/* количество ребер в скелете */
	unsigned int TimeTrace;		/* время трассировки */
	unsigned int TimeTree;		/* время построения дерева смежности контуров*/
	unsigned int TimeSkelet;	/* время скелетизации */
	unsigned int TimePrun;		/* время стрижки*/
	unsigned int TimeSpectrum;	/* время вычисления спектра */
	unsigned int Total;			/* время общее*/
	ResultTable() {
		BMPm = BMPn = Polygons = ConnectComp = Points = Vertex = Edges = 0;
		TimeTrace = TimeTree = TimeSkelet = TimePrun = TimeSpectrum = Total = 0;
	}
};


class TConnected : public LinkedListElement < TConnected >
	/*Связная компонента многоугольной фигуры */ {
	friend class TPolFigure;
public:
	TContour* Border;					/* Внешний контур компоненты */

	list<TBone*> Bones;
	list<TNode*> Nodes;
	/* скелет компоненты - только то, что лежит внутри нее */
	
	vector<TContour*> HoleList;			/* Список дыр - внутренних контуров компоненты */
	TConnected() { CompCount++; };
	virtual ~TConnected() { HoleList.clear(); CompCount--; };
};


class TPolFigure : public Domain
	/*Граничное и скелетное представление бинарной области */
	/*Описание скелета области в виде списков узлов и костей*/ {
	typedef Domain inherited;
public:
	LinkedListTail<TConnected>* Components;		/*Список компонент многоуг.фигуры*/
	ResultTable RTab;							/*Статистика вычислений*/
	virtual ~TPolFigure();
	TPolFigure(BitRaster* PM, double Amin);		/*Построение компонент для битмапа с отбрасыванием малых контуров*/
	TPolFigure(QString filepath);				/*Чтение компонент из текстового файла*/
	void CreateContours(ContourTracer* BinIm);	/*Построение контуров для бинарного образа BinIm*/
	void MakeComponents();						/*Формирование компонент из контуров*/
	void Invert();								/*Инверсия фигуры*/
	void RestoreInversion();					/*Восстановление фигуры после инверсии*/
	void ClearAll();							/*Чистка всех списков*/
	TPolFigure() : Components(NULL) {}
};