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

/*������� ��� ������������ ���������� ����������*/
struct ResultTable {
	int BMPm, BMPn;				/* ������� ������ ������ - ������� */
	int Polygons;				/* ���������� ������������� �������� */
	int ConnectComp;			/* ���������� ������� �������� */
	int Points;					/* ���������� ������ � ��������� */
	int Vertex;					/* ���������� ������ � ������� */
	int Edges;					/* ���������� ����� � ������� */
	unsigned int TimeTrace;		/* ����� ����������� */
	unsigned int TimeTree;		/* ����� ���������� ������ ��������� ��������*/
	unsigned int TimeSkelet;	/* ����� ������������ */
	unsigned int TimePrun;		/* ����� �������*/
	unsigned int TimeSpectrum;	/* ����� ���������� ������� */
	unsigned int Total;			/* ����� �����*/
	ResultTable() {
		BMPm = BMPn = Polygons = ConnectComp = Points = Vertex = Edges = 0;
		TimeTrace = TimeTree = TimeSkelet = TimePrun = TimeSpectrum = Total = 0;
	}
};


class TConnected : public LinkedListElement < TConnected >
	/*������� ���������� ������������� ������ */ {
	friend class TPolFigure;
public:
	TContour* Border;					/* ������� ������ ���������� */

	list<TBone*> Bones;
	list<TNode*> Nodes;
	/* ������ ���������� - ������ ��, ��� ����� ������ ��� */
	
	vector<TContour*> HoleList;			/* ������ ��� - ���������� �������� ���������� */
	TConnected() { CompCount++; };
	virtual ~TConnected() { HoleList.clear(); CompCount--; };
};


class TPolFigure : public Domain
	/*��������� � ��������� ������������� �������� ������� */
	/*�������� ������� ������� � ���� ������� ����� � ������*/ {
	typedef Domain inherited;
public:
	LinkedListTail<TConnected>* Components;		/*������ ��������� �������.������*/
	ResultTable RTab;							/*���������� ����������*/
	virtual ~TPolFigure();
	TPolFigure(BitRaster* PM, double Amin);		/*���������� ��������� ��� ������� � ������������� ����� ��������*/
	TPolFigure(QString filepath);				/*������ ��������� �� ���������� �����*/
	void CreateContours(ContourTracer* BinIm);	/*���������� �������� ��� ��������� ������ BinIm*/
	void MakeComponents();						/*������������ ��������� �� ��������*/
	void Invert();								/*�������� ������*/
	void RestoreInversion();					/*�������������� ������ ����� ��������*/
	void ClearAll();							/*������ ���� �������*/
	TPolFigure() : Components(NULL) {}
};