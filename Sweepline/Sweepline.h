#include <QString>
#include "SkeletonLib/BitRaster.h"
#include "SkeletonLib/TPolFigure.h"

#define CROP_BISECTORS		false

using namespace std;

class TBone;

class TNode {
	friend class TBone;
	TBone* bones[3];
	Element* sites[3];
	Point_3 p;

	bool addBone(TBone* bone) {
		int i = 0;
		for (; i < 3 && bones[i] != NULL; i++);
		if (i == 3) {
			return false;
		}
		else {
			bones[i] = bone;
			return true;
		}
	}

public:
	TNode() {
		for (int i = 0; i < 3; i++) {
			bones[i] = NULL;
			sites[i] = NULL;
		}
		p = Point_3(0, 0, 0);
	};
	TNode(Element* A, Element* B, Element* C) {
		sites[0] = A;
		sites[1] = B;
		sites[2] = C;
		bones[0] = NULL;
		bones[1] = NULL;
		bones[2] = NULL;
		p = Point_3(0, 0, 0);
	};
	void setPoint(Point_3 p_) { p = p_; };
	Point_3 pt() { return p; };
	const TBone* bone(int i) const {
		return bones[i];
	}
};

class TBone {
	Element* sites[2];
	TNode* org_;
	TNode* dest_;
	list<Point_3> path;
public:
	TBone(TNode* org, TNode* dest) : org_(org), dest_(dest) {

		sites[0] = NULL;
		for (int i = 0; i < 3 && org_->sites[i] != NULL; i++) {
			for (int j = 0; j < 3 && dest_->sites[j] != NULL; j++) {
				if (org_->sites[i] == dest_->sites[j]) {
					if (sites[0] == NULL) {
						sites[0] = org_->sites[i];
					}
					else {
						sites[1] = org_->sites[i];
					}
				}
			}
		}
		org->addBone(this);
		dest->addBone(this);

	};
	void setPath(const list<Point_3>& arg) {
		path = arg;
	}
	const list<Point_3>& getPath() const { return path; };
	const TNode* org() const { return org_; };
	const TNode* dest() const { return dest_; };
	const Element* site(int i) const { return sites[i]; };
};

int ComparePEvents(void* Item1, void* Item2);

int CompareBisectors(void* Item1, void* Item2);

enum EventType { POINT, CROSS };

struct Event;

enum HalfPlane { LOWER, UPPER };

struct Bisector {
	list<Point_3> vertices;
	Element* bot_site;
	Element* top_site;
	Event* cross;
	Point_3 A;
	Point_3 B;
	bool lower;

	Bisector(Element* ts, Element* bs, Polygon_2& strel, bool low);
	pair<Gmpq, Vector_3> verticalLevel();
	void cutInHalf(HalfPlane side);
	void cutByStart(Element* prev, Polygon_2& strel);
	void display() {
		for (auto v : vertices) {
			printf("%f %f %f\n", v.x().exact().to_double(), v.y().exact().to_double(), v.z().exact().to_double());
		}
		printf("\n");
	}
};

struct Event {
	Gmpq x;
	Gmpq y;
	Gmpq z;
	Vertex* p;
	Bisector* bot_bis;
	Bisector* top_bis;
	EventType type;
	Event(Gmpq x_, Gmpq y_, Gmpq z_, Vertex* p_) : x(x_), y(y_), z(z_), p(p_), bot_bis(NULL), top_bis(NULL), type(POINT) {};
	Event(Gmpq x_, Gmpq y_, Gmpq z_, Bisector* b1_, Bisector* b2_) : x(x_), y(y_), z(z_), p(NULL), bot_bis(b1_), top_bis(b2_), type(CROSS) {};
};

pair<Gmpq, Gmpq> distanceAndProjection(Point_2 p, Element* site, Polygon_2& strel);
list<Point_3> medialPolyline(Element* elem1, Element* elem2, Polygon_2& strel);
list<Point_3> medialPolyline(Edge* edge1, Edge* edge2, Polygon_2& strel);
list<Point_3> medialPolyline(Edge* edge, Vertex* vertex, Polygon_2& strel);
list<Point_3> medialPolyline(Vertex* vertex1, Vertex* vertex2, Polygon_2& strel);
int compareProjections(Gmpq a, Gmpq b, Polygon_2& strel);
int getDirection(Gmpq dist[3][5], Gmpq proj[3][5], Element sites[3], Polygon_2& strel);

pair<list<Point_3>, list<Segment_3>> polylineIntersection(const list<Point_3>& A, const list<Point_3>& B);

Gmpq paramByCoords(const list<Point_3>& pline, Point_3 a);
list<Point_3> cutByEnds(const list<Point_3>& arg, Point_3 a, Point_3 b);
Point_3 levelPoint(const list<Point_3>& arg, Point_2 level);


class Sweepline
{
	TPolFigure* figure;
	Polygon_2 strel;
	TAVL* events;
	list<Event*> points;
	TAVL* status;
	Event* last_event;

	void readScene(QString filepath);
	void readStrel(QString filepath);
	void normalizeStrel();

	int numberOfRays(Vertex* p);

	void collectPoints();
	void findIntersections(Bisector* bis);
	bool checkStatus(bool display);
	void displayEvents();

	pair<bool, Point_3> bisectorIntersection(Bisector* A, Bisector* B);

	void processBackward(Bisector* bis, Point_3& B);
	void processConcave(Bisector* bis);
	void processCross(Event* e);
	inline void updateStatus(Bisector* bis);
	TKnot* findBisector(Bisector* bis);
	pair<Element*,bool> currentSite(Bisector* bis, Element* test);

public:

	Sweepline();
	list<Bisector*> bisectors;
	void constructSkeleton(QString scene_name, QString strel_name);

};