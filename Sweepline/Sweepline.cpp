#include "Sweepline.h"
#include <QImage>


pair<bool, Point_3> point2point_point2line(const Point_2& p1, const Point_2& q, const Point_2& p2, const Segment_2& e) {
	// Вычислить преобразование, переводящее первую точку в точку, а вторую -- помещающее на линию
	// Ответ в формате [смещение по x, смещение по y, коэффициент растяжения]
	Vector_2 v(-e.target().y() + e.source().y(), e.target().x() - e.source().x());

	CGAL_Matrix A(3, 3);
	A[0][0] = 1;
	A[0][1] = 0;
	A[0][2] = p1.x();
	A[1][0] = 0;
	A[1][1] = 1;
	A[1][2] = p1.y();
	A[2][0] = v.x();
	A[2][1] = v.y();
	A[2][2] = p2.x() * v.x() + p2.y() * v.y();

	if (CGAL::determinant(A[0][0], A[0][1], A[0][2], A[1][0], A[1][1], A[1][2], A[2][0], A[2][1], A[2][2]) == 0) {
		return { false, Point_3() };
	}

	CGAL_Matrix b(3, 1);
	b[0][0] = q.x();
	b[1][0] = q.y();
	b[2][0] = e.source().x() * v.x() + e.source().y() * v.y();

	Gmpq d;
	auto x = (CGAL::Linear_algebraHd<Gmpq>::inverse(A, d) * b).column(0);
	x = 1 / d * x;
	return { true, Point_3(x[0], x[1], x[2]) };
}


Point_3 point2point_point2point(const Point_2& p1, const Point_2& q1, const Point_2& p2, const Point_2& q2) {
	Gmpq z = (q2.x() != q1.x()) ? (q2.x() - q1.x()) / (p2.x() - p1.x()) : (q2.y() - q1.y()) / (p2.y() - p1.y());
	Gmpq x = q1.x() - z * p1.x();
	Gmpq y = q1.y() - z * p1.y();
	return Point_3(x, y, z);
}


Point_2 point2line_point2line(const Point_2& p1, const Segment_2& e1, const Point_2& p2, const Segment_2& e2) {
	Gmpq a11 = e1.target().y() - e1.source().y();
	Gmpq a12 = e1.source().x() - e1.target().x();
	Gmpq a21 = e2.target().y() - e2.source().y();
	Gmpq a22 = e2.source().x() - e2.target().x();
	Gmpq b1 = (e1.source().x() - p1.x()) * a11 + (e1.source().y() - p1.y()) * a12;
	Gmpq b2 = (e2.source().x() - p2.x()) * a21 + (e2.source().y() - p2.y()) * a22;
	Gmpq d = a11 * a22 - a12 * a21;
	
	Gmpq x = a22 * b1 - a12 * b2;
	Gmpq y = a11 * b2 - a21 * b1;
	return Point_2(x / d, y / d);
}


Sweepline::Sweepline() : figure(NULL), events(new TAVL(ComparePEvents)), status(new TAVL(CompareBisectors)), last_event(NULL) {}


void Sweepline::readScene(QString filepath) {
	if (filepath.endsWith(".txt")) {
		figure = new TPolFigure(filepath);
	}
	else {
		QImage image(filepath);
		BitRaster* srcimg = new BitRaster(image.width(), image.height());
		for (int i = 0; i < image.height(); i++) {
			for (int j = 0; j < image.width(); j++) {
				double alpha = qAlpha(image.pixel(j, i)) / 255.0;
				int gray = qGray(image.pixel(j, i)) * alpha + 255 * (1 - alpha);
				bool isBlack = gray < 128;
				srcimg->setBit(j, i, isBlack);
			}
		}
		int AreaIgnore = 0;
		figure = new TPolFigure(srcimg, AreaIgnore);
	}
}


void Sweepline::readStrel(QString filepath) {
	strel.clear();
	ifstream in_file(filepath.toStdString());
	in_file >> strel;
	in_file.close();
}


void Sweepline::normalizeStrel() {
	Point_2 p = *strel.right_vertex();
	for (auto v = strel.vertices_begin(); v != strel.vertices_end(); v++) {
		*v = Point_2(v->x() / p.x(), (v->y() - p.y()) / p.x());
	}
	auto right = strel.vertices_circulator();
	while (*right != *strel.right_vertex()) {
		right++;
	}
}


int ComparePEvents(void* Item1, void* Item2) {
	Event* e1 = (Event*)Item1;
	Event* e2 = (Event*)Item2;
	
	if (e1->x < e2->x)
		return -1;
	if (e2->x < e1->x)
		return 1;

	if (e1->y < e2->y)
		return -1;
	if (e2->y < e1->y)
		return 1;

	if (e1->x - e1->z < e2->x - e2->z)
		return -1;
	if (e1->x - e1->z > e2->x - e2->z)
		return 1;

	if (e1->type == CROSS && e2->type == POINT)
		return -1;
	if (e1->type == POINT && e2->type == CROSS)
		return 1;

	if (e1->type == POINT) {
		if (e1->p->NEl < e2->p->NEl) {
			return -1;
		}
		if (e1->p->NEl > e2->p->NEl) {
			return 1;
		}
		return 0;
	}
	/*
	int bot_comp = CompareBisectors(e1->bot_bis, e2->bot_bis);
	if (bot_comp != 0) {
		return bot_comp;
	}
	int top_comp = CompareBisectors(e1->top_bis, e2->top_bis);
	if (top_comp != 0) {
		return top_comp;
	}
	*/
	if (e1->bot_bis->bot_site->NEl < e2->bot_bis->bot_site->NEl) {
		return -1;
	}
	if (e1->bot_bis->bot_site->NEl > e2->bot_bis->bot_site->NEl) {
		return 1;
	}

	if (e1->bot_bis->top_site->NEl < e2->bot_bis->top_site->NEl) {
		return -1;
	}
	if (e1->bot_bis->top_site->NEl > e2->bot_bis->top_site->NEl) {
		return 1;
	}

	if (e1->top_bis->top_site->NEl < e2->top_bis->top_site->NEl) {
		return -1;
	}
	if (e1->top_bis->top_site->NEl > e2->top_bis->top_site->NEl) {
		return 1;
	}
	return 0;
}


Point_2 SWEEPLINE_LEVEL;
Bisector* BELOW;
Bisector* ABOVE;


int CompareBisectors(void* Item1, void* Item2) {
	Bisector* b1 = (Bisector*)Item1;
	Bisector* b2 = (Bisector*)Item2;

	// Один и тот же бисектор
	if (b1->bot_site == b2->bot_site && b1->top_site == b2->top_site && b1->lower == b2->lower) {
		return 0;
	}
	if (b1->bot_site == b2->top_site && b1->top_site == b2->bot_site && b1->lower == b2->lower) {
		return 0;
	}

	pair<Gmpq, Vector_3> l1 = b1->verticalLevel();
	pair<Gmpq, Vector_3> l2 = b2->verticalLevel();
	if (l1.first < l2.first) {
		return -1;
	}
	if (l1.first > l2.first) {
		return 1;
	}

	if (b1 == BELOW) {
		return -1;
	}
	if (b1 == ABOVE) {
		return 1;
	}
	if (b2 == BELOW) {
		return 1;
	}
	if (b2 == ABOVE) {
		return -1;
	}

	Vector_3 v1 = l1.second;
	Vector_3 v2 = l2.second;

	bool on_the_right = Point_2(b1->A.x(), b1->A.y()) == SWEEPLINE_LEVEL || Point_2(b2->A.x(), b2->A.y()) == SWEEPLINE_LEVEL;

	if (v1.x() + v1.z() == 0 && v1.y() != 0 && v2.x() + v2.z() != 0) {
		return on_the_right ? 1 : -1;
	}
	if (v2.x() + v2.z() == 0 && v2.y() != 0 && v1.x() + v1.z() != 0) {
		return on_the_right ? -1 : 1;
	}
	if (v1.x() + v1.z() != 0 && v2.x() + v2.z() != 0) {
		Gmpq h1 = v1.y() / (v1.x() + v1.z());
		Gmpq h2 = v2.y() / (v2.x() + v2.z());
		if (h1 > h2) {
			return on_the_right ? 1 : -1;
		}
		if (h1 < h2) {
			return on_the_right ? -1 : 1;
		}
		if (b1->top_site->NEl == b2->bot_site->NEl) {
			return -1;
		}
		if (b2->top_site->NEl == b1->bot_site->NEl) {
			return 1;
		}

		if (b1->A.x() + b1->A.z() < b2->A.x() + b2->A.z()) {
			return h1 < 0 ? -1 : 1;
		}
		if (b1->A.x() + b1->A.z() > b2->A.x() + b2->A.z()) {
			return h1 < 0 ? 1 : -1;
		}
		if (b1->A.y() < b2->A.y()) {
			return h1 < 0 ? -1 : 1;
		}
		if (b1->A.y() > b2->A.y()) {
			return h1 < 0 ? 1 : -1;
		}
	}

	if (b1->top_site->NEl == b2->bot_site->NEl) {
		return -1;
	}
	if (b2->top_site->NEl == b1->bot_site->NEl) {
		return 1;
	}

	if (b1->bot_site->NEl < b2->bot_site->NEl) {
		return -1;
	}
	if (b1->bot_site->NEl > b2->bot_site->NEl) {
		return 1;
	}
	if (b1->top_site->NEl < b2->top_site->NEl) {
		return -1;
	}
	if (b1->top_site->NEl > b2->top_site->NEl) {
		return 1;
	}
	return 0;
}


void Sweepline::collectPoints() {
	int num = 0;
	for (auto iComp = figure->Components->first(); iComp != NULL; iComp = iComp->getNext()) {
		for (int i = 0; i < iComp->HoleList.size() + 1; i++) {
			TContour* iCont = (i == 0) ? iComp->Border : iComp->HoleList[i - 1];
			for (auto iElem = iCont->ListElements->first(); iElem != NULL; iElem = iElem->getNext()) {
				iElem->NEl = num++;
				if (iElem->isVertex) {
					Vertex* q = (Vertex*)iElem;
					Event* e = new Event(q->pos.x(), q->pos.y(), 0, q);
					events->Insert(e);
				}
				else {
					((Edge*)iElem)->setBottom(strel);
				}
			}
		}
	}
	points.sort([](Event* a, Event* b) { return a->p->pos < b->p->pos; });
}


Bisector::Bisector(Element* bs, Element* ts, Polygon_2& strel, bool low = false) :
	top_site(ts), bot_site(bs), cross(NULL), lower(low), A(-1,-1,-1), B(-1,-1,-1)
{
	vertices = medialPolyline(bot_site, top_site, strel);

	if (vertices.size() > 1) {

		Point_3 a = *vertices.begin();
		if (a.z() < 0) {
			Point_3 a1(a.x(), a.y(), -a.z());
			Point_3 a2 = *++vertices.begin();
			if (a2.z() < 0) {
				a2 = Point_3(a2.x(), a2.y(), -a2.z());
			}
			Gmpq t = INT_MAX / max(abs(a1.x() - a2.x()), abs(a1.y() - a2.y()));
			a = a2 + t * (a1 - a2);
		}
		
		Point_3 b = *vertices.rbegin();
		if (b.z() < 0) {
			Point_3 b1(b.x(), b.y(), -b.z());
			Point_3 b2 = *++vertices.rbegin();
			if (b2.z() < 0) {
				b2 = Point_3(b2.x(), b2.y(), -b2.z());
			}
			Gmpq t = INT_MAX / max(abs(b1.x() - b2.x()), abs(b1.y() - b2.y()));
			b = b2 + t * (b1 - b2);
		}

		*vertices.begin() = a;
		*vertices.rbegin() = b;
	}

}


void Bisector::cutInHalf(HalfPlane side) {
	if (!bot_site->isVertex && !top_site->isVertex) {
		return;
	}
	if (vertices.size() < 3) {
		return;
	}

	auto left = vertices.begin();
	for (auto v = ++vertices.begin(); v != vertices.end(); v++) {
		if (v->x() + abs(v->z()) < left->x() + abs(left->z()) || (v->x() + abs(v->z()) == left->x() + abs(left->z()) && v->y() < left->y())) {
			left = v;
		}
	}
	auto prev = left;
	if (left != vertices.begin()) {
		prev--;
	}
	else {
		prev++;
	}
	if (prev->y() < left->y() && side == UPPER || prev->y() > left->y() && !(side == UPPER)) {
		// Должны взять часть списка после left
		prev = vertices.begin();
		while (prev != left) {
			prev = vertices.erase(prev);
		}
		return;
	}
	else if (prev->y() > left->y() && side == UPPER || prev->y() < left->y() && !(side == UPPER)){
		// Должны взять часть списка до left
		left++;
		while (left != vertices.end()) {
			left = vertices.erase(left);
		}
		return;
	}

	auto p1 = vertices.begin();
	auto p2 = p1; p2++;
	while (p2 != vertices.end()) {
		if ((p1->y() - SWEEPLINE_LEVEL.y()) * (p2->y() - SWEEPLINE_LEVEL.y()) <= 0) {

			if (p1->y() == p2->y()) {
				return;
			}

			Point_3 q = *p1 + (SWEEPLINE_LEVEL.y() - p1->y()) / (p2->y() - p1->y()) * (*p2 - *p1) ;
			if ((p1->y() < p2->y()) == (side == UPPER)) {
				// Убираем начало
				auto it = vertices.begin();
				while (it != p2) {
					it = vertices.erase(it);
				}
				if (q != vertices.front()) {
					vertices.push_front(q);
				}
				return;
			}
			else {
				// Убираем конец
				while (p2 != vertices.end()) {
					p2 = vertices.erase(p2);
				}
				if (q != vertices.back()) {
					vertices.push_back(q);
				}
				return;
			}
		}
		p1++; p2++;
	}

	printf("I cannot cut in half!\n");
	std::exit(EXIT_FAILURE);
}


pair<Gmpq, Gmpq> distanceAndProjection(Point_2 p, Element* site, Polygon_2& strel) {
	if (!site->isVertex) {
		Edge* edge = (Edge*)site;
		Point_3 res = point2point_point2line(Point_2(0, 0), p, edge->bot, edge->pos).second;
		Gmpq t;
		if (edge->pos.source().x() != edge->pos.target().x()) {
			t = (edge->bot.x() * res.z() + res.x() - edge->pos.source().x()) / (edge->pos.target().x() - edge->pos.source().x());
		}
		else {
			t = (edge->bot.y() * res.z() + res.y() - edge->pos.source().y()) / (edge->pos.target().y() - edge->pos.source().y());
		}
		return { res.z(), t * Gmpq(strel.size()) };
	}

	Point_2 q = ((Vertex*)site)->pos;
	if (q == p) {
		return { 0, 0 };
	}
	Segment_2 seg1(Point_2(0, 0), Point_2(0, 0) + INT_MAX / max(abs(q.x() - p.x()), abs(q.y() - p.y())) * (q - p));
	auto a = strel.vertices_circulator();
	int i = 0;
	do {
		auto b = a; b++;
		Segment_2 seg2(*a, *b);
		if (CGAL::cpp11::result_of<Intersect_2(Segment_2, Segment_2)>::type cross = CGAL::intersection(seg1, seg2)) {
			// Нашли нужный отрезок
			Point_3 res = point2point_point2line(p, Point_2(0, 0), q, seg2).second;
			Gmpq t = (a->x() != b->x()) ? (q.x() * res.z() + res.x() - a->x()) / (b->x() - a->x()) :
				(q.y() * res.z() + res.y() - a->y()) / (b->y() - a->y());
			return { 1 / res.z(), i + t };
		}
		a++; b++; i++;
	} while (a != strel.vertices_circulator());
	
	return { INT_MIN, 0 };
}


int getDirection(Gmpq dist[3][5], Gmpq proj[3][5], Element* sites[3], Polygon_2& strel) {

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 5; j++) {
			printf("%.3f (%.3f)%s", dist[i][j].exact().to_double(), proj[i][j].exact().to_double(), j == 4 ? " \n " : " | ");
		}
	}

	int grades[2][5];
	for (int j = 0; j < 5; j++) {
		Gmpq a = (proj[0][j] >= 0 && proj[0][j] <= Gmpq(strel.size())) ? dist[0][j] : LLONG_MAX;
		for (int i = 1; i < 3; i++) {
			Gmpq b = (proj[i][j] >= 0 && proj[i][j] <= Gmpq(strel.size())) ? dist[i][j] : LLONG_MAX;
			if (b < a) {
				grades[i - 1][j] = (a < Gmpq(LLONG_MAX)) ? -1 : -3;
			}
			else if (b > a) {
				grades[i - 1][j] = (b < Gmpq(LLONG_MAX)) ? 1 : 3;
			}
			else if (a < Gmpq(LLONG_MAX) && sites[0]->isVertex && sites[i]->isVertex) {
				grades[i - 1][j] = 2 * compareProjections(proj[i][j], proj[0][j], strel);
			}
			else {
				grades[i - 1][j] = 0;
			}
		}
	}

	for (int type = 1; type <= 3; type++) {
		// Один из смежных ближе без сброса (-1) или со сбросом (-2)
		for (int i : {-1, +1}) {
			if (grades[0][2 + i] <= 0 && grades[1][2 + i] <= 0 && min(grades[0][2 + i], grades[1][2 + i]) == -type) {
				return i;
			}
		}
		// Один из следующих за смежным ближе без сброса (-1) или со сбросом (-2)
		for (int i : {-1, +1}) {
			if (grades[0][2 + i] == 0 && grades[1][2 + i] == 0 &&
				grades[0][2 + 2 * i] <= 0 && grades[1][2 + 2 * i] <= 0 && min(grades[0][2 + 2 * i], grades[1][2 + 2 * i]) == -type)
			{
				return i;
			}
		}
	}

	for (int type = 1; type <= 3; type++) {
		// Один из смежных дальше без сброса (1) или со сбросом (2)
		for (int i : {-1, +1}) {
			if (grades[0][2 + i] >= 0 && grades[1][2 + i] >= 0 && max(grades[0][2 + i], grades[1][2 + i]) == type) {
				return -i;
			}
		}
		// Один из следующих за смежным дальше без сброса (-1) или со сбросом (2)
		for (int i : {-1, +1}) {
			if (grades[0][2 + i] == 0 && grades[1][2 + i] == 0 &&
				grades[0][2 + 2 * i] >= 0 && grades[1][2 + 2 * i] >= 0 && max(grades[0][2 + 2 * i], grades[1][2 + 2 * i]) == type)
			{
				return -i;
			}
		}
	}

	for (int type = 1; type <= 3; type++) {
		// В одну сторону, по крайней мере, можем идти с уверенностью
		for (int i : {-1, +1}) {
			if (grades[0][2 - i] == 0 && grades[1][2 - i] == 0 && 
				(grades[0][2 + i] > 0 || grades[1][2 + i] > 0) && max(grades[0][2 + i], grades[1][2 + i]) == type)
			{
				return -i;
			}
		}
		// В одну сторону, по крайней мере, можем идти с уверенностью
		for (int i : {-1, +1}) {
			if (grades[0][2 - i] == 0 && grades[1][2 - i] == 0 && grades[0][2 + i] == 0 && grades[1][2 + i] == 0 &&
				grades[0][2 - 2 * i] == 0 && grades[1][2 - 2 * i] == 0 &&
				(grades[0][2 + 2 * i] > 0 || grades[1][2 + 2 * i] > 0) && max(grades[0][2 + 2 * i], grades[1][2 + 2 * i]) == type)
			{
				return -i;
			}
		}
	}

	return 0;
}


void Bisector::cutByStart(Element* prev, Polygon_2& strel) {
	if (vertices.size() < 2) {
		return;
	}
	if (vertices.size() == 2) {
		return;
	}

	Point_2 dir1(0, 0, 0);
	Point_2 dir2(0, 0, 0);
	auto p1 = vertices.begin();
	auto p2 = p1; p2++;
	auto r = vertices.end();
	while (p2 != vertices.end()) {
		if (Segment_3(*p1, *p2).has_on(A)) {
			if (r == vertices.end()) {
				r = p1;
			}
			if (A != *p1) {
				dir1 = Point_2(p2->x() + p2->z() - p1->x() - p1->z(), p2->y() - p1->y());
			}
			if (A != *p2) {
				dir2 = Point_2(p2->x() + p2->z() - p1->x() - p1->z(), p2->y() - p1->y());
			}
		}
		p1++; p2++;
	}
	int dir = 0;
	if (dir1 > Point_2(0, 0) && dir2 >= Point_2(0, 0) || dir1 >= Point_2(0, 0) && dir2 > Point_2(0, 0)) {
		dir = 1;
	}
	else if (dir1 < Point_2(0, 0) && dir2 <= Point_2(0, 0) || dir1 <= Point_2(0, 0) && dir2 < Point_2(0, 0)) {
		dir = -1;
	}

	if (dir == 0) {
		list<Point_3> part;
		auto p1 = vertices.begin();
		auto p2 = p1; p2++;
		while (part.empty()) {
			Segment_3 seg(*p1, *p2);
			if (seg.has_on(A)) {
				part.push_back(*p1);
				if (A != seg.source() && A != seg.target()) {
					part.push_back(A);
				}
				part.push_back(*p2);
				if (p1 != vertices.begin()) {
					part.push_front(*--p1);
				}
				if (p1 != vertices.begin() && A != seg.source()) {
					part.push_front(*--p1);
				}
				p2++;
				if (p2 != vertices.end()) {
					part.push_back(*p2++);
				}
				if (p2 != vertices.end() && A != seg.target()) {
					part.push_back(*p2++);
				}
			}
			else {
				p1++; p2++;
			}
		}

		Element* sites[3] = { prev, bot_site, top_site };
		list<Point_3> crosses;
		for (auto site : sites) {
			list<list<Point_3>> lines;
			if (!site->isVertex) {
				lines = { medialPolyline(site, site->getPrevLooped(), strel), medialPolyline(site->getNextLooped(), site, strel) };
				for (auto& l : lines) {
					Gmpq t = INT_MAX / max(abs(l.back().x() - l.front().x()), abs(l.back().y() - l.front().y()));
					l = list<Point_3>({ l.back() + t * (l.front() - l.back()), l.front() + t * (l.back() - l.front()) });	
				}
			}
			else if (site == prev) {
				Point_3 p(((Vertex*)site)->pos.x(), ((Vertex*)site)->pos.y(), 0);
				for (auto q = strel.vertices_begin(); q != strel.vertices_end(); q++) {
					Gmpq t = INT_MAX / max(abs(q->x()), abs(q->y()));
					Point_3 r = p + t * Vector_3(-q->x(), -q->y(), 1);
					lines.push_back(list<Point_3>({ p, r }));
				}
			}

			for (auto& p : part) {
				p = Point_3(p.x(), p.y(), 0);
			}
			for (auto line : lines) {
				for (auto& p : line) {
					p = Point_3(p.x(), p.y(), 0);
				}
				auto res = polylineIntersection(part, line);
				for (auto p : res.first) {
					crosses.push_back(p);
				}
				for (auto seg : res.second) {
					crosses.push_back(seg.source());
					crosses.push_back(seg.target());
				}
			}
		}

		for (auto q : crosses) {
			auto p1 = part.begin();
			auto p2 = p1; p2++;
			while (!Segment_3(*p1, *p2).has_on(q)) {
				p1++; p2++;
			}
			if (q != *p1 && q != *p2) {
				part.insert(p2, q);
			}
		}

		Point_2 tests[5];
		tests[2] = Point_2(A.x(), A.y());
		auto p = part.begin();
		while (*p != Point_3(A.x(), A.y(), 0)) {
			p++;
		}
		auto q = p;
		for (int i = 0; i < 2; i++) {
			if (q != part.begin()) {
				q--;
				tests[1 - i] = Point_2(q->x(), q->y());
			}
			else {
				tests[1 - i] = tests[2 - i];
			}
		}
		q = p; q++;
		for (int i = 0; i < 2; i++) {
			if (q != part.end()) {
				tests[3 + i] = Point_2(q->x(), q->y());
				q++;
			}
			else {
				tests[3 + i] = tests[2 + i];
			}
		}

		Gmpq dist[3][5];
		Gmpq proj[3][5];
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 5; j++) {
				auto res = distanceAndProjection(tests[j], sites[i], strel);
				dist[i][j] = res.first;
				proj[i][j] = res.second;
			}
		}
		dir = getDirection(dist, proj, sites, strel);

		if (dir == 0) {
			printf("I cannot determine direction!\n");
			std::exit(EXIT_FAILURE);
		}
	}

	// ======================================================================
	
	
	if (dir == 1) {
		r++;
		auto p = vertices.begin();
		while (p != r) {
			p = vertices.erase(p);
		}
		if (A != *r) {
			vertices.push_front(A);
		}
		return;
	}
	else if (dir == -1) {
		auto p = r; p++;
		while (p != vertices.end()) {
			p = vertices.erase(p);
		}
		if (A != *r) {
			vertices.push_back(A);
		}
		return;
	}

}


list<Point_3> medialPolyline(Element* elem1, Element* elem2, Polygon_2& strel) {
	list<Point_3> res;
	if (!elem1->isVertex && !elem2->isVertex) {
		res = medialPolyline((Edge*)elem1, (Edge*)elem2, strel);
	}
	else if (!elem1->isVertex && elem2->isVertex) {
		res = medialPolyline((Edge*)elem1, (Vertex*)elem2, strel);
	}
	else if (elem1->isVertex && !elem2->isVertex) {
		res = medialPolyline((Edge*)elem2, (Vertex*)elem1, strel);
	}
	else{
		res = medialPolyline((Vertex*)elem1, (Vertex*)elem2, strel);
	}
	return res;
}


list<Point_3> medialPolyline(Edge* edge1, Edge* edge2, Polygon_2& strel) {

	list<Point_3> res;
	Segment_2 seg1 = edge1->pos;
	Segment_2 seg2 = edge2->pos;

	bool par = CGAL::parallel(seg1, seg2);
	if (par) {
		if (CGAL::collinear(seg1.source(), seg1.target(), seg2.source())) {
			// Вектора лежат на одной прямой
			return res;
		}
		Vector_2 dir = seg2.source() - seg1.source();
		Vector_2 v1(-seg1.target().y() + seg1.source().y(), seg1.target().x() - seg1.source().x());
		Vector_2 v2(-seg2.target().y() + seg2.source().y(), seg2.target().x() - seg2.source().x());
		if (!(v1 * dir > 0 && v2 * dir < 0)) {
			// Пересечение плоскостей -- это не полоса
			return res;
		}
		else if (!CROP_BISECTORS) {
			Point_2 bot1 = edge1->bot;
			Point_2 bot2 = edge2->bot;
			Point_3 a = point2point_point2line(bot1, seg1.source(), bot2, seg2).second;
			res.push_back(Point_3(a.x(), a.y(), -a.z()));
			res.push_back(Point_3(a.x() + seg1.target().x() - seg1.source().x(), a.y() + seg1.target().y() - seg1.source().y(), -a.z()));
			return res;
		}
	}

	Point_2 bot1 = edge1->bot;
	Point_2 bot2 = edge2->bot;

	if (!CROP_BISECTORS) {
		CGAL::cpp11::result_of<Intersect_2(Line_2, Line_2)>::type cross = CGAL::intersection(Line_2(seg1), Line_2(seg2));
		Point_2 p = boost::get<Point_2>(*cross);
		res.push_back(Point_3(p.x(), p.y(), 0));
		if (bot1 == bot2) {
			res.push_back(Point_3(p.x() - bot1.x(), p.y() - bot1.y(), -1));
		}
		else {
			Point_2 p = point2line_point2line(bot1, seg1, bot2, seg2);
			res.push_back(Point_3(p.x(), p.y(), -1));
		}

		return res;
	}

	CGAL::cpp11::result_of<Intersect_2(Segment_2, Segment_2)>::type cross = CGAL::intersection(seg1, seg2);
	if (cross && bot1 == bot2) {
		Point_2 p = boost::get<Point_2>(*cross);
		res.push_back(Point_3(p.x(), p.y(), 0));
		res.push_back(Point_3(p.x() - bot1.x(), p.y() - bot1.y(), -1));
		return res;
	}

	pair<bool, Point_3> ans;
	Point_3 a;
	Point_2 p;
	bool flag;

	ans = point2point_point2line(bot1, seg1.source(), bot2, seg2);
	a = ans.second;
	p = Point_2(bot2.x() * a.z() + a.x(), bot2.y() * a.z() + a.y());
	if (ans.first && a.z() >= 0 && seg2.collinear_has_on(p)) {
		res.push_back(a);
	}

	if (par) {
		Vector_2 u = seg1.to_vector();
		p += u;
		a = Point_3(a.x() + u.x(), a.y() + u.y(), a.z());
	}
	else {
		ans = point2point_point2line(bot1, seg1.target(), bot2, seg2);
		a = ans.second;
		p = Point_2(bot2.x() * a.z() + a.x(), bot2.y() * a.z() + a.y());
	}
	if (ans.first && a.z() >= 0 && seg2.collinear_has_on(p)) {
		res.push_back(a);
		if (res.size() == 2) {
			return res;
		}
	}

	ans = point2point_point2line(bot2, seg2.source(), bot1, seg1);
	a = ans.second;
	p = Point_2(bot1.x() * a.z() + a.x(), bot1.y() * a.z() + a.y());
	if (ans.first && a.z() >= 0 && seg1.collinear_has_on(p)) {
		if (res.empty() || res.front() != a) {
			res.push_back(a);
			if (res.size() == 2) {
				return res;
			}
		}
	}

	if (par) {
		Vector_2 u = seg2.to_vector();
		p += u;
		a = Point_3(a.x() + u.x(), a.y() + u.y(), a.z());
	}
	else {
		ans = point2point_point2line(bot2, seg2.target(), bot1, seg1);
		a = ans.second;
		p = Point_2(bot1.x() * a.z() + a.x(), bot1.y() * a.z() + a.y());
	}
	if (ans.first && a.z() >= 0 && seg1.collinear_has_on(p)) {
		if (res.empty() || res.front() != a) {
			res.push_back(a);
			if (res.size() == 2) {
				return res;
			}
		}
	}

	return res;
}


list<Point_3> medialPolyline(Edge* edge, Vertex* vertex, Polygon_2& strel) {

	list<Point_3> res;
	Segment_2 seg = edge->pos;
	Point_2 pt = vertex->pos;

	Vector_2 v(-seg.target().y() + seg.source().y(), seg.target().x() - seg.source().x());
	Gmpq val = (pt - seg.source()) * v;
	if (val < 0) {
		return res;
	}
	/*
	if (val == 0 && !seg.collinear_has_on(pt)) {
		return res;
	}
	*/

	Gmpq botval = INT_MAX;
	CGAL::Polygon_circulator<vector<Point_2>> botcir, s;
	auto q = strel.vertices_circulator();
	do {
		Gmpq curval = v * Vector_2(q->x(), q->y());
		if (curval < botval) {
			botval = curval;
			botcir = q;
		}
		q++;
	} while (q != strel.vertices_circulator());

	q = botcir; q++;
	s = botcir; s--;
	Point_2 bot;
	if (v * Vector_2(q->x(), q->y()) == botval) {
		bot = CGAL::midpoint(*botcir, *q);
		s = botcir;
		botcir = q;
	}
	else if (v * Vector_2(s->x(), s->y()) == botval) {
		bot = CGAL::midpoint(*botcir, *s);
	}
	else {
		bot = *botcir;
		s = botcir;
	}

	double bot_[2] = { bot.x().exact().to_double(), bot.y().exact().to_double() };

	if (val == 0) {
		res.push_back(Point_3(pt.x(), pt.y(), 0));
		res.push_back(Point_3(pt.x() - bot.x(), pt.y() - bot.y(), -1));
		return res;
	}

	Point_3 a;
	q = botcir; q++;
	Point_2 p = CGAL::midpoint(*botcir, *q);
	a = point2point_point2line(p, pt, *botcir, seg).second;
	res.push_back(a);

	while (q != s) {
		a = point2point_point2line(*q, pt, *botcir, seg).second;
		res.push_back(a);
		q++;
	}
	q--;
	p = CGAL::midpoint(*q, *s);
	a = point2point_point2line(p, pt, *botcir, seg).second;
	res.push_back(a);

	if (!CROP_BISECTORS) {
		res.front() = Point_3(res.front().x(), res.front().y(), -res.front().z());
		res.back() = Point_3(res.back().x(), res.back().y(), -res.back().z());
		return res;
	}

	pair <bool, Point_3> ans;
	if (res.size() > 1) {

		Point_2 z(0, 0);

		Gmpq tmin = INT_MAX;
		Gmpq tmax = INT_MIN;
		for (int i = 0; i < 2; i++) {

			Gmpq t;
			int k = 0;
			auto c = res.begin();
			auto d = c; d++;
			do {
				auto s = Segment_2(Point_2(c->x(), c->y()), Point_2(d->x(), d->y()));
				ans = point2point_point2line(bot, (i == 0) ? seg.source() : seg.target(), z, s);
				if (ans.first == true) {
					a = ans.second;
					// Положение на прямой, если отрезок [c,d] параметризован сегментом [0,1]
					t = (d->x() != c->x()) ? (a[0] - c->x()) / (d->x() - c->x()) : (a[1] - c->y()) / (d->y() - c->y());
					double t_ = t.exact().to_double();
					if (k == 0 && t <= 1 || k == res.size() - 2 && t >= 0 || (t >= 0 && t < 1)) {
						tmin = min(tmin, k + t);
						tmax = max(tmax, k + t);
					}
				}
				k++; c++; d++;
			} while (d != res.end());

		}

		for (int i = 0; i < 2; i++) {

			auto c = i == 0 ? res.begin() : --res.end();
			auto d = c;
			if (i == 0) {
				d++;
			}
			else {
				d--;
			}

			auto ans1 = point2point_point2line(z, Point_2(c->x(), c->y()), bot, seg);
			Point_2 p(bot.x() * ans1.second[2] + ans1.second[0], bot.y() * ans1.second[2] + ans1.second[1]);
			auto ans2 = point2point_point2line(z, Point_2(d->x(), d->y()), bot, seg);
			Point_2 q(bot.x() * ans2.second[2] + ans2.second[0], bot.y() * ans2.second[2] + ans2.second[1]);

			if (ans1.first && ans2.first && p == q && seg.collinear_has_on(p)) {
				if (i == 0) {
					tmin = INT_MIN;
					tmax = max(tmax, Gmpq(1));
				}
				else {
					tmin = min(tmin, Gmpq(res.size() - 2));
					tmax = INT_MAX;
				}
			}
		}

		double t_[2] = { tmin.exact().to_double(), tmax.exact().to_double() };

		if (tmin > tmax) {
			return list<Point_3>();
		}

		list<Point_3> pline;
		int k = 0;
		auto c = res.begin();
		auto d = c; d++;
		while (k + 1 <= tmin && k < res.size() - 2) {
			c++; d++; k++;
		}
		if (tmin <= INT_MIN) {
			pline.push_back(Point_3(c->x(), c->y(), -c->z()));
		}
		else {
			Point_3 v = *c + (tmin - k) * (*d - *c);
			pline.push_back(v);
		}
		if (tmin == tmax) {
			return pline;
		}

		while (k + 1 < tmax && k < res.size() - 2) {
			pline.push_back(*d);
			c++; d++; k++;
		}

		if (tmax >= INT_MAX) {
			pline.push_back(Point_3(d->x(), d->y(), -d->z()));
		}
		else {
			Point_3 v = *c + (tmax - k) * (*d - *c);
			pline.push_back(v);
		}

		return pline;
	}
	else {
		return list<Point_3>();
	}
}


list<Point_3> medialPolyline(Vertex* vertex1, Vertex* vertex2, Polygon_2& strel) {

	Point_2 point1 = vertex1->pos;
	Point_2 point2 = vertex2->pos;
	list<Point_3> res;

	if (vertex1 == vertex2) {
		return res;
	}

	Vector_2 v(-point2.y() + point1.y(), point2.x() - point1.x());
	CGAL::Polygon_circulator<vector<Point_2>> topcir, botcir;
	Gmpq topval = INT_MIN;
	Gmpq botval = INT_MAX;
	auto q = strel.vertices_circulator();
	do {
		Gmpq curval = v * Vector_2(q->x(), q->y());
		if (curval <= botval) {
			botval = curval;
			botcir = q;
		}
		if (curval >= topval) {
			topval = curval;
			topcir = q;
		}
		q++;
	} while (q != strel.vertices_circulator());

	q = botcir; q++;
	if (v * Vector_2(q->x(), q->y()) == botval) {
		botcir++;
	}
	q = topcir; q++;
	if (v * Vector_2(q->x(), q->y()) == topval) {
		topcir++;
	}

	CGAL::Polygon_circulator<vector<Point_2>> left = botcir;
	left--;
	CGAL::Polygon_circulator<vector<Point_2>> right = botcir;
	right++;

	Point_2 p;
	Point_3 a;
	if (v * Vector_2(left->x(), left->y()) == botval) {
		// Нижняя сторона параллельна отрезку, соединяющему сайты
		Point_2 p1((3 * left->x() + botcir->x()) / 4, (3 * left->y() + botcir->y()) / 4);
		Point_2 p2((left->x() + 3 * botcir->x()) / 4, (left->y() + 3 * botcir->y()) / 4);

		Point_3 a = point2point_point2point(p1, point1, p2, point2);
		res.push_back(Point_3(a[0], a[1], -a[2]));
	}
	else {
		if (v * Vector_2(left->x(), left->y()) < v * Vector_2(right->x(), right->y())) {
			// Левая точка ниже
			p = CGAL::midpoint(*botcir, *left);
			a = point2point_point2line(point1, p, point2, Segment_2(*botcir, *right)).second;
		}
		else {
			p = CGAL::midpoint(*botcir, *right);
			a = point2point_point2line(point2, p, point1, Segment_2(*botcir, *left)).second;
		}
		res.push_back(Point_3(-a[0] / a[2], -a[1] / a[2], -1 / a[2]));
	}

	Gmpq curval = INT_MIN;
	while (left != topcir || right != topcir) {
		Gmpq lval = v * Vector_2(left->x(), left->y());
		Gmpq rval = v * Vector_2(right->x(), right->y());
		if (lval < rval) {
			q = right; q--;
			if (lval > curval) {
				a = point2point_point2line(point1, *left, point2, Segment_2(*q, *right)).second;
				res.push_back(Point_3(-a[0] / a[2], -a[1] / a[2], 1 / a[2]));
				curval = lval;
			}
			left--;
		}
		else {
			q = left; q++;
			if (rval > curval) {
				a = point2point_point2line(point2, *right, point1, Segment_2(*q, *left)).second;
				res.push_back(Point_3(-a[0] / a[2], -a[1] / a[2], 1 / a[2]));
				curval = rval;
			}
			right++;
		}
	}

	left++;
	right--;
	if (v * Vector_2(right->x(), right->y()) == topval) {
		// Верхняя сторона параллельна отрезку, соединяющему сайты
		Point_2 p1((3 * topcir->x() + right->x()) / 4, (3 * topcir->y() + right->y()) / 4);
		Point_2 p2((topcir->x() + 3 * right->x()) / 4, (topcir->y() + 3 * right->y()) / 4);
		a = point2point_point2point(p1, point1, p2, point2);
		res.push_back(Point_3(a[0], a[1], -a[2]));
	}
	else {
		if (v * Vector_2(left->x(), left->y()) > v * Vector_2(right->x(), right->y())) {
			// Левая точка выше
			p = CGAL::midpoint(*topcir, *left);
			a = point2point_point2line(point1, p, point2, Segment_2(*topcir, *right)).second;
		}
		else {
			p = CGAL::midpoint(*topcir, *right);
			a = point2point_point2line(point2, p, point1, Segment_2(*topcir, *left)).second;
		}
		res.push_back(Point_3(-a[0] / a[2], -a[1] / a[2], -1 / a[2]));
	}

	return res;
}


pair<Gmpq, Vector_3> Bisector::verticalLevel(){
	if (vertices.size() == 0) {
		return { INT_MIN, Vector_3() };
	}
	if (vertices.size() == 1) {
		if (vertices.front().x() + vertices.front().z() == SWEEPLINE_LEVEL.x()) {
			return { vertices.front().y(), Vector_3() };
		}
	}
	if (vertices.front().x() + vertices.front().z() == SWEEPLINE_LEVEL.x() && vertices.back().x() + vertices.back().z() == SWEEPLINE_LEVEL.x() &&
		(vertices.front().y() - SWEEPLINE_LEVEL.y()) * (vertices.back().y() - SWEEPLINE_LEVEL.y()) <= 0)
	{
		return { SWEEPLINE_LEVEL.y(), vertices.back() - vertices.front() };
	}

	auto p = vertices.begin();
	auto q = p; q++;
	while (q != vertices.end()) {
		if (p->x() + p->z() == SWEEPLINE_LEVEL.x() && q->x() + q->z() == SWEEPLINE_LEVEL.x()) {
			if (p->y() == q->y()) {
				return { p->y(), *q - *p };
			}
			if ((p->y() - SWEEPLINE_LEVEL.y())*(q->y() - SWEEPLINE_LEVEL.y()) <= 0) {
				return { SWEEPLINE_LEVEL.y(), *q - *p };
			}
		}
		if (p->x() + p->z() != q->x() + q->z() && (p->x() + p->z() - SWEEPLINE_LEVEL.x()) * (q->x() + q->z() - SWEEPLINE_LEVEL.x()) <= 0) {
			Gmpq t = (SWEEPLINE_LEVEL.x() - p->x() - p->z()) / (q->x() + q->z() - p->x() - p->z());
			return { p->y() + t * (q->y() - p->y()), *q - *p };
		}
		p++; q++;
	}

	return { INT_MIN, Vector_3() };
}


int compareProjections(Gmpq a, Gmpq b, Polygon_2& strel) {
	if (max(a, b) - min(a, b) > Gmpq(strel.size() - 1)) {
		if (a < b) {
			a += Gmpq(strel.size());
		}
		else {
			b += Gmpq(strel.size());
		}
	}
	if (a == b || abs(b - a) == 1) {
		return 0;
	}
	if (abs(b - a) < 1) {
		Gmpq i = 0;
		while (i <= min(a, b)) {
			i += 1;
		}
		if (a >= i || b >= i) {
			return 0;
		}
		i -= 1;
		a -= i;
		b -= i;
		return abs(a - 0.5) < abs(b - 0.5) ? -1 : 1;
	}
	return 0;
}


pair<bool, Point_3> Sweepline::bisectorIntersection(Bisector* A, Bisector* B) {

	auto res = polylineIntersection(A->vertices, B->vertices);
	Bisector* C = new Bisector(A->bot_site, B->top_site, strel);

	if (res.first.empty() && !res.second.empty()) {
		auto res1 = polylineIntersection(A->vertices, C->vertices);

		if (res.second.front().squared_length() > INT_MAX && !res1.first.empty()) {
			
			for (auto line : { A, B }) {

				Point_3 q1 = line->vertices.front();
				Point_3 q2 = line->vertices.back();
				auto ans1 = distanceAndProjection(Point_2(q1.x(), q1.y()), line == A ? C->top_site : C->bot_site, strel);
				auto ans2 = distanceAndProjection(Point_2(q2.x(), q2.y()), line == A ? C->top_site : C->bot_site, strel);

				int s1 = ans1.first - q1.z() > 0 ? 1 : (ans1.first - q1.z() < 0 ? -1 : 0);
				if (s1 == 0 && C->bot_site->isVertex && C->top_site->isVertex) {
					auto ans0 = distanceAndProjection(Point_2(q1.x(), q1.y()), line == A ? C->bot_site : C->top_site, strel);
					s1 = compareProjections(ans1.second, ans0.second, strel);
				}
				int s2 = ans2.first - q2.z() > 0 ? 1 : (ans2.first - q2.z() < 0 ? -1 : 0);
				if (s2 == 0 && C->bot_site->isVertex && C->top_site->isVertex) {
					auto ans0 = distanceAndProjection(Point_2(q2.x(), q2.y()), line == A ? C->bot_site : C->top_site, strel);
					s2 = compareProjections(ans1.second, ans0.second, strel);
				}
				if (s1 > 0 && s2 > 0) {
					return { false, Point_3() };
				}
			}
		}

		if (res1.first.empty() && !res1.second.empty()) {
			res.second.clear();
			CGAL::cpp11::result_of<Intersect_3(Segment_3, Segment_3)>::type cross = CGAL::intersection(res1.second.front(), res1.second.front());
			if (cross) {
				if (const Segment_3* w = boost::get<Segment_3>(&*cross)) {
					if (Point_2(w->source().x() + w->source().z(), w->source().y()) <
						Point_2(w->target().x() + w->target().z(), w->target().y()))
					{
						res.first.push_back(w->source());
					}
					else {
						res.first.push_back(w->target());
					}
				}
				else if (const Point_3* p = boost::get<Point_3>(&*cross)) {
					res.first.push_back(*p);
				}
			}
		}
		else {
			res = res1;
		}
	}

	auto p = res.first.begin();
	Element* sites[3] = { A->bot_site, A->top_site, B->top_site };
	while (p != res.first.end()) {
		bool valid = true;
		for (int i = 0; valid && i < 3; i++) {
			if (!(sites[i]->isVertex)) {
				Edge* edge = (Edge*)sites[i];
				Point_2 q = edge->bot;
				Point_2 r(p->x() + p->z() * q.x(), p->y() + p->z() * q.y());
				if (!edge->pos.collinear_has_on(r)) {
					valid = false;
				}
			}
		}

		if (valid) {
			valid = false;
			auto q1 = C->vertices.begin();
			auto q2 = q1; q2++;
			while (!valid && q2 != C->vertices.end()) {
				if (Segment_3(*q1, *q2).has_on(*p)) {
					valid = true;
				}
				q1++; q2++;
			}
		}
		if (valid) {
			p++;
		}
		else {
			p = res.first.erase(p);
		}
	}

	delete C;
	if (res.first.empty()) {
		return { false, Point_3() };
	}
	else if (res.first.size() == 1) {
		return { true, res.first.front() };
	}
	else {
		printf("Two candidates for intersection found!");
		std::exit(EXIT_FAILURE);
	}

}


pair<list<Point_3>, list<Segment_3>> polylineIntersection(const list<Point_3>& A, const list<Point_3>& B) {

	list<Point_3> dots;
	list<Segment_3> segs;

	if (A.size() == 0 || B.size() == 0) {
		return { dots, segs };
	}
	if (A.size() == 1 && B.size() == 1) {
		if (A.front() == B.front()) {
			dots.push_back(A.front());
		}
		return { dots, segs };
	}
	if (A.size() == 1 || B.size() == 1) {
		Point_3 p = (A.size() == 1) ? A.front() : B.front();
		const list<Point_3>& L = (A.size() == 1) ? B : A;
		auto p1 = L.begin();
		for (int i = 0; i + 1 < L.size(); i++, p1++) {
			auto p2 = p1; p2++;
			if (Segment_2(Point_2(p1->x(), p1->y()), Point_2(p2->x(), p2->y())).has_on(Point_2(p.x(), p.y()))) {
				dots.push_back(p);
				return { dots, segs };
			}
		}
		return { dots, segs };
	}

	auto p1 = A.begin();
	for (int i = 0; i + 1 < A.size(); i++, p1++)
	{
		auto p2 = p1; p2++;
		Segment_3 a(*p1, *p2);
		auto q1 = B.begin();
		for (int j = 0; j + 1 < B.size(); j++, q1++) {
			auto q2 = q1; q2++;
			Segment_3 b(*q1, *q2);
			CGAL::cpp11::result_of<Intersect_3(Segment_3, Segment_3)>::type cross = CGAL::intersection(a, b);
			if (cross) {
				if (const Segment_3* w = boost::get<Segment_3>(&*cross)) {
					segs.push_back(*w);
				}
				else {
					const Point_3* p = boost::get<Point_3>(&*cross);
					if (dots.empty() || dots.back() != *p) {
						dots.push_back(*p);
					}
				}
			}
		}
	}

	for (auto s : segs) {
		auto p = dots.begin();
		while (p != dots.end()) {
			if (*p == s.source() || *p == s.target()) {
				p = dots.erase(p);
			}
			else {
				p++;
			}
		}
	}

	return { dots, segs };
}


void Sweepline::findIntersections(Bisector* bis) {

	TKnot* knot = findBisector(bis);
	TKnot* prev = knot->getPrev();
	TKnot* next = knot->getNext();

	Event* bot_event = NULL;
	Event* top_event = NULL;

	if (prev != NULL) {
		Bisector* below = (Bisector*)prev->key;
		if (bis->bot_site == below->top_site && bis->top_site != below->bot_site &&
			!(bis->top_site == bis->bot_site->getNextLooped() && below->bot_site == below->top_site->getPrevLooped()) &&
			!(bis->top_site == bis->bot_site->getPrevLooped() && below->bot_site == below->top_site->getNextLooped()))
		{

			if (bis->bot_site->NEl == 562 && bis->top_site->NEl == 241)
				below->display();

			pair<bool, Point_3> ans = bisectorIntersection(below, bis);
			if (ans.first) {
				Point_3 p = ans.second;
				Event* e = new Event(p.x() + p.z(), p.y(), p.z(), below, bis);
				if (!(e->x < last_event->x || e->x == last_event->x && e->y < last_event->y) &&
					(bot_event == NULL || ComparePEvents(e, bot_event) < 0))
				{
					bot_event = e;
				}
				else {
					delete e;
				}
			}
		}
	}

	if (next != NULL) {
		Bisector* above = (Bisector*)next->key;

		if (bis->top_site == above->bot_site && bis->bot_site != above->top_site &&
			!(bis->bot_site == bis->top_site->getNextLooped() && above->top_site == above->bot_site->getPrevLooped()) &&
			!(bis->bot_site == bis->top_site->getPrevLooped() && above->top_site == above->bot_site->getNextLooped()))
		{
			pair<bool, Point_3> ans = bisectorIntersection(bis, above);
			if (ans.first) {
				Point_3 p = ans.second;
				Event* e = new Event(p.x() + p.z(), p.y(), p.z(), bis, above);
				if (!(e->x < last_event->x || e->x == last_event->x && e->y < last_event->y) &&
					 (top_event == NULL || ComparePEvents(e, top_event) < 0))
				{
					top_event = e;
				}
				else {
					delete e;
				}
			}
		}
	}

	Event* new_event = NULL;
	if (bot_event != NULL) {
		if (top_event != NULL) {
			auto p1 = bis->vertices.begin();
			auto p2 = p1; p2++;
			int i = 0;
			Gmpq t1, t2;
			while (p2 != bis->vertices.end()) {
				Segment_3 seg(*p1, *p2);
				if (seg.has_on(Point_3(bot_event->x - bot_event->z, bot_event->y, bot_event->z))) {
					if (p1->x() != p2->x()) {
						t1 = i + (bot_event->x - bot_event->z - p1->x()) / (p2->x() - p1->x());
					}
					else {
						t1 = i + (bot_event->y - p1->y()) / (p2->y() - p1->y());
					}
				}
				if (seg.has_on(Point_3(top_event->x - top_event->z, top_event->y, top_event->z))) {
					if (p1->x() != p2->x()) {
						t2 = i + (top_event->x - top_event->z - p1->x()) / (p2->x() - p1->x());
					}
					else {
						t2 = i + (top_event->y - p1->y()) / (p2->y() - p1->y());
					}
				}
				i++;
				p1++; p2++;
			}

			if ((bis->bot_site->isVertex || bis->top_site->isVertex) && bis->vertices.back() == bis->A) {
				t1 = -t1;
				t2 = -t2;
			}
			if (!(bis->bot_site->isVertex || bis->top_site->isVertex)) {
				Vector_3 v = bis->vertices.back() - bis->vertices.front();
				if (v.x() + v.z() <= 0 || v.x() + v.z() == 0 && v.y() <= 0) {
					t1 = -t1;
					t2 = -t2;
				}
			}

			if (ComparePEvents(bot_event, top_event) < 0) {
			//if (t1 <= t2) {
				new_event = bot_event;
				delete top_event;
			}
			else {
				new_event = top_event;
				delete bot_event;
			}
		}
		else {
			new_event = bot_event;
		}
	}
	else {
		if (top_event != NULL) {
			new_event = top_event;
		}
	}

	if (new_event != NULL) {

		if ((new_event->bot_bis->cross == NULL || ComparePEvents(new_event, new_event->bot_bis->cross) < 0) &&
			(new_event->top_bis->cross == NULL || ComparePEvents(new_event, new_event->top_bis->cross) < 0))
		{
			Event* e = new_event->bot_bis->cross;
			if (e != NULL) {
				e->bot_bis->cross = NULL;
				e->top_bis->cross = NULL;
				events->Remove(e);
			}
			e = new_event->top_bis->cross;
			if (e != NULL) {
				e->bot_bis->cross = NULL;
				e->top_bis->cross = NULL;
				events->Remove(e);
			}
			new_event->bot_bis->cross = new_event;
			new_event->top_bis->cross = new_event;
			events->Insert(new_event);
		}
		else {
			delete new_event;
		}
	}

}


int Sweepline::numberOfRays(Vertex* p) {
	Edge* e1 = (Edge*)p->getPrevLooped();
	Edge* e2 = (Edge*)p->getNextLooped();
	Vector_2 v1 = e1->pos.to_vector();
	Vector_2 v2 = e2->pos.to_vector();
	if (v1.x() * v2.y() - v2.x() * v1.y() < 0) {
		if (e1->bot != e2->bot) {
			return 2;
		}
	}
	return 1;
}


void Sweepline::processBackward(Bisector* bis, Point_3& B) {
	TKnot* knot = findBisector(bis);
	Bisector* twin = (Bisector*)knot->key;
	if (twin->cross != NULL) {
		Event* e = twin->cross;
		events->Remove(e);
		e->bot_bis->cross = NULL;
		e->top_bis->cross = NULL;
	}
	twin->B = B;
	status->Remove(twin);
	delete bis;
}


pair<Element*,bool> Sweepline::currentSite(Bisector* bis, Element* test) {
	
	TKnot* knot = status->FindNearest(bis);
	if (((Bisector*)knot->key)->verticalLevel().first > SWEEPLINE_LEVEL.y()) {
		knot = knot->getPrev();
	}
	Bisector* below, *above;
	if (knot == NULL) {
		below = NULL;
		above = (Bisector*)(status->Ruler->first())->key;
	}
	else {
		while (knot != NULL && ((Bisector*)knot->key)->verticalLevel().first <= SWEEPLINE_LEVEL.y()) {
			knot = knot->getNext();
		}
		knot = (knot == NULL) ? status->Ruler->last() : knot->getPrev();
		if (((Bisector*)knot->key)->verticalLevel().first == SWEEPLINE_LEVEL.y()) {
			
			TKnot* right = knot;
			TKnot* left = knot->getPrev();
			while (left != NULL && ((Bisector*)left->key)->verticalLevel().first == SWEEPLINE_LEVEL.y()) {
				left = left->getPrev();
			}
			left = (left == NULL) ? status->Ruler->first() : left->getNext();

			Bisector* left_bis = (Bisector*)left->key;
			Bisector* right_bis = (Bisector*)right->key;

			int dir = 0;
			Vector_3 v = right_bis->verticalLevel().second;
			if (v.x() + v.z() != 0) {
				Gmpq dy = v.y() / (v.x() + v.z());
				if (dy < 0) {
					dir = 1;
				}
				else if (dy > 0) {
					dir = -1;
				}
			}
			if (dir == 0 && test != NULL) {
				Bisector* test_bis = new Bisector(right_bis->top_site, test, strel);
				test_bis->cutInHalf(LOWER);
				if (!bisectorIntersection(right_bis, test_bis).first) {
					dir = -1;
				}
				else {
					dir = 1;
				}
				delete test_bis;
			}

			if (dir == -1) {
				return { ((Bisector*)left->key)->bot_site, false };
			}
			else {
				return { ((Bisector*)right->key)->top_site, true };
			}
		}
		else {
			below = (Bisector*)knot->key;
			above = knot->getNext() == NULL ? NULL : (Bisector*)knot->getNext()->key;
		}
	}

	bool from_below;
	if (below == NULL) {
		from_below = false;
	}
	else if (above == NULL || below->top_site == above->bot_site) {
		from_below = true;
	}
	else {
		Gmpq pos1 = distanceAndProjection(SWEEPLINE_LEVEL, below->top_site, strel).first;
		Gmpq pos2 = distanceAndProjection(SWEEPLINE_LEVEL, above->bot_site, strel).first;
		if (pos1 < 0) {
			from_below = false;
		}
		else if (pos2 < 0) {
			from_below = true;
		}
		else {
			from_below = pos1 < pos2;
		}
	}
	return { from_below ? below->top_site : above->bot_site, from_below };
}


void Sweepline::processConcave(Bisector* bis) {
	auto res = currentSite(bis, bis->top_site);
	Element* curr = res.first;
	Vertex* p = (Vertex*)bis->bot_site->getPrevLooped();
	Bisector* bot_bis = new Bisector(curr, p->getPrevLooped(), strel);
	Bisector* top_bis = new Bisector(p->getNextLooped(), curr, strel);
	bot_bis->cutInHalf(LOWER);
	top_bis->cutInHalf(UPPER);

	Point_3 cross = bisectorIntersection(bis, bot_bis).second;
	bis->B = cross;
	bot_bis->A = cross;
	top_bis->A = cross;

	if (res.second) {
		updateStatus(ABOVE = bot_bis);
		updateStatus(ABOVE = top_bis);
	}
	else {
		updateStatus(BELOW = top_bis);
		updateStatus(BELOW = bot_bis);
	}
}


void Sweepline::processCross(Event* e) {
	TKnot* bot = findBisector(e->bot_bis);
	if (bot != NULL) {
		TKnot* top = findBisector(e->top_bis);
		
		//if (true) {
		if (top != NULL && top->getPrev() == bot) {

			Bisector* bis = new Bisector(e->bot_bis->bot_site, e->top_bis->top_site, strel);
			Point_3 cut(e->x - e->z, e->y, e->z);
			((Bisector*)bot->key)->B = cut;
			((Bisector*)top->key)->B = cut;
			bis->A = cut;
			
			if (bis->bot_site->NEl == 300 && bis->top_site->NEl == 139)
				int a = 2;

			bis->cutByStart(e->bot_bis->top_site, strel);

			if (status->FindKey(top->key) != NULL) {
				status->Remove(top->key);
			}
			else {
				TKnot* it = top->getNext();
				while (it != NULL && ((Bisector*)it->key)->verticalLevel().first == SWEEPLINE_LEVEL.y()) {
					it = it->getNext();
				}
				it = (it == NULL) ? status->Ruler->last() : it->getPrev();
				Bisector* above = ABOVE; ABOVE = (Bisector*)it->key;

				if (it == top) {
					status->Remove(ABOVE);
				}
				else {
					Bisector* last = (Bisector*)it->key;
					status->Remove(ABOVE);
					top = findBisector(e->top_bis);
					it = top->getNext();
					while (it != NULL && ((Bisector*)it->key)->verticalLevel().first == SWEEPLINE_LEVEL.y()) {
						it = it->getNext();
					}
					do {
						it = (it == NULL) ? status->Ruler->last() : it->getPrev();
						Bisector* temp = (Bisector*)it->key;
						it->key = last;
						last = temp;
					} while (it != top);
				}
				ABOVE = above;
			}

			TKnot* bot = findBisector(e->bot_bis);
			bot->key = bis;
			bisectors.push_back(bis);
			findIntersections(bis);
		}
	}
	else {
		printf("I can't find this event!\n");
		//std::exit(EXIT_FAILURE);
	}
}


void Sweepline::updateStatus(Bisector* bis) {
	bisectors.push_back(bis);
	status->Insert(bis);
	findIntersections(bis);
}


TKnot* Sweepline::findBisector(Bisector* bis) {
	TKnot* knot = status->FindNearest(bis);
	if (knot == NULL) {
		return NULL;
	}
	Bisector* cand = (Bisector*)knot->key;
	if (cand->bot_site == bis->bot_site && cand->top_site == bis->top_site && cand->lower == bis->lower) {
		return knot;
	}
	if (cand->bot_site == bis->top_site && cand->top_site == bis->bot_site && cand->lower == bis->lower) {
		return knot;
	}

	Gmpq level = bis->verticalLevel().first;
	TKnot* pair = knot->getPrev();
	while (pair != NULL && ((Bisector*)pair->key)->verticalLevel().first == level) {
		Bisector* cand = (Bisector*)pair->key;
		if (cand->bot_site == bis->bot_site && cand->top_site == bis->top_site && cand->lower == bis->lower) {
			return pair;
		}
		if (cand->bot_site == bis->top_site && cand->top_site == bis->bot_site && cand->lower == bis->lower) {
			return pair;
		}
		pair = pair->getPrev();
	}
	pair = knot->getNext();
	while (pair != NULL && ((Bisector*)pair->key)->verticalLevel().first == level) {
		Bisector* cand = (Bisector*)pair->key;
		if (cand->bot_site == bis->bot_site && cand->top_site == bis->top_site && cand->lower == bis->lower) {
			return pair;
		}
		if (cand->bot_site == bis->top_site && cand->top_site == bis->bot_site && cand->lower == bis->lower) {
			return pair;
		}
		pair = pair->getNext();
	}

	return NULL;
}


void Sweepline::constructSkeleton(QString scene_name, QString strel_name) {

	readScene(scene_name);
	readStrel(strel_name);
	normalizeStrel();

	if (figure != NULL && !strel.is_empty()) {

		int num_event = 0;
		collectPoints();
		while (!points.empty() || !events->empty()) {

			Event* e;
			if (!points.empty() && (events->empty() || points.front()->p->pos <
				Point_2(((Event*)events->Ruler->first()->key)->x, ((Event*)events->Ruler->first()->key)->y)))
			{
				e = points.front();
				points.pop_front();
			}
			else {
				e = (Event*)events->MinKey();
			}
			last_event = e;
			SWEEPLINE_LEVEL = Point_2(e->x, e->y);
			BELOW = NULL;
			ABOVE = NULL;

			if (e->type == POINT) {
				if (e->x.exact().to_double() == 16 && e->y.exact().to_double() == 16)
					int a = 1;

				printf("(%f, %f, %f, %d)\n", e->x.exact().to_double(), e->y.exact().to_double(), e->z.exact().to_double(), e->p->NEl);
			}
			else {
				printf("(%f, %f, %f, %d, %d, %d)\n", e->x.exact().to_double(), e->y.exact().to_double(),
					e->z.exact().to_double(), e->bot_bis->bot_site->NEl, e->bot_bis->top_site->NEl, e->top_bis->top_site->NEl);
			}

			if (e->type == POINT) {

				int num_rays = numberOfRays(e->p);
				if (num_rays == 1) {

					Bisector* bis = new Bisector(e->p->getNextLooped(), e->p->getPrevLooped(), strel);
					bis->A = Point_3(e->x, e->y, 0);
					Vector_3 v(bis->vertices.back().x() + bis->vertices.back().z() - e->x, bis->vertices.back().y() - e->y, bis->vertices.back().z());

					if (v.x() > 0 || v.x() == 0 && v.y() > 0) {
						updateStatus(bis);
					}
					else {
						Point_2 prev_pt = ((Vertex*)e->p->getPrevLooped()->getPrevLooped())->pos;
						Point_2 next_pt = ((Vertex*)e->p->getNextLooped()->getNextLooped())->pos;
						if (prev_pt < SWEEPLINE_LEVEL) {
							if (next_pt < SWEEPLINE_LEVEL) {
								Point_3 p(e->x, e->y, 0);
								processBackward(bis, p);
							}
							else {
								// Вставляем первым среди равных
								updateStatus(BELOW = bis);
							}
						}
						else if (next_pt < SWEEPLINE_LEVEL) {
							// Вставляем последним среди равных
							updateStatus(ABOVE = bis);
						}
						else{
							processConcave(bis);
						}
					}

				}
				else {
					Bisector* prev_bis = new Bisector(e->p, e->p->getPrevLooped(), strel);
					Bisector* next_bis = new Bisector(e->p->getNextLooped(), e->p, strel);
					prev_bis->A = Point_3(e->x, e->y, 0);
					next_bis->A = Point_3(e->x, e->y, 0);
					Vector_3 v1(prev_bis->vertices.back().x() + prev_bis->vertices.back().z() - e->x, prev_bis->vertices.back().y() - e->y, prev_bis->vertices.back().z());
					Vector_3 v2(next_bis->vertices.back().x() + next_bis->vertices.back().z() - e->x, next_bis->vertices.back().y() - e->y, next_bis->vertices.back().z());

					if (v1.y() > 0 || v2.y() < 0) {
						Point_2 prev_pt = ((Vertex*)e->p->getPrevLooped()->getPrevLooped())->pos;
						if (prev_pt < SWEEPLINE_LEVEL) {
							updateStatus(BELOW = prev_bis);
							updateStatus(BELOW = next_bis);
						}
						else {
							updateStatus(ABOVE = next_bis);
							updateStatus(ABOVE = prev_bis);
						}
					}
					else if (v1.y() < 0 && v2.y() == 0 || v1.y() == 0 && v1.x() > 0 && v2.y() == 0 ) {
						Point_2 next_pt = ((Vertex*)e->p->getNextLooped()->getNextLooped())->pos;
						if (next_pt < SWEEPLINE_LEVEL) {
							updateStatus(ABOVE = next_bis);
							updateStatus(ABOVE = prev_bis);
						}
						else{
							auto res = currentSite(next_bis, e->p);
							Element* curr = res.first;
							Bisector* bot_bis = new Bisector(curr, e->p, strel);
							Bisector* top_bis = new Bisector(e->p->getNextLooped(), curr, strel);
							
							Point_3 cross = bisectorIntersection(next_bis, bot_bis).second;
							next_bis->B = cross;
							bot_bis->cutInHalf(LOWER);
							bot_bis->A = cross;

							if (res.second) {
								updateStatus(ABOVE = bot_bis);
								updateStatus(ABOVE = prev_bis);
								bisectors.push_back(next_bis);
							}
							if (findBisector(top_bis) != NULL) {
								processBackward(top_bis, cross);
							}
							else {
								top_bis->cutInHalf(UPPER);
								top_bis->A = cross;
								if (res.second) {
									updateStatus(ABOVE = top_bis);
								}
								else {
									updateStatus(BELOW = top_bis);
								}
							}
							if (!res.second) {
								updateStatus(BELOW = prev_bis);
								updateStatus(BELOW = bot_bis);
							}
						}

					}
					else if (v1.y() == 0 && v2.y() > 0 || v1.y() == 0 && v2.y() == 0 && v2.x() > 0) {
						Point_2 prev_pt = ((Vertex*)e->p->getPrevLooped()->getPrevLooped())->pos;
						if (prev_pt < SWEEPLINE_LEVEL) {
							updateStatus(BELOW = prev_bis);
							updateStatus(BELOW = next_bis);
						}
						else {
							auto res = currentSite(prev_bis, e->p->getPrevLooped());
							Element* curr = res.first;
							Bisector* bot_bis = new Bisector(curr, e->p->getPrevLooped(), strel);
							Bisector* top_bis = new Bisector(e->p, curr, strel);

							Point_3 cross = bisectorIntersection(prev_bis, bot_bis).second;
							prev_bis->B = cross;
							top_bis->cutInHalf(UPPER);
							top_bis->A = cross;

							if (!res.second) {
								updateStatus(BELOW = top_bis);
								updateStatus(BELOW = next_bis);
								bisectors.push_back(prev_bis);
							}
							if (findBisector(bot_bis) != NULL) {
								processBackward(bot_bis, cross);
							}
							else {
								bot_bis->cutInHalf(LOWER);
								bot_bis->A = cross;
								if (res.second) {
									updateStatus(ABOVE = bot_bis);
								}
								else {
									updateStatus(BELOW = bot_bis);
								}
							}
							if (res.second) {
								updateStatus(ABOVE = next_bis);
								updateStatus(ABOVE = top_bis);
							}
						}
					}
					else{
						auto res = currentSite(prev_bis, e->p);
						Element* curr = res.first;
						Bisector* bot_bis = new Bisector(curr, e->p, strel, true);
						bot_bis->cutInHalf(LOWER);
						Bisector* top_bis = new Bisector(e->p, curr, strel, false);
						top_bis->cutInHalf(UPPER);
						Point_3 cross = top_bis->vertices.front() == bot_bis->vertices.front() ||
							top_bis->vertices.front() == bot_bis->vertices.back() ?
							top_bis->vertices.front() : top_bis->vertices.back();

						bot_bis->A = cross;
						top_bis->A = cross;

						if (res.second) {
							updateStatus(ABOVE = bot_bis);
							updateStatus(ABOVE = prev_bis);
							updateStatus(ABOVE = next_bis);
							updateStatus(ABOVE = top_bis);
						}
						else {
							updateStatus(BELOW = top_bis);
							updateStatus(BELOW = next_bis);
							updateStatus(BELOW = prev_bis);
							updateStatus(BELOW = bot_bis);
						}

					}
				}
			}
			else if (e->type == CROSS) {
				processCross(e);
			}

			/*
			FILE* fid = fopen(("status" + QString::number(num_event++) + ".txt").toStdString().c_str(), "w");
			int ready = 0;
			for (auto bis : bisectors) {
				if (bis->B != Point_3(-1, -1, -1)) {
					ready++;
				}
			}
			fprintf(fid, "%d\n", ready);
			for (auto bis : bisectors) {
				if (bis->B != Point_3(-1, -1, -1)) {
					list<Point_3> pline = cutByEnds(bis->vertices, bis->A, bis->B);
					fprintf(fid, "%d\n", pline.size());
					for (auto p : pline) {
						fprintf(fid, "%f %f %f\n", p.x().exact().to_double(), p.y().exact().to_double(), p.z().exact().to_double());
					}
				}
			}
			fprintf(fid, "%d\n", status->Count);
			for (auto knot = status->Ruler->first(); knot != NULL; knot = knot->getNext()) {
				Bisector* bis = (Bisector*)knot->key;
				list<Point_3> pline = cutByEnds(bis->vertices, bis->A, levelPoint(bis->vertices, SWEEPLINE_LEVEL));
				fprintf(fid, "%d\n", pline.size());
				for (auto p : pline) {
					fprintf(fid, "%f %f %f\n", p.x().exact().to_double(), p.y().exact().to_double(), p.z().exact().to_double());
				}
			}
			fclose(fid);
			
			if (points.empty() && events->empty()) {
				FILE* fid = fopen((scene_name.left(scene_name.lastIndexOf('.')) + ".poly").toStdString().c_str(), "wb");
				int k = bisectors.size();
				fwrite(&k, sizeof(int), 1, fid);
				for (auto bis : bisectors) {
					bis->vertices = cutByEnds(bis->vertices, bis->A, bis->B);
					double pos[8];
					Element* bot = bis->bot_site;
					Element* top = bis->top_site;
					int count = 0;
					if (!bot->isVertex) {
						pos[0] = ((Edge*)bot)->org->X;
						pos[1] = ((Edge*)bot)->org->Y;
						pos[2] = ((Edge*)bot)->dest->X;
						pos[3] = ((Edge*)bot)->dest->Y;
						if (!top->isVertex) {
							pos[4] = ((Edge*)top)->org->X;
							pos[5] = ((Edge*)top)->org->Y;
							pos[6] = ((Edge*)top)->dest->X;
							pos[7] = ((Edge*)top)->dest->Y;
							count = 8;
						}
						else {
							pos[4] = ((Vertex*)top)->p->X;
							pos[5] = ((Vertex*)top)->p->Y;
							count = 6;
						}
					}
					else {
						if (!top->isVertex) {
							pos[0] = ((Edge*)top)->org->X;
							pos[1] = ((Edge*)top)->org->Y;
							pos[2] = ((Edge*)top)->dest->X;
							pos[3] = ((Edge*)top)->dest->Y;
							pos[4] = ((Vertex*)bot)->p->X;
							pos[5] = ((Vertex*)bot)->p->Y;
							count = 6;
						}
						else {
							pos[0] = ((Vertex*)bot)->p->X;
							pos[1] = ((Vertex*)bot)->p->Y;
							pos[2] = ((Vertex*)top)->p->X;
							pos[3] = ((Vertex*)top)->p->Y;
							count = 4;
						}
					}
					fwrite(&count, sizeof(int), 1, fid);
					fwrite(&pos, sizeof(double), count, fid);

					int k = bis->vertices.size();
					fwrite(&k, sizeof(int), 1, fid);
					for (auto v : bis->vertices) {
						double coords[3] = { v.x().exact().to_double(), v.y().exact().to_double(), v.z().exact().to_double() };
						fwrite(&coords, sizeof(double), 3, fid);
					}
				}
				fclose(fid);
			}
			*/


#ifdef DISPLAY
			bool ans = checkStatus(false);
			//displayEvents();
#else
			//bool ans = checkStatus(false);
			bool ans = true;
#endif
			if (!ans) {
#ifndef DISPLAY
				checkStatus(true);
#endif
				printf("Status is corrupted\n");
				return;
			}
			
		}

	}

}


Gmpq paramByCoords(const list<Point_3>& pline, Point_3 a) {
	if (pline.size() == 1) {
		return 0;
	}
	Gmpq res;

	auto p = pline.begin();
	auto q = p; q++;
	int k = 0;

	while (k < pline.size() - 1) {

		if (CGAL::collinear(Point_2(p->x(), p->y()), Point_2(q->x(), q->y()), Point_2(a.x(), a.y()))) {
			if (p->z() < 0) {
				// Бесконечный луч сначала
				Gmpq val = (p->x() != q->x()) ? ((a.x() - q->x()) * (p->x() - q->x())) : ((a.y() - q->y()) * (p->y() - q->y()));
				if (val >= 0) {
					res = (p->x() != q->x()) ? (a.x() - p->x()) / (q->x() - p->x()) : (a.y() - p->y()) / (q->y() - p->y());
					return k + res;
				}
			}
			else if (q->z() < 0) {
				// Бесконечный луч с конца
				Gmpq val = (p->x() != q->x()) ? ((a.x() - p->x()) * (q->x() - p->x())) : ((a.y() - p->y()) * (q->y() - p->y()));
				if (val >= 0) {
					res = (p->x() != q->x()) ? (a.x() - p->x()) / (q->x() - p->x()) : (a.y() - p->y()) / (q->y() - p->y());
					return k + res;
				}
			}
			else {
				if (Segment_2(Point_2(p->x(), p->y()), Point_2(q->x(), q->y())).collinear_has_on(Point_2(a.x(), a.y()))) {
					res = (p->x() != q->x()) ? (a.x() - p->x()) / (q->x() - p->x()) : (a.y() - p->y()) / (q->y() - p->y());
					return k + res;
				}
			}
		}
		p++; q++; k++;
	}
	return res;
}


Point_3 levelPoint(const list<Point_3>& arg, Point_2 level) {

	if (arg.size() == 0) {
		return Point_3();
	}
	if (arg.size() == 1) {
		if (arg.front().x() + arg.front().z() == SWEEPLINE_LEVEL.x()) {
			return arg.front();
		}
	}
	if (arg.front().x() + arg.front().z() == SWEEPLINE_LEVEL.x() && arg.back().x() + arg.back().z() == SWEEPLINE_LEVEL.x() &&
		(arg.front().y() - SWEEPLINE_LEVEL.y()) * (arg.back().y() - SWEEPLINE_LEVEL.y()) <= 0)
	{
		if (arg.front().y() == arg.back().y()) {
			return arg.back();
		}
		else {
			Gmpq t = (SWEEPLINE_LEVEL.y() - arg.front().x()) / (arg.back().x() - arg.front().x());
			return arg.front() + t * (arg.back() - arg.front());
		}
	}

	auto p = arg.begin();
	auto q = p; q++;
	while (q != arg.end()) {
		if (p->x() + p->z() == SWEEPLINE_LEVEL.x() && q->x() + q->z() == SWEEPLINE_LEVEL.x()) {
			if (p->y() == q->y()) {
				return *q;
			}
			if ((p->y() - SWEEPLINE_LEVEL.y()) * (q->y() - SWEEPLINE_LEVEL.y()) <= 0) {
				Gmpq t = (SWEEPLINE_LEVEL.y() - p->y()) / (q->y() - p->y());
				return *p + t * (*q - *p);
			}
		}
		if (p->x() + p->z() != q->x() + q->z() && (p->x() + p->z() - SWEEPLINE_LEVEL.x()) * (q->x() + q->z() - SWEEPLINE_LEVEL.x()) <= 0) {
			Gmpq t = (SWEEPLINE_LEVEL.x() - p->x() - p->z()) / (q->x() + q->z() - p->x() - p->z());
			return *p + t * (*q - *p);
		}
		p++; q++;
	}

	return Point_3();
}


list<Point_3> cutByEnds(const list<Point_3>& arg, Point_3 a, Point_3 b) {

	Gmpq t1 = paramByCoords(arg, a);
	Gmpq t2 = paramByCoords(arg, b);

	list<Point_3> res;
	if (arg.size() == 1) {
		res.push_back(arg.front());
		return res;
	}

	Point_3 q;
	if (t1 <= t2) {
		int k = 0;
		auto p1 = arg.begin();
		while (k + 1 <= t1 && k < arg.size() - 2) {
			p1++;
			k++;
		}
		auto p2 = p1; p2++;
		q = *p1 + (t1 - k) * (*p2 - *p1);
		res.push_back(Point_3(q.x(), q.y(), (1 - t1 + k) * abs(p1->z()) + (t1 - k) * abs(p2->z())));

		while (k + 1 < t2 && k < arg.size() - 2) {
			res.push_back(*p2);
			p1++; p2++; k++;
		}
		q = (*p1 + (t2 - k) * (*p2 - *p1));
		res.push_back(Point_3(q.x(), q.y(), (1 - t2 + k) * abs(p1->z()) + (t2 - k) * abs(p2->z())));
	}
	else {
		int k = arg.size() - 1;
		auto p1 = arg.rbegin();
		while (k - 1 >= t1 && k > 1) {
			p1++;
			k--;
		}
		auto p2 = p1; p2++;
		q = (*p1 + (k - t1) * (*p2 - *p1));
		res.push_back(Point_3(q.x(), q.y(), (1 - k + t1) * abs(p1->z()) + (k - t1) * abs(p2->z())));

		while (k - 1 > t2 && k > 1) {
			res.push_back(*p2);
			p1++; p2++; k--;
		}
		q = (*p1 + (k - t2) * (*p2 - *p1));
		res.push_back(Point_3(q.x(), q.y(), (1 - k + t2) * abs(p1->z()) + (k - t2) * abs(p2->z())));
	}

	return res;
}


bool Sweepline::checkStatus(bool display) {
	bool ans = true;
	Gmpq y_curr = INT_MIN;
	auto r = status->Ruler->first();
	while (r != NULL) {
		Bisector* bis = (Bisector*)r->key;
		Gmpq y_prev = y_curr;
		y_curr = bis->verticalLevel().first;
		if (y_curr < y_prev) {
			ans = false;
		}
		if (display) {
			printf("[%d, %d, %f]%s\n", bis->bot_site->NEl, bis->top_site->NEl, y_curr.exact().to_double(), y_curr < y_prev ? " - DEFECT" : "");
		}
		r = r->getNext();
	}
	if (display) {
		printf("\n");
	}
	return ans;
}


void Sweepline::displayEvents() {
	auto r = events->Ruler->first();
	while (r != NULL) {
		Event* e = (Event*)r->key;
		if (e->type == POINT) {
			//printf("(%f, %f, %f, %d)\n", e->x.exact().to_double(), e->y.exact().to_double(), e->p->NEl);
		}
		else {
			printf("(%f, %f, %f, %d, %d, %d)\n", e->x.exact().to_double(), e->y.exact().to_double(),
				   e->z.exact().to_double(), e->bot_bis->bot_site->NEl, e->bot_bis->top_site->NEl, e->top_bis->top_site->NEl);
		}
		r = r->getNext();
	}
	printf("\n");
}

