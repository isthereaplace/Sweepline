#include "Sweepline.h"

int main(int argc, char** argv) {

	QString scene_name(argv[1]);
	QString strel_name(argv[2]);

	Sweepline SL;
	int start = clock();
	SL.constructSkeleton(scene_name, strel_name);
		
	printf("Elapsed time: %d ms\n", clock() - start);

	return 0;
}