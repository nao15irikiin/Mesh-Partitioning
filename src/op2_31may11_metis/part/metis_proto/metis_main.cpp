#include <string.h>
#include <iostream>
#include <time.h>
#include <unistd.h>

#include "input.h"
#include "../op_part_metis.h"
#include "../op_part.h"
#include "../../common/simple_timer.h"

// main program
int main(int argc, char **argv){

    string filename;
    string mode;
    int npart=2048;


    // Read input coordinate data
    //filename = "../../../airfoil-input/new_grid-60-40.dat";
    filename = "../../../airfoil-input/new_grid-600-400.dat";


    // Get arguments
    for (int n=1; n<argc; n++) {
        if (strncmp(argv[n],"FILE=",5)==0) {
            filename = argv[n]+5;
            cout << "Mesh file=" << filename << endl;
        }

        if (strncmp(argv[n],"PART=",5)==0) {
            npart = atoi(argv[n]+5);
            cout << "The number of partition=" << npart <<endl;
        }

        if (strncmp(argv[n],"MODE=",5)==0) {
            mode = argv[n]+5;
            cout << "RCB Mode=" << mode <<endl;
        }
    }


int loop=10;
SimpleTimer part_timer("partition");
for(int i=0; i<loop; i++){

    InputData inputs(filename);

    // get pointer of parameter
    struct params *p = inputs.getParam(); 

    point *partnode  = new point[p->nnode];
    point *partedge  = new point[p->nedge];
    point *partbedge = new point[p->nbedge];
    point *partcell  = new point[p->ncell];

    part_timer.start();
    part_metis(METIS_PART_D_MESH, npart, CELL_TYPE_QUADRILATERAL,
               p->nnode, p->nedge, p->nbedge, p->ncell,
               partnode, partedge, partbedge, partcell,
               p->cell, p->ecell, p->becell);
    part_timer.stop_and_add_to_total();
}

    printf("part time=%lf\n", part_timer.total_time()/loop);

	// Write profile data
	char proffile[256];
	strcpy(proffile, "profile.dat");
	bool isFileExist=true;
	if(access(proffile, R_OK)==-1) isFileExist=false;

	FILE *fp;
	if ( (fp = fopen(proffile,"a")) == NULL) {
		printf("can't generate/find file %s\n", proffile); exit(-1);
	}

	// Write header
	if(isFileExist==false) {
		fprintf(fp, "<filename>, <mode>, <numofpart>, <RCB exec time>, <timestamp>\n");
	}

	time_t timer;
	time(&timer);
	tm *ct = localtime(&timer);
	fprintf(fp, "%s, %s, %d, %lf, %s",filename.c_str(), mode.c_str(), npart, part_timer.total_time()/loop, asctime(ct));
	fclose(fp);

}
