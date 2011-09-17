
// standard headers
#include <stdlib.h>
#include <stdio.h>
#include <vector>

#include "input.h"

using namespace std;

InputData::InputData(string fname){

  // read in grid
  printf("reading in grid \n");

  FILE *fp;
  if ( (fp = fopen(fname.c_str(),"r")) == NULL) {
    printf("can't open file new_grid.dat\n"); exit(-1);
  }

  p.dims=2;

  if (fscanf(fp,"%d %d %d %d \n",&(p.nnode), &(p.ncell), &(p.nedge), &(p.nbedge)) != 4) {
    printf("error reading from new_grid.dat\n"); exit(-1);
  }

  printf("number of node	:%d\n", p.nnode);
  printf("number of cell	:%d\n", p.ncell);
  printf("number of edge	:%d\n", p.nedge);
  printf("number of bedge	:%d\n", p.nbedge);

  // cell   - mapping between cell and node (1 cell shares 4 nodes)
  // edge   - mapping between edge and node (1 edge shares 2 nodes)
  // ecell  - mapping between edge and cell (1 edge shares 2 cells)
  // bedge  - mapping between ???? and node
  // becell - mapping between ???? and cell
  // bound  - mapping between ???? and ???
	
  p.cell   = (int *) malloc(4*p.ncell*sizeof(int));
  p.edge   = (int *) malloc(2*p.nedge*sizeof(int));
  p.ecell  = (int *) malloc(2*p.nedge*sizeof(int));
  p.bedge  = (int *) malloc(2*p.nbedge*sizeof(int));
  p.becell = (int *) malloc(  p.nbedge*sizeof(int));

  p.bound  = (int *) malloc(  p.nbedge*sizeof(int));
  p.x      = (float *) malloc(2*p.nnode*sizeof(float));

  // Read x coordinate
  for (int n=0; n<(p.nnode); n++) {
    if (fscanf(fp,"%f %f \n",&(p.x[2*n]), &(p.x[2*n+1])) != 2) {
      printf("error reading x data from new_grid.dat\n"); exit(-1);
    }
  }
  printf("Completed reading x data from new_grid.dat\n");

  // Read cell data (cell-node) 
  for (int n=0; n<(p.ncell); n++) {
    if (fscanf(fp,"%d %d %d %d \n",&(p.cell[4*n  ]), &(p.cell[4*n+1]),
                                   &(p.cell[4*n+2]), &(p.cell[4*n+3])) != 4) {
      printf("error reading cell data from new_grid.dat\n"); exit(-1);
    }
  }
  printf("Completed reading cell data from new_grid.dat\n");

  // Read edge data (edge-node) and ecell data (edge-cell) 
  for (int n=0; n<(p.nedge); n++) {
    if (fscanf(fp,"%d %d %d %d \n",&(p.edge[2*n]), &(p.edge[2*n+1]),
                                   &(p.ecell[2*n]),&(p.ecell[2*n+1])) != 4) {
      printf("error reading edge and ecell data from new_grid.dat\n"); exit(-1);
    }
  }
  printf("Completed reading edge and ecell data from new_grid.dat\n");

  for (int n=0; n<(p.nbedge); n++) {
    if (fscanf(fp,"%d %d %d %d \n",&(p.bedge[2*n]),&(p.bedge[2*n+1]),
                                   &(p.becell[n]), &(p.bound[n])) != 4) {
      printf("error reading bedge, becell, bound data from new_grid.dat\n"); exit(-1);
    }
  }
  printf("Completed reading bedge, becell and bound data from new_grid.dat\n");

  fclose(fp);

}

InputData::~InputData(){
  free(p.cell);
  free(p.edge);
  free(p.ecell);
  free(p.bedge);
  free(p.becell);
  free(p.bound);
  free(p.x);
}


struct params* InputData::getParam(){
  return &p;
} 

