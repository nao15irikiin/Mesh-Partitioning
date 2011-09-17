
// standard headers
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

enum eCellType
{
  CELL_TYPE_TRIANGLE = 1,
  CELL_TYPE_TETRAHEDRAL = 2,
  CELL_TYPE_HEXAHEDRAL = 3,
  CELL_TYPE_QUADRILATERAL = 4
};

void readfile(const char* name, 
              int &nnode,
	      int &ncell,
	      int &nedge,
	      int &nbedge,
	      int *&cell,
              int *&edge,
	      int *&ecell,
	      int *&bedge,
	      int *&becell,
	      int *&bound,
	      float *&x,
	      float *&q,
	      float *&qold,
	      float *&res,
	      float *&adt){

  // read in grid
  printf("reading in grid \n");

  FILE *fp;
  if ( (fp = fopen(name,"r")) == NULL) {
    printf("can't open file new_grid.dat\n"); exit(-1);
  }

  if (fscanf(fp,"%d %d %d %d \n",&nnode, &ncell, &nedge, &nbedge) != 4) {
    printf("error reading from new_grid.dat\n"); exit(-1);
  }

  printf("number of node	:%d\n", nnode);
  printf("number of cell	:%d\n", ncell);
  printf("number of edge	:%d\n", nedge);
  printf("number of bedge	:%d\n", nbedge);

  // cell   - mapping between cell and node (1 cell shares 4 nodes)
  // edge   - mapping between edge and node (1 edge shares 2 nodes)
  // ecell  - mapping between edge and cell (1 edge shares 2 cells)
  // bedge  - mapping between ???? and node
  // becell - mapping between ???? and cell
  // bound  - mapping between ???? and ???
	
  cell   = (int *) malloc(4*ncell*sizeof(int));
  edge   = (int *) malloc(2*nedge*sizeof(int));
  ecell  = (int *) malloc(2*nedge*sizeof(int));
  bedge  = (int *) malloc(2*nbedge*sizeof(int));
  becell = (int *) malloc(  nbedge*sizeof(int));
  bound  = (int *) malloc(  nbedge*sizeof(int));

  x      = (float *) malloc(2*nnode*sizeof(float));
  q      = (float *) malloc(4*ncell*sizeof(float));
  qold   = (float *) malloc(4*ncell*sizeof(float));
  res    = (float *) malloc(4*ncell*sizeof(float));
  adt    = (float *) malloc(  ncell*sizeof(float));

  // Read x coordinate
  for (int n=0; n<nnode; n++) {
    if (fscanf(fp,"%f %f \n",&x[2*n], &x[2*n+1]) != 2) {
      printf("error reading x data from new_grid.dat\n"); exit(-1);
    }
  }
  printf("Completed reading x data from new_grid.dat\n");

  // Read cell data (cell-node) 
  for (int n=0; n<ncell; n++) {
    if (fscanf(fp,"%d %d %d %d \n",&cell[4*n  ], &cell[4*n+1],
                                   &cell[4*n+2], &cell[4*n+3]) != 4) {
      printf("error reading cell data from new_grid.dat\n"); exit(-1);
    }
  }
  printf("Completed reading cell data from new_grid.dat\n");

  // Read edge data (edge-node) and ecell data (edge-cell) 
  for (int n=0; n<nedge; n++) {
    if (fscanf(fp,"%d %d %d %d \n",&edge[2*n], &edge[2*n+1],
                                   &ecell[2*n],&ecell[2*n+1]) != 4) {
      printf("error reading edge and ecell data from new_grid.dat\n"); exit(-1);
    }
  }
  printf("Completed reading edge and ecell data from new_grid.dat\n");

  for (int n=0; n<nbedge; n++) {
    if (fscanf(fp,"%d %d %d %d \n",&bedge[2*n],&bedge[2*n+1],
                                   &becell[n], &bound[n]) != 4) {
      printf("error reading bedge, becell, bound data from new_grid.dat\n"); exit(-1);
    }
  }
  printf("Completed reading bedge, becell and bound data from new_grid.dat\n");

  fclose(fp);

}

void output_metis(const char* name, int ncell, int* cell, int type){

  int celldim = 0;

  switch(type){
  case 1:
    celldim = 3;
    break;
  case 2:
  case 4:
    celldim = 4;
    break;
  case 3:
    celldim = 8;
    break;
  default:
    return;
  }

  FILE *fp;
  fp = fopen(name,"w");

  fprintf(fp, "%d %d\n", ncell, type);
  for(int i=0; i<ncell; i++){
    for(int j=0; j<celldim; j++){
      // need +1 since Metis expect node start from 1.
      fprintf(fp, "%d ", cell[i*celldim + j]+1);
    }
    fprintf(fp, "\n");
  }

  fclose(fp);

  /*
  // Export bash script
  int bsize = 64;
  int gridsize = (ncell - 1) / bsize + 1;

  fp = fopen("../dparts/part.sh","w");

  fprintf(fp, "#!/bin/bash\n");
  fprintf(fp, "../../../../../metis-4.0.3/partdmesh ./%s %d", name, gridsize);

  fclose(fp);
  */
}

void output_int(const char* name, 
		int num_element,
		int num_subelement, 
		int* array){

  int tmp;
  int min=99999, max =-1;

  FILE *fp;
  fp = fopen(name,"w");

  fprintf(fp, "%d %d \n", num_element, num_subelement);

  for(int i=0; i<num_element; i++){
    for(int j=0; j<num_subelement; j++){
      tmp =array[i*num_subelement + j];
      if(min>tmp) min=tmp;
      if(max<tmp) max=tmp;
      //fprintf(fp, "%d ", tmp);
    }
    //fprintf(fp, "\n");
  }

  fprintf(fp, "min:%d max:%d \n", min,max);
  fclose(fp);

  printf("Completed writing %s\n", name);

}

void output_float(const char* name, 
		  int num_element,
		  int num_subelement, 
		  float* array){

  float tmp;
  float min=99999.9, max =-1.0;

  FILE *fp;
  fp = fopen(name,"w");

  fprintf(fp, "%d %d \n", num_element, num_subelement);

  for(int i=0; i<num_element; i++){
    for(int j=0; j<num_subelement; j++){
      tmp =array[i*num_subelement + j];
      if(min>tmp) min=tmp;
      if(max<tmp) max=tmp;
      fprintf(fp, "%f ", tmp);
    }
    fprintf(fp, "\n");
  }

  fprintf(fp, "min:%f max:%f \n", min,max);
  fclose(fp);

  printf("Completed writing %s\n", name);

}

void output_vtk(const char* name,
                int nnode, 
                int ncell, 
                float *x, 
                int *cell){

  FILE *fp;
  fp = fopen(name,"w");

  // (1) Header
  fprintf(fp, "# vtk DataFile Version 2.0\n");

  // (2) Title
  fprintf(fp, "OP2 Airfoil Mesh\n");

  // (3) Format
  fprintf(fp, "ASCII\n");

  // (4) Data set
  fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

  // (5.1) Data attributes: Nodes Coordinate
  fprintf(fp, "POINTS %d float\n", nnode);

  for (int i=0; i<nnode; i++){
      fprintf(fp, "%lf %lf 0.0 \n", x[2*i], x[2*i+1]);
  }

  // (5.2) Data attributes: cell-node mapping
  fprintf(fp, "CELLS %d %d\n", ncell, ncell*5);
  for (int i=0; i<ncell; i++){
    fprintf(fp, "4 %d %d %d %d\n", cell[4*i], cell[4*i+1], cell[4*i+2], cell[4*i+3]);
  }

  // (5.3) Data attributes: cell type
  fprintf(fp, "CELL_TYPES %d\n", ncell);
  for (int i=0; i<ncell; i++){
    fprintf(fp, "9\n");
  }

  fclose(fp);

  printf("Completed writing %s\n", name);

}

// main program
int main(int argc, char **argv){

  int    *becell, *ecell,  *bound, *bedge, *edge, *cell;
  float  *x, *q, *qold, *adt, *res;
  int    nnode,ncell,nedge,nbedge,niter;

  readfile("../../../airfoil-input/new_grid.dat",
           nnode, ncell, nedge, nbedge,
	   cell, edge, ecell, bedge, becell, bound,
	   x, q, qold, res, adt);

  //output_metis("../dparts/metis.mesh", ncell, cell, CELL_TYPE_QUADRILATERAL);

  //output_float("grid_x.dat", nnode, 2, x);
  //output_int("grid_cell.dat", ncell, 4, cell);
  //output_int("grid_edge.dat", nedge, 2, edge);
  //output_int("grid_ecell.dat", nedge, 2, ecell);
  //output_int("grid_bedge.dat", nbedge, 2, bedge);
  //output_int("grid_becell.dat", nbedge, 1, becell);
  //output_int("grid_bound.dat", nbedge, 1, bound);

  // This is for paraview file format
  output_vtk("../../../airfoil-input/grid.vtk", nnode, ncell, x, cell);

  free(cell);
  free(edge);
  free(ecell);
  free(bedge);
  free(becell);
  free(bound);

  free(x);
  free(q);
  free(qold);
  free(res);
  free(adt);


}

