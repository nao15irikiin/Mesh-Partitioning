
//
// standard headers
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h> // Added by Nao for C99 boolean


void shuffle(int *array, int count, int dim){

	int *tmp = (int*)malloc(dim*sizeof(int));
	int target;

	for (int i = 0; i<count; i++) {
		target=rand()%count;
//		target=i;
		for (int j=0; j<dim; j++) {
			tmp[j]=array[2*i+j];
			array[2*i+j] = array[2*target+j];
			array[2*target+j] = tmp[j];
		}
	}

	free(tmp);
}

void edgeshuffle(int *array1, int *array2, int count, int dim){

	int *tmp1 = (int*)malloc(dim*sizeof(int));
	int *tmp2 = (int*)malloc(dim*sizeof(int));
	int target;

	for (int i = 0; i<count; i++) {
		target=rand()%count;
		for (int j=0; j<dim; j++) {
			tmp1[j]=array1[dim*i+j];
			array1[dim*i+j] = array1[dim*target+j];
			array1[dim*target+j] = tmp1[j];

			tmp2[j]=array2[dim*i+j];
			array2[dim*i+j] = array2[dim*target+j];
			array2[dim*target+j] = tmp2[j];

		}
	}

	free(tmp1);
	free(tmp2);
}

// main program
int main(int argc, char **argv){

	int    *becell, *ecell,  *bound, *bedge, *edge, *cell;
	float  *x, *q, *qold, *adt, *res;

	int    nnode,ncell,nedge,nbedge;


	printf("reading in grid \n");

	// generate input file name
	char *mesh_data_path;
	mesh_data_path = getenv("OP2_MESH_PATH");

	char infile[256];
	if(mesh_data_path == NULL)
		strcpy(infile, "../../../airfoil-input");
	else
		strcpy(infile, mesh_data_path);

	strcat(infile, "/new_grid-600-400.dat");

	FILE *fp;
	if ( (fp = fopen(infile,"r")) == NULL) {
		printf("can't open file %s\n", infile); exit(-1);
	}

	if (fscanf(fp,"%d %d %d %d \n",&nnode, &ncell, &nedge, &nbedge) != 4) {
		printf("error reading from new_grid.dat\n"); exit(-1);
	}

	// <op_map memory allocation>
	// cell   : mapping between cell and node (1 cell shares 4 nodes)
	// edge   : mapping between edge and node (1 edge shares 2 nodes)
	// ecell  : mapping between edge and cell (1 edge shares 2 cells)
	// bedge  : mapping between boundary edge and node
	// becell : mapping between boundary cell and cell	
	cell   = (int *) malloc(4*ncell*sizeof(int));
	edge   = (int *) malloc(2*nedge*sizeof(int));
	ecell  = (int *) malloc(2*nedge*sizeof(int));
	bedge  = (int *) malloc(2*nbedge*sizeof(int));
	becell = (int *) malloc(  nbedge*sizeof(int));

	// <op_data memory allocation>
	bound  = (int *)   malloc(  nbedge*sizeof(int));
	x      = (float *) malloc(2*nnode*sizeof(float));
	q      = (float *) malloc(4*ncell*sizeof(float));
	qold   = (float *) malloc(4*ncell*sizeof(float));
	res    = (float *) malloc(4*ncell*sizeof(float));
	adt    = (float *) malloc(  ncell*sizeof(float));

	// Read x coordinate
	for (int n=0; n<nnode; n++){
		if (fscanf(fp,"%f %f \n",&x[2*n], &x[2*n+1]) != 2) {
			printf("error reading from new_grid.dat\n"); exit(-1);
   		}
	}

	// Read cell data (cell-node) 
	for (int n=0; n<ncell; n++){
		if (fscanf(fp,"%d %d %d %d \n",&cell[4*n  ], &cell[4*n+1],
											&cell[4*n+2], &cell[4*n+3]) != 4) {
			printf("error reading from new_grid.dat\n"); exit(-1);
		}
	}

	// Read edge data (edge-node) and ecell data (edge-cell) 
	for (int n=0; n<nedge; n++){
		if (fscanf(fp,"%d %d %d %d \n",&edge[2*n], &edge[2*n+1],
											&ecell[2*n],&ecell[2*n+1]) != 4) {
			printf("error reading from new_grid.dat\n"); exit(-1);
		}
	}

	// Read boundary edge data (edge-node), boundary edge cell data (edge-cell) 
	for (int n=0; n<nbedge; n++){
		if (fscanf(fp,"%d %d %d %d \n",&bedge[2*n],&bedge[2*n+1],
											&becell[n], &bound[n]) != 4) {
			printf("error reading from new_grid.dat\n"); exit(-1);
		}
	}

	fclose(fp);


	// Randomization
	printf("randomizing data \n");
	edgeshuffle(edge, ecell, nedge, 2);	
	//shuffle(edge, nedge, 2);
	//shuffle(ecell, nedge, 2);	


	// generate output file name
	printf("writing mesh data \n");

	char outfile[256];
	if(mesh_data_path == NULL)
		strcpy(outfile, "../../../airfoil-input");
	else
		strcpy(outfile, mesh_data_path);

	strcat(outfile, "/new_grid_rand-600-400.dat");

	if ( (fp = fopen(outfile,"w")) == NULL) {
		printf("can't generate file %s\n", outfile); exit(-1);
	}

	fprintf(fp,"%d %d %d %d \n",nnode, ncell, nedge, nbedge);

	// Read x coordinate
	for (int n=0; n<nnode; n++)
		fprintf(fp,"%f %f \n",x[2*n], x[2*n+1]);


	// Read cell data (cell-node) 
	for (int n=0; n<ncell; n++)
		fprintf(fp,"%d %d %d %d \n", cell[4*n  ], cell[4*n+1],
										  cell[4*n+2], cell[4*n+3]);

	// Read edge data (edge-node) and ecell data (edge-cell) 
	for (int n=0; n<nedge; n++)
		fprintf(fp,"%d %d %d %d \n", edge[2*n], edge[2*n+1],
										  ecell[2*n],ecell[2*n+1]);

	// Read boundary edge data (edge-node), boundary edge cell data (edge-cell) 
	for (int n=0; n<nbedge; n++)
		fprintf(fp,"%d %d %d %d \n", bedge[2*n], bedge[2*n+1],
											becell[n], bound[n]);

	fclose(fp);

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

