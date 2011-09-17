// include files
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "op_part_metis.h"

extern "C"{
#include <metis.h>
}

void part_metisFile_NodeCell(int numofpart, 
			 int nnode, int ncell, 
			 point *partnode, point *partcell){
  
	char number[256];
	sprintf(number, "%d", numofpart);

	printf("metis partition: num of part=%s\n",number);

	char *mesh_data_path;
	mesh_data_path = getenv("OP2_MESH_PATH");

	// Read node partition
	char filepath[128];
	if(mesh_data_path == NULL)
		strcpy(filepath, "../../airfoil-input");
	else{
		strcpy(filepath, mesh_data_path);
	}
	strcat(filepath, "/part/metis");

	char file1[256];
	strcpy(file1, filepath);
	strcat(file1, "/metis.mesh.npart.");
	strcat(file1, number);

	FILE *fp;
	if((fp = fopen(file1,"r")) == NULL){
		printf("can't open file %s\n", file1); exit(-1);
	}

	for(int i=0; i<nnode; i++){
		fscanf(fp, "%d\n", &(partnode[i].part));
		partnode[i].index = i;
	}
	fclose(fp);
	
	// Read cell partition
	char file2[256];
	strcpy(file2, filepath);
	strcat(file2, "/metis.mesh.epart.");
	strcat(file2, number);

	if((fp = fopen(file2,"r")) == NULL){
		printf("can't open file %s\n", file2); exit(-1);
	}
	for(int i=0; i<ncell; i++){
		fscanf(fp, "%d\n", &(partcell[i].part));
		partcell[i].index = i;
	}
	fclose(fp);

}


void part_metisAPI_NodeCell(int metistype,
                            int nparts, int etype, 
                            int nnode, int ncell, 
                            int* cell, 
                            point *partnode, point *partcell){

	int celldim = 0;
	int numflag=0, edgecut;

	switch(etype){
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

	idxtype *elmnts;
	elmnts = (idxtype *)cell;
	//elmnts = idxmalloc(celldim*(ncell), "ReadMesh: elmnts");
	//elmnts = (idxtype*) malloc(celldim*(ncell)*sizeof(int));
	//for (int i=0; i<celldim*ncell; i++) elmnts[i]=cell[i];
 
	idxtype *epart, *npart;
	//epart  = idxmalloc(ncell, "main: epart");
	//npart  = idxmalloc(nnode, "main: npart");
	epart  = (idxtype*) malloc(ncell*sizeof(int));
	npart  = (idxtype*) malloc(nnode*sizeof(int));
	char etypestr[4][5] = {"TRI", "TET", "HEX", "QUAD"};

	printf("**********************************************************************\n");
	printf("%s", METISTITLE);
	printf("Mesh Information ----------------------------------------------------\n");
	printf("  #Elements: %d, #Nodes: %d, Etype: %s\n\n", ncell, nnode, etypestr[etype-1]);

	switch(metistype){
		case METIS_PART_N_MESH:
			printf("Partitioning Nodal Graph... -----------------------------------------\n");
			METIS_PartMeshNodal(&ncell, &nnode, elmnts, &etype, &numflag, &nparts, &edgecut, epart, npart);
			break;

		case METIS_PART_D_MESH:
			printf("Partitioning Dual Graph... ------------------------------------------\n");
			METIS_PartMeshDual(&ncell, &nnode, elmnts, &etype, &numflag, &nparts, &edgecut, epart, npart);
			break;

		default:		
			printf(" \n");
			exit(-1);
	}
	printf("  %d-way Edge-Cut: %7d, Balance: %5.2f\n", nparts, edgecut, ComputeElementBalance(ncell, nparts, epart));
	printf("**********************************************************************\n");

	// Copy partition result data
	for (int i=0; i<ncell; i++){
		partcell[i].part=epart[i];
		partcell[i].index = i;
	}

	for(int i=0; i<nnode; i++){
		partnode[i].part=npart[i];
		partnode[i].index = i;
	}

	// free temporary memory for metis
	//GKfree(&elmnts, &epart, &npart, LTERM);
	//free(elmnts);
	free(epart);
	free(npart);
}

void part_metis_Edge(int nedge, 
                     point *partedge,
                     point *partcell,
                     int *ecell){

	//Generate edge partition information from edge and cell partition
	for(int i=0; i<nedge; i++){
		//int part1 = partcell[ ecell[i*2] ].part;
		int part2 = partcell[ ecell[i*2+1] ].part;
		int target_part = part2;
		
		partedge[i].part = target_part;
		partedge[i].index = i;
	}

}

void part_metis_Bedge(int nbedge, 
                      point *partbedge,
                      point *partcell,
                      int *becell){

	//Generate bedge partition information from edge and cell partition
	int part;
	for(int i=0; i<nbedge; i++){
		part = partcell[ becell[i] ].part;
		
		partbedge[i].part = part;
		partbedge[i].index = i;
	}

}

void part_metis(int metistype,
				int numofpart, int etype, 
				int nnode, int nedge,  int nbedge, int ncell,
				point *partnode, point *partedge, point *partbedge, point *partcell, 
				int *cell, int *ecell, int *becell){

	//part_metisFile_NodeCell(numofpart, nnode, ncell, partnode, partcell);
	part_metisAPI_NodeCell(metistype, numofpart, etype, nnode, ncell, cell, partnode, partcell);

	part_metis_Edge(nedge, partedge, partcell, ecell);

	part_metis_Bedge(nbedge, partbedge, partcell, becell);

/*
	// debug
	for (int i=0; i<nnode; i++){
		if (partnode[i].part!=partnodeAPI[i].part){
			printf("partnode[%d].part=%d partnodeAPI[%d].part=%d\n", i, partnode[i].part, i, partnodeAPI[i].part);
		}
		if (partnode[i].index!=partnodeAPI[i].index){
			printf("partnode[%d].index=%d partnodeAPI[%d].index=%d\n", i, partnode[i].index, i, partnodeAPI[i].index);
		}
	}
	for (int i=0; i<ncell; i++){
		if (partcell[i].part!=partcellAPI[i].part){
			printf("partcell[%d].part=%d partcellAPI[%d].part=%d \n", i, partcell[i].part, i, partcellAPI[i].part);
		}
		if (partcell[i].index!=partcellAPI[i].index){
			printf("partcell[%d].index=%d partcellAPI[%d].index=%d \n", i, partcell[i].index, i, partcellAPI[i].index);
		}
	}

	for (int i=0; i<nedge; i++){
		if (partedge[i].part!=partedgeAPI[i].part){
			printf("partedge[%d].part=%d partedgeAPI[%d].part=%d\n", i, partedge[i].part, i, partedgeAPI[i].part);
		}
		if (partedge[i].index!=partedgeAPI[i].index){
			printf("partedge[%d].index=%d partedgeAPI[%d].index=%d\n", i, partedge[i].index, i, partedgeAPI[i].index);
		}
	}
 */

}


