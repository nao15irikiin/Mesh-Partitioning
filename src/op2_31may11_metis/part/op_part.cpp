#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <vector>
using namespace std;

#include "op_part.h"

vector<int> *prev_parts = NULL;

int comp2(const void *a2, const void *b2) {

	point *a = (point *)a2;
	point *b = (point *)b2;

	if ( prev_parts->at(a->part) == prev_parts->at(b->part) )
		return 0;
	else
		if ( prev_parts->at(a->part) < prev_parts->at(b->part) )
			return -1;
		else
			return 1;
}


void part_MemAlignment(int npart,
							int nnode, int nedge, int nbedge, int ncell,
							point *partnode, point *partedge, point *partbedge, point *partcell, 
							int *node_part_info, int *edge_part_info, int *bedge_part_info, int *cell_part_info,
							int *cell, int *edge, int *ecell, int *bedge, int *becell,
							int *boun, float *x, float *q, float *qold, float *res, float *adt){

	// 1.1 fixup backtrack info - N/A

	// 1.2 sort partedge (In the order of appearance : Nao)
	//     prev_parts[] stores the appearance order of edge
	//     sorted partedge[] stores partition data in the same chunck of block
	int counter = 0;
	prev_parts = new vector<int>(npart, -1);

	for(int i=0; i<nedge; i++){

		int part = partedge[i].part;

		if(prev_parts->at(part) < 0)
			prev_parts->at(part) = counter++;
	}

	qsort(partedge,nedge,sizeof(point),comp2);
	delete prev_parts;

	// 1.3 fixup edge_part_info
	// edge_part_info[] stores the number of cells in each block
	counter = -1;
	int old_part = -1;
	for(int i=0; i<nedge; i++){

		int part = partedge[i].part;

		if(old_part != part){
			counter++;
		}

		old_part = part;
		edge_part_info[counter]++;
	}

	// 1.4 fixup edge watchers i.e. N/A
	
	// 1.5 sort actual edge pointers
	int* edge2 = new int[2*nedge];
	int* ecell2 = new int[2*nedge];

	for(int i=0; i<nedge; i++){

		int index = partedge[i].index;

		for(int j=0; j<2; j++){
			edge2[i*2+j] = edge[index*2+j];
			ecell2[i*2+j] = ecell[index*2+j];
		}
	}

	memcpy(edge, edge2, sizeof(int)*2*nedge);
	memcpy(ecell, ecell2, sizeof(int)*2*nedge); 
	delete[] edge2;
	delete[] ecell2;

/*
	// 2.1.1 Fixup ecell backtrack info
	for(int i=0; i<nedge; i++)
	{
		int index1 = ecell[2*i];
		int index2 = ecell[2*i+1];
		partcell[index1].watcher[0].push_back(2*i);
		partcell[index2].watcher[0].push_back(2*i+1);
	}

	// 2.1.2 Fixup becell backtrack info
	for(int i=0; i<nbedge; i++)
	{
		int index = becell[i];
		partcell[index].watcher[1].push_back(i);
	}

	// 2.2 Sort cells
	counter = 0;
	prev_parts = new vector<int>(npart, -1);
	for(int i=0; i<ncell; i++)
	{
		int part = partcell[i].part;
		if(prev_parts->at(part) < 0)
				prev_parts->at(part) = counter++;
	}
	qsort(partcell,ncell,sizeof(point),comp2);
	delete prev_parts;

	
	// 2.3 fixup cell_part_info
	counter = -1;
	old_part = -1;
	for(int i=0; i<ncell; i++)
	{
		int part = partcell[i].part;
		if(old_part != part)
		{
			counter++;
		}
		old_part = part;
		cell_part_info[counter]++;
	}

	// 2.4 Fixup cell watchers i.e. ecell[] and becell[]
	for(int i=0; i<ncell; i++)
	{
		for(int j=0; j<partcell[i].watcher[0].size(); j++)
		{
			int index = partcell[i].watcher[0][j];
			ecell[index] = i;
		}

		for(int j=0; j<partcell[i].watcher[1].size(); j++)
		{
			int index = partcell[i].watcher[1][j];
			becell[index] = i;
		}

	}
	
	// 2.5 sort cell data and pointers
	int* cell2 = new int[4*ncell];
	float* q2 = new float[4*ncell];
	float* qold2 = new float[4*ncell];
	float* adt2 = new float[1*ncell];
	float* res2 = new float[4*ncell];
	for(int i=0; i<ncell; i++)
	{
		int index = partcell[i].index;
		for(int j=0; j<4; j++)
		{
			cell2[4*i+j] = cell[4*index+j];
//			   q2[4*i+j] =    q[4*index+j];
//			qold2[4*i+j] = qold[4*index+j];
//			 res2[4*i+j] =  res[4*index+j];
		}
//		adt2[i] = adt[index];
	}
	memcpy(cell, cell2, sizeof(int)*4*ncell);
	//memcpy(q, q2, sizeof(float)*4*ncell);
	//memcpy(qold, qold2, sizeof(float)*4*ncell);
	//memcpy(adt, adt2, sizeof(float)*ncell);
	//memcpy(res, res2, sizeof(float)*4*ncell);
	delete[] cell2;
	delete[] q2;
	delete[] qold2;
	delete[] adt2;
	delete[] res2;
*/

/*
	// 3.1.1 Fixup edge backtrack info
	for(int i=0; i<nedge; i++)
	{
		int index1 = edge[2*i];
		int index2 = edge[2*i+1];
		partnode[index1].watcher[0].push_back(2*i);
		partnode[index2].watcher[0].push_back(2*i+1);
	}
	// 3.1.2 Fixup cell backtrack info
	for(int i=0; i<ncell; i++)
	{
		int index1 = cell[4*i];
		int index2 = cell[4*i+1];
		int index3 = cell[4*i+2];
		int index4 = cell[4*i+3];
		partnode[index1].watcher[1].push_back(4*i);
		partnode[index2].watcher[1].push_back(4*i+1);
		partnode[index3].watcher[1].push_back(4*i+2);
		partnode[index4].watcher[1].push_back(4*i+3);
	}

	// 3.1.3 Fixup bedge backtrack info
	for(int i=0; i<nbedge; i++)
	{
		int index1 = bedge[2*i];
		int index2 = bedge[2*i+1];
		partnode[index1].watcher[2].push_back(2*i);
		partnode[index2].watcher[2].push_back(2*i+1);
	}

	// 3.2 Sort nodes
	counter = 0;
	prev_parts = new vector<int>(npart, -1);
	for(int i=0; i<nnode; i++)
	{
		int part = partnode[i].part;
		if(prev_parts->at(part) < 0)
				prev_parts->at(part) = counter++;
	}
	qsort(partnode,nnode,sizeof(point),comp2);
	delete prev_parts;

	// 3.3 Fixup node_part_info
	counter = -1;
	old_part = -1;
	for(int i=0; i<nnode; i++)
	{
		int part = partnode[i].part;
		if(old_part != part)
		{
			counter++;
		}
		old_part = part;
		node_part_info[counter]++;
	}

	// 3.4 Fixup node watchers i.e. edge[], cell[] and bedge[]
	for(int i=0; i<nnode; i++)
	{
		for(int j=0; j<partnode[i].watcher[0].size(); j++)
		{
			int index = partnode[i].watcher[0][j];
			edge[index] = i;
		}

		for(int j=0; j<partnode[i].watcher[1].size(); j++)
		{
			int index = partnode[i].watcher[1][j];
			cell[index] = i;
		}

		for(int j=0; j<partnode[i].watcher[2].size(); j++)
		{
			int index = partnode[i].watcher[2][j];
			bedge[index] = i;
		}
	}

	// 3.5 sort node data
	float* x2 = new float[2*nnode];
	for(int i=0; i<nnode; i++)
	{
		int index = partnode[i].index;
		for(int j=0; j<2; j++)
		{
			x2[2*i+j] = x[2*index+j];
		}
	}
	memcpy(x, x2, sizeof(float)*2*nnode);
	delete[] x2;

  
	// 4.1 fixup backtrack info - N/A

	// 4.2 sort partedge (In the order of appearance : Nao)
	counter = 0;
	prev_parts = new vector<int>(npart, -1);

	for(int i=0; i<nbedge; i++){

		int part = partbedge[i].part;

		if(prev_parts->at(part) < 0)
			prev_parts->at(part) = counter++;
	}

	qsort(partbedge,nbedge,sizeof(point),comp2);
	delete prev_parts;

	// 4.3 fixup bedge_part_info
	counter = -1;
	old_part = -1;
	for(int i=0; i<nbedge; i++){

		int part = partbedge[i].part;

		if(old_part != part){
			counter++;
		}

		old_part = part;
		bedge_part_info[counter]++;
	}

	// 4.4 fixup bedge watchers i.e. N/A
	
	// 4.5 sort actual bedge pointers
	int* bedge2 = new int[2*nbedge];
	int* becell2 = new int[nbedge];
	int* boun2 = new int[nbedge];

	for(int i=0; i<nbedge; i++){

		int index = partbedge[i].index;

		for(int j=0; j<2; j++)
			bedge2[i*2+j] = bedge[index*2+j];

		becell2[i] = becell[index];
		boun2[i] = boun[index];
	}

	memcpy(bedge, bedge2, sizeof(int)*2*nbedge);
	memcpy(becell, becell2, sizeof(int)*nbedge); 
	memcpy(boun, boun2, sizeof(int)*nbedge); 
	delete[] bedge2;
	delete[] becell2;
	delete[] boun2;
*/
}


int part_CalcEdgecut(int nedge, int nnode, 
                     int *edge, point *partnode){

	int nodePart1, nodePart2;
	int edge_cut=0;
	for(int idx_edge=0; idx_edge<nedge; idx_edge++){
		nodePart1 = partnode[edge[idx_edge*2  ]].part;
		nodePart2 = partnode[edge[idx_edge*2+1]].part;
		if(nodePart1!=nodePart2) edge_cut++;
	} 

	return edge_cut;
}






