#ifndef OP_PART_H
#define OP_PART_H

#include <stdio.h>
#include <vector>


using namespace std;

enum eSMARTPartType
{
  METIS_PART_N_MESH = 1,
  METIS_PART_D_MESH = 2,
  RCB_PART = 10
};


// struct declaration
struct point
{
	int index;
	int part;
	vector<int> watcher[3];

	point(){
		index = 0;
		part = 0;
	}
	
	point(const point& ref){
		index = ref.index;
		part = ref.part;
	}
};

void part_MemAlignment( int npart,
							int nnode, int nedge, int nbedge, int ncell,
							point *partnode, point *partedge, point *partbedge, point *partcell, 
							int *node_part_info, int *edge_part_info, int *bedge_part_info, int *cell_part_info,
							int *cell, int *edge, int *ecell, int *bedge, int *becell,
							int *boun, float *x, float *q, float *qold, float *res, float *adt);


int part_CalcEdgecut(int nedge, int nnode, 
                     int *edge, point *partnode);

#endif

