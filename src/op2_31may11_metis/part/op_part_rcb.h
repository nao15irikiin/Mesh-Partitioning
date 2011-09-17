#ifndef OP_PART_RCB_H
#define OP_PART_RCB_H

#include <string>
#include "op_part.h"

void part_rcb(int numoflevel, int dims,
              int nnode, int nedge, int nbedge, int ncell,
              point *partnode, point *partedge, point *partbedge, point *partcell,
              int *cell, int *ecell, int *becell, int *edge,
              float *x);


void part_rcb_writefile(int nnode, int nedge, int nbedge, int ncell,
                        string fname_node,  point *partnode,
                        string fname_edge,  point *partedge,
                        string fname_bedge, point *partbedge,
                        string fname_cell,  point *partcell);


#endif

