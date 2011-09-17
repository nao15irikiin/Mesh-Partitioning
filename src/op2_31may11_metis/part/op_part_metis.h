#ifndef OP_PART_METIS_H
#define OP_PART_METIS_H

#include "op_part.h"

enum eCellType
{
  CELL_TYPE_TRIANGLE = 1,
  CELL_TYPE_TETRAHEDRAL = 2,
  CELL_TYPE_HEXAHEDRAL = 3,
  CELL_TYPE_QUADRILATERAL = 4
};


void part_metis(int metistype, int numofpart, int etype,
                int nnode, int nedge, int nbedge, int ncell,
                point *partnode, point *partedge, point *partbedge, point *partcell,
                int *cell, int *ecell, int *becell);
  
#endif

