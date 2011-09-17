
enum eCellType
{
  CELL_TYPE_TRIANGLE = 1,
  CELL_TYPE_TETRAHEDRAL = 2,
  CELL_TYPE_HEXAHEDRAL = 3,
  CELL_TYPE_QUADRILATERAL = 4
};

void readfile(int &nnode, int &ncell, int &nedge, int &nbedge,
              int *&cell, int *&edge, int *&ecell, int *&bedge, int *&becell,
              int *&bound, float *&x, float *&q, float *&qold, float *&res, float *&adt);

void generate_metisdata(int ncell, int nnode, int* cell, int etype, int nparts);

