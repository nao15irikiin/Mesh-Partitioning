
#ifndef INPUT_H
#define INPUT_H

// standard headers
#include <stdlib.h>
#include <string>

using namespace std;

struct params {
  // dimension size
  int dims;

  // number of sets
  int nnode;
  int ncell;
  int nedge;
  int nbedge;

  // map data
  int *cell;
  int *edge;
  int *ecell;
  int *bedge;
  int *becell;

  // data
  int *bound;
  float *x;
};

class InputData{
private:
  struct params p;

public:
  InputData(string fname);
  ~InputData();
  struct params* getParam();

};


#endif
