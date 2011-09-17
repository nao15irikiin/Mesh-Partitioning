#ifndef UNPICKLE_H
#define UNPICKLE_H

#include <assert.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

//datastructure from serialized data (input and expected_output)
struct params {
  //constants
  double dt;
  double nktv2p;
  int ntype;
  double *yeff;
  double *geff;
  double *betaeff;
  double *coeffFrict;

  //node data
  int nnode;
  double *x;
  double *v;
  double *omega;
  double *radius;
  double *mass;
  int    *type;
  double *force;
  double *torque;

  //edge data
  int nedge;
  int *edge;
  double *shear;

  //partition data (OP2 only)
  int npartition;
  std::vector<int> partition_length;

  //expected results
  double *expected_force;
  double *expected_torque;
  double *expected_shear;

  //to pass additinal argument from main:Nao
  int argc;
  char **argv;

};

void print_params(struct params *p);

//unpickle array
template<class T>
inline void fill_array(std::ifstream &file, T *array, int num_elements);

struct params *parse_file(std::string fname);

void parse_partition_file(struct params *input, std::string fname);

void delete_params(struct params *p);

#endif
