
//
// header file (includes op_lib_core.h and various system header files)
//

#include "op_lib.h"

//
// routines called by user code and kernels
// these wrappers are used by non-CUDA versions
// op_lib.cu provides wrappers for CUDA version
//


void op_init(int argc, char **argv, int diags){
  op_init_core(argc, argv, diags);
}

op_dat op_decl_dat_char(op_set set, int dim, char const *type,
                        int size, char *data, char const *name){
  return op_decl_dat_core(set, dim, type, size, data, name);
}

void op_fetch_data(op_dat dat) {}
void op_decl_const_char(int, char const*, int, char*, char const*){}

op_plan * op_plan_get(char const *name, op_set set, int part_size,
                      int nargs, op_arg *args, int ninds, int *inds){
  return op_plan_core(name, set, part_size, nargs, args, ninds, inds);
}

void op_exit(){
  for(int ip=0; ip<OP_plan_index; ip++) {
    for (int m=0; m<OP_plans[ip].nargs; m++)
      if (OP_plans[ip].loc_maps[m] != NULL)
        free(OP_plans[ip].loc_maps[m]);
    for (int m=0; m<OP_plans[ip].ninds; m++)
      free(OP_plans[ip].ind_maps[m]);
    free(OP_plans[ip].ind_offs);
    free(OP_plans[ip].ind_sizes);
    free(OP_plans[ip].nthrcol);
    free(OP_plans[ip].thrcol);
    free(OP_plans[ip].offset);
    free(OP_plans[ip].nelems);
    free(OP_plans[ip].blkmap);
  }

  op_exit_core();
}
