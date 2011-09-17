#include "user_defined_types.h"
#include "op_datatypes.h"
#include "kernels.h"
__device__
#include <update.h>
__global__

void op_cuda_update(double *arg0,double *arg1,double *arg2,double *arg3,double *arg4,int set_size)
{
  for (int n = threadIdx.x + blockIdx.x * blockDim.x; n < set_size; n += blockDim.x * gridDim.x) {
    update(arg0 + n * 3,arg1 + n * 3,arg2 + n * 3,arg3 + n * 3,arg4 + n * 3);
  }
}

float op_par_loop_update(const char *name,op_set set,struct op_dat<void> *arg0,int idx0,op_ptr *ptr0,enum op_access acc0,struct op_dat<void> *arg1,int idx1,op_ptr *ptr1,enum op_access acc1,struct op_dat<void> *arg2,int idx2,op_ptr *ptr2,enum op_access acc2,struct op_dat<void> *arg3,int idx3,op_ptr *ptr3,enum op_access acc3,struct op_dat<void> *arg4,int idx4,op_ptr *ptr4,enum op_access acc4)
{
  int bsize = BSIZE;
  int gridsize = (set.size - 1) / bsize + 1;
  int reduct_bytes = 0;
  int reduct_size = 0;
  int reduct_shared = reduct_size * (BSIZE / 2);
  int const_bytes = 0;
cudaEvent_t start, stop;
  float elapsed_time_ms = 0.00000F;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start,0);
  op_cuda_update<<<gridsize,bsize,reduct_shared>>>(((double *)arg0->dat_d),((double *)arg1->dat_d),((double *)arg2->dat_d),((double *)arg3->dat_d),((double *)arg4->dat_d),set.size);
  cudaEventRecord(stop,0);
  cudaThreadSynchronize();
  cudaEventElapsedTime(&elapsed_time_ms,start,stop);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);
  return elapsed_time_ms;
}
