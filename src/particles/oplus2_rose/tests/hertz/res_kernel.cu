#include "user_defined_types.h"
#include "op_datatypes.h"
#include "kernels.h"
__device__
#include <res.h>
__global__

void op_cuda_res(double *ind_arg0,int *ind_arg0_ptrs,int *ind_arg0_sizes,int *ind_arg0_offset,double *ind_arg1,int *ind_arg1_ptrs,int *ind_arg1_sizes,int *ind_arg1_offset,double *ind_arg2,int *ind_arg2_ptrs,int *ind_arg2_sizes,int *ind_arg2_offset,struct contact *arg0_d,double *arg1_d,int *arg2_ptrs,int *arg3_ptrs,int *arg4_ptrs,int *arg5_ptrs,int block_offset,int *blkmap,int *offset,int *nelems,int *ncolors,int *colors)
{
  struct contact arg0_l[1];
  double arg1_l[3];
  double arg2_l[3];
  double arg3_l[3];
  double arg4_l[3];
  double arg5_l[3];
  extern __shared__ 
  char shared[];
  __shared__ 
  int *ind_arg0_ptr;
  __shared__ 
  int *ind_arg1_ptr;
  __shared__ 
  int *ind_arg2_ptr;
  __shared__ 
  int ind_arg0_size;
  __shared__ 
  int ind_arg1_size;
  __shared__ 
  int ind_arg2_size;
  __shared__ 
  double *ind_arg0_s;
  __shared__ 
  double *ind_arg1_s;
  __shared__ 
  double *ind_arg2_s;
  __shared__ 
  struct contact *arg0;
  __shared__ 
  double *arg1;
  __shared__ 
  int *arg2_ptr;
  __shared__ 
  int *arg3_ptr;
  __shared__ 
  int *arg4_ptr;
  __shared__ 
  int *arg5_ptr;
  __shared__ 
  int nelem2;
  __shared__ 
  int ncolor;
  __shared__ 
  int *color;
  __shared__ 
  int blockId;
  __shared__ 
  int nelem;
  if (threadIdx.x == 0) {
    blockId = blkmap[blockIdx.x + block_offset];
    nelem = nelems[blockId];
    ncolor = ncolors[blockId];
    int cur_offset = offset[blockId];
    color = colors + cur_offset;
    nelem2 = blockDim.x * (1 + (nelem - 1) / blockDim.x);
    ind_arg0_size = ind_arg0_sizes[blockId];
    ind_arg1_size = ind_arg1_sizes[blockId];
    ind_arg2_size = ind_arg2_sizes[blockId];
    ind_arg0_ptr = ind_arg0_ptrs + ind_arg0_offset[blockId];
    ind_arg1_ptr = ind_arg1_ptrs + ind_arg1_offset[blockId];
    ind_arg2_ptr = ind_arg2_ptrs + ind_arg2_offset[blockId];
    arg0 = arg0_d + cur_offset * 1;
    arg1 = arg1_d + cur_offset * 3;
    arg2_ptr = arg2_ptrs + cur_offset;
    arg3_ptr = arg3_ptrs + cur_offset;
    arg4_ptr = arg4_ptrs + cur_offset;
    arg5_ptr = arg5_ptrs + cur_offset;
    int nbytes = 0;
    ind_arg0_s = ((double *)(&shared[nbytes]));
    nbytes += ROUND_UP(ind_arg0_size * (sizeof(float ) * 3));
    ind_arg1_s = ((double *)(&shared[nbytes]));
    nbytes += ROUND_UP(ind_arg1_size * (sizeof(float ) * 3));
    ind_arg2_s = ((double *)(&shared[nbytes]));
  }
  __syncthreads();
  for (int n = threadIdx.x; n < ind_arg0_size; n += blockDim.x) {
    ind_arg0_s[0+n*3] = 0;
    ind_arg0_s[1+n*3] = 0;
    ind_arg0_s[2+n*3] = 0;
  }
  for (int n = threadIdx.x; n < ind_arg1_size; n += blockDim.x) {
    ind_arg1_s[0+n*3] = 0;
    ind_arg1_s[1+n*3] = 0;
    ind_arg1_s[2+n*3] = 0;
  }
  for (int n = threadIdx.x; n < ind_arg2_size; n += blockDim.x) {
    ind_arg2_s[0+n*3] = 0;
    ind_arg2_s[1+n*3] = 0;
    ind_arg2_s[2+n*3] = 0;
  }
  __syncthreads();
  for (int n = threadIdx.x; n < nelem2; n += blockDim.x) {
    int col2 = -1;
    if (n < nelem) {
      arg2_l[0] = 0;
      arg2_l[1] = 0;
      arg2_l[2] = 0;
      arg3_l[0] = 0;
      arg3_l[1] = 0;
      arg3_l[2] = 0;
      arg4_l[0] = 0;
      arg4_l[1] = 0;
      arg4_l[2] = 0;
      arg5_l[0] = 0;
      arg5_l[1] = 0;
      arg5_l[2] = 0;
      arg0_l[0] =  *(arg0 + (n * 1 + 0));
      arg1_l[0] =  *(arg1 + (n * 3 + 0));
      arg1_l[1] =  *(arg1 + (n * 3 + 1));
      arg1_l[2] =  *(arg1 + (n * 3 + 2));
      res(arg0_l,arg1_l,arg2_l,arg3_l,arg4_l,arg5_l);
       *(arg1 + (n * 3 + 0)) = arg1_l[0];
       *(arg1 + (n * 3 + 1)) = arg1_l[1];
       *(arg1 + (n * 3 + 2)) = arg1_l[2];
      col2 = color[n];
    }
    for (int col = 0; col < ncolor; ++col) {
      if (col == col2) {
        int ind_index = arg2_ptr[n];
        ind_arg0_s[0+ind_index*3] += arg2_l[0];
        ind_arg0_s[1+ind_index*3] += arg2_l[1];
        ind_arg0_s[2+ind_index*3] += arg2_l[2];
        ind_index = arg3_ptr[n];
        ind_arg1_s[0+ind_index*3] += arg3_l[0];
        ind_arg1_s[1+ind_index*3] += arg3_l[1];
        ind_arg1_s[2+ind_index*3] += arg3_l[2];
        ind_index = arg4_ptr[n];
        ind_arg2_s[0+ind_index*3] += arg4_l[0];
        ind_arg2_s[1+ind_index*3] += arg4_l[1];
        ind_arg2_s[2+ind_index*3] += arg4_l[2];
        ind_index = arg5_ptr[n];
        ind_arg2_s[0+ind_index*3] += arg5_l[0];
        ind_arg2_s[1+ind_index*3] += arg5_l[1];
        ind_arg2_s[2+ind_index*3] += arg5_l[2];
      }
      __syncthreads();
    }
  }
  for (int n = threadIdx.x; n < ind_arg0_size; n += blockDim.x) {
    int ind_index = ind_arg0_ptr[n];
    ind_arg0[0+ind_index*3] += ind_arg0_s[0+n*3];
    ind_arg0[1+ind_index*3] += ind_arg0_s[1+n*3];
    ind_arg0[2+ind_index*3] += ind_arg0_s[2+n*3];
  }
  for (int n = threadIdx.x; n < ind_arg1_size; n += blockDim.x) {
    int ind_index = ind_arg1_ptr[n];
    ind_arg1[0+ind_index*3] += ind_arg1_s[0+n*3];
    ind_arg1[1+ind_index*3] += ind_arg1_s[1+n*3];
    ind_arg1[2+ind_index*3] += ind_arg1_s[2+n*3];
  }
  for (int n = threadIdx.x; n < ind_arg2_size; n += blockDim.x) {
    int ind_index = ind_arg2_ptr[n];
    ind_arg2[0+ind_index*3] += ind_arg2_s[0+n*3];
    ind_arg2[1+ind_index*3] += ind_arg2_s[1+n*3];
    ind_arg2[2+ind_index*3] += ind_arg2_s[2+n*3];
  }
}

float op_par_loop_res(const char *name,op_set set,struct op_dat<void> *arg0,int idx0,op_ptr *ptr0,enum op_access acc0,struct op_dat<void> *arg1,int idx1,op_ptr *ptr1,enum op_access acc1,struct op_dat<void> *arg2,int idx2,op_ptr *ptr2,enum op_access acc2,struct op_dat<void> *arg3,int idx3,op_ptr *ptr3,enum op_access acc3,struct op_dat<void> *arg4,int idx4,op_ptr *ptr4,enum op_access acc4,struct op_dat<void> *arg5,int idx5,op_ptr *ptr5,enum op_access acc5)
{
  int nargs = 6;
  int ninds = 3;
  int gridsize = (set.size - 1) / BSIZE + 1;
  struct op_dat<void> args[6] = { *arg0,  *arg1,  *arg2,  *arg3,  *arg4,  *arg5};
  int idxs[6] = {-1, -1, idx2, idx3, idx4, idx5};
  op_ptr ptrs[6] = {OP_ID, OP_ID,  *ptr2,  *ptr3,  *ptr4,  *ptr5};
  int dims[6] = {arg0->dim, arg1->dim, arg2->dim, arg3->dim, arg4->dim, arg5->dim};
  enum op_access accs[6] = {acc0, acc1, acc2, acc3, acc4, acc5};
  int inds[6] = {-1, -1, 0, 1, 2, 2};
  op_plan *Plan = plan(name,set,nargs,args,idxs,ptrs,dims,accs,ninds,inds);
  int block_offset = 0;
  int reduct_bytes = 0;
  int reduct_size = 0;
  int reduct_shared = reduct_size * (BSIZE / 2);
  int const_bytes = 0;
  float total_time = 0.00000F;
  for (int col = 0; col < Plan->ncolors; ++col) {
    int nblocks = Plan->ncolblk[col];
    int nshared = Plan->nshared;
cudaEvent_t start, stop;
    float elapsed_time_ms = 0.00000F;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start,0);
    op_cuda_res<<<nblocks,BSIZE,nshared>>>(((double *)arg2->dat_d),Plan->ind_ptrs[0],Plan->ind_sizes[0],Plan->ind_offs[0],((double *)arg3->dat_d),Plan->ind_ptrs[1],Plan->ind_sizes[1],Plan->ind_offs[1],((double *)arg4->dat_d),Plan->ind_ptrs[2],Plan->ind_sizes[2],Plan->ind_offs[2],((struct contact *)arg0->dat_d),((double *)arg1->dat_d),Plan->ptrs[2],Plan->ptrs[3],Plan->ptrs[4],Plan->ptrs[5],block_offset,Plan->blkmap,Plan->offset,Plan->nelems,Plan->nthrcol,Plan->thrcol);
    cudaEventRecord(stop,0);
    cudaThreadSynchronize();
    cudaEventElapsedTime(&elapsed_time_ms,start,stop);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    total_time += elapsed_time_ms;
    cudaThreadSynchronize();
    block_offset += nblocks;
  }
  return total_time;
}
