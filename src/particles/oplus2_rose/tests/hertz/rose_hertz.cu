/*
 * OP2 instantiation of the hertz pairwise kernel.
 *
 * We take two datasets (from file) inputs (per-particle and contact data) and
 * expected output following a pairwise calculation.
 *
 * We use OP2 to express the problem as an iteration over contacts.
 *
 */
//#define AOS_EXTRA_DEBUG //< add (i,j) index information to struct
#ifdef GPU_TIMER
  #include "cuda_timer.h"
#elif POSIX_TIMER
  #include "posix_timer.h"
#else
  #include "simple_timer.h"
#endif
#include "check_result_vector.h"
#include "user_defined_types.h"
#include "op_seq.h"
#include "res.h"
#include "update.h"
#include "framework.h"
#include <sstream>
#include "kernels.h" 
using namespace std;

void run(struct params *input,int num_iter)
{
//--------------------
// One-time only costs
//--------------------
// Added by Nao for passing optional arguments
  int c;
  extern int OP_part_size;
  while((c = getopt((input -> params::argc),(input -> params::argv),"b:")) != (-1)){
    switch(c){
      case 'b':
{
        printf("found -b !!!\n");
        OP_part_size = atoi(optarg);
        break; 
      }
    }
  }
  printf("OP_part_size=%d\n",OP_part_size);
//AoS generate
  one_time. push_back (::SimpleTimer::SimpleTimer(("aos_gen")));
  one_time. back (). start ();
  struct contact *aos = new contact [(input -> params::nedge)];
  for (int e = 0; e < (input -> params::nedge); e++) {
    int i = (input -> params::edge)[e * 2];
    int j = (input -> params::edge)[(e * 2) + 1];
#ifdef AOS_EXTRA_DEBUG
#endif
    aos[e].contact::xi[0] = (input -> params::x)[i * 3];
    aos[e].contact::xi[1] = (input -> params::x)[(i * 3) + 1];
    aos[e].contact::xi[2] = (input -> params::x)[(i * 3) + 2];
    aos[e].contact::xj[0] = (input -> params::x)[j * 3];
    aos[e].contact::xj[1] = (input -> params::x)[(j * 3) + 1];
    aos[e].contact::xj[2] = (input -> params::x)[(j * 3) + 2];
    aos[e].contact::vi[0] = (input -> params::v)[i * 3];
    aos[e].contact::vi[1] = (input -> params::v)[(i * 3) + 1];
    aos[e].contact::vi[2] = (input -> params::v)[(i * 3) + 2];
    aos[e].contact::vj[0] = (input -> params::v)[j * 3];
    aos[e].contact::vj[1] = (input -> params::v)[(j * 3) + 1];
    aos[e].contact::vj[2] = (input -> params::v)[(j * 3) + 2];
    aos[e].contact::omegai[0] = (input -> params::omega)[i * 3];
    aos[e].contact::omegai[1] = (input -> params::omega)[(i * 3) + 1];
    aos[e].contact::omegai[2] = (input -> params::omega)[(i * 3) + 2];
    aos[e].contact::omegaj[0] = (input -> params::omega)[j * 3];
    aos[e].contact::omegaj[1] = (input -> params::omega)[(j * 3) + 1];
    aos[e].contact::omegaj[2] = (input -> params::omega)[(j * 3) + 2];
    aos[e].contact::radiusi = (input -> params::radius)[i];
    aos[e].contact::radiusj = (input -> params::radius)[j];
    aos[e].contact::massi = (input -> params::mass)[i];
    aos[e].contact::massj = (input -> params::mass)[j];
    aos[e].contact::typei = (input -> params::type)[i];
    aos[e].contact::typej = (input -> params::type)[j];
  }
  double *force_inc = new double [((input -> params::nnode) * 3)];
  double *force_dec = new double [((input -> params::nnode) * 3)];
  double *torque_dec = new double [((input -> params::nnode) * 3)];
  one_time. back (). stop_and_add_to_total ();
// OP2 init
  const char *fake_argv = "foo";
  one_time. push_back (::SimpleTimer::SimpleTimer(("op2_init")));
  one_time. back (). start ();
  op_init(1,((char **)fake_argv));
  one_time. back (). stop_and_add_to_total ();
// OP2 declarations of sets, ptrs and datasets
  one_time. push_back (::SimpleTimer::SimpleTimer(("op2_decl")));
  one_time. back (). start ();
  op_set nodes((input -> params::nnode),0L,"nodes");
  op_set edges((input -> params::nedge),(((input -> params::npartition) == 1)?0L : &input -> params::partition_length),"edges");
  op_ptr edge_map((edges),(nodes),2,(input -> params::edge),"edge_map");
  struct op_dat< contact  > p_aos((edges),1,aos,"p_aos");
  struct op_dat< double  > p_force_inc((nodes),3,force_inc,"p_force_inc");
  struct op_dat< double  > p_force_dec((nodes),3,force_dec,"p_force_dec");
  struct op_dat< double  > p_torque_dec((nodes),3,torque_dec,"p_torque_dec");
  struct op_dat< double  > p_shear((edges),3,(input -> params::shear),"p_shear");
  struct op_dat< double  > p_force((nodes),3,(input -> params::force),"p_force");
  struct op_dat< double  > p_torque((nodes),3,(input -> params::torque),"p_torque");
  one_time. back (). stop_and_add_to_total ();
// OP2 plan (see HACK below)
  one_time. push_back (::SimpleTimer::SimpleTimer(("op2_plan")));
  op_diagnostic_output();
//--------------------
// Per-iteration costs
//--------------------
  double res_time[2UL];
  per_iter. push_back (::SimpleTimer::SimpleTimer(("compute_kernel")));
  per_iter. push_back (::SimpleTimer::SimpleTimer(("add_kernel")));
  per_iter. push_back (::SimpleTimer::SimpleTimer(("result_fetch")));
  for (int run = 0; run < num_iter; run++) {
//-----------------------------------------------------------------------
//KERNEL INVOCATION
//res kernel computes delta values (force/torque) and shear
//update kernel produces force/torque results
    per_iter[0]. start ();
    op_par_loop_res("res",(edges),(struct op_dat<void> *)(&p_aos),0,0L,OP_READ,(struct op_dat<void> *)(&p_shear),0,0L,OP_RW,(struct op_dat<void> *)(&p_force_inc),0,&edge_map,OP_INC,(struct op_dat<void> *)(&p_force_dec),1,&edge_map,OP_INC,(struct op_dat<void> *)(&p_torque_dec),0,&edge_map,OP_INC,(struct op_dat<void> *)(&p_torque_dec),1,&edge_map,OP_INC);
    per_iter[0]. stop_and_add_to_total ();
//HACK for op2_plan one-time cost
//This cost is wrapped up into the first invocation of the res kernel 
//(run 0). We record the runtime for the first two invocations of the
//kernel and use the difference to approximate the plan cost. We also fixup
//the per_iter cost for the res kernel.
    if (run == 0) {
      res_time[0] = per_iter[0]. total_time ();
    }
    else if (run == 1) {
//total_time is a cumulative total so we subtract the run 0 time
      res_time[1] = (per_iter[0]. total_time () - res_time[0]);
      double op2_plan_time = (res_time[0] - res_time[1]);
      one_time. back (). set_total_time (op2_plan_time);
      per_iter[0]. set_total_time ((2 * res_time[1]));
    }
    per_iter[1]. start ();
    op_par_loop_update("update",(nodes),(struct op_dat<void> *)(&p_force_inc),0,0L,OP_RW,(struct op_dat<void> *)(&p_force_dec),0,0L,OP_RW,(struct op_dat<void> *)(&p_torque_dec),0,0L,OP_RW,(struct op_dat<void> *)(&p_force),0,0L,OP_RW,(struct op_dat<void> *)(&p_torque),0,0L,OP_RW);
    per_iter[1]. stop_and_add_to_total ();
//-----------------------------------------------------------------------
//POSTPROCESSING
//memcpy data back to host
    per_iter[2]. start ();
    op_fetch_data(&p_shear);
    op_fetch_data(&p_force);
    op_fetch_data(&p_torque);
    per_iter[2]. stop_and_add_to_total ();
//-----------------------------------------------------------------------
#if 0
//CHECKING
//only check results the first time around
#endif
  }
//cleanup
  :: delete []aos;
  delete []force_inc;
  delete []force_dec;
  delete []torque_dec;
}
