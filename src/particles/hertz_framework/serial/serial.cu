/*
 * Serial implementation of the hertz pairwise kernel.
 *
 */

#ifdef GPU_TIMER
  #include "cuda_timer.h"
#elif POSIX_TIMER
  #include "posix_timer.h"
#else
  #include "simple_timer.h"
#endif

#include "check_result_vector.h"
#include "framework.h"
#include "hertz_constants.h"
#include <math.h>
#include <sstream>

using namespace std;

double dt;
double nktv2p;
double yeff;
double geff;
double betaeff;
double coeffFrict;

inline void pair_interaction(
    int i, int j,
  //inputs
    double *xi, double *xj,           //position
    double *vi, double *vj,           //velocity
    double *omegai, double *omegaj,   //rotational velocity
    double radi, double radj,         //radius
    double massi, double massj,       //mass
    int typei, int typej,             //type
  //inouts
    double *shear,
    double *forcei,
    double *forcej,
    double *torquei,
    double *torquej) {

  // del is the vector from j to i
  double delx = xi[0] - xj[0];
  double dely = xi[1] - xj[1];
  double delz = xi[2] - xj[2];

  double rsq = delx*delx + dely*dely + delz*delz;
  double radsum = radi + radj;
  if (rsq >= radsum*radsum) {
    //unset non-touching atoms
    shear[0] = 0.0;
    shear[1] = 0.0;
    shear[2] = 0.0;
  } else {
    //distance between centres of atoms i and j
    //or, magnitude of del vector
    double r = sqrt(rsq);
    double rinv = 1.0/r;
    double rsqinv = 1.0/rsq;

    // relative translational velocity
    double vr1 = vi[0] - vj[0];
    double vr2 = vi[1] - vj[1];
    double vr3 = vi[2] - vj[2];

    // normal component
    double vnnr = vr1*delx + vr2*dely + vr3*delz;
    double vn1 = delx*vnnr * rsqinv;
    double vn2 = dely*vnnr * rsqinv;
    double vn3 = delz*vnnr * rsqinv;

    // tangential component
    double vt1 = vr1 - vn1;
    double vt2 = vr2 - vn2;
    double vt3 = vr3 - vn3;

    // relative rotational velocity
    double deltan = radsum-r;
    double cri = radi-0.5*deltan;
    double crj = radj-0.5*deltan;
    double wr1 = (cri*omegai[0] + crj*omegaj[0]) * rinv;
    double wr2 = (cri*omegai[1] + crj*omegaj[1]) * rinv;
    double wr3 = (cri*omegai[2] + crj*omegaj[2]) * rinv;

    // normal forces = Hookian contact + normal velocity damping
    double meff = massi*massj/(massi+massj);
    //not-implemented: freeze_group_bit

    //derive contact model parameters (inlined)
    //Yeff, Geff, betaeff, coeffFrict are lookup tables
    double reff = radi * radj / (radi + radj);
    double sqrtval = sqrt(reff * deltan);
    double Sn = 2.    * yeff * sqrtval;
    double St = 8.    * geff * sqrtval;
    double kn = 4./3. * yeff * sqrtval;
    double kt = St;
    double gamman=-2.*sqrtFiveOverSix*betaeff*sqrt(Sn*meff);
    double gammat=-2.*sqrtFiveOverSix*betaeff*sqrt(St*meff);
    double xmu=coeffFrict;
    //not-implemented if (dampflag == 0) gammat = 0;
    kn /= nktv2p;
    kt /= nktv2p;

    double damp = gamman*vnnr*rsqinv;
	  double ccel = kn*(radsum-r)*rinv - damp;

    //not-implemented cohesionflag

    // relative velocities
    double vtr1 = vt1 - (delz*wr2-dely*wr3);
    double vtr2 = vt2 - (delx*wr3-delz*wr1);
    double vtr3 = vt3 - (dely*wr1-delx*wr2);

    // shear history effects
    shear[0] += vtr1 * dt;
    shear[1] += vtr2 * dt;
    shear[2] += vtr3 * dt;

    // rotate shear displacements
    double rsht = shear[0]*delx + shear[1]*dely + shear[2]*delz;
    rsht *= rsqinv;

    shear[0] -= rsht*delx;
    shear[1] -= rsht*dely;
    shear[2] -= rsht*delz;

    // tangential forces = shear + tangential velocity damping
    double fs1 = - (kt*shear[0]);
    double fs2 = - (kt*shear[1]);
    double fs3 = - (kt*shear[2]);

    // rescale frictional displacements and forces if needed
    double fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
    double fn = xmu * fabs(ccel*r);
    double shrmag = sqrt(shear[0]*shear[0] +
                         shear[1]*shear[1] +
                         shear[2]*shear[2]);
    if (fs > fn) {
      if (shrmag != 0.0) {
        fs1 *= fn/fs;
        fs2 *= fn/fs;
        fs3 *= fn/fs;
        shear[0] = -fs1/kt;
        shear[1] = -fs2/kt;
        shear[2] = -fs3/kt;
      } else {
        fs1 = 0.0;
        fs2 = 0.0;
        fs3 = 0.0;
      }
    } else {
      fs1 -= (gammat*vtr1);
      fs2 -= (gammat*vtr2);
      fs3 -= (gammat*vtr3);
    }

    double fx = delx*ccel + fs1;
    double fy = dely*ccel + fs2;
    double fz = delz*ccel + fs3;

    double tor1 = rinv * (dely*fs3 - delz*fs2);
    double tor2 = rinv * (delz*fs1 - delx*fs3);
    double tor3 = rinv * (delx*fs2 - dely*fs1);

    // this is what we've been working up to!
    forcei[0] += fx;
    forcei[1] += fy;
    forcei[2] += fz;

    forcej[0] -= fx;
    forcej[1] -= fy;
    forcej[2] -= fz;

    torquei[0] -= cri*tor1;
    torquei[1] -= cri*tor2;
    torquei[2] -= cri*tor3;

    torquej[0] -= crj*tor1;
    torquej[1] -= crj*tor2;
    torquej[2] -= crj*tor3;
  }
}

void run(struct params *input, int num_iter) {

  //--------------------
  // Per-iteration costs
  //--------------------

  per_iter.push_back(SimpleTimer("kernel"));

  //internal copies of outputs
  double *force = new double[input->nnode*3];
  double *torque = new double[input->nnode*3];
  double *shear = new double[input->nedge*3];

  for (int run=0; run<num_iter; run++) {
    //setup constants
    dt = input->dt;
    nktv2p = input->nktv2p;
    yeff = input->yeff[3];
    geff = input->geff[3];
    betaeff = input->betaeff[3];
    coeffFrict = input->coeffFrict[3];

    //make some internal copies
    double *force = new double[input->nnode*3];
    double *torque = new double[input->nnode*3];
    double *shear = new double[input->nedge*3];

    //make copies
    for (int n=0; n<input->nnode*3; n++) {
      force[n] = input->force[n];
      torque[n] = input->torque[n];
    }
    for (int e=0; e<input->nedge*3; e++) {
      shear[e] = input->shear[e];
    }

    per_iter[0].start();
    for (int e=0; e<input->nedge; e++) {
      int i = input->edge[(e*2)];
      int j = input->edge[(e*2)+1];
      pair_interaction(
        i, j,
        &input->x[(i*3)],     &input->x[(j*3)],
        &input->v[(i*3)],     &input->v[(j*3)],
        &input->omega[(i*3)], &input->omega[(j*3)],
         input->radius[i],     input->radius[j],
         input->mass[i],       input->mass[j],
         input->type[i],       input->type[j],
        &shear[(e*3)],
        &force[(i*3)],        &force[(j*3)],
        &torque[(i*3)],       &torque[(j*3)]
      );
    }
    per_iter[0].stop_and_add_to_total();

    //only check results the first time around
    if (run == 0) {
      for (int n=0; n<input->nnode; n++) {
        std::stringstream out;
        out << "force[" << n << "]";
        check_result_vector(
            out.str().c_str(),
            &input->expected_force[(n*3)], &force[(n*3)]);
        out.str("");

        out << "torque[" << n << "]";
        check_result_vector(
            out.str().c_str(),
            &input->expected_torque[(n*3)], &torque[(n*3)]);
      }
      for (int n=0; n<input->nedge; n++) {
        stringstream out;
        out << "shear[" << n << "]";
        check_result_vector(
            out.str().c_str(),
            &input->expected_shear[(n*3)], &shear[(n*3)]);
      }
    }

  }

  delete[] force;
  delete[] torque;
  delete[] shear;
}
