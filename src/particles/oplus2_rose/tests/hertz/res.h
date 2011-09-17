#include <math.h>
#define sqrtFiveOverSix 0.91287092917527685576161630466800355658790782499663875

static void res(struct contact *aos, 
  double *shear, 
  double *force_inc, double *force_dec, 
  double *torquei, double *torquej) {
  //TODO: don't hardcode, push these into constant memory
  double dt = 0.00001;
  double nktv2p = 1;
  double yeff = 3134796.2382445144467056;
  double geff = 556173.5261401557363570;
  double betaeff = -0.3578571305033167;
  double coeffFrict = 0.5;

  // del is the vector from j to i
  double delx = aos->xi[0] - aos->xj[0];
  double dely = aos->xi[1] - aos->xj[1];
  double delz = aos->xi[2] - aos->xj[2];

  double rsq = delx*delx + dely*dely + delz*delz;
  double radsum = aos->radiusi + aos->radiusj;
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
    double vr1 = aos->vi[0] - aos->vj[0];
    double vr2 = aos->vi[1] - aos->vj[1];
    double vr3 = aos->vi[2] - aos->vj[2];

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
    double wr1 = (aos->radiusi*aos->omegai[0] + aos->radiusj*aos->omegaj[0]) * rinv;
    double wr2 = (aos->radiusi*aos->omegai[1] + aos->radiusj*aos->omegaj[1]) * rinv;
    double wr3 = (aos->radiusi*aos->omegai[2] + aos->radiusj*aos->omegaj[2]) * rinv;

    // normal forces = Hookian contact + normal velocity damping
    double meff = aos->massi*aos->massj/(aos->massi+aos->massj);
    //not-implemented: freeze_group_bit

    double deltan = radsum-r;

    //derive contact model parameters (inlined)
    //yeff, geff, betaeff, coeffFrict are constant lookup tables
    double reff = aos->radiusi * aos->radiusj / (aos->radiusi + aos->radiusj);
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
    double fs1 = - (kt*shear[0] + gammat*vtr1);
    double fs2 = - (kt*shear[1] + gammat*vtr2);
    double fs3 = - (kt*shear[2] + gammat*vtr3);

    // rescale frictional displacements and forces if needed
    double fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
    double fn = xmu * fabs(ccel*r);
    double shrmag = 0;
    if (fs > fn) {
      shrmag = sqrt(shear[0]*shear[0] +
          shear[1]*shear[1] +
          shear[2]*shear[2]);
      if (shrmag != 0.0) {
        shear[0] = (fn/fs) * (shear[0] + gammat*vtr1/kt) - gammat*vtr1/kt;
        shear[1] = (fn/fs) * (shear[1] + gammat*vtr2/kt) - gammat*vtr2/kt;
        shear[2] = (fn/fs) * (shear[2] + gammat*vtr3/kt) - gammat*vtr3/kt;
        fs1 *= fn/fs;
        fs2 *= fn/fs;
        fs3 *= fn/fs;
      } else {
        fs1 = fs2 = fs3 = 0.0;
      }
    }

    double fx = delx*ccel + fs1;
    double fy = dely*ccel + fs2;
    double fz = delz*ccel + fs3;

    double tor1 = rinv * (dely*fs3 - delz*fs2);
    double tor2 = rinv * (delz*fs1 - delx*fs3);
    double tor3 = rinv * (delx*fs2 - dely*fs1);

    // this is what we've been working up to!
    force_inc[0] += fx;
    force_inc[1] += fy;
    force_inc[2] += fz;

    force_dec[0] += fx;
    force_dec[1] += fy;
    force_dec[2] += fz;

    torquei[0] += -(aos->radiusi*tor1);
    torquei[1] += -(aos->radiusi*tor2);
    torquei[2] += -(aos->radiusi*tor3);

    torquej[0] += -(aos->radiusj*tor1);
    torquej[1] += -(aos->radiusj*tor2);
    torquej[2] += -(aos->radiusj*tor3);
  }
}
