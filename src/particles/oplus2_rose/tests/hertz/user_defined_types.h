// This is where we put the user defined types;

#ifndef _USER_DEFINED_TYPES_
#define _USER_DEFINED_TYPES_

// AoS read-dataset for device
struct contact {
  int i; int j;
  double xi[3];     double xj[3];
  double vi[3];     double vj[3];
  double omegai[3]; double omegaj[3];
  double radiusi;   double radiusj;
  double massi;     double massj;
  int    typei;     int typej;
};

#endif
