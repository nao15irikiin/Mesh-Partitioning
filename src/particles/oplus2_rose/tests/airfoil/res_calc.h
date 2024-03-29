void res_calc(float *x1,  float *x2,  float *q1,  float *q2,
              float *adt1,float *adt2,float *res1,float *res2, int *boun) {
  float dx,dy,mu, ri, p1,vol1, p2,vol2, f, fac;

  if (*boun)
    fac = 0.0;       // boundary face
  else
    fac = 1.0;       // regular face

  dx = x1[0] - x2[0];
  dy = x1[1] - x2[1];

  ri   =  1.0/q1[0];
  p1   = gm1*(q1[3]-0.5*ri*(q1[1]*q1[1]+q1[2]*q1[2]));
  vol1 =  fac*ri*(q1[1]*dy - q1[2]*dx);

  ri   =  1.0/q2[0];
  p2   = gm1*(q2[3]-0.5*ri*(q2[1]*q2[1]+q2[2]*q2[2]));
  vol2 =  fac*ri*(q2[1]*dy - q2[2]*dx);

  if (*boun) {
    vol1 = 0.0;
    vol2 = 0.0;
  }

  mu = 0.5*((*adt1)+(*adt2))*eps;

  f = 0.5*( vol1* q1[0]         + vol2* q2[0]         ) + mu*(q1[0]-q2[0]);
  res1[0] += f*fac;
  res2[0] -= f;
  f = 0.5*( vol1* q1[1] + p1*dy + vol2* q2[1] + p2*dy ) + mu*(q1[1]-q2[1]);
  res1[1] += f*fac;
  res2[1] -= f;
  f = 0.5*( vol1* q1[2] - p1*dx + vol2* q2[2] - p2*dx ) + mu*(q1[2]-q2[2]);
  res1[2] += f*fac;
  res2[2] -= f;
  f = 0.5*( vol1*(q1[3]+p1)     + vol2*(q2[3]+p2)     ) + mu*(q1[3]-q2[3]);
  res1[3] += f*fac;
  res2[3] -= f;
}
