static void update(
  double *force_inc, double *force_dec,
  double *torque_dec,
  double *force, 
  double *torque) {
  //add/subtract in delta values
  force[0] += force_inc[0];
  force[1] += force_inc[1];
  force[2] += force_inc[2];

  force[0] -= force_dec[0];
  force[1] -= force_dec[1];
  force[2] -= force_dec[2];

  torque[0] -= torque_dec[0];
  torque[1] -= torque_dec[1];
  torque[2] -= torque_dec[2];

  //zero out intermediate delta values
  force_inc[0] = 0.0;
  force_inc[1] = 0.0;
  force_inc[2] = 0.0;

  force_dec[0] = 0.0;
  force_dec[1] = 0.0;
  force_dec[2] = 0.0;

  torque_dec[0] = 0.0;
  torque_dec[1] = 0.0;
  torque_dec[2] = 0.0;
}
