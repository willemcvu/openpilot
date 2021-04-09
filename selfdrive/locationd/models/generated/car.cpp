
extern "C"{

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}

}
extern "C" {
#include <math.h>
/******************************************************************************
 *                      Code generated with sympy 1.7.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_7194484907591779792) {
   out_7194484907591779792[0] = delta_x[0] + nom_x[0];
   out_7194484907591779792[1] = delta_x[1] + nom_x[1];
   out_7194484907591779792[2] = delta_x[2] + nom_x[2];
   out_7194484907591779792[3] = delta_x[3] + nom_x[3];
   out_7194484907591779792[4] = delta_x[4] + nom_x[4];
   out_7194484907591779792[5] = delta_x[5] + nom_x[5];
   out_7194484907591779792[6] = delta_x[6] + nom_x[6];
   out_7194484907591779792[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_2576233854144261535) {
   out_2576233854144261535[0] = -nom_x[0] + true_x[0];
   out_2576233854144261535[1] = -nom_x[1] + true_x[1];
   out_2576233854144261535[2] = -nom_x[2] + true_x[2];
   out_2576233854144261535[3] = -nom_x[3] + true_x[3];
   out_2576233854144261535[4] = -nom_x[4] + true_x[4];
   out_2576233854144261535[5] = -nom_x[5] + true_x[5];
   out_2576233854144261535[6] = -nom_x[6] + true_x[6];
   out_2576233854144261535[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_1210464128241612005) {
   out_1210464128241612005[0] = 1.0;
   out_1210464128241612005[1] = 0.0;
   out_1210464128241612005[2] = 0.0;
   out_1210464128241612005[3] = 0.0;
   out_1210464128241612005[4] = 0.0;
   out_1210464128241612005[5] = 0.0;
   out_1210464128241612005[6] = 0.0;
   out_1210464128241612005[7] = 0.0;
   out_1210464128241612005[8] = 0.0;
   out_1210464128241612005[9] = 1.0;
   out_1210464128241612005[10] = 0.0;
   out_1210464128241612005[11] = 0.0;
   out_1210464128241612005[12] = 0.0;
   out_1210464128241612005[13] = 0.0;
   out_1210464128241612005[14] = 0.0;
   out_1210464128241612005[15] = 0.0;
   out_1210464128241612005[16] = 0.0;
   out_1210464128241612005[17] = 0.0;
   out_1210464128241612005[18] = 1.0;
   out_1210464128241612005[19] = 0.0;
   out_1210464128241612005[20] = 0.0;
   out_1210464128241612005[21] = 0.0;
   out_1210464128241612005[22] = 0.0;
   out_1210464128241612005[23] = 0.0;
   out_1210464128241612005[24] = 0.0;
   out_1210464128241612005[25] = 0.0;
   out_1210464128241612005[26] = 0.0;
   out_1210464128241612005[27] = 1.0;
   out_1210464128241612005[28] = 0.0;
   out_1210464128241612005[29] = 0.0;
   out_1210464128241612005[30] = 0.0;
   out_1210464128241612005[31] = 0.0;
   out_1210464128241612005[32] = 0.0;
   out_1210464128241612005[33] = 0.0;
   out_1210464128241612005[34] = 0.0;
   out_1210464128241612005[35] = 0.0;
   out_1210464128241612005[36] = 1.0;
   out_1210464128241612005[37] = 0.0;
   out_1210464128241612005[38] = 0.0;
   out_1210464128241612005[39] = 0.0;
   out_1210464128241612005[40] = 0.0;
   out_1210464128241612005[41] = 0.0;
   out_1210464128241612005[42] = 0.0;
   out_1210464128241612005[43] = 0.0;
   out_1210464128241612005[44] = 0.0;
   out_1210464128241612005[45] = 1.0;
   out_1210464128241612005[46] = 0.0;
   out_1210464128241612005[47] = 0.0;
   out_1210464128241612005[48] = 0.0;
   out_1210464128241612005[49] = 0.0;
   out_1210464128241612005[50] = 0.0;
   out_1210464128241612005[51] = 0.0;
   out_1210464128241612005[52] = 0.0;
   out_1210464128241612005[53] = 0.0;
   out_1210464128241612005[54] = 1.0;
   out_1210464128241612005[55] = 0.0;
   out_1210464128241612005[56] = 0.0;
   out_1210464128241612005[57] = 0.0;
   out_1210464128241612005[58] = 0.0;
   out_1210464128241612005[59] = 0.0;
   out_1210464128241612005[60] = 0.0;
   out_1210464128241612005[61] = 0.0;
   out_1210464128241612005[62] = 0.0;
   out_1210464128241612005[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_4667158935505101972) {
   out_4667158935505101972[0] = state[0];
   out_4667158935505101972[1] = state[1];
   out_4667158935505101972[2] = state[2];
   out_4667158935505101972[3] = state[3];
   out_4667158935505101972[4] = state[4];
   out_4667158935505101972[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_4667158935505101972[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_4667158935505101972[7] = state[7];
}
void F_fun(double *state, double dt, double *out_8846265460316692692) {
   out_8846265460316692692[0] = 1;
   out_8846265460316692692[1] = 0;
   out_8846265460316692692[2] = 0;
   out_8846265460316692692[3] = 0;
   out_8846265460316692692[4] = 0;
   out_8846265460316692692[5] = 0;
   out_8846265460316692692[6] = 0;
   out_8846265460316692692[7] = 0;
   out_8846265460316692692[8] = 0;
   out_8846265460316692692[9] = 1;
   out_8846265460316692692[10] = 0;
   out_8846265460316692692[11] = 0;
   out_8846265460316692692[12] = 0;
   out_8846265460316692692[13] = 0;
   out_8846265460316692692[14] = 0;
   out_8846265460316692692[15] = 0;
   out_8846265460316692692[16] = 0;
   out_8846265460316692692[17] = 0;
   out_8846265460316692692[18] = 1;
   out_8846265460316692692[19] = 0;
   out_8846265460316692692[20] = 0;
   out_8846265460316692692[21] = 0;
   out_8846265460316692692[22] = 0;
   out_8846265460316692692[23] = 0;
   out_8846265460316692692[24] = 0;
   out_8846265460316692692[25] = 0;
   out_8846265460316692692[26] = 0;
   out_8846265460316692692[27] = 1;
   out_8846265460316692692[28] = 0;
   out_8846265460316692692[29] = 0;
   out_8846265460316692692[30] = 0;
   out_8846265460316692692[31] = 0;
   out_8846265460316692692[32] = 0;
   out_8846265460316692692[33] = 0;
   out_8846265460316692692[34] = 0;
   out_8846265460316692692[35] = 0;
   out_8846265460316692692[36] = 1;
   out_8846265460316692692[37] = 0;
   out_8846265460316692692[38] = 0;
   out_8846265460316692692[39] = 0;
   out_8846265460316692692[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_8846265460316692692[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_8846265460316692692[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8846265460316692692[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8846265460316692692[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_8846265460316692692[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_8846265460316692692[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_8846265460316692692[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_8846265460316692692[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_8846265460316692692[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_8846265460316692692[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8846265460316692692[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8846265460316692692[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_8846265460316692692[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_8846265460316692692[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_8846265460316692692[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8846265460316692692[56] = 0;
   out_8846265460316692692[57] = 0;
   out_8846265460316692692[58] = 0;
   out_8846265460316692692[59] = 0;
   out_8846265460316692692[60] = 0;
   out_8846265460316692692[61] = 0;
   out_8846265460316692692[62] = 0;
   out_8846265460316692692[63] = 1;
}
void h_25(double *state, double *unused, double *out_2358645640862431631) {
   out_2358645640862431631[0] = state[6];
}
void H_25(double *state, double *unused, double *out_6075262642728191814) {
   out_6075262642728191814[0] = 0;
   out_6075262642728191814[1] = 0;
   out_6075262642728191814[2] = 0;
   out_6075262642728191814[3] = 0;
   out_6075262642728191814[4] = 0;
   out_6075262642728191814[5] = 0;
   out_6075262642728191814[6] = 1;
   out_6075262642728191814[7] = 0;
}
void h_24(double *state, double *unused, double *out_7700269130267956316) {
   out_7700269130267956316[0] = state[4];
   out_7700269130267956316[1] = state[5];
}
void H_24(double *state, double *unused, double *out_8855632334565165841) {
   out_8855632334565165841[0] = 0;
   out_8855632334565165841[1] = 0;
   out_8855632334565165841[2] = 0;
   out_8855632334565165841[3] = 0;
   out_8855632334565165841[4] = 1;
   out_8855632334565165841[5] = 0;
   out_8855632334565165841[6] = 0;
   out_8855632334565165841[7] = 0;
   out_8855632334565165841[8] = 0;
   out_8855632334565165841[9] = 0;
   out_8855632334565165841[10] = 0;
   out_8855632334565165841[11] = 0;
   out_8855632334565165841[12] = 0;
   out_8855632334565165841[13] = 1;
   out_8855632334565165841[14] = 0;
   out_8855632334565165841[15] = 0;
}
void h_30(double *state, double *unused, double *out_4412189585487919305) {
   out_4412189585487919305[0] = state[4];
}
void H_30(double *state, double *unused, double *out_3538636268220991466) {
   out_3538636268220991466[0] = 0;
   out_3538636268220991466[1] = 0;
   out_3538636268220991466[2] = 0;
   out_3538636268220991466[3] = 0;
   out_3538636268220991466[4] = 1;
   out_3538636268220991466[5] = 0;
   out_3538636268220991466[6] = 0;
   out_3538636268220991466[7] = 0;
}
void h_26(double *state, double *unused, double *out_9032881220560786666) {
   out_9032881220560786666[0] = state[7];
}
void H_26(double *state, double *unused, double *out_6666888276042150529) {
   out_6666888276042150529[0] = 0;
   out_6666888276042150529[1] = 0;
   out_6666888276042150529[2] = 0;
   out_6666888276042150529[3] = 0;
   out_6666888276042150529[4] = 0;
   out_6666888276042150529[5] = 0;
   out_6666888276042150529[6] = 0;
   out_6666888276042150529[7] = 1;
}
void h_27(double *state, double *unused, double *out_668197806126750671) {
   out_668197806126750671[0] = state[3];
}
void H_27(double *state, double *unused, double *out_4826218256057616778) {
   out_4826218256057616778[0] = 0;
   out_4826218256057616778[1] = 0;
   out_4826218256057616778[2] = 0;
   out_4826218256057616778[3] = 1;
   out_4826218256057616778[4] = 0;
   out_4826218256057616778[5] = 0;
   out_4826218256057616778[6] = 0;
   out_4826218256057616778[7] = 0;
}
void h_29(double *state, double *unused, double *out_6823645001967410800) {
   out_6823645001967410800[0] = state[1];
}
void H_29(double *state, double *unused, double *out_7362435037349823419) {
   out_7362435037349823419[0] = 0;
   out_7362435037349823419[1] = 1;
   out_7362435037349823419[2] = 0;
   out_7362435037349823419[3] = 0;
   out_7362435037349823419[4] = 0;
   out_7362435037349823419[5] = 0;
   out_7362435037349823419[6] = 0;
   out_7362435037349823419[7] = 0;
}
void h_28(double *state, double *unused, double *out_4820792134792247001) {
   out_4820792134792247001[0] = state[5];
   out_4820792134792247001[1] = state[6];
}
void H_28(double *state, double *unused, double *out_1837799856328670975) {
   out_1837799856328670975[0] = 0;
   out_1837799856328670975[1] = 0;
   out_1837799856328670975[2] = 0;
   out_1837799856328670975[3] = 0;
   out_1837799856328670975[4] = 0;
   out_1837799856328670975[5] = 1;
   out_1837799856328670975[6] = 0;
   out_1837799856328670975[7] = 0;
   out_1837799856328670975[8] = 0;
   out_1837799856328670975[9] = 0;
   out_1837799856328670975[10] = 0;
   out_1837799856328670975[11] = 0;
   out_1837799856328670975[12] = 0;
   out_1837799856328670975[13] = 0;
   out_1837799856328670975[14] = 1;
   out_1837799856328670975[15] = 0;
}
}

extern "C"{
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_25 = 3.841459;
void update_25(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_24 = 5.991465;
void update_24(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_30 = 3.841459;
void update_30(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_26 = 3.841459;
void update_26(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_27 = 3.841459;
void update_27(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_29 = 3.841459;
void update_29(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_28 = 5.991465;
void update_28(double *, double *, double *, double *, double *);
}

#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}



extern "C"{

      void update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
      }
    
      void update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
      }
    
      void update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
      }
    
      void update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
      }
    
      void update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
      }
    
      void update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
      }
    
      void update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
      }
    
}
