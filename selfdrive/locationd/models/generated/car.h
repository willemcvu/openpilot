/******************************************************************************
 *                      Code generated with sympy 1.7.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_7194484907591779792);
void inv_err_fun(double *nom_x, double *true_x, double *out_2576233854144261535);
void H_mod_fun(double *state, double *out_1210464128241612005);
void f_fun(double *state, double dt, double *out_4667158935505101972);
void F_fun(double *state, double dt, double *out_8846265460316692692);
void h_25(double *state, double *unused, double *out_2358645640862431631);
void H_25(double *state, double *unused, double *out_6075262642728191814);
void h_24(double *state, double *unused, double *out_7700269130267956316);
void H_24(double *state, double *unused, double *out_8855632334565165841);
void h_30(double *state, double *unused, double *out_4412189585487919305);
void H_30(double *state, double *unused, double *out_3538636268220991466);
void h_26(double *state, double *unused, double *out_9032881220560786666);
void H_26(double *state, double *unused, double *out_6666888276042150529);
void h_27(double *state, double *unused, double *out_668197806126750671);
void H_27(double *state, double *unused, double *out_4826218256057616778);
void h_29(double *state, double *unused, double *out_6823645001967410800);
void H_29(double *state, double *unused, double *out_7362435037349823419);
void h_28(double *state, double *unused, double *out_4820792134792247001);
void H_28(double *state, double *unused, double *out_1837799856328670975);
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
void set_mass(double x);

void set_rotational_inertia(double x);

void set_center_to_front(double x);

void set_center_to_rear(double x);

void set_stiffness_front(double x);

void set_stiffness_rear(double x);
