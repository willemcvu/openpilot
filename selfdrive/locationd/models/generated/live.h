/******************************************************************************
 *                      Code generated with sympy 1.7.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_4823297603783166018);
void inv_err_fun(double *nom_x, double *true_x, double *out_1209893589484191652);
void H_mod_fun(double *state, double *out_5446294350774946742);
void f_fun(double *state, double dt, double *out_7890352782137742207);
void F_fun(double *state, double dt, double *out_6641875001366453990);
void h_3(double *state, double *unused, double *out_1847321279850234523);
void H_3(double *state, double *unused, double *out_7082165555105143204);
void h_4(double *state, double *unused, double *out_8536922137666975393);
void H_4(double *state, double *unused, double *out_8292369790379984057);
void h_9(double *state, double *unused, double *out_3828234516502945534);
void H_9(double *state, double *unused, double *out_3612517932160513364);
void h_10(double *state, double *unused, double *out_4156621338336708487);
void H_10(double *state, double *unused, double *out_5922521417932686327);
void h_12(double *state, double *unused, double *out_8247874486266393421);
void H_12(double *state, double *unused, double *out_3649461388605989862);
void h_31(double *state, double *unused, double *out_4590369118245091644);
void H_31(double *state, double *unused, double *out_6324551163749409683);
void h_32(double *state, double *unused, double *out_7471663611537064466);
void H_32(double *state, double *unused, double *out_1150389945245636849);
void h_13(double *state, double *unused, double *out_3210587242573528524);
void H_13(double *state, double *unused, double *out_6067144793180650559);
void h_14(double *state, double *unused, double *out_3828234516502945534);
void H_14(double *state, double *unused, double *out_3612517932160513364);
void h_19(double *state, double *unused, double *out_4278315153228771932);
void H_19(double *state, double *unused, double *out_8045197783262944545);
#define DIM 23
#define EDIM 22
#define MEDIM 22
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_3 = 3.841459;
void update_3(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_4 = 7.814728;
void update_4(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_9 = 7.814728;
void update_9(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_10 = 7.814728;
void update_10(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_12 = 7.814728;
void update_12(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_31 = 7.814728;
void update_31(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_32 = 9.487729;
void update_32(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_13 = 7.814728;
void update_13(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_14 = 7.814728;
void update_14(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_19 = 7.814728;
void update_19(double *, double *, double *, double *, double *);