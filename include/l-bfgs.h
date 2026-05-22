#ifndef __LBFGS_H
#define __LBFGS_H

extern int lbfgs_step_num;
extern double * prev_grad_delta;
extern double * prev_grad;
extern double * prev_step;
extern double * prev_skal_prod;
extern double max_lenth;
extern int lbfgs_vec_num;
int lbfgs_init(int , int);
int lbfgs_finalize();
double lbfgs_step(double *, double *, double *, int (*g_calc)(int, double *, double * ),int (*inv_h_calc)(int, double *, double * ),int,double);
int do_nothing(int, double *, double * );
#endif
