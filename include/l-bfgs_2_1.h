#ifndef __LBFGS_H
#define __LBFGS_H


class l_bfgs_engine;

class l_bfgs_engine 
{
    public:
        int lbfgs_step_num;
        double * prev_grad_delta;
        double * prev_grad;
        double * prev_step;
        double * prev_skal_prod;
        double max_lenth;
        int lbfgs_vec_num;
        
        l_bfgs_engine();
        int init(int , int);
        double step(double *, double *, double *, int (*g_calc)(int, double *, double * ),int (*inv_h_calc)(int, double *, double * ),int,double);
        ~l_bfgs_engine ();
        
};

class bfgs_engine;

class bfgs_engine 
{
    public:
        int dim;
        int bfgs_step_num;
        double * grad_delta;
        double * prev_grad;
        double * prev_step;
        double * inv_hess;
        double * I_m_rsy;
        double * B;
        
        double max_lenth;
//         double * prev_skal_prod;
//         int lbfgs_vec_num;
        
        bfgs_engine();
        int init(int);
        int set_eq_diag(double c);
        int set_diag(double * c);
        
        double step(double *, double *, int (*g_calc)(int, double *, double * ));
        ~bfgs_engine ();
        
};

int do_nothing(int, double *, double * );



#endif
