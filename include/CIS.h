#ifndef __cis
#define __cis

//# include "molecule.h"
//# include "molecule2.h"
//# include "davidson.h"
# include "inp_par_read.h"


int CIS(molecule * M, cis_par * cis, char * job_name);

// double two_el_int(const double* B, int N, int n_aux, int p, int q, int r, int s);

class CIS_engine;

class CIS_engine{
    public:
      
      molecule* M;
      long n_ao;
      int method;
      int n_f;
      long n_c;
      long n_v;
      long n_aux;
      int n_s;
      int n_b;
      
      double* H;
      double* DM;
      
      int* nums;
      
      //double* A;
      double* B;
      double* F;
      double* H_CIS; //новое
      double* A1;
      double* C;
      double* d;
      double* E_app;
      double* E_p;
      double* H_p;
      double C1;
      
      CIS_engine();
      int set_par(molecule* ext_M, int ext_method, int ext_n_f, long ext_n_c, long ext_n_v, long ext_n_aux, int ext_n_s, int ext_n_b);
      int calc_H_CIS(); //новое
      int H_mult(int);
      int calc_MO_INTS();
      int calc_F();
      int gen_bf();
      int gen_H();
      int calc_evec();
      int print_excited_states_info();
      int gen_NO(std::vector <double> w_state, char* job_name);
      ~CIS_engine();
};

#endif
