#ifndef RI_H
#define RI_H
#include "molecule.h"


int gen_RI_AA(molecule * A);

int regen_RI_AA(molecule * A);

int RI_core_realloc(long n_c, long n_ao);

int DM_to_F_transform_RI(double * J, double * K, double * DM, double * MO, long n_c, long n_ao);

int calc_2el_MO_INTS_RI(long n_ao, double * DM_D, double * J, double * K, double * MO_INTS, 
                        double * MO_VEC, double * MO_VEC1, double * MO_VEC2, 
                        long n_c, long n_a, long t);

int calc_2el_MO_INTS_AAAG_RI(double * MO_INTS, double * MO_VEC, long n_c, long n_a, long n_ao);

int calc_SOSCF_hess_INTS_RI (double * aa_J_ints, 
                             double * aa_K_ints,
                             double * cv_J_ints,
                             double * cv_K_ints, double * MO_VEC, long n_c, long n_a, long n_ao);

class RI_data;

class RI_data 
{
    public:
        long n_cor;
        long n_act;
        long n_vac;
        long n_ao;
        long aux_n_ao;
        
        double * VC_RI_M;
        double * CA_RI_M;
        double * CC_RI_M;
        double * VA_RI_M;
        double * AA_RI_M;
        
        double * COR_CVEC;
        double * ACT_CVEC;
        double * VAC_CVEC;

        std::vector <Shell> s;
        
        RI_data();
        int set_par(molecule * A, long ext_n_cor, long ext_n_act, long ext_n_vac);
        int MO_calc(double * VEC, long ext_n_ao);
        int VAAA_calc(double * VAAA);
        int CAAA_calc(double * CAAA);
        int VACC_calc(double * VACC_J, double * VACC_K);
        int CACC_calc(double * CACC_J, double * CACC_K);
        
        int M_calc(double ** M, std::vector<Shell> shells,
                   double * V1, double * V2, long n1, long n2, char * f_name);
                   
        std::vector<Shell> aux_shells;
//         double* AUX_Vm05;
        ~RI_data();
        
    private:
//         double* AUX_V;
//         double* AUX_Vp05;
        
        
};

class RI_complex_data;

class RI_complex_data 
{
    public:
        long n_cor;
        long n_act;
        long n_vac;
        long n_ao;
        long aux_n_ao;
        
        double * VC_RI_M_r;
        double * CA_RI_M_r;
        double * CC_RI_M_r;
        double * VA_RI_M_r;
        double * AA_RI_M_r;
        
        double * VC_RI_M_i;
        double * CA_RI_M_i;
        double * CC_RI_M_i;
        double * VA_RI_M_i;
        double * AA_RI_M_i;
        
        
        RI_complex_data();
        int set_par(molecule * A, long ext_n_cor, long ext_n_act, long ext_n_vac);
        int MO_calc(double * VEC_r, double * VEC_i, long ext_n_ao);
//         int VAAA_calc(double * VAAA);
//         int CAAA_calc(double * CAAA);
//         int VACC_calc(double * VACC_J, double * VACC_K);
//         int CACC_calc(double * CACC_J, double * CACC_K);
//         
//         int M_calc(double ** M, std::vector<Shell> shells,
//                    double * V1, double * V2, long n1, long n2, char * f_name);
//                    
//         std::vector<Shell> aux_shells;
// //         double* AUX_Vm05;
        ~RI_complex_data();
        
    private:
//         double* AUX_V;
//         double* AUX_Vp05;
        
        
};



#endif
