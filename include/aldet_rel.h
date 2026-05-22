#ifndef ALDET_REL_H
#define ALDET_REL_H

class aldet_rel_data 
{
    public:
        //CI coefficients
        double ** coef_r;//el_s -- electron-state real part
        double ** coef_i;//el_s -- electron-state imaginary part
        
        
        //for storing multiple data
        int n_sets;
        int * n_states;
        
        //active space
        int n_act;
        int n_el;
        int print_number;
        
        //calculated variables
        long Nd;//number of all configurations
        int * f;
        int * vec;
        int * bit;
        int * buf;
        
        //bufer for energy
        double ** E_states;
        double ** S2;
        double ** L2;
        double ** P ;
        //sum of actent orbital energies
        double ** E_act;
        
        // arrays for storing results of multiplication of H and CI vectors
        double * Hv_r;    //actual
        double * Hv_i;    //actual
        double * Hv_buf_r;//bufer
        double * Hv_buf_i;//bufer
        
        // diagonal elements of H
        double * H_diag;
        double * H_diag_appr;
        
        int import_done;
        double   E_core;    
        double * F_act_r;    
        double * F_act_i;    
        double * act_INTS_AA;
        double * act_INTS_AB;
        double * act_INTS_BB;
        
        int do_PT;
        double * T3_AAA;
        double * T3_AAB;
        double * T3_BBA;
        double * T3_BBB;
        double * T2_AA ;
        double * T2_AB ;
        double * T2_BB ;
        double * T1_A  ;
        double * T1_B  ;
        double   T0    ;
        
        double   Lambda_core;
        double * Lambda_act;
        int * act_rep_num;
        
        
        double * J_act;
        double * K_act;
        
        
        
        int n_e1;
        int * __restrict__ e1_ind ;  
        int * __restrict__ e1_orbs;  
        
        double * __restrict__ e1_sign;  
                
        int n_e1_so;
        int * __restrict__ e1_ind_so;  
        double * __restrict__ e1_h_so;  
        
        int n_e2;
        int * __restrict__ e2_orbs;  
        int * __restrict__ e2_ind ;  
        double * __restrict__ e2_V;  
        double * __restrict__ e2_sign;  
        
        
        int n_e3;
        int * __restrict__ e3_ind ;  
        double * __restrict__ e3_V;  
        
        
        double* spin_sign;
        
        aldet_rel_data(); 
        int get_dim(int ext_n_act, int ext_na, int ext_nb, int ext_n_sets, int ext_print_number);
        
        int simple_import_data(double * ext_act_INTS,
                               double * ext_act_INTS_AB,
                               double * ext_F_act,
                               double   ext_E_core);
        
        int U_simple_import_data(double * ext_act_INTS_AA,
                                 double * ext_act_INTS_AB,
                                 double * ext_act_INTS_BB,
                                 double * ext_F_act_A,
                                 double * ext_F_act_B,
                                 double   ext_E_core);
        
        int PT2_alloc();
        
        int PT2_import_data(double * ext_T3,
                            double * ext_T3_AB,
                            double * ext_T2,
                            double * ext_T2_AB,
                            double * ext_T1,
                            double   ext_T0);
        
        int UPT2_import_data(double * ext_T3_AAA,
                             double * ext_T3_AAB,
                             double * ext_T3_BBA,
                             double * ext_T3_BBB,
                             double * ext_T2_AA,
                             double * ext_T2_AB,
                             double * ext_T2_BB,
                             double * ext_T1_A,
                             double * ext_T1_B,
                             double   ext_T0);
        
        
        
        int gen_ext_ind();
        int gen_ext_ind_PT();
        int gen_ext_ind_SO();
        
        int PT_update();
        int SO_update(double * Sx, double * Sy, double * Sz);
        
        int calc_spin_matr(double * S2, double * Sz, double * C_r, double * C_i,int n_s, int ld);
        int init_zero_vec(int n_s, int i_set);
        
        int printf_occ(int i_CI);
        int print_states(int i_s, int n_s, int print);        
        
        int H_full_calc_and_diag(int i_set);
        int H_full_calc  (double * H );
        int H_full_calc_i(double * Hi);
        
        ~aldet_rel_data();

        
};

int compress_H_full_to_H_ab(double * H, double * Hf, aldet_data * CI, aldet_rel_data * CI_r);

#endif
 
