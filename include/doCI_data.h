#ifndef dCd_H
#define dCd_H

# include "aldet.h"
# include "CI.h"

double extern svd_eps;

class doCI_data;

class doCI_data 
{
    public:
        
        molecule * A;
        molecule * B;
        int * n_ao      ;
        int * n_alp_el  ;
        int * n_bet_el  ;
        int * n_act     ;
        int **n_alp_el_f;
        int **n_bet_el_f;
        int **n_act_f   ;
        int * n_cor     ;
        
        int n_svd;
        int n_zero_svd;
        int n_frag;
        
        double **VEC;
        
        double * S_DD ;
        double * L_MO ;
        double * R_MO ;
        double * L_COR;
        double * R_COR;
        double * L_ACT;
        double * R_ACT;
        double * L_COR_1st_order;
        double * R_COR_1st_order;
        double * L_ACT_1st_order;
        double * R_ACT_1st_order;
        double * SVD  ;
        double * BUF  ;
        double * S_SV ;
        double ** L_ACT_F;
        double ** R_ACT_F;
        double * S_FV ;
        double * U_ACT;
        double * V_ACT;
        double * S_SA ;
        double * S_SB ;
        double * DM_C ;//core
        double * DM_R_1st;//core
        double * DM_L_1st;//core
        double ** DM_C_F ;//core by fragment
        double * DM_A;//alp active only
        double * DM_B;//bet active only
        double * DM_T ;//total for 1-el ints
        double * DM_TA;//total for 1-el ints alp
        double * DM_TB;//total for 1-el ints bet
        
        double  E_core;
        double Dx_core;
        double Dy_core;
        double Dz_core;        
        
        double *  H_ACT;
        double * Dx_ACT;
        double * Dy_ACT;
        double * Dz_ACT;        
        
            
        int * n_st_total;
        int ** n_states;
        ci_map ** ci1_f;
        ci_map ** ci2_f;
        ci_map_arr ci1;
        ci_map_arr ci2;
        
        double * ci_aldet1;int ld_ci1;long n_det1;//delete
        double * ci_aldet2;int ld_ci2;long n_det2;//delete
        
        aldet_data aldet;
        int two_charge_states;
        aldet_data aldet_2;
        
        
        double * h_2e;
        double * ACT_DET_A;
        double * ACT_ADJ_A;
        double * ACT_ADJ_B;
        double * J;
        double * K;
        double * K_R_1st;
        double * K_L_1st;
        double * act_INTS;
        double * act_INTS_1st_order;
        
        double p_SVD;
        
        ci_pr_vec_aldet ci1_lt_aldet;
        ci_pr_vec_aldet ci2_lt_aldet;
    
        double * ACT_DMT_A;
        double * ACT_DMT_B;
        double * ACT_DMT_A_1st_order;
        double * ACT_DMT_B_1st_order;
        double * S_new;
        double * H_2el;
    
        
        int first_alloc(molecule *A, molecule * B);
        
        int gen_aldet_data();
        
        int calc_S_MO(molecule * D);///to be deleted 
        
        int AO_to_MO(double * M);
        
        int cor_svd(molecule * D);
        
        int cor_svd_PT(molecule * D);
        
//         int calc_DM_D_MO(double * D);
        
        int DM_gen_ortogonal(int s, double coef);
        
        int second_alloc();
        
        int act_ort(molecule * D);
        
        int AO_to_act(double * M_act, double *M_AO);
        
        int cpy_VEC_to_ACT();
        
        int cpy_ACT_VEC_to_ACT(double * ACT);
        
        int calc_S_SV(molecule * D);
        
        int ci_ortogonalization(double * S_MO);
        
        int read_act_ints(char * prefix);
        
        int write_act_ints(char * prefix);
        
        int calc_CI_for_dec_wf();
        
        int calc_CI_for_PT();
        
        int calc_DM_for_dec_wf();
        
        int calc_DM_for_PT();
        
        int calc_S_for_dec_wf();
        
        int ci_link();
        
        int link_and_transform_for_separate_spaces(molecule * cD);
        
        int link_and_transform_for_PT(molecule * cD);
        
        int S_calc();
        
        int calc_1el_DM();
        
        int average_DM(std::vector<double> avecoe);
        
        int nat_orb_calc_1mol(int i, int renorm);
        
        int transform_1el_DM(int N);
        
        int V_2el_calc();
        
        double H_2el_calc(double S);
        
        int clear();// to be changed to the destructor!
        

    private:
        
};


#endif
