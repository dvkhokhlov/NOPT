#ifndef RES_FIT_H
#define RES_FIT_H

# include "RI.h"
# include "PT_tensors.h"


class UPT_tensors;

class UPT_tensors 
{
    public:
        
        RI_data * RI_A;
        RI_data * RI_B;
        
        double * e_c_A;
        double * e_a_A;
        double * e_IP_A;
        double * e_EA_A;
        double * e_v_A;
        
        double * e_c_B;
        double * e_a_B;
        double * e_IP_B;
        double * e_EA_B;
        double * e_v_B;
        
        double * IP_U_A;
        double * EA_U_A;
        double * IP_U_B;
        double * EA_U_B;
        double * IP_A_Um;
        double * EA_A_Um;
        double * IP_A_Up;
        double * EA_A_Up;

        double * e_IP_AA;
        double * e_EA_AA;
        double * IP_AA_U;
        double * EA_AA_U;
        double * IP_AA_Um;
        double * EA_AA_Um;
        double * IP_AA_Up;
        double * EA_AA_Up;
        
        double * e_IP_BB;
        double * e_EA_BB;
        double * IP_BB_U;
        double * EA_BB_U;
        double * IP_BB_Um;
        double * EA_BB_Um;
        double * IP_BB_Up;
        double * EA_BB_Up;
        
        double * e_IP_AB;
        double * e_EA_AB;
        double * IP_AB_U;
        double * EA_AB_U;
        double * IP_AB_Um;
        double * EA_AB_Um;
        double * IP_AB_Up;
        double * EA_AB_Up;
        
        
        int n_c; 
        int n_a;
        int n_v;
        int aux_n_ao;
        double * H_AV_A;
        double * H_CA_A;
        double * H_CV_A;
        double * H_AV_B;
        double * H_CA_B;
        double * H_CV_B;
        double edshift;
        
        
        
        // RF operator tables
//         double * RF_PV   ;
        double * RF_PH_A   ;
        double * RF_PH_B   ;

        double RF_PS;
        double * RF_PV_JK_A;
        double * RF_PV_JK_B;
        double * RF_PV_AB;
        double * RF_PV_BA;
        
        double * RF_P3_JK_A;
        double * RF_P3_JK_B;
        double * RF_P3_AB;
        double * RF_P3_BA;
        
        
        //functions
        
        //constructor
        UPT_tensors();
        
        //initialization
        int set_par(RI_data * ext_RI_A, RI_data * ext_RI_B,
                    double * ext_eps, double * ext_eps_B,
                    int ext_n_c, int ext_n_a, int ext_n_v, 
                    double * ext_H_AV_A, double * ext_H_CA_A, double * ext_H_CV_A,
                    double * ext_H_AV_B, double * ext_H_CA_B, double * ext_H_CV_B,
                    double ext_edshift);
//         int IPEA(aldet_data * I, int i_set, std::vector<double> avecoe);
        int MPPT(aldet_data * I, int i_set, std::vector<double> avecoe);
        int set_zero();
        int symm();
        
        //1-el properties table generators
//         int calc_IPEA_1_AV(double * P_AV);
//         int calc_IPEA_1_CA(double * P_CA);
//         int calc_IPEA_1_CV(double * P_CV);

//         int calc_EE_1_AV(double * P_AV);
//         int calc_EE_1_CA(double * P_CA);
//         int calc_EE_1_CV(double * P_CV);

        //energy table generators
//         int calc_IPEA_2_CV();        
//         int calc_IPEA_2_AV();        
//         int calc_IPEA_2_CA();
//         int calc_IPEA_2_CCVV();
//         int calc_IPEA_2_CAVV();
//         int calc_IPEA_2_AAVV();
//         int calc_IPEA_2_CCAV();
//         int calc_IPEA_2_CCAA();
        
        int calc_EE_2_CV();        
        int calc_EE_2_AV();        
        int calc_EE_2_CA();
        int calc_EE_2_CCVV();
        int calc_EE_2_CAVV();
        int calc_EE_2_AAVV();
        int calc_EE_2_CCAV();
        int calc_EE_2_CCAA();

//         int calc_U_EE_2_CV();        
//         int calc_U_EE_2_AV();        
//         int calc_U_EE_2_CA();
//         int calc_U_EE_2_CCVV();
//         int calc_U_EE_2_CAVV();
//         int calc_U_EE_2_AAVV();
//         int calc_U_EE_2_CCAV();
//         int calc_U_EE_2_CCAA();
        
//         int calc_2_AV_bc();        
//         int calc_2_CA_bc();
//         int calc_2_CAVV_bc();
//         int calc_2_AAVV_bc();
//         int calc_2_CCAV_bc();
//         int calc_2_CCAA_bc();
        
        
        //main functions
//         int P1_calc_IPEA(double * P1, aldet_data * CI, double * P_AV, double * P_CA, double * P_CV, int AV, int CV, int CA);
//         int P1_calc_EE  (double * P1, aldet_data * CI, double * P_AV, double * P_CA, double * P_CV, int AV, int CV, int CA);
//         int E2_calc_IPEA();
        int E2_calc_EE  ();
//         int E2_calc_U_EE  ();
        
        
//         int E2_calc_bc();
//         int E2_calc_orb_der(int orb);
        
        //destructor
        ~UPT_tensors();
       
       
};

#endif
