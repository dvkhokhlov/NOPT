#ifndef RES_FIT_H
#define RES_FIT_H

# include "RI.h"

class PT_tensors_rel;

class PT_tensors_rel 
{
    public:
        
        RI_complex_data * RI; 

        double * e_c;
        double * e_a;
        double * e_IP;
        double * e_EA;
        double * e_v;
        double * IP_U;
        double * EA_U;
        double * IP_Um;
        double * EA_Um;
        double * IP_Up;
        double * EA_Up;

        double * e_IP_2;
        double * e_EA_2;
        double * IP_2_U;
        double * EA_2_U;
        double * IP_2_Um;
        double * EA_2_Um;
        double * IP_2_Up;
        double * EA_2_Up;
        
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
        double * H_AV_r;
        double * H_CA_r;
        double * H_CV_r;
        
        double * H_AV_i;
        double * H_CA_i;
        double * H_CV_i;
        
        double edshift;
        
        
        
        // RF operator tables
        double * RF_PV   ;
        double * RF_PH_r ;
        double * RF_PH_i ;

        double RF_PS;
        double * RF_PV_r;
        double * RF_PV_i;
        
        double * RF_P3_r;
        double * RF_P3_i;
        
        
        //functions
        
        //constructor
        PT_tensors_rel();
        
        //initialization
        int set_par(RI_complex_data * ext_RI, 
                    double * ext_eps_c,
                    int ext_n_c, int ext_n_a, int ext_n_v, 
                    double * ext_H_AV, double * ext_H_CA, double * ext_H_CV,
                    double ext_edshift);
        int IPEA(aldet_data * I, int i_set, std::vector<double> avecoe);
        int set_zero();
        int symm();
        
        //1-el properties table generators
        int calc_IPEA_1_AV(double * P_AV);
        int calc_IPEA_1_CA(double * P_CA);
        int calc_IPEA_1_CV(double * P_CV);

        int calc_EE_1_AV(double * P_AV);
        int calc_EE_1_CA(double * P_CA);
        int calc_EE_1_CV(double * P_CV);

        //energy table generators
        int calc_IPEA_2_CV();        
        int calc_IPEA_2_AV();        
        int calc_IPEA_2_CA();
        int calc_IPEA_2_CCVV();
        int calc_IPEA_2_CAVV();
        int calc_IPEA_2_AAVV();
        int calc_IPEA_2_CCAV();
        int calc_IPEA_2_CCAA();
        
        int calc_EE_2_CV();        
        int calc_EE_2_AV();        
        int calc_EE_2_CA();
        int calc_EE_2_CCVV();
        int calc_EE_2_CAVV();
        int calc_EE_2_AAVV();
        int calc_EE_2_CCAV();
        int calc_EE_2_CCAA();

        int calc_2_AV_bc();        
        int calc_2_CA_bc();
        int calc_2_CAVV_bc();
        int calc_2_AAVV_bc();
        int calc_2_CCAV_bc();
        int calc_2_CCAA_bc();
        
        
        //main functions
        int P1_calc_IPEA(double * P1, aldet_data * CI, double * P_AV, double * P_CA, double * P_CV, int AV, int CV, int CA);
        int P1_calc_EE  (double * P1, aldet_data * CI, double * P_AV, double * P_CA, double * P_CV, int AV, int CV, int CA);
        int E2_calc_IPEA();
        int E2_calc_EE  ();
        
        int E2_calc_bc();
        int E2_calc_orb_der(int orb);
        
        //destructor
        ~PT_tensors_rel();
       
       
};

int sym_JK_T(double *RF_PV_JK, double * RF_PV,int n_a);

int sym_JK_CA_T(double *RF_PV_JK, double * RF_PV,int n_a);

int sym_AB_T(double * RF_PV_AB, double * RF_PV, int n_a);

int sym_AB_CA_T(double * RF_PV_AB, double * RF_PV, int n_a);

#endif
