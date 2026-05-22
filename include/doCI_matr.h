#ifndef NOCI_MATR_H
#define NOCI_MATR_H

# include "inp_par_read.h"

double E_1el_calc(double * H, double * DM, int n1, int n2);

int gen_S_DD(double *D, double *S, int n_1, int d_1, int n_2, int d_2);

int act_to_cor_ort(double * O,double * S,int d_o,int d_n, int v, int n, double * T,double * V, double * L, double * D,int trans);

int act_to_act_ort(double * O,double * S, int n_b, int n_t, int ld, int trans);

int act_to_act_so_matr(double * O, double * S, double * U, int n_b, int n_t, int ld, int trans);

double E_1el_2MO_calc(double * P, double * B, double * K, double * BUF, int n);

int d_o_DM_gen  (double * DM, double coef,
                 double * L, int dim_l,
                 double * R, int dim_r,
                 double * sigma, int n_el);
int d_o_DMV_gen (double * DM,
                 double * L, int dim_l,
                 double * R, int dim_r,
                 double * TH, int n_v , double * B);

int d_o_DMDE_gen(double * DM,
                 double * D, int dim,
                 double * E,
                 double * TH, int n_d, int n_e,
                 int a, int b);


int print_CI_results(FILE * target, int i, int j, double S, double V_nuc, double H_1el, double H_2el, double d_x, double d_y, double d_z);

int gen_S_DD(double *D, double *S, int n_1, int d_1, int n_2, int d_2);

int make_doCI_matr(molecule * A, molecule * B, molecule * EF, int have_efrag, FILE * out_file, int dec);

int make_doCI_fake(molecule * A, int n, double * H, double * S, molecule * EF, double ext_en, int non_diag);

int make_CI_SR_matr(molecule * A, int n, double * H, double * S, molecule * EF, int have_efrag, char * out_folder);

int write_ef(molecule * A, int n, double * H, double * S, molecule * EF, int non_diag);

int Sorb_calc(molecule * A, molecule * Ap);

double gen_doCI_DM(double ** DM, double * S_MO, molecule * A, molecule * B, int n_ao,double ** C, int index1, int index2, int n_s, int dec);

double gen_doCI_DM_1block(double ** DM, double * S_WF, double * S_MO, molecule * A, molecule * B, int dec);
double gen_doCI_DM_1block_PT(double ** DM, double * S_WF, double * S_MO, molecule * A, molecule * B, int dec);

// int gen_S_AO(double * S_MO, molecule * A);

int make_doCI_Smatr(molecule * A, molecule * B, molecule * EF, int have_efrag, FILE * out_file, int dec);

int make_doCI_Vmatr(molecule * A, molecule * B, molecule * EF, int have_efrag, FILE * out_file, int dec);

int make_doCI_H1el(molecule * A, molecule * B, molecule * EF, int have_efrag);

int make_doCI_PT(molecule * A, molecule * B, molecule * EF, int have_efrag, FILE * out_file, int dec, int nshA, int nshB);

#endif
