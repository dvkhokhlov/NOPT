#include <stdio.h>
#include <cstring>

inline int set_zero_matr(double *  M, long n){
//     for(int i=0;i<n;i++)M[i]=0.0;
    memset(M,0,n*sizeof(double));
    return 0;
}

inline int set_zero_matr_int(int *  M, long n){
//     for(int i=0;i<n;i++)M[i]=0.0;
    memset(M,0,n*sizeof(int));
    return 0;
}

int PrintMatr(double * M, int n1, int n2, int g);

int fPrintMatr(FILE * stream, double * M, int n1, int n2, int g);

int PrintMatr10(double * M, int n1, int n2, int g);

int fPrintMatr10(FILE * stream, double * M, int n1, int n2, int g);

int PrintMatr_int(int * M, int n1, int n2, int g);

// int PrintCMatr(double _Complex * M, int n1, int n2, int g);

// int copy_with_realloc(double * Out, double * Inp, int dim);


int sort_int_array(int * A, int n);

int transpose(double * M, int n1, int n2);
int transpose_A_to_B(double * Out, double * In, int n1, int n2, int ld1, int ld2);

int transpose_3d_abc_to_bac(double * Out, double * In, int na, int nb, int nc);
int transpose_3d_abc_to_acb(double * Out, double * In, int na, int nb, int nc);
int transpose_3d_abc_to_acb_part_b_sum(double * Out, double * In, int na, int nb0, int nbf, int nb, int nc);
int transpose_3d_abc_to_bca(double * Out, double * In, int na, int nb, int nc);
int transpose_3d_abc_to_cab_sum(double * Out, double * In, int na, int nb, int nc);
int transpose_3d_abc_to_cab_part_b_sum(double * Out, double * In, int na, int nb0, int nbf, int nb, int nc);
int transpose_3d_abc_to_cab_sum_with_b_mult(double * Out, double * In, int na, int nb, int nc, double m);
int transpose_3d_abc_to_abc_part_c_sum(double * Out, double * In, int na, int nb, int nc0, int ncf, int nc);

int transpose_3d_abc_to_bac_int(int * Out, int * In, int na, int nb, int nc);

int lapack_diag(double * mat, double * eval, int n);
int normalize_rotation_rows(double * U, int n);
extern "C"{
    int lapack_herm_diag(double* Mr, double* Mi, double * eval, int n);
}
int lapack_right_eig(double * M, double * eval, int n);

int tr_and_diag(double * H, double * M, double * B, double * V, double * E, int dim);

int S05_calc(double * S, double * M, double * P, int n);

int S05_calc_p(double * S, double * M, double * P, int n, double eps);

int adj_constr_by_svd(double * O, double * L, double * R, double * S, double * B, int n, double sc);

int adj2_constr_by_svd(double * O, const double * L, const double * R, double * S, double * B, int n, double sc);

int matr_constr_by_svd(double * O, double * L, double * R, double * S, double * B, int n);

int zero_svd_ortogonalization(double * L, double * R, double * S, double * B, double * B2, int n);

int lapack_svd(double * mat, double * L, double * R, double * svd, int n);

int magma_svd(double * mat, double * L, double * R, double * svd, int n);

double minor_det(double* MATR, int dim, int a, int b);

double mat_det_calc(double * mat, int N);

double double_minor_det(double* MATR, int dim, int a, int b, int c, int d);

double mat_det_calc_lapack(double * mat, int N);

double adj_calc(double* MATR, int dim, int a, int b);

double mat_det_calc_lapack_no_change(double * mat, int N);

double double_adj_calc(double* MATR, int dim, int a, int b, int c, int d);

int double_adj_mat_constr(double* MATR, int dim, int a, int b, double* OUT);

double adj_matr_constr(double * IN, int N);

int inv_matr_constr(double * IN, int N);

double mat_det_calc_2(double * mat, int N);

double minor_matr_constr(double * matr,int dim);

double qtr_mat_det_calc(double * mat, int dim, int incr);

int triang_minor_matr_constr(double * mat, double * inv_mat,int N);

int ortogonalization(double * matr, int n);

int ortogonalization_no_sq(double * matr, int n, int m);

int ortogonalization_no_sq_row_with_base(double * matr, int n, int m, int ld);

int ortogonalization_no_sq_by_ref(double ** matr, int n, int m);

int ortonormalization_no_sq(double * matr, int n, int m);

double make_ort(double *a, int start_a, double *b, int start_b,int dim);

double make_norm(double *a,int start, int dim);

int matr_cpy(double * Q, int l_Q, double * O, int l_O, int l_col, int l_row);

int matr_cpy_vf(double * Q, int l_Q, double * O, int l_O, int l_col, int l_row, int * VF); //filtration of vertical(v) lines

int matr_cpy_gf(double * Q, int l_Q, double * O, int l_O, int l_col, int l_row, int * GF); //filtration of horizontal(g) lines

int matr_cpy_gvf(double * Q, int l_Q, double * O, int l_O, int l_col, int l_row, int * GF, int * VF); //filtration of all lines

int find_max_abs_el(double * A, int n);

double find_max_abs_value(const double * A, int n, int m, int lda);

int find_N_min_els(int * nums, int N, double * A, int n);

int find_max_coef(int * nums, int N, double * A, int n, int a, int ld);

int find_max_coef_in_row(int * nums, int N, double * A, int n, int a);

int compare_int_ar(int * A, int * B, int n);

int transform(double * O, double * In, double * T, double * buf, int dim);

int transform_back(double * O, double * In, double * T, double * buf, int dim);

int transform_back_nm(double * O, double * In, double * T, double * buf, int dim_i, int dim_o);

int HC_SCE(double * H, double * S, double * C, double * B, double * ev, int n);

int HC_SCE_test(double * H, double * S, double * C, double * B, double * ev, int n);
int HPC_SPCE(double * H, double * S, double * P, double * C, double * B, double * ev, int n);
int HC_SCE_p(double * H, double * S, double * C, double * B, double * ev, int n, double eps);
int HC_SCE_low(double * H, double * S, double * C, double * B, double * ev, int n, double eps);

int HSC_CE(double * H, double * S, double * C, double * B, double * ev, int n);
int HPSPC_PCE(double * H, double * S, double * P, double * C, double * B, double * ev, int n);
int HSC_CE_p(double * H, double * S, double * C, double * B, double * ev, int n, double eps);

int symmetrization(double * A, int n);

int symmetrization_with_scaling(double * A, int n, double c);

int anti_symmetrization(double * A, int n);

int anti_symmetrization_with_scaling(double * A, int n, double c);

double sym_check_for_tensor(double * A, int N, int n, int ld, double f);

int tensor_symmetrization(double * A, double * s, int N, int n, int ld);

int Cholesky(double * L, double * M, int n, double eps);
