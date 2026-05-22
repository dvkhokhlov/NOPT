// int n_ext_calc(int *o1, int * o2, int n_v);

// #define GOOGLE_D_M

#define STD_U_M

// #define NO_2EL_DEBUG

// #define DEBUG_NO_EXCH

#include <bitset>
// #include "molecule.h"

#define CI_MAX_SPACE 25

using ci_key = std::bitset<2*CI_MAX_SPACE>;
using ci_key_long = std::bitset<4*CI_MAX_SPACE>;


#ifdef GOOGLE_D_M

#include <google/dense_hash_map>
using ci_map = google::dense_hash_map< ci_key, double> ;
using ci_map_long = google::dense_hash_map< ci_key_long, double> ;
// using ci_map_arr = google::dense_hash_map< ci_key, double *> ;
using ci_map_arr = google::dense_hash_map< ci_key,double *>;
#endif

#ifdef TSL_R_M

#include </home/ilya/src/tsl_robin_map/robin-map/include/tsl/robin_map.h>
using ci_map = tsl::robin_map< ci_key, double> ;
using ci_map_long = tsl::robin_map< ci_key_long, double> ;
// using ci_map_arr = google::dense_hash_map< ci_key, double *> ;
using ci_map_arr = tsl::robin_map< ci_key, double *>;
#endif

#ifdef STD_U_M

#include <unordered_map>
using ci_map = std::unordered_map< ci_key, double> ;
using ci_map_long = std::unordered_map< ci_key_long, double> ;
// using ci_map_arr = google::dense_hash_map< ci_key, double *> ;
using ci_map_arr = std::unordered_map< ci_key,double *>;
#endif


extern double ci_t_cutoff;
extern double ci_v_cutoff;
extern double ci_l_cutoff;

int init_ci_map(ci_map * m);

int init_ci_map_arr(ci_map_arr * m);

// int key_to_oc(long long key, int n);

int no_eq_occ(ci_key k1,ci_key k2,int n, int n2);

int printf_ci_key(ci_key k,int n, int p);

int ci_link_sign(ci_key k1,ci_key k2,int n, int n2);

ci_key det_key(int * a, int n, int z);

long long det_id_1ext(int * a, int x, int n);

double adj_sign(long long key, double * SVD, int n_v);

double det_F_calc(ci_key key,double *E,int n);

int gen_1_el_matr(double * O,ci_map * A,ci_map * B, int n_act, int f);

int gen_1_el_matr_arr(double * O,ci_map_arr * A, int nA, ci_map_arr * B, int nB, int n_act, int f);

int gen_2_el_matr_AA(double * O,ci_map * A,ci_map * B, int n_act, int f);

int gen_2_el_matr_AB(double * O,ci_map * A,ci_map * B, int n_act);

int gen_1_ext_data(ci_map_arr * O, ci_map * A, int n_act, int f);

int gen_2_el_tensor_AA(double * O,ci_map_arr * A,ci_map_arr * B, int n_act, int f);

int gen_2_el_tensor_AB(double * O,ci_map_arr * A,ci_map_arr * B, int n_act);

int ci_skalar_product(double * O, ci_map_arr *I1, int n1, ci_map_arr *I2, int n2);

// ci_map CI_map(molecule * A, int s, int n_act, int f, int i_S);

// ci_map transform_CI(molecule * A, int s, int n_act_i, int n_act_o, int f, int i_S, double * U);

int ci_link(ci_map * O, ci_map ci1, int n1, ci_map ci2, int n2);

int ci_add(ci_map_arr * O, ci_map * In, int i, int i_max);

int ci_arr_clear(ci_map_arr * O);

int ci_array_link(ci_map * O,ci_map ** In,int f,int * s, int * v);

int ci_array_link_to_array(ci_map_arr * O,ci_map ** In,int f, int * S, int * v, int i, int i_max);

int calc_CI_2el_all(double * O, ci_map_arr * I1, ci_map_arr * I2, double * V, int n, int ld, int s1, int s2);

int print_ci_map(ci_map * ci, int n);

int print_ci_map_arr(ci_map_arr * ci, int s, int n);

double ci_den_sum(ci_map_arr * ci, double E, double * eps, int n);

double V_AB_22_J_calc(double * V, double * D_A, double * D_B, int a0, int am, int b0, int bm, int ld);

double V_AB_22_K_calc(double * V, double * D_A, double * D_B, int a0, int am, int b0, int bm, int ld);


ci_map gen_1e_e_a(ci_map_arr * ci,double * e, int n);

ci_map gen_1e_e(ci_map * ci,double * e, int n);

// ci_map gen_2e_e(ci_map_arr * ci,double * e, int n);

ci_map gen_1e_d_a(ci_map_arr * ci,double * e, int n);

ci_map gen_1e_d(ci_map * ci,double * e, int n);

int gen_arr_by_map(ci_map_arr * O, ci_map * ci, int n);
// ci_map gen_2e_d(ci_map_arr * ci,double * e, int n);

// int calc_CI_2el_AB(ci_map * O, ci_map * In, double * V, int n/*, int f*/);

int ci_map_vec_lin_tr(ci_map_arr * ci,double * U, int n);

int ci_vec_lin_tr(double * ci,double * U, int n, int Na, int Nb, int ld);

int ci_arr_norm(ci_map_arr * ci, int n);
