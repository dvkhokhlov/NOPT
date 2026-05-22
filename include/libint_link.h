#ifndef LI_LINK
#define LI_LINK
# include <vector>
# include <libint2/shell.h>
# include <libint2/engine.h>

using libint2::Shell;

int AO_1el_from_2shells(double * M, std::vector<Shell> s1, std::vector<Shell> s2, int dim1, int dim2, char E, int r);

int AO_1el_from_shell(double * M, std::vector<Shell> s, int dim, char E, int r);

int DM_to_F_transform(double ** J, double ** K,
                      double ** DM_J, int n_j,
                      double ** DM_K, int n_k, 
                      std::vector<Shell> shells, int n_ao);


double calc_2el_MO_INTS(std::vector<Shell> shells, int n_ao, const double * __restrict__ DM_D, double * J, double * K, double * MO_INTS, const double * __restrict__ R_VEC, const double * __restrict__ L_VEC, int n_v);

double calc_2el_MO_INTS_multi_JK(std::vector<Shell> shells, int n_ao, double ** __restrict__ DM_J, double ** J, int n_j, double ** __restrict__ DM_K, double ** K, int n_k, double * MO_INTS, const double * __restrict__ R_VEC, const double * __restrict__ L_VEC, int n_v, std::vector<libint2::Engine> * engines);

double calc_2el_MO_INTS_4sets(std::vector<Shell> shells, int n_ao, double ** __restrict__ DM_J, double ** J, int n_j, double ** __restrict__ DM_K, double ** K, int n_k, double * MO_INTS, const double * __restrict__ R_VEC1, const double * __restrict__ R_VEC2, const double * __restrict__ L_VEC1, const double * __restrict__ L_VEC2, int n_v, std::vector<libint2::Engine> * engines);

double calc_2el_MO_INTS_4sets_2shells(std::vector<Shell> shells1,std::vector<Shell> shells2, int n_ao1, int n_ao2, double ** __restrict__ DM_J, double ** J, int n_j, double ** __restrict__ DM_K, double ** K, int n_k, double * MO_INTS, const double * __restrict__ R_VEC1, const double * __restrict__ R_VEC2, const double * __restrict__ L_VEC1, const double * __restrict__ L_VEC2, int n_v, std::vector<libint2::Engine> * engines);

double calc_2el_MO_INTS_XYXY(std::vector<Shell>  shells, int n_ao, double * XYXY_INTS, const double * __restrict__ X_VEC, const double * __restrict__ Y_VEC, int n_x, int n_y);

double calc_2el_MO_INTS_XXXY(std::vector<Shell> shells, int n_ao, double * XXXY_INTS, const double * __restrict__ X_VEC, const double * __restrict__ Y_VEC, int n_x, int n_y);

int calc_2el_MP2(std::vector<Shell> shells, int n_ao, int n_mo, double * T, double * V, std::vector<libint2::Engine> * engines);

int RI_V_calc(std::vector<Shell> shells, int n_ao, double * V);

double calc_2el_MO_INTS_XY_RI(std::vector<Shell> shells, int n_ao, std::vector<Shell> aux_shells, int aux_n_ao,double * XY_RI_INTS, const double * __restrict__ X_VEC, const double * __restrict__ Y_VEC, int n_x, int n_y, std::vector<libint2::Engine> * engines);

int calc_2el_AO_INTS_RI(std::vector<Shell> shells, int n_ao, std::vector<Shell> aux_shells, int aux_n_ao, double * AO_RI_INTS);

int calc_2el_AO_INTS_RI_1atom(std::vector<Shell> shells, int n_ao, std::vector<Shell> aux_shells, int aux_n_ao, double * AO_RI_INTS);

int calc_2el_AO_INTS_RI_aux_norm(std::vector<Shell> shells, int n_ao, std::vector<Shell> * aux_shells, int aux_n_ao);

int calc_2el_AO_INTS_RI_aux_transform(std::vector<Shell> shells, int n_ao, std::vector<Shell> aux_shells, int aux_n_ao, double * T, double * AO_RI_INTS);
#endif
