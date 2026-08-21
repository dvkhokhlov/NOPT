// localized_dmrg — localized DMRG-CASSCF helpers
// See include/localized_dmrg.h for conventions

#include "localized_dmrg.h"

#include "blas_link.h"   // cblas_dgemm

void build_loc_orbitals(const double* C_act, const double* U, int n_ao, int n_act, double* C_loc) {
    // C_loc = C_act * U : (n_ao x n_act) = (n_ao x n_act)(n_act x n_act)
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_ao, n_act, n_act, 1.0, C_act, n_act,
                U, n_act, 0.0, C_loc, n_act);
}
