#ifndef CDAS_SF_TENSORS_H
#define CDAS_SF_TENSORS_H

// Public interface for the spin-free EE dressing-tensor engine. Builds
// (E0, g1, g2, g3) for the EE family from DF B-tensors and orbital energies.
// Raw pointers at the NOPT boundary, RAII (std::vector) inside the builders.

#include <vector>

class RI_data;

// Resolvent kernel. e_IP/e_EA are size n_a; under EE every entry equals the
// single active energy eps_A. deriv==0 is the bare 1/D resolvent; deriv!=0
// (derivative kernel) is not implemented and rejected by the builders.
struct cdas_sf_kernel {
    std::vector<double> e_IP, e_EA; // size n_a; EE: all = eps_A
    int deriv = 0;                  // 0: 1/D ; 1: derivative kernel (unimplemented)
};

// Spin-free dressing tensors in the conventions.md §4 storage layouts. raw_av/
// raw_ca hold the bare rank-3 both-double representatives, kept alongside the
// final tensors; empty when the class is absent (n_v==0 kills AV, n_c==0 CA).
struct cdas_sf_tensors {
    int n_a = 0;
    double E0 = 0.0;             // all class constants summed
    std::vector<double> g1;      // n_a^2   g1[t*n_a+u]
    std::vector<double> g2;      // n_a^4
    std::vector<double> g3;      // n_a^6   final (gauge-fixed, Hermitized)
    std::vector<double> raw_av;  // n_a^6
    std::vector<double> raw_ca;  // n_a^6
    void alloc(int n_a);
};

// Builds the aldet RF spin-case arrays from the final g1/g2 plus the raw rank-3
// tables. Defined in a separate translation unit.
int cdas_sf_to_rf(const cdas_sf_tensors& t,
                  double* RF_PS, double* RF_PH,
                  double* RF_PV_JK, double* RF_PV_AB,
                  double* RF_P3_JK, double* RF_P3_AB);

// eps: full orbital-energy array (n_c core | n_a active | n_v virt); ACTIVE
// entries are NOT read (actives come from the kernel). H_AV[t*n_v+a],
// H_CA[i*n_a+t], H_CV[i*n_v+a] are the F^c coupling blocks.
int cdas_sf_build(const RI_data& R, const double* eps,
                  int n_c, int n_a, int n_v,
                  const double* H_AV, const double* H_CA, const double* H_CV,
                  const cdas_sf_kernel& K, cdas_sf_tensors& out);

// Active-basis rotation. Q is n_a×n_a row-major; the new active orbital p is
// phi'_p = Σ_q Q[p*n_a+q] phi_q. Rotates every open active index of g1,g2,g3,
// raw_av,raw_ca by one GEMM per index (2,4,6,6,6 GEMMs; empty raw tables
// skipped). E0 is scalar-invariant. conventions.md §6.
int cdas_sf_rotate(const double* Q, cdas_sf_tensors& t);

// Per-class builders (uniform signature). Each ACCUMULATES into out and never
// zeroes it; individually callable for the per-class gates.
namespace cdas_sf_detail {
    int build_ccvv(const RI_data& R, const double* eps, int n_c, int n_a, int n_v,
                   const double* H_AV, const double* H_CA, const double* H_CV,
                   const cdas_sf_kernel& K, cdas_sf_tensors& out);
    int build_cavv(const RI_data& R, const double* eps, int n_c, int n_a, int n_v,
                   const double* H_AV, const double* H_CA, const double* H_CV,
                   const cdas_sf_kernel& K, cdas_sf_tensors& out);
    int build_ccav(const RI_data& R, const double* eps, int n_c, int n_a, int n_v,
                   const double* H_AV, const double* H_CA, const double* H_CV,
                   const cdas_sf_kernel& K, cdas_sf_tensors& out);
    int build_aavv(const RI_data& R, const double* eps, int n_c, int n_a, int n_v,
                   const double* H_AV, const double* H_CA, const double* H_CV,
                   const cdas_sf_kernel& K, cdas_sf_tensors& out);
    int build_ccaa(const RI_data& R, const double* eps, int n_c, int n_a, int n_v,
                   const double* H_AV, const double* H_CA, const double* H_CV,
                   const cdas_sf_kernel& K, cdas_sf_tensors& out);
    int build_cv(const RI_data& R, const double* eps, int n_c, int n_a, int n_v,
                 const double* H_AV, const double* H_CA, const double* H_CV,
                 const cdas_sf_kernel& K, cdas_sf_tensors& out);
    int build_av(const RI_data& R, const double* eps, int n_c, int n_a, int n_v,
                 const double* H_AV, const double* H_CA, const double* H_CV,
                 const cdas_sf_kernel& K, cdas_sf_tensors& out);
    int build_ca(const RI_data& R, const double* eps, int n_c, int n_a, int n_v,
                 const double* H_AV, const double* H_CA, const double* H_CV,
                 const cdas_sf_kernel& K, cdas_sf_tensors& out);
}

#endif
