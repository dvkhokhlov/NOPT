#ifndef DSRG_SF_BATCH_H
#define DSRG_SF_BATCH_H

// Batched, never-materialized <0|[Vr,T2]|0> energy terms over RI B-tensors.
// B-block layout is aux fastest, [(p*n_q+q)*naux + k], matching RI_data::MO_calc.
// e_c/e_a/e_v are the semicanonical Fock eigenvalues (core/active/virtual);
// s is the DSRG flow parameter. Each call pins BLAS to one thread and parallelises
// over the outer batch index (OMP), restoring the BLAS thread count on return.

// CCVV: E = sum_{i<=j core} f * [scaled (ie|jf)] * [2 (ie|jf) - (if|je)], f=(i==j)?1:2.
// ccvv_zero=true -> bare MP2 (T2 = 1/D, source dressing dropped): CCVV_SOURCE=ZERO.
double dsrg_sf_E_CCVV(const double* VC_RI_M,
                      long n_cor, long n_vac, long naux,
                      const double* e_c, const double* e_v,
                      double s, bool ccvv_zero);

// CAVV: active Hbar1V intermediate (batch over core, virtual inner) contracted
// with the spin-summed 1-RDM L1 (na*na, layout [u*n_act+v], symmetric).
double dsrg_sf_E_CAVV(const double* VC_RI_M, const double* VA_RI_M,
                      long n_cor, long n_act, long n_vac, long naux,
                      const double* e_c, const double* e_a, const double* e_v,
                      const double* L1, double s);

// CCAV: active Hbar1C intermediate (batch over virtual) contracted with
// Eta1 = 2*delta - L1 (na*na, layout [u*n_act+v], symmetric).
double dsrg_sf_E_CCAV(const double* VC_RI_M, const double* CA_RI_M,
                      long n_cor, long n_act, long n_vac, long naux,
                      const double* e_c, const double* e_a, const double* e_v,
                      const double* Eta1, double s);

// CAVV Hbar1 correction: the same virtual-branch intermediate C[v,u] as
// dsrg_sf_E_CAVV, folded in place of the L1 trace. Accumulates the Hermitized
// half-commutator Hbar1[u,v] += 0.5*(C[u,v]+C[v,u]) (na*na, layout [u*n_act+v]).
void dsrg_sf_Hbar1_CAVV(const double* VC_RI_M, const double* VA_RI_M,
                        long n_cor, long n_act, long n_vac, long naux,
                        const double* e_c, const double* e_a, const double* e_v,
                        double s, double* Hbar1);

// CCAV Hbar1 correction: the core-branch intermediate C[v,u] of dsrg_sf_E_CCAV,
// folded with opposite sign: Hbar1[u,v] -= 0.5*(C[u,v]+C[v,u]) (+= into the
// caller buffer, na*na, layout [u*n_act+v]).
void dsrg_sf_Hbar1_CCAV(const double* VC_RI_M, const double* CA_RI_M,
                        long n_cor, long n_act, long n_vac, long naux,
                        const double* e_c, const double* e_a, const double* e_v,
                        double s, double* Hbar1);

#endif
