// Spin-free EE dressing-tensor engine: the rank-<=1 class builders (CCVV,
// CAVV, CCAV) plus assembly. Formulas transcribed verbatim from S1_tensor_
// engine.md §5.1; integral slabs from the RI_data DF B-tensors (layouts in
// conventions.md §1); bare 1/D kernel (deriv==0). No Hermitization here.

#include "cdas_sf_tensors.h"

#include "blas_link.h"
#include "RI.h"

#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>

extern int num_threads;

#ifndef NDEBUG
namespace {
// True in the EE regime: every e_IP/e_EA entry equals one shared scalar, the
// precondition under which each class denominator must collapse to the
// conventions.md §5 table. Debug guard only; never gates production numerics.
bool sf_ee_flat(const cdas_sf_kernel& K){
    if(K.e_IP.empty()) return false;
    const double c = K.e_IP[0];
    for(double x : K.e_IP) if(std::fabs(x-c) > 1e-12) return false;
    for(double x : K.e_EA) if(std::fabs(x-c) > 1e-12) return false;
    return true;
}
}
#endif

void cdas_sf_tensors::alloc(int na){
    n_a = na;
    const size_t n2 = (size_t)na*na;
    g1.assign(n2, 0.0);
    g2.assign(n2*n2, 0.0);
    g3.assign(n2*n2*n2, 0.0);
    raw_av.clear();
    raw_ca.clear();
    E0 = 0.0;
}

// CCVV (S1 §5.1, tex e0spatial):
//   E0 += Σ_{ij,ab} (ia|jb)[2(ia|jb) − (ib|ja)] / (e_i+e_j−e_a−e_b).
// Per core pair (i,j): G[a][b]=(ia|jb)=Σ_K B^{ia}_K B^{jb}_K via one dgemm of
// the VC B-tensors; (ib|ja)=G[b][a]. D carries no active labels, so the EE
// collapse is the identity — no runtime assertion needed.
int cdas_sf_detail::build_ccvv(const RI_data& R, const double* eps,
                               int n_c, int n_a, int n_v,
                               const double* H_AV, const double* H_CA,
                               const double* H_CV, const cdas_sf_kernel& K,
                               cdas_sf_tensors& out){
    (void)H_AV; (void)H_CA; (void)H_CV; (void)K;
    if(n_c==0 || n_v==0) return 0;
    const long aux = R.aux_n_ao;
    const double* e_c = eps;
    const double* e_v = eps + n_c + n_a;

    std::vector<double> E0_th(num_threads, 0.0);

#pragma omp parallel
    {
        const int nt = omp_get_thread_num();
        std::vector<double> G((size_t)n_v*n_v);
        double e0 = 0.0;
        for(int ij = nt; ij < n_c*n_c; ij += num_threads){
            const int i = ij / n_c;
            const int j = ij % n_c;
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                        n_v, n_v, (int)aux, 1.0,
                        R.VC_RI_M + i*aux, (int)(n_c*aux),
                        R.VC_RI_M + j*aux, (int)(n_c*aux), 0.0,
                        G.data(), n_v);
            for(int a=0; a<n_v; a++)
            for(int b=0; b<n_v; b++){
                const double g = G[(size_t)a*n_v+b];
                const double D = e_c[i]+e_c[j]-e_v[a]-e_v[b];
                e0 += g*(2.0*g - G[(size_t)b*n_v+a]) / D;
            }
        }
        E0_th[nt] = e0;
    }
    for(int th=0; th<num_threads; th++) out.E0 += E0_th[th];
    return 0;
}

// CAVV (S1 §5.1, tex cavvspatial):
//   g1_{tu} += Σ_{i,ab} (ia|tb)[2(ia|ub) − (ib|ua)] / (e_i+e_IP[t]−e_a−e_b).
// D uses the REMOVED active label e_IP[t]; the e_IP[u] twin is supplied by the
// later global Hermitization, so no t<->u symmetrization here. Per core i:
// slab C[a][b][t]=(ia|tb)=VC_i·VAᵀ over aux (one dgemm); (ib|ua)=C[b][a][u].
int cdas_sf_detail::build_cavv(const RI_data& R, const double* eps,
                               int n_c, int n_a, int n_v,
                               const double* H_AV, const double* H_CA,
                               const double* H_CV, const cdas_sf_kernel& K,
                               cdas_sf_tensors& out){
    (void)H_AV; (void)H_CA; (void)H_CV;
    if(n_c==0 || n_v==0 || n_a==0) return 0;
    const long aux = R.aux_n_ao;
    const double* e_c = eps;
    const double* e_v = eps + n_c + n_a;
#ifndef NDEBUG
    const bool ee = sf_ee_flat(K);
    const double epsA = K.e_IP.empty() ? 0.0 : K.e_IP[0];
#endif
    std::vector<std::vector<double>> PH_th(
        num_threads, std::vector<double>((size_t)n_a*n_a, 0.0));

#pragma omp parallel
    {
        const int nt = omp_get_thread_num();
        std::vector<double> C((size_t)n_v*n_v*n_a);
        double* ph = PH_th[nt].data();
        for(int i = nt; i < n_c; i += num_threads){
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                        n_v, n_v*n_a, (int)aux, 1.0,
                        R.VC_RI_M + i*aux, (int)(n_c*aux),
                        R.VA_RI_M, (int)aux, 0.0,
                        C.data(), n_v*n_a);
            for(int a=0; a<n_v; a++)
            for(int b=0; b<n_v; b++){
                const double* Cab = &C[((size_t)a*n_v+b)*n_a]; // (ia|·b)
                const double* Cba = &C[((size_t)b*n_v+a)*n_a]; // (ib|·a)
                for(int t=0; t<n_a; t++){
                    const double D = e_c[i]+K.e_IP[t]-e_v[a]-e_v[b];
#ifndef NDEBUG
                    if(ee) assert(std::fabs(D-(e_c[i]+epsA-e_v[a]-e_v[b])) < 1e-8);
#endif
                    const double pref = Cab[t] / D;
                    for(int u=0; u<n_a; u++)
                        ph[t*n_a+u] += pref*(2.0*Cab[u] - Cba[u]);
                }
            }
        }
    }
    for(int th=0; th<num_threads; th++){
        const double* ph = PH_th[th].data();
        for(int k=0; k<n_a*n_a; k++) out.g1[k] += ph[k];
    }
    return 0;
}

// CCAV (S1 §5.1, tex ccavspatial; overall sign −1; certified KET label u in
// the g1 denominator — NOT e_t):
//   g1_{tu} −= Σ_{ij,a} (it|ja)[2(iu|ja) − (ia|ju)] / (e_i+e_j−e_EA[u]−e_a)
//   E0  += 2 Σ_{ija,t} (it|ja)[2(it|ja) − (ia|jt)] / (e_i+e_j−e_EA[t]−e_a)
// Slabs per (i,j): TA[t][a]=(it|ja)=CA_i·VC_jᵀ; TB[a][u]=(ia|ju)=VC_i·CA_jᵀ.
int cdas_sf_detail::build_ccav(const RI_data& R, const double* eps,
                               int n_c, int n_a, int n_v,
                               const double* H_AV, const double* H_CA,
                               const double* H_CV, const cdas_sf_kernel& K,
                               cdas_sf_tensors& out){
    (void)H_AV; (void)H_CA; (void)H_CV;
    if(n_c==0 || n_v==0 || n_a==0) return 0;
    const long aux = R.aux_n_ao;
    const double* e_c = eps;
    const double* e_v = eps + n_c + n_a;
#ifndef NDEBUG
    const bool ee = sf_ee_flat(K);
    const double epsA = K.e_EA.empty() ? 0.0 : K.e_EA[0];
#endif
    std::vector<std::vector<double>> PH_th(
        num_threads, std::vector<double>((size_t)n_a*n_a, 0.0));
    std::vector<double> E0_th(num_threads, 0.0);

#pragma omp parallel
    {
        const int nt = omp_get_thread_num();
        std::vector<double> TA((size_t)n_a*n_v);
        std::vector<double> TB((size_t)n_v*n_a);
        double* ph = PH_th[nt].data();
        double e0 = 0.0;
        for(int ij = nt; ij < n_c*n_c; ij += num_threads){
            const int i = ij / n_c;
            const int j = ij % n_c;
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                        n_a, n_v, (int)aux, 1.0,
                        R.CA_RI_M + (long)i*n_a*aux, (int)aux,
                        R.VC_RI_M + j*aux, (int)(n_c*aux), 0.0,
                        TA.data(), n_v);
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                        n_v, n_a, (int)aux, 1.0,
                        R.VC_RI_M + i*aux, (int)(n_c*aux),
                        R.CA_RI_M + (long)j*n_a*aux, (int)aux, 0.0,
                        TB.data(), n_a);
            for(int a=0; a<n_v; a++){
                const double* TBa = &TB[(size_t)a*n_a]; // (ia|j·)
                for(int u=0; u<n_a; u++){
                    const double D = e_c[i]+e_c[j]-K.e_EA[u]-e_v[a];
#ifndef NDEBUG
                    if(ee) assert(std::fabs(D-(e_c[i]+e_c[j]-epsA-e_v[a])) < 1e-8);
#endif
                    const double fac = (2.0*TA[(size_t)u*n_v+a] - TBa[u]) / D;
                    for(int t=0; t<n_a; t++)
                        ph[t*n_a+u] -= TA[(size_t)t*n_v+a]*fac;
                }
                for(int t=0; t<n_a; t++){ // E0 constant spill, KET label t
                    const double Dt  = e_c[i]+e_c[j]-K.e_EA[t]-e_v[a];
                    const double TAt = TA[(size_t)t*n_v+a];
                    e0 += 2.0*TAt*(2.0*TAt - TBa[t]) / Dt;
                }
            }
        }
        E0_th[nt] = e0;
    }
    for(int th=0; th<num_threads; th++){
        const double* ph = PH_th[th].data();
        for(int k=0; k<n_a*n_a; k++) out.g1[k] += ph[k];
        out.E0 += E0_th[th];
    }
    return 0;
}

// The remaining classes are not implemented in this increment — fail loudly
// rather than silently accumulate nothing.
int cdas_sf_detail::build_aavv(const RI_data&, const double*, int, int, int,
                               const double*, const double*, const double*,
                               const cdas_sf_kernel&, cdas_sf_tensors&){
    fprintf(stderr, "cdas_sf: build_aavv is not implemented\n");
    abort();
}
int cdas_sf_detail::build_ccaa(const RI_data&, const double*, int, int, int,
                               const double*, const double*, const double*,
                               const cdas_sf_kernel&, cdas_sf_tensors&){
    fprintf(stderr, "cdas_sf: build_ccaa is not implemented\n");
    abort();
}
int cdas_sf_detail::build_cv(const RI_data&, const double*, int, int, int,
                             const double*, const double*, const double*,
                             const cdas_sf_kernel&, cdas_sf_tensors&){
    fprintf(stderr, "cdas_sf: build_cv is not implemented\n");
    abort();
}
int cdas_sf_detail::build_av(const RI_data&, const double*, int, int, int,
                             const double*, const double*, const double*,
                             const cdas_sf_kernel&, cdas_sf_tensors&){
    fprintf(stderr, "cdas_sf: build_av is not implemented\n");
    abort();
}
int cdas_sf_detail::build_ca(const RI_data&, const double*, int, int, int,
                             const double*, const double*, const double*,
                             const cdas_sf_kernel&, cdas_sf_tensors&){
    fprintf(stderr, "cdas_sf: build_ca is not implemented\n");
    abort();
}

int cdas_sf_build(const RI_data& R, const double* eps,
                  int n_c, int n_a, int n_v,
                  const double* H_AV, const double* H_CA, const double* H_CV,
                  const cdas_sf_kernel& K, cdas_sf_tensors& out){
    if((int)K.e_IP.size()!=n_a || (int)K.e_EA.size()!=n_a){
        fprintf(stderr, "cdas_sf_build: kernel e_IP/e_EA must have size n_a=%d\n", n_a);
        abort();
    }
    if(K.deriv!=0){
        fprintf(stderr, "cdas_sf_build: only deriv==0 (bare resolvent) is implemented\n");
        abort();
    }
    out.alloc(n_a);
    // All eight classes, always: an unimplemented builder aborts loudly, so a
    // partial engine can never return a silently incomplete dressing.
    cdas_sf_detail::build_ccvv(R, eps, n_c, n_a, n_v, H_AV, H_CA, H_CV, K, out);
    cdas_sf_detail::build_cavv(R, eps, n_c, n_a, n_v, H_AV, H_CA, H_CV, K, out);
    cdas_sf_detail::build_ccav(R, eps, n_c, n_a, n_v, H_AV, H_CA, H_CV, K, out);
    cdas_sf_detail::build_aavv(R, eps, n_c, n_a, n_v, H_AV, H_CA, H_CV, K, out);
    cdas_sf_detail::build_ccaa(R, eps, n_c, n_a, n_v, H_AV, H_CA, H_CV, K, out);
    cdas_sf_detail::build_cv  (R, eps, n_c, n_a, n_v, H_AV, H_CA, H_CV, K, out);
    cdas_sf_detail::build_av  (R, eps, n_c, n_a, n_v, H_AV, H_CA, H_CV, K, out);
    cdas_sf_detail::build_ca  (R, eps, n_c, n_a, n_v, H_AV, H_CA, H_CV, K, out);
    return 0;
}
