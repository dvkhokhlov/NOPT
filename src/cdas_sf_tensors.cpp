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
#include <algorithm>

extern int num_threads;

namespace {
// True in the EE regime: every e_IP/e_EA entry equals one shared scalar, the
// precondition under which each class denominator collapses to the
// conventions.md §5 table AND under which the rank-3 gauge projection is
// certified. cdas_sf_build gates the aggregate on it; the debug asserts reuse it.
bool sf_ee_flat(const cdas_sf_kernel& K){
    if(K.e_IP.empty()) return false;
    const double c = K.e_IP[0];
    for(double x : K.e_IP) if(std::fabs(x-c) > 1e-12) return false;
    for(double x : K.e_EA) if(std::fabs(x-c) > 1e-12) return false;
    return true;
}
}

void cdas_sf_tensors::alloc(int na){
    n_a = na;
    const size_t n2 = (size_t)na*na;
    g1.assign(n2, 0.0);
    g2.assign(n2*n2, 0.0);
    g3.clear();          // sized in sf_finalize_g3 (keeps the build peak at 3·n_a^6)
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
    if(K.deriv!=0){ fprintf(stderr,"build_ccvv: only deriv==0 is implemented\n"); abort(); }
    if(n_c==0 || n_v==0) return 0;
    const long aux = R.aux_n_ao;
    const double* e_c = eps;
    const double* e_v = eps + n_c + n_a;

    std::vector<double> E0_th(num_threads, 0.0);

    // In-region per-(ij) slab GEMMs: pin the BLAS pool to one thread for the
    // whole region so outer OpenMP owns the parallelism (blas_link.h fence).
#ifdef _OPENBLAS
    int ntb = openblas_get_num_threads(); openblas_set_num_threads(1);
#endif
#ifdef _MKL
    int ntb = mkl_get_max_threads(); mkl_set_num_threads(1);
#endif
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
#ifdef _OPENBLAS
    openblas_set_num_threads(ntb);
#endif
#ifdef _MKL
    mkl_set_num_threads(ntb);
#endif
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
    if(K.deriv!=0){ fprintf(stderr,"build_cavv: only deriv==0 is implemented\n"); abort(); }
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

#ifdef _OPENBLAS
    int ntb = openblas_get_num_threads(); openblas_set_num_threads(1);
#endif
#ifdef _MKL
    int ntb = mkl_get_max_threads(); mkl_set_num_threads(1);
#endif
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
#ifdef _OPENBLAS
    openblas_set_num_threads(ntb);
#endif
#ifdef _MKL
    mkl_set_num_threads(ntb);
#endif
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
    if(K.deriv!=0){ fprintf(stderr,"build_ccav: only deriv==0 is implemented\n"); abort(); }
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

#ifdef _OPENBLAS
    int ntb = openblas_get_num_threads(); openblas_set_num_threads(1);
#endif
#ifdef _MKL
    int ntb = mkl_get_max_threads(); mkl_set_num_threads(1);
#endif
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
#ifdef _OPENBLAS
    openblas_set_num_threads(ntb);
#endif
#ifdef _MKL
    mkl_set_num_threads(ntb);
#endif
    for(int th=0; th<num_threads; th++){
        const double* ph = PH_th[th].data();
        for(int k=0; k<n_a*n_a; k++) out.g1[k] += ph[k];
        out.E0 += E0_th[th];
    }
    return 0;
}

// AAVV (S1 §5.2, tex aavvspatial):
//   g2_{tu,vw} += Σ_{ab} (ta|vb)(ua|wb) / (e_IP[t]+e_IP[v]−e_a−e_b).
// Kernel labels are the CREATION labels t,v. Chunk external virtual a; per a one
// dgemm builds the slab (ta|vb)=Σ_K VA[a,t]·VA[b,v], the kernel scales one
// factor, one dgemm contracts b into the working [(t,v)|(u,w)] tensor.
int cdas_sf_detail::build_aavv(const RI_data& R, const double* eps,
                               int n_c, int n_a, int n_v,
                               const double* H_AV, const double* H_CA,
                               const double* H_CV, const cdas_sf_kernel& K,
                               cdas_sf_tensors& out){
    (void)H_AV; (void)H_CA; (void)H_CV;
    if(K.deriv!=0){ fprintf(stderr,"build_aavv: only deriv==0 is implemented\n"); abort(); }
    if(n_a==0 || n_v==0) return 0;
    const long aux = R.aux_n_ao;
    const double* e_v = eps + n_c + n_a;
    const size_t na2 = (size_t)n_a*n_a;
#ifndef NDEBUG
    const bool ee = sf_ee_flat(K);
    const double epsA = K.e_IP.empty() ? 0.0 : K.e_IP[0];
#endif
    std::vector<std::vector<double>> W_th(
        num_threads, std::vector<double>(na2*na2, 0.0));

#ifdef _OPENBLAS
    int ntb = openblas_get_num_threads(); openblas_set_num_threads(1);
#endif
#ifdef _MKL
    int ntb = mkl_get_max_threads(); mkl_set_num_threads(1);
#endif
#pragma omp parallel
    {
        const int nt = omp_get_thread_num();
        std::vector<double> M((size_t)n_a*n_v*n_a);  // M[t][b·v]=(ta|vb)
        std::vector<double> base(na2*n_v);           // base[(t,v)][b]=(ta|vb)
        std::vector<double> P(na2*n_v);              // base scaled by 1/D
        double* W = W_th[nt].data();
        for(int a = nt; a < n_v; a += num_threads){
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                        n_a, n_v*n_a, (int)aux, 1.0,
                        R.VA_RI_M + (long)a*n_a*aux, (int)aux,
                        R.VA_RI_M, (int)aux, 0.0,
                        M.data(), n_v*n_a);
            for(int t=0; t<n_a; t++)
            for(int b=0; b<n_v; b++)
            for(int v=0; v<n_a; v++)
                base[((size_t)t*n_a+v)*n_v + b] = M[((size_t)t*n_v+b)*n_a + v];
            for(int t=0; t<n_a; t++)
            for(int v=0; v<n_a; v++)
            for(int b=0; b<n_v; b++){
                const double D = K.e_IP[t]+K.e_IP[v]-e_v[a]-e_v[b];
#ifndef NDEBUG
                if(ee) assert(std::fabs(D-(2.0*epsA-e_v[a]-e_v[b])) < 1e-8);
#endif
                P[((size_t)t*n_a+v)*n_v + b] = base[((size_t)t*n_a+v)*n_v + b] / D;
            }
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                        (int)na2, (int)na2, n_v, 1.0,
                        P.data(), n_v,
                        base.data(), n_v, 1.0,
                        W, (int)na2);
        }
    }
#ifdef _OPENBLAS
    openblas_set_num_threads(ntb);
#endif
#ifdef _MKL
    mkl_set_num_threads(ntb);
#endif
    // Repack working->storage (middle labels u,v swap; creators (t,v) grouped in
    // the working tensor, pairs (t,u),(v,w) in storage):
    //   g2[((t·u)v)w] += W[(t·n_a+v)·n_a² + (u·n_a+w)].
    for(int th=0; th<num_threads; th++){
        const double* W = W_th[th].data();
        for(int t=0; t<n_a; t++)
        for(int u=0; u<n_a; u++)
        for(int v=0; v<n_a; v++)
        for(int w=0; w<n_a; w++)
            out.g2[(((size_t)t*n_a+u)*n_a+v)*n_a+w] +=
                W[((size_t)t*n_a+v)*na2 + ((size_t)u*n_a+w)];
    }
    return 0;
}

// CCAA (S1 §5.2, tex ccaaspatial):
//   g2_{tu,vw} += Σ_{ij} (it|jv)(iu|jw) / (e_i+e_j−e_EA[t]−e_EA[v]).
// Creation labels t,v carry e_EA; GEMM shape as AAVV, external core i, integrals
// (it|jv)=Σ_K CA[i,t]·CA[j,v]. Built in a LOCAL buffer, Hermitized, then the two
// certified spills are contracted from it BEFORE it merges into out.g2.
int cdas_sf_detail::build_ccaa(const RI_data& R, const double* eps,
                               int n_c, int n_a, int n_v,
                               const double* H_AV, const double* H_CA,
                               const double* H_CV, const cdas_sf_kernel& K,
                               cdas_sf_tensors& out){
    (void)H_AV; (void)H_CA; (void)H_CV; (void)n_v;
    if(K.deriv!=0){ fprintf(stderr,"build_ccaa: only deriv==0 is implemented\n"); abort(); }
    if(n_c==0 || n_a==0) return 0;
    const long aux = R.aux_n_ao;
    const double* e_c = eps;
    const size_t na2 = (size_t)n_a*n_a;
#ifndef NDEBUG
    const bool ee = sf_ee_flat(K);
    const double epsA = K.e_EA.empty() ? 0.0 : K.e_EA[0];
#endif
    std::vector<std::vector<double>> W_th(
        num_threads, std::vector<double>(na2*na2, 0.0));

#ifdef _OPENBLAS
    int ntb = openblas_get_num_threads(); openblas_set_num_threads(1);
#endif
#ifdef _MKL
    int ntb = mkl_get_max_threads(); mkl_set_num_threads(1);
#endif
#pragma omp parallel
    {
        const int nt = omp_get_thread_num();
        std::vector<double> M((size_t)n_a*n_c*n_a);  // M[t][j·v]=(it|jv)
        std::vector<double> base(na2*n_c);           // base[(t,v)][j]=(it|jv)
        std::vector<double> P(na2*n_c);              // base scaled by 1/D
        double* W = W_th[nt].data();
        for(int i = nt; i < n_c; i += num_threads){
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                        n_a, n_c*n_a, (int)aux, 1.0,
                        R.CA_RI_M + (long)i*n_a*aux, (int)aux,
                        R.CA_RI_M, (int)aux, 0.0,
                        M.data(), n_c*n_a);
            for(int t=0; t<n_a; t++)
            for(int j=0; j<n_c; j++)
            for(int v=0; v<n_a; v++)
                base[((size_t)t*n_a+v)*n_c + j] = M[((size_t)t*n_c+j)*n_a + v];
            for(int t=0; t<n_a; t++)
            for(int v=0; v<n_a; v++)
            for(int j=0; j<n_c; j++){
                const double D = e_c[i]+e_c[j]-K.e_EA[t]-K.e_EA[v];
#ifndef NDEBUG
                if(ee) assert(std::fabs(D-(e_c[i]+e_c[j]-2.0*epsA)) < 1e-8);
#endif
                P[((size_t)t*n_a+v)*n_c + j] = base[((size_t)t*n_a+v)*n_c + j] / D;
            }
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                        (int)na2, (int)na2, n_c, 1.0,
                        P.data(), n_c,
                        base.data(), n_c, 1.0,
                        W, (int)na2);
        }
    }
#ifdef _OPENBLAS
    openblas_set_num_threads(ntb);
#endif
#ifdef _MKL
    mkl_set_num_threads(ntb);
#endif
    // Reduce + repack into a LOCAL class buffer (storage g2[((t·u)v)w]; middle
    // labels u,v swap, as in AAVV).
    std::vector<double> g2C(na2*na2, 0.0);
    for(int th=0; th<num_threads; th++){
        const double* W = W_th[th].data();
        for(int t=0; t<n_a; t++)
        for(int u=0; u<n_a; u++)
        for(int v=0; v<n_a; v++)
        for(int w=0; w<n_a; w++)
            g2C[(((size_t)t*n_a+u)*n_a+v)*n_a+w] +=
                W[((size_t)t*n_a+v)*na2 + ((size_t)u*n_a+w)];
    }
    auto idx = [n_a](int t,int u,int v,int w){
        return (((size_t)t*n_a+u)*n_a+v)*n_a+w; };
    // Hermitize the class tensor: g2C_{tu,vw} = ½(g2C_{tu,vw}+g2C_{ut,wv}).
    {
        std::vector<double> s(g2C);
        for(int t=0; t<n_a; t++)
        for(int u=0; u<n_a; u++)
        for(int v=0; v<n_a; v++)
        for(int w=0; w<n_a; w++)
            g2C[idx(t,u,v,w)] = 0.5*(s[idx(t,u,v,w)] + s[idx(u,t,w,v)]);
    }
    // Certified spills from the Hermitized buffer (S1 §5.2):
    //   g1_{tu} += −2 Σ_w g2C_{tu,ww} + Σ_w g2C_{tw,wu}
    //   E0      += 2 Σ_{tw} g2C_{tt,ww} − Σ_{tw} g2C_{tw,wt}
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++){
        double g = 0.0;
        for(int w=0; w<n_a; w++)
            g += -2.0*g2C[idx(t,u,w,w)] + g2C[idx(t,w,w,u)];
        out.g1[t*n_a+u] += g;
    }
    double e0 = 0.0;
    for(int t=0; t<n_a; t++)
    for(int w=0; w<n_a; w++)
        e0 += 2.0*g2C[idx(t,t,w,w)] - g2C[idx(t,w,w,t)];
    out.E0 += e0;
    for(size_t k=0; k<na2*na2; k++) out.g2[k] += g2C[k];
    return 0;
}

// Both-double resolvent g2 (cert CV.2 and its CA (2e)²-δ mirror — identical
// 8-product form). J,P,invD are n_a²×ne (row = ordered active pair, col =
// external index e); each term is (kernel-scaled slab)·(bare slab) over e.
// Six products pack straight to (t,u)(v,w), two are crossed reads.
// Accumulates into g2[((t·u)v)w] (already g2's storage layout).
static void sf_bothdouble_g2(int n_a, int ne,
                             const std::vector<double>& J,
                             const std::vector<double>& P,
                             const std::vector<double>& invD, double* g2){
    const size_t na2=(size_t)n_a*n_a, NE=(size_t)ne, na4=na2*na2;
    // SJd[(a,b)]=J·invD; M1[(a,b)]=−P[(b,a)]·invD[(a,b)]; M2[(a,b)]=−P·invD[(b,a)].
    std::vector<double> SJd(na2*NE), M1(na2*NE), M2(na2*NE);
#pragma omp parallel for collapse(2) schedule(static)
    for(int p=0;p<n_a;p++) for(int q=0;q<n_a;q++)
        for(size_t e=0;e<NE;e++){
            const size_t pq=((size_t)p*n_a+q)*NE+e, qp=((size_t)q*n_a+p)*NE+e;
            SJd[pq]=J[pq]*invD[pq];
            M1[pq]=-P[qp]*invD[pq];
            M2[pq]=-P[pq]*invD[qp];
        }
    std::vector<double> W(na4);
    // six (t,u)(v,w) products, accumulated into W by successive GEMMs.
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,(int)na2,(int)na2,ne, 2.0,SJd.data(),ne,J.data(),ne,   0.0,W.data(),(int)na2);
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,(int)na2,(int)na2,ne, 1.0,M1.data(), ne,J.data(),ne,   1.0,W.data(),(int)na2);
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,(int)na2,(int)na2,ne,-1.0,SJd.data(),ne,P.data(),ne,   1.0,W.data(),(int)na2);
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,(int)na2,(int)na2,ne, 2.0,J.data(),  ne,SJd.data(),ne, 1.0,W.data(),(int)na2);
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,(int)na2,(int)na2,ne,-1.0,P.data(),  ne,SJd.data(),ne, 1.0,W.data(),(int)na2);
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,(int)na2,(int)na2,ne, 1.0,J.data(),  ne,M1.data(),ne,  1.0,W.data(),(int)na2);
    // two crossed products, then scatter into W (P7→(t,w)(u,v); P8→(w,t)(v,u)).
    std::vector<double> C7(na4), C8(na4);
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,(int)na2,(int)na2,ne, 1.0,P.data(), ne,M2.data(),ne, 0.0,C7.data(),(int)na2);
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,(int)na2,(int)na2,ne, 1.0,M2.data(),ne,P.data(), ne, 0.0,C8.data(),(int)na2);
#pragma omp parallel for collapse(2) schedule(static)
    for(int t=0;t<n_a;t++) for(int u=0;u<n_a;u++)
    for(int v=0;v<n_a;v++) for(int w=0;w<n_a;w++)
        g2[(((size_t)t*n_a+u)*n_a+v)*n_a+w] +=
              W[((size_t)t*n_a+u)*na2 + (size_t)v*n_a+w]
            + C7[((size_t)t*n_a+w)*na2 + (size_t)u*n_a+v]
            + C8[((size_t)w*n_a+t)*na2 + (size_t)v*n_a+u];
}

// CV (S1 §5.3 both-double g2 + F^c-coupled g1 + constant; §5.4 (2e)² δ-image g1).
// Assembled slabs J[(t,u)][(a,i)]=(ia|tu) and P[(x,y)][(a,i)]=(ix|ay) over the
// combined external index e=(a,i). The both-double g2 is realized as §4 GEMMs via
// sf_bothdouble_g2; the rank-1 g1 and the E0 scalar keep the certified scalar form
// (their denominators couple the OUTPUT active label with e — no dense GEMM).
int cdas_sf_detail::build_cv(const RI_data& R, const double* eps,
                             int n_c, int n_a, int n_v,
                             const double* H_AV, const double* H_CA,
                             const double* H_CV, const cdas_sf_kernel& K,
                             cdas_sf_tensors& out){
    (void)H_AV; (void)H_CA;
    if(K.deriv!=0){ fprintf(stderr,"build_cv: only deriv==0 is implemented\n"); abort(); }
    if(n_c==0 || n_v==0 || n_a==0) return 0;
    const long aux = R.aux_n_ao;
    const double* e_c = eps;
    const double* e_v = eps + n_c + n_a;
    const size_t na2=(size_t)n_a*n_a;
    const int ne = n_v*n_c;                 // e = a*n_c + i
#ifndef NDEBUG
    const bool ee = sf_ee_flat(K);
#endif
    // J[(t,u)][(a,i)] = (ia|tu) = AA·VCᵀ.
    std::vector<double> J(na2*(size_t)ne);
    nopt_par_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, (int)na2, ne, (int)aux,
                   1.0, R.AA_RI_M, (int)aux, R.VC_RI_M, (int)aux, 0.0, J.data(), ne);
    // Pbig[(i,x)][(a,y)] = (ix|ay) = CA·VAᵀ, gathered to P[(x,y)][(a,i)].
    std::vector<double> Pbig((size_t)n_c*n_a*n_v*n_a);
    nopt_par_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n_c*n_a, n_v*n_a, (int)aux,
                   1.0, R.CA_RI_M, (int)aux, R.VA_RI_M, (int)aux, 0.0, Pbig.data(), n_v*n_a);
    std::vector<double> P(na2*(size_t)ne);
#pragma omp parallel for collapse(2) schedule(static)
    for(int x=0;x<n_a;x++) for(int y=0;y<n_a;y++)
    for(int a=0;a<n_v;a++) for(int i=0;i<n_c;i++)
        P[((size_t)x*n_a+y)*ne + (size_t)a*n_c+i] =
            Pbig[((size_t)(i*n_a+x)*n_v+a)*n_a+y];
    // invD[(p,q)][(a,i)] = 1/D(p,q); D = e_i+e_IP[p]−e_EA[q]−e_a, p==q → e_i−e_a.
    std::vector<double> invD(na2*(size_t)ne);
#pragma omp parallel for collapse(2) schedule(static)
    for(int p=0;p<n_a;p++) for(int q=0;q<n_a;q++)
    for(int a=0;a<n_v;a++) for(int i=0;i<n_c;i++){
        double D = e_c[i]+K.e_IP[p]-K.e_EA[q]-e_v[a];
        if(p==q) D = e_c[i]-e_v[a];
#ifndef NDEBUG
        if(ee) assert(std::fabs(D-(e_c[i]-e_v[a])) < 1e-8);
#endif
        invD[((size_t)p*n_a+q)*ne + (size_t)a*n_c+i] = 1.0/D;
    }
    sf_bothdouble_g2(n_a, ne, J, P, invD, out.g2.data());

    // rank-1 g1 (F^c-coupled + F^c² constant + (2e)² δ-image) and the E0 scalar.
    std::vector<std::vector<double>> PH_th(num_threads, std::vector<double>(na2,0.0));
    std::vector<double> E0_th(num_threads,0.0);
#pragma omp parallel
    {
        const int nt=omp_get_thread_num();
        double* ph=PH_th[nt].data();
        double e0=0.0;
        for(int e=nt;e<ne;e+=num_threads){
            const int a=e/n_c, i=e%n_c;
            const double Dc=e_c[i]-e_v[a], Fia=H_CV[i*n_v+a];
            auto Je=[&](int p,int q){ return J[((size_t)p*n_a+q)*ne+e]; };
            auto Pe=[&](int p,int q){ return P[((size_t)p*n_a+q)*ne+e]; };
            auto iDe=[&](int p,int q){ return invD[((size_t)p*n_a+q)*ne+e]; };
            for(int t=0;t<n_a;t++) for(int u=0;u<n_a;u++)
                ph[t*n_a+u] += Fia*( (2.0*Je(t,u)-Pe(t,u))/Dc
                                   + (2.0*Je(t,u)-Pe(u,t))*iDe(t,u) );
            e0 += 2.0*Fia*Fia/Dc;
            for(int t=0;t<n_a;t++) for(int u=0;u<n_a;u++) for(int v=0;v<n_a;v++)
                ph[t*n_a+u] += ( 2.0*Je(t,v)*Je(u,v) + 2.0*Pe(v,t)*Pe(v,u)
                               - Je(t,v)*Pe(v,u) - Pe(v,t)*Je(u,v) )*iDe(t,v);
        }
        E0_th[nt]=e0;
    }
    for(int th=0;th<num_threads;th++){
        const double* ph=PH_th[th].data();
        for(size_t k=0;k<na2;k++) out.g1[k]+=ph[k];
        out.E0+=E0_th[th];
    }
    return 0;
}

// AV (S1 §5.3 both-double M table -> raw_av; collapsed g2 + F^c² g1; §5.4 (2e)²
// δ-image g2). External virtual a; slabs Va[a][p][q][r]=(ap|qr), S[p][q][a][r]=
// (pq|ra). GEMM realization: each term's D depends only on ONE factor's own
// labels (D_av symmetric in its +active args), so ONE scaled copy P=Va/D_av
// serves the collapsed g2 and the whole M table. Every D collapses to e_A−e_a.
int cdas_sf_detail::build_av(const RI_data& R, const double* eps,
                             int n_c, int n_a, int n_v,
                             const double* H_AV, const double* H_CA,
                             const double* H_CV, const cdas_sf_kernel& K,
                             cdas_sf_tensors& out){
    (void)H_CA; (void)H_CV;
    if(K.deriv!=0){ fprintf(stderr,"build_av: only deriv==0 is implemented\n"); abort(); }
    if(n_a==0 || n_v==0) return 0;
    const long aux = R.aux_n_ao;
    const double* e_v = eps + n_c + n_a;
    const size_t na2=(size_t)n_a*n_a, na3=na2*n_a, na6=na3*na3;
    if(out.raw_av.empty()) out.raw_av.assign(na6, 0.0);   // accumulate on reuse
#ifndef NDEBUG
    const bool ee = sf_ee_flat(K);
    const double epsA = K.e_IP.empty() ? 0.0 : K.e_IP[0];
#endif
    // D(p,q;r)=e_IP[p]+e_IP[q]−e_EA[r]−e_a; a removed==added coincidence drops one
    // pair (multiplicity: r matching both p and q still drops only one).
    auto Dav=[&](int p,int q,int r,double ea)->double{
        double d = K.e_IP[p]+K.e_IP[q]-K.e_EA[r]-ea;
        if(r==p) d = K.e_IP[q]-ea;
        if(r==q) d = K.e_IP[p]-ea;
#ifndef NDEBUG
        if(ee) assert(std::fabs(d-(epsA-ea)) < 1e-8);
#endif
        return d;
    };
    // Va[((a·n_a+p)·n_a+q)·n_a+r]=(ap|qr); S[((p·n_a+q)·n_v+a)·n_a+r]=(pq|ra).
    std::vector<double> Va((size_t)n_v*na3);
    nopt_par_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                   n_v*n_a, n_a*n_a, (int)aux, 1.0,
                   R.VA_RI_M, (int)aux, R.AA_RI_M, (int)aux, 0.0,
                   Va.data(), n_a*n_a);
    std::vector<double> S(na2*(size_t)n_v*n_a);
    nopt_par_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                   n_a*n_a, n_v*n_a, (int)aux, 1.0,
                   R.AA_RI_M, (int)aux, R.VA_RI_M, (int)aux, 0.0,
                   S.data(), n_v*n_a);
    // P[a][p][q][r] = Va[a][p][q][r] / D_av(p,q;r): the single kernel-scaled slab.
    std::vector<double> P((size_t)n_v*na3);
#pragma omp parallel for schedule(static)
    for(int a=0; a<n_v; a++){
        const double ea=e_v[a];
        for(int p=0;p<n_a;p++) for(int q=0;q<n_a;q++) for(int r=0;r<n_a;r++){
            const size_t k=(size_t)a*na3 + ((size_t)p*n_a+q)*n_a+r;
            P[k] = Va[k]/Dav(p,q,r,ea);
        }
    }

    // --- F^c² g1: g1_{tu} += Σ_a F_ta F_ua / (e_IP[t]−e_a). SFav = F/(e_IP−e). ---
    std::vector<double> SFav((size_t)n_a*n_v);
#pragma omp parallel for schedule(static)
    for(int t=0;t<n_a;t++) for(int a=0;a<n_v;a++)
        SFav[(size_t)t*n_v+a] = H_AV[t*n_v+a]/(K.e_IP[t]-e_v[a]);
    std::vector<double> Wg1(na2);
    nopt_par_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n_a, n_a, n_v, 1.0,
                   SFav.data(), n_v, H_AV, n_v, 0.0, Wg1.data(), n_a);
    for(size_t k=0;k<na2;k++) out.g1[k] += Wg1[k];

    // --- collapsed g2 (4 terms): Cf=SFav·Va, Cp=H_AV·P (both n_a×n_a^3), then
    // g2_{tu,vw} += Cf[t][(u,v,w)] + Cf[v][(w,t,u)] + Cp[u][(t,v,w)] + Cp[w][(v,t,u)]. ---
    {
        std::vector<double> Cf((size_t)n_a*na3), Cp((size_t)n_a*na3);
        nopt_par_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_a, (int)na3, n_v,
                       1.0, SFav.data(), n_v, Va.data(), (int)na3, 0.0, Cf.data(), (int)na3);
        nopt_par_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_a, (int)na3, n_v,
                       1.0, H_AV, n_v, P.data(), (int)na3, 0.0, Cp.data(), (int)na3);
#pragma omp parallel for collapse(2) schedule(static)
        for(int t=0;t<n_a;t++) for(int u=0;u<n_a;u++)
        for(int v=0;v<n_a;v++) for(int w=0;w<n_a;w++)
            out.g2[(((size_t)t*n_a+u)*n_a+v)*n_a+w] +=
                  Cf[(size_t)t*na3 + ((size_t)u*n_a+v)*n_a+w]
                + Cf[(size_t)v*na3 + ((size_t)w*n_a+t)*n_a+u]
                + Cp[(size_t)u*na3 + ((size_t)t*n_a+v)*n_a+w]
                + Cp[(size_t)w*na3 + ((size_t)v*n_a+t)*n_a+u];
    }

    // --- (2e)² δ-image g2 (§5.4): contract the combined (a,m) index. Factor1
    // scaled by 1/Dm(t,v;m,a) [collapse m==t/m==v]; two products, mirror pack. ---
    {
        const int nam = n_v*n_a;                   // e = a*n_a + m
        std::vector<double> A1(na2*(size_t)nam), B1((size_t)nam*na2);
        std::vector<double> A2(na2*(size_t)nam), B2((size_t)nam*na2);
#pragma omp parallel for collapse(2) schedule(static)
        for(int a=0;a<n_v;a++) for(int m=0;m<n_a;m++){
            const double ea=e_v[a];
            for(int t=0;t<n_a;t++) for(int v=0;v<n_a;v++){
                double Dm=K.e_IP[t]+K.e_IP[v]-K.e_EA[m]-ea;
                if(m==t) Dm=K.e_IP[v]-ea;
                if(m==v) Dm=K.e_IP[t]-ea;
#ifndef NDEBUG
                if(ee) assert(std::fabs(Dm-(epsA-ea))<1e-8);
#endif
                const size_t e=(size_t)a*n_a+m;
                A1[((size_t)t*n_a+v)*nam + e] = S[((size_t)(t*n_a+m)*n_v+a)*n_a+v]/Dm; // (tm|va)
                A2[((size_t)t*n_a+v)*nam + e] = Va[(size_t)a*na3+((size_t)t*n_a+v)*n_a+m]/Dm; // (at|vm)
            }
            for(int u=0;u<n_a;u++) for(int w=0;w<n_a;w++){
                const size_t e=(size_t)a*n_a+m;
                B1[e*na2 + (size_t)u*n_a+w] = S[((size_t)(m*n_a+u)*n_v+a)*n_a+w]; // (mu|aw)
                B2[e*na2 + (size_t)w*n_a+u] = S[((size_t)(m*n_a+w)*n_v+a)*n_a+u]; // (mw|ua)
            }
        }
        std::vector<double> C1(na2*na2), C2(na2*na2);
        nopt_par_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, (int)na2, (int)na2, nam,
                       1.0, A1.data(), nam, B1.data(), (int)na2, 0.0, C1.data(), (int)na2);
        nopt_par_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, (int)na2, (int)na2, nam,
                       1.0, A2.data(), nam, B2.data(), (int)na2, 0.0, C2.data(), (int)na2);
#pragma omp parallel for collapse(2) schedule(static)
        for(int t=0;t<n_a;t++) for(int u=0;u<n_a;u++)
        for(int v=0;v<n_a;v++) for(int w=0;w<n_a;w++)
            out.g2[(((size_t)t*n_a+u)*n_a+v)*n_a+w] +=
                  C1[((size_t)t*n_a+v)*na2 + (size_t)u*n_a+w]
                + C2[((size_t)t*n_a+v)*na2 + (size_t)w*n_a+u];
    }

    // --- both-double M table (6 terms, all +1): ONE GEMM C=Pᵀ·Va (n_a^3×n_a^3),
    // then raw_av[t,u,v,w,x,y] += Σ of the 6 certified reads of C[(α)][(β)]. ---
    {
        std::vector<double> C(na6);
        nopt_par_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, (int)na3, (int)na3, n_v,
                       1.0, P.data(), (int)na3, Va.data(), (int)na3, 0.0, C.data(), (int)na3);
        auto Cx=[&](int a1,int a2,int a3,int b1,int b2,int b3)->double{
            return C[(((size_t)a1*n_a+a2)*n_a+a3)*na3 + ((size_t)b1*n_a+b2)*n_a+b3]; };
#pragma omp parallel for collapse(2) schedule(static)
        for(int t=0;t<n_a;t++) for(int u=0;u<n_a;u++)
        for(int v=0;v<n_a;v++) for(int w=0;w<n_a;w++)
        for(int x=0;x<n_a;x++) for(int y=0;y<n_a;y++)
            out.raw_av[(((((size_t)t*n_a+u)*n_a+v)*n_a+w)*n_a+x)*n_a+y] +=
                  Cx(t,v,w, u,x,y) + Cx(t,x,y, u,v,w) + Cx(v,t,u, w,x,y)
                + Cx(v,x,y, w,t,u) + Cx(x,t,u, y,v,w) + Cx(x,v,w, y,t,u);
    }
    return 0;
}

// CA (S1 §5.3 both-double M table -> raw_ca; collapsed g2 with OVERALL −1; §5.4
// (2e)² δ-image g2, double-δ g1 with the mandatory ½, F^c-linear δ g1, F^c²
// g1 with the KET label u, scalar E0). External core i; slab CAAA[i][p][q][r]=
// (ip|qr). M table needs TWO scaled copies (Pq, Pr — +active label in slot q or
// r); Pq also serves the collapsed g2; δ-image via sf_bothdouble_g2 over (i,x).
int cdas_sf_detail::build_ca(const RI_data& R, const double* eps,
                             int n_c, int n_a, int n_v,
                             const double* H_AV, const double* H_CA,
                             const double* H_CV, const cdas_sf_kernel& K,
                             cdas_sf_tensors& out){
    (void)H_AV; (void)H_CV; (void)n_v;
    if(K.deriv!=0){ fprintf(stderr,"build_ca: only deriv==0 is implemented\n"); abort(); }
    if(n_a==0 || n_c==0) return 0;
    const long aux = R.aux_n_ao;
    const double* e_c = eps;
    const size_t na2=(size_t)n_a*n_a, na3=na2*n_a, na6=na3*na3;
    if(out.raw_ca.empty()) out.raw_ca.assign(na6, 0.0);   // accumulate on reuse
#ifndef NDEBUG
    const bool ee = sf_ee_flat(K);
    const double epsA = K.e_EA.empty() ? 0.0 : K.e_EA[0];
#endif
    // CAAA[((i·n_a+p)·n_a+q)·n_a+r]=(ip|qr).
    std::vector<double> CAAA((size_t)n_c*na3);
    nopt_par_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                   n_c*n_a, n_a*n_a, (int)aux, 1.0,
                   R.CA_RI_M, (int)aux, R.AA_RI_M, (int)aux, 0.0,
                   CAAA.data(), n_a*n_a);
    // D(p;q,r)=e_i+e_IP[p]−e_EA[q]−e_EA[r] (p the +active label) with the single-
    // removed collapse; D1(q)=e_i−e_EA[q]. eci passed per call.
    auto Dca=[&](double eci,int p,int q,int r)->double{
        double d=eci+K.e_IP[p]-K.e_EA[q]-K.e_EA[r];
        if(p==q) d=eci-K.e_EA[r];
        if(p==r) d=eci-K.e_EA[q];
#ifndef NDEBUG
        if(ee) assert(std::fabs(d-(eci-epsA))<1e-8);
#endif
        return d;
    };
    auto D1=[&](double eci,int q)->double{
        double d=eci-K.e_EA[q];
#ifndef NDEBUG
        if(ee) assert(std::fabs(d-(eci-epsA))<1e-8);
#endif
        return d;
    };
    // scalar E0 (both-single, KET label t): +2 Σ_{i,t} F^c_{it}²/(e_i−e_EA[t]).
    {
        double e0=0.0;
        for(int i=0;i<n_c;i++) for(int t=0;t<n_a;t++){
            const double F=H_CA[i*n_a+t];
            e0 += 2.0*F*F/D1(e_c[i],t);
        }
        out.E0 += e0;
    }

    // rank-1 g1 (F^c², F^c-linear δ, double-δ): kept scalar — the denominators
    // couple the OUTPUT label with the contracted ones. (t,u) striped, i inner.
#pragma omp parallel
    {
        const int nt=omp_get_thread_num();
        for(int tu=nt; tu<n_a*n_a; tu+=num_threads){
            const int t=tu/n_a, u=tu%n_a;
            for(int i=0;i<n_c;i++){
                const double eci=e_c[i];
                const double* Ci=CAAA.data()+(size_t)i*na3;
                const double Fit=H_CA[i*n_a+t], Fiu=H_CA[i*n_a+u];
                out.g1[t*n_a+u] += -Fit*Fiu/D1(eci,u);      // F^c² g1 (KET label u)
                double g1acc=0.0;
                for(int v=0;v<n_a;v++){
                    const double ivtu=Ci[(v*n_a+t)*n_a+u];  // (iv|tu)
                    g1acc += H_CA[i*n_a+v]*( (2.0*ivtu-Ci[(t*n_a+v)*n_a+u])/D1(eci,v)
                                           + (2.0*ivtu-Ci[(u*n_a+v)*n_a+t])/Dca(eci,t,u,v) );
                    for(int w=0;w<n_a;w++){                  // double-δ g1 (½ mandatory)
                        const double ivtw=Ci[(v*n_a+t)*n_a+w], ivuw=Ci[(v*n_a+u)*n_a+w];
                        const double iwtv=Ci[(w*n_a+t)*n_a+v], iwuv=Ci[(w*n_a+u)*n_a+v];
                        g1acc += 0.5*(2.0*ivtw*ivuw + 2.0*iwtv*iwuv
                                      - ivtw*iwuv - iwtv*ivuw)/Dca(eci,t,v,w);
                    }
                }
                out.g1[t*n_a+u] += g1acc;
            }
        }
    }

    // Pq[i][(p,q,r)] = CAAA / D(q; p,r); Pr[i][(p,q,r)] = CAAA / D(r; p,q).
    std::vector<double> Pq((size_t)n_c*na3), Pr((size_t)n_c*na3);
#pragma omp parallel for schedule(static)
    for(int i=0;i<n_c;i++){
        const double eci=e_c[i];
        for(int p=0;p<n_a;p++) for(int q=0;q<n_a;q++) for(int r=0;r<n_a;r++){
            const size_t k=(size_t)i*na3 + ((size_t)p*n_a+q)*n_a+r;
            Pq[k]=CAAA[k]/Dca(eci,q,p,r);
            Pr[k]=CAAA[k]/Dca(eci,r,p,q);
        }
    }

    // collapsed g2 (4 terms, OVERALL −1): C1=SFcaᵀ·CAAA, C3=H_CAᵀ·Pq (n_a×n_a^3);
    // g2 −= C1[u][(t,v,w)] + C1[w][(v,t,u)] + C3[t][(u,v,w)] + C3[v][(w,t,u)].
    {
        std::vector<double> SFca((size_t)n_c*n_a);
#pragma omp parallel for schedule(static)
        for(int i=0;i<n_c;i++) for(int q=0;q<n_a;q++)
            SFca[(size_t)i*n_a+q] = H_CA[i*n_a+q]/(e_c[i]-K.e_EA[q]);
        std::vector<double> C1((size_t)n_a*na3), C3((size_t)n_a*na3);
        nopt_par_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n_a, (int)na3, n_c,
                       1.0, SFca.data(), n_a, CAAA.data(), (int)na3, 0.0, C1.data(), (int)na3);
        nopt_par_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n_a, (int)na3, n_c,
                       1.0, H_CA, n_a, Pq.data(), (int)na3, 0.0, C3.data(), (int)na3);
#pragma omp parallel for collapse(2) schedule(static)
        for(int t=0;t<n_a;t++) for(int u=0;u<n_a;u++)
        for(int v=0;v<n_a;v++) for(int w=0;w<n_a;w++)
            out.g2[(((size_t)t*n_a+u)*n_a+v)*n_a+w] -=
                  C1[(size_t)u*na3 + ((size_t)t*n_a+v)*n_a+w]
                + C1[(size_t)w*na3 + ((size_t)v*n_a+t)*n_a+u]
                + C3[(size_t)t*na3 + ((size_t)u*n_a+v)*n_a+w]
                + C3[(size_t)v*na3 + ((size_t)w*n_a+t)*n_a+u];
    }

    // (2e)² δ-image g2: CA (2e)² mirror of CV.2 over e=(i,x). Jix[(t,u)]=(ix|tu),
    // Qix[(p,q)]=(ip|xq), invDx[(p,q)]=1/D(p,q;x). Reuses sf_bothdouble_g2 (+).
    {
        const int ne=n_c*n_a;                   // e = i*n_a + x
        std::vector<double> Jix(na2*(size_t)ne), Qix(na2*(size_t)ne), invDx(na2*(size_t)ne);
#pragma omp parallel for collapse(2) schedule(static)
        for(int i=0;i<n_c;i++) for(int x=0;x<n_a;x++){
            const double eci=e_c[i];
            const size_t e=(size_t)i*n_a+x;
            for(int p=0;p<n_a;p++) for(int q=0;q<n_a;q++){
                Jix[((size_t)p*n_a+q)*ne+e] = CAAA[(size_t)i*na3+((size_t)x*n_a+p)*n_a+q];
                Qix[((size_t)p*n_a+q)*ne+e] = CAAA[(size_t)i*na3+((size_t)p*n_a+x)*n_a+q];
                invDx[((size_t)p*n_a+q)*ne+e] = 1.0/Dca(eci,p,q,x);
            }
        }
        sf_bothdouble_g2(n_a, ne, Jix, Qix, invDx, out.g2.data());
    }

    // both-double M table (6 terms, all +1): Cq=Pqᵀ·CAAA and Cr=Prᵀ·CAAA share ONE
    // n_a^6 scratch; raw_ca gathers the 4 Cq-reads then the 2 Cr-reads (cert §4).
    {
        std::vector<double> C(na6);
        nopt_par_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, (int)na3, (int)na3, n_c,
                       1.0, Pq.data(), (int)na3, CAAA.data(), (int)na3, 0.0, C.data(), (int)na3);
        auto Cx=[&](int a1,int a2,int a3,int b1,int b2,int b3)->double{
            return C[(((size_t)a1*n_a+a2)*n_a+a3)*na3 + ((size_t)b1*n_a+b2)*n_a+b3]; };
#pragma omp parallel for collapse(2) schedule(static)
        for(int t=0;t<n_a;t++) for(int u=0;u<n_a;u++)
        for(int v=0;v<n_a;v++) for(int w=0;w<n_a;w++)
        for(int x=0;x<n_a;x++) for(int y=0;y<n_a;y++)
            out.raw_ca[(((((size_t)t*n_a+u)*n_a+v)*n_a+w)*n_a+x)*n_a+y] +=
                  Cx(w,x,y, t,u,v) + Cx(u,t,w, v,x,y) + Cx(u,x,y, v,t,w) + Cx(y,t,w, x,u,v);
        nopt_par_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, (int)na3, (int)na3, n_c,
                       1.0, Pr.data(), (int)na3, CAAA.data(), (int)na3, 0.0, C.data(), (int)na3);
#pragma omp parallel for collapse(2) schedule(static)
        for(int t=0;t<n_a;t++) for(int u=0;u<n_a;u++)
        for(int v=0;v<n_a;v++) for(int w=0;w<n_a;w++)
        for(int x=0;x<n_a;x++) for(int y=0;y<n_a;y++)
            out.raw_ca[(((((size_t)t*n_a+u)*n_a+v)*n_a+w)*n_a+x)*n_a+y] +=
                  Cx(w,u,v, t,x,y) + Cx(y,u,v, x,t,w);
    }
    return 0;
}

// §6 finalization, ranks <= 2 (S1 §6 steps 1-2 / conventions.md §3). g1
// Hermitized ½(g1+g1ᵀ); g2 Hermitized ½(g2_{tu,vw}+g2_{ut,wv}) then pair-
// symmetrized ½(g2_{tu,vw}+g2_{vw,tu}). g3 stays zero and raw_av/raw_ca stay
// raw — the §6.3 g3 gauge projection is a later increment.
static void sf_finalize_rank12(cdas_sf_tensors& out){
    const int n = out.n_a;
    auto id2 = [n](int t,int u,int v,int w){ return (((size_t)t*n+u)*n+v)*n+w; };
    {   std::vector<double> s(out.g1);
        for(int t=0;t<n;t++) for(int u=0;u<n;u++)
            out.g1[(size_t)t*n+u] = 0.5*(s[(size_t)t*n+u] + s[(size_t)u*n+t]); }
    {   std::vector<double> s(out.g2);              // Hermitize
        for(int t=0;t<n;t++) for(int u=0;u<n;u++)
        for(int v=0;v<n;v++) for(int w=0;w<n;w++)
            out.g2[id2(t,u,v,w)] = 0.5*(s[id2(t,u,v,w)] + s[id2(u,t,w,v)]); }
    {   std::vector<double> s(out.g2);              // pair-symmetrize
        for(int t=0;t<n;t++) for(int u=0;u<n;u++)
        for(int v=0;v<n;v++) for(int w=0;w<n;w++)
            out.g2[id2(t,u,v,w)] = 0.5*(s[id2(t,u,v,w)] + s[id2(v,w,t,u)]); }
}

#ifndef NDEBUG
// §6.4 assertions (ranks <= 2): the finalized g1/g2 are Hermitian and g2 is
// pair-symmetric. The per-class EE-collapse asserts are live in the builders;
// the CCAA spill is closed against its own Hermitized class tensor in build_ccaa.
static void sf_assert_rank12(const cdas_sf_tensors& out){
    const int n = out.n_a;
    auto id2 = [n](int t,int u,int v,int w){ return (((size_t)t*n+u)*n+v)*n+w; };
    double sc = 1.0;
    for(double x : out.g2) sc = std::max(sc, std::fabs(x));
    for(int t=0;t<n;t++) for(int u=0;u<n;u++){
        assert(std::fabs(out.g1[(size_t)t*n+u]-out.g1[(size_t)u*n+t]) <= 1e-10*sc);
        for(int v=0;v<n;v++) for(int w=0;w<n;w++){
            assert(std::fabs(out.g2[id2(t,u,v,w)]-out.g2[id2(u,t,w,v)]) <= 1e-10*sc);
            assert(std::fabs(out.g2[id2(t,u,v,w)]-out.g2[id2(v,w,t,u)]) <= 1e-10*sc);
        }
    }
}
#endif

// §6.3 g3 gauge projection (S1 §6.3 / cert_a2_g3_gauge.md §9.1). ACCUMULATES the
// UN-Hermitized projection sym(D)−½asymC(D) of ONE raw both-double table into
// G_accum; the Hermitization lives in sf_finalize_g3. D(·)=T(·)−T(·, u↔w) is
// read straight from T inline (no full-size D/G temporary). One code path for
// AV and CA. T read-only.
namespace cdas_sf_detail {
void sf_gauge_project_g3(const std::vector<double>& T, int n_a,
                         std::vector<double>& G_accum){
    const int n = n_a;
    auto ix = [n](int t,int u,int v,int w,int x,int y)->size_t {
        return (((((size_t)t*n+u)*n+v)*n+w)*n+x)*n+y; };
    // D on a pair-tuple: swap the two annihilator slots (positions 2 and 4).
    auto Dv = [&](int a,int b,int c,int d,int e,int f)->double {
        return T[ix(a,b,c,d,e,f)] - T[ix(a,d,c,b,e,f)]; };
    // sym: 6 whole-pair orderings; asymC: 6 signed creator-slot permutations,
    // annihilators fixed (id/3-cycles +, swaps −). Accumulate s−½c into g3.
#pragma omp parallel for collapse(2) schedule(static)
    for(int t=0;t<n;t++) for(int u=0;u<n;u++)
    for(int v=0;v<n;v++) for(int w=0;w<n;w++)
    for(int x=0;x<n;x++) for(int y=0;y<n;y++){
        const double s = ( Dv(t,u,v,w,x,y) + Dv(t,u,x,y,v,w)
                         + Dv(v,w,t,u,x,y) + Dv(v,w,x,y,t,u)
                         + Dv(x,y,t,u,v,w) + Dv(x,y,v,w,t,u) ) / 6.0;
        const double c = ( Dv(t,u,v,w,x,y) - Dv(v,u,t,w,x,y)
                         - Dv(x,u,v,w,t,y) - Dv(t,u,x,w,v,y)
                         + Dv(v,u,x,w,t,y) + Dv(x,u,t,w,v,y) ) / 6.0;
        G_accum[ix(t,u,v,w,x,y)] += s - 0.5*c;
    }
}
}

// Fill g3 from the raw both-double tables: allocate g3 here (kept out of alloc()
// so the build peaks at 3·n_a^6), accumulate the UN-Hermitized projection of each
// class, then Hermitize ONCE in place. Hermitization is linear, so ½(G+G†) of the
// AV+CA sum equals the per-class Hermitizations summed. raw_av/raw_ca read only.
static void sf_finalize_g3(cdas_sf_tensors& out){
    const int n = out.n_a;
    const size_t na6 = (size_t)n*n*n*n*n*n;
    out.g3.assign(na6, 0.0);
    if(!out.raw_av.empty())
        cdas_sf_detail::sf_gauge_project_g3(out.raw_av, n, out.g3);
    if(!out.raw_ca.empty())
        cdas_sf_detail::sf_gauge_project_g3(out.raw_ca, n, out.g3);
    // In-place Hermitization: pair (t,u,v,w,x,y) with its dagger (u,t,w,v,y,x);
    // average each unordered pair once (self-conjugate elements untouched).
    auto ix = [n](int t,int u,int v,int w,int x,int y)->size_t {
        return (((((size_t)t*n+u)*n+v)*n+w)*n+x)*n+y; };
    double* g3 = out.g3.data();
#pragma omp parallel for collapse(2) schedule(static)
    for(int t=0;t<n;t++) for(int u=0;u<n;u++)
    for(int v=0;v<n;v++) for(int w=0;w<n;w++)
    for(int x=0;x<n;x++) for(int y=0;y<n;y++){
        const size_t i = ix(t,u,v,w,x,y), id = ix(u,t,w,v,y,x);
        if(i < id){ double s = 0.5*(g3[i]+g3[id]); g3[i] = g3[id] = s; }
    }
}

#ifndef NDEBUG
// §6.4 assertion (rank 3): finalized g3 is pair-S3 symmetric (invariant under the
// two generating pair swaps) and Hermitian under g3_{tu,vw,xy}=g3_{ut,wv,yx}.
static void sf_assert_g3(const cdas_sf_tensors& out){
    const int n = out.n_a;
    auto ix = [n](int t,int u,int v,int w,int x,int y)->size_t {
        return (((((size_t)t*n+u)*n+v)*n+w)*n+x)*n+y; };
    double sc = 1.0;
    for(double z : out.g3) sc = std::max(sc, std::fabs(z));
    const double tol = 1e-10*sc;
    for(int t=0;t<n;t++) for(int u=0;u<n;u++)
    for(int v=0;v<n;v++) for(int w=0;w<n;w++)
    for(int x=0;x<n;x++) for(int y=0;y<n;y++){
        const double g = out.g3[ix(t,u,v,w,x,y)];
        assert(std::fabs(g - out.g3[ix(v,w,t,u,x,y)]) <= tol); // pair0<->pair1
        assert(std::fabs(g - out.g3[ix(t,u,x,y,v,w)]) <= tol); // pair1<->pair2
        assert(std::fabs(g - out.g3[ix(u,t,w,v,y,x)]) <= tol); // Hermitian
    }
}
#endif

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
    if(!sf_ee_flat(K)){
        fprintf(stderr, "cdas_sf_build: the aggregate's rank-3 gauge projection is "
                        "certified for a single EE scalar only; non-EE callers must use "
                        "the cdas_sf_detail builders (rank<=2 only)\n");
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
    sf_finalize_rank12(out);
    sf_finalize_g3(out);
#ifndef NDEBUG
    sf_assert_rank12(out);
    sf_assert_g3(out);
#endif
    return 0;
}

// Rotate one open active index of a rank-r tensor by Q: view it as
// (n_a)×(n_a^{r-1}) row-major; ONE GEMM computes out[j,p]=Σ_q Q[p,q] in[q,j],
// rotating the leading index AND cyclically shifting it to the back, so r passes
// rotate every index once and restore the original axis order. cur/scratch swap.
static void sf_rotate_indices(const double* Q, int n_a, int rank,
                              std::vector<double>& T){
    if(T.empty()) return;
    size_t total = 1; for(int k=0; k<rank; k++) total *= (size_t)n_a;
    const int R = (int)(total / (size_t)n_a);   // n_a^{rank-1}
    std::vector<double> scratch(total);
    double* cur = T.data();
    double* nxt = scratch.data();
    for(int pass=0; pass<rank; ++pass){
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans,
                    R, n_a, n_a, 1.0,
                    cur, R, Q, n_a, 0.0, nxt, n_a);
        std::swap(cur, nxt);
    }
    if(cur != T.data()) std::copy(cur, cur+total, T.data());
}

// Rotate every open active index of each tensor by Q (2,4,6,6,6 GEMMs; empty raw
// tables skipped). E0 is a rotation-invariant scalar contraction, left untouched.
// Rotation commutes with the §6 Hermitization/gauge (uniform per-index Q), so a
// rotated build and a build-then-rotate agree elementwise (G1c). conventions.md §6.
int cdas_sf_rotate(const double* Q, cdas_sf_tensors& t){
    const int n = t.n_a;
    if(n <= 0) return 0;
    sf_rotate_indices(Q, n, 2, t.g1);
    sf_rotate_indices(Q, n, 4, t.g2);
    sf_rotate_indices(Q, n, 6, t.g3);
    if(!t.raw_av.empty()) sf_rotate_indices(Q, n, 6, t.raw_av);
    if(!t.raw_ca.empty()) sf_rotate_indices(Q, n, 6, t.raw_ca);
    return 0;
}
