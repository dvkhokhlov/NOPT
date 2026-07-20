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

// CV (S1 §5.3 both-double g2 + F^c-coupled g1 + constant; §5.4 (2e)² δ-image
// g1). Per core i: slabs J[a][t][u]=(ia|tu) (VC_i·AAᵀ) and P[x][a][y]=(ix|ay)
// (CA_i·VAᵀ). Each g2 assignment carries its OWN denominator D(p,q)=e_i+e_IP[p]
// −e_EA[q]−e_a with the p==q collapse → e_i−e_a; g2 is written straight to the
// storage layout (no BLAS outer product, so no repack). External i is striped.
int cdas_sf_detail::build_cv(const RI_data& R, const double* eps,
                             int n_c, int n_a, int n_v,
                             const double* H_AV, const double* H_CA,
                             const double* H_CV, const cdas_sf_kernel& K,
                             cdas_sf_tensors& out){
    (void)H_AV; (void)H_CA;
    if(n_c==0 || n_v==0 || n_a==0) return 0;
    const long aux = R.aux_n_ao;
    const double* e_c = eps;
    const double* e_v = eps + n_c + n_a;
    const size_t na2 = (size_t)n_a*n_a;
#ifndef NDEBUG
    const bool ee = sf_ee_flat(K);
#endif
    std::vector<std::vector<double>> G2_th(
        num_threads, std::vector<double>(na2*na2, 0.0));
    std::vector<std::vector<double>> PH_th(
        num_threads, std::vector<double>(na2, 0.0));
    std::vector<double> E0_th(num_threads, 0.0);

#pragma omp parallel
    {
        const int nt = omp_get_thread_num();
        std::vector<double> J((size_t)n_v*na2);      // J[(a·n_a+t)n_a+u]=(ia|tu)
        std::vector<double> P((size_t)n_a*n_v*n_a);  // P[(x·n_v+a)n_a+y]=(ix|ay)
        std::vector<double> Pa(na2);                 // Pa[x·n_a+y]=(ix|ay), fixed a
        std::vector<double> invD(na2);               // 1/D(p,q), fixed (i,a)
        double* g2 = G2_th[nt].data();
        double* ph = PH_th[nt].data();
        double e0 = 0.0;
        for(int i = nt; i < n_c; i += num_threads){
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                        n_v, (int)na2, (int)aux, 1.0,
                        R.VC_RI_M + i*aux, (int)(n_c*aux),
                        R.AA_RI_M, (int)aux, 0.0,
                        J.data(), (int)na2);
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                        n_a, n_v*n_a, (int)aux, 1.0,
                        R.CA_RI_M + (long)i*n_a*aux, (int)aux,
                        R.VA_RI_M, (int)aux, 0.0,
                        P.data(), n_v*n_a);
            const double eci = e_c[i];
            for(int a=0; a<n_v; a++){
                const double Dc = eci - e_v[a];          // EE-collapsed denom
                const double* Ja = &J[(size_t)a*na2];    // Ja[t·n_a+u]=(ia|tu)
                for(int x=0; x<n_a; x++)
                for(int y=0; y<n_a; y++)
                    Pa[x*n_a+y] = P[((size_t)x*n_v+a)*n_a+y];
                for(int p=0; p<n_a; p++)
                for(int q=0; q<n_a; q++){
                    double D = eci + K.e_IP[p] - e_v[a] - K.e_EA[q];
                    if(p==q) D = Dc;                     // removed==added collapse
#ifndef NDEBUG
                    if(ee) assert(std::fabs(D-Dc) < 1e-8);
#endif
                    invD[p*n_a+q] = 1.0/D;
                }
                // g2, 7 assignments — each keeps its OWN D(p,q); the (ia|tu)(ia|vw)
                // diagonal carries factor 2 split over D(t,u) and D(v,w).
                for(int t=0; t<n_a; t++)
                for(int u=0; u<n_a; u++){
                    const double Jtu = Ja[t*n_a+u];
                    const double id_tu = invD[t*n_a+u];
                    for(int v=0; v<n_a; v++)
                    for(int w=0; w<n_a; w++){
                        const double Jvw = Ja[v*n_a+w];
                        const double g =
                            (2.0*Jtu*Jvw - Pa[u*n_a+t]*Jvw - Jtu*Pa[v*n_a+w])*id_tu
                          + (2.0*Jtu*Jvw - Pa[t*n_a+u]*Jvw - Jtu*Pa[w*n_a+v])*invD[v*n_a+w]
                          - Pa[t*n_a+w]*Pa[u*n_a+v]*invD[v*n_a+u]
                          - Pa[w*n_a+t]*Pa[v*n_a+u]*invD[t*n_a+w];
                        g2[((size_t)(t*n_a+u)*n_a+v)*n_a+w] += g;
                    }
                }
                // g1 F^c-coupled (§5.3) + the F^c² constant.
                const double Fia = H_CV[i*n_v+a];
                for(int t=0; t<n_a; t++)
                for(int u=0; u<n_a; u++){
                    const double Jtu = Ja[t*n_a+u];
                    ph[t*n_a+u] += Fia*( (2.0*Jtu - Pa[t*n_a+u])/Dc
                                       + (2.0*Jtu - Pa[u*n_a+t])*invD[t*n_a+u] );
                }
                e0 += 2.0*Fia*Fia/Dc;
                // g1 (2e)² δ-image (§5.4); denominator D(t,v)=invD[t·n_a+v].
                for(int t=0; t<n_a; t++)
                for(int u=0; u<n_a; u++)
                for(int v=0; v<n_a; v++)
                    ph[t*n_a+u] += ( 2.0*Ja[t*n_a+v]*Ja[u*n_a+v]
                                   + 2.0*Pa[v*n_a+t]*Pa[v*n_a+u]
                                   - Ja[t*n_a+v]*Pa[v*n_a+u]
                                   - Pa[v*n_a+t]*Ja[u*n_a+v] )*invD[t*n_a+v];
            }
        }
        E0_th[nt] = e0;
    }
    for(int th=0; th<num_threads; th++){
        const double* g2 = G2_th[th].data();
        const double* ph = PH_th[th].data();
        for(size_t k=0; k<na2*na2; k++) out.g2[k] += g2[k];
        for(size_t k=0; k<na2;     k++) out.g1[k] += ph[k];
        out.E0 += E0_th[th];
    }
    return 0;
}

// AV (S1 §5.3 both-double M table -> raw_av; collapsed g2/g1; §5.4 (2e)² δ-image
// g2). External virtual a. Shared slabs Va[a][p][q][r]=(ap|qr) (VA·AA over aux)
// and S[p][q][a][r]=(pq|ra) (AA·VA over aux). The OUTPUT pair (t,u) is striped
// over threads — each owns disjoint raw_av/g2/g1 rows (as the stock 3-body
// builder stripes its active pair), so no n_a^6 per-thread copy is needed. The
// M table is a straight loop nest with hoisted slab pointers: literal
// transcription over GEMM-shaping at this stage. Every D collapses to e_A−e_a.
int cdas_sf_detail::build_av(const RI_data& R, const double* eps,
                             int n_c, int n_a, int n_v,
                             const double* H_AV, const double* H_CA,
                             const double* H_CV, const cdas_sf_kernel& K,
                             cdas_sf_tensors& out){
    (void)H_CA; (void)H_CV;
    if(n_a==0 || n_v==0) return 0;
    const long aux = R.aux_n_ao;
    const double* e_v = eps + n_c + n_a;
    const size_t na2=(size_t)n_a*n_a, na3=na2*n_a, na6=na3*na3;
    out.raw_av.assign(na6, 0.0);
#ifndef NDEBUG
    const bool ee = sf_ee_flat(K);
    const double epsA = K.e_IP.empty() ? 0.0 : K.e_IP[0];
#endif
    // Va[((a·n_a+p)·n_a+q)·n_a+r]=(ap|qr); S[((p·n_a+q)·n_v+a)·n_a+r]=(pq|ra).
    std::vector<double> Va((size_t)n_v*na3);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                n_v*n_a, n_a*n_a, (int)aux, 1.0,
                R.VA_RI_M, (int)aux, R.AA_RI_M, (int)aux, 0.0,
                Va.data(), n_a*n_a);
    std::vector<double> S(na2*(size_t)n_v*n_a);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                n_a*n_a, n_v*n_a, (int)aux, 1.0,
                R.AA_RI_M, (int)aux, R.VA_RI_M, (int)aux, 0.0,
                S.data(), n_v*n_a);
    auto Sidx=[&](int p,int q,int a,int r){ return ((size_t)(p*n_a+q)*n_v+a)*n_a+r; };
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
    auto D1=[&](int p,double ea)->double{             // single removed active
        double d = K.e_IP[p]-ea;
#ifndef NDEBUG
        if(ee) assert(std::fabs(d-(epsA-ea)) < 1e-8);
#endif
        return d;
    };

#pragma omp parallel
    {
        const int nt = omp_get_thread_num();
        for(int tu=nt; tu<n_a*n_a; tu+=num_threads){
            const int t=tu/n_a, u=tu%n_a;
            for(int a=0; a<n_v; a++){
                const double ea=e_v[a];
                const double* Vp = Va.data() + (size_t)a*na3; // Vp[(p·n_a+q)·n_a+r]=(ap|qr)
                const double Fta=H_AV[t*n_v+a], Fua=H_AV[u*n_v+a];
                out.g1[t*n_a+u] += Fta*Fua/D1(t,ea);        // collapsed F^c² g1
                for(int v=0; v<n_a; v++)
                for(int w=0; w<n_a; w++){
                    // collapsed g2 (4 terms, all +) + (2e)² δ-image (sum over m)
                    double g = Fta*Vp[(u*n_a+v)*n_a+w]/D1(t,ea)
                             + H_AV[v*n_v+a]*Vp[(w*n_a+t)*n_a+u]/D1(v,ea)
                             + Fua*Vp[(t*n_a+v)*n_a+w]/Dav(t,v,w,ea)
                             + H_AV[w*n_v+a]*Vp[(v*n_a+t)*n_a+u]/Dav(t,v,u,ea);
                    for(int m=0; m<n_a; m++){
                        double Dm=K.e_IP[t]+K.e_IP[v]-K.e_EA[m]-ea;
                        if(m==t) Dm=K.e_IP[v]-ea;
                        if(m==v) Dm=K.e_IP[t]-ea;
#ifndef NDEBUG
                        if(ee) assert(std::fabs(Dm-(epsA-ea))<1e-8);
#endif
                        g += (S[Sidx(t,m,a,v)]*S[Sidx(m,u,a,w)]
                            + Vp[(t*n_a+v)*n_a+m]*S[Sidx(m,w,a,u)])/Dm;
                    }
                    out.g2[(((size_t)t*n_a+u)*n_a+v)*n_a+w] += g;
                }
                // both-double M table (6 terms, all +1) -> raw_av (RAW, un-Hermitized)
                for(int v=0; v<n_a; v++)
                for(int w=0; w<n_a; w++)
                for(int x=0; x<n_a; x++)
                for(int y=0; y<n_a; y++){
                    double m = Vp[(t*n_a+v)*n_a+w]*Vp[(u*n_a+x)*n_a+y]/Dav(t,v,w,ea)
                             + Vp[(t*n_a+x)*n_a+y]*Vp[(u*n_a+v)*n_a+w]/Dav(t,x,y,ea)
                             + Vp[(v*n_a+t)*n_a+u]*Vp[(w*n_a+x)*n_a+y]/Dav(t,v,u,ea)
                             + Vp[(w*n_a+t)*n_a+u]*Vp[(v*n_a+x)*n_a+y]/Dav(v,x,y,ea)
                             + Vp[(x*n_a+t)*n_a+u]*Vp[(y*n_a+v)*n_a+w]/Dav(t,x,u,ea)
                             + Vp[(y*n_a+t)*n_a+u]*Vp[(x*n_a+v)*n_a+w]/Dav(v,x,w,ea);
                    out.raw_av[(((((size_t)t*n_a+u)*n_a+v)*n_a+w)*n_a+x)*n_a+y] += m;
                }
            }
        }
    }
    return 0;
}

// CA (S1 §5.3 both-double M table -> raw_ca; collapsed g2 with OVERALL −1; §5.4
// (2e)² δ-image g2, double-δ g1 with the mandatory ½, F^c-linear δ g1, F^c²
// g1 with the KET label u, scalar E0). External core i. Shared slab
// CAAA[i][p][q][r]=(ip|qr) (CA·AA over aux). OUTPUT pair (t,u) striped over
// threads (disjoint rows); E0 (no integrals) summed serially. Every D collapses
// to e_i−e_A. M table stored RAW — NOT operator-consistent, so its rank-3 gate
// is deferred to the read-off/gauge steps (cert_a1 §1).
int cdas_sf_detail::build_ca(const RI_data& R, const double* eps,
                             int n_c, int n_a, int n_v,
                             const double* H_AV, const double* H_CA,
                             const double* H_CV, const cdas_sf_kernel& K,
                             cdas_sf_tensors& out){
    (void)H_AV; (void)H_CV; (void)n_v;
    if(n_a==0 || n_c==0) return 0;
    const long aux = R.aux_n_ao;
    const double* e_c = eps;
    const size_t na2=(size_t)n_a*n_a, na3=na2*n_a, na6=na3*na3;
    out.raw_ca.assign(na6, 0.0);
#ifndef NDEBUG
    const bool ee = sf_ee_flat(K);
    const double epsA = K.e_EA.empty() ? 0.0 : K.e_EA[0];
#endif
    // CAAA[((i·n_a+p)·n_a+q)·n_a+r]=(ip|qr).
    std::vector<double> CAAA((size_t)n_c*na3);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                n_c*n_a, n_a*n_a, (int)aux, 1.0,
                R.CA_RI_M, (int)aux, R.AA_RI_M, (int)aux, 0.0,
                CAAA.data(), n_a*n_a);
    // D(p;q,r)=e_i+e_IP[p]−e_EA[q]−e_EA[r] and the single-removed D1(q)=e_i−e_EA[q]
    // (kept in a q-only closure over the current core i via the eci argument).
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

#pragma omp parallel
    {
        const int nt=omp_get_thread_num();
        for(int tu=nt; tu<n_a*n_a; tu+=num_threads){
            const int t=tu/n_a, u=tu%n_a;
            for(int i=0;i<n_c;i++){
                const double eci=e_c[i];
                const double* Ci=CAAA.data()+(size_t)i*na3; // Ci[(p·n_a+q)·n_a+r]=(ip|qr)
                const double Fit=H_CA[i*n_a+t], Fiu=H_CA[i*n_a+u];
                out.g1[t*n_a+u] += -Fit*Fiu/D1(eci,u);      // F^c² g1 (KET label u)
                double g1acc=0.0;
                for(int v=0;v<n_a;v++){
                    const double ivtu=Ci[(v*n_a+t)*n_a+u];  // (iv|tu)
                    // F^c-linear δ g1: 2/−1 mirror of CV.3
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
                for(int v=0;v<n_a;v++)
                for(int w=0;w<n_a;w++){
                    // collapsed g2 (4 terms) — accumulated then subtracted (OVERALL −1)
                    double gc = H_CA[i*n_a+u]*Ci[(t*n_a+v)*n_a+w]/D1(eci,u)
                              + H_CA[i*n_a+w]*Ci[(v*n_a+t)*n_a+u]/D1(eci,w)
                              + Fit*Ci[(u*n_a+v)*n_a+w]/Dca(eci,v,u,w)
                              + H_CA[i*n_a+v]*Ci[(w*n_a+t)*n_a+u]/Dca(eci,t,u,w);
                    // (2e)² δ-image g2 (+, sum over active x), mirror of CV g2
                    double gd=0.0;
                    for(int x=0;x<n_a;x++){
                        const double ixtu=Ci[(x*n_a+t)*n_a+u], ixvw=Ci[(x*n_a+v)*n_a+w];
                        gd += (2.0*ixtu*ixvw - Ci[(u*n_a+x)*n_a+t]*ixvw
                               - ixtu*Ci[(v*n_a+x)*n_a+w])/Dca(eci,t,u,x)
                            + (2.0*ixtu*ixvw - Ci[(t*n_a+x)*n_a+u]*ixvw
                               - ixtu*Ci[(w*n_a+x)*n_a+v])/Dca(eci,v,w,x)
                            - Ci[(t*n_a+x)*n_a+w]*Ci[(u*n_a+x)*n_a+v]/Dca(eci,v,u,x)
                            - Ci[(w*n_a+x)*n_a+t]*Ci[(v*n_a+x)*n_a+u]/Dca(eci,t,w,x);
                    }
                    out.g2[(((size_t)t*n_a+u)*n_a+v)*n_a+w] += gd - gc;
                }
                // both-double M table (6 terms, all +1) -> raw_ca (RAW, un-Hermitized)
                for(int v=0;v<n_a;v++)
                for(int w=0;w<n_a;w++)
                for(int x=0;x<n_a;x++)
                for(int y=0;y<n_a;y++){
                    double m = Ci[(t*n_a+u)*n_a+v]*Ci[(w*n_a+x)*n_a+y]/Dca(eci,x,w,y)
                             + Ci[(t*n_a+x)*n_a+y]*Ci[(w*n_a+u)*n_a+v]/Dca(eci,v,u,w)
                             + Ci[(u*n_a+t)*n_a+w]*Ci[(v*n_a+x)*n_a+y]/Dca(eci,t,u,w)
                             + Ci[(v*n_a+t)*n_a+w]*Ci[(u*n_a+x)*n_a+y]/Dca(eci,x,u,y)
                             + Ci[(x*n_a+t)*n_a+w]*Ci[(y*n_a+u)*n_a+v]/Dca(eci,v,u,y)
                             + Ci[(y*n_a+t)*n_a+w]*Ci[(x*n_a+u)*n_a+v]/Dca(eci,t,w,y);
                    out.raw_ca[(((((size_t)t*n_a+u)*n_a+v)*n_a+w)*n_a+x)*n_a+y] += m;
                }
            }
        }
    }
    return 0;
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
