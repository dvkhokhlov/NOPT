# include "blas_link.h"
# include "dsrg_sf_batch.h"
# include <omp.h>
# include <vector>
# include <utility>
# include <cmath>

// STD_SOURCE regularizer (Forte dsrg_source.h). Kept file-local so this TU stays
// independent of the in-core algebra TU. Taylor guard: small = 1e-3, order 4
// (TAYLOR_THRESHOLD = 3), the Forte form the pinned gate uses.
namespace {

inline double dsrg_renorm(double D, double s){ return std::exp(-s*D*D); }

// [1 - exp(-Z^2)] / Z series, Z = sqrt(s) D, k = 1: value = Z - Z^3/2 + Z^5/6 - ...
inline double dsrg_taylor_exp(double Z, int n){
    if(n < 0) return 0.0;
    double value = Z, tmp = Z;
    for(int x=0; x<n-1; ++x){ tmp *= -Z*Z/(x+2); value += tmp; }
    return value;
}

// [1 - exp(-s D^2)] / D with the small-|Z| Taylor guard.
inline double dsrg_renorm_denom(double D, double s){
    const double small = 1e-3;   // 10^-TAYLOR_THRESHOLD, tt = 3
    const int    order = 4;      // floor(0.5*(15/tt + 1)) + 1, tt = 3
    double S = std::sqrt(s), Z = S*D;
    if(std::fabs(Z) < small) return dsrg_taylor_exp(Z, order) * S;
    return (1.0 - std::exp(-Z*Z)) / D;
}

} // namespace

double dsrg_sf_E_CCVV(const double* VC_RI_M,
                      long n_cor, long n_vac, long naux,
                      const double* e_c, const double* e_v,
                      double s, bool ccvv_zero){
    if(n_cor==0 || n_vac==0) return 0.0;
    const long nc=n_cor, nv=n_vac, nQ=naux;

    // triangular core-pair batch (i<=j); OMP over the flat pair list
    std::vector<std::pair<int,int>> ij;
    ij.reserve((size_t)nc*(nc+1)/2);
    for(int i=0;i<nc;i++) for(int j=i;j<nc;j++) ij.emplace_back(i,j);
    const long npair = (long)ij.size();

    std::vector<std::vector<double>> Jab(num_threads, std::vector<double>((size_t)nv*nv));

    double Eout = 0.0;
#ifdef _OPENBLAS
    int ntb = openblas_get_num_threads();
    openblas_set_num_threads(1);
#endif
#ifdef _MKL
    int ntb = mkl_get_max_threads();
    mkl_set_num_threads(1);
#endif
    omp_set_num_threads(num_threads);

#pragma omp parallel for reduction(+:Eout)
    for(long p=0;p<npair;p++){
        int th = omp_get_thread_num();
        double* J = Jab[th].data();
        int i = ij[p].first, j = ij[p].second;

        // J[e,f] = sum_k B(virt e, core i, k) B(virt f, core j, k) = (ie|jf).
        // VC_RI_M is virt-major; per fixed core the (e,k) view is strided by nc*naux.
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                    nv,nv,nQ,1.0,
                    VC_RI_M+(long)i*naux, nc*naux,
                    VC_RI_M+(long)j*naux, nc*naux,
                    0.0, J, nv);

        double loc = 0.0;               // 2J-K formed inline; J stays raw here
        for(long e=0;e<nv;e++)
        for(long f=0;f<nv;f++){
            double D = e_c[i]+e_c[j]-e_v[e]-e_v[f];
            double val = J[e*nv+f];
            if(ccvv_zero) val /= D;
            else          val *= (1.0 + dsrg_renorm(D,s)) * dsrg_renorm_denom(D,s);
            loc += val * (2.0*J[e*nv+f] - J[f*nv+e]);
        }
        Eout += ((i==j) ? 1.0 : 2.0) * loc;
    }

#ifdef _OPENBLAS
    openblas_set_num_threads(ntb);
#endif
#ifdef _MKL
    mkl_set_num_threads(ntb);
#endif
    return Eout;
}

double dsrg_sf_E_CAVV(const double* VC_RI_M, const double* VA_RI_M,
                      long n_cor, long n_act, long n_vac, long naux,
                      const double* e_c, const double* e_a, const double* e_v,
                      const double* L1, double s){
    if(n_cor==0 || n_act==0 || n_vac==0) return 0.0;
    const long nc=n_cor, na=n_act, nv=n_vac, nQ=naux;

    // batch over (core i, virtual c); the free GEMM axes are (virt e, act u)
    std::vector<std::pair<int,int>> ic;
    ic.reserve((size_t)nc*nv);
    for(int i=0;i<nc;i++) for(int c=0;c<nv;c++) ic.emplace_back(i,c);
    const long npair = (long)ic.size();

    std::vector<std::vector<double>> J1(num_threads, std::vector<double>((size_t)nv*na));
    std::vector<std::vector<double>> J2(num_threads, std::vector<double>((size_t)nv*na));
    std::vector<std::vector<double>> JK(num_threads, std::vector<double>((size_t)nv*na));
    std::vector<std::vector<double>> C1(num_threads, std::vector<double>((size_t)na*na, 0.0));

#ifdef _OPENBLAS
    int ntb = openblas_get_num_threads();
    openblas_set_num_threads(1);
#endif
#ifdef _MKL
    int ntb = mkl_get_max_threads();
    mkl_set_num_threads(1);
#endif
    omp_set_num_threads(num_threads);

#pragma omp parallel for
    for(long p=0;p<npair;p++){
        int th = omp_get_thread_num();
        double* j1 = J1[th].data();
        double* j2 = J2[th].data();
        double* jk = JK[th].data();
        double* c1 = C1[th].data();
        int i = ic[p].first, c = ic[p].second;

        // J1[e,u] = (ie|cu): VC strided per fixed core i, VA contiguous per fixed virt c
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                    nv,na,nQ,1.0,
                    VC_RI_M+(long)i*naux,    nc*naux,
                    VA_RI_M+(long)c*na*naux, naux,
                    0.0, j1, na);
        // J2[e,u] = (ic|ue) exchange: full VA block times the B(virt c, core i) column
        cblas_dgemv(CblasRowMajor,CblasNoTrans,
                    nv*na,nQ,1.0,
                    VA_RI_M, naux,
                    VC_RI_M+((long)c*nc+i)*naux, 1,
                    0.0, j2, 1);

        for(long e=0;e<nv;e++)
        for(long u=0;u<na;u++){
            double D = e_c[i]+e_a[u]-e_v[c]-e_v[e];
            jk[e*na+u]  = (2.0*j1[e*na+u] - j2[e*na+u]) * dsrg_renorm_denom(D,s);
            j1[e*na+u] *= (1.0 + dsrg_renorm(D,s));
        }
        // C1[v,u] += sum_e JK[e,v] J1[e,u]
        cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                    na,na,nv,1.0,
                    jk, na, j1, na, 1.0, c1, na);
    }

#ifdef _OPENBLAS
    openblas_set_num_threads(ntb);
#endif
#ifdef _MKL
    mkl_set_num_threads(ntb);
#endif

    std::vector<double> C((size_t)na*na, 0.0);
    for(int th=0;th<num_threads;th++)
        for(long k=0;k<na*na;k++) C[k]+=C1[th][k];

    double Eout=0.0;                          // E = C[v,u] L1[u,v]
    for(long v=0;v<na;v++)
    for(long u=0;u<na;u++)
        Eout += C[v*na+u]*L1[u*na+v];
    return Eout;
}

double dsrg_sf_E_CCAV(const double* VC_RI_M, const double* CA_RI_M,
                      long n_cor, long n_act, long n_vac, long naux,
                      const double* e_c, const double* e_a, const double* e_v,
                      const double* Eta1, double s){
    if(n_cor==0 || n_act==0 || n_vac==0) return 0.0;
    const long nc=n_cor, na=n_act, nv=n_vac, nQ=naux;

    // "ac" block once (shared, read-only): Bac[(u*nc+m)*naux+k] = B(act u, core m, k),
    // a repack of the core-major CA_RI_M so the per-virtual GEMM yields J[u,m,n].
    std::vector<double> Bac((size_t)na*nc*naux);
    for(long u=0;u<na;u++)
    for(long m=0;m<nc;m++)
    for(long k=0;k<naux;k++)
        Bac[(u*nc+m)*naux+k] = CA_RI_M[(m*na+u)*naux+k];

    std::vector<std::vector<double>>  J(num_threads, std::vector<double>((size_t)na*nc*nc));
    std::vector<std::vector<double>> JK(num_threads, std::vector<double>((size_t)na*nc*nc));
    std::vector<std::vector<double>> C1(num_threads, std::vector<double>((size_t)na*na, 0.0));

#ifdef _OPENBLAS
    int ntb = openblas_get_num_threads();
    openblas_set_num_threads(1);
#endif
#ifdef _MKL
    int ntb = mkl_get_max_threads();
    mkl_set_num_threads(1);
#endif
    omp_set_num_threads(num_threads);

#pragma omp parallel for
    for(long c=0;c<nv;c++){
        int th = omp_get_thread_num();
        double* j  =  J[th].data();
        double* jk = JK[th].data();
        double* c1 = C1[th].data();

        // J[u,m,n] = (um|cn): Bac [na*nc x aux] times B(virt c, core n) [nc x aux]
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                    na*nc,nc,nQ,1.0,
                    Bac.data(), naux,
                    VC_RI_M+(long)c*nc*naux, naux,
                    0.0, j, nc);

        for(long u=0;u<na;u++)          // JK[u,m,n] = 2 J[u,m,n] - J[u,n,m] (raw J)
        for(long m=0;m<nc;m++)
        for(long n=0;n<nc;n++)
            jk[(u*nc+m)*nc+n] = 2.0*j[(u*nc+m)*nc+n] - j[(u*nc+n)*nc+m];

        for(long u=0;u<na;u++)          // scale J by (1+exp), JK by [1-exp]/D
        for(long m=0;m<nc;m++)
        for(long n=0;n<nc;n++){
            double D = e_c[m]+e_c[n]-e_v[c]-e_a[u];
            j [(u*nc+m)*nc+n] *= (1.0 + dsrg_renorm(D,s));
            jk[(u*nc+m)*nc+n] *= dsrg_renorm_denom(D,s);
        }
        // C1[v,u] += sum_mn J[v,m,n] JK[u,m,n]
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                    na,na,nc*nc,1.0,
                    j, nc*nc, jk, nc*nc, 1.0, c1, na);
    }

#ifdef _OPENBLAS
    openblas_set_num_threads(ntb);
#endif
#ifdef _MKL
    mkl_set_num_threads(ntb);
#endif

    std::vector<double> C((size_t)na*na, 0.0);
    for(int th=0;th<num_threads;th++)
        for(long k=0;k<na*na;k++) C[k]+=C1[th][k];

    double Eout=0.0;                          // E = C[v,u] Eta1[u,v]
    for(long v=0;v<na;v++)
    for(long u=0;u<na;u++)
        Eout += C[v*na+u]*Eta1[u*na+v];
    return Eout;
}

void dsrg_sf_Hbar1_CAVV(const double* VC_RI_M, const double* VA_RI_M,
                        long n_cor, long n_act, long n_vac, long naux,
                        const double* e_c, const double* e_a, const double* e_v,
                        double s, double* Hbar1){
    if(n_cor==0 || n_act==0 || n_vac==0) return;
    const long nc=n_cor, na=n_act, nv=n_vac, nQ=naux;

    // batch over (core i, virtual c); the free GEMM axes are (virt e, act u)
    std::vector<std::pair<int,int>> ic;
    ic.reserve((size_t)nc*nv);
    for(int i=0;i<nc;i++) for(int c=0;c<nv;c++) ic.emplace_back(i,c);
    const long npair = (long)ic.size();

    std::vector<std::vector<double>> J1(num_threads, std::vector<double>((size_t)nv*na));
    std::vector<std::vector<double>> J2(num_threads, std::vector<double>((size_t)nv*na));
    std::vector<std::vector<double>> JK(num_threads, std::vector<double>((size_t)nv*na));
    std::vector<std::vector<double>> C1(num_threads, std::vector<double>((size_t)na*na, 0.0));

#ifdef _OPENBLAS
    int ntb = openblas_get_num_threads();
    openblas_set_num_threads(1);
#endif
#ifdef _MKL
    int ntb = mkl_get_max_threads();
    mkl_set_num_threads(1);
#endif
    omp_set_num_threads(num_threads);

#pragma omp parallel for
    for(long p=0;p<npair;p++){
        int th = omp_get_thread_num();
        double* j1 = J1[th].data();
        double* j2 = J2[th].data();
        double* jk = JK[th].data();
        double* c1 = C1[th].data();
        int i = ic[p].first, c = ic[p].second;

        // J1[e,u] = (ie|cu): VC strided per fixed core i, VA contiguous per fixed virt c
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                    nv,na,nQ,1.0,
                    VC_RI_M+(long)i*naux,    nc*naux,
                    VA_RI_M+(long)c*na*naux, naux,
                    0.0, j1, na);
        // J2[e,u] = (ic|ue) exchange: full VA block times the B(virt c, core i) column
        cblas_dgemv(CblasRowMajor,CblasNoTrans,
                    nv*na,nQ,1.0,
                    VA_RI_M, naux,
                    VC_RI_M+((long)c*nc+i)*naux, 1,
                    0.0, j2, 1);

        for(long e=0;e<nv;e++)
        for(long u=0;u<na;u++){
            double D = e_c[i]+e_a[u]-e_v[c]-e_v[e];
            jk[e*na+u]  = (2.0*j1[e*na+u] - j2[e*na+u]) * dsrg_renorm_denom(D,s);
            j1[e*na+u] *= (1.0 + dsrg_renorm(D,s));
        }
        // C1[v,u] += sum_e JK[e,v] J1[e,u]
        cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                    na,na,nv,1.0,
                    jk, na, j1, na, 1.0, c1, na);
    }

#ifdef _OPENBLAS
    openblas_set_num_threads(ntb);
#endif
#ifdef _MKL
    mkl_set_num_threads(ntb);
#endif

    std::vector<double> C((size_t)na*na, 0.0);
    for(int th=0;th<num_threads;th++)
        for(long k=0;k<na*na;k++) C[k]+=C1[th][k];

    // Hermitized half-commutator fold: Hbar1[a,b] += 0.5*(C[a,b] + C[b,a])
    for(long a=0;a<na;a++)
    for(long b=0;b<na;b++)
        Hbar1[a*na+b] += 0.5*(C[a*na+b] + C[b*na+a]);
}

void dsrg_sf_Hbar1_CCAV(const double* VC_RI_M, const double* CA_RI_M,
                        long n_cor, long n_act, long n_vac, long naux,
                        const double* e_c, const double* e_a, const double* e_v,
                        double s, double* Hbar1){
    if(n_cor==0 || n_act==0 || n_vac==0) return;
    const long nc=n_cor, na=n_act, nv=n_vac, nQ=naux;

    // "ac" block once (shared, read-only): Bac[(u*nc+m)*naux+k] = B(act u, core m, k),
    // a repack of the core-major CA_RI_M so the per-virtual GEMM yields J[u,m,n].
    std::vector<double> Bac((size_t)na*nc*naux);
    for(long u=0;u<na;u++)
    for(long m=0;m<nc;m++)
    for(long k=0;k<naux;k++)
        Bac[(u*nc+m)*naux+k] = CA_RI_M[(m*na+u)*naux+k];

    std::vector<std::vector<double>>  J(num_threads, std::vector<double>((size_t)na*nc*nc));
    std::vector<std::vector<double>> JK(num_threads, std::vector<double>((size_t)na*nc*nc));
    std::vector<std::vector<double>> C1(num_threads, std::vector<double>((size_t)na*na, 0.0));

#ifdef _OPENBLAS
    int ntb = openblas_get_num_threads();
    openblas_set_num_threads(1);
#endif
#ifdef _MKL
    int ntb = mkl_get_max_threads();
    mkl_set_num_threads(1);
#endif
    omp_set_num_threads(num_threads);

#pragma omp parallel for
    for(long c=0;c<nv;c++){
        int th = omp_get_thread_num();
        double* j  =  J[th].data();
        double* jk = JK[th].data();
        double* c1 = C1[th].data();

        // J[u,m,n] = (um|cn): Bac [na*nc x aux] times B(virt c, core n) [nc x aux]
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                    na*nc,nc,nQ,1.0,
                    Bac.data(), naux,
                    VC_RI_M+(long)c*nc*naux, naux,
                    0.0, j, nc);

        for(long u=0;u<na;u++)          // JK[u,m,n] = 2 J[u,m,n] - J[u,n,m] (raw J)
        for(long m=0;m<nc;m++)
        for(long n=0;n<nc;n++)
            jk[(u*nc+m)*nc+n] = 2.0*j[(u*nc+m)*nc+n] - j[(u*nc+n)*nc+m];

        for(long u=0;u<na;u++)          // scale J by (1+exp), JK by [1-exp]/D
        for(long m=0;m<nc;m++)
        for(long n=0;n<nc;n++){
            double D = e_c[m]+e_c[n]-e_v[c]-e_a[u];
            j [(u*nc+m)*nc+n] *= (1.0 + dsrg_renorm(D,s));
            jk[(u*nc+m)*nc+n] *= dsrg_renorm_denom(D,s);
        }
        // C1[v,u] += sum_mn J[v,m,n] JK[u,m,n]
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                    na,na,nc*nc,1.0,
                    j, nc*nc, jk, nc*nc, 1.0, c1, na);
    }

#ifdef _OPENBLAS
    openblas_set_num_threads(ntb);
#endif
#ifdef _MKL
    mkl_set_num_threads(ntb);
#endif

    std::vector<double> C((size_t)na*na, 0.0);
    for(int th=0;th<num_threads;th++)
        for(long k=0;k<na*na;k++) C[k]+=C1[th][k];

    // opposite-sign Hermitized fold: Hbar1[a,b] -= 0.5*(C[a,b] + C[b,a])
    for(long a=0;a<na;a++)
    for(long b=0;b<na;b++)
        Hbar1[a*na+b] -= 0.5*(C[a*na+b] + C[b*na+a]);
}
