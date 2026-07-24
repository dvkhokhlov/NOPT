// State-specific unrelaxed spin-free DSRG-PT2 algebra core. Assembles the v-tilde
// integral blocks and T1/T2 amplitudes from the RI B-tensors, then the per-class
// <0|[Hbar1,T]|0> energy: Forte's spin-adapted [F,T1]/[F,T2]/[V,T1] closures
// (sadsrg_comm.cc) plus the certified in-core [V,T2] einsum (cert D1) and the DIRECT
// lambda3 assembly (cert D3). Plain integrals, Eta1 = 2I-gamma, S2 = 2T-K.
//
// Block storage is physicist row-major: a 2-body block <p1 p2|p3 p4> is [p1][p2][p3][p4]
// and an amplitude t^{i j}_{a b} is [i][j][a][b]. All hot contractions go through the
// RI aux GEMM / nopt_par_dgemm; small active folds are scalar (OMP over the fat axis).

#include "dsrg_sf_tensors.h"

#include "blas_link.h"
#include "RI.h"
#include "casci_solver.h"

#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cmath>

extern int num_threads;

namespace {

// Raw aux contraction raw[(a),(b)] = sum_k A[a*aux+k] B[b*aux+k] over the RI auxiliary
// index (aux fastest in both slabs). raw is rA x rB row-major (ld = rB).
void contract_aux(const double* A, int rA, const double* B, int rB, long aux, double* raw){
    nopt_par_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, rA, rB, (int)aux,
                   1.0, A, (int)aux, B, (int)aux, 0.0, raw, rB);
}

// out[p][q] = sum_e A[e*dp+p] B[e*dq+q]  (contract the leading external index e).
void contract_ext(const double* A, const double* B, int next, int dp, int dq, double* out){
    nopt_par_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, dp, dq, next,
                   1.0, A, dp, B, dq, 0.0, out, dq);
}

} // namespace

// ---- regularizer (Forte STD source; TAYLOR_THRESHOLD=3 -> small=1e-3, order 4) ----

double dsrg_sf_tensors::rs(double D) const {
    const double sq = std::sqrt(s_flow);
    const double Z = sq * D;
    if(std::fabs(Z) < 1e-3){                 // small-denominator Taylor branch
        double value = Z, tmp = Z;           // (1 - e^{-Z^2})/Z = Z - Z^3/2 + Z^5/6 - ...
        for(int x = 0; x < 3; x++){ tmp *= -Z*Z/(x+2); value += tmp; }
        return value * sq;
    }
    return (1.0 - std::exp(-s_flow*D*D)) / D;
}

double dsrg_sf_tensors::rfac(double D) const {
    return std::exp(-s_flow*D*D);
}

// ---------------------------------------------------------------------------------

void dsrg_sf_tensors::set_problem(int nc, int na, int nv, double s, bool ccvv_source_zero,
                                  const double* eps_c, const double* eps_a, const double* eps_v,
                                  const double* fca, const double* fcv, const double* fav,
                                  const RI_data* Rin){
    n_c = nc; n_a = na; n_v = nv;
    s_flow = s; ccvv_zero = ccvv_source_zero;
    e_c = eps_c; e_a = eps_a; e_v = eps_v;
    f_ca = fca; f_cv = fcv; f_av = fav;
    R = Rin;
}

void dsrg_sf_tensors::set_densities(const double* L1in, const double* GAMMA){
    const int na = n_a;
    const size_t na2 = (size_t)na*na, na4 = na2*na2;
    GAMMA_in = GAMMA;
    L1.assign(L1in, L1in + na2);
    Eta1.assign(na2, 0.0);
    for(int p=0;p<na;p++) for(int q=0;q<na;q++)
        Eta1[p*na+q] = (p==q ? 2.0 : 0.0) - L1[p*na+q];
    // Physicist G2[p][q][r][s] = <p+ q+ s r> from the seam layout GAMMA[t,u,v,w]
    // (which pairs with chemist (tu|vw)): G2_phys[p,q,r,s] = GAMMA[p,r,q,s].
    G2.assign(na4, 0.0);
    for(int p=0;p<na;p++) for(int q=0;q<na;q++)
    for(int r=0;r<na;r++) for(int s=0;s<na;s++)
        G2[((size_t)(p*na+q)*na+r)*na+s] = GAMMA[((size_t)(p*na+r)*na+q)*na+s];
    // SF_L2 cumulant (cert 2.5): SF_L2^{pq}_{rs} = G2 - L1^p_r L1^q_s + 1/2 L1^p_s L1^q_r.
    SF_L2.assign(na4, 0.0);
    for(int p=0;p<na;p++) for(int q=0;q<na;q++)
    for(int r=0;r<na;r++) for(int s=0;s<na;s++)
        SF_L2[((size_t)(p*na+q)*na+r)*na+s] =
              G2[((size_t)(p*na+q)*na+r)*na+s]
            - L1[p*na+r]*L1[q*na+s]
            + 0.5*L1[p*na+s]*L1[q*na+r];
}

// ---- renormalized v-tilde blocks: vt = (1 + e^{-s D^2}) <p1 p2|p3 p4> from B ----

void dsrg_sf_tensors::build_vt(){
    const int na=n_a, nc=n_c, nv=n_v;
    const long aux = R->aux_n_ao;
    const double* VA = R->VA_RI_M;   // row = v*na + a  (virt,act)
    const double* CA = R->CA_RI_M;   // row = c*na + a  (core,act)
    const double* AA = R->AA_RI_M;   // row = a1*na + a2  (only VA/CA/AA feed the vt blocks)

    std::vector<double> raw;

    // vvaa  <e f|x y> = (ex|fy):  slabA=VA(e,x), slabB=VA(f,y)
    Vt_vvaa.assign((size_t)nv*nv*na*na, 0.0);
    raw.assign((size_t)nv*na*nv*na, 0.0);
    contract_aux(VA, nv*na, VA, nv*na, aux, raw.data());
    for(int e=0;e<nv;e++) for(int f=0;f<nv;f++)
    for(int x=0;x<na;x++) for(int y=0;y<na;y++){
        const double D = e_v[e]+e_v[f]-e_a[x]-e_a[y];
        Vt_vvaa[(((size_t)e*nv+f)*na+x)*na+y] =
            raw[(size_t)(e*na+x)*(nv*na)+(f*na+y)]*(1.0+rfac(D));
    }

    // aacc  <w x|i j> = (wi|xj):  slabA=CA(w,i)->row i*na+w, slabB=CA(x,j)->row j*na+x
    Vt_aacc.assign((size_t)na*na*nc*nc, 0.0);
    raw.assign((size_t)nc*na*nc*na, 0.0);
    contract_aux(CA, nc*na, CA, nc*na, aux, raw.data());
    for(int w=0;w<na;w++) for(int x=0;x<na;x++)
    for(int i=0;i<nc;i++) for(int j=0;j<nc;j++){
        const double D = e_a[w]+e_a[x]-e_c[i]-e_c[j];
        Vt_aacc[(((size_t)w*na+x)*nc+i)*nc+j] =
            raw[(size_t)(i*na+w)*(nc*na)+(j*na+x)]*(1.0+rfac(D));
    }

    // avca  <v e|m x> = (vm|ex):  slabA=CA(v,m)->row m*na+v, slabB=VA(e,x)->row e*na+x
    Vt_avca.assign((size_t)na*nv*nc*na, 0.0);
    raw.assign((size_t)nc*na*nv*na, 0.0);
    contract_aux(CA, nc*na, VA, nv*na, aux, raw.data());
    for(int v=0;v<na;v++) for(int e=0;e<nv;e++)
    for(int m=0;m<nc;m++) for(int x=0;x<na;x++){
        const double D = e_a[v]+e_v[e]-e_c[m]-e_a[x];
        Vt_avca[(((size_t)v*nv+e)*nc+m)*na+x] =
            raw[(size_t)(m*na+v)*(nv*na)+(e*na+x)]*(1.0+rfac(D));
    }

    // aaca  <u v|m y> = (um|vy):  slabA=CA(u,m)->row m*na+u, slabB=AA(v,y)->row v*na+y
    Vt_aaca.assign((size_t)na*na*nc*na, 0.0);
    raw.assign((size_t)nc*na*na*na, 0.0);
    contract_aux(CA, nc*na, AA, na*na, aux, raw.data());
    for(int u=0;u<na;u++) for(int v=0;v<na;v++)
    for(int m=0;m<nc;m++) for(int y=0;y<na;y++){
        const double D = e_a[u]+e_a[v]-e_c[m]-e_a[y];
        Vt_aaca[(((size_t)u*na+v)*nc+m)*na+y] =
            raw[(size_t)(m*na+u)*(na*na)+(v*na+y)]*(1.0+rfac(D));
    }

    // avaa  <z a|x y> = (zx|ay):  slabA=AA(z,x)->row z*na+x, slabB=VA(a,y)->row a*na+y
    Vt_avaa.assign((size_t)na*nv*na*na, 0.0);
    raw.assign((size_t)na*na*nv*na, 0.0);
    contract_aux(AA, na*na, VA, nv*na, aux, raw.data());
    for(int z=0;z<na;z++) for(int a=0;a<nv;a++)
    for(int x=0;x<na;x++) for(int y=0;y<na;y++){
        const double D = e_a[z]+e_v[a]-e_a[x]-e_a[y];
        Vt_avaa[(((size_t)z*nv+a)*na+x)*na+y] =
            raw[(size_t)(z*na+x)*(nv*na)+(a*na+y)]*(1.0+rfac(D));
    }

    // vaaa  <e v|w x> = (ew|vx):  slabA=VA(e,w)->row e*na+w, slabB=AA(v,x)->row v*na+x
    Vt_vaaa.assign((size_t)nv*na*na*na, 0.0);
    raw.assign((size_t)nv*na*na*na, 0.0);
    contract_aux(VA, nv*na, AA, na*na, aux, raw.data());
    for(int e=0;e<nv;e++) for(int v=0;v<na;v++)
    for(int w=0;w<na;w++) for(int x=0;x<na;x++){
        const double D = e_v[e]+e_a[v]-e_a[w]-e_a[x];
        Vt_vaaa[(((size_t)e*na+v)*na+w)*na+x] =
            raw[(size_t)(e*na+w)*(na*na)+(v*na+x)]*(1.0+rfac(D));
    }
}

// ---- amplitudes: T2 = <ab|ij> R(D_ij^ab), 7 stored blocks (aaaa internal, zeroed);
//      T1 with the gamma-coupled source; Ft = f + source * e^{-sD^2} (Hermitian). ----

void dsrg_sf_tensors::build_amplitudes(){
    const int na=n_a, nc=n_c, nv=n_v;
    const long aux = R->aux_n_ao;
    const double* VA = R->VA_RI_M;
    const double* CA = R->CA_RI_M;
    const double* VC = R->VC_RI_M;
    const double* AA = R->AA_RI_M;
    std::vector<double> raw;

    // aavv  t^{uv}_{ef} = <ef|uv> R,  <ef|uv>=(eu|fv):  VA(e,u), VA(f,v)
    T2_aavv.assign((size_t)na*na*nv*nv, 0.0);
    raw.assign((size_t)nv*na*nv*na, 0.0);
    contract_aux(VA, nv*na, VA, nv*na, aux, raw.data());
    for(int u=0;u<na;u++) for(int v=0;v<na;v++)
    for(int e=0;e<nv;e++) for(int f=0;f<nv;f++){
        const double D = e_a[u]+e_a[v]-e_v[e]-e_v[f];
        T2_aavv[(((size_t)u*na+v)*nv+e)*nv+f] =
            raw[(size_t)(e*na+u)*(nv*na)+(f*na+v)]*rs(D);
    }

    // ccaa  t^{ij}_{uv} = <uv|ij> R,  <uv|ij>=(ui|vj):  CA(u,i)->row i*na+u, CA(v,j)->j*na+v
    T2_ccaa.assign((size_t)nc*nc*na*na, 0.0);
    raw.assign((size_t)nc*na*nc*na, 0.0);
    contract_aux(CA, nc*na, CA, nc*na, aux, raw.data());
    for(int i=0;i<nc;i++) for(int j=0;j<nc;j++)
    for(int u=0;u<na;u++) for(int v=0;v<na;v++){
        const double D = e_c[i]+e_c[j]-e_a[u]-e_a[v];
        T2_ccaa[(((size_t)i*nc+j)*na+u)*na+v] =
            raw[(size_t)(i*na+u)*(nc*na)+(j*na+v)]*rs(D);
    }

    // caav t^{iv}_{ua} = <ua|iv> R, <ua|iv>=(ui|av): CA(u,i)->i*na+u, VA(a,v)->a*na+v
    T2_caav.assign((size_t)nc*na*na*nv, 0.0);
    raw.assign((size_t)nc*na*nv*na, 0.0);
    contract_aux(CA, nc*na, VA, nv*na, aux, raw.data());
    for(int i=0;i<nc;i++) for(int v=0;v<na;v++)
    for(int u=0;u<na;u++) for(int a=0;a<nv;a++){
        const double D = e_c[i]+e_a[v]-e_a[u]-e_v[a];
        T2_caav[(((size_t)i*na+v)*na+u)*nv+a] =
            raw[(size_t)(i*na+u)*(nv*na)+(a*na+v)]*rs(D);
    }

    // acav t^{vi}_{ua} = <ua|vi> R, <ua|vi>=(uv|ai): AA(u,v)->u*na+v, VC(a,i)->a*nc+i
    T2_acav.assign((size_t)na*nc*na*nv, 0.0);
    raw.assign((size_t)na*na*nv*nc, 0.0);
    contract_aux(AA, na*na, VC, nv*nc, aux, raw.data());
    for(int v=0;v<na;v++) for(int i=0;i<nc;i++)
    for(int u=0;u<na;u++) for(int a=0;a<nv;a++){
        const double D = e_a[v]+e_c[i]-e_a[u]-e_v[a];
        T2_acav[(((size_t)v*nc+i)*na+u)*nv+a] =
            raw[(size_t)(u*na+v)*(nv*nc)+(a*nc+i)]*rs(D);
    }

    // aava t^{uv}_{ey} = <ey|uv> R, <ey|uv>=(eu|yv): VA(e,u)->e*na+u, AA(y,v)->y*na+v
    T2_aava.assign((size_t)na*na*nv*na, 0.0);
    raw.assign((size_t)nv*na*na*na, 0.0);
    contract_aux(VA, nv*na, AA, na*na, aux, raw.data());
    for(int u=0;u<na;u++) for(int v=0;v<na;v++)
    for(int e=0;e<nv;e++) for(int y=0;y<na;y++){
        const double D = e_a[u]+e_a[v]-e_v[e]-e_a[y];
        T2_aava[(((size_t)u*na+v)*nv+e)*na+y] =
            raw[(size_t)(e*na+u)*(na*na)+(y*na+v)]*rs(D);
    }

    // caaa t^{iw}_{uv} = <uv|iw> R, <uv|iw>=(ui|vw): CA(u,i)->i*na+u, AA(v,w)->v*na+w
    T2_caaa.assign((size_t)nc*na*na*na, 0.0);
    raw.assign((size_t)nc*na*na*na, 0.0);
    contract_aux(CA, nc*na, AA, na*na, aux, raw.data());
    for(int i=0;i<nc;i++) for(int w=0;w<na;w++)
    for(int u=0;u<na;u++) for(int v=0;v<na;v++){
        const double D = e_c[i]+e_a[w]-e_a[u]-e_a[v];
        T2_caaa[(((size_t)i*na+w)*na+u)*na+v] =
            raw[(size_t)(i*na+u)*(na*na)+(v*na+w)]*rs(D);
    }
    // aaaa (internal) never formed -> zero by construction.

    // ---- T1: source s^i_a = f^i_a + sum_{ux}(e_x-e_u) gamma_{xu} t^{iu}_{ax};
    //          t^i_a = s^i_a R(e_i-e_a); Ft^i_a = f^i_a + s^i_a e^{-sD^2}. ----
    T1_cv.assign((size_t)nc*nv, 0.0); Ft_cv.assign((size_t)nc*nv, 0.0);
    T1_ca.assign((size_t)nc*na, 0.0); Ft_ca.assign((size_t)nc*na, 0.0);
    T1_av.assign((size_t)na*nv, 0.0); Ft_av.assign((size_t)na*nv, 0.0);

    // The gamma-coupled source uses S2 = 2T-K: 0.5(e_x-e_u) L1[x,u] (2 t - t_swap).
    // cv:  t^{iu}_{ax}=T2_acav[u,i,x,a] (2J); t^{iu}_{xa}=T2_caav[i,u,x,a] (K)
#pragma omp parallel for collapse(2) schedule(static)
    for(int i=0;i<nc;i++) for(int a=0;a<nv;a++){
        double src = f_cv[i*nv+a];
        for(int u=0;u<na;u++) for(int x=0;x<na;x++)
            src += 0.5*(e_a[x]-e_a[u])*L1[x*na+u]*(
                     2.0*T2_acav[(((size_t)u*nc+i)*na+x)*nv+a]
                       - T2_caav[(((size_t)i*na+u)*na+x)*nv+a]);
        const double D = e_c[i]-e_v[a];
        T1_cv[i*nv+a] = ccvv_zero ? src/D : src*rs(D);
        Ft_cv[i*nv+a] = f_cv[i*nv+a] + src*rfac(D);
    }
    // ca:  t^{iu}_{wx}=T2_caaa[i,u,w,x] (2J); t^{iu}_{xw}=T2_caaa[i,u,x,w] (K)
#pragma omp parallel for collapse(2) schedule(static)
    for(int i=0;i<nc;i++) for(int w=0;w<na;w++){
        double src = f_ca[i*na+w];
        for(int u=0;u<na;u++) for(int x=0;x<na;x++)
            src += 0.5*(e_a[x]-e_a[u])*L1[x*na+u]*(
                     2.0*T2_caaa[(((size_t)i*na+u)*na+w)*na+x]
                       - T2_caaa[(((size_t)i*na+u)*na+x)*na+w]);
        const double D = e_c[i]-e_a[w];
        T1_ca[i*na+w] = src*rs(D);
        Ft_ca[i*na+w] = f_ca[i*na+w] + src*rfac(D);
    }
    // av:  t^{vu}_{ax}=T2_aava[v,u,a,x] (2J); t^{vu}_{xa}=t^{uv}_{ax}=T2_aava[u,v,a,x] (K)
#pragma omp parallel for collapse(2) schedule(static)
    for(int v=0;v<na;v++) for(int a=0;a<nv;a++){
        double src = f_av[v*nv+a];
        for(int u=0;u<na;u++) for(int x=0;x<na;x++)
            src += 0.5*(e_a[x]-e_a[u])*L1[x*na+u]*(
                     2.0*T2_aava[(((size_t)v*na+u)*nv+a)*na+x]
                       - T2_aava[(((size_t)u*na+v)*nv+a)*na+x]);
        const double D = e_a[v]-e_v[a];
        T1_av[v*nv+a] = src*rs(D);
        Ft_av[v*nv+a] = f_av[v*nv+a] + src*rfac(D);
    }
}

// ---- E(2): Forte [F,T1]/[F,T2]/[V,T1] closures + certified in-core [V,T2] ----

void dsrg_sf_tensors::compute_e2(){
    const int na=n_a, nc=n_c, nv=n_v;
    const size_t na2=(size_t)na*na, na4=na2*na2;
    auto id4 = [na](int a,int b,int c,int d){ return (((size_t)a*na+b)*na+c)*na+d; };

    // ---- [F,T1]  E = 2 Ft[am] T1[ma] + L1[vu](Ft[ev]T1[ue] - Ft[um]T1[mv]) ----
    // The 2 Ft[am]T1[ma] particle-hole sum runs over BOTH virtual-core and active-core.
    {
        double e=0.0;
        for(int i=0;i<nc;i++) for(int a=0;a<nv;a++)
            e += 2.0*Ft_cv[i*nv+a]*T1_cv[i*nv+a];
        for(int m=0;m<nc;m++) for(int u=0;u<na;u++)
            e += 2.0*Ft_ca[m*na+u]*T1_ca[m*na+u];
        for(int u=0;u<na;u++) for(int v=0;v<na;v++){
            double t1=0.0, t2=0.0;
            for(int a=0;a<nv;a++) t1 += Ft_av[v*nv+a]*T1_av[u*nv+a];   // Ft[ev]T1[ue]
            for(int m=0;m<nc;m++) t2 += Ft_ca[m*na+u]*T1_ca[m*na+v];   // Ft[um]T1[mv]
            e += L1[v*na+u]*(t1 - t2);
        }
        ledger.E_FT1 = e;
    }

    // ---- [F,T2]  E = L2[xyuv](Ft[ex]T2[uvey] - Ft[vm]T2[muyx]) ----
    {
        double e=0.0;
        for(int x=0;x<na;x++) for(int y=0;y<na;y++)
        for(int u=0;u<na;u++) for(int v=0;v<na;v++){
            double a1=0.0, a2=0.0;
            for(int ee=0;ee<nv;ee++)                                  // Ft[ex]T2_aava[u,v,e,y]
                a1 += Ft_av[x*nv+ee]*T2_aava[(((size_t)u*na+v)*nv+ee)*na+y];
            for(int m=0;m<nc;m++)                                     // Ft[vm]T2_caaa[m,u,y,x]
                a2 += Ft_ca[m*na+v]*T2_caaa[(((size_t)m*na+u)*na+y)*na+x];
            e += SF_L2[id4(x,y,u,v)]*(a1 - a2);
        }
        ledger.E_FT2 = e;
    }

    // ---- [V,T1]  E = L2[xyuv](Vt[evxy]T1[ue] - Vt[uvmy]T1[mx]) ----
    {
        double e=0.0;
        for(int x=0;x<na;x++) for(int y=0;y<na;y++)
        for(int u=0;u<na;u++) for(int v=0;v<na;v++){
            double a1=0.0, a2=0.0;
            for(int ee=0;ee<nv;ee++)                                  // Vt_vaaa[e,v,x,y]T1_av[u,e]
                a1 += Vt_vaaa[(((size_t)ee*na+v)*na+x)*na+y]*T1_av[u*nv+ee];
            for(int m=0;m<nc;m++)                                     // Vt_aaca[u,v,m,y]T1_ca[m,x]
                a2 += Vt_aaca[(((size_t)u*na+v)*nc+m)*na+y]*T1_ca[m*na+x];
            e += SF_L2[id4(x,y,u,v)]*(a1 - a2);
        }
        ledger.E_VT1 = e;
    }

    // ------------------------- in-core [V,T2], per class -------------------------
    // Forte H2_T2_C0_T2small spin-free reduction: each class is
    //   0.25 (S2.Vt) (density-products) + 0.5 (T2.Vt) L2,  S2 = 2J-K = 2T2 - T2_swap.
    // The fat (core^2/virt^2/core*virt) index is contracted by GEMM to pure-active
    // tensors; the small active density folds and closings are scalar.
    const long aux = R->aux_n_ao;
    // C[(row)][(col)] = sum_k A[(row),k] B[k,(col)], A is na^2 x K, B is K x na^2.
    auto gemm_aa = [&](const std::vector<double>& A, const std::vector<double>& B, int K,
                       std::vector<double>& C){
        nopt_par_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, (int)na2, (int)na2, K,
                       1.0, A.data(), K, B.data(), (int)na2, 0.0, C.data(), (int)na2);
    };

    // AAVV:  0.25 (S2.Vt) L1 L1 + 0.5 (T2.Vt) L2.  Fat index (e,f) over virt^2.
    {
        std::vector<double> MT(na4, 0.0), MS(na4, 0.0);      // [(u,v)][(w,x)]
        gemm_aa(T2_aavv, Vt_vvaa, nv*nv, MT);                // sum_ef T2[u,v,e,f] Vt[e,f,w,x]
        std::vector<double> S2((size_t)na2*nv*nv);           // S2_aavv[u,v,e,f]=2T2-T2(e<->f)
        for(int u=0;u<na;u++) for(int v=0;v<na;v++)
        for(int e=0;e<nv;e++) for(int f=0;f<nv;f++)
            S2[(((size_t)u*na+v)*nv+e)*nv+f] =
                2.0*T2_aavv[(((size_t)u*na+v)*nv+e)*nv+f] - T2_aavv[(((size_t)u*na+v)*nv+f)*nv+e];
        gemm_aa(S2, Vt_vvaa, nv*nv, MS);
        double e1=0.0, e2=0.0;
        for(int u=0;u<na;u++) for(int v=0;v<na;v++)
        for(int w=0;w<na;w++) for(int x=0;x<na;x++){
            e1 += MS[id4(u,v,w,x)]*L1[x*na+v]*L1[w*na+u];    // 0.25 xv,wu (S2)
            e2 += MT[id4(u,v,w,x)]*SF_L2[id4(u,v,w,x)];      // 0.5  uvwx (T2)
        }
        ledger.aavv_L1 = 0.25*e1;
        ledger.aavv_L2 = 0.5*e2;
    }

    // CCAA:  0.25 (Vt.S2) Eta Eta + 0.5 (Vt.T2) L2.  Fat index (i,j) over core^2.
    {
        std::vector<double> NT(na4, 0.0), NS(na4, 0.0);      // [(w,x)][(u,v)]
        nopt_par_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, (int)na2, (int)na2, nc*nc,
                       1.0, Vt_aacc.data(), nc*nc, T2_ccaa.data(), (int)na2, 0.0, NT.data(), (int)na2);
        std::vector<double> S2((size_t)nc*nc*na2);           // S2_ccaa[i,j,u,v]=2T2-T2(u<->v)
        for(int i=0;i<nc;i++) for(int j=0;j<nc;j++)
        for(int u=0;u<na;u++) for(int v=0;v<na;v++)
            S2[(((size_t)i*nc+j)*na+u)*na+v] =
                2.0*T2_ccaa[(((size_t)i*nc+j)*na+u)*na+v] - T2_ccaa[(((size_t)i*nc+j)*na+v)*na+u];
        nopt_par_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, (int)na2, (int)na2, nc*nc,
                       1.0, Vt_aacc.data(), nc*nc, S2.data(), (int)na2, 0.0, NS.data(), (int)na2);
        double e1=0.0, e2=0.0;
        for(int u=0;u<na;u++) for(int v=0;v<na;v++)
        for(int w=0;w<na;w++) for(int x=0;x<na;x++){
            const size_t wxuv = (size_t)(w*na+x)*na2 + (u*na+v);
            e1 += NS[wxuv]*Eta1[v*na+x]*Eta1[u*na+w];        // 0.25 vx,uw (S2)
            e2 += NT[wxuv]*SF_L2[id4(u,v,w,x)];              // 0.5  uvwx (T2)
        }
        ledger.ccaa_L1 = 0.25*e1;
        ledger.ccaa_L2 = 0.5*e2;
    }

    // CAAV:  v_l1 = 0.5 temp[uxyv] Eta L1 ; v_l2 = temp[uvxy] L2. Needs the avac H2 block
    // (<ue|xm>=(ux|em)) beside the stored avca; fat index (e,m) over virt*core.
    {
        // Vt_avac[u,e,x,m] = <ue|xm> = (ux|em)(1+e^{-sD^2}), D=e_a[u]+e_v[e]-e_a[x]-e_c[m].
        std::vector<double> Vt_avac((size_t)na*nv*na*nc, 0.0);
        {
            std::vector<double> raw((size_t)na*na*nv*nc, 0.0);
            contract_aux(R->AA_RI_M, na*na, R->VC_RI_M, nv*nc, aux, raw.data());   // (ux|em)
            for(int u=0;u<na;u++) for(int e=0;e<nv;e++)
            for(int x=0;x<na;x++) for(int m=0;m<nc;m++){
                const double D = e_a[u]+e_v[e]-e_a[x]-e_c[m];
                Vt_avac[(((size_t)u*nv+e)*na+x)*nc+m] =
                    raw[(size_t)(u*na+x)*(nv*nc)+(e*nc+m)]*(1.0+rfac(D));
            }
        }
        auto Vavca = [&](int p,int e,int m,int x){ return Vt_avca[(((size_t)p*nv+e)*nc+m)*na+x]; };
        auto Vavac = [&](int p,int e,int x,int m){ return Vt_avac[(((size_t)p*nv+e)*na+x)*nc+m]; };
        auto Tcaav = [&](int i,int v,int u,int a){ return T2_caav[(((size_t)i*na+v)*na+u)*nv+a]; };
        auto Tacav = [&](int v,int i,int u,int a){ return T2_acav[(((size_t)v*nc+i)*na+u)*nv+a]; };
        const int nem = nv*nc;
        // Vt gathers laid fat-first [(e,m)][(p,x)]: Wca(p)=<pe|mx>, Wac(p)=<pe|xm>.
        std::vector<double> Wca((size_t)nem*na2, 0.0), Wac((size_t)nem*na2, 0.0);
        for(int e=0;e<nv;e++) for(int m=0;m<nc;m++)
        for(int p=0;p<na;p++) for(int x=0;x<na;x++){
            Wca[(size_t)(e*nc+m)*na2 + (p*na+x)] = Vavca(p,e,m,x);
            Wac[(size_t)(e*nc+m)*na2 + (p*na+x)] = Vavac(p,e,x,m);
        }
        // v_l1: C[(v,x)][(y,u)] from S2_caav / S2_acav.
        std::vector<double> R1((size_t)nem*na2, 0.0), R2((size_t)nem*na2, 0.0);
        for(int e=0;e<nv;e++) for(int m=0;m<nc;m++)
        for(int y=0;y<na;y++) for(int u=0;u<na;u++){
            R1[(size_t)(e*nc+m)*na2 + (y*na+u)] = 2.0*Tcaav(m,y,u,e) - Tacav(y,m,u,e);   // S2_caav[m,y,u,e]
            R2[(size_t)(e*nc+m)*na2 + (y*na+u)] = 2.0*Tacav(y,m,u,e) - Tcaav(m,y,u,e);   // S2_acav[y,m,u,e]
        }
        std::vector<double> C1(na4, 0.0), C2(na4, 0.0);
        contract_ext(Wca.data(), R1.data(), nem, (int)na2, (int)na2, C1.data());
        contract_ext(Wac.data(), R2.data(), nem, (int)na2, (int)na2, C2.data());
        double vl1=0.0;
        for(int v=0;v<na;v++) for(int x=0;x<na;x++)
        for(int y=0;y<na;y++) for(int u=0;u<na;u++){
            const double c = C1[(size_t)(v*na+x)*na2+(y*na+u)] + C2[(size_t)(v*na+x)*na2+(y*na+u)];
            vl1 += c*Eta1[u*na+v]*L1[x*na+y];                // E1[u,v] L1[x,y]
        }
        ledger.caav_L1 = 0.5*vl1;
        // v_l2: temp[uvxy] = RA[(u,x)(v,y)] - RB[(u,x)(v,y)] - RC[(v,x)(u,y)].
        std::vector<double> QA((size_t)nem*na2, 0.0), QB((size_t)nem*na2, 0.0), QC((size_t)nem*na2, 0.0);
        for(int e=0;e<nv;e++) for(int m=0;m<nc;m++)
        for(int a1=0;a1<na;a1++) for(int a2=0;a2<na;a2++){
            QA[(size_t)(e*nc+m)*na2 + (a1*na+a2)] = 2.0*Tacav(a1,m,a2,e) - Tcaav(m,a1,a2,e); // S2_acav[v,m,y,e]
            QB[(size_t)(e*nc+m)*na2 + (a1*na+a2)] = Tacav(a1,m,a2,e);                        // T2_acav[v,m,y,e]
            QC[(size_t)(e*nc+m)*na2 + (a1*na+a2)] = Tcaav(m,a1,a2,e);                        // T2_caav[m,u,y,e]
        }
        std::vector<double> RA(na4, 0.0), RB(na4, 0.0), RC(na4, 0.0);
        contract_ext(Wac.data(), QA.data(), nem, (int)na2, (int)na2, RA.data());   // [(u,x)][(v,y)]
        contract_ext(Wca.data(), QB.data(), nem, (int)na2, (int)na2, RB.data());   // [(u,x)][(v,y)]
        contract_ext(Wca.data(), QC.data(), nem, (int)na2, (int)na2, RC.data());   // [(v,x)][(u,y)]
        double vl2=0.0;
        for(int u=0;u<na;u++) for(int v=0;v<na;v++)
        for(int x=0;x<na;x++) for(int y=0;y<na;y++)
            vl2 += ( RA[(size_t)(u*na+x)*na2+(v*na+y)] - RB[(size_t)(u*na+x)*na2+(v*na+y)]
                     - RC[(size_t)(v*na+x)*na2+(u*na+y)] )*SF_L2[id4(u,v,x,y)];
        ledger.caav_L2 = vl2;
    }

    // CAAA:  fat index m (core). S2_caaa = 2T2 - T2(x<->y). E1 folds the active (w,z) sums.
    {
        auto Tcaaa = [&](int m,int w,int x,int y){ return T2_caaa[(((size_t)m*na+w)*na+x)*na+y]; };
        auto Vaaca = [&](int u,int v,int m,int w){ return Vt_aaca[(((size_t)u*na+v)*nc+m)*na+w]; };
        const int nmz = nc*na;
        std::vector<double> tmp_l2(na4, 0.0), G(na4, 0.0);
        // --- L1 part: tmp_l1[u,x,y,v] = 0.25 sum Vt_aaca[v,z,m,x] S2_caaa[m,y,u,w] E1[w,z] ---
        {
            std::vector<double> L((size_t)na2*nmz, 0.0), Rr((size_t)nmz*na2, 0.0);  // [(v,x)][(m,w)],[(m,w)][(y,u)]
            for(int v=0;v<na;v++) for(int x=0;x<na;x++)
            for(int m=0;m<nc;m++) for(int w=0;w<na;w++){
                double acc=0.0; for(int z=0;z<na;z++) acc += Vaaca(v,z,m,x)*Eta1[w*na+z];
                L[(size_t)(v*na+x)*nmz + (m*na+w)] = acc;                            // sum_z Vt[v,z,m,x] E1[w,z]
            }
            for(int m=0;m<nc;m++) for(int w=0;w<na;w++)
            for(int y=0;y<na;y++) for(int u=0;u<na;u++)
                Rr[(size_t)(m*na+w)*na2 + (y*na+u)] = 2.0*Tcaaa(m,y,u,w) - Tcaaa(m,y,w,u); // S2_caaa[m,y,u,w]
            gemm_aa(L, Rr, nmz, G);                                                  // [(v,x)][(y,u)]
            double e1=0.0;
            for(int v=0;v<na;v++) for(int x=0;x<na;x++)
            for(int y=0;y<na;y++) for(int u=0;u<na;u++)
                e1 += G[(size_t)(v*na+x)*na2+(y*na+u)]*Eta1[u*na+v]*L1[x*na+y];
            ledger.caaa_L1 = 0.25*e1;
        }
        // --- L2 part: four terms accumulate into tmp_l2[u,v,x,y], closed with L2 ---
        // T1: +0.5 (Vt_aaca[u,v,m,w] L1[w,z]) T2_caaa[m,z,x,y]
        {
            std::vector<double> L((size_t)na2*nmz, 0.0);   // [(u,v)][(m,z)]
            for(int u=0;u<na;u++) for(int v=0;v<na;v++)
            for(int m=0;m<nc;m++) for(int z=0;z<na;z++){
                double acc=0.0; for(int w=0;w<na;w++) acc += Vaaca(u,v,m,w)*L1[w*na+z];
                L[(size_t)(u*na+v)*nmz + (m*na+z)] = acc;
            }
            gemm_aa(L, T2_caaa, nmz, G);                    // [(u,v)][(x,y)] (T2_caaa is [(m,z)][(x,y)])
            for(int u=0;u<na;u++) for(int v=0;v<na;v++)
            for(int x=0;x<na;x++) for(int y=0;y<na;y++)
                tmp_l2[id4(u,v,x,y)] += 0.5*G[(size_t)(u*na+v)*na2+(x*na+y)];
        }
        // T2: +0.5 (Vt_aaca[w,u,m,x] E1[w,z]) S2_caaa[m,v,z,y]  -> G[(u,x)][(v,y)]
        // T3: -0.5 (Vt_aaca[u,w,m,x] E1[w,z]) T2_caaa[m,v,z,y]  -> G[(u,x)][(v,y)]
        // T4: -0.5 (Vt_aaca[v,w,m,x] E1[w,z]) T2_caaa[m,u,y,z]  -> G[(v,x)][(u,y)]
        std::vector<double> Rs2((size_t)nmz*na2, 0.0), Rt3((size_t)nmz*na2, 0.0), Rt4((size_t)nmz*na2, 0.0);
        for(int m=0;m<nc;m++) for(int z=0;z<na;z++)
        for(int a1=0;a1<na;a1++) for(int a2=0;a2<na;a2++){
            Rs2[(size_t)(m*na+z)*na2 + (a1*na+a2)] = 2.0*Tcaaa(m,a1,z,a2) - Tcaaa(m,a1,a2,z); // S2_caaa[m,v,z,y]
            Rt3[(size_t)(m*na+z)*na2 + (a1*na+a2)] = Tcaaa(m,a1,z,a2);                        // T2_caaa[m,v,z,y]
            Rt4[(size_t)(m*na+z)*na2 + (a1*na+a2)] = Tcaaa(m,a1,a2,z);                        // T2_caaa[m,u,y,z]
        }
        {
            std::vector<double> L2f((size_t)na2*nmz, 0.0), L3f((size_t)na2*nmz, 0.0), L4f((size_t)na2*nmz, 0.0);
            for(int u=0;u<na;u++) for(int x=0;x<na;x++)
            for(int m=0;m<nc;m++) for(int z=0;z<na;z++){
                double a2=0.0, a3=0.0;
                for(int w=0;w<na;w++){ a2 += Vaaca(w,u,m,x)*Eta1[w*na+z]; a3 += Vaaca(u,w,m,x)*Eta1[w*na+z]; }
                L2f[(size_t)(u*na+x)*nmz + (m*na+z)] = a2;
                L3f[(size_t)(u*na+x)*nmz + (m*na+z)] = a3;
            }
            for(int v=0;v<na;v++) for(int x=0;x<na;x++)
            for(int m=0;m<nc;m++) for(int z=0;z<na;z++){
                double a4=0.0; for(int w=0;w<na;w++) a4 += Vaaca(v,w,m,x)*Eta1[w*na+z];
                L4f[(size_t)(v*na+x)*nmz + (m*na+z)] = a4;
            }
            std::vector<double> G2t(na4,0.0), G3t(na4,0.0), G4t(na4,0.0);
            gemm_aa(L2f, Rs2, nmz, G2t);   // [(u,x)][(v,y)]
            gemm_aa(L3f, Rt3, nmz, G3t);   // [(u,x)][(v,y)]
            gemm_aa(L4f, Rt4, nmz, G4t);   // [(v,x)][(u,y)]
            for(int u=0;u<na;u++) for(int v=0;v<na;v++)
            for(int x=0;x<na;x++) for(int y=0;y<na;y++)
                tmp_l2[id4(u,v,x,y)] += 0.5*G2t[(size_t)(u*na+x)*na2+(v*na+y)]
                                      - 0.5*G3t[(size_t)(u*na+x)*na2+(v*na+y)]
                                      - 0.5*G4t[(size_t)(v*na+x)*na2+(u*na+y)];
        }
        double e2=0.0;
        for(int u=0;u<na;u++) for(int v=0;v<na;v++)
        for(int x=0;x<na;x++) for(int y=0;y<na;y++)
            e2 += tmp_l2[id4(u,v,x,y)]*SF_L2[id4(u,v,x,y)];
        ledger.caaa_L2 = e2;
    }

    // AAAV:  fat index e (virt). S2_aava = 2T2 - T2(u<->v). L1/E1 fold the active (w,z) sums.
    {
        auto Taava = [&](int u,int v,int e,int y){ return T2_aava[(((size_t)u*na+v)*nv+e)*na+y]; };
        auto Vvaaa = [&](int e,int v,int w,int x){ return Vt_vaaa[(((size_t)e*na+v)*na+w)*na+x]; };
        const int nez = nv*na;
        std::vector<double> tmp_l2(na4, 0.0), G(na4, 0.0);
        // --- L1 part: tmp_l1[u,x,y,v] = 0.25 sum Vt_vaaa[e,v,w,x] S2_aava[z,y,e,u] L1[w,z] ---
        {
            std::vector<double> L((size_t)na2*nez, 0.0), Rr((size_t)nez*na2, 0.0);  // [(v,x)][(e,z)],[(e,z)][(y,u)]
            for(int v=0;v<na;v++) for(int x=0;x<na;x++)
            for(int e=0;e<nv;e++) for(int z=0;z<na;z++){
                double acc=0.0; for(int w=0;w<na;w++) acc += Vvaaa(e,v,w,x)*L1[w*na+z];
                L[(size_t)(v*na+x)*nez + (e*na+z)] = acc;
            }
            for(int e=0;e<nv;e++) for(int z=0;z<na;z++)
            for(int y=0;y<na;y++) for(int u=0;u<na;u++)
                Rr[(size_t)(e*na+z)*na2 + (y*na+u)] = 2.0*Taava(z,y,e,u) - Taava(y,z,e,u); // S2_aava[z,y,e,u]
            gemm_aa(L, Rr, nez, G);                                                  // [(v,x)][(y,u)]
            double e1=0.0;
            for(int v=0;v<na;v++) for(int x=0;x<na;x++)
            for(int y=0;y<na;y++) for(int u=0;u<na;u++)
                e1 += G[(size_t)(v*na+x)*na2+(y*na+u)]*Eta1[u*na+v]*L1[x*na+y];
            ledger.aaav_L1 = 0.25*e1;
        }
        // --- L2 part: four terms accumulate into tmp_l2[u,v,x,y], closed with L2 ---
        // T1: +0.5 T2_aava[u,v,e,w] (Vt_vaaa[e,z,x,y] E1[w,z])  -> G[(u,v)][(x,y)]
        {
            std::vector<double> Rr((size_t)nez*na2, 0.0);   // [(e,w)][(x,y)]
            for(int e=0;e<nv;e++) for(int w=0;w<na;w++)
            for(int x=0;x<na;x++) for(int y=0;y<na;y++){
                double acc=0.0; for(int z=0;z<na;z++) acc += Vvaaa(e,z,x,y)*Eta1[w*na+z];
                Rr[(size_t)(e*na+w)*na2 + (x*na+y)] = acc;
            }
            gemm_aa(T2_aava, Rr, nez, G);                   // T2_aava is [(u,v)][(e,w)]
            for(int u=0;u<na;u++) for(int v=0;v<na;v++)
            for(int x=0;x<na;x++) for(int y=0;y<na;y++)
                tmp_l2[id4(u,v,x,y)] += 0.5*G[(size_t)(u*na+v)*na2+(x*na+y)];
        }
        // T2: +0.5 (Vt_vaaa[e,u,w,x] L1[w,z]) S2_aava[z,v,e,y]  -> G[(u,x)][(v,y)]
        // T3: -0.5 (Vt_vaaa[e,u,x,w] L1[w,z]) T2_aava[z,v,e,y]  -> G[(u,x)][(v,y)]
        // T4: -0.5 (Vt_vaaa[e,v,x,w] L1[w,z]) T2_aava[u,z,e,y]  -> G[(v,x)][(u,y)]
        std::vector<double> Rs2((size_t)nez*na2, 0.0), Rt3((size_t)nez*na2, 0.0), Rt4((size_t)nez*na2, 0.0);
        for(int e=0;e<nv;e++) for(int z=0;z<na;z++)
        for(int a1=0;a1<na;a1++) for(int a2=0;a2<na;a2++){
            Rs2[(size_t)(e*na+z)*na2 + (a1*na+a2)] = 2.0*Taava(z,a1,e,a2) - Taava(a1,z,e,a2); // S2_aava[z,v,e,y]
            Rt3[(size_t)(e*na+z)*na2 + (a1*na+a2)] = Taava(z,a1,e,a2);                        // T2_aava[z,v,e,y]
            Rt4[(size_t)(e*na+z)*na2 + (a1*na+a2)] = Taava(a1,z,e,a2);                        // T2_aava[u,z,e,y]
        }
        {
            std::vector<double> L2f((size_t)na2*nez, 0.0), L3f((size_t)na2*nez, 0.0), L4f((size_t)na2*nez, 0.0);
            for(int u=0;u<na;u++) for(int x=0;x<na;x++)
            for(int e=0;e<nv;e++) for(int z=0;z<na;z++){
                double a2=0.0, a3=0.0;
                for(int w=0;w<na;w++){ a2 += Vvaaa(e,u,w,x)*L1[w*na+z]; a3 += Vvaaa(e,u,x,w)*L1[w*na+z]; }
                L2f[(size_t)(u*na+x)*nez + (e*na+z)] = a2;
                L3f[(size_t)(u*na+x)*nez + (e*na+z)] = a3;
            }
            for(int v=0;v<na;v++) for(int x=0;x<na;x++)
            for(int e=0;e<nv;e++) for(int z=0;z<na;z++){
                double a4=0.0; for(int w=0;w<na;w++) a4 += Vvaaa(e,v,x,w)*L1[w*na+z];
                L4f[(size_t)(v*na+x)*nez + (e*na+z)] = a4;
            }
            std::vector<double> G2t(na4,0.0), G3t(na4,0.0), G4t(na4,0.0);
            gemm_aa(L2f, Rs2, nez, G2t);   // [(u,x)][(v,y)]
            gemm_aa(L3f, Rt3, nez, G3t);   // [(u,x)][(v,y)]
            gemm_aa(L4f, Rt4, nez, G4t);   // [(v,x)][(u,y)]
            for(int u=0;u<na;u++) for(int v=0;v<na;v++)
            for(int x=0;x<na;x++) for(int y=0;y<na;y++)
                tmp_l2[id4(u,v,x,y)] += 0.5*G2t[(size_t)(u*na+x)*na2+(v*na+y)]
                                      - 0.5*G3t[(size_t)(u*na+x)*na2+(v*na+y)]
                                      - 0.5*G4t[(size_t)(v*na+x)*na2+(u*na+y)];
        }
        double e2=0.0;
        for(int u=0;u<na;u++) for(int v=0;v<na;v++)
        for(int x=0;x<na;x++) for(int y=0;y<na;y++)
            e2 += tmp_l2[id4(u,v,x,y)]*SF_L2[id4(u,v,x,y)];
        ledger.aaav_L2 = e2;
    }
}

// ---- DIRECT lambda3 (cert 2.6). E3v (virtual) / E3c (core) as separate ledger
//      entries. Solver overlap Omega via the h2caa seam; the remaining pieces are
//      -t*v*D2 plus <=2-body cumulant completions. Hazard H1: one E3c completion
//      line uses S2 = 2T-K, not T2 (marked). ----

void dsrg_sf_tensors::compute_e3(casci_solver* CI, int root){
    const int na=n_a, nc=n_c, nv=n_v;
    const size_t na2=(size_t)na*na, na3=na2*na;
    auto id4 = [na](int a,int b,int c,int d){ return (((size_t)a*na+b)*na+c)*na+d; };
    const int ns = CI->n_states();

    // ---- Omega via the seam. Tensors are [np][n_a^3], flat [(p*na+i1)*na2+i2*na+i3].
    //      virtual leg map: Tbra[e,w,x,y]=Vt_vaaa[e,w,x,y], Tket[e,z,u,v]=T2_aava[u,v,e,z].
    //      core leg map:    Tbra[m,w,x,y]=Vt_aaca[x,y,m,w],  Tket[m,z,u,v]=T2_caaa[m,z,u,v]. ----
    om_v.assign(ns, 0.0);
    om_c.assign(ns, 0.0);
    {
        std::vector<double> Tbra((size_t)nv*na3, 0.0), Tket((size_t)nv*na3, 0.0);
        for(int e=0;e<nv;e++) for(int w=0;w<na;w++) for(int x=0;x<na;x++) for(int y=0;y<na;y++)
            Tbra[((size_t)e*na+w)*na2 + x*na+y] = Vt_vaaa[(((size_t)e*na+w)*na+x)*na+y];
        for(int e=0;e<nv;e++) for(int z=0;z<na;z++) for(int u=0;u<na;u++) for(int v=0;v<na;v++)
            Tket[((size_t)e*na+z)*na2 + u*na+v] = T2_aava[(((size_t)u*na+v)*nv+e)*na+z];
        CI->h2caa_overlap(Tbra.data(), Tket.data(), nv, om_v.data());
    }
    {
        std::vector<double> Tbra((size_t)nc*na3, 0.0), Tket((size_t)nc*na3, 0.0);
        for(int m=0;m<nc;m++) for(int w=0;w<na;w++) for(int x=0;x<na;x++) for(int y=0;y<na;y++)
            Tbra[((size_t)m*na+w)*na2 + x*na+y] = Vt_aaca[(((size_t)x*na+y)*nc+m)*na+w];
        for(int m=0;m<nc;m++) for(int z=0;z<na;z++) for(int u=0;u<na;u++) for(int v=0;v<na;v++)
            Tket[((size_t)m*na+z)*na2 + u*na+v] = T2_caaa[(((size_t)m*na+z)*na+u)*na+v];
        CI->h2caa_overlap(Tbra.data(), Tket.data(), nc, om_c.data());
    }

    // S2 blocks needed by the completions: S2[i,j,a,b] = 2 T2[i,j,a,b] - T2[j,i,a,b].
    // aava: S2_aava[u,v,e,z] = 2 T2_aava[u,v,e,z] - T2_aava[v,u,e,z].
    // caaa: S2_caaa[m,w,x,y] = 2 T2_caaa[m,w,x,y] - T2_caaa[m,w,y,x] (particles both active).
    std::vector<double> S2_aava(T2_aava.size()), S2_caaa(T2_caaa.size());
    for(int u=0;u<na;u++) for(int v=0;v<na;v++) for(int e=0;e<nv;e++) for(int z=0;z<na;z++)
        S2_aava[(((size_t)u*na+v)*nv+e)*na+z] =
            2.0*T2_aava[(((size_t)u*na+v)*nv+e)*na+z] - T2_aava[(((size_t)v*na+u)*nv+e)*na+z];
    for(int m=0;m<nc;m++) for(int w=0;w<na;w++) for(int x=0;x<na;x++) for(int y=0;y<na;y++)
        S2_caaa[(((size_t)m*na+w)*na+x)*na+y] =
            2.0*T2_caaa[(((size_t)m*na+w)*na+x)*na+y] - T2_caaa[(((size_t)m*na+w)*na+y)*na+x];

    // Convenience accessors (physicist / amplitude storage).
    auto H2v = [&](int e,int a1,int a2,int a3)->double{     // Vt_vaaa[e,v,w,x]
        return Vt_vaaa[(((size_t)e*na+a1)*na+a2)*na+a3]; };
    auto Ta  = [&](int u,int v,int e,int z)->double{        // T2_aava
        return T2_aava[(((size_t)u*na+v)*nv+e)*na+z]; };
    auto Sa  = [&](int u,int v,int e,int z)->double{        // S2_aava
        return S2_aava[(((size_t)u*na+v)*nv+e)*na+z]; };
    auto H2c = [&](int a1,int a2,int m,int a3)->double{     // Vt_aaca[u,v,m,y]
        return Vt_aaca[(((size_t)a1*na+a2)*nc+m)*na+a3]; };
    auto Tc  = [&](int m,int w,int x,int y)->double{        // T2_caaa
        return T2_caaa[(((size_t)m*na+w)*na+x)*na+y]; };
    auto Sc  = [&](int m,int w,int x,int y)->double{        // S2_caaa
        return S2_caaa[(((size_t)m*na+w)*na+x)*na+y]; };
    auto g2  = [&](int x,int y,int u,int v)->double{ return G2[id4(x,y,u,v)]; };
    auto l2  = [&](int x,int y,int u,int v)->double{ return SF_L2[id4(x,y,u,v)]; };
    auto l1  = [&](int p,int q)->double{ return L1[p*na+q]; };

    // ---------------------------- E3v (virtual branch) ----------------------------
    // All completion terms first fold L1 into H2v over its two ket actives, contract e
    // against T2/S2 into a pure-active tensor, then close with G2/L2. n_a^4 intermediates.
    double E3v = om_v[root];

    // term1: - H2v[e,z,x,y] T2[u,v,e,z] G2[x,y,u,v]
    {
        std::vector<double> A((size_t)nv*na*na2, 0.0), B((size_t)nv*na*na2, 0.0);  // [(e,z)][(x,y)],[(e,z)][(u,v)]
        for(int e=0;e<nv;e++) for(int z=0;z<na;z++)
        for(int x=0;x<na;x++) for(int y=0;y<na;y++)
            A[(size_t)(e*na+z)*na2 + (x*na+y)] = H2v(e,z,x,y);
        for(int e=0;e<nv;e++) for(int z=0;z<na;z++)
        for(int u=0;u<na;u++) for(int v=0;v<na;v++)
            B[(size_t)(e*na+z)*na2 + (u*na+v)] = Ta(u,v,e,z);
        std::vector<double> M(na2*na2, 0.0);
        contract_ext(A.data(), B.data(), nv*na, (int)na2, (int)na2, M.data());
        double e=0.0;
        for(int x=0;x<na;x++) for(int y=0;y<na;y++)
        for(int u=0;u<na;u++) for(int v=0;v<na;v++)
            e += M[(size_t)(x*na+y)*na2+(u*na+v)]*g2(x,y,u,v);
        E3v -= e;
    }
    // Fold helper: HH[e][x] = sum_{w,y} (H2v[e,w,x,y] c1 + H2v[e,w,y,x] c2) l1[y,w] then
    // contract e with a T2/S2 block -> N[x][(z,u,v)], close with L2/G2. Because the six
    // completion terms differ in which H2v legs carry L1 and which density closes, they
    // are written out; each keeps intermediates at n_a^4.
    auto contract_e = [&](const std::vector<double>& HH_ex,
                          const double* amp, bool amp_is_S)->std::vector<double>{
        // amp is T2_aava(amp_is_S=false) or S2_aava(true): amp[u,v,e,z]. Returns
        // N[x][(z,u,v)] = sum_e HH_ex[e,x] amp[u,v,e,z].
        std::vector<double> Ar((size_t)nv*na, 0.0), Br((size_t)nv*na3, 0.0);
        for(int e=0;e<nv;e++) for(int x=0;x<na;x++)
            Ar[(size_t)e*na + x] = HH_ex[(size_t)e*na + x];
        for(int e=0;e<nv;e++) for(int z=0;z<na;z++)
        for(int u=0;u<na;u++) for(int v=0;v<na;v++)
            Br[(size_t)e*na3 + ((z*na+u)*na+v)] = amp[(((size_t)u*na+v)*nv+e)*na+z];
        std::vector<double> N((size_t)na*na3, 0.0);   // [x][(z,u,v)]
        contract_ext(Ar.data(), Br.data(), nv, na, (int)na3, N.data());
        (void)amp_is_S; return N;
    };

    // term3: - (H2v[e,w,x,y] l1[y,w] - 1/2 H2v[e,w,y,x] l1[y,w]) T2[u,v,e,z] G2[x,z,u,v]
    {
        std::vector<double> HH((size_t)nv*na, 0.0);   // [e][x]
        for(int e=0;e<nv;e++) for(int x=0;x<na;x++){
            double acc=0.0;
            for(int w=0;w<na;w++) for(int y=0;y<na;y++)
                acc += (H2v(e,w,x,y) - 0.5*H2v(e,w,y,x))*l1(y,w);
            HH[(size_t)e*na+x]=acc;
        }
        std::vector<double> N = contract_e(HH, T2_aava.data(), false);   // N[x][(z,u,v)]
        double e=0.0;
        for(int x=0;x<na;x++) for(int z=0;z<na;z++)
        for(int u=0;u<na;u++) for(int v=0;v<na;v++)
            e += N[(size_t)x*na3+((z*na+u)*na+v)]*g2(x,z,u,v);
        E3v -= e;
    }
    // term4: - H2v[e,w,x,y] (1/2 S2[u,v,e,z] l1[z,v]) L2[x,y,u,w]
    {
        // fold l1[z,v] into S2 over z,v -> SF[e][u] (leading e for the contract_e GEMM);
        // SF[e,u] = sum_{v,z} S2[u,v,e,z] l1[z,v].
        std::vector<double> SF((size_t)nv*na, 0.0);   // [e][u]
        for(int e=0;e<nv;e++) for(int u=0;u<na;u++){
            double acc=0.0;
            for(int v=0;v<na;v++) for(int z=0;z<na;z++) acc += Sa(u,v,e,z)*l1(z,v);
            SF[(size_t)e*na+u]=acc;
        }
        // C[u][w,x,y] = sum_e H2v[e,w,x,y] SF[e,u]; close with L2[x,y,u,w] (0.5 factor).
        std::vector<double> Hr((size_t)nv*na3, 0.0);
        for(int e=0;e<nv;e++) for(int w=0;w<na;w++) for(int x=0;x<na;x++) for(int y=0;y<na;y++)
            Hr[(size_t)e*na3 + ((w*na+x)*na+y)] = H2v(e,w,x,y);
        std::vector<double> C((size_t)na*na3, 0.0);   // [u][(w,x,y)]
        contract_ext(SF.data(), Hr.data(), nv, na, (int)na3, C.data());
        double e=0.0;
        for(int u=0;u<na;u++) for(int w=0;w<na;w++)
        for(int x=0;x<na;x++) for(int y=0;y<na;y++)
            e += C[(size_t)u*na3+((w*na+x)*na+y)]*l2(x,y,u,w);
        E3v -= 0.5*e;
    }
    // term5: - 1/2 (H2v[e,w,x,y] l1[x,u]) S2[u,v,e,z] L2[y,z,w,v]
    {
        std::vector<double> HH((size_t)nv*na, 0.0);   // [e][(w,y)] folded? need x->u via l1
        // Fold l1[x,u] into H2v over x -> HL[e,w,u,y]; then contract e with S2, close L2.
        std::vector<double> HL((size_t)nv*na3, 0.0);  // [e][(w,u,y)]
        for(int e=0;e<nv;e++) for(int w=0;w<na;w++) for(int u=0;u<na;u++) for(int y=0;y<na;y++){
            double acc=0.0; for(int x=0;x<na;x++) acc += H2v(e,w,x,y)*l1(x,u);
            HL[(size_t)e*na3 + ((w*na+u)*na+y)] = acc;
        }
        // D[(w,u,y)][(v,z)] = sum_e HL[e,w,u,y] S2[u',v,e,z]? u is shared: sum over e only,
        // keep u,v,z,w,y. Build M[w,u,y | z,v] = sum_e HL[e,w,u,y] S2[u,v,e,z].
        double e=0.0;
        for(int u=0;u<na;u++){
            // slab over fixed u: HLu[e][(w,y)], Su[e][(v,z)]
            std::vector<double> HLu((size_t)nv*na2, 0.0), Su((size_t)nv*na2, 0.0);
            for(int ee=0;ee<nv;ee++) for(int w=0;w<na;w++) for(int y=0;y<na;y++)
                HLu[(size_t)ee*na2 + (w*na+y)] = HL[(size_t)ee*na3 + ((w*na+u)*na+y)];
            for(int ee=0;ee<nv;ee++) for(int v=0;v<na;v++) for(int z=0;z<na;z++)
                Su[(size_t)ee*na2 + (v*na+z)] = Sa(u,v,ee,z);
            std::vector<double> Mu(na2*na2, 0.0);   // [(w,y)][(v,z)]
            contract_ext(HLu.data(), Su.data(), nv, (int)na2, (int)na2, Mu.data());
            for(int w=0;w<na;w++) for(int y=0;y<na;y++)
            for(int v=0;v<na;v++) for(int z=0;z<na;z++)
                e += Mu[(size_t)(w*na+y)*na2+(v*na+z)]*l2(y,z,w,v);
        }
        E3v -= 0.5*e;
    }
    // term6: + 1/2 (H2v[e,w,x,y] l1[y,u]) T2[u,v,e,z] L2[x,z,w,v]
    {
        std::vector<double> HL((size_t)nv*na3, 0.0);  // [e][(w,u,x)] = sum_y H2v[e,w,x,y] l1[y,u]
        for(int e=0;e<nv;e++) for(int w=0;w<na;w++) for(int u=0;u<na;u++) for(int x=0;x<na;x++){
            double acc=0.0; for(int y=0;y<na;y++) acc += H2v(e,w,x,y)*l1(y,u);
            HL[(size_t)e*na3 + ((w*na+u)*na+x)] = acc;
        }
        double e=0.0;
        for(int u=0;u<na;u++){
            std::vector<double> HLu((size_t)nv*na2, 0.0), Tu((size_t)nv*na2, 0.0);
            for(int ee=0;ee<nv;ee++) for(int w=0;w<na;w++) for(int x=0;x<na;x++)
                HLu[(size_t)ee*na2 + (w*na+x)] = HL[(size_t)ee*na3 + ((w*na+u)*na+x)];
            for(int ee=0;ee<nv;ee++) for(int v=0;v<na;v++) for(int z=0;z<na;z++)
                Tu[(size_t)ee*na2 + (v*na+z)] = Ta(u,v,ee,z);
            std::vector<double> Mu(na2*na2, 0.0);   // [(w,x)][(v,z)]
            contract_ext(HLu.data(), Tu.data(), nv, (int)na2, (int)na2, Mu.data());
            for(int w=0;w<na;w++) for(int x=0;x<na;x++)
            for(int v=0;v<na;v++) for(int z=0;z<na;z++)
                e += Mu[(size_t)(w*na+x)*na2+(v*na+z)]*l2(x,z,w,v);
        }
        E3v += 0.5*e;
    }
    // term7: + 1/2 (H2v[e,w,x,y] l1[y,v]) T2[u,v,e,z] L2[x,z,u,w]
    {
        std::vector<double> HL((size_t)nv*na3, 0.0);  // [e][(w,v,x)] = sum_y H2v[e,w,x,y] l1[y,v]
        for(int e=0;e<nv;e++) for(int w=0;w<na;w++) for(int v=0;v<na;v++) for(int x=0;x<na;x++){
            double acc=0.0; for(int y=0;y<na;y++) acc += H2v(e,w,x,y)*l1(y,v);
            HL[(size_t)e*na3 + ((w*na+v)*na+x)] = acc;
        }
        double e=0.0;
        for(int v=0;v<na;v++){
            std::vector<double> HLv((size_t)nv*na2, 0.0), Tv((size_t)nv*na2, 0.0);
            for(int ee=0;ee<nv;ee++) for(int w=0;w<na;w++) for(int x=0;x<na;x++)
                HLv[(size_t)ee*na2 + (w*na+x)] = HL[(size_t)ee*na3 + ((w*na+v)*na+x)];
            for(int ee=0;ee<nv;ee++) for(int u=0;u<na;u++) for(int z=0;z<na;z++)
                Tv[(size_t)ee*na2 + (u*na+z)] = Ta(u,v,ee,z);
            std::vector<double> Mv(na2*na2, 0.0);   // [(w,x)][(u,z)]
            contract_ext(HLv.data(), Tv.data(), nv, (int)na2, (int)na2, Mv.data());
            for(int w=0;w<na;w++) for(int x=0;x<na;x++)
            for(int u=0;u<na;u++) for(int z=0;z<na;z++)
                e += Mv[(size_t)(w*na+x)*na2+(u*na+z)]*l2(x,z,u,w);
        }
        E3v += 0.5*e;
    }
    // term8: + 1/2 (H2v[e,w,x,y] l1[z,w]) T2[u,v,e,z] G2[x,y,u,v]
    {
        // fold l1[z,w] into T2 over w? w is on H2v; z on T2. Fold l1 into T2 over z -> TL[u,v,e,w]
        std::vector<double> TL((size_t)na2*nv*na, 0.0);  // [u][v][e][w]
        for(int u=0;u<na;u++) for(int v=0;v<na;v++) for(int e=0;e<nv;e++) for(int w=0;w<na;w++){
            double acc=0.0; for(int z=0;z<na;z++) acc += Ta(u,v,e,z)*l1(z,w);
            TL[(((size_t)u*na+v)*nv+e)*na+w] = acc;
        }
        // M[(x,y)][(u,v)] = sum_{e,w} H2v[e,w,x,y] TL[u,v,e,w]
        std::vector<double> Hr((size_t)nv*na*na2, 0.0), Tr((size_t)nv*na*na2, 0.0);
        for(int e=0;e<nv;e++) for(int w=0;w<na;w++) for(int x=0;x<na;x++) for(int y=0;y<na;y++)
            Hr[(size_t)(e*na+w)*na2 + (x*na+y)] = H2v(e,w,x,y);
        for(int e=0;e<nv;e++) for(int w=0;w<na;w++) for(int u=0;u<na;u++) for(int v=0;v<na;v++)
            Tr[(size_t)(e*na+w)*na2 + (u*na+v)] = TL[(((size_t)u*na+v)*nv+e)*na+w];
        std::vector<double> M(na2*na2, 0.0);
        contract_ext(Hr.data(), Tr.data(), nv*na, (int)na2, (int)na2, M.data());
        double e=0.0;
        for(int x=0;x<na;x++) for(int y=0;y<na;y++)
        for(int u=0;u<na;u++) for(int v=0;v<na;v++)
            e += M[(size_t)(x*na+y)*na2+(u*na+v)]*g2(x,y,u,v);
        E3v += 0.5*e;
    }
    ledger.E3v = E3v;

    // ---------------------------- E3c (core branch) ----------------------------
    // H2c = Vt_aaca[u,v,m,y]; Tc/Sc = T2_caaa/S2_caaa[m,w,x,y]; external m (core).
    double E3c = -om_c[root];

    // term1: + H2c[u,v,m,z] T2[m,z,x,y] G2[x,y,u,v]
    {
        std::vector<double> A((size_t)nc*na*na2, 0.0), B((size_t)nc*na*na2, 0.0);  // [(m,z)][(u,v)],[(m,z)][(x,y)]
        for(int u=0;u<na;u++) for(int v=0;v<na;v++) for(int m=0;m<nc;m++) for(int z=0;z<na;z++)
            A[(size_t)(m*na+z)*na2 + (u*na+v)] = H2c(u,v,m,z);
        for(int m=0;m<nc;m++) for(int z=0;z<na;z++) for(int x=0;x<na;x++) for(int y=0;y<na;y++)
            B[(size_t)(m*na+z)*na2 + (x*na+y)] = Tc(m,z,x,y);
        std::vector<double> M(na2*na2, 0.0);   // [(u,v)][(x,y)]
        contract_ext(A.data(), B.data(), nc*na, (int)na2, (int)na2, M.data());
        double e=0.0;
        for(int u=0;u<na;u++) for(int v=0;v<na;v++)
        for(int x=0;x<na;x++) for(int y=0;y<na;y++)
            e += M[(size_t)(u*na+v)*na2+(x*na+y)]*g2(x,y,u,v);
        E3c += e;
    }
    // term3: + (H2c[u,v,m,z] l1[z,v] - 1/2 H2c[v,u,m,z] l1[z,v]) T2[m,w,x,y] L2[x,y,u,w]
    {
        std::vector<double> HH((size_t)nc*na, 0.0);   // [m][u] (leading m for the GEMM)
        for(int u=0;u<na;u++) for(int m=0;m<nc;m++){
            double acc=0.0;
            for(int v=0;v<na;v++) for(int z=0;z<na;z++)
                acc += (H2c(u,v,m,z) - 0.5*H2c(v,u,m,z))*l1(z,v);
            HH[(size_t)m*na+u]=acc;
        }
        // C[u][(w,x,y)] = sum_m HH[m,u] T2_caaa[m,w,x,y]; close L2[x,y,u,w].
        std::vector<double> Tr((size_t)nc*na3, 0.0);
        for(int m=0;m<nc;m++) for(int w=0;w<na;w++) for(int x=0;x<na;x++) for(int y=0;y<na;y++)
            Tr[(size_t)m*na3 + ((w*na+x)*na+y)] = Tc(m,w,x,y);
        std::vector<double> C((size_t)na*na3, 0.0);
        contract_ext(HH.data(), Tr.data(), nc, na, (int)na3, C.data());
        double e=0.0;
        for(int u=0;u<na;u++) for(int w=0;w<na;w++)
        for(int x=0;x<na;x++) for(int y=0;y<na;y++)
            e += C[(size_t)u*na3+((w*na+x)*na+y)]*l2(x,y,u,w);
        E3c += e;
    }
    // term4: + 1/2 H2c[u,v,m,z] (S2[m,w,x,y] l1[y,w]) G2[x,z,u,v]
    {
        std::vector<double> SF((size_t)nc*na, 0.0);   // [m][x] = sum_{w,y} S2_caaa[m,w,x,y] l1[y,w]
        for(int m=0;m<nc;m++) for(int x=0;x<na;x++){
            double acc=0.0; for(int w=0;w<na;w++) for(int y=0;y<na;y++) acc += Sc(m,w,x,y)*l1(y,w);
            SF[(size_t)m*na+x]=acc;
        }
        // C[(u,v,z)][x] = sum_m H2c[u,v,m,z] SF[m,x]; close G2[x,z,u,v].
        std::vector<double> Hr((size_t)nc*na3, 0.0);
        for(int u=0;u<na;u++) for(int v=0;v<na;v++) for(int m=0;m<nc;m++) for(int z=0;z<na;z++)
            Hr[(size_t)m*na3 + ((u*na+v)*na+z)] = H2c(u,v,m,z);
        std::vector<double> C((size_t)na3*na, 0.0);   // [(u,v,z)][x]
        contract_ext(Hr.data(), SF.data(), nc, (int)na3, na, C.data());
        double e=0.0;
        for(int u=0;u<na;u++) for(int v=0;v<na;v++) for(int z=0;z<na;z++)
        for(int x=0;x<na;x++)
            e += C[(size_t)((u*na+v)*na+z)*na + x]*g2(x,z,u,v);
        E3c += 0.5*e;
    }
    // term5: + 1/2 (H2c[u,v,m,z] l1[x,u]) S2[m,w,x,y] L2[y,z,w,v]   <-- S2, NOT T2 (H1)
    {
        std::vector<double> HL((size_t)nc*na3, 0.0);  // [m][(v,z,x)] = sum_u H2c[u,v,m,z] l1[x,u]
        for(int m=0;m<nc;m++) for(int v=0;v<na;v++) for(int z=0;z<na;z++) for(int x=0;x<na;x++){
            double acc=0.0; for(int u=0;u<na;u++) acc += H2c(u,v,m,z)*l1(x,u);
            HL[(size_t)m*na3 + ((v*na+z)*na+x)] = acc;
        }
        double e=0.0;
        for(int x=0;x<na;x++){
            std::vector<double> HLx((size_t)nc*na2, 0.0), Sx((size_t)nc*na2, 0.0);
            for(int m=0;m<nc;m++) for(int v=0;v<na;v++) for(int z=0;z<na;z++)
                HLx[(size_t)m*na2 + (v*na+z)] = HL[(size_t)m*na3 + ((v*na+z)*na+x)];
            for(int m=0;m<nc;m++) for(int w=0;w<na;w++) for(int y=0;y<na;y++)
                Sx[(size_t)m*na2 + (w*na+y)] = Sc(m,w,x,y);
            std::vector<double> Mx(na2*na2, 0.0);   // [(v,z)][(w,y)]
            contract_ext(HLx.data(), Sx.data(), nc, (int)na2, (int)na2, Mx.data());
            for(int v=0;v<na;v++) for(int z=0;z<na;z++)
            for(int w=0;w<na;w++) for(int y=0;y<na;y++)
                e += Mx[(size_t)(v*na+z)*na2+(w*na+y)]*l2(y,z,w,v);
        }
        E3c += 0.5*e;
    }
    // term6: -0.5 (H2c[u,v,m,z] l1[y,v]) T2[m,w,x,y] L2[x,z,u,w]  (fold vs T2 4th index y)
    // term7: -0.5 (H2c[u,v,m,z] l1[x,v]) T2[m,w,x,y] L2[y,z,w,u]  (fold vs T2 3rd index x)
    // Pr[(m,a)][(u,z)] = sum_v H2c[u,v,m,z] l1[a,v] serves both; only the T2 leg contracted
    // and the L2 closing differ.
    {
        std::vector<double> Pr((size_t)nc*na*na2, 0.0);   // [(m,a)][(u,z)]
        for(int m=0;m<nc;m++) for(int a=0;a<na;a++)
        for(int u=0;u<na;u++) for(int z=0;z<na;z++){
            double acc=0.0; for(int v=0;v<na;v++) acc += H2c(u,v,m,z)*l1(a,v);
            Pr[(size_t)(m*na+a)*na2 + (u*na+z)] = acc;
        }
        std::vector<double> T6((size_t)nc*na*na2, 0.0), T7((size_t)nc*na*na2, 0.0);
        for(int m=0;m<nc;m++) for(int w=0;w<na;w++) for(int x=0;x<na;x++) for(int y=0;y<na;y++){
            T6[(size_t)(m*na+y)*na2 + (w*na+x)] = Tc(m,w,x,y);   // Pr's a<->y, free (w,x)
            T7[(size_t)(m*na+x)*na2 + (w*na+y)] = Tc(m,w,x,y);   // Pr's a<->x, free (w,y)
        }
        std::vector<double> t6(na2*na2, 0.0), t7(na2*na2, 0.0);  // [(u,z)][(w,x)] / [(u,z)][(w,y)]
        contract_ext(Pr.data(), T6.data(), nc*na, (int)na2, (int)na2, t6.data());
        contract_ext(Pr.data(), T7.data(), nc*na, (int)na2, (int)na2, t7.data());
        double e6=0.0, e7=0.0;
        for(int u=0;u<na;u++) for(int z=0;z<na;z++)
        for(int w=0;w<na;w++) for(int x=0;x<na;x++)
            e6 += t6[(size_t)(u*na+z)*na2+(w*na+x)]*l2(x,z,u,w);
        for(int u=0;u<na;u++) for(int z=0;z<na;z++)
        for(int w=0;w<na;w++) for(int y=0;y<na;y++)
            e7 += t7[(size_t)(u*na+z)*na2+(w*na+y)]*l2(y,z,w,u);
        E3c -= 0.5*e6;
        E3c -= 0.5*e7;
    }
    // term8: - 1/2 (H2c[u,v,m,z] l1[z,w]) T2[m,w,x,y] G2[x,y,u,v]
    {
        std::vector<double> HL((size_t)nc*na3, 0.0);  // [m][(u,v,w)] = sum_z H2c[u,v,m,z] l1[z,w]
        for(int m=0;m<nc;m++) for(int u=0;u<na;u++) for(int v=0;v<na;v++) for(int w=0;w<na;w++){
            double acc=0.0; for(int z=0;z<na;z++) acc += H2c(u,v,m,z)*l1(z,w);
            HL[(size_t)m*na3 + ((u*na+v)*na+w)] = acc;
        }
        // M[(u,v)][(x,y)] = sum_{m,w} HL[m,u,v,w] T2_caaa[m,w,x,y]
        std::vector<double> Hr((size_t)nc*na*na2, 0.0), Tr((size_t)nc*na*na2, 0.0);
        for(int m=0;m<nc;m++) for(int w=0;w<na;w++) for(int u=0;u<na;u++) for(int v=0;v<na;v++)
            Hr[(size_t)(m*na+w)*na2 + (u*na+v)] = HL[(size_t)m*na3 + ((u*na+v)*na+w)];
        for(int m=0;m<nc;m++) for(int w=0;w<na;w++) for(int x=0;x<na;x++) for(int y=0;y<na;y++)
            Tr[(size_t)(m*na+w)*na2 + (x*na+y)] = Tc(m,w,x,y);
        std::vector<double> M(na2*na2, 0.0);
        contract_ext(Hr.data(), Tr.data(), nc*na, (int)na2, (int)na2, M.data());
        double e=0.0;
        for(int u=0;u<na;u++) for(int v=0;v<na;v++)
        for(int x=0;x<na;x++) for(int y=0;y<na;y++)
            e += M[(size_t)(u*na+v)*na2+(x*na+y)]*g2(x,y,u,v);
        E3c -= 0.5*e;
    }
    ledger.E3c = E3c;
}

// ---- reference energy (cert 2.7): bare (non-antisymmetrized) active integrals (H3) ----

double dsrg_sf_tensors::compute_eref(const double* h_core_diag, const double* h_active,
                                     double e_scalar){
    const int na=n_a, nc=n_c;
    const long aux = R->aux_n_ao;
    auto id4 = [na](int a,int b,int c,int d){ return (((size_t)a*na+b)*na+c)*na+d; };

    dump_h_core.assign(h_core_diag, h_core_diag + nc);
    dump_h_act.assign(h_active, h_active + (size_t)na*na);
    dump_e_scalar = e_scalar;

    double E0 = e_scalar;
    for(int m=0;m<nc;m++) E0 += h_core_diag[m] + e_c[m];          // core: (h + f) diagonal

    double e1=0.0;                                               // 0.5 (h + f)^u_v L1^v_u
    for(int u=0;u<na;u++) for(int v=0;v<na;v++){
        const double F_uv = (u==v ? e_a[u] : 0.0);              // semicanonical active Fock diag
        e1 += (h_active[u*na+v] + F_uv)*L1[v*na+u];
    }
    E0 += 0.5*e1;

    // 0.5 g^{uv}_{xy} SF_L2^{xy}_{uv} with BARE g = <uv|xy> = (ux|vy) = sum_k AA(u,x) AA(v,y).
    std::vector<double> graw((size_t)na*na*na*na, 0.0);          // [(u,x)][(v,y)]
    contract_aux(R->AA_RI_M, na*na, R->AA_RI_M, na*na, aux, graw.data());
    double e2=0.0;
    for(int u=0;u<na;u++) for(int v=0;v<na;v++)
    for(int x=0;x<na;x++) for(int y=0;y<na;y++){
        const double g = graw[(size_t)(u*na+x)*(na*na)+(v*na+y)];
        e2 += g*SF_L2[id4(x,y,u,v)];
    }
    E0 += 0.5*e2;

    ledger.Eref = E0;
    return E0;
}

// ---- dressed active operator Hbar = H_seed + 1/2[H~1,A]_{aa,aaaa} (Forte H_A_Ca_small).
//      H2 = renormalized tei (Vt blocks), H1 = renormalized Fock (Ft cross blocks); the
//      commutator temp is accumulated then Hermitized (C[uv]+=t; C[vu]+=t, alpha=1/2).
//      The two batched CAVV/CCAV Hbar1 pieces are left for the driver. Fat legs GEMM,
//      active density folds/closings are scalar -- same idiom as compute_e2. ----

void dsrg_sf_tensors::build_hbar(){
    const int na=n_a, nc=n_c, nv=n_v;
    const size_t na2=(size_t)na*na, na3=na2*na, na4=na2*na2;
    const long aux = R->aux_n_ao;
    auto id4=[na](int a,int b,int c,int d){ return (((size_t)a*na+b)*na+c)*na+d; };

    // C[na2 x na2] = A[na2 x K] B[K x na2]
    auto gemm_aa=[&](const std::vector<double>&A,const std::vector<double>&B,int K,std::vector<double>&C){
        C.assign(na4,0.0);
        if(K<=0) return;
        nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,(int)na2,(int)na2,K,
                       1.0,A.data(),K,B.data(),(int)na2,0.0,C.data(),(int)na2);
    };
    // M[dA x dB] = sum_fat gA[fat x dA] gB[fat x dB] (fat leading); M must be pre-sized.
    auto cx=[&](const std::vector<double>&gA,const std::vector<double>&gB,int nfat,int dA,int dB,std::vector<double>&M){
        for(size_t i=0;i<(size_t)dA*dB;i++) M[i]=0.0;
        if(nfat<=0) return;
        contract_ext(gA.data(),gB.data(),nfat,dA,dB,M.data());
    };

    // renormalized H2 block accessors (physicist <p1 p2|p3 p4>)
    auto Vvvaa=[&](int e,int f,int x,int y){ return Vt_vvaa[(((size_t)e*nv+f)*na+x)*na+y]; };
    auto Vaacc=[&](int w,int x,int i,int j){ return Vt_aacc[(((size_t)w*na+x)*nc+i)*nc+j]; };
    auto Vavca=[&](int v,int e,int m,int x){ return Vt_avca[(((size_t)v*nv+e)*nc+m)*na+x]; };
    auto Vvaaa=[&](int e,int v,int w,int x){ return Vt_vaaa[(((size_t)e*na+v)*na+w)*na+x]; };
    auto Vaaca=[&](int u,int v,int m,int y){ return Vt_aaca[(((size_t)u*na+v)*nc+m)*na+y]; };
    auto Vavaa=[&](int z,int a,int x,int y){ return Vt_avaa[(((size_t)z*nv+a)*na+x)*na+y]; };

    // avac <ue|xm> = (ux|em); not stored by the S2 pipeline, so built here like compute_e2.
    std::vector<double> Vt_avac((size_t)na*nv*na*nc, 0.0);
    {
        std::vector<double> raw((size_t)na*na*nv*nc, 0.0);
        contract_aux(R->AA_RI_M, na*na, R->VC_RI_M, nv*nc, aux, raw.data());
        for(int u=0;u<na;u++) for(int e=0;e<nv;e++)
        for(int x=0;x<na;x++) for(int m=0;m<nc;m++){
            const double D = e_a[u]+e_v[e]-e_a[x]-e_c[m];
            Vt_avac[(((size_t)u*nv+e)*na+x)*nc+m] =
                raw[(size_t)(u*na+x)*(nv*nc)+(e*nc+m)]*(1.0+rfac(D));
        }
    }
    auto Vavac=[&](int u,int e,int x,int m){ return Vt_avac[(((size_t)u*nv+e)*na+x)*nc+m]; };

    // amplitude accessors
    auto Taavv=[&](int u,int v,int e,int f){ return T2_aavv[(((size_t)u*na+v)*nv+e)*nv+f]; };
    auto Tccaa=[&](int i,int j,int u,int v){ return T2_ccaa[(((size_t)i*nc+j)*na+u)*na+v]; };
    auto Tcaav=[&](int i,int v,int u,int a){ return T2_caav[(((size_t)i*na+v)*na+u)*nv+a]; };
    auto Tacav=[&](int v,int i,int u,int a){ return T2_acav[(((size_t)v*nc+i)*na+u)*nv+a]; };
    auto Taava=[&](int u,int v,int e,int y){ return T2_aava[(((size_t)u*na+v)*nv+e)*na+y]; };
    auto Tcaaa=[&](int i,int w,int u,int v){ return T2_caaa[(((size_t)i*na+w)*na+u)*na+v]; };

    // S2 = 2J - K (Forte sa_mrpt2.cc:394-397): ijab->ijba is the particle swap, except
    // caav/acav/aava whose swapped partner is the other stored block.
    std::vector<double> S2_aavv(T2_aavv.size()), S2_ccaa(T2_ccaa.size()), S2_caaa(T2_caaa.size());
    std::vector<double> S2_caav(T2_caav.size()), S2_acav(T2_acav.size()), S2_aava(T2_aava.size());
    for(int u=0;u<na;u++) for(int v=0;v<na;v++) for(int e=0;e<nv;e++) for(int f=0;f<nv;f++)
        S2_aavv[(((size_t)u*na+v)*nv+e)*nv+f]=2.0*Taavv(u,v,e,f)-Taavv(u,v,f,e);
    for(int i=0;i<nc;i++) for(int j=0;j<nc;j++) for(int u=0;u<na;u++) for(int v=0;v<na;v++)
        S2_ccaa[(((size_t)i*nc+j)*na+u)*na+v]=2.0*Tccaa(i,j,u,v)-Tccaa(i,j,v,u);
    for(int i=0;i<nc;i++) for(int w=0;w<na;w++) for(int u=0;u<na;u++) for(int v=0;v<na;v++)
        S2_caaa[(((size_t)i*na+w)*na+u)*na+v]=2.0*Tcaaa(i,w,u,v)-Tcaaa(i,w,v,u);
    for(int i=0;i<nc;i++) for(int v=0;v<na;v++) for(int u=0;u<na;u++) for(int a=0;a<nv;a++)
        S2_caav[(((size_t)i*na+v)*na+u)*nv+a]=2.0*Tcaav(i,v,u,a)-Tacav(v,i,u,a);
    for(int v=0;v<na;v++) for(int i=0;i<nc;i++) for(int u=0;u<na;u++) for(int a=0;a<nv;a++)
        S2_acav[(((size_t)v*nc+i)*na+u)*nv+a]=2.0*Tacav(v,i,u,a)-Tcaav(i,v,u,a);
    for(int u=0;u<na;u++) for(int v=0;v<na;v++) for(int e=0;e<nv;e++) for(int y=0;y<na;y++)
        S2_aava[(((size_t)u*na+v)*nv+e)*na+y]=2.0*Taava(u,v,e,y)-Taava(v,u,e,y);
    auto Saavv=[&](int u,int v,int e,int f){ return S2_aavv[(((size_t)u*na+v)*nv+e)*nv+f]; };
    auto Sccaa=[&](int i,int j,int u,int v){ return S2_ccaa[(((size_t)i*nc+j)*na+u)*na+v]; };
    auto Scaav=[&](int i,int v,int u,int a){ return S2_caav[(((size_t)i*na+v)*na+u)*nv+a]; };
    auto Sacav=[&](int v,int i,int u,int a){ return S2_acav[(((size_t)v*nc+i)*na+u)*nv+a]; };
    auto Saava=[&](int u,int v,int e,int y){ return S2_aava[(((size_t)u*na+v)*nv+e)*na+y]; };
    auto Scaaa=[&](int i,int w,int u,int v){ return S2_caaa[(((size_t)i*na+w)*na+u)*na+v]; };

    // G2 = 2*H2 - H2(swap ket) on {avac, avaa, aaac}
    std::vector<double> G2_avac((size_t)na*nv*na*nc,0.0), G2_avaa((size_t)na*nv*na*na,0.0),
                        G2_aaac((size_t)na*na*na*nc,0.0);
    for(int u=0;u<na;u++) for(int e=0;e<nv;e++) for(int v=0;v<na;v++) for(int m=0;m<nc;m++)
        G2_avac[(((size_t)u*nv+e)*na+v)*nc+m]=2.0*Vavac(u,e,v,m)-Vavca(u,e,m,v);
    for(int u=0;u<na;u++) for(int e=0;e<nv;e++) for(int x=0;x<na;x++) for(int y=0;y<na;y++)
        G2_avaa[(((size_t)u*nv+e)*na+x)*na+y]=2.0*Vavaa(u,e,x,y)-Vavaa(u,e,y,x);   // 2<ue|xy>-<ue|yx>
    for(int u=0;u<na;u++) for(int v=0;v<na;v++) for(int w=0;w<na;w++) for(int m=0;m<nc;m++)
        G2_aaac[(((size_t)u*na+v)*na+w)*nc+m]=2.0*Vaaca(v,u,m,w)-Vaaca(u,v,m,w);
    auto Gavac=[&](int u,int e,int v,int m){ return G2_avac[(((size_t)u*nv+e)*na+v)*nc+m]; };
    auto Gavaa=[&](int u,int e,int x,int y){ return G2_avaa[(((size_t)u*nv+e)*na+x)*na+y]; };
    auto Gaaac=[&](int u,int v,int w,int m){ return G2_aaac[(((size_t)u*na+v)*na+w)*nc+m]; };

    auto l1=[&](int p,int q){ return L1[p*na+q]; };
    auto e1=[&](int p,int q){ return Eta1[p*na+q]; };
    auto l2=[&](int x,int y,int u,int v){ return SF_L2[id4(x,y,u,v)]; };
    auto Fav=[&](int t,int a){ return Ft_av[t*nv+a]; };
    auto Fca=[&](int i,int t){ return Ft_ca[i*na+t]; };
    auto Fcv=[&](int i,int a){ return Ft_cv[i*nv+a]; };
    auto T1av=[&](int t,int a){ return T1_av[t*nv+a]; };
    auto T1ca=[&](int i,int t){ return T1_ca[i*na+t]; };
    auto T1cv=[&](int i,int a){ return T1_cv[i*nv+a]; };

    std::vector<double> temp1(na2,0.0);      // C1 accumulator
    std::vector<double> M(na4,0.0), gA, gB;  // reusable GEMM output / gathers

    // ============================ C1 : H1_T_C1a_smallS ============================
    // C1[uv] += H1[ev]T1[ue] - H1[um]T1[mv] + H1[em]S2[umve] + H1[xm]S2[muxv]
    //         + 0.5 H1[ex]S2[yuev]L1[xy] - 0.5 H1[ym]S2[muxv]L1[xy]  (aaaa amps drop out)
    for(int u=0;u<na;u++) for(int v=0;v<na;v++){
        double t=0.0;
        for(int e=0;e<nv;e++) t += Fav(v,e)*T1av(u,e);
        for(int m=0;m<nc;m++) t -= Fca(m,u)*T1ca(m,v);
        for(int e=0;e<nv;e++) for(int m=0;m<nc;m++) t += Fcv(m,e)*Sacav(u,m,v,e);
        for(int x=0;x<na;x++) for(int m=0;m<nc;m++) t += Fca(m,x)*Scaaa(m,u,x,v);
        for(int e=0;e<nv;e++) for(int x=0;x<na;x++) for(int y=0;y<na;y++)
            t += 0.5*Fav(x,e)*Saava(y,u,e,v)*l1(x,y);
        for(int x=0;x<na;x++) for(int y=0;y<na;y++) for(int m=0;m<nc;m++)
            t -= 0.5*Fca(m,y)*Scaaa(m,u,x,v)*l1(x,y);
        temp1[u*na+v]+=t;
    }

    // ============================ C1 : H2_T_C1a_smallG ============================
    // C1[uv] += T1[ma]G2[uavm] + 0.5 T1[xe]L1[yx]G2[uevy] - 0.5 T1[mx]L1[xy]G2[uyvm]
    // The first term's a is a PARTICLE composite: both the virtual and active branches.
    {   // fold T1 with L1 into virtual/core-indexed intermediates, then contract.
        std::vector<double> TLe((size_t)nv*na,0.0), TLm((size_t)nc*na,0.0);   // [e][y],[m][y]
        for(int e=0;e<nv;e++) for(int y=0;y<na;y++){ double a=0; for(int x=0;x<na;x++) a+=T1av(x,e)*l1(y,x); TLe[(size_t)e*na+y]=a; }
        for(int m=0;m<nc;m++) for(int y=0;y<na;y++){ double a=0; for(int x=0;x<na;x++) a+=T1ca(m,x)*l1(x,y); TLm[(size_t)m*na+y]=a; }
        for(int u=0;u<na;u++) for(int v=0;v<na;v++){
            double t=0.0;
            for(int a=0;a<nv;a++) for(int m=0;m<nc;m++) t += T1cv(m,a)*Gavac(u,a,v,m);
            for(int a=0;a<na;a++) for(int m=0;m<nc;m++) t += T1ca(m,a)*Gaaac(u,a,v,m);
            for(int e=0;e<nv;e++) for(int y=0;y<na;y++) t += 0.5*TLe[(size_t)e*na+y]*Gavaa(u,e,v,y);
            for(int m=0;m<nc;m++) for(int y=0;y<na;y++) t -= 0.5*TLm[(size_t)m*na+y]*Gaaac(u,y,v,m);
            temp1[u*na+v]+=t;
        }
    }
    // C1[wz] += 0.5 G2[wezx]T2[uvey]L2[xyuv] - 0.5 G2[wuzm]T2[mvxy]L2[xyuv]
    {   // contract the fat leg (e or m) into an n_a^6 tensor, then close with L2.
        std::vector<double> P((size_t)na3*na3);
        auto Pat=[&](int a1,int a2,int a3,int b1,int b2,int b3){
            return (((size_t)a1*na+a2)*na+a3)*na3+((size_t)b1*na+b2)*na+b3; };
        gA.assign((size_t)nv*na3,0.0); gB.assign((size_t)nv*na3,0.0);   // [e][(w,z,x)],[e][(u,v,y)]
        for(int e=0;e<nv;e++) for(int w=0;w<na;w++) for(int z=0;z<na;z++) for(int x=0;x<na;x++)
            gA[(size_t)e*na3+((w*na+z)*na+x)]=Gavaa(w,e,z,x);
        for(int e=0;e<nv;e++) for(int u=0;u<na;u++) for(int v=0;v<na;v++) for(int y=0;y<na;y++)
            gB[(size_t)e*na3+((u*na+v)*na+y)]=Taava(u,v,e,y);
        if(nv>0) contract_ext(gA.data(),gB.data(),nv,(int)na3,(int)na3,P.data());
        else for(size_t i=0;i<(size_t)na3*na3;i++) P[i]=0.0;
        for(int w=0;w<na;w++) for(int z=0;z<na;z++){
            double t=0.0;
            for(int x=0;x<na;x++) for(int u=0;u<na;u++) for(int v=0;v<na;v++) for(int y=0;y<na;y++)
                t += P[Pat(w,z,x,u,v,y)]*l2(x,y,u,v);
            temp1[w*na+z]+=0.5*t;
        }
        gA.assign((size_t)nc*na3,0.0); gB.assign((size_t)nc*na3,0.0);   // [m][(w,z,u)],[m][(v,x,y)]
        for(int m=0;m<nc;m++) for(int w=0;w<na;w++) for(int z=0;z<na;z++) for(int u=0;u<na;u++)
            gA[(size_t)m*na3+((w*na+z)*na+u)]=Gaaac(w,u,z,m);
        for(int m=0;m<nc;m++) for(int v=0;v<na;v++) for(int x=0;x<na;x++) for(int y=0;y<na;y++)
            gB[(size_t)m*na3+((v*na+x)*na+y)]=Tcaaa(m,v,x,y);
        if(nc>0) contract_ext(gA.data(),gB.data(),nc,(int)na3,(int)na3,P.data());
        else for(size_t i=0;i<(size_t)na3*na3;i++) P[i]=0.0;
        for(int w=0;w<na;w++) for(int z=0;z<na;z++){
            double t=0.0;
            for(int u=0;u<na;u++) for(int v=0;v<na;v++) for(int x=0;x<na;x++) for(int y=0;y<na;y++)
                t += P[Pat(w,z,u,v,x,y)]*l2(x,y,u,v);
            temp1[w*na+z]-=0.5*t;
        }
    }

    // ===================== C1 : H2_T_C1a_smallS direct terms =====================
    // C1[wz] += H2[uemz]S2[mwue] + H2[uezm]S2[wmue] + H2[vumz]S2[mwvu]
    //         - H2[wemu]S2[muze] - H2[weum]S2[umze] - H2[ewvu]S2[vuez]
    {
        const int nuem=na*nv*nc;
        gA.assign((size_t)nuem*na,0.0); gB.assign((size_t)nuem*na,0.0);   // fat (u,e,m)
        for(int u=0;u<na;u++) for(int e=0;e<nv;e++) for(int m=0;m<nc;m++){
            const size_t f=((size_t)(u*nv+e)*nc+m)*na;
            for(int z=0;z<na;z++) gA[f+z]=Vavca(u,e,m,z);
            for(int w=0;w<na;w++) gB[f+w]=Scaav(m,w,u,e);
        }
        cx(gA,gB,nuem,na,na,M);
        for(int w=0;w<na;w++) for(int z=0;z<na;z++) temp1[w*na+z]+=M[(size_t)z*na+w];   // + H2[uemz]S2[mwue]
        for(int u=0;u<na;u++) for(int e=0;e<nv;e++) for(int m=0;m<nc;m++){
            const size_t f=((size_t)(u*nv+e)*nc+m)*na;
            for(int z=0;z<na;z++) gA[f+z]=Vavac(u,e,z,m);
            for(int w=0;w<na;w++) gB[f+w]=Sacav(w,m,u,e);
        }
        cx(gA,gB,nuem,na,na,M);
        for(int w=0;w<na;w++) for(int z=0;z<na;z++) temp1[w*na+z]+=M[(size_t)z*na+w];   // + H2[uezm]S2[wmue]
        for(int e=0;e<nv;e++) for(int m=0;m<nc;m++) for(int u=0;u<na;u++){
            const size_t f=((size_t)(e*nc+m)*na+u)*na;                                  // fat (e,m,u)
            for(int w=0;w<na;w++) gA[f+w]=Vavca(w,e,m,u);
            for(int z=0;z<na;z++) gB[f+z]=Scaav(m,u,z,e);
        }
        cx(gA,gB,nuem,na,na,M);
        for(int w=0;w<na;w++) for(int z=0;z<na;z++) temp1[w*na+z]-=M[(size_t)w*na+z];   // - H2[wemu]S2[muze]
        for(int e=0;e<nv;e++) for(int m=0;m<nc;m++) for(int u=0;u<na;u++){
            const size_t f=((size_t)(e*nc+m)*na+u)*na;
            for(int w=0;w<na;w++) gA[f+w]=Vavac(w,e,u,m);
            for(int z=0;z<na;z++) gB[f+z]=Sacav(u,m,z,e);
        }
        cx(gA,gB,nuem,na,na,M);
        for(int w=0;w<na;w++) for(int z=0;z<na;z++) temp1[w*na+z]-=M[(size_t)w*na+z];   // - H2[weum]S2[umze]
    }
    {   // H2[vumz]S2[mwvu]  fat (v,u,m);  H2[ewvu]S2[vuez]  fat (e,v,u)
        const int nvum=na*na*nc;
        gA.assign((size_t)nvum*na,0.0); gB.assign((size_t)nvum*na,0.0);
        for(int v=0;v<na;v++) for(int u=0;u<na;u++) for(int m=0;m<nc;m++){
            const size_t f=((size_t)(v*na+u)*nc+m)*na;
            for(int z=0;z<na;z++) gA[f+z]=Vaaca(v,u,m,z);
            for(int w=0;w<na;w++) gB[f+w]=Scaaa(m,w,v,u);
        }
        cx(gA,gB,nvum,na,na,M);
        for(int w=0;w<na;w++) for(int z=0;z<na;z++) temp1[w*na+z]+=M[(size_t)z*na+w];   // + H2[vumz]S2[mwvu]
        const int nevu=nv*na*na;
        gA.assign((size_t)nevu*na,0.0); gB.assign((size_t)nevu*na,0.0);
        for(int e=0;e<nv;e++) for(int v=0;v<na;v++) for(int u=0;u<na;u++){
            const size_t f=((size_t)(e*na+v)*na+u)*na;
            for(int w=0;w<na;w++) gA[f+w]=Vvaaa(e,w,v,u);
            for(int z=0;z<na;z++) gB[f+z]=Saava(v,u,e,z);
        }
        cx(gA,gB,nevu,na,na,M);
        for(int w=0;w<na;w++) for(int z=0;z<na;z++) temp1[w*na+z]-=M[(size_t)w*na+z];   // - H2[ewvu]S2[vuez]
    }

    // ============= C1 : H2_T_C1a_smallS temp[wzuv] closed with L1[uv] =============
    {
        std::vector<double> twz(na4,0.0);
        auto addM_wuvz=[&](double c){ for(int w=0;w<na;w++) for(int u=0;u<na;u++) for(int v=0;v<na;v++) for(int z=0;z<na;z++)
            twz[id4(w,z,u,v)]+=c*M[(size_t)(w*na+u)*na2+(v*na+z)]; };
        // +0.5 S2[wvef]H2[efzu]  (fat ef): [(w,v)][(z,u)]
        gA.assign((size_t)nv*nv*na2,0.0); gB.assign((size_t)nv*nv*na2,0.0);
        for(int e=0;e<nv;e++) for(int f=0;f<nv;f++) for(int w=0;w<na;w++) for(int v=0;v<na;v++)
            gA[(size_t)(e*nv+f)*na2+(w*na+v)]=Saavv(w,v,e,f);
        for(int e=0;e<nv;e++) for(int f=0;f<nv;f++) for(int z=0;z<na;z++) for(int u=0;u<na;u++)
            gB[(size_t)(e*nv+f)*na2+(z*na+u)]=Vvvaa(e,f,z,u);
        cx(gA,gB,nv*nv,(int)na2,(int)na2,M);
        for(int w=0;w<na;w++) for(int v=0;v<na;v++) for(int z=0;z<na;z++) for(int u=0;u<na;u++)
            twz[id4(w,z,u,v)]+=0.5*M[(size_t)(w*na+v)*na2+(z*na+u)];
        // +0.5 S2[wvex]H2[exzu]  (fat e,x): [(w,v)][(z,u)]
        gA.assign((size_t)nv*na*na2,0.0); gB.assign((size_t)nv*na*na2,0.0);
        for(int e=0;e<nv;e++) for(int x=0;x<na;x++) for(int w=0;w<na;w++) for(int v=0;v<na;v++)
            gA[(size_t)(e*na+x)*na2+(w*na+v)]=Saava(w,v,e,x);
        for(int e=0;e<nv;e++) for(int x=0;x<na;x++) for(int z=0;z<na;z++) for(int u=0;u<na;u++)
            gB[(size_t)(e*na+x)*na2+(z*na+u)]=Vvaaa(e,x,z,u);
        cx(gA,gB,nv*na,(int)na2,(int)na2,M);
        for(int w=0;w<na;w++) for(int v=0;v<na;v++) for(int z=0;z<na;z++) for(int u=0;u<na;u++)
            twz[id4(w,z,u,v)]+=0.5*M[(size_t)(w*na+v)*na2+(z*na+u)];
        // +0.5 S2[vwex]H2[exuz]  (fat e,x): [(v,w)][(u,z)]
        for(int e=0;e<nv;e++) for(int x=0;x<na;x++) for(int v=0;v<na;v++) for(int w=0;w<na;w++)
            gA[(size_t)(e*na+x)*na2+(v*na+w)]=Saava(v,w,e,x);
        for(int e=0;e<nv;e++) for(int x=0;x<na;x++) for(int u=0;u<na;u++) for(int z=0;z<na;z++)
            gB[(size_t)(e*na+x)*na2+(u*na+z)]=Vvaaa(e,x,u,z);
        cx(gA,gB,nv*na,(int)na2,(int)na2,M);
        for(int v=0;v<na;v++) for(int w=0;w<na;w++) for(int u=0;u<na;u++) for(int z=0;z<na;z++)
            twz[id4(w,z,u,v)]+=0.5*M[(size_t)(v*na+w)*na2+(u*na+z)];
        // -0.5 S2[wmue]H2[vezm]  (fat e,m): [(w,u)][(v,z)]
        gA.assign((size_t)nv*nc*na2,0.0); gB.assign((size_t)nv*nc*na2,0.0);
        for(int e=0;e<nv;e++) for(int m=0;m<nc;m++) for(int w=0;w<na;w++) for(int u=0;u<na;u++)
            gA[(size_t)(e*nc+m)*na2+(w*na+u)]=Sacav(w,m,u,e);
        for(int e=0;e<nv;e++) for(int m=0;m<nc;m++) for(int v=0;v<na;v++) for(int z=0;z<na;z++)
            gB[(size_t)(e*nc+m)*na2+(v*na+z)]=Vavac(v,e,z,m);
        cx(gA,gB,nv*nc,(int)na2,(int)na2,M); addM_wuvz(-0.5);
        // -0.5 S2[mwue]H2[vemz]  (fat m,e): [(w,u)][(v,z)]
        for(int m=0;m<nc;m++) for(int e=0;e<nv;e++) for(int w=0;w<na;w++) for(int u=0;u<na;u++)
            gA[(size_t)(m*nv+e)*na2+(w*na+u)]=Scaav(m,w,u,e);
        for(int m=0;m<nc;m++) for(int e=0;e<nv;e++) for(int v=0;v<na;v++) for(int z=0;z<na;z++)
            gB[(size_t)(m*nv+e)*na2+(v*na+z)]=Vavca(v,e,m,z);
        cx(gA,gB,nc*nv,(int)na2,(int)na2,M); addM_wuvz(-0.5);
        // -0.5 S2[mwxu]H2[xvmz]  (fat m,x): [(w,u)][(v,z)]
        gA.assign((size_t)nc*na*na2,0.0); gB.assign((size_t)nc*na*na2,0.0);
        for(int m=0;m<nc;m++) for(int x=0;x<na;x++) for(int w=0;w<na;w++) for(int u=0;u<na;u++)
            gA[(size_t)(m*na+x)*na2+(w*na+u)]=Scaaa(m,w,x,u);
        for(int m=0;m<nc;m++) for(int x=0;x<na;x++) for(int v=0;v<na;v++) for(int z=0;z<na;z++)
            gB[(size_t)(m*na+x)*na2+(v*na+z)]=Vaaca(x,v,m,z);
        cx(gA,gB,nc*na,(int)na2,(int)na2,M); addM_wuvz(-0.5);
        // -0.5 S2[mwux]H2[vxmz]  (fat m,x): [(w,u)][(v,z)]
        for(int m=0;m<nc;m++) for(int x=0;x<na;x++) for(int w=0;w<na;w++) for(int u=0;u<na;u++)
            gA[(size_t)(m*na+x)*na2+(w*na+u)]=Scaaa(m,w,u,x);
        for(int m=0;m<nc;m++) for(int x=0;x<na;x++) for(int v=0;v<na;v++) for(int z=0;z<na;z++)
            gB[(size_t)(m*na+x)*na2+(v*na+z)]=Vaaca(v,x,m,z);
        cx(gA,gB,nc*na,(int)na2,(int)na2,M); addM_wuvz(-0.5);
        // +0.25 (S2[jwxu]L1[xy])H2[yvjz]  (fat j,y): [(w,u)][(v,z)]
        {
            std::vector<double> SL((size_t)nc*na3,0.0);   // [j][w][y][u] = sum_x S2[jwxu]L1[xy]
            for(int j=0;j<nc;j++) for(int w=0;w<na;w++) for(int y=0;y<na;y++) for(int u=0;u<na;u++){
                double a=0; for(int x=0;x<na;x++) a+=Scaaa(j,w,x,u)*l1(x,y);
                SL[(size_t)j*na3+((w*na+y)*na+u)]=a;
            }
            gA.assign((size_t)nc*na*na2,0.0); gB.assign((size_t)nc*na*na2,0.0);
            for(int j=0;j<nc;j++) for(int y=0;y<na;y++) for(int w=0;w<na;w++) for(int u=0;u<na;u++)
                gA[(size_t)(j*na+y)*na2+(w*na+u)]=SL[(size_t)j*na3+((w*na+y)*na+u)];
            for(int j=0;j<nc;j++) for(int y=0;y<na;y++) for(int v=0;v<na;v++) for(int z=0;z<na;z++)
                gB[(size_t)(j*na+y)*na2+(v*na+z)]=Vaaca(y,v,j,z);
            cx(gA,gB,nc*na,(int)na2,(int)na2,M); addM_wuvz(0.25);
        }
        // -0.25 (S2[ywbu]L1[xy])H2[bvxz]  and  -0.25 (S2[wybu]L1[xy])H2[bvzx]  (fat b,x)
        {
            std::vector<double> SF1((size_t)nv*na3,0.0), SF2((size_t)nv*na3,0.0);  // [b][x][w][u]
            for(int b=0;b<nv;b++) for(int x=0;x<na;x++) for(int w=0;w<na;w++) for(int u=0;u<na;u++){
                double a1=0,a2=0; for(int y=0;y<na;y++){ a1+=Saava(y,w,b,u)*l1(x,y); a2+=Saava(w,y,b,u)*l1(x,y); }
                SF1[(size_t)b*na3+((x*na+w)*na+u)]=a1;
                SF2[(size_t)b*na3+((x*na+w)*na+u)]=a2;
            }
            gA.assign((size_t)nv*na*na2,0.0); gB.assign((size_t)nv*na*na2,0.0);
            for(int b=0;b<nv;b++) for(int x=0;x<na;x++) for(int w=0;w<na;w++) for(int u=0;u<na;u++)
                gA[(size_t)(b*na+x)*na2+(w*na+u)]=SF1[(size_t)b*na3+((x*na+w)*na+u)];
            for(int b=0;b<nv;b++) for(int x=0;x<na;x++) for(int v=0;v<na;v++) for(int z=0;z<na;z++)
                gB[(size_t)(b*na+x)*na2+(v*na+z)]=Vvaaa(b,v,x,z);
            cx(gA,gB,nv*na,(int)na2,(int)na2,M); addM_wuvz(-0.25);
            for(int b=0;b<nv;b++) for(int x=0;x<na;x++) for(int w=0;w<na;w++) for(int u=0;u<na;u++)
                gA[(size_t)(b*na+x)*na2+(w*na+u)]=SF2[(size_t)b*na3+((x*na+w)*na+u)];
            for(int b=0;b<nv;b++) for(int x=0;x<na;x++) for(int v=0;v<na;v++) for(int z=0;z<na;z++)
                gB[(size_t)(b*na+x)*na2+(v*na+z)]=Vvaaa(b,v,z,x);
            cx(gA,gB,nv*na,(int)na2,(int)na2,M); addM_wuvz(-0.25);
        }
        for(int w=0;w<na;w++) for(int z=0;z<na;z++){
            double t=0.0; for(int u=0;u<na;u++) for(int v=0;v<na;v++) t+=twz[id4(w,z,u,v)]*l1(u,v);
            temp1[w*na+z]+=t;
        }
    }

    // ============= C1 : H2_T_C1a_smallS temp[wzuv] closed with Eta1[uv] =============
    {
        std::vector<double> twz(na4,0.0);
        auto addM_vzwu=[&](double c){ for(int v=0;v<na;v++) for(int z=0;z<na;z++) for(int w=0;w<na;w++) for(int u=0;u<na;u++)
            twz[id4(w,z,u,v)]+=c*M[(size_t)(v*na+z)*na2+(w*na+u)]; };
        // -0.5 S2[mnzu]H2[wvmn]  (fat m,n): [(z,u)][(w,v)]
        gA.assign((size_t)nc*nc*na2,0.0); gB.assign((size_t)nc*nc*na2,0.0);
        for(int m=0;m<nc;m++) for(int n=0;n<nc;n++) for(int z=0;z<na;z++) for(int u=0;u<na;u++)
            gA[(size_t)(m*nc+n)*na2+(z*na+u)]=Sccaa(m,n,z,u);
        for(int m=0;m<nc;m++) for(int n=0;n<nc;n++) for(int w=0;w<na;w++) for(int v=0;v<na;v++)
            gB[(size_t)(m*nc+n)*na2+(w*na+v)]=Vaacc(w,v,m,n);
        cx(gA,gB,nc*nc,(int)na2,(int)na2,M);
        for(int z=0;z<na;z++) for(int u=0;u<na;u++) for(int w=0;w<na;w++) for(int v=0;v<na;v++)
            twz[id4(w,z,u,v)]-=0.5*M[(size_t)(z*na+u)*na2+(w*na+v)];
        // -0.5 S2[mxzu]H2[wvmx]  (fat m,x): [(z,u)][(w,v)]
        gA.assign((size_t)nc*na*na2,0.0); gB.assign((size_t)nc*na*na2,0.0);
        for(int m=0;m<nc;m++) for(int x=0;x<na;x++) for(int z=0;z<na;z++) for(int u=0;u<na;u++)
            gA[(size_t)(m*na+x)*na2+(z*na+u)]=Scaaa(m,x,z,u);
        for(int m=0;m<nc;m++) for(int x=0;x<na;x++) for(int w=0;w<na;w++) for(int v=0;v<na;v++)
            gB[(size_t)(m*na+x)*na2+(w*na+v)]=Vaaca(w,v,m,x);
        cx(gA,gB,nc*na,(int)na2,(int)na2,M);
        for(int z=0;z<na;z++) for(int u=0;u<na;u++) for(int w=0;w<na;w++) for(int v=0;v<na;v++)
            twz[id4(w,z,u,v)]-=0.5*M[(size_t)(z*na+u)*na2+(w*na+v)];
        // -0.5 S2[mxuz]H2[vwmx]  (fat m,x): [(u,z)][(v,w)]
        for(int m=0;m<nc;m++) for(int x=0;x<na;x++) for(int u=0;u<na;u++) for(int z=0;z<na;z++)
            gA[(size_t)(m*na+x)*na2+(u*na+z)]=Scaaa(m,x,u,z);
        for(int m=0;m<nc;m++) for(int x=0;x<na;x++) for(int v=0;v<na;v++) for(int w=0;w<na;w++)
            gB[(size_t)(m*na+x)*na2+(v*na+w)]=Vaaca(v,w,m,x);
        cx(gA,gB,nc*na,(int)na2,(int)na2,M);
        for(int u=0;u<na;u++) for(int z=0;z<na;z++) for(int v=0;v<na;v++) for(int w=0;w<na;w++)
            twz[id4(w,z,u,v)]-=0.5*M[(size_t)(u*na+z)*na2+(v*na+w)];
        // +0.5 S2[vmze]H2[weum]  (fat m,e): [(v,z)][(w,u)]
        gA.assign((size_t)nc*nv*na2,0.0); gB.assign((size_t)nc*nv*na2,0.0);
        for(int m=0;m<nc;m++) for(int e=0;e<nv;e++) for(int v=0;v<na;v++) for(int z=0;z<na;z++)
            gA[(size_t)(m*nv+e)*na2+(v*na+z)]=Sacav(v,m,z,e);
        for(int m=0;m<nc;m++) for(int e=0;e<nv;e++) for(int w=0;w<na;w++) for(int u=0;u<na;u++)
            gB[(size_t)(m*nv+e)*na2+(w*na+u)]=Vavac(w,e,u,m);
        cx(gA,gB,nc*nv,(int)na2,(int)na2,M); addM_vzwu(0.5);
        // +0.5 S2[mvze]H2[wemu]  (fat m,e): [(v,z)][(w,u)]
        for(int m=0;m<nc;m++) for(int e=0;e<nv;e++) for(int v=0;v<na;v++) for(int z=0;z<na;z++)
            gA[(size_t)(m*nv+e)*na2+(v*na+z)]=Scaav(m,v,z,e);
        for(int m=0;m<nc;m++) for(int e=0;e<nv;e++) for(int w=0;w<na;w++) for(int u=0;u<na;u++)
            gB[(size_t)(m*nv+e)*na2+(w*na+u)]=Vavca(w,e,m,u);
        cx(gA,gB,nc*nv,(int)na2,(int)na2,M); addM_vzwu(0.5);
        // +0.5 S2[xvez]H2[ewxu]  (fat e,x): [(v,z)][(w,u)]
        gA.assign((size_t)nv*na*na2,0.0); gB.assign((size_t)nv*na*na2,0.0);
        for(int e=0;e<nv;e++) for(int x=0;x<na;x++) for(int v=0;v<na;v++) for(int z=0;z<na;z++)
            gA[(size_t)(e*na+x)*na2+(v*na+z)]=Saava(x,v,e,z);
        for(int e=0;e<nv;e++) for(int x=0;x<na;x++) for(int w=0;w<na;w++) for(int u=0;u<na;u++)
            gB[(size_t)(e*na+x)*na2+(w*na+u)]=Vvaaa(e,w,x,u);
        cx(gA,gB,nv*na,(int)na2,(int)na2,M); addM_vzwu(0.5);
        // +0.5 S2[vxez]H2[ewux]  (fat e,x): [(v,z)][(w,u)]
        for(int e=0;e<nv;e++) for(int x=0;x<na;x++) for(int v=0;v<na;v++) for(int z=0;z<na;z++)
            gA[(size_t)(e*na+x)*na2+(v*na+z)]=Saava(v,x,e,z);
        for(int e=0;e<nv;e++) for(int x=0;x<na;x++) for(int w=0;w<na;w++) for(int u=0;u<na;u++)
            gB[(size_t)(e*na+x)*na2+(w*na+u)]=Vvaaa(e,w,u,x);
        cx(gA,gB,nv*na,(int)na2,(int)na2,M); addM_vzwu(0.5);
        // -0.25 (S2[yvbz]Eta1[xy])H2[bwxu]  (fat b,x): [(v,z)][(w,u)]
        {
            std::vector<double> SF((size_t)nv*na3,0.0);   // [b][x][v][z] = sum_y S2[yvbz]Eta1[xy]
            for(int b=0;b<nv;b++) for(int x=0;x<na;x++) for(int v=0;v<na;v++) for(int z=0;z<na;z++){
                double a=0; for(int y=0;y<na;y++) a+=Saava(y,v,b,z)*e1(x,y);
                SF[(size_t)b*na3+((x*na+v)*na+z)]=a;
            }
            gA.assign((size_t)nv*na*na2,0.0); gB.assign((size_t)nv*na*na2,0.0);
            for(int b=0;b<nv;b++) for(int x=0;x<na;x++) for(int v=0;v<na;v++) for(int z=0;z<na;z++)
                gA[(size_t)(b*na+x)*na2+(v*na+z)]=SF[(size_t)b*na3+((x*na+v)*na+z)];
            for(int b=0;b<nv;b++) for(int x=0;x<na;x++) for(int w=0;w<na;w++) for(int u=0;u<na;u++)
                gB[(size_t)(b*na+x)*na2+(w*na+u)]=Vvaaa(b,w,x,u);
            cx(gA,gB,nv*na,(int)na2,(int)na2,M); addM_vzwu(-0.25);
        }
        // +0.25 (Eta1[xy]H2[ywju])S2[jvxz]  and  +0.25 (Eta1[xy]H2[wyju])S2[jvzx]  (fat j,x)
        {
            std::vector<double> HF1((size_t)nc*na3,0.0), HF2((size_t)nc*na3,0.0);  // [j][x][w][u]
            for(int j=0;j<nc;j++) for(int x=0;x<na;x++) for(int w=0;w<na;w++) for(int u=0;u<na;u++){
                double a1=0,a2=0; for(int y=0;y<na;y++){ a1+=Vaaca(y,w,j,u)*e1(x,y); a2+=Vaaca(w,y,j,u)*e1(x,y); }
                HF1[(size_t)j*na3+((x*na+w)*na+u)]=a1;
                HF2[(size_t)j*na3+((x*na+w)*na+u)]=a2;
            }
            gA.assign((size_t)nc*na*na2,0.0); gB.assign((size_t)nc*na*na2,0.0);
            for(int j=0;j<nc;j++) for(int x=0;x<na;x++) for(int v=0;v<na;v++) for(int z=0;z<na;z++)
                gA[(size_t)(j*na+x)*na2+(v*na+z)]=Scaaa(j,v,x,z);
            for(int j=0;j<nc;j++) for(int x=0;x<na;x++) for(int w=0;w<na;w++) for(int u=0;u<na;u++)
                gB[(size_t)(j*na+x)*na2+(w*na+u)]=HF1[(size_t)j*na3+((x*na+w)*na+u)];
            cx(gA,gB,nc*na,(int)na2,(int)na2,M); addM_vzwu(0.25);
            for(int j=0;j<nc;j++) for(int x=0;x<na;x++) for(int v=0;v<na;v++) for(int z=0;z<na;z++)
                gA[(size_t)(j*na+x)*na2+(v*na+z)]=Scaaa(j,v,z,x);
            for(int j=0;j<nc;j++) for(int x=0;x<na;x++) for(int w=0;w<na;w++) for(int u=0;u<na;u++)
                gB[(size_t)(j*na+x)*na2+(w*na+u)]=HF2[(size_t)j*na3+((x*na+w)*na+u)];
            cx(gA,gB,nc*na,(int)na2,(int)na2,M); addM_vzwu(0.25);
        }
        for(int w=0;w<na;w++) for(int z=0;z<na;z++){
            double t=0.0; for(int u=0;u<na;u++) for(int v=0;v<na;v++) t+=twz[id4(w,z,u,v)]*e1(u,v);
            temp1[w*na+z]+=t;
        }
    }

    // ================== C1 : H2_T_C1a_smallS L2 ring terms ==================
    // Each contracts one fat leg (j/a/i/b) into an n_a^6 tensor, then closes with L2.
    {
        std::vector<double> P((size_t)na3*na3);
        auto idx6=[&](int a1,int a2,int a3,int b1,int b2,int b3){
            return (((size_t)a1*na+a2)*na+a3)*na3+((size_t)b1*na+b2)*na+b3; };
        auto contractFat=[&](int nfat, auto Afn, auto Bfn){
            gA.assign((size_t)nfat*na3,0.0); gB.assign((size_t)nfat*na3,0.0);
            for(int p=0;p<nfat;p++) for(int a1=0;a1<na;a1++) for(int a2=0;a2<na;a2++) for(int a3=0;a3<na;a3++)
                gA[(size_t)p*na3+((a1*na+a2)*na+a3)]=Afn(p,a1,a2,a3);
            for(int p=0;p<nfat;p++) for(int b1=0;b1<na;b1++) for(int b2=0;b2<na;b2++) for(int b3=0;b3<na;b3++)
                gB[(size_t)p*na3+((b1*na+b2)*na+b3)]=Bfn(p,b1,b2,b3);
            if(nfat>0) contract_ext(gA.data(),gB.data(),nfat,(int)na3,(int)na3,P.data());
            else for(size_t i=0;i<(size_t)na3*na3;i++) P[i]=0.0;
        };
        // r1: +0.5 H2[vujz]T2[jwyx]L2[xyuv]  fat j: A[j](v,u,z) B[j](w,y,x)
        contractFat(nc,[&](int j,int v,int u,int z){return Vaaca(v,u,j,z);},
                       [&](int j,int w,int y,int x){return Tcaaa(j,w,y,x);});
        for(int v=0;v<na;v++) for(int u=0;u<na;u++) for(int z=0;z<na;z++)
        for(int w=0;w<na;w++) for(int y=0;y<na;y++) for(int x=0;x<na;x++)
            temp1[w*na+z]+=0.5*P[idx6(v,u,z,w,y,x)]*l2(x,y,u,v);
        // r2: +0.5 H2[auzx]S2[wvay]L2[xyuv]  fat a: A[a](u,z,x) B[a](w,v,y)
        contractFat(nv,[&](int a,int u,int z,int x){return Vvaaa(a,u,z,x);},
                       [&](int a,int w,int v,int y){return Saava(w,v,a,y);});
        for(int u=0;u<na;u++) for(int z=0;z<na;z++) for(int x=0;x<na;x++)
        for(int w=0;w<na;w++) for(int v=0;v<na;v++) for(int y=0;y<na;y++)
            temp1[w*na+z]+=0.5*P[idx6(u,z,x,w,v,y)]*l2(x,y,u,v);
        // r3: -0.5 H2[auxz]T2[wvay]L2[xyuv]  fat a: A[a](u,x,z) B[a](w,v,y)
        contractFat(nv,[&](int a,int u,int x,int z){return Vvaaa(a,u,x,z);},
                       [&](int a,int w,int v,int y){return Taava(w,v,a,y);});
        for(int u=0;u<na;u++) for(int x=0;x<na;x++) for(int z=0;z<na;z++)
        for(int w=0;w<na;w++) for(int v=0;v<na;v++) for(int y=0;y<na;y++)
            temp1[w*na+z]-=0.5*P[idx6(u,x,z,w,v,y)]*l2(x,y,u,v);
        // r4: -0.5 H2[auxz]T2[vway]L2[xyvu]  fat a: A[a](u,x,z) B[a](v,w,y)
        contractFat(nv,[&](int a,int u,int x,int z){return Vvaaa(a,u,x,z);},
                       [&](int a,int v,int w,int y){return Taava(v,w,a,y);});
        for(int u=0;u<na;u++) for(int x=0;x<na;x++) for(int z=0;z<na;z++)
        for(int v=0;v<na;v++) for(int w=0;w<na;w++) for(int y=0;y<na;y++)
            temp1[w*na+z]-=0.5*P[idx6(u,x,z,v,w,y)]*l2(x,y,v,u);
        // r5: -0.5 H2[bwyx]T2[vubz]L2[xyuv]  fat b: A[b](w,y,x) B[b](v,u,z)
        contractFat(nv,[&](int b,int w,int y,int x){return Vvaaa(b,w,y,x);},
                       [&](int b,int v,int u,int z){return Taava(v,u,b,z);});
        for(int w=0;w<na;w++) for(int y=0;y<na;y++) for(int x=0;x<na;x++)
        for(int v=0;v<na;v++) for(int u=0;u<na;u++) for(int z=0;z<na;z++)
            temp1[w*na+z]-=0.5*P[idx6(w,y,x,v,u,z)]*l2(x,y,u,v);
        // r6: -0.5 H2[wuix]S2[ivzy]L2[xyuv]  fat i: A[i](w,u,x) B[i](v,z,y)
        contractFat(nc,[&](int i,int w,int u,int x){return Vaaca(w,u,i,x);},
                       [&](int i,int v,int z,int y){return Scaaa(i,v,z,y);});
        for(int w=0;w<na;w++) for(int u=0;u<na;u++) for(int x=0;x<na;x++)
        for(int v=0;v<na;v++) for(int z=0;z<na;z++) for(int y=0;y<na;y++)
            temp1[w*na+z]-=0.5*P[idx6(w,u,x,v,z,y)]*l2(x,y,u,v);
        // r7: +0.5 H2[uwix]T2[ivzy]L2[xyuv]  fat i: A[i](u,w,x) B[i](v,z,y)
        contractFat(nc,[&](int i,int u,int w,int x){return Vaaca(u,w,i,x);},
                       [&](int i,int v,int z,int y){return Tcaaa(i,v,z,y);});
        for(int u=0;u<na;u++) for(int w=0;w<na;w++) for(int x=0;x<na;x++)
        for(int v=0;v<na;v++) for(int z=0;z<na;z++) for(int y=0;y<na;y++)
            temp1[w*na+z]+=0.5*P[idx6(u,w,x,v,z,y)]*l2(x,y,u,v);
        // r8: +0.5 H2[uwix]T2[ivyz]L2[xyvu]  fat i: A[i](u,w,x) B[i](v,y,z)
        contractFat(nc,[&](int i,int u,int w,int x){return Vaaca(u,w,i,x);},
                       [&](int i,int v,int y,int z){return Tcaaa(i,v,y,z);});
        for(int u=0;u<na;u++) for(int w=0;w<na;w++) for(int x=0;x<na;x++)
        for(int v=0;v<na;v++) for(int y=0;y<na;y++) for(int z=0;z<na;z++)
            temp1[w*na+z]+=0.5*P[idx6(u,w,x,v,y,z)]*l2(x,y,v,u);
        // r9: +0.5 H2[avxy]S2[uwaz]L2[xyuv]  fat a: A[a](v,x,y) B[a](u,w,z)
        contractFat(nv,[&](int a,int v,int x,int y){return Vvaaa(a,v,x,y);},
                       [&](int a,int u,int w,int z){return Saava(u,w,a,z);});
        for(int v=0;v<na;v++) for(int x=0;x<na;x++) for(int y=0;y<na;y++)
        for(int u=0;u<na;u++) for(int w=0;w<na;w++) for(int z=0;z<na;z++)
            temp1[w*na+z]+=0.5*P[idx6(v,x,y,u,w,z)]*l2(x,y,u,v);
        // r10: -0.5 H2[uviy]S2[iwxz]L2[xyuv]  fat i: A[i](u,v,y) B[i](w,x,z)
        contractFat(nc,[&](int i,int u,int v,int y){return Vaaca(u,v,i,y);},
                       [&](int i,int w,int x,int z){return Scaaa(i,w,x,z);});
        for(int u=0;u<na;u++) for(int v=0;v<na;v++) for(int y=0;y<na;y++)
        for(int w=0;w<na;w++) for(int x=0;x<na;x++) for(int z=0;z<na;z++)
            temp1[w*na+z]-=0.5*P[idx6(u,v,y,w,x,z)]*l2(x,y,u,v);
    }

    // ================= C2 : H_T_C2a_smallS ladder/ring (direct into C2t) =================
    std::vector<double> C2t(na4,0.0);
    {
        gemm_aa(T2_aavv, Vt_vvaa, nv*nv, M);            // += H2[efxy]T2[uvef]
        for(size_t i=0;i<na4;i++) C2t[i]+=M[i];
        gemm_aa(Vt_aacc, T2_ccaa, nc*nc, M);            // += H2[uvmn]T2[mnxy]
        for(size_t i=0;i<na4;i++) C2t[i]+=M[i];
        // += H2[ewxy]T2[uvew] and H2[ewyx]T2[vuew]  (fat e,w) -> [(u,v)][(x,y)]
        gA.assign((size_t)nv*na*na2,0.0); gB.assign((size_t)nv*na*na2,0.0);
        for(int e=0;e<nv;e++) for(int w=0;w<na;w++) for(int u=0;u<na;u++) for(int v=0;v<na;v++)
            gA[(size_t)(e*na+w)*na2+(u*na+v)]=Taava(u,v,e,w);
        for(int e=0;e<nv;e++) for(int w=0;w<na;w++) for(int x=0;x<na;x++) for(int y=0;y<na;y++)
            gB[(size_t)(e*na+w)*na2+(x*na+y)]=Vvaaa(e,w,x,y);
        cx(gA,gB,nv*na,(int)na2,(int)na2,M);
        for(size_t i=0;i<na4;i++) C2t[i]+=M[i];
        for(int e=0;e<nv;e++) for(int w=0;w<na;w++) for(int u=0;u<na;u++) for(int v=0;v<na;v++)
            gA[(size_t)(e*na+w)*na2+(u*na+v)]=Taava(v,u,e,w);
        for(int e=0;e<nv;e++) for(int w=0;w<na;w++) for(int x=0;x<na;x++) for(int y=0;y<na;y++)
            gB[(size_t)(e*na+w)*na2+(x*na+y)]=Vvaaa(e,w,y,x);
        cx(gA,gB,nv*na,(int)na2,(int)na2,M);
        for(size_t i=0;i<na4;i++) C2t[i]+=M[i];
        // += H2[vumw]T2[mwyx] and H2[uvmw]T2[mwxy]  (fat m,w) -> [(u,v)][(x,y)]
        gA.assign((size_t)nc*na*na2,0.0); gB.assign((size_t)nc*na*na2,0.0);
        for(int m=0;m<nc;m++) for(int w=0;w<na;w++) for(int u=0;u<na;u++) for(int v=0;v<na;v++)
            gA[(size_t)(m*na+w)*na2+(u*na+v)]=Vaaca(v,u,m,w);
        for(int m=0;m<nc;m++) for(int w=0;w<na;w++) for(int x=0;x<na;x++) for(int y=0;y<na;y++)
            gB[(size_t)(m*na+w)*na2+(x*na+y)]=Tcaaa(m,w,y,x);
        cx(gA,gB,nc*na,(int)na2,(int)na2,M);
        for(size_t i=0;i<na4;i++) C2t[i]+=M[i];
        for(int m=0;m<nc;m++) for(int w=0;w<na;w++) for(int u=0;u<na;u++) for(int v=0;v<na;v++)
            gA[(size_t)(m*na+w)*na2+(u*na+v)]=Vaaca(u,v,m,w);
        for(int m=0;m<nc;m++) for(int w=0;w<na;w++) for(int x=0;x<na;x++) for(int y=0;y<na;y++)
            gB[(size_t)(m*na+w)*na2+(x*na+y)]=Tcaaa(m,w,x,y);
        cx(gA,gB,nc*na,(int)na2,(int)na2,M);
        for(size_t i=0;i<na4;i++) C2t[i]+=M[i];
    }

    // ============ C2 : particle/hole inner temp2, then C2t[uvxy]+=t; C2t[vuyx]+=t ============
    {
        std::vector<double> t2(na4,0.0);
        auto addM_uxvy=[&](double c){ for(int u=0;u<na;u++) for(int x=0;x<na;x++) for(int v=0;v<na;v++) for(int y=0;y<na;y++)
            t2[id4(u,v,x,y)]+=c*M[(size_t)(u*na+x)*na2+(v*na+y)]; };
        auto addM_vyux=[&](double c){ for(int v=0;v<na;v++) for(int y=0;y<na;y++) for(int u=0;u<na;u++) for(int x=0;x<na;x++)
            t2[id4(u,v,x,y)]+=c*M[(size_t)(v*na+y)*na2+(u*na+x)]; };
        // + H1[ax]T2[uvay]  (fat a) -> [(u,v,y)][x]
        gA.assign((size_t)nv*na3,0.0); gB.assign((size_t)nv*na,0.0);
        for(int a=0;a<nv;a++) for(int u=0;u<na;u++) for(int v=0;v<na;v++) for(int y=0;y<na;y++)
            gA[(size_t)a*na3+((u*na+v)*na+y)]=Taava(u,v,a,y);
        for(int a=0;a<nv;a++) for(int x=0;x<na;x++) gB[(size_t)a*na+x]=Fav(x,a);
        cx(gA,gB,nv,(int)na3,na,M);
        for(int u=0;u<na;u++) for(int v=0;v<na;v++) for(int x=0;x<na;x++) for(int y=0;y<na;y++)
            t2[id4(u,v,x,y)]+=M[(size_t)((u*na+v)*na+y)*na+x];
        // - H1[ui]T2[ivxy]  (fat i) -> [(v,x,y)][u]
        gA.assign((size_t)nc*na3,0.0); gB.assign((size_t)nc*na,0.0);
        for(int i=0;i<nc;i++) for(int v=0;v<na;v++) for(int x=0;x<na;x++) for(int y=0;y<na;y++)
            gA[(size_t)i*na3+((v*na+x)*na+y)]=Tcaaa(i,v,x,y);
        for(int i=0;i<nc;i++) for(int u=0;u<na;u++) gB[(size_t)i*na+u]=Fca(i,u);
        cx(gA,gB,nc,(int)na3,na,M);
        for(int u=0;u<na;u++) for(int v=0;v<na;v++) for(int x=0;x<na;x++) for(int y=0;y<na;y++)
            t2[id4(u,v,x,y)]-=M[(size_t)((v*na+x)*na+y)*na+u];
        // + T1[ua]H2[avxy]  (fat a) -> [u][(v,x,y)]
        gA.assign((size_t)nv*na,0.0); gB.assign((size_t)nv*na3,0.0);
        for(int a=0;a<nv;a++) for(int u=0;u<na;u++) gA[(size_t)a*na+u]=T1av(u,a);
        for(int a=0;a<nv;a++) for(int v=0;v<na;v++) for(int x=0;x<na;x++) for(int y=0;y<na;y++)
            gB[(size_t)a*na3+((v*na+x)*na+y)]=Vvaaa(a,v,x,y);
        cx(gA,gB,nv,na,(int)na3,M);
        for(int u=0;u<na;u++) for(int v=0;v<na;v++) for(int x=0;x<na;x++) for(int y=0;y<na;y++)
            t2[id4(u,v,x,y)]+=M[(size_t)u*na3+((v*na+x)*na+y)];
        // - T1[ix]H2[uviy]  (fat i) -> [x][(u,v,y)]
        gA.assign((size_t)nc*na,0.0); gB.assign((size_t)nc*na3,0.0);
        for(int i=0;i<nc;i++) for(int x=0;x<na;x++) gA[(size_t)i*na+x]=T1ca(i,x);
        for(int i=0;i<nc;i++) for(int u=0;u<na;u++) for(int v=0;v<na;v++) for(int y=0;y<na;y++)
            gB[(size_t)i*na3+((u*na+v)*na+y)]=Vaaca(u,v,i,y);
        cx(gA,gB,nc,na,(int)na3,M);
        for(int u=0;u<na;u++) for(int v=0;v<na;v++) for(int x=0;x<na;x++) for(int y=0;y<na;y++)
            t2[id4(u,v,x,y)]-=M[(size_t)x*na3+((u*na+v)*na+y)];
        // - 0.5 L1[wz]T2[vuaw]H2[azyx]  (fold L1 into T2 over w) -> [(a,z)][(u,v)]x[(a,z)][(x,y)]
        {
            std::vector<double> TL((size_t)na2*nv*na,0.0);   // [v][u][a][z]
            for(int v=0;v<na;v++) for(int u=0;u<na;u++) for(int a=0;a<nv;a++) for(int z=0;z<na;z++){
                double s=0; for(int w=0;w<na;w++) s+=Taava(v,u,a,w)*l1(w,z); TL[(((size_t)v*na+u)*nv+a)*na+z]=s;
            }
            gA.assign((size_t)nv*na*na2,0.0); gB.assign((size_t)nv*na*na2,0.0);
            for(int a=0;a<nv;a++) for(int z=0;z<na;z++) for(int u=0;u<na;u++) for(int v=0;v<na;v++)
                gA[(size_t)(a*na+z)*na2+(u*na+v)]=TL[(((size_t)v*na+u)*nv+a)*na+z];
            for(int a=0;a<nv;a++) for(int z=0;z<na;z++) for(int x=0;x<na;x++) for(int y=0;y<na;y++)
                gB[(size_t)(a*na+z)*na2+(x*na+y)]=Vvaaa(a,z,y,x);
            cx(gA,gB,nv*na,(int)na2,(int)na2,M);
            for(int u=0;u<na;u++) for(int v=0;v<na;v++) for(int x=0;x<na;x++) for(int y=0;y<na;y++)
                t2[id4(u,v,x,y)]-=0.5*M[(size_t)(u*na+v)*na2+(x*na+y)];
        }
        // - 0.5 Eta1[wz]T2[izyx]H2[vuiw]  (fold Eta1 into H2 over w) -> [(i,z)][(x,y)]x[(i,z)][(u,v)]
        {
            std::vector<double> HF((size_t)na2*nc*na,0.0);   // [v][u][i][z]
            for(int v=0;v<na;v++) for(int u=0;u<na;u++) for(int i=0;i<nc;i++) for(int z=0;z<na;z++){
                double s=0; for(int w=0;w<na;w++) s+=Vaaca(v,u,i,w)*e1(w,z); HF[(((size_t)v*na+u)*nc+i)*na+z]=s;
            }
            gA.assign((size_t)nc*na*na2,0.0); gB.assign((size_t)nc*na*na2,0.0);
            for(int i=0;i<nc;i++) for(int z=0;z<na;z++) for(int x=0;x<na;x++) for(int y=0;y<na;y++)
                gA[(size_t)(i*na+z)*na2+(x*na+y)]=Tcaaa(i,z,y,x);
            for(int i=0;i<nc;i++) for(int z=0;z<na;z++) for(int u=0;u<na;u++) for(int v=0;v<na;v++)
                gB[(size_t)(i*na+z)*na2+(u*na+v)]=HF[(((size_t)v*na+u)*nc+i)*na+z];
            cx(gA,gB,nc*na,(int)na2,(int)na2,M);
            for(int u=0;u<na;u++) for(int v=0;v<na;v++) for(int x=0;x<na;x++) for(int y=0;y<na;y++)
                t2[id4(u,v,x,y)]-=0.5*M[(size_t)(x*na+y)*na2+(u*na+v)];
        }
        // + H2[uexm]S2[vmye]  (fat e,m) -> [(u,x)][(v,y)]
        gA.assign((size_t)nv*nc*na2,0.0); gB.assign((size_t)nv*nc*na2,0.0);
        for(int e=0;e<nv;e++) for(int m=0;m<nc;m++) for(int u=0;u<na;u++) for(int x=0;x<na;x++)
            gA[(size_t)(e*nc+m)*na2+(u*na+x)]=Vavac(u,e,x,m);
        for(int e=0;e<nv;e++) for(int m=0;m<nc;m++) for(int v=0;v<na;v++) for(int y=0;y<na;y++)
            gB[(size_t)(e*nc+m)*na2+(v*na+y)]=Sacav(v,m,y,e);
        cx(gA,gB,nv*nc,(int)na2,(int)na2,M); addM_uxvy(1.0);
        // + H2[wumx]S2[mvwy]  (fat w,m) -> [(u,x)][(v,y)]
        gA.assign((size_t)na*nc*na2,0.0); gB.assign((size_t)na*nc*na2,0.0);
        for(int w=0;w<na;w++) for(int m=0;m<nc;m++) for(int u=0;u<na;u++) for(int x=0;x<na;x++)
            gA[(size_t)(w*nc+m)*na2+(u*na+x)]=Vaaca(w,u,m,x);
        for(int w=0;w<na;w++) for(int m=0;m<nc;m++) for(int v=0;v<na;v++) for(int y=0;y<na;y++)
            gB[(size_t)(w*nc+m)*na2+(v*na+y)]=Scaaa(m,v,w,y);
        cx(gA,gB,na*nc,(int)na2,(int)na2,M); addM_uxvy(1.0);
        // + 0.5 (L1[wz]S2[zvay])H2[auwx]  (fold L1 into S2 over z) -> [(a,w)][(u,x)]x[(a,w)][(v,y)]
        {
            std::vector<double> SL((size_t)na*na*nv*na,0.0);   // [w][v][a][y]
            for(int w=0;w<na;w++) for(int v=0;v<na;v++) for(int a=0;a<nv;a++) for(int y=0;y<na;y++){
                double s=0; for(int z=0;z<na;z++) s+=Saava(z,v,a,y)*l1(w,z); SL[(((size_t)w*na+v)*nv+a)*na+y]=s;
            }
            gA.assign((size_t)nv*na*na2,0.0); gB.assign((size_t)nv*na*na2,0.0);
            for(int a=0;a<nv;a++) for(int w=0;w<na;w++) for(int u=0;u<na;u++) for(int x=0;x<na;x++)
                gA[(size_t)(a*na+w)*na2+(u*na+x)]=Vvaaa(a,u,w,x);
            for(int a=0;a<nv;a++) for(int w=0;w<na;w++) for(int v=0;v<na;v++) for(int y=0;y<na;y++)
                gB[(size_t)(a*na+w)*na2+(v*na+y)]=SL[(((size_t)w*na+v)*nv+a)*na+y];
            cx(gA,gB,nv*na,(int)na2,(int)na2,M); addM_uxvy(0.5);
        }
        // - 0.5 (L1[wz]S2[ivwy])H2[zuix]  (fold L1 into H2 over z) -> [(i,w)][(v,y)]x[(i,w)][(u,x)]
        {
            std::vector<double> HL((size_t)na*na*nc*na,0.0);   // [w][u][i][x]
            for(int w=0;w<na;w++) for(int u=0;u<na;u++) for(int i=0;i<nc;i++) for(int x=0;x<na;x++){
                double s=0; for(int z=0;z<na;z++) s+=Vaaca(z,u,i,x)*l1(w,z); HL[(((size_t)w*na+u)*nc+i)*na+x]=s;
            }
            gA.assign((size_t)nc*na*na2,0.0); gB.assign((size_t)nc*na*na2,0.0);
            for(int i=0;i<nc;i++) for(int w=0;w<na;w++) for(int v=0;v<na;v++) for(int y=0;y<na;y++)
                gA[(size_t)(i*na+w)*na2+(v*na+y)]=Scaaa(i,v,w,y);
            for(int i=0;i<nc;i++) for(int w=0;w<na;w++) for(int u=0;u<na;u++) for(int x=0;x<na;x++)
                gB[(size_t)(i*na+w)*na2+(u*na+x)]=HL[(((size_t)w*na+u)*nc+i)*na+x];
            cx(gA,gB,nc*na,(int)na2,(int)na2,M); addM_vyux(-0.5);
        }
        // - H2[uemx]T2[vmye]  (fat e,m) -> [(u,x)][(v,y)]
        gA.assign((size_t)nv*nc*na2,0.0); gB.assign((size_t)nv*nc*na2,0.0);
        for(int e=0;e<nv;e++) for(int m=0;m<nc;m++) for(int u=0;u<na;u++) for(int x=0;x<na;x++)
            gA[(size_t)(e*nc+m)*na2+(u*na+x)]=Vavca(u,e,m,x);
        for(int e=0;e<nv;e++) for(int m=0;m<nc;m++) for(int v=0;v<na;v++) for(int y=0;y<na;y++)
            gB[(size_t)(e*nc+m)*na2+(v*na+y)]=Tacav(v,m,y,e);
        cx(gA,gB,nv*nc,(int)na2,(int)na2,M); addM_uxvy(-1.0);
        // - H2[uwmx]T2[mvwy]  (fat w,m) -> [(u,x)][(v,y)]
        gA.assign((size_t)na*nc*na2,0.0); gB.assign((size_t)na*nc*na2,0.0);
        for(int w=0;w<na;w++) for(int m=0;m<nc;m++) for(int u=0;u<na;u++) for(int x=0;x<na;x++)
            gA[(size_t)(w*nc+m)*na2+(u*na+x)]=Vaaca(u,w,m,x);
        for(int w=0;w<na;w++) for(int m=0;m<nc;m++) for(int v=0;v<na;v++) for(int y=0;y<na;y++)
            gB[(size_t)(w*nc+m)*na2+(v*na+y)]=Tcaaa(m,v,w,y);
        cx(gA,gB,na*nc,(int)na2,(int)na2,M); addM_uxvy(-1.0);
        // - 0.5 (L1[wz]T2[zvay])H2[auxw]  (fold L1 into T2 over z) -> [(a,w)][(u,x)]x[(a,w)][(v,y)]
        {
            std::vector<double> TL((size_t)na*na*nv*na,0.0);   // [w][v][a][y]
            for(int w=0;w<na;w++) for(int v=0;v<na;v++) for(int a=0;a<nv;a++) for(int y=0;y<na;y++){
                double s=0; for(int z=0;z<na;z++) s+=Taava(z,v,a,y)*l1(w,z); TL[(((size_t)w*na+v)*nv+a)*na+y]=s;
            }
            gA.assign((size_t)nv*na*na2,0.0); gB.assign((size_t)nv*na*na2,0.0);
            for(int a=0;a<nv;a++) for(int w=0;w<na;w++) for(int u=0;u<na;u++) for(int x=0;x<na;x++)
                gA[(size_t)(a*na+w)*na2+(u*na+x)]=Vvaaa(a,u,x,w);
            for(int a=0;a<nv;a++) for(int w=0;w<na;w++) for(int v=0;v<na;v++) for(int y=0;y<na;y++)
                gB[(size_t)(a*na+w)*na2+(v*na+y)]=TL[(((size_t)w*na+v)*nv+a)*na+y];
            cx(gA,gB,nv*na,(int)na2,(int)na2,M); addM_uxvy(-0.5);
        }
        // + 0.5 (L1[wz]T2[ivwy])H2[uzix]  (fold L1 into H2 over z) -> [(i,w)][(v,y)]x[(i,w)][(u,x)]
        {
            std::vector<double> HL((size_t)na*na*nc*na,0.0);   // [w][u][i][x]
            for(int w=0;w<na;w++) for(int u=0;u<na;u++) for(int i=0;i<nc;i++) for(int x=0;x<na;x++){
                double s=0; for(int z=0;z<na;z++) s+=Vaaca(u,z,i,x)*l1(w,z); HL[(((size_t)w*na+u)*nc+i)*na+x]=s;
            }
            gA.assign((size_t)nc*na*na2,0.0); gB.assign((size_t)nc*na*na2,0.0);
            for(int i=0;i<nc;i++) for(int w=0;w<na;w++) for(int v=0;v<na;v++) for(int y=0;y<na;y++)
                gA[(size_t)(i*na+w)*na2+(v*na+y)]=Tcaaa(i,v,w,y);
            for(int i=0;i<nc;i++) for(int w=0;w<na;w++) for(int u=0;u<na;u++) for(int x=0;x<na;x++)
                gB[(size_t)(i*na+w)*na2+(u*na+x)]=HL[(((size_t)w*na+u)*nc+i)*na+x];
            cx(gA,gB,nc*na,(int)na2,(int)na2,M); addM_vyux(0.5);
        }
        // - H2[vemx]T2[muye]  (fat e,m) -> [(v,x)][(u,y)]
        gA.assign((size_t)nv*nc*na2,0.0); gB.assign((size_t)nv*nc*na2,0.0);
        for(int e=0;e<nv;e++) for(int m=0;m<nc;m++) for(int v=0;v<na;v++) for(int x=0;x<na;x++)
            gA[(size_t)(e*nc+m)*na2+(v*na+x)]=Vavca(v,e,m,x);
        for(int e=0;e<nv;e++) for(int m=0;m<nc;m++) for(int u=0;u<na;u++) for(int y=0;y<na;y++)
            gB[(size_t)(e*nc+m)*na2+(u*na+y)]=Tcaav(m,u,y,e);
        cx(gA,gB,nv*nc,(int)na2,(int)na2,M);
        for(int v=0;v<na;v++) for(int x=0;x<na;x++) for(int u=0;u<na;u++) for(int y=0;y<na;y++)
            t2[id4(u,v,x,y)]-=M[(size_t)(v*na+x)*na2+(u*na+y)];
        // - H2[vwmx]T2[muyw]  (fat w,m) -> [(v,x)][(u,y)]
        gA.assign((size_t)na*nc*na2,0.0); gB.assign((size_t)na*nc*na2,0.0);
        for(int w=0;w<na;w++) for(int m=0;m<nc;m++) for(int v=0;v<na;v++) for(int x=0;x<na;x++)
            gA[(size_t)(w*nc+m)*na2+(v*na+x)]=Vaaca(v,w,m,x);
        for(int w=0;w<na;w++) for(int m=0;m<nc;m++) for(int u=0;u<na;u++) for(int y=0;y<na;y++)
            gB[(size_t)(w*nc+m)*na2+(u*na+y)]=Tcaaa(m,u,y,w);
        cx(gA,gB,na*nc,(int)na2,(int)na2,M);
        for(int v=0;v<na;v++) for(int x=0;x<na;x++) for(int u=0;u<na;u++) for(int y=0;y<na;y++)
            t2[id4(u,v,x,y)]-=M[(size_t)(v*na+x)*na2+(u*na+y)];
        // - 0.5 (L1[wz]T2[uzay])H2[avxw]  (fold L1 into T2 over z) -> [(a,w)][(v,x)]x[(a,w)][(u,y)]
        {
            std::vector<double> TL((size_t)na*na*nv*na,0.0);   // [u][w][a][y]
            for(int u=0;u<na;u++) for(int w=0;w<na;w++) for(int a=0;a<nv;a++) for(int y=0;y<na;y++){
                double s=0; for(int z=0;z<na;z++) s+=Taava(u,z,a,y)*l1(w,z); TL[(((size_t)u*na+w)*nv+a)*na+y]=s;
            }
            gA.assign((size_t)nv*na*na2,0.0); gB.assign((size_t)nv*na*na2,0.0);
            for(int a=0;a<nv;a++) for(int w=0;w<na;w++) for(int v=0;v<na;v++) for(int x=0;x<na;x++)
                gA[(size_t)(a*na+w)*na2+(v*na+x)]=Vvaaa(a,v,x,w);
            for(int a=0;a<nv;a++) for(int w=0;w<na;w++) for(int u=0;u<na;u++) for(int y=0;y<na;y++)
                gB[(size_t)(a*na+w)*na2+(u*na+y)]=TL[(((size_t)u*na+w)*nv+a)*na+y];
            cx(gA,gB,nv*na,(int)na2,(int)na2,M);
            for(int v=0;v<na;v++) for(int x=0;x<na;x++) for(int u=0;u<na;u++) for(int y=0;y<na;y++)
                t2[id4(u,v,x,y)]-=0.5*M[(size_t)(v*na+x)*na2+(u*na+y)];
        }
        // + 0.5 (L1[wz]T2[iuyw])H2[vzix]  (fold L1 into H2 over z) -> [(i,w)][(v,x)]x[(i,w)][(u,y)]
        {
            std::vector<double> HL((size_t)na*na*nc*na,0.0);   // [w][v][i][x]
            for(int w=0;w<na;w++) for(int v=0;v<na;v++) for(int i=0;i<nc;i++) for(int x=0;x<na;x++){
                double s=0; for(int z=0;z<na;z++) s+=Vaaca(v,z,i,x)*l1(w,z); HL[(((size_t)w*na+v)*nc+i)*na+x]=s;
            }
            gA.assign((size_t)nc*na*na2,0.0); gB.assign((size_t)nc*na*na2,0.0);
            for(int i=0;i<nc;i++) for(int w=0;w<na;w++) for(int u=0;u<na;u++) for(int y=0;y<na;y++)
                gA[(size_t)(i*na+w)*na2+(u*na+y)]=Tcaaa(i,u,y,w);
            for(int i=0;i<nc;i++) for(int w=0;w<na;w++) for(int v=0;v<na;v++) for(int x=0;x<na;x++)
                gB[(size_t)(i*na+w)*na2+(v*na+x)]=HL[(((size_t)w*na+v)*nc+i)*na+x];
            cx(gA,gB,nc*na,(int)na2,(int)na2,M);
            for(int u=0;u<na;u++) for(int y=0;y<na;y++) for(int v=0;v<na;v++) for(int x=0;x<na;x++)
                t2[id4(u,v,x,y)]+=0.5*M[(size_t)(u*na+y)*na2+(v*na+x)];
        }
        // symmetrize the particle/hole block: C2t[uvxy]+=t2[uvxy]; C2t[vuyx]+=t2[uvxy]
        for(int u=0;u<na;u++) for(int v=0;v<na;v++) for(int x=0;x<na;x++) for(int y=0;y<na;y++){
            const double t=t2[id4(u,v,x,y)];
            C2t[id4(u,v,x,y)]+=t; C2t[id4(v,u,y,x)]+=t;
        }
    }

    // ---- seed (bare F_act diag, bare <uv|xy>) then symmetric-accumulation Hermitize ----
    Hbar1.assign(na2,0.0);
    for(int u=0;u<na;u++) Hbar1[u*na+u]=e_a[u];
    Hbar2.assign(na4,0.0);
    {
        std::vector<double> graw(na4,0.0);
        contract_aux(R->AA_RI_M, na*na, R->AA_RI_M, na*na, aux, graw.data());
        for(int u=0;u<na;u++) for(int v=0;v<na;v++) for(int x=0;x<na;x++) for(int y=0;y<na;y++)
            Hbar2[id4(u,v,x,y)]=graw[(size_t)(u*na+x)*na2+(v*na+y)];
    }
    for(int u=0;u<na;u++) for(int v=0;v<na;v++){
        const double t=temp1[u*na+v];
        Hbar1[u*na+v]+=0.5*t; Hbar1[v*na+u]+=0.5*t;
    }
    for(int u=0;u<na;u++) for(int v=0;v<na;v++) for(int x=0;x<na;x++) for(int y=0;y<na;y++){
        const double t=C2t[id4(u,v,x,y)];
        Hbar2[id4(u,v,x,y)]+=0.5*t; Hbar2[id4(x,y,u,v)]+=0.5*t;
    }
}

// ---- de-normal-order the active operator to a bare 0/1/2-body operator (Forte sadsrg.cc
//      :584, spin-free). Pair-symmetrize Hbar2 first so it carries the 4-fold symmetry,
//      then scalar -H1.L1 + 1/4 L1.hbar.L1 - 1/2 H2.L2 and 1-body -1/2 hbar.L1 with
//      hbar = 2*H2 - H2^swap. Hbar0 = E_corr (asymmetric <[H~1,T]>) from the ledger. ----

void dsrg_sf_tensors::degno(){
    const int na=n_a;
    const size_t na2=(size_t)na*na, na4=na2*na2;
    auto id4=[na](int a,int b,int c,int d){ return (((size_t)a*na+b)*na+c)*na+d; };

    // pair-symmetrize: <uv|xy> <- 1/2(<uv|xy> + <vu|yx>)
    {
        std::vector<double> t(Hbar2);
        for(int u=0;u<na;u++) for(int v=0;v<na;v++) for(int x=0;x<na;x++) for(int y=0;y<na;y++)
            Hbar2[id4(u,v,x,y)]=0.5*(t[id4(u,v,x,y)]+t[id4(v,u,y,x)]);
    }

    Hbar0_val = ledger.E_corr();

    // hbar[pqrs] = 2 H2[pqrs] - H2[pqsr]
    std::vector<double> hbar(na4,0.0);
    for(int p=0;p<na;p++) for(int q=0;q<na;q++) for(int r=0;r<na;r++) for(int s=0;s<na;s++)
        hbar[id4(p,q,r,s)]=2.0*Hbar2[id4(p,q,r,s)]-Hbar2[id4(p,q,s,r)];

    // scalar: -H1[vu]L1[uv] + 0.25 L1[uv]hbar[vyux]L1[xy] - 0.5 H2[xyuv]L2[uvxy]
    double s1=0.0;
    for(int u=0;u<na;u++) for(int v=0;v<na;v++) s1 -= Hbar1[v*na+u]*L1[u*na+v];
    double s2=0.0;
    for(int u=0;u<na;u++) for(int v=0;v<na;v++) for(int x=0;x<na;x++) for(int y=0;y<na;y++)
        s2 += 0.25*L1[u*na+v]*hbar[id4(v,y,u,x)]*L1[x*na+y];
    for(int x=0;x<na;x++) for(int y=0;y<na;y++) for(int u=0;u<na;u++) for(int v=0;v<na;v++)
        s2 -= 0.5*Hbar2[id4(x,y,u,v)]*SF_L2[id4(u,v,x,y)];
    e0_d = ledger.Eref + Hbar0_val + s1 + s2;

    // 1-body: h1_d[uv] = Hbar1[uv] - 0.5 hbar[uxvy]L1[yx]  (Hbar1 stays pre-deGNO)
    h1_d.assign(Hbar1.begin(), Hbar1.end());
    for(int u=0;u<na;u++) for(int v=0;v<na;v++){
        double a=0.0; for(int x=0;x<na;x++) for(int y=0;y<na;y++) a += hbar[id4(u,x,v,y)]*L1[y*na+x];
        h1_d[u*na+v] -= 0.5*a;
    }

    // chemist re-pack for the aldet import: (ab|cd)_chem = <ac|bd>_phys (inverse unpack)
    h2_d_chem.assign(na4,0.0);
    for(int a=0;a<na;a++) for(int b=0;b<na;b++) for(int c=0;c<na;c++) for(int d=0;d<na;d++)
        h2_d_chem[id4(a,b,c,d)]=Hbar2[id4(a,c,b,d)];
}

void dsrg_sf_tensors::set_root_data(const double* w, const double* E_bare,
                                    const double* E_dressed, const int* rmap, int ns){
    root_w.assign(w, w+ns);
    root_Ebare.assign(E_bare, E_bare+ns);
    root_Edressed.assign(E_dressed, E_dressed+ns);
    root_map.assign(rmap, rmap+ns);
    root_data_set = true;
}

// ---- DSRG_PILOT_DUMP: extends the DSRG_PRIM_DUMP whitespace-token grammar (cert 4.2).
//      One %.17e token per line; a stream error aborts (a truncated dump silently
//      poisons the pilot gate). ----

void dsrg_sf_tensors::pilot_dump(const char* path, int na_el, int nb_el, int ns, int root) const {
    const int na=n_a, nc=n_c, nv=n_v;
    const long aux = R->aux_n_ao;
    const size_t na2=(size_t)na*na, na4=na2*na2;

    FILE* f = fopen(path, "w");
    if(!f){ fprintf(out_stream, "ERROR: DSRG_PILOT_DUMP cannot open %s\n", path); abort(); }

    fprintf(f, "%d %d %d %d %d %d %d\n", nc, na, nv, na_el, nb_el, ns, root);
    fprintf(f, "%.17e\n", s_flow);
    fprintf(f, "%.17e\n", ledger.Eref);
    for(int m=0;m<nc;m++) fprintf(f, "%.17e\n", e_c[m]);
    for(int t=0;t<na;t++) fprintf(f, "%.17e\n", e_a[t]);
    for(int a=0;a<nv;a++) fprintf(f, "%.17e\n", e_v[a]);
    for(long i=0;i<(long)nc*na;i++) fprintf(f, "%.17e\n", f_ca[i]);
    for(long i=0;i<(long)nc*nv;i++) fprintf(f, "%.17e\n", f_cv[i]);
    for(long i=0;i<(long)na*nv;i++) fprintf(f, "%.17e\n", f_av[i]);
    for(long i=0;i<(long)nc;i++)   fprintf(f, "%.17e\n", dump_h_core.empty()?0.0:dump_h_core[i]);
    for(long i=0;i<(long)na*na;i++) fprintf(f, "%.17e\n", dump_h_act.empty()?0.0:dump_h_act[i]);
    fprintf(f, "%.17e\n", dump_e_scalar);
    fprintf(f, "%ld\n", aux);
    // RI B slabs (aux fastest) in NOPT layout; the pilot re-forms any (pq|rs) block.
    for(long i=0;i<(long)nv*nc*aux;i++) fprintf(f, "%.17e\n", R->VC_RI_M[i]);
    for(long i=0;i<(long)nc*na*aux;i++) fprintf(f, "%.17e\n", R->CA_RI_M[i]);
    for(long i=0;i<(long)nv*na*aux;i++) fprintf(f, "%.17e\n", R->VA_RI_M[i]);
    for(long i=0;i<(long)nc*nc*aux;i++) fprintf(f, "%.17e\n", R->CC_RI_M[i]);
    for(long i=0;i<(long)na*na*aux;i++) fprintf(f, "%.17e\n", R->AA_RI_M[i]);
    for(size_t i=0;i<na2;i++)  fprintf(f, "%.17e\n", L1[i]);
    for(size_t i=0;i<na4;i++)  fprintf(f, "%.17e\n", SF_L2[i]);
    for(size_t i=0;i<na4;i++)  fprintf(f, "%.17e\n", GAMMA_in[i]);
    for(int r=0;r<ns;r++) fprintf(f, "%.17e\n", om_v.empty()?0.0:om_v[r]);
    for(int r=0;r<ns;r++) fprintf(f, "%.17e\n", om_c.empty()?0.0:om_c[r]);
    // per-class energy ledger (Layer-A/Forte cross-check line)
    fprintf(f, "%.17e %.17e %.17e\n", ledger.E_FT1, ledger.E_FT2, ledger.E_VT1);
    fprintf(f, "%.17e %.17e %.17e %.17e %.17e\n",
            ledger.aavv(), ledger.ccaa(), ledger.caav(), ledger.caaa(), ledger.aaav());
    fprintf(f, "%.17e %.17e\n", ledger.VT2_L1_incore(), ledger.VT2_L2());
    fprintf(f, "%.17e %.17e %.17e\n", ledger.ccvv, ledger.cavv, ledger.ccav);
    fprintf(f, "%.17e %.17e %.17e\n", ledger.E3v, ledger.E3c, ledger.Eref);

    // ---- relaxation-rung tail (append-only after the 16 ledger tokens). has_hbar==0
    //      leaves the S2 parser section untouched; ==1 emits the dressed operator and
    //      the driver-set per-root data, aborting if that data is missing (partial dump). ----
    const int has_hbar = Hbar1.empty() ? 0 : 1;
    fprintf(f, "%d\n", has_hbar);
    if(has_hbar){
        if(!root_data_set || (int)root_w.size()!=ns || h2_d_chem.empty() || h1_d.empty()){
            fprintf(out_stream, "ERROR: DSRG_PILOT_DUMP has_hbar set but dressed/per-root "
                                "data missing on %s\n", path); abort();
        }
        for(int r=0;r<ns;r++) fprintf(f, "%.17e\n", root_w[r]);
        for(int r=0;r<ns;r++) fprintf(f, "%.17e\n", root_Ebare[r]);
        fprintf(f, "%.17e\n", Hbar0_val);
        for(size_t i=0;i<na2;i++) fprintf(f, "%.17e\n", Hbar1[i]);
        for(size_t i=0;i<na4;i++) fprintf(f, "%.17e\n", h2_d_chem[i]);
        fprintf(f, "%.17e\n", e0_d);
        for(size_t i=0;i<na2;i++) fprintf(f, "%.17e\n", h1_d[i]);
        for(int r=0;r<ns;r++) fprintf(f, "%.17e\n", root_Edressed[r]);
        for(int r=0;r<ns;r++) fprintf(f, "%d\n", root_map[r]);
    }

    if(ferror(f)){ fprintf(out_stream, "ERROR: DSRG_PILOT_DUMP write error on %s\n", path); abort(); }
    if(fclose(f)!=0){ fprintf(out_stream, "ERROR: DSRG_PILOT_DUMP close error on %s\n", path); abort(); }
    fprintf(out_stream, "  [dump] DSRG_PILOT_DUMP -> %s (nC=%d nA=%d nV=%d naux=%ld)\n",
            path, nc, na, nv, aux);
}
