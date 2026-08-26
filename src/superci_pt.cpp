// Super-CI-PT orbital converger: the first-order amplitudes of the Dyall H0 over the
// singly excited space, taken as the orbital rotation parameters. NOPT's diagonal B is
// twice the diagonal of the same zeroth-order matrix, so this is that preconditioner
// with the diagonal approximation removed.
# include <cmath>
# include <cstdlib>
# include <cstdio>
# include <cstring>
# include <algorithm>

# include "blas_link.h"
# include "matr.h"
# include "molecule.h"
# include "CAS.h"
# include "superci_pt.h"
# include "common_vars.h"

namespace {

// A metric direction below TAU_DROP is a linearly dependent SX state and carries no
// rotation. TAU_ADMIT re-admits it only well clear of the drop, so the kept set cannot
// flip between macro-iterations and make the amplitude jump.
const double TAU_DROP  = 1e-10;
const double TAU_ADMIT = 1e-9;

// max|offdiag| the canonicalization is allowed to leave, relative to the block's spread
const double CANON_TOL = 1e-9;

double max_abs(const double * x, long n){
    double m=0.0;
    for(long i=0;i<n;i++) m = std::max(m, std::fabs(x[i]));
    return m;
}

}   // namespace


int superci_pt_engine::init(int ext_n_c, int ext_n_a, int ext_n_v, int ext_n_ao,
                            const int * ext_rep_num, int ext_n_rep, double ext_x_max){

    n_c   = ext_n_c;
    n_a   = ext_n_a;
    n_v   = ext_n_v;
    n_ao  = ext_n_ao;
    n_mo  = n_c+n_a+n_v;
    n_rep = ext_n_rep;
    rep_num = ext_rep_num;
    x_max = ext_x_max;
    drop_reported = false;

    keep_p.assign(std::max(n_rep,1), -1);
    keep_h.assign(std::max(n_rep,1), -1);

    return 0;
}

double superci_pt_engine::calc(const double * G){

    return max_abs(G, (long)n_c*n_a + (long)n_c*n_v + (long)n_a*n_v);
}

// Window-local indices carrying irrep i_r. With symmetry off every rep_num is -1, so the
// caller asks for i_r=-1 and gets the whole window in one block.
void superci_pt_engine::members(int base, int dim, int i_r, std::vector<int>& mem) const{

    mem.clear();
    for(int p=0;p<dim;p++) if(rep_num[base+p]==i_r) mem.push_back(p);
}

// Per-irrep block diagonalization of F[n0..n0+dim), touching nothing: F, MO_VEC,
// rep_num and orb_energy are all left as they were.
void superci_pt_engine::canonicalize_block(const double * F, int n0, int dim,
                                           std::vector<double>& V, std::vector<double>& eps){

    V.assign((size_t)dim*dim, 0.0);
    eps.assign(std::max(dim,1), 0.0);
    if(dim==0) return;

    std::vector<int> reps;
    if(IS_SYM==0) reps.push_back(-1);
    else{
        for(int p=0;p<dim;p++){
            const int r = rep_num[n0+p];
            if(r<0 || r>=n_rep){
                fprintf(out_stream,"ERROR: converger=sxpt needs an irrep label on every optimized orbital,"
                                   " but MO %d carries rep_num=%d; use converger=soscf\n", n0+p, r);
                exit(EXIT_FAILURE);
            }
        }
        for(int r=0;r<n_rep;r++) reps.push_back(r);
    }

    std::vector<int> mem;
    std::vector<double> blk, ev;
    for(size_t ir=0; ir<reps.size(); ir++){
        members(n0, dim, reps[ir], mem);
        const int m = mem.size();
        if(m==0) continue;
        blk.assign((size_t)m*m, 0.0);
        ev .assign(m, 0.0);
        for(int i=0;i<m;i++)
        for(int j=0;j<m;j++)
            blk[(size_t)i*m+j] = F[(size_t)(n0+mem[i])*n_ao + n0+mem[j]];

        lapack_diag(blk.data(), ev.data(), m);   // rows of blk are the eigenvectors

        for(int mu=0;mu<m;mu++){
            eps[mem[mu]] = ev[mu];
            for(int j=0;j<m;j++)
                V[(size_t)mem[j]*dim + mem[mu]] = blk[(size_t)mu*m + j];
        }
    }
}

// max|offdiag(V^T F_blk V)|. Everything downstream assumes the two blocks are diagonal,
// and the transformation is silent when it is not.
double superci_pt_engine::canonical_residual(const double * F, int n0, int dim,
                                             const std::vector<double>& V){

    if(dim<2) return 0.0;

    buf1.assign((size_t)dim*dim, 0.0);
    buf2.assign((size_t)dim*dim, 0.0);
    for(int i=0;i<dim;i++)
    for(int j=0;j<dim;j++)
        buf1[(size_t)i*dim+j] = F[(size_t)(n0+i)*n_ao + n0+j];

    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, dim,dim,dim, 1.0,
                   buf1.data(),dim, V.data(),dim, 0.0, buf2.data(),dim);
    nopt_par_dgemm(CblasRowMajor,CblasTrans  ,CblasNoTrans, dim,dim,dim, 1.0,
                   V.data(),dim, buf2.data(),dim, 0.0, buf1.data(),dim);

    double off=0.0;
    for(int i=0;i<dim;i++)
    for(int j=0;j<dim;j++)
        if(i!=j) off = std::max(off, std::fabs(buf1[(size_t)i*dim+j]));

    return off;
}

// K[t,u] = -sum_w F^I_uw gamma_tw - sum_wxy (uw|xy) GAMMA_tw,xy, and Ktilde = K + 2 F_tot,
// both on the active-active block and both from the state-averaged RDMs.
void superci_pt_engine::build_koopmans(CAS_engine * CAS){

    const long na2 = (long)n_a*n_a;
    const long na3 = na2*n_a;

    K  .assign((size_t)na2, 0.0);
    K_t.assign((size_t)na2, 0.0);
    if(n_a==0) return;

    // F^I on the active block, laid out as [u*n_a+w] so the GEMM sees a contiguous operand
    buf1.assign((size_t)na2, 0.0);
    for(int u=0;u<n_a;u++)
    for(int w=0;w<n_a;w++)
        buf1[(size_t)u*n_a+w] = CAS->F_core_MO[(size_t)(n_c+u)*n_ao + n_c+w];

    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans, n_a,n_a,n_a, -1.0,
                   CAS->gamma,n_a, buf1.data(),n_a, 0.0, K.data(),n_a);
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans, n_a,n_a,na3, -1.0,
                   CAS->GAMMA,na3, CAS->aaaa_ints,na3, 1.0, K.data(),n_a);

    for(int t=0;t<n_a;t++)
    for(int u=0;u<n_a;u++)
        K_t[(size_t)t*n_a+u] = K[(size_t)t*n_a+u] + 2.0*CAS->F_tot[(size_t)(n_c+t)*n_ao + n_c+u];
}

// M C_mu = sign * eig_mu * metric * C_mu with C^T metric C = 1, solved per active irrep
// through a canonical orthogonalization of the metric. Rows of C are the C_mu.
void superci_pt_engine::solve_pencil(const double * M_in, const double * metric, double sign,
                                     std::vector<int>& keep, std::vector<double>& C,
                                     std::vector<double>& eig, double & min_metric,
                                     const char * what){

    C.clear();
    eig.clear();
    min_metric = 2.0;
    if(n_a==0) return;

    std::vector<int> reps;
    if(IS_SYM==0) reps.push_back(-1);
    else          for(int r=0;r<n_rep;r++) reps.push_back(r);

    std::vector<int> mem;
    std::vector<double> g_r, n_r, X, sub, red, lam;
    for(size_t ir=0; ir<reps.size(); ir++){
        members(n_c, n_a, reps[ir], mem);
        const int m = mem.size();
        if(m==0) continue;

        g_r.assign((size_t)m*m, 0.0);
        n_r.assign(m, 0.0);
        for(int i=0;i<m;i++)
        for(int j=0;j<m;j++)
            g_r[(size_t)i*m+j] = metric[(size_t)mem[i]*n_a + mem[j]];
        lapack_diag(g_r.data(), n_r.data(), m);   // ascending
        min_metric = std::min(min_metric, n_r[0]);

        int n_lo=0, n_hi=0;
        for(int k=0;k<m;k++){
            if(n_r[k]>TAU_DROP ) n_lo++;
            if(n_r[k]>TAU_ADMIT) n_hi++;
        }
        int & prev = keep[IS_SYM==0 ? 0 : reps[ir]];
        int nk = n_lo;
        if(prev>=0 && n_lo>prev) nk = std::max(prev, n_hi);
        prev = nk;
        if(nk<m && !drop_reported){
            drop_reported = true;
            fprintf(out_stream,"NOTE: super-CI-PT drops %d %s direction(s) of the active metric"
                               " (smallest occupation %.2e); those rotations are frozen\n",
                               m-nk, what, n_r[0]);
        }
        if(nk==0) continue;

        // X[p*nk+j] = U_keep n_keep^-1/2, the kept metric eigenvectors as columns
        X.assign((size_t)m*nk, 0.0);
        for(int j=0;j<nk;j++){
            const int k = m-nk+j;
            const double s = 1.0/std::sqrt(n_r[k]);
            for(int p=0;p<m;p++) X[(size_t)p*nk+j] = g_r[(size_t)k*m+p]*s;
        }

        sub.assign((size_t)m*m, 0.0);
        for(int i=0;i<m;i++)
        for(int j=0;j<m;j++)
            sub[(size_t)i*m+j] = M_in[(size_t)mem[i]*n_a + mem[j]];

        std::vector<double> tmp((size_t)m*nk, 0.0);
        red.assign((size_t)nk*nk, 0.0);
        lam.assign(nk, 0.0);
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, m,nk,m, 1.0,
                    sub.data(),m, X.data(),nk, 0.0, tmp.data(),nk);
        cblas_dgemm(CblasRowMajor,CblasTrans  ,CblasNoTrans, nk,nk,m, 1.0,
                    X.data(),nk, tmp.data(),nk, 0.0, red.data(),nk);
        lapack_diag(red.data(), lam.data(), nk);   // rows of red are the eigenvectors

        const size_t base = C.size();
        C.resize(base + (size_t)nk*n_a, 0.0);
        for(int mu=0;mu<nk;mu++){
            eig.push_back(sign*lam[mu]);
            double * row = C.data() + base + (size_t)mu*n_a;
            for(int p=0;p<m;p++){
                double s=0.0;
                for(int j=0;j<nk;j++) s += X[(size_t)p*nk+j]*red[(size_t)mu*nk+j];
                row[mem[p]] = s;
            }
        }
    }
}

// Cayley: U = (I - kappa/2)^-1 (I + kappa/2), exactly orthogonal for antisymmetric kappa.
// MO_VEC holds orbitals as rows, so C_new = U C_old; only the n_mo optimized rows exist.
void superci_pt_engine::apply_rotation(CAS_engine * CAS){

    const size_t nn = (size_t)n_mo*n_mo;
    buf1.assign(nn, 0.0);
    buf2.assign(nn, 0.0);
    for(size_t i=0;i<nn;i++){
        buf1[i] = -0.5*kappa[i];
        buf2[i] =  0.5*kappa[i];
    }
    for(int i=0;i<n_mo;i++){
        buf1[(size_t)i*n_mo+i] += 1.0;
        buf2[(size_t)i*n_mo+i] += 1.0;
    }

    // A row-major buffer read column-major is its own transpose, and inversion commutes
    // with that flip, so buf1 comes back holding (I - kappa/2)^-1 row-major.
    lapack_int N = n_mo, info = 0, lwork = 2*N*N+6*N+1;
    std::vector<lapack_int> piv(std::max(n_mo,1));
    std::vector<double> work(lwork);
#ifdef _OPENBLAS
    LAPACK_dgetrf(&N,&N,buf1.data(),&N,piv.data(),&info);
    if(info==0) LAPACK_dgetri(&N,buf1.data(),&N,piv.data(),work.data(),&lwork,&info);
#endif
#ifdef _MKL
    DGETRF(&N,&N,buf1.data(),&N,piv.data(),&info);
    if(info==0) DGETRI(&N,buf1.data(),&N,piv.data(),work.data(),&lwork,&info);
#endif
    if(info!=0){
        fprintf(out_stream,"ERROR: super-CI-PT Cayley transform failed (LAPACK info=%d)\n",(int)info);
        exit(EXIT_FAILURE);
    }

    std::vector<double> U(nn, 0.0);
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, n_mo,n_mo,n_mo, 1.0,
                   buf1.data(),n_mo, buf2.data(),n_mo, 0.0, U.data(),n_mo);

    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, n_mo,n_ao,n_mo, 1.0,
                   U.data(),n_mo, CAS->MO_VEC,n_ao, 0.0, CAS->MO_BUF,n_ao);
    memcpy(CAS->MO_VEC, CAS->MO_BUF, (size_t)n_mo*n_ao*sizeof(double));
}

double superci_pt_engine::step(CAS_engine * CAS){

    const double * g_it = CAS->G;
    const double * g_ia = CAS->G + (size_t)n_c*n_a;
    const double * g_ta = CAS->G + (size_t)n_c*n_a + (size_t)n_c*n_v;

    // --- canonical frame -------------------------------------------------------------
    canonicalize_block(CAS->F_tot, 0       , n_c, V_c, eps_c);
    canonicalize_block(CAS->F_tot, n_c+n_a , n_v, V_v, eps_v);

    const double off_c = canonical_residual(CAS->F_tot, 0      , n_c, V_c);
    const double off_v = canonical_residual(CAS->F_tot, n_c+n_a, n_v, V_v);
    const double span  = std::max(1.0, std::max(max_abs(eps_c.data(),n_c),
                                                max_abs(eps_v.data(),n_v)));
    if(std::max(off_c,off_v) > CANON_TOL*span){
        fprintf(out_stream,"ERROR: super-CI-PT canonicalization left the Fock blocks non-diagonal"
                           " (core %.3e, virtual %.3e, tolerance %.3e)\n",
                           off_c, off_v, CANON_TOL*span);
        exit(EXIT_FAILURE);
    }

    // --- gradient into that frame ----------------------------------------------------
    gc_it.assign((size_t)n_c*n_a, 0.0);
    gc_ia.assign((size_t)n_c*n_v, 0.0);
    gc_ta.assign((size_t)n_a*n_v, 0.0);
    if(n_c&&n_a)
        nopt_par_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans, n_c,n_a,n_c, 1.0,
                       V_c.data(),n_c, g_it,n_a, 0.0, gc_it.data(),n_a);
    if(n_a&&n_v)
        nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, n_a,n_v,n_v, 1.0,
                       g_ta,n_v, V_v.data(),n_v, 0.0, gc_ta.data(),n_v);
    if(n_c&&n_v){
        buf1.assign((size_t)n_c*n_v, 0.0);
        nopt_par_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans, n_c,n_v,n_c, 1.0,
                       V_c.data(),n_c, g_ia,n_v, 0.0, buf1.data(),n_v);
        nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, n_c,n_v,n_v, 1.0,
                       buf1.data(),n_v, V_v.data(),n_v, 0.0, gc_ia.data(),n_v);
    }

    // --- Koopmans matrices and the two pencils ---------------------------------------
    build_koopmans(CAS);

    std::vector<double> gam_h((size_t)n_a*n_a, 0.0);
    for(int t=0;t<n_a;t++)
    for(int u=0;u<n_a;u++)
        gam_h[(size_t)t*n_a+u] = (t==u?2.0:0.0) - CAS->gamma[(size_t)t*n_a+u];

    double min_p=2.0, min_h=2.0;
    solve_pencil(K  .data(), CAS->gamma  , -1.0, keep_p, C_p, e_p, min_p, "particle");
    solve_pencil(K_t.data(), gam_h.data(), +1.0, keep_h, C_h, e_h, min_h, "hole");

    const int nk_p = e_p.size();
    const int nk_h = e_h.size();

    // --- amplitudes ------------------------------------------------------------------
    T_it.assign((size_t)n_c*n_a, 0.0);
    T_ia.assign((size_t)n_c*n_v, 0.0);
    T_ta.assign((size_t)n_a*n_v, 0.0);
    double t_dir = 0.0;         // largest amplitude any single pencil direction carries

    std::vector<double> Tc_it((size_t)n_c*n_a, 0.0);
    std::vector<double> Tc_ia((size_t)n_c*n_v, 0.0);
    std::vector<double> Tc_ta((size_t)n_a*n_v, 0.0);

    for(int i=0;i<n_c;i++)
    for(int a=0;a<n_v;a++)
        Tc_ia[(size_t)i*n_v+a] = -gc_ia[(size_t)i*n_v+a]/(4.0*(eps_v[a]-eps_c[i]));

    if(nk_p&&n_v){
        std::vector<double> Y((size_t)nk_p*n_v, 0.0);
        nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, nk_p,n_v,n_a, 1.0,
                       C_p.data(),n_a, gc_ta.data(),n_v, 0.0, Y.data(),n_v);
        for(int mu=0;mu<nk_p;mu++)
        for(int a=0;a<n_v;a++)
            Y[(size_t)mu*n_v+a] *= -0.5/(eps_v[a]-e_p[mu]);
        nopt_par_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans, n_a,n_v,nk_p, 1.0,
                       C_p.data(),n_a, Y.data(),n_v, 0.0, Tc_ta.data(),n_v);
        for(int mu=0;mu<nk_p;mu++)
            t_dir = std::max(t_dir, max_abs(C_p.data()+(size_t)mu*n_a,n_a)
                                   *max_abs(Y  .data()+(size_t)mu*n_v,n_v));
    }

    if(nk_h&&n_c){
        std::vector<double> Y((size_t)nk_h*n_c, 0.0);
        nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans, nk_h,n_c,n_a, 1.0,
                       C_h.data(),n_a, gc_it.data(),n_a, 0.0, Y.data(),n_c);
        for(int mu=0;mu<nk_h;mu++)
        for(int i=0;i<n_c;i++)
            Y[(size_t)mu*n_c+i] *= -0.5/(e_h[mu]-eps_c[i]);
        nopt_par_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans, n_c,n_a,nk_h, 1.0,
                       Y.data(),n_c, C_h.data(),n_a, 0.0, Tc_it.data(),n_a);
        for(int mu=0;mu<nk_h;mu++)
            t_dir = std::max(t_dir, max_abs(C_h.data()+(size_t)mu*n_a,n_a)
                                   *max_abs(Y  .data()+(size_t)mu*n_c,n_c));
    }

    // --- back out of the canonical frame ---------------------------------------------
    if(n_c&&n_a)
        nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, n_c,n_a,n_c, 1.0,
                       V_c.data(),n_c, Tc_it.data(),n_a, 0.0, T_it.data(),n_a);
    if(n_a&&n_v)
        nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans, n_a,n_v,n_v, 1.0,
                       Tc_ta.data(),n_v, V_v.data(),n_v, 0.0, T_ta.data(),n_v);
    if(n_c&&n_v){
        buf1.assign((size_t)n_c*n_v, 0.0);
        nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, n_c,n_v,n_c, 1.0,
                       V_c.data(),n_c, Tc_ia.data(),n_v, 0.0, buf1.data(),n_v);
        nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans, n_c,n_v,n_v, 1.0,
                       buf1.data(),n_v, V_v.data(),n_v, 0.0, T_ia.data(),n_v);
    }

    // symmetry-forbidden amplitudes, mirroring calc_grad's masks
    for(int i=0;i<n_c;i++)
    for(int t=0;t<n_a;t++)
        if(rep_num[i]!=rep_num[n_c+t]) T_it[(size_t)i*n_a+t]=0.0;
    for(int i=0;i<n_c;i++)
    for(int a=0;a<n_v;a++)
        if(rep_num[i]!=rep_num[n_c+n_a+a]) T_ia[(size_t)i*n_v+a]=0.0;
    for(int t=0;t<n_a;t++)
    for(int a=0;a<n_v;a++)
        if(rep_num[n_c+t]!=rep_num[n_c+n_a+a]) T_ta[(size_t)t*n_v+a]=0.0;

    // --- trust cap, then the rotation ------------------------------------------------
    double mx = std::max(max_abs(T_it.data(),(long)n_c*n_a),
                std::max(max_abs(T_ia.data(),(long)n_c*n_v),
                         max_abs(T_ta.data(),(long)n_a*n_v)));
    if(mx>x_max){
        const double s = x_max/mx;
        for(size_t i=0;i<T_it.size();i++) T_it[i]*=s;
        for(size_t i=0;i<T_ia.size();i++) T_ia[i]*=s;
        for(size_t i=0;i<T_ta.size();i++) T_ta[i]*=s;
        fprintf(out_stream," SX-PT is scaling rotation angle matrix Xmax=%.5e                         |\n", mx);
        mx = x_max;
    }

    fprintf(out_stream," SX-PT: metric min %.1e/%.1e  kept %3d+%3d/%3d  max dir |T| %.1e       |\n",
                       min_p, min_h, nk_p, nk_h, n_a, t_dir);

    kappa.assign((size_t)n_mo*n_mo, 0.0);
    for(int i=0;i<n_c;i++)
    for(int t=0;t<n_a;t++){
        kappa[(size_t)i*n_mo + n_c+t] =  T_it[(size_t)i*n_a+t];
        kappa[(size_t)(n_c+t)*n_mo+i] = -T_it[(size_t)i*n_a+t];
    }
    for(int i=0;i<n_c;i++)
    for(int a=0;a<n_v;a++){
        kappa[(size_t)i*n_mo + n_c+n_a+a] =  T_ia[(size_t)i*n_v+a];
        kappa[(size_t)(n_c+n_a+a)*n_mo+i] = -T_ia[(size_t)i*n_v+a];
    }
    for(int t=0;t<n_a;t++)
    for(int a=0;a<n_v;a++){
        kappa[(size_t)(n_c+t)*n_mo + n_c+n_a+a] =  T_ta[(size_t)t*n_v+a];
        kappa[(size_t)(n_c+n_a+a)*n_mo + n_c+t] = -T_ta[(size_t)t*n_v+a];
    }

    apply_rotation(CAS);

    return mx;
}
