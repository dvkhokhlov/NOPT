// Kollmar's orbital DIIS: the error vector is the amplitude vector itself, which goes to
// zero at convergence, and the history is carried in one frame by accumulating Theta.
# include <cmath>
# include <cstdio>
# include <algorithm>
# include <limits>

# include "blas_link.h"
# include "matr.h"
# include "orbital_diis.h"


void orbital_diis::init(int ext_depth, size_t ext_n_rot, double ext_x_max){

    depth = ext_depth;
    n_rot = ext_n_rot;
    x_max = ext_x_max;
    restart();
}

void orbital_diis::restart(){

    kap  .clear();
    theta.clear();
    theta_cur.assign(n_rot, 0.0);
    n_keep = 0;
}

// tau minimizes <kappa_bar|kappa_bar> under sum(tau) = 1. Substitution and elimination:
// tau_p = 1 - sum_{k!=p} tau_k turns that into an unconstrained least-squares fit of -kappa^p
// by the differences kappa^k - kappa^p, whose solution inverts sigma once instead of twice.
// The constraint holds by construction, so there is no normalization that can vanish.
void orbital_diis::solve_pulay(){

    const lapack_int n = (lapack_int)kap.size();
    tau.assign(n, 0.0);
    n_keep = 1;
    if(n==1){ tau[0] = 1.0; return; }

    // the eliminated coefficient becomes the right-hand side, so it must be the smallest
    // column: the most recent amplitude. Eliminating the largest loses 3-6 digits.
    const lapack_int p = n-1;

    // column-major n_rot x (n-1) difference history, column j = kappa^{col[j]} - kappa^p;
    // dgesvd overwrites it
    lapack_int nc = n-1;
    std::vector<double> Et((size_t)n_rot*nc);
    std::vector<lapack_int> col(nc);
    for(lapack_int k=0,j=0;k<n;k++){
        if(k==p) continue;
        col[j] = k;
        std::copy(kap[k].begin(), kap[k].end(), Et.begin()+(size_t)j*n_rot);
        cblas_daxpy((lapack_int)n_rot, -1.0, kap[p].data(),1, Et.data()+(size_t)j*n_rot,1);
        j++;
    }

    lapack_int m = (lapack_int)n_rot, info = 0;
    lapack_int lwork = std::max(5*nc, 3*nc+m) + 64;
    std::vector<double> work((size_t)lwork), U((size_t)n_rot*nc), VT((size_t)nc*nc);
    sig.assign(nc, 0.0);
#ifdef _OPENBLAS
    LAPACK_dgesvd("S","S",&m,&nc,Et.data(),&m,sig.data(),U.data(),&m,
                  VT.data(),&nc,work.data(),&lwork,&info);
#endif
#ifdef _MKL
    DGESVD("S","S",&m,&nc,Et.data(),&m,sig.data(),U.data(),&m,
           VT.data(),&nc,work.data(),&lwork,&info);
#endif
    if(info!=0){
        fprintf(out_stream,"NOTE: super-CI-PT DIIS history SVD failed (info=%d);"
                           " taking the bare step this iteration\n",(int)info);
        tau[n-1] = 1.0;
        return;
    }

    // ctilde = -V S^-1 U^T kappa^p over the kept directions; V[j,mu] = VT[mu + nc*j]
    // (column-major, sigma descending). Cut at Etilde's numerical rank: below
    // max(m,n)*eps*sigma_max a direction is zero at working precision. Always keeps
    // sigma_max itself, so the fit never empties.
    const double cut = (double)std::max(m,nc)*std::numeric_limits<double>::epsilon()*sig[0];
    std::vector<double> y(nc, 0.0);
    cblas_dgemv(CblasColMajor, CblasTrans, m, nc, 1.0, U.data(), m,
                kap[p].data(),1, 0.0, y.data(),1);
    n_keep = 1;
    for(lapack_int mu=0;mu<nc;mu++){
        if(sig[mu]<=cut) continue;
        n_keep++;
        const double s = -y[mu]/sig[mu];
        for(lapack_int j=0;j<nc;j++) tau[col[j]] += VT[mu+(size_t)nc*j]*s;
    }
    double sum_c = 0.0;
    for(lapack_int j=0;j<nc;j++) sum_c += tau[col[j]];
    tau[p] = 1.0 - sum_c;
}

double orbital_diis::extrapolate(const std::vector<double>& kappa, std::vector<double>& applied){

    if((int)kap.size()==depth){
        kap  .pop_front();
        theta.pop_front();
    }
    kap  .push_back(kappa);
    theta.push_back(theta_cur);

    if(kap.size()>1) solve_pulay();

    // a single history entry forces tau = 1, so the extrapolation is exactly the identity
    if(kap.size()==1){ applied = kappa; n_keep = 1; }
    else{
        const int n = (int)kap.size();
        std::vector<double> tbar(n_rot, 0.0), kbar(n_rot, 0.0);
        for(int i=0;i<n;i++){
            cblas_daxpy((lapack_int)n_rot, tau[i], theta[i].data(),1, tbar.data(),1);
            cblas_daxpy((lapack_int)n_rot, tau[i], kap  [i].data(),1, kbar.data(),1);
        }
        const std::vector<double>& th = theta.back();
        applied.assign(n_rot, 0.0);
        for(size_t p=0;p<n_rot;p++) applied[p] = -th[p] + tbar[p] + kbar[p];
    }

    // one trust region for the whole converger: the extrapolated rotation obeys the same
    // x_max the bare amplitude does
    double mx = 0.0;
    for(size_t p=0;p<n_rot;p++) mx = std::max(mx, std::fabs(applied[p]));
    if(mx>x_max){
        const double s = x_max/mx;
        for(size_t p=0;p<n_rot;p++) applied[p] *= s;
    }

    // Theta carries the rotation that was applied -- identical to Theta_bar + kappa_bar
    // when nothing was scaled, and still exact when it was
    const std::vector<double>& th_k = theta.back();
    for(size_t p=0;p<n_rot;p++) theta_cur[p] = th_k[p] + applied[p];

    return mx;
}
