// Kollmar's orbital DIIS: the error vector is the amplitude vector itself, which goes to
// zero at convergence, and the history is carried in one frame by accumulating Theta.
# include <cmath>
# include <algorithm>

# include "blas_link.h"
# include "matr.h"
# include "orbital_diis.h"

namespace {

// Smallest eigenvalue the Gram matrix, scaled to unit largest diagonal, may have before
// the oldest vector is dropped: below it the history is linearly dependent and tau keeps
// no significant digits.
const double DIIS_MIN_EIG = 1e-12;

}   // namespace


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
}

// tau minimizes <kappa_bar|kappa_bar> under sum(tau) = 1, i.e. tau = B^-1 1 / (1^T B^-1 1),
// taken from an eigendecomposition of B so the same factorization reports its conditioning.
void orbital_diis::solve_pulay(){

    while(true){
        const int n = (int)kap.size();
        B.assign((size_t)n*n, 0.0);
        for(int i=0;i<n;i++)
        for(int j=0;j<=i;j++){
            const double b = 2.0*cblas_ddot((lapack_int)n_rot, kap[i].data(),1, kap[j].data(),1);
            B[(size_t)i*n+j] = b;
            B[(size_t)j*n+i] = b;
        }
        double s = 0.0;
        for(int i=0;i<n;i++) s = std::max(s, B[(size_t)i*n+i]);
        if(s>0.0) for(size_t p=0;p<B.size();p++) B[p] /= s;

        ev.assign(n, 0.0);
        lapack_diag(B.data(), ev.data(), n);   // rows of B are the eigenvectors, ascending
        if(n==1 || ev[0]>DIIS_MIN_EIG) break;
        kap  .pop_front();
        theta.pop_front();
    }

    const int n = (int)kap.size();
    tau.assign(n, 0.0);
    if(n==1){ tau[0] = 1.0; return; }
    for(int mu=0;mu<n;mu++){
        double s = 0.0;
        for(int j=0;j<n;j++) s += B[(size_t)mu*n+j];
        s /= ev[mu];
        for(int i=0;i<n;i++) tau[i] += B[(size_t)mu*n+i]*s;
    }
    double denom = 0.0;
    for(int i=0;i<n;i++) denom += tau[i];
    for(int i=0;i<n;i++) tau[i] /= denom;
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
    if(kap.size()==1) applied = kappa;
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
