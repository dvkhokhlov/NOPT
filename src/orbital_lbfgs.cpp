// Two-loop L-BFGS (Nocedal & Wright Alg. 7.4) on the orbital rotation parameters. The pair
// skip s.y > eps*y.y is the L-BFGS-B guard of Byrd, Lu, Nocedal and Zhu (1995).
# include <cmath>
# include <cstdio>
# include <cstdlib>
# include <algorithm>
# include <limits>
# include <utility>

# include "blas_link.h"
# include "orbital_lbfgs.h"


void orbital_lbfgs::init(int ext_depth, size_t ext_n_rot, double ext_x_max){

    depth = ext_depth;
    n_rot = ext_n_rot;
    x_max = ext_x_max;
    diag  = getenv("NOPT_SXPT_DIAG")!=nullptr;
    restart();
}

void orbital_lbfgs::restart(){

    S  .clear();
    Y  .clear();
    rho.clear();
    g_prev.clear();
    s_prev.clear();
    have_prev = false;
}

double orbital_lbfgs::direction(const std::vector<double>& g, const apply_fn& h0inv,
                                const project_fn& project, std::vector<double>& applied,
                                double t_max){

    const lapack_int n = (lapack_int)n_rot;

    // the pair of the step applied last, measured in the gradient it produced
    double sy = 0.0, yy = 0.0, cs = 0.0;
    if(have_prev){
        std::vector<double> y(g);
        cblas_daxpy(n, -1.0, g_prev.data(),1, y.data(),1);
        sy = cblas_ddot(n, s_prev.data(),1, y.data(),1);
        yy = cblas_ddot(n, y.data(),1, y.data(),1);
        const double ss = cblas_ddot(n, s_prev.data(),1, s_prev.data(),1);
        if(ss*yy>0.0) cs = sy/std::sqrt(ss*yy);
        if(sy > std::numeric_limits<double>::epsilon()*yy){
            if((int)S.size()==depth){
                S  .pop_front();
                Y  .pop_front();
                rho.pop_front();
            }
            S  .push_back(s_prev);
            Y  .push_back(std::move(y));
            rho.push_back(1.0/sy);
        }
    }

    // two-loop recursion; h0inv writes r = H0^-1 q into applied, which then becomes -H^-1 g
    const int m = (int)S.size();
    q = g;
    alpha.assign(m, 0.0);
    for(int i=m-1;i>=0;i--){
        alpha[i] = rho[i]*cblas_ddot(n, S[i].data(),1, q.data(),1);
        cblas_daxpy(n, -alpha[i], Y[i].data(),1, q.data(),1);
    }
    h0inv(q, applied);
    for(int i=0;i<m;i++){
        const double beta = rho[i]*cblas_ddot(n, Y[i].data(),1, applied.data(),1);
        cblas_daxpy(n, alpha[i]-beta, S[i].data(),1, applied.data(),1);
    }
    cblas_dscal(n, -1.0, applied.data(),1);

    // with no pair the direction is h0inv's own output, already in the caller's space
    if(m>0 && project) project(applied);

    // one trust region for the whole converger: the L-BFGS step obeys the same x_max the
    // bare amplitude does
    double mx = 0.0;
    for(size_t p=0;p<n_rot;p++) mx = std::max(mx, std::fabs(applied[p]));
    if(mx>x_max) cblas_dscal(n, x_max/mx, applied.data(),1);
    const double gdot = cblas_ddot(n, applied.data(),1, g.data(),1);

    s_prev = applied;
    g_prev = g;
    have_prev = true;

    if(diag)
        fprintf(out_stream," SX-PT LBFGS n %2d sy %.2e yy %.2e cos %.2e gdot %.2e pre %.2e r %.2e\n",
                m, sy, yy, cs, gdot, mx, mx/t_max);

    return mx;
}
