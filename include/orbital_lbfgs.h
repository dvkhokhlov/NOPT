#ifndef __orbital_lbfgs
#define __orbital_lbfgs

# include <deque>
# include <functional>
# include <vector>

// Orbital L-BFGS over the rotation parameters with the caller's inverse zeroth-order Hessian
// at the centre of the two-loop recursion. A pair is the applied step and the gradient change
// over it, admitted only when its curvature s.y is positive.
class orbital_lbfgs{
    public:
        using apply_fn   = std::function<void(const std::vector<double>&, std::vector<double>&)>;
        using project_fn = std::function<void(std::vector<double>&)>;

        void   init(int ext_depth, size_t ext_n_rot, double ext_x_max);
        bool   active() const { return depth>0; }
        bool   traced() const { return diag; }      // NOPT_SXPT_DIAG set: one trace line per step

        // Admits the pair (last applied step, g - previous g) when its curvature is positive,
        // forms d = -H^-1 g by the two-loop recursion with h0inv at the centre, projects it
        // (skipped while no pair is stored), scales it into the trust region and keeps it as
        // the next pair's step. applied must not alias g; t_max = max|T| feeds the trace only.
        double direction(const std::vector<double>& g, const apply_fn& h0inv,
                         const project_fn& project, std::vector<double>& applied, double t_max);
        void   restart();                            // drop the pairs and the previous point

    private:
        int    depth;
        size_t n_rot;
        double x_max;
        bool   diag;

        std::deque< std::vector<double> > S, Y;
        std::deque<double> rho;
        std::vector<double> g_prev, s_prev;
        bool   have_prev;
        std::vector<double> q, alpha;
};

#endif
