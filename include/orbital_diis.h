#ifndef __orbital_diis
#define __orbital_diis

# include <deque>
# include <vector>

// Orbital DIIS over the rotation parameters. The history pairs every amplitude vector
// kappa^i with the accumulated rotation Theta^i of the orbitals it was computed in, and
// the applied rotation is -Theta^k + Theta_bar + kappa_bar, scaled into the trust region.
// kappa must be the raw amplitude: capping it first makes the Gram matrix singular.
class orbital_diis{
    public:
        void init(int ext_depth, size_t ext_n_rot, double ext_x_max);
        bool active() const { return depth>0; }
        int  in_use() const { return (int)kap.size(); }

        // applied = -Theta^k + Theta_bar + kappa_bar, scaled to x_max; returns its
        // magnitude before that scaling. Must not alias kappa.
        double extrapolate(const std::vector<double>& kappa, std::vector<double>& applied);

    private:
        int    depth;
        size_t n_rot;
        double x_max;

        std::deque< std::vector<double> > kap, theta;
        std::vector<double> theta_cur, tau, B, ev;

        void solve_pulay();
        void restart();
};

#endif
