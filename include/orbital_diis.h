#ifndef __orbital_diis
#define __orbital_diis

# include <deque>
# include <vector>

// Orbital DIIS over the rotation parameters. The history pairs every amplitude vector
// kappa^i with the accumulated rotation Theta^i of the orbitals it was computed in, and
// the applied rotation is -Theta^k + Theta_bar + kappa_bar, scaled into the trust region.
// kappa must be the raw amplitude: capping it first makes the history near-degenerate.
class orbital_diis{
    public:
        void init(int ext_depth, size_t ext_n_rot, double ext_x_max);
        bool active() const { return depth>0; }
        int  in_use() const { return n_keep; }      // history directions the fit kept

        // applied = -Theta^k + Theta_bar + kappa_bar, scaled to x_max; returns its
        // magnitude before that scaling. Must not alias kappa.
        double extrapolate(const std::vector<double>& kappa, std::vector<double>& applied);

        // Drop the history: the surface it was accumulated on is gone.
        void restart();

    private:
        int    depth;
        size_t n_rot;
        double x_max;
        int    n_keep;

        std::deque< std::vector<double> > kap, theta;
        std::vector<double> theta_cur, tau, sig;

        void solve_pulay();
};

#endif
