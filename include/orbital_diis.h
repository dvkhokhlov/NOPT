#ifndef __orbital_diis
#define __orbital_diis

# include <deque>
# include <vector>

// Orbital DIIS over the rotation parameters. The history pairs every amplitude vector
// kappa^i with the accumulated rotation Theta^i of the orbitals it was computed in, and
// the applied rotation is -Theta^k + Theta_bar + kappa_bar. kappa must already carry the
// trust cap, or the stored Theta stops equalling the rotation that was applied.
class orbital_diis{
    public:
        void init(int ext_depth, size_t ext_n_rot, double ext_x_max);
        bool active() const { return depth>0; }
        int  in_use() const { return (int)kap.size(); }

        // applied = -Theta^k + Theta_bar + kappa_bar, or kappa itself with the history
        // restarted when the extrapolation leaves the trust region; must not alias kappa
        void extrapolate(const std::vector<double>& kappa, std::vector<double>& applied);

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
