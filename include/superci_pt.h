#ifndef __superci_pt
#define __superci_pt

# include <vector>

# include "orbital_diis.h"

class CAS_engine;

// Perturbative super-CI orbital converger ($CAS converger=sxpt). The first-order
// amplitudes of the Dyall zeroth-order Hamiltonian are the orbital rotation
// parameters; two n_act x n_act pencils serve every core and every virtual orbital.
// State-averaged only -- the Koopmans matrices are built from the SA 1- and 2-RDM.
class superci_pt_engine{
    public:
        int    init(int ext_n_c, int ext_n_a, int ext_n_v, int ext_n_ao,
                    const int * ext_rep_num, int ext_n_rep, double ext_x_max,
                    int ext_diis);
        double calc(const double * G);      // max|g| over the three rotatable blocks
        double step(CAS_engine * CAS);      // build kappa, apply the rotation; returns max|kappa|
        double applied() const { return app_max; }   // max|rotation actually applied|
        void   reset_history();             // drop everything carried across macro-iterations

    private:
        int n_c, n_a, n_v, n_ao, n_mo, n_rep;
        const int * rep_num;
        double x_max;
        double app_max;
        bool   drop_reported;
        orbital_diis diis;

        // canonical frame: V[p*dim+mu] holds eigenvectors as columns, mu reusing the
        // window's own slots (no permutation), eps the matching orbital energies
        std::vector<double> V_c, V_v, eps_c, eps_v;
        std::vector<double> K, K_t;                  // Koopmans matrices, n_a x n_a
        std::vector<double> C_p, C_h, e_p, e_h;      // pencil solutions, rows are vectors
        std::vector<int>    keep_p, keep_h;          // kept counts per irrep, previous macro-iter
        std::vector<double> gc_it, gc_ia, gc_ta;     // gradient in the canonical frame
        std::vector<double> T_it, T_ia, T_ta;
        std::vector<double> kappa, buf1, buf2;
        std::vector<double> kap_vec, app_vec;

        void members(int base, int dim, int i_r, std::vector<int>& mem) const;
        void canonicalize_block(const double * F, int n0, int dim,
                                std::vector<double>& V, std::vector<double>& eps);
        double canonical_residual(const double * F, int n0, int dim,
                                  const std::vector<double>& V);
        double orthogonality_defect(int dim, const std::vector<double>& V);
        void build_koopmans(CAS_engine * CAS);
        void solve_pencil(const double * M_in, const double * metric, double sign,
                          std::vector<int>& keep, std::vector<double>& C,
                          std::vector<double>& eig, double & min_metric, const char * what);
        void apply_rotation(CAS_engine * CAS);
};

#endif
