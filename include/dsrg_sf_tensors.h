#ifndef DSRG_SF_TENSORS_H
#define DSRG_SF_TENSORS_H

// State-specific unrelaxed spin-free DSRG-PT2 algebra core. Semicanonical only:
// the driver diagonalizes the SA generalized Fock, so eps are the diagonal and the
// three Fock cross blocks are the Hermitian off-diagonals. Plain (non-antisymmetrized)
// integrals, Eta1 = 2I-gamma, S2 = 2T-K -- Forte's spin-adapted realization. All 2-e
// data enters through the RI B-tensors (non-owning); densities through the seam GAMMA.
// The batched CCVV/CAVV/CCAV terms live in dsrg_sf_batch.{h,cpp} (I2); the driver adds
// them into the [Vr,T2] L1 sub-line.

#include <vector>

class RI_data;
class casci_solver;

// Forte six-line breakdown, per contraction class (state-specific, unrelaxed). The
// in-core [Vr,T2] is split both per amplitude-block class (Layer-A pilot rows) and
// per density factor (Forte L1/L2 sub-lines); L3 == E3v+E3c. The batched
// CCVV/CAVV/CCAV energies are computed in dsrg_sf_batch and stored by the driver.
struct dsrg_pt2_ledger {
    double E_FT1 = 0.0, E_FT2 = 0.0, E_VT1 = 0.0;
    double ccvv = 0.0, cavv = 0.0, ccav = 0.0;
    double aavv_L1 = 0.0, aavv_L2 = 0.0;
    double ccaa_L1 = 0.0, ccaa_L2 = 0.0;
    double caav_L1 = 0.0, caav_L2 = 0.0;
    double caaa_L1 = 0.0, caaa_L2 = 0.0;
    double aaav_L1 = 0.0, aaav_L2 = 0.0;
    double E3v = 0.0, E3c = 0.0;
    double Eref = 0.0;

    double aavv() const { return aavv_L1 + aavv_L2; }
    double ccaa() const { return ccaa_L1 + ccaa_L2; }
    double caav() const { return caav_L1 + caav_L2; }
    double caaa() const { return caaa_L1 + caaa_L2; }
    double aaav() const { return aaav_L1 + aaav_L2; }
    double VT2_L1_incore() const { return aavv_L1+ccaa_L1+caav_L1+caaa_L1+aaav_L1; }
    double VT2_L2() const { return aavv_L2+ccaa_L2+caav_L2+caaa_L2+aaav_L2; }
    // In-core correlation energy (adds the six-line minus the batched CCVV/CAVV/CCAV,
    // which the driver contributes separately).
    double E2_incore() const {
        return E_FT1+E_FT2+E_VT1 + VT2_L1_incore()+VT2_L2() + E3v+E3c;
    }
    // Full correlation energy = in-core plus the three batched DF terms (this is Hbar0).
    double E_corr() const { return E2_incore()+ccvv+cavv+ccav; }
};

class dsrg_sf_tensors {
public:
    dsrg_sf_tensors() = default;

    // One-shot problem setup (non-owning pointers; the buffers must outlive the engine).
    // eps_c/eps_a/eps_v: semicanonical Fock eigenvalues. f_ca[i*n_a+t], f_cv[i*n_v+a],
    // f_av[t*n_v+a]: the Hermitian Fock cross blocks. R: RI B-tensors in the same basis.
    void set_problem(int n_c, int n_a, int n_v, double s, bool ccvv_source_zero,
                     const double* eps_c, const double* eps_a, const double* eps_v,
                     const double* f_ca, const double* f_cv, const double* f_av,
                     const RI_data* R);
    // Densities in the semicanonical active basis. L1: spin-summed 1-RDM, trace =
    // N_act_el, layout [t*n_a+u]. GAMMA: seam 2-RDM, GAMMA[((t*n_a+u)*n_a+v)*n_a+w],
    // spin-summed. Builds Eta1, the physicist 2-RDM G2, and the SF_L2 cumulant.
    void set_densities(const double* L1, const double* GAMMA);

    // Build pipeline -- call in order.
    void build_vt();          // 6 renormalized v-tilde blocks from B
    void build_amplitudes();  // 7 T2 blocks (internal aaaa zeroed), 3 T1 blocks, Ft
    void compute_e2();        // Forte [F,T1]/[F,T2]/[V,T1] + certified in-core [V,T2]
    void compute_e3(casci_solver* CI, int root);        // DIRECT lambda3 (cert 2.6)
    // Reference energy (cert 2.7). h_core_diag: bare 1-e h on core (diagonal, n_c).
    // h_active: bare 1-e h active block, n_a*n_a. e_scalar: frozen-core + nuclear.
    double compute_eref(const double* h_core_diag, const double* h_active, double e_scalar);

    // ---- dressed active operator + de-normal-ordering (the relaxation rung) ----
    // Hbar1 = F_act(seed) + 1/2[H~1,A]_aa, Hbar2 = V_act(seed) + 1/2[H~1,A]_aaaa built
    // from the in-core term list (symmetric-accumulation Hermitization, alpha=1/2),
    // EXCLUDING the two batched CAVV/CCAV corrections the driver folds into hbar1_ref()
    // before degno(). Semicanonical active basis. Needs build_vt/build_amplitudes first.
    void build_hbar();

    // De-normal-order the MK-GNO active operator to a bare 0/1/2-body operator. Pair-
    // symmetrizes Hbar2 first, then the spin-free fold; fills e0_d, h1_d, and the chemist
    // re-pack h2_d_chem for the determinant re-solve. Hbar0 = E_corr from the ledger.
    void degno();

    // Mutable Hbar1 so the driver folds in the batched CAVV/CCAV corrections; the
    // pre-deGNO dressed blocks and the deGNO hand-off for import_integrals.
    std::vector<double>&       hbar1_ref()            { return Hbar1; }
    const std::vector<double>& hbar1()          const { return Hbar1; }
    const std::vector<double>& hbar2()          const { return Hbar2; }
    double                     hbar0()          const { return Hbar0_val; }
    double                     dressed_e0()     const { return e0_d; }
    const std::vector<double>& dressed_h1()     const { return h1_d; }
    const std::vector<double>& dressed_h2_chem()const { return h2_d_chem; }

    // Per-root data the driver snapshots for the pilot tail (reference weights, bare and
    // dressed per-root energies, overlap root map). Must be set before pilot_dump once
    // build_hbar has run, else the dump aborts (a partial dump poisons the gate).
    void set_root_data(const double* w, const double* E_bare,
                       const double* E_dressed, const int* root_map, int ns);

    // Env-gated DSRG_PILOT_DUMP writer (extends the DSRG_PRIM_DUMP grammar, cert 4.2).
    // Aborts if the stream errors; a truncated dump would silently poison the gate.
    void pilot_dump(const char* path, int na_el, int nb_el, int ns, int root) const;

    dsrg_pt2_ledger ledger;
    const std::vector<double>& omega_v() const { return om_v; }
    const std::vector<double>& omega_c() const { return om_c; }

private:
    int n_c = 0, n_a = 0, n_v = 0;
    double s_flow = 0.5;
    bool ccvv_zero = false;
    const double* e_c = nullptr;
    const double* e_a = nullptr;
    const double* e_v = nullptr;
    const double* f_ca = nullptr;  // [i*n_a+t]
    const double* f_cv = nullptr;  // [i*n_v+a]
    const double* f_av = nullptr;  // [t*n_v+a]
    const RI_data* R = nullptr;

    // densities (owned)
    std::vector<double> L1, Eta1;   // n_a^2
    std::vector<double> G2;         // n_a^4, physicist G2[p][q][r][s]
    std::vector<double> SF_L2;      // n_a^4, cumulant SF_L2[p][q][r][s]
    const double* GAMMA_in = nullptr;   // seam layout, kept for the dump

    // renormalized v-tilde blocks (physicist [bra1][bra2][ket1][ket2])
    std::vector<double> Vt_vvaa, Vt_aacc, Vt_avca, Vt_aaca, Vt_avaa, Vt_vaaa;
    // bare amplitude blocks (T2[hole1][hole2][part1][part2]); aaaa is internal (zero)
    std::vector<double> T2_aavv, T2_ccaa, T2_caav, T2_acav, T2_aava, T2_caaa;
    // singles + renormalized Fock cross blocks
    std::vector<double> T1_cv, T1_ca, T1_av;
    std::vector<double> Ft_cv, Ft_ca, Ft_av;
    // lambda3 solver overlaps per root (virtual/core branch)
    std::vector<double> om_v, om_c;

    // stashed for the dump (bare 1-e + scalar, from compute_eref)
    std::vector<double> dump_h_core, dump_h_act;
    double dump_e_scalar = 0.0;

    // dressed active operator (physicist) and the de-normal-ordered hand-off. Hbar2 is
    // pair-symmetrized in place by degno(); h1_d/h2_d_chem are the post-deGNO operators.
    std::vector<double> Hbar1, Hbar2;
    std::vector<double> h1_d, h2_d_chem;
    double Hbar0_val = 0.0, e0_d = 0.0;

    // driver-owned per-root data for the pilot tail (empty until set_root_data)
    std::vector<double> root_w, root_Ebare, root_Edressed;
    std::vector<int>    root_map;
    bool root_data_set = false;

    // file-local regularizer (Forte STD source; Taylor guard small=1e-3, order 4)
    double rs(double D) const;      // (1 - e^{-sD^2})/D  with the small-D Taylor branch
    double rfac(double D) const;    // e^{-sD^2}
};

#endif
