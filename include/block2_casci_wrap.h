#pragma once
//
// block2_casci_wrap — casci_solver backend driving an external block2 DMRG-CI

#include <memory>

#include "casci_solver.h"
#include "inp_par_read.h"   // dmrg_par

struct dmrgci_engine;       // opaque; defined in the .cpp (holds all block2 state)

class block2_casci_wrap final : public casci_solver {
    std::unique_ptr<dmrgci_engine> impl_;   // pimpl; complete type only in the .cpp

public:
    // active-space dims (from M->CI[0]) + DMRG config; inits the block2 runtime once.
    block2_casci_wrap(int n_act, int na, int nb, int mult, int n_s, int print_number,
                      const dmrg_par& cfg);
    ~block2_casci_wrap() override;

    // --- configuration / lifecycle ---
    void init_state_storage(int n_s, int i_set) override;
    bool has_coef(int i_set) const override;
    void set_act_rep_num(int* rep_num) override;
    void set_localization(const double* U) override;
    void set_active_rotation(const double* R) override;
    void set_report_rotation(const double* U) override;
    void set_state_weights(const double* w, int n_s) override;
    void import_integrals(double* aaaa, double* f_act, double e_core) override;
    // Encode a TOTAL dressed active-space operator (F_act+g1, (tu|vw)+g2, g3, E_core+E0) as one
    // spin-adapted GeneralFCIDUMP -> GeneralMPO and swap it into the solve. Tensors arrive in the
    // frozen lattice (localized) basis; only the frozen Fiedler reorder is applied. h3 may be null.
    bool supports_dressed_import() const override { return true; }
    void import_dressed_operator(const double* h1_total, const double* h2_total,
                                 const double* h3_total, double const_total) override;

    // --- solve ---
    int solve(int primary, int read, bool use_prev_guess) override;

    // --- reduced density matrices ---
    void calc_DM_diag(double* gamma, int a) override;
    void G_calc(double* GAMMA) override;
    void calc_DMA(double* g, int a, int b) override;
    void calc_DMB(double* g, int a, int b) override;

    // --- queries ---
    int    n_states()      const override;
    int    mult()          const override;
    double E_state(int i)   const override;
    double S2_state(int i)  const override;
    double L2_state(int i)  const override;
    double P_state(int i)   const override;
    double* E_states_ptr()  const override;
    double last_solve_resid() const override;
    bool last_solve_hit_max() const override;

    // --- IO / diagnostics ---
    void gen_ext_ind() override;
    void print_states(int a, int n_s, int print) override;
    void write_civec(int i_s, char* name) override;

    // supports_civec_rotation() stays false (base default): DMRG re-solves rather than
    // rotating CI vectors; tracking/malmqvist/calc_S therefore stay no-ops. as_aldet()
    // stays nullptr
};
