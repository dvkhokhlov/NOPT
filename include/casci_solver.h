#pragma once
//
// casci_solver — abstract interface between the CAS-SCF orbital optimizer (CAS_engine)
// and a concrete active-space CI solver (built-in aldet determinant CI, DMRG, ...).
//
// The optimizer consumes only the RDMs (gamma/GAMMA); it never reads a CI internal.
// Any backend that can import the active-space Hamiltonian, solve, and hand back 1- and
// 2-RDMs in NOPT's convention can drive CAS-SCF through this one type.
//
// Method names deliberately mirror the existing aldet_data routines so the CAS_engine
// touch-point conversion is a mechanical CI[0].foo() -> CI->foo() rename (keeps the
// behavioural diff small and reviewable). Determinant-specific operations that do not
// generalise to an MPS backend are capability-gated via supports_civec_rotation().
//
// Pointer parameters stay raw (double*/int*) for zero-copy interop with the existing
// flat-array code; this interface introduces no ownership.

class aldet_data;  // opaque to consumers; only the aldet adapter dereferences it

class casci_solver {
public:
    virtual ~casci_solver();

    // --- configuration / lifecycle ---
    virtual void init_state_storage(int n_s, int i_set) = 0;   // allocate coef/E_states/... (aldet: init_zero_vec)
    virtual bool has_coef(int i_set) const = 0;                // is the CI vector storage allocated?
    virtual void set_act_rep_num(int* rep_num) = 0;            // per-active-orbital irrep numbers
    virtual void import_integrals(double* aaaa,                // active 2-e (tu|vw), chemist, n_act^4
                                  double* f_act,               // embedded 1-e h_tu (core folded in), n_act^2
                                  double e_core) = 0;          // inactive + nuclear scalar
    virtual void import_lambda(double* lambda_act,             // linear-molecule Lambda machinery (optional;
                               double lambda_core) {}          //   no-op unless the backend supports it)
    // Active-space localizing rotation U (n_act x n_act, [a*n_act+p], C_loc=C*U, U^T U=I). The
    // backend solves in the rotated basis and reports RDMs back in the original basis; nullptr or
    // never-called means solve in the supplied basis. aldet ignores it.
    virtual void set_localization(const double* U) {}

    // --- solve ---
    // Encapsulates the full diagonalisation (aldet: copy_coef -> set_par -> H_diag_calc -> run).
    // use_prev_guess: warm-start from the previous CI vector. Returns an iteration/convergence count.
    virtual int solve(int primary, int read, bool use_prev_guess) = 0;

    // --- reduced density matrices (the actual contract) ---
    // 1-RDM diagonal-form: spin-summed, symmetric; trace = N_active_el; each diagonal
    // entry (orbital occupation) satisfies 0 <= gamma_tt <= 2.
    virtual void calc_DM_diag(double* gamma, int a) = 0;
    // 2-RDM, NOPT convention: E2 = 1/2 * sum_{tuvw} GAMMA_{tuvw} (tu|vw),
    // packed GAMMA[((t*n_act+u)*n_act+v)*n_act+w].
    virtual void G_calc(double* GAMMA) = 0;
    // spin-resolved / transition 1-RDM blocks (properties): alpha and beta.
    virtual void calc_DMA(double* g, int a, int b) = 0;
    virtual void calc_DMB(double* g, int a, int b) = 0;

    // --- queries ---
    virtual int    n_states()       const = 0;
    virtual int    mult()           const = 0;
    virtual double E_state(int i)    const = 0;
    virtual double S2_state(int i)   const = 0;
    virtual double L2_state(int i)   const = 0;                // linear molecules
    virtual double P_state(int i)    const = 0;                // parity (linear molecules)
    virtual double* E_states_ptr() const = 0;                 // contiguous state-energy block (raw view, for PrintEnergy)
    // Energy convergence actually achieved by the last solve (DMRG: |dE| between the final two
    // sweeps; iterative CI: final residual). 0 if not tracked. Reported in the CAS-SCF table.
    virtual double last_solve_resid() const { return 0.0; }

    // --- relating the wavefunction across an active-orbital-basis change (capability-gated) ---
    // All three operations need the same thing: representing/comparing the wavefunction
    // when the active orbitals are rotated. A determinant CI has explicit, index-comparable
    // CI vectors, so it rotates them (malmqvist), follows roots by overlap against the
    // previous iteration's vectors (calc_S), and does the linear-molecule pi-pair rotation.
    // An MPS backend CAN do these but only at real cost: block2 rotates an MPS across orbital
    // bases via exp(kappa) time evolution (logm(U) -> anti-Hermitian 1-body MPO -> td_dmrg;
    // see block2-preview docs/.../orbital-rotation.rst), and cross-basis overlap then needs
    // that rotation plus an identity-MPO sweep. Cheap for small near-converged steps, costly
    // for large rotations (e.g. localization). So such a backend typically advertises false
    // and skips all three -- re-solving from a warm-started MPS instead -- but may opt in.
    virtual bool supports_civec_rotation() const { return false; }
    virtual void malmqvist(int i_set, double* U) {}                          // rotate CI vector by active-block U
    virtual void rotate_pi_pair(int i_set, double s, double c,               // linear-molecule pi-pair rotation
                                int pair, int* ind_pi) {}
    virtual void calc_S(double* S_track, int a, int b) {}                    // overlap vs previous iter's CI vectors

    // --- IO / diagnostics ---
    virtual void gen_ext_ind() = 0;
    virtual void print_states(int a, int n_s, int print) = 0;
    virtual void write_civec(int i_s, char* name) = 0;

    // --- escape hatch ---
    // Determinant-coupled paths outside the SCF loop (PT/XMCQDPT) keep using aldet_data
    // directly. The aldet adapter returns its wrapped object; other backends return nullptr.
    virtual aldet_data* as_aldet() { return nullptr; }
};
