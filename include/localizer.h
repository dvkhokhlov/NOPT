#pragma once
//
// localizer — orbital localization producing the unitary U only (orbitals never mutated).
//
//   C_loc = C * U ;   phi_loc[p] = sum_a phi_deloc[a] * U[a,p] ;   U^T U = I.
//
// orbital_localizer is an abstract base: the non-virtual
// localize() handles subset / irrep plumbing and delegates the actual sweep to the protected
// localize_block(), which each scheme (Pipek-Mezey now, Boys later) implements. 

#include <vector>

class molecule;

// Iteration / convergence controls for a localization run.
struct loc_options {
    int    max_sweeps = 200;     // hard cap on Jacobi sweeps
    double grad_tol   = 1e-10;   // converged when max|B_st| < grad_tol ...
    double ftol       = 1e-12;   // ... OR the functional increment over a sweep < ftol
    bool   by_irrep   = false;   // PROVISIONED, NOT IMPLEMENTED YET.
                                 //   false: one block over ALL columns -> maximal locality,
                                 //          symmetry deliberately broken.
                                 //   true : per-rep_num blocks (symmetry preserved). Currently
                                 //          ERRORS LOUDLY (returns converged=false); reserved.
    bool   verbose    = false;   // trace L and max|B_st| per sweep to out_stream
    const double* U0  = nullptr; // optional n_sub^2 orthogonal seed; nullptr => identity
};

// Outcome of a localization run.
struct loc_result {
    bool   converged  = false;   // false => U is the U=I fallback (caller MUST check this)
    int    n_sweeps   = 0;
    double final_grad = 0.0;     // max|B_st| at the last sweep
    double L          = 0.0;     // final localization functional
};

// Abstract localizer. Subclasses implement only localize_block().
class orbital_localizer {
public:
    virtual ~orbital_localizer() = default;

    // Localize the orbital block C (n_ao x n_sub, AO x orbital, row-major [ao*n_sub+orb]).
    // Writes the orthogonal U (n_sub x n_sub, row-major [a*n_sub+p]); C_loc = C*U. C is never
    // modified.
    // orb_irrep: per-column irrep label (== rep_num+n_core slice), length n_sub. RESERVED —
    //   only consulted once by_irrep is implemented; pass nullptr for now.
    loc_result localize(const double* C, int n_ao, int n_sub,
                        const int* orb_irrep, double* U,
                        const loc_options& opt = loc_options());

protected:
    // Localize one orbital block. Ug enters pre-seeded orthogonal (identity or U0); the
    // implementation accumulates its rotations into Ug and must leave it orthogonal.
    virtual loc_result localize_block(const double* C_blk, int n_ao, int n_g,
                                      double* Ug, const loc_options& opt) = 0;
};

// Pipek-Mezey localizer. Constructed from a molecule (for S_AO + shells); caches only the AO->atom
// map (C-independent and fixed across macro-iterations).
// Maximizes L(U) = sum_A sum_i (Q^A_ii)^2, with the symmetrized Mulliken population matrix
//   Q^A_pq = 1/2 sum_{mu in A} ( C[mu,p] (SC)[mu,q] + (SC)[mu,p] C[mu,q] ),   SC = S_AO C,
// built directly in the active space — no dense n_ao x n_ao S^A per atom.
class pm_localizer : public orbital_localizer {
public:
    explicit pm_localizer(molecule& mol);

protected:
    loc_result localize_block(const double* C_blk, int n_ao, int n_g,
                              double* Ug, const loc_options& opt) override;

private:
    int n_ao_   = 0;
    int n_atoms_ = 0;
    const double* S_AO_ = nullptr;   // molecule-owned AO overlap (n_ao*n_ao), not owned
    std::vector<int> ao_atom_;       // AO -> atom, from the shell centers
};

// Pure 2x2 Jacobi sweep
loc_result pm_jacobi_sweep(double* Qpops, int n_atoms, int n_blk,
                           double* U, const loc_options& opt);
