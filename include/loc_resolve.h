#pragma once
// loc_resolve — pins the rotation a localizer leaves undetermined. A converged Pipek-Mezey pair
// has curvature A_st = sum_A ( Q^A_st^2 - (Q^A_ss - Q^A_tt)^2/4 ) <= 0, and A_st ~ 0 is a
// direction the functional cannot see: an active space built from atomic shells leaves one such
// block per atom. Each flat block is canonicalized in a one-body operator instead.

struct loc_resolve_report {
    int    n_blocks  = 0;    // flat blocks larger than one orbital
    int    max_block = 0;    // orbitals in the largest of them
    double gap_in    = 0.0;  // largest -A_st inside a block
    double gap_out   = 0.0;  // smallest -A_st between blocks; the two bracket the threshold
    double dL        = 0.0;  // functional change, ~0 confirming the blocks were flat
    double ev_gap    = 0.0;  // smallest neighbouring eigenvalue gap of F inside a block; a gap at
                             // zero is a rotation F does not pin either (symmetry-degenerate
                             // shells), and is only free where the operator truly is
};

// Qpops: n_atoms population matrices (n x n) in the localized frame, as pm_jacobi_sweep leaves
// them. F: one-body operator (n x n) in the basis U starts from. U (n x n, [a*n+p], C_loc = C U)
// is updated in place, and left untouched when no block is flat.
loc_resolve_report resolve_flat_blocks(const double* Qpops, int n_atoms, int n,
                                       const double* F, double* U);
