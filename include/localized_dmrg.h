#pragma once
// localized_dmrg — active-space basis transforms for localized DMRG-CASSCF.
// U is the localizing rotation (n_act x n_act, [a*n_act+p]), C_loc = C*U, U^T U = I:
//   forward  (deloc -> loc, integrals):  X_loc   = U^T X U
//   backward (loc -> deloc, RDMs):       X_deloc = U   X U^T
// Buffers are flat row-major in NOPT's CAS layout: 1-body [t*n+u], 2-body [((t*n+u)*n+v)*n+w].

// Two-index transform of an n x n matrix. X and out must not alias.
void rotate1(const double* X, const double* U, int n, double* out, bool forward);

// Four-index transform of G[((a*n+b)*n+c)*n+d]: four quarter-transforms, O(n^5). out doubles as
// scratch, so G and out must not alias.
void rotate2(const double* G, const double* U, int n, double* out, bool forward);

// Build the localized active orbitals C_loc = C_act * U for the diagnostic orbital dump.
// C_act, C_loc: n_ao x n_act, [ao*n_act+orb] (== CAS_engine::ACT_CVEC). They must not alias.
void build_loc_orbitals(const double* C_act, const double* U, int n_ao, int n_act, double* C_loc);
