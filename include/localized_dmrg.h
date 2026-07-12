#pragma once
//
// localized_dmrg — active-space basis transforms for localized DMRG-CASSCF.
//
// Given the localizing rotation U (n_act x n_act, row-major [a*n_act+p]) with C_loc = C*U and
// U^T U = I, these rotate the active-space integrals into the localized basis (forward) and the
// localized-basis RDMs back into the delocalized basis (backward). U is real & orthogonal, so the
// two directions are exact mutual inverses.
//
//   forward  (phi_deloc -> phi_loc,   integrals):  X_loc   = U^T X U
//   backward (phi_loc   -> phi_deloc, RDMs):        X_deloc = U   X U^T
//
// All buffers are flat row-major in
// NOPT's CAS layout: 1-body [t*n+u]; 2-body [((t*n+u)*n+v)*n+w] (chemist (tu|vw) / 2-RDM).

// Two-index transform of an n x n matrix: out = U^T X U (forward) or U X U^T (backward).
// X and out must not alias.
void rotate1(const double* X, const double* U, int n, double* out, bool forward);

// Four-index transform of the 2-body tensor G[((a*n+b)*n+c)*n+d] (chemist (ab|cd) integral or
// 2-RDM): four sequential quarter-transforms, O(n^5). out doubles as quarter-transform scratch
// (it is overwritten from the first pass on), so G and out must not alias.
void rotate2(const double* G, const double* U, int n, double* out, bool forward);

// Build the localized active orbitals C_loc = C_act * U for the diagnostic orbital dump.
//   C_act, C_loc: n_ao x n_act, row-major [ao*n_act+orb] (== CAS_engine::ACT_CVEC layout).
//   U: n_act x n_act, [a*n_act+p]. C_act and C_loc must not alias.
void build_loc_orbitals(const double* C_act, const double* U, int n_ao, int n_act, double* C_loc);
