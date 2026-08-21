#pragma once
// localized_dmrg — localized DMRG-CASSCF helpers. U is the localizing rotation
// (n_act x n_act, [a*n_act+p]), C_loc = C*U, U^T U = I; the tensor transforms it
// drives live in tensor_rotate.h.

// Build the localized active orbitals C_loc = C_act * U for the diagnostic orbital dump.
// C_act, C_loc: n_ao x n_act, [ao*n_act+orb] (== CAS_engine::ACT_CVEC). They must not alias.
void build_loc_orbitals(const double* C_act, const double* U, int n_ao, int n_act, double* C_loc);
