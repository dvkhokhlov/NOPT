#pragma once
//
// mps_rotation — rotate a state-averaged MultiMPS by a one-body unitary via an externally-driven
// multi-root TangentSpace time-evolution sweep.
//
// block2's built-in TimeEvolution on a MultiMPS is complex-only: it treats the two center
// wavefunctions as the real/imag parts of ONE complex vector and hard-asserts ket.size()==2, so it
// cannot evolve an N>2 state-averaged MultiMPS's physical roots. But a one-body propagator carries
// no cross-root coupling, so the multi-root sweep can be driven with public block2 primitives
// (multi_eff_ham -> per-root real Krylov expo -> shared multi-target density-matrix truncation).
// This is the only new capability the SA warm-start needs; everything else uses block2 as-is.

#include <memory>

#include "block2_core.hpp"
#include "block2_dmrg.hpp"

// Rotate `mps` in place by the one-body propagator encoded in `mpo_rot` (the anti-Hermitian NC MPO
// exp(-kappa t) built from kappa = log(U)), applied over t in [0,1] as `n_steps` TangentSpace TE
// sweeps at bond dimension `rot_m`. Works for any nroots (unlike block2's own MultiMPS TE). Returns
// the mean per-root norm^2 after the sweeps (~1 for a faithful unitary rotation; a large drift
// signals the basis change was too big to carry, and the caller should cold-restart). Switches the
// block2 threading seq_type to Simple internally (Tasked corrupts time evolution) and restores it.
double evolve_sa_multimps(
    const std::shared_ptr<block2::MultiMPS<block2::SU2, double>> &mps,
    const std::shared_ptr<block2::MPO<block2::SU2, double>> &mpo_rot,
    block2::ubond_t rot_m, double dt, int n_steps);
