#pragma once
// mps_rotation — rotate a state-averaged MultiMPS by a one-body unitary, as an externally-driven
// multi-root TangentSpace TE sweep. block2's own MultiMPS TimeEvolution is complex-only (reads the
// two center wavefunctions as one complex vector, asserts ket.size()==2) and cannot evolve N>2
// physical roots. A one-body propagator has no cross-root coupling, so the sweep is driven from
// public primitives: multi_eff_ham -> per-root real Krylov expo -> multi-target DM truncation.

#include <cstdint>
#include <memory>
#include <vector>

#include "block2_core.hpp"
#include "block2_dmrg.hpp"

// Rotate `mps` in place by the one-body propagator in `mpo_rot` (the anti-Hermitian NC MPO
// exp(-kappa t), kappa = log(U)) over t in [0,1] as `n_steps` TangentSpace TE sweeps at bond
// dimension `rot_m`. Works for any nroots. Returns the mean per-root norm^2 after the sweeps (~1 for
// a faithful rotation; large drift means the basis change was too big to carry).
double evolve_sa_multimps(
    const std::shared_ptr<block2::MultiMPS<block2::SU2, double>> &mps,
    const std::shared_ptr<block2::MPO<block2::SU2, double>> &mpo_rot,
    block2::ubond_t rot_m, double dt, int n_steps);

// Outcome of apply_orbital_rotation_mps (below). Each caller applies its own thresholds/messages.
struct mps_rotation_result {
    double norm2 = 1.0;              // mean per-root norm^2 after the sweeps (drift => change too large)
    double im_norm = 0.0;           // |Im(log U)|; large => U was not a proper rotation
    bool complex_generator = false; // log U came back complex -> caller cold-starts / degrades
    bool skipped = false;           // |kappa| below threshold: exp(-kappa) ~ I, MPS left unrotated
};

// Rotate `mps` in place by a proper one-body active-orbital unitary U (n x n, [a*n+p]): kappa =
// log(U), reindexed into the frozen lattice (reorder_perm) order, applied as the anti-Hermitian NC
// MPO exp(-kappa) over `rot_steps` TE sweeps at bond dim `rot_m`. Shared by the warm-start MPS reuse
// and the determinant read-out. A near-identity U (kappa below threshold) is a no-op.
mps_rotation_result apply_orbital_rotation_mps(
    const std::shared_ptr<block2::MultiMPS<block2::SU2, double>> &mps,
    const double *U, int n, int n_elec, int twos,
    const std::vector<uint8_t> &orbsym, const std::vector<uint16_t> &reorder_perm,
    int rot_m, int rot_steps);
