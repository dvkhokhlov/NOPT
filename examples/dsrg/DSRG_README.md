# DSRG-PT2 on a CAS-SCF reference

NOPT can compute second-order driven-similarity-renormalization-group perturbation theory
(DSRG-PT2, Li & Evangelista) on top of a converged CAS-SCF reference. The method is
spin-free and density-fitted; the multi-state step is *uncontracted* — the dressed active
Hamiltonian is re-diagonalized in the full CAS determinant space.

Enable it with `DSRG=1` in `$PAR` plus a `$DSRG` group. `$CAS` is switched on automatically.

Minimal input:

```
$PAR RHF=1 DSRG=1 RI=1 NAME=job $PAREND
...
$act_space n_alp=1 n_bet=1 n_val=2 mult=1 $end
$CAS n_s=1 cisolver=aldet $END
$dsrg s=0.5 $end
```

## 1. Requirements

- **`RI=1` is mandatory.** The amplitudes and the batched CCVV/CAVV/CCAV terms are built
  from RI B-tensors; the driver exits if RI is off. Supply `_RI_BASIS` in `$MOL` (a
  `-jkfit` set for the orbital basis).
- **Either CI backend, energies and `relax=once` alike.** The engine reads the per-state
  spin-summed 2-RDM (`G2_calc_diag`) and 3-body moment (`G3_calc_diag`) through the CI seam,
  and both aldet and `cisolver=dmrg` provide them. The λ3 term rotates those tensors into the
  semicanonical active basis rather than the wavefunction, so nothing is asked of the solver's
  own basis; the relaxation step likewise moves the dressed *operator* to the solver (see §3).
- **All occupied orbitals are correlated.** There is no frozen-core option on this path;
  a Forte-style `frozen_docc` has no counterpart here.

## 2. `$DSRG` group keywords

Default shown in parentheses.

### Flow parameter and amplitudes

- **s** *(0.5)* — DSRG flow parameter in Eh⁻², must be > 0. The regularizer is
  `(1-e^{-sD²})/D`, so `s` sets a suppression edge at `1/sqrt(s)`: denominators well below
  it are damped, well above it recover the bare `1/D`. The useful zone is roughly
  `[0.1, 1.0]`; this is not enforced (see §3).
- **ccvv_source** *(normal)* — `normal | zero`. `zero` sets the DSRG source operator to zero
  for the amplitude blocks whose denominators carry core and virtual labels only (T1
  core→virtual, T2 CCVV), so those amplitudes take the bare `1/D` form instead of the
  regularized one. Equivalent to Forte's `CCVV_SOURCE ZERO`.

### Reference and states

- **root** *(0)* — zero-based state-specific root, must be `< $CAS n_s`. Perturbation is
  built on this root's densities.
- **sa** *(0)* — `0 | 1`. `1` builds the reference from the weighted ensemble of all `$CAS`
  states (weights taken from `$CAS w_state`, normalized): ensemble 1-/2-RDM → SA
  generalized Fock → ensemble amplitudes. Mutually exclusive with `root`. On its own it
  yields only an ensemble-averaged energy — pair it with `relax=once` (see §3).
- **relax** *(none)* — `none | once`. `once` performs the uncontracted multi-state step: the
  DSRG dressing is folded into a bare 0/1/2-body active operator, re-diagonalized in the
  CAS determinant space, and the dressed roots are matched back to the bare roots by
  CI-vector overlap. This is what produces **per-state** energies.

### Output

- **print** *(1)* — verbosity. `>= 2` appends the per-class ledger: in-core `[Vr,T2]`
  L1/L2 by excitation class (AAVV, CCAA, CAAV, CAAA, AAAV), the two λ3 contributions
  (E3 v/c), and the `[Fr,T1]`/`[Fr,T2]`/`[Vr,T1]` breakdown. Useful for localizing a
  disagreement to a single contraction class.

### Guards

Bad values exit before any work: `s <= 0`, unknown `ccvv_source` or `relax` value,
`root < 0` or `root >= n_s`, `sa` outside `{0,1}`, `sa=1` together with an explicit `root`,
and a zero total reference weight under `sa=1`. During the run, the semicanonical
block-diagonalization of the Fock matrix is verified and aborts if the residual
off-diagonal exceeds `1e-8`.

## 3. Hints & known problems

### Hints

- **`sa=1` without `relax=once` is not a chemically meaningful result.** It prints one
  number — the energy of the ensemble averaged over the reference states — which is neither
  a state energy nor an observable, and yields no excitation energies. It exists as a
  well-defined intermediate (and matches what Forte reports as its unrelaxed SA summary),
  but for a state-averaged calculation you want `sa=1 relax=once`.
- **DSRG-PT2 is more sensitive to CAS-SCF convergence than the reference energy is.** The
  denominators are orbital energies, so they respond at *first* order to a residual orbital
  rotation, while `E_ref` responds at second order. A reference that looks perfectly
  converged by its energy can still move the correlation terms: tightening `$CAS grad` from
  `1e-6` to `1e-8` shifts individual terms at the `1e-8` level and costs nothing. Do it for
  any comparison work.
- **Per-state energies inherit the CI convergence.** With the aldet default Davidson stop,
  the dressed per-root energies carry the CI residual; `$dav dr=1e-14` sharpens them when
  you are comparing states rather than looking at a single total.

### Known problems

- **Large `s` silently leaves the trustworthy envelope.** Nothing breaks loudly — over
  `s = 0.1 … 1e7` the method never produced a NaN, never diverged and never failed to
  converge. But the regularizer's peak gain grows as `0.638·sqrt(s)`, so at large `s` any
  difference in the reference orbitals is amplified accordingly, and past roughly `s ~ 1e5`
  the result becomes sensitive enough to reference convergence that two correct
  implementations of the same equations disagree far above their normal agreement level.
  The disagreement concentrates in the 2- and 3-body cumulant terms while the total still
  looks reasonable. The documented `[0.1, 1.0]` zone sits four decades inside this, so
  ordinary calculations are unaffected — but "turn the regularizer off by cranking `s`"
  gives a number whose last digits are not meaningful, and the output does not say so.
- **`relax=once` on a DMRG reference re-solves on the frozen lattice.** The bare states are
  snapshotted, the dressed operator is rotated to the solver's native basis and imported as a
  general MPO, and the re-solve warm-starts from the converged bare MPS; dressed roots are
  matched back by MPS overlap with the same `|overlap| > 0.9` ambiguity rule as aldet.
  `$DMRG warm_start=off` makes the dressed re-solve cold by design. `$DMRG localize=pm` needs
  no special handling anywhere on this path.
- **The λ3 term holds the 3-body moment explicitly.** Two `n_act⁶` buffers (the averaged
  moment and its semicanonical rotation) plus two more inside the contraction: 2 MB each at
  `n_act = 8`, 512 MB each at `n_act = 20`. On a DMRG reference that, and not the determinant
  space, is what caps the active space on this path.
