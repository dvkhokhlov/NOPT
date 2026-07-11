# DMRG-CASSCF backend (block2)

NOPT can drive the CAS-SCF active-space CI with a DMRG solver (the block2 backend)
instead of the native `aldet` full-CI. Select it per calculation with `cisolver=dmrg`
in `$CAS` plus a `$DMRG` group; the backend is opt-in at build time (off by default).

Minimal input:

```
$CAS n_s=1 w_state=all_1; cisolver=dmrg $END
$DMRG m=300 $END
```

See `examples/dmrg/{4ene,8ene,naphtalene,pBQ}` for worked RHF → CAS(aldet) → DMRG
inputs with reference logs.

## 1. Installing block2 compatible with NOPT

block2 must be built as a C++ shared library (`libblock2.so`) whose ABI matches the
flags NOPT compiles its block2 sources with. A mismatch either fails to link or
silently corrupts numerics.

### Build `libblock2.so`

From the block2-preview root:

```
mkdir build_clib && cd build_clib
cmake .. -DBUILD_LIB=OFF -DBUILD_CLIB=ON -DUSE_MKL=OFF -DLARGE_BOND=ON \
         -DF77UNDERSCORE=ON -DUSE_OPENBLAS=ON \
         -DCMAKE_INSTALL_PREFIX=/opt/block2
make -j$(nproc) && make install
```

- `BUILD_CLIB=ON` builds the C++ library + headers + CMake config (`BUILD_LIB` is the
  *python* extension — leave OFF). `cmake` configure needs a `python3` on PATH; runtime
  does not.
- `F77UNDERSCORE=ON` is mandatory with OpenBLAS: without it the library calls `dgemm`,
  OpenBLAS exports `dgemm_`, and NOPT won't link.
- `LARGE_BOND=ON` and `USE_OPENBLAS=ON` match NOPT's build flags. The library is
  real-`double` (`USE_COMPLEX=OFF`), which is all CASSCF needs.
- Installs to `/opt/block2`.

### Enable in NOPT

In your per-machine `Makefile.vars` (template: `Makefile_example.vars`):

```
USE_BLOCK2:=yes
BLOCK2_DIR:=/opt/block2
```

Off (`no`) by default: a normal build links no block2, and `cisolver=dmrg` then exits.
The Makefile applies the matching ABI macros (including the by-hand
`-D_USE_GLOBAL_VARIABLE`) to the block2 sources automatically. 

## 2. `$DMRG` group keywords

Default shown in parentheses.

### Core solve

- **m** *(required, > 0)* — MPS bond dimension; No default — must be
  set when `cisolver=dmrg`.
- **sweeps** *(40)* — maximum DMRG sweeps per CI solve.
- **sweep_tol** *(1e-8)* — sweep-to-sweep energy convergence threshold.
- **hf_occ** *(integral)* — initial MPS occupancy scheme. Only `integral` accepted.
- **schedule** *(default)* — bond-dimension / noise sweep schedule. Only `default` accepted.
- **save_dir** *(/dev/shm)* — block2 scratch root (renormalized operators, MPS tensors).
  RAM-backed by default.
- **memory** *(1.0)* — size of block2's double stack, in GB. This is what a large active space
  or bond dimension exhausts; raise it if block2 aborts on a stack overflow.

### Localization & ordering

- **localize** *(off)* — localize the active orbitals before the solve. `off | pm`
  (Pipek-Mezey). `boys` parses but is not implemented.
- **loc_order** *(fiedler)* — orbital (site) ordering for the MPS lattice.
  `fiedler | none`. `gaopt` parses but is not implemented.
- **dump_loc_orbs** *(off)* — flag: presence writes the localized orbitals (GAMESS format)
  at iteration 0, then continues. The `=on`/`=off` value is ignored — the keyword's
  presence alone enables it.

### Warm start (reuse the converged MPS across CAS-SCF macro-iterations)

- **warm_start** *(on)* — reuse the previous macro-iteration's MPS as the next guess.
  `on | off`.
- **warm_sweeps** *(0 = auto)* — max sweeps for a warm re-solve. `0` = auto = `sweeps/2`.
- **warm_start_after** *(1)* — cold macro-iterations before the localized frame is frozen
  and warm start begins.
- **warm_rotate** *(on)* — rotate the reused MPS into the new orbital basis.
  `on | off` (off = reuse tensors without rotation).
- **rot_m** *(0 = use m)* — bond dimension for the MPS-rotation time evolution.
- **rot_steps** *(1)* — time-evolution steps for the MPS rotation (dt = 1/rot_steps).

### Determinant read-out (leading CI configurations after convergence)

- **print_dets** *(on)* — report the leading determinants after convergence. `on | off`.
- **extract_m** *(0 = no compression)* — bond dimension the canonical MPS is compressed to
  before determinant enumeration. Compression is lossy (see below); the default extracts from
  the uncompressed MPS.
- **extract_cutoff** *(1e-3)* — magnitude cutoff for the determinant search; drops
  determinants below it.
- **det_rot_m** *(0 = auto)* — bond dimension for the localized→canonical rotation applied
  before read-out. auto = min(2m, 1500).
- **det_rot_steps** *(10)* — time-evolution steps for that read-out rotation.

## 3. Hints & known problems

### Hints

- **`localize=pm` is effectively required for convergence at moderate `m`.** PM-localizing
  the active space lowers MPS entanglement, so a modest bond dimension reaches FCI energy; 
- **`save_dir` defaults to `/dev/shm` (RAM).** Fast, but a large-`m` or long run can exhaust
  RAM; point it at a disk path for big active spaces.

- **CAS-SCF must be state-averaged (`method=1`, the default).** Separate minimization needs a
  2-RDM per state; the DMRG backend only forms their average. It is also not well defined here:
  block2 solves the state-averaged roots in one shared renormalized basis, so the individual
  roots are not variational and minimizing orbitals against them separately is not a stationary
  point of anything. `method != 1` with `cisolver=dmrg` is rejected.

### Known problems

- **Small `extract_m` corrupts the determinant dump.** Setting `extract_m` > 0 compresses the
  converged MPS before enumerating determinants; if `extract_m` ≪ `m` the compression is lossy
  and, in certain systems, can trigger a segfault on the block2 side. The default (`0`) does not
  compress. If you do set it, watch the captured-weight table printed under the determinants —
  a note is emitted when too little of the CI vector is captured for the leading-configuration
  picture to mean anything.
