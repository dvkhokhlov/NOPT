# AVAS active-space steering

AVAS (Atomic Valence Active Space) picks *which* orbitals land in the CAS active window by
their overlap with a set of atomic reference shells, instead of by orbital energy or by a
hand-written `reorder=1 orbitals=` list. It is opt-in: without an `$AVAS` group nothing in a
run changes.

Minimal input:

```
$PAR RHF=1 CAS=1 D5=1 RI=1 NAME=cr2 $PAREND
...
$act_space n_alp=6 n_bet=6 n_val=12 mult=1 $end
$AVAS atoms=1 2; shells=4s 3d; $END
```

## What it does

After the reference orbitals are available (RHF, or `MO_orth` when `RHF=0`), AVAS builds the
projector of the requested atomic shells onto the occupied and the virtual orbital block
separately, diagonalizes each one, and rotates the two blocks so that

The requested shells are first orthogonalized against the shells of the reference basis that lie
below them, over all selected atoms at once. In a segmented basis a shell overlaps the ones under
it (def2-SVP Cr: `<3s|4s>` = 0.40, `<3p|4p>` = 0.32), so without this the target span holds
semicore that the virtual block cannot reach and part of it is stranded in the core. How many
functions were removed is printed; if a requested shell turns out to lie inside them, the run
stops rather than inverting a singular target overlap.

- the σ-largest occupied orbitals become the **last** occupied orbitals, and
- the σ-largest virtual orbitals become the **first** virtual orbitals,

which is exactly the window `[n_core, n_core+n_val)` that `$act_space` defines. The eigenvalues
σ ∈ [0,1] measure how much of each rotated orbital lies in the reference span; they are printed
for both blocks with the selection boundary marked, and are stored in the orbital-energy field
so the dumped orbital files carry them.

The counts are **not** chosen by AVAS: `$act_space` stays authoritative. AVAS fills
`k_occ = (n_alp+n_bet)/2` occupied and `n_val - k_occ` virtual slots. If the forced counts cut
across a σ tier rather than at the largest gap of the spectrum, a `NOTE:` line says so and the
run proceeds as asked.

σ sums to the number of reference functions over both blocks, so the weight the window leaves
behind is reported too:

```
Target weight in the active space: 27.824 of 28
          left in occupied orbitals: 0.002
          left in virtual  orbitals: 0.076
```

One target direction can also be shared between an occupied and a virtual orbital, in which case
the two σ sum to one and the fixed window takes only the virtual half. Such pairs are counted in
a `NOTE:` line. They are not always a defect — a covalent target splits this way legitimately —
but they say the window `n_alp`/`n_bet` fixes is worth re-reading.

The rotated orbitals are always written as `<NAME>_AVAS.orb`, `<NAME>_AVAS.orb_GAMESS` and
`<NAME>_AVAS.out`, so a steered run can be inspected and restarted from its window.

## Keywords

- **atoms=** *(required, no default)* — 1-based indices of the atoms carrying the target
  shells, `;`-terminated: `atoms=1 2;`.
- **shells=** *(required, no default)* — nl labels applied to every selected atom,
  `;`-terminated: `shells=4s 3d;`. Within an atom the k-th reference shell of angular
  momentum l is the principal number n = k+l+1, so for a 3d metal `4s` and `3d` are the
  valence labels. A label that the reference basis does not carry for an element is an error.
- **ref_basis=** *(cc-pvtz-minao)* — the minimal basis the reference shells are taken from
  (H–Kr in the shipped library; it is also the SAD-guess basis).

Double-shell spaces come out on their own: the virtual σ tier contains both the antibonding
partners of the target shell and the radially-similar next shell (4d-like for a 3d reference),
so no extra reference functions are needed. The printed spectrum shows the tier structure.

## Restrictions

AVAS is rejected loudly, not silently ignored, when

- there is no `CAS=1` — nothing downstream consumes the steered window;
- `MP2=1` or `CIS=1` is set in the same run — both need canonical orbitals;
- `$act_space reorder=1` is set — two contradictory steering mechanisms;
- the point group is not C1 — the rotation mixes irreps;
- an ECP is in use — the reference shells are not all in the calculation basis.

Localization and DMRG orbital ordering (`localize=pm`, `loc_order=`) are unaffected: they take
the active window as given, so they compose with AVAS normally.
