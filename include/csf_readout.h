#pragma once
//
// CSF-native leading-configuration read-out for the DMRG backend. The DMRG solve's MPS is rotated by
// the proper (continuous) part of the solve->canonical basis change and sampled as spin-adapted CSFs
// (block2 DeterminantTRIE<SU2>, step-vector encoded). The discrete part of the basis change -- the
// lattice->canonical orbital permutation and the per-orbital sign the rotation cannot carry -- is not
// applied to the MPS here; it is reported alongside the CSFs (report_csf_orbital_map) so the consumer
// maps the step vectors (lattice order) onto the canonical orbitals exactly (fermionic reorder +
// per-orbital signs) when expanding to determinants.

#include <cstdio>
#include <cstdint>
#include <vector>

struct leading_csf {
    std::vector<uint8_t> step;  // n_act entries, 0/1/2/3 = empty/up/down/doubly, DMRG lattice order
    double coef;                // <CSF | Psi_state>
};

// Print the lattice->canonical orbital map and the per-orbital sign once, before the state blocks.
// canon_of_lattice[k] = canonical orbital index reported at lattice site k; lattice_sign[k] = +/-1 the
// sign to reinstate on singly-occupied site k.
void report_csf_orbital_map(std::FILE *out, int n_act,
                            const std::vector<int> &canon_of_lattice,
                            const std::vector<double> &lattice_sign);

// Print one state's leading CSFs in the native print_states header format:
//   "State <s>  E  = <E> S^2 = <S2>:"  then  "<coef>  | <step vector>"  rows (lattice order), where the
// step vector uses 0/u/d/2 = empty/up/down/doubly per site; then a captured-weight diagnostic. coef is
// the CSF (spin-adapted) amplitude, NOT a determinant amplitude -- expand each CSF to determinants
// (using the orbital map) to compare against aldet.
void report_leading_csfs(std::FILE *out, int state_idx, double E, double S2, int n_act,
                         const std::vector<leading_csf> &csfs,
                         int print_number, double captured_weight);
