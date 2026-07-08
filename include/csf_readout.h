#pragma once
//
// CSF-native leading-configuration read-out for the DMRG backend. The DMRG solve's MPS is rotated by
// the proper (continuous) part of the solve->canonical basis change. The retained vector is extracted
// with block2's DeterminantTRIE<SU2> as step-vector CSFs, then locally recast into canonical orbital
// order/sign gauge. Determinants are still printed as a temporary diagnostic from that retained SU2
// vector.

#include <cstdio>
#include <cstdint>
#include <vector>

struct leading_csf {
    std::vector<uint8_t> step;  // n_act entries, 0/1/2/3 = empty/up/down/doubly
    double coef;                // <CSF | Psi_state>
};

// Print the lattice->canonical orbital map and the peeled residual-rotation sign diagnostic once,
// before the state blocks. canon_of_lattice[k] = canonical orbital index reported at lattice site k.
void report_csf_orbital_map(std::FILE *out, int n_act,
                            const std::vector<int> &canon_of_lattice,
                            const std::vector<double> &lattice_sign);

// Print one state's leading canonical CSFs. The step vector uses 0/u/d/2 =
// empty/up/down/doubly per canonical orbital. coef is the CSF (spin-adapted) amplitude, NOT a
// determinant amplitude.
void report_leading_csfs(std::FILE *out, int state_idx, double E, double S2, int n_act,
                         const std::vector<leading_csf> &csfs,
                         int print_number, double captured_weight);
