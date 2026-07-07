#pragma once
//
// Native-style leading-determinant read-out for the DMRG backend. The determinants are already
// in the canonical (delocalized MO) basis: the DMRG solve's localized MPS is rotated to the
// canonical reporting basis before extraction (see block2_casci_wrap::print_states), so the
// configurations are reported as-is, sorted by weight, in the aldet print_states format.

#include <cstdio>
#include <cstdint>
#include <vector>

struct leading_config {
    std::vector<uint8_t> occ;   // n_act entries, 0/1/2/3 (empty/alpha/beta/double), input orbital order
    double coef;
};

// Print one state's leading determinants in the native print_states format:
//   "State <s>  E  = <E> S^2 = <S2>:"  then the top-<print_number> rows  "<coef>  | <alpha> | <beta>",
// followed by a captured-weight diagnostic line. captured_weight is the fraction of the state norm
// carried by the read-out (rotation + compression + extraction-cutoff truncation); ~1 is faithful,
// a low value flags that the localized->canonical rotation bled significant weight past the bond dim.
void report_leading_configs(std::FILE *out, int state_idx, double E, double S2, int n_act,
                            const std::vector<leading_config> &configs,
                            int print_number, double captured_weight);
