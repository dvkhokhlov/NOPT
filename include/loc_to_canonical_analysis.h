#pragma once
//
// Leading-configuration analysis for the DMRG backend: report a state's dominant
// determinants in the canonical (delocalized) MO basis. The DMRG solve runs in a
// localized + Fiedler-ordered basis where individual configurations are not
// interpretable; the physically meaningful leading configs are the canonical-basis ones.
//
// Input is a truncated set of determinants extracted from the converged MPS, already in
// input orbital order with physical signs applied, each as a per-orbital occupation code
// (0 empty, 1 alpha, 2 beta, 3 doubly) and its coefficient. When the run is not localized
// the configs are already canonical and are reported as-is; when localized, the truncated
// CI vector is rotated localized->canonical before reporting (see the .cpp).

#include <cstdio>
#include <cstdint>
#include <vector>

struct leading_config {
    std::vector<uint8_t> occ;   // n_act entries, 0/1/2/3, in input orbital order
    double coef;
};

// Print one state's leading determinants in the native print_states format:
//   "State <s>  E  = <E> S^2 = <S2>:"  then rows  "<coef>  | <alpha> | <beta>".
// U (n_act x n_act, [new*n_act+old]) is the rotation from the extracted-determinant basis into
// the canonical reporting basis (localization + active canonicalization composed by the caller);
// pass nullptr for no rotation (configs already in the reporting basis). Coefficients use the
// standard spin-blocked determinant sign convention (<J|I> = det(U_alpha) det(U_beta)); magnitudes
// and configuration identities match the native aldet read-out, signs follow this convention.
void report_leading_configs(std::FILE *out, int state_idx, double E, double S2,
                            int n_act, int na, int nb,
                            const std::vector<leading_config> &configs,
                            const double *U, int print_number);
