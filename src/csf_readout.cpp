// CSF-native leading-configuration read-out for the DMRG backend -- see csf_readout.h.

#include "csf_readout.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <string>

void report_csf_orbital_map(std::FILE *out, int n_act,
                            const std::vector<int> &canon_of_lattice,
                            const std::vector<double> &lattice_sign) {
    std::fprintf(out, "CSF read-out orbital map (lattice k -> canonical p, sign):\n");
    std::fprintf(out, "  canonical:");
    for (int k = 0; k < n_act; k++) std::fprintf(out, " %d", canon_of_lattice[k]);
    std::fprintf(out, "\n  sign     :");
    for (int k = 0; k < n_act; k++) std::fprintf(out, " %c", lattice_sign[k] < 0.0 ? '-' : '+');
    std::fprintf(out, "\n");
}

void report_leading_csfs(std::FILE *out, int state_idx, double E, double S2, int n_act,
                         const std::vector<leading_csf> &csfs,
                         int print_number, double captured_weight) {
    std::fprintf(out, "State %d  E  = % 18.10f S^2 = %.2f:\n", state_idx, E, S2);

    const int top = (int)std::min<size_t>(print_number, csfs.size());
    std::vector<int> ord(csfs.size());
    std::iota(ord.begin(), ord.end(), 0);
    std::partial_sort(ord.begin(), ord.begin() + top, ord.end(),
                      [&](int i, int j) { return std::fabs(csfs[i].coef) > std::fabs(csfs[j].coef); });
    static const char code[4] = {'0', 'u', 'd', '2'}; // empty / up / down / doubly
    std::string s;
    for (int k = 0; k < top; k++) {
        const leading_csf &c = csfs[ord[k]];
        s.assign(n_act, '0');
        for (int o = 0; o < n_act; o++) s[o] = code[c.step[o] & 3];
        std::fprintf(out, "% .10e  | %s\n", c.coef, s.c_str());
    }
    std::fprintf(out, "  (captured weight %.6f)\n", captured_weight);
}
