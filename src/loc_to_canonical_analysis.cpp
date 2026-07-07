// Native-style leading-determinant read-out for the DMRG backend — see loc_to_canonical_analysis.h.

#include "loc_to_canonical_analysis.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <string>

namespace {

// Occupation codes (0/1/2/3, input order) -> alpha/beta "1010" strings, orbital 1 leftmost
// (matching aldet printf_occ_a/b).
void occ_codes_to_strings(const std::vector<uint8_t> &occ, int n_act, std::string &a, std::string &b) {
    a.assign(n_act, '0');
    b.assign(n_act, '0');
    for (int o = 0; o < n_act; o++) {
        if (occ[o] == 1 || occ[o] == 3) a[o] = '1';
        if (occ[o] == 2 || occ[o] == 3) b[o] = '1';
    }
}

} // namespace

void report_leading_configs(std::FILE *out, int state_idx, double E, double S2, int n_act,
                            const std::vector<leading_config> &configs,
                            int print_number, double captured_weight) {
    std::fprintf(out, "State %d  E  = % 18.10f S^2 = %.2f:\n", state_idx, E, S2);

    const int top = (int)std::min<size_t>(print_number, configs.size());
    std::vector<int> ord(configs.size());
    std::iota(ord.begin(), ord.end(), 0);
    std::partial_sort(ord.begin(), ord.begin() + top, ord.end(),
                      [&](int i, int j) { return std::fabs(configs[i].coef) > std::fabs(configs[j].coef); });
    std::string a, b;
    for (int k = 0; k < top; k++) {
        occ_codes_to_strings(configs[ord[k]].occ, n_act, a, b);
        std::fprintf(out, "% .10e  | %s | %s\n", configs[ord[k]].coef, a.c_str(), b.c_str());
    }
    std::fprintf(out, "  (captured weight %.6f)\n", captured_weight);
}
