// Leading-configuration analysis for the DMRG backend — see loc_to_canonical_analysis.h.

#include "loc_to_canonical_analysis.h"
#include "from_hash.h" // get_factorials / get_ind_from_ON / get_vec

#include <algorithm>
#include <cmath>
#include <numeric>
#include <string>

namespace {

long binomial(int n, int k) {
    if (k < 0 || k > n) return 0;
    return (long)std::lround(tgammal(n + 1) / tgammal(k + 1) / tgammal(n - k + 1));
}

// Determinant of a small n x n row-major matrix (destroys m); LU with partial pivoting.
double det_small(double *m, int n) {
    double det = 1.0;
    for (int k = 0; k < n; k++) {
        int piv = k;
        double mx = std::fabs(m[k * n + k]);
        for (int i = k + 1; i < n; i++) {
            double v = std::fabs(m[i * n + k]);
            if (v > mx) { mx = v; piv = i; }
        }
        if (mx == 0.0) return 0.0;
        if (piv != k) {
            for (int j = 0; j < n; j++) std::swap(m[k * n + j], m[piv * n + j]);
            det = -det;
        }
        det *= m[k * n + k];
        for (int i = k + 1; i < n; i++) {
            double f = m[i * n + k] / m[k * n + k];
            for (int j = k; j < n; j++) m[i * n + j] -= f * m[k * n + j];
        }
    }
    return det;
}

// Compound (nel-th exterior power) matrix of the orbital rotation U ([new*n_act+old]):
//   A[J,I] = det( U[ J_occ , I_occ ] )  over occupied-orbital sets of strings J (new) and I (old),
// i.e. the Slater-determinant overlap <J_new|I_old>. Row/col strings are the get_vec enumeration.
std::vector<double> compound_matrix(const double *U, int n_act, int nel, long Nx,
                                    const std::vector<int> &vecx) {
    std::vector<double> A((size_t)Nx * Nx, 0.0);
    std::vector<double> sub((size_t)nel * nel);
    if (nel == 0) { A[0] = 1.0; return A; } // single empty string, overlap 1
    for (long J = 0; J < Nx; J++)
        for (long I = 0; I < Nx; I++) {
            for (int r = 0; r < nel; r++) {
                const int jo = vecx[(size_t)J * (nel + 1) + r + 1] - 1;
                for (int c = 0; c < nel; c++) {
                    const int io = vecx[(size_t)I * (nel + 1) + c + 1] - 1;
                    sub[r * nel + c] = U[(size_t)jo * n_act + io];
                }
            }
            A[(size_t)J * Nx + I] = det_small(sub.data(), nel);
        }
    return A;
}

// Occupied-orbital list (get_vec order, 1-indexed ascending) -> "1010"-style 0/1 string,
// orbital 1 leftmost (matching aldet printf_occ_a/b).
void occ_string(long i_CI, int n_act, int n_el, const std::vector<int> &vecx, std::string &s) {
    s.assign(n_act, '0');
    for (int k = 1; k <= n_el; k++) s[vecx[(size_t)i_CI * (n_el + 1) + k] - 1] = '1';
}

} // namespace

void report_leading_configs(std::FILE *out, int state_idx, double E, double S2,
                            int n_act, int na, int nb,
                            const std::vector<leading_config> &configs,
                            const double *U, int print_number) {
    const long Na = binomial(n_act, na), Nb = binomial(n_act, nb);

    // Embed the extracted determinants into the dense CI grid ci[i_alpha*Nb + j_beta]
    // (alpha-major, beta fastest). The set is truncated; its captured weight (~1 for a
    // well-captured MPS) is reported as a diagnostic.
    std::vector<double> ci((size_t)Na * Nb, 0.0);
    std::vector<int> fa((size_t)na * n_act), fb((size_t)nb * n_act);
    if (na > 0) get_factorials(na, n_act, fa.data());
    if (nb > 0) get_factorials(nb, n_act, fb.data());
    std::vector<int> bit_a(n_act), bit_b(n_act), buf(std::max(na, nb) + 2);
    double w = 0.0;
    for (const auto &c : configs) {
        for (int o = 0; o < n_act; o++) {
            bit_a[o] = (c.occ[o] == 1 || c.occ[o] == 3);
            bit_b[o] = (c.occ[o] == 2 || c.occ[o] == 3);
        }
        buf[0] = 0;
        int i_CI = get_ind_from_ON(bit_a.data(), n_act, na, fa.data(), buf.data());
        buf[0] = 0;
        int j_CI = get_ind_from_ON(bit_b.data(), n_act, nb, fb.data(), buf.data());
        ci[(size_t)i_CI * Nb + j_CI] = c.coef;
        w += c.coef * c.coef;
    }

    std::vector<int> vec_a((size_t)Na * (na + 1)), vec_b((size_t)Nb * (nb + 1));
    get_vec(na, n_act, (int)Na, vec_a.data());
    get_vec(nb, n_act, (int)Nb, vec_b.data());

    // Rotate the (truncated) CI vector into the canonical reporting basis by the exact
    // string-factorized transform d[Ja,Jb] = sum_{Ia,Ib} A[Ja,Ia] ci[Ia,Ib] B[Jb,Ib], with
    // A, B the alpha/beta compound matrices of U. Norm-preserving for a unitary U, so the
    // reported amplitudes keep the extracted-MPS scale.
    if (U != nullptr) {
        std::vector<double> A = compound_matrix(U, n_act, na, Na, vec_a);
        std::vector<double> B = compound_matrix(U, n_act, nb, Nb, vec_b);
        std::vector<double> tmp((size_t)Na * Nb, 0.0); // tmp[Ja,Ib] = sum_Ia A[Ja,Ia] ci[Ia,Ib]
        for (long Ja = 0; Ja < Na; Ja++)
            for (long Ia = 0; Ia < Na; Ia++) {
                const double a = A[(size_t)Ja * Na + Ia];
                if (a == 0.0) continue;
                for (long Ib = 0; Ib < Nb; Ib++)
                    tmp[(size_t)Ja * Nb + Ib] += a * ci[(size_t)Ia * Nb + Ib];
            }
        std::fill(ci.begin(), ci.end(), 0.0); // d[Ja,Jb] = sum_Ib tmp[Ja,Ib] B[Jb,Ib]
        for (long Ja = 0; Ja < Na; Ja++)
            for (long Jb = 0; Jb < Nb; Jb++) {
                double d = 0.0;
                for (long Ib = 0; Ib < Nb; Ib++)
                    d += tmp[(size_t)Ja * Nb + Ib] * B[(size_t)Jb * Nb + Ib];
                ci[(size_t)Ja * Nb + Jb] = d;
            }
    }

    std::fprintf(out, "State %d  E  = % 18.10f S^2 = %.2f:\n", state_idx, E, S2);

    // Order the grid by descending |coef| and print the top print_number.
    const int top = (int)std::min<long>(print_number, Na * Nb);
    std::vector<int> ord((size_t)Na * Nb);
    std::iota(ord.begin(), ord.end(), 0);
    std::partial_sort(ord.begin(), ord.begin() + top, ord.end(),
                      [&](int i, int j) { return std::fabs(ci[i]) > std::fabs(ci[j]); });
    std::string a, b;
    for (int k = 0; k < top; k++) {
        const int idx = ord[k];
        occ_string(idx / (int)Nb, n_act, na, vec_a, a);
        occ_string(idx % (int)Nb, n_act, nb, vec_b, b);
        std::fprintf(out, "% .10e  | %s | %s\n", ci[idx], a.c_str(), b.c_str());
    }
    std::fprintf(out, "  (%d configs, captured weight %.6f)\n", (int)configs.size(), w);
}
