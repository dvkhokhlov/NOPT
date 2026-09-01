#include "loc_resolve.h"

#include "blas_link.h"        // cblas_dgemm
#include "localized_dmrg.h"   // rotate1
#include "matr.h"             // lapack_diag, normalize_rotation_rows

#include <algorithm>
#include <cmath>
#include <vector>

// Two population profiles differing by d electrons give -A_st ~ d^2/2, so this admits a block
// whose profiles agree to about 0.14 e.
#define LOC_FLAT_TOL 1.0e-2

namespace {

// sum_A sum_k (Q^A_kk)^2 over one block, its orbitals rotated by V ([j*n_b+k], columns are the
// new orbitals). V = nullptr leaves the block as it stands.
double block_functional(const double* Qpops, int n_atoms, int n,
                        const std::vector<int>& idx, const double* V) {
    const int n_b = (int)idx.size();
    double L = 0.0;
    for (int at = 0; at < n_atoms; ++at) {
        const double* Q = Qpops + (size_t)at * n * n;
        for (int k = 0; k < n_b; ++k) {
            double q = 0.0;
            if (V == nullptr)
                q = Q[(size_t)idx[k] * n + idx[k]];
            else
                for (int i = 0; i < n_b; ++i)
                    for (int j = 0; j < n_b; ++j)
                        q += V[i * n_b + k] * Q[(size_t)idx[i] * n + idx[j]] * V[j * n_b + k];
            L += q * q;
        }
    }
    return L;
}

} // namespace

loc_resolve_report resolve_flat_blocks(const double* Qpops, int n_atoms, int n,
                                       const double* F, double* U) {
    loc_resolve_report rep;
    if (n < 2) return rep;

    // pair curvature of the converged functional
    std::vector<double> A((size_t)n * n, 0.0);
    for (int s = 0; s < n; ++s)
        for (int t = s + 1; t < n; ++t) {
            double a = 0.0;
            for (int at = 0; at < n_atoms; ++at) {
                const double* Q = Qpops + (size_t)at * n * n;
                const double d    = Q[(size_t)s * n + t];
                const double diff = Q[(size_t)s * n + s] - Q[(size_t)t * n + t];
                a += d * d - 0.25 * diff * diff;
            }
            A[(size_t)s * n + t] = A[(size_t)t * n + s] = a;
        }

    // blocks = connected components of the flat pairs
    std::vector<int> blk((size_t)n, -1);
    int n_comp = 0;
    for (int s = 0; s < n; ++s) {
        if (blk[s] >= 0) continue;
        blk[s] = n_comp;
        std::vector<int> stack(1, s);
        while (!stack.empty()) {
            const int p = stack.back();
            stack.pop_back();
            for (int q = 0; q < n; ++q)
                if (blk[q] < 0 && -A[(size_t)p * n + q] < LOC_FLAT_TOL) {
                    blk[q] = n_comp;
                    stack.push_back(q);
                }
        }
        n_comp++;
    }

    rep.gap_out = 1.0e300;
    for (int s = 0; s < n; ++s)
        for (int t = s + 1; t < n; ++t) {
            const double c = -A[(size_t)s * n + t];
            if (blk[s] == blk[t]) rep.gap_in  = std::max(rep.gap_in, c);
            else                  rep.gap_out = std::min(rep.gap_out, c);
        }
    if (rep.gap_out > 1.0e299) rep.gap_out = 0.0;   // a single block spans everything

    std::vector<double> F_loc((size_t)n * n);
    rotate1(F, U, n, F_loc.data(), /*forward=*/true);

    rep.ev_gap = 1.0e300;
    std::vector<double> Fb, Ub, Un, ev;
    for (int c = 0; c < n_comp; ++c) {
        std::vector<int> idx;
        for (int p = 0; p < n; ++p)
            if (blk[p] == c) idx.push_back(p);
        const int n_b = (int)idx.size();
        if (n_b < 2) continue;

        rep.n_blocks++;
        rep.max_block = std::max(rep.max_block, n_b);

        Fb.resize((size_t)n_b * n_b);
        ev.resize((size_t)n_b);
        for (int i = 0; i < n_b; ++i)
            for (int j = 0; j < n_b; ++j)
                Fb[(size_t)i * n_b + j] = F_loc[(size_t)idx[i] * n + idx[j]];

        // rows of Fb are the eigenvectors; the sign gauge keeps them reproducible across
        // macro-iterations, and V[j,k] = Fb[k,j] carries the block's orbitals into them
        lapack_diag(Fb.data(), ev.data(), n_b);
        normalize_rotation_rows(Fb.data(), n_b);
        for (int k = 0; k + 1 < n_b; ++k)
            rep.ev_gap = std::min(rep.ev_gap, ev[k + 1] - ev[k]);

        std::vector<double> V((size_t)n_b * n_b);
        for (int k = 0; k < n_b; ++k)
            for (int j = 0; j < n_b; ++j)
                V[(size_t)j * n_b + k] = Fb[(size_t)k * n_b + j];

        rep.dL += block_functional(Qpops, n_atoms, n, idx, V.data())
                - block_functional(Qpops, n_atoms, n, idx, nullptr);

        // U[:,idx] <- U[:,idx] V
        Ub.resize((size_t)n * n_b);
        Un.resize((size_t)n * n_b);
        for (int a = 0; a < n; ++a)
            for (int j = 0; j < n_b; ++j)
                Ub[(size_t)a * n_b + j] = U[(size_t)a * n + idx[j]];
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n_b, n_b, 1.0,
                    Ub.data(), n_b, V.data(), n_b, 0.0, Un.data(), n_b);
        for (int a = 0; a < n; ++a)
            for (int k = 0; k < n_b; ++k)
                U[(size_t)a * n + idx[k]] = Un[(size_t)a * n_b + k];
    }

    if (rep.ev_gap > 1.0e299) rep.ev_gap = 0.0;
    return rep;
}
