// localized_dmrg — active-space basis transforms for the localized DMRG-CASSCF
// See include/localized_dmrg.h for conventions

#include "localized_dmrg.h"

#include "blas_link.h"   // cblas_dgemm

#include <cmath>
#include <vector>

namespace {

// Cyclic index move on a 4-index tensor laid out [first][rest], rest = n^3 contiguous:
// transpose the (n x n^3) view to (n^3 x n), i.e. send the leading index to the back
// (abcd -> bcda). Mirrors matr.cpp's transpose_A_to_B
void cycle_first_to_last(const double* in, double* out, int n) {
    const size_t n3 = (size_t)n * n * n;
    for (size_t r = 0; r < n3; ++r)
        for (int p = 0; p < n; ++p)
            out[r * n + p] = in[(size_t)p * n3 + r];
}

} // namespace

void rotate1(const double* X, const double* U, int n, double* out, bool forward) {
    std::vector<double> T((size_t)n * n);
    if (forward) { // out = U^T X U
        // T = U^T X
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, n, n, 1.0, U, n, X, n, 0.0,
                    T.data(), n);
        // out = T U
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, T.data(), n, U, n,
                    0.0, out, n);
    } else { // out = U X U^T
        // T = U X
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, U, n, X, n, 0.0,
                    T.data(), n);
        // out = T U^T
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, n, n, 1.0, T.data(), n, U, n, 0.0,
                    out, n);
    }
}

void rotate2(const double* G, const double* U, int n, double* out, bool forward) {
    const size_t n3 = (size_t)n * n * n, n4 = n3 * n;
    std::vector<double> a(n4), b(n4);
    const double* src = G; // pass 0 reads the input directly (no initial copy)
    for (int pass = 0; pass < 4; ++pass) {
        // Quarter-transform the leading index of the (n x n^3) view.
        //   forward : b = U^T src  -> contracts U's first (deloc) index
        //   backward: b = U   src  -> contracts U's second (loc) index
        cblas_dgemm(CblasRowMajor, forward ? CblasTrans : CblasNoTrans, CblasNoTrans, n, (int)n3,
                    n, 1.0, U, n, src, (int)n3, 0.0, b.data(), (int)n3);
        // Move the just-transformed index to the back; after 4 passes order is restored.
        cycle_first_to_last(b.data(), a.data(), n);
        src = a.data();
    }
    std::copy(a.begin(), a.end(), out);
}

double orthogonality_error(const double* U, int n) {
    double maxdev = 0.0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            double s = 0.0;
            for (int k = 0; k < n; ++k)
                s += U[k * n + i] * U[k * n + j]; // (U^T U)[i,j]
            double d = std::fabs(s - (i == j ? 1.0 : 0.0));
            if (d > maxdev) maxdev = d;
        }
    return maxdev;
}

void build_loc_orbitals(const double* C_act, const double* U, int n_ao, int n_act, double* C_loc) {
    // C_loc = C_act * U : (n_ao x n_act) = (n_ao x n_act)(n_act x n_act)
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_ao, n_act, n_act, 1.0, C_act, n_act,
                U, n_act, 0.0, C_loc, n_act);
}
