// tensor_rotate — orthogonal basis transforms of flat active-space tensors
// See include/tensor_rotate.h for conventions

#include "tensor_rotate.h"

#include "blas_link.h"   // cblas_dgemm

#include <algorithm>
#include <vector>

namespace {

// Quarter-transform scratch of rotate2/rotate3, grown on demand and kept across calls: the caller's
// out is the other half of a two-buffer ping-pong, so no n^4 (n^6) buffer is allocated per call.
thread_local std::vector<double> qt_scratch;

// Cyclic index move on a tensor laid out [first][rest], rest = ncols contiguous: transpose the
// (n x ncols) view to (ncols x n), i.e. send the leading index to the back. Mirrors matr.cpp's
// transpose_A_to_B. Blocked over the long index so the gathered reads and the n-wide write band
// both stay resident; pure data movement.
void cycle_first_to_last(const double* in, double* out, int n, size_t ncols) {
    const size_t tile = 64;
    const long n_tiles = (long)((ncols + tile - 1) / tile);
    // Small tensors stay in cache and the fork/join costs more than the move (crossover ~n=12).
#pragma omp parallel for schedule(static) if (ncols * n >= 16384)
    for (long t = 0; t < n_tiles; ++t) {
        const size_t r0 = (size_t)t * tile, r1 = std::min(r0 + tile, ncols);
        for (int p = 0; p < n; ++p)
            for (size_t r = r0; r < r1; ++r)
                out[r * n + p] = in[(size_t)p * ncols + r];
    }
}

void cycle_first_to_last(const double* in, double* out, int n) {
    cycle_first_to_last(in, out, n, (size_t)n * n * n);
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
    if (qt_scratch.size() < n4) qt_scratch.resize(n4);
    double* b = qt_scratch.data();
    const double* src = G; // pass 0 reads the input directly (no initial copy)
    for (int pass = 0; pass < 4; ++pass) {
        // Quarter-transform the leading index of the (n x n^3) view.
        //   forward : b = U^T src  -> contracts U's first (deloc) index
        //   backward: b = U   src  -> contracts U's second (loc) index
        cblas_dgemm(CblasRowMajor, forward ? CblasTrans : CblasNoTrans, CblasNoTrans, n, (int)n3,
                    n, 1.0, U, n, src, (int)n3, 0.0, b, (int)n3);
        // Move the just-transformed index to the back; after 4 passes order is restored. out is the
        // scratch's ping-pong partner, so the last pass lands the result in place.
        cycle_first_to_last(b, out, n);
        src = out;
    }
}

void rotate3(const double* G, const double* U, int n, double* out, bool forward) {
    const size_t n5 = (size_t)n * n * n * n * n, n6 = n5 * n;
    if (qt_scratch.size() < n6) qt_scratch.resize(n6);
    double* b = qt_scratch.data();
    const double* src = G; // pass 0 reads the input directly (no initial copy)
    for (int pass = 0; pass < 6; ++pass) {
        // Quarter-transform the leading index of the (n x n^5) view.
        //   forward : b = U^T src  -> contracts U's first (deloc) index
        //   backward: b = U   src  -> contracts U's second (loc) index
        cblas_dgemm(CblasRowMajor, forward ? CblasTrans : CblasNoTrans, CblasNoTrans, n, (int)n5,
                    n, 1.0, U, n, src, (int)n5, 0.0, b, (int)n5);
        // Move the just-transformed index to the back; after 6 passes order is restored. out is the
        // scratch's ping-pong partner, so the last pass lands the result in place.
        cycle_first_to_last(b, out, n, n5);
        src = out;
    }
}
