#pragma once
// tensor_rotate — orthogonal basis transforms of flat active-space tensors.
// U is n x n with U^T U = I, stored [a*n+p]:
//   forward  (integrals): X' = U^T X U, contracting U's first index
//   backward (densities): X' = U   X U^T, contracting U's second index
// Buffers are flat row-major in NOPT's CAS layout: 1-body [t*n+u], 2-body [((t*n+u)*n+v)*n+w].

// Two-index transform of an n x n matrix. X and out must not alias.
void rotate1(const double* X, const double* U, int n, double* out, bool forward);

// Four-index transform of G[((a*n+b)*n+c)*n+d]: four quarter-transforms, O(n^5). out doubles as
// scratch, so G and out must not alias.
void rotate2(const double* G, const double* U, int n, double* out, bool forward);

// Six-index transform of G[(((((a*n+b)*n+c)*n+d)*n+e)*n+f]: six quarter-transforms, O(n^7). out
// doubles as scratch, so G and out must not alias.
void rotate3(const double* G, const double* U, int n, double* out, bool forward);
