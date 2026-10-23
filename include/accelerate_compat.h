#ifndef _ACCELERATE_COMPAT_H
#define _ACCELERATE_COMPAT_H

// Apple Accelerate standing in for OpenBLAS (BLAS_LIB:=accelerate defines _OPENBLAS and _ACCELERATE).
// The CBLAS calls are identical; this maps the remaining OpenBLAS / netlib-lapack.h spellings that
// the _OPENBLAS code paths use onto Accelerate. The new LAPACK interface (macOS 13.3+) takes no
// hidden Fortran string lengths and uses __LAPACK_int (int unless ACCELERATE_LAPACK_ILP64), which
// matches the int arguments NOPT passes.

#define ACCELERATE_NEW_LAPACK
#include <omp.h>

// The vecLib pieces only, never the <Accelerate/Accelerate.h> umbrella: that drags in MacTypes.h,
// whose `wide`, `Point`, ... typedefs collide with NOPT globals (keywords.h `wide`). The rename
// guards against MacTypes.h arriving through lapack.h anyway.
#define wide nopt_mactypes_wide
#ifdef NOPT_BLAS_CBLAS_ONLY
// block2 TUs: CBLAS + threading only (see blas_link.h) -- lapack.h would clash with block2's own.
# include <vecLib/cblas_new.h>
# include <vecLib/thread_api.h>
#else
# include <vecLib/cblas_new.h>
# include <vecLib/lapack.h>
# include <vecLib/thread_api.h>
# define lapack_int __LAPACK_int
# define LAPACK_dsyevd dsyevd_
# define LAPACK_dgetrf dgetrf_
# define LAPACK_dgetri dgetri_
# define LAPACK_dgeev  dgeev_
# define LAPACK_dgesvd dgesvd_
# define LAPACK_zheevd zheevd_
#endif
#undef wide

// Accelerate threads its BLAS internally; NOPT asks for one thread inside its own OpenMP loops
// (BLAS_PAR=nopt) and restores afterwards. BLASSetThreading/BLASGetThreading: macOS 15 SDK.
// static inline: this header is also included from C (complex_diag.c).
static inline void openblas_set_num_threads(int n) {
    BLASSetThreading(n == 1 ? BLAS_THREADING_SINGLE_THREADED : BLAS_THREADING_MULTI_THREADED);
}
static inline int openblas_get_num_threads(void) {
    return BLASGetThreading() == BLAS_THREADING_SINGLE_THREADED ? 1 : omp_get_max_threads();
}

#endif
