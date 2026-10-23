#ifndef _BLAS_LINK_H
#define _BLAS_LINK_H

#include "common_vars.h"
#include "omp.h"
// NOPT_BLAS_CBLAS_ONLY (set for the block2 TUs by the Makefile): CBLAS + thread control only.
// block2 declares its own Fortran LAPACK prototypes, which clang rejects as conflicting with
// lapack.h's (GCC lets them pass); those TUs call no LAPACK themselves.
#ifdef _OPENBLAS
# ifdef _ACCELERATE
#  include "accelerate_compat.h"
# else
#  include "cblas.h"
#  ifndef NOPT_BLAS_CBLAS_ONLY
#   include "lapacke.h"
#  endif
# endif
# ifndef lapack_int
#  define lapack_int int
# endif
#endif
#ifdef _MKL
# include "mkl_cblas.h"
#define lapack_int              MKL_INT
# include "mkl_lapack.h"
# include "mkl_service.h"
#endif


inline void set_blas_non_par(){
#ifdef _OPENBLAS
    // int ntb = openblas_get_num_threads();
    openblas_set_num_threads(1);
#endif
#ifdef _MKL
    // int ntb = mkl_get_max_threads();
    mkl_set_num_threads(1);
#endif
    omp_set_num_threads(num_threads);

}


inline void set_blas_par(){
#ifdef _OPENBLAS
    openblas_set_num_threads(num_threads);
#endif
#ifdef _MKL
    mkl_set_num_threads(num_threads);
#endif

}

inline void nopt_par_dgemm(const CBLAS_ORDER Layout, const CBLAS_TRANSPOSE TransA,
                           const CBLAS_TRANSPOSE TransB, const lapack_int M, const lapack_int N,
                           const lapack_int K, const double alpha, const double *A,
                           const lapack_int lda, const double *B, const lapack_int ldb,
                           const double beta, double *C, const lapack_int ldc){
    
    if(M==0)
        return;
    if(N==0)
        return;
    if(K==0){
        for(int i=0;i<M;i++)
        for(int j=0;j<N;j++)
            C[i*ldc+j]=C[i*ldc+j]*beta;
        return;
    }
#ifdef _SELF_BLAS_PAR
    cblas_dgemm(Layout, TransA, TransB, 
                M, N, K, alpha,
                A, lda, 
                B, ldb, beta, 
                C,ldc);
#endif
    
#ifdef _NOPT_BLAS_PAR
    for(int i=0;i<num_threads+1;i++)
        nopt_parallel_index[i]=(i*M)/num_threads;
#ifdef _OPENBLAS
    int ntb = openblas_get_num_threads();
    openblas_set_num_threads(1);
#endif
#ifdef _MKL
    int ntb = mkl_get_max_threads();
    mkl_set_num_threads(1);
#endif
    omp_set_num_threads(num_threads);
#pragma omp parallel
    {
        int i = omp_get_thread_num();
//         fprintf(stderr,"%d\n",i);
        int incr_a=1;
        if(TransA==CblasNoTrans){
            incr_a=lda;
        }
        cblas_dgemm(Layout, TransA, TransB, 
                    nopt_parallel_index[i+1]-nopt_parallel_index[i], N, K, alpha,
                    A+(long)nopt_parallel_index[i]*incr_a, lda, 
                    B, ldb, beta,
                    C+(long)nopt_parallel_index[i]*ldc,ldc);
    }
#ifdef _OPENBLAS
    openblas_set_num_threads(ntb);
#endif
#ifdef _MKL
    mkl_set_num_threads(ntb);
#endif

    
#endif
    
    
    
}

#endif
