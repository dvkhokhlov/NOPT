# include <stdio.h>
# include <stdlib.h>
// #include "lapacke.h"
#include <math.h>

#ifdef _OPENBLAS
# include "cblas.h"
# include "lapacke.h"
#endif
#ifdef _MKL
# include "mkl_cblas.h"
#define lapack_int              MKL_INT
# include "mkl_lapack.h"
# include "mkl_service.h"
#endif
# include "complex.h"

// int PrintCMatr(NOPT_complex * M, int n1, int n2, int g){
    // for(int i=0;i<n1;i++){
    // for(int j=0;j<n2;j++)printf(" % .4e+% .4e*i ",M[i*n2+j]);
        // printf("\n");
        // }
    // if (g==1)getchar();
    // return 0;
// }
#ifdef _OPENBLAS
#define NOPT_complex double _Complex
#endif
#ifdef _MKL
#define NOPT_complex MKL_Complex16
#endif



int lapack_herm_diag(double* Mr, double* Mi, double * eval, int n){
    
    
    NOPT_complex * M;
    M = (NOPT_complex *) malloc ( n*n * sizeof ( NOPT_complex));

#ifdef _OPENBLAS
    for(int i=0;i<n*n;i++)M[i]=Mr[i]+Mi[i]*I;
#endif
#ifdef _MKL
    for(int i=0;i<n*n;i++){M[i].real=Mr[i];M[i].imag=Mi[i];}
#endif
    
    
//     printf("H_ci:\n");
//     PrintCMatr(M,n,n,0);
            
    NOPT_complex * h_work;

    lapack_int lwork;
    lapack_int * iwork;
    lapack_int liwork;
    double * rwork;
    lapack_int lrwork;

    lapack_int info;

    double aux_work [1];
    lapack_int aux_iwork [1];

    //memory query
    lwork = n*n+2*n;//2*n*n+6*n+1;
    liwork = 5*n+3;
    lrwork = 2*n*n+5*n+1;
    h_work = (NOPT_complex *) malloc ( lwork * sizeof ( NOPT_complex));
    iwork = (lapack_int *) malloc ( liwork * sizeof ( lapack_int ));
    rwork = (double *) malloc ( lrwork * sizeof ( double ));
    
    if(h_work==NULL){
        if(lwork>1024*1024*1024)printf("ERROR: eigensystem() could not allocate %.0f GB of memory\n\n",lwork/1024/1024/1024);
        else                      printf("ERROR: eigensystem() could not allocate %.0f mB of memory\n\n",lwork/1024/1024);
        exit(EXIT_FAILURE);
    }

#ifdef _OPENBLAS
    LAPACK_zheevd("V","L",&n,M,&n,eval,h_work,&lwork,rwork,&lrwork,iwork,&liwork,&info);
    for(int i=0;i<n*n;i++){
        Mr[i]=creal(M[i]);
        Mi[i]=cimag(M[i]);
    }
#endif
#ifdef _MKL
//     printf("not supported\n");
//     exit(0);
    ZHEEVD("V","L",&n,M,&n,eval,h_work,&lwork,rwork,&lrwork,iwork,&liwork,&info);
    for(int i=0;i<n*n;i++){
        Mr[i]=M[i].real;
        Mi[i]=M[i].imag;
    }
#endif
    
        
    free(h_work);
    free(iwork);
    free(rwork);
    free(M);
    
    return 0;
}
 
