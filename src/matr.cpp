# include <stdio.h>
# include <stdlib.h>
// #include "lapacke.h"
#include <math.h>
#include <cstring>
# include "matr.h"

# include "blas_link.h"

# define max(a,b)  (((a)<(b))?(b):(a))



int PrintMatr(double * M, int n1, int n2, int g){
    for(int i=0;i<n1;i++){
    for(int j=0;j<n2;j++)fprintf(out_stream," % .4e ",M[i*n2+j]);
        fprintf(out_stream,"\n");
        }
    if (g==1)getchar();
    return 0;
}

int fPrintMatr(FILE * stream, double * M, int n1, int n2, int g){
    for(int i=0;i<n1;i++){
    for(int j=0;j<n2;j++)fprintf(stream," % .4e ",M[i*n2+j]);
        fprintf(stream,"\n");
        }
    if (g==1)getchar();
    return 0;
}

int PrintMatr10(double * M, int n1, int n2, int g){
    for(int i=0;i<n1;i++){
    for(int j=0;j<n2;j++)fprintf(out_stream," % .10e ",M[i*n2+j]);
        fprintf(out_stream,"\n");
        }
    if (g)getchar();
    return 0;
}

int fPrintMatr10(FILE * stream, double * M, int n1, int n2, int g){
    
//     FILE * out =fopen(out_name,"w");
    
//     if(out==NULL){
//         printf("ERROR: fPrintMatr10 could not open %s\n",out_name);
//         exit(0);
//     }
    
    for(int i=0;i<n1;i++){
    for(int j=0;j<n2;j++)fprintf(stream," % .10e ",M[i*n2+j]);
        fprintf(stream,"\n");
    }    
//     fclose(out);
    if (g)getchar();
    
    return 0;
}


int PrintMatr_int(int * M, int n1, int n2, int g){
    for(int i=0;i<n1;i++){
    for(int j=0;j<n2;j++)printf(" %d ",M[i*n2+j]);
        printf("\n");
        }
    if (g)getchar();
    return 0;
}





// int copy_with_realloc(double * Out, double * Inp, int dim){
//     
//     if(Out!=nullptr) delete[] Out;
//     
//     Out= new double[dim];
//     
//     memcpy(Out, Inp, sizeof(double)*dim);
//     
//     return 0;
// }
    



int transpose(double * M, int n1, int n2){
    double c;
    for(int i=0;i<n1;i++){
    for(int j=i+1;j<n2;j++){
        c=M[i*n2+j];
        M[i*n2+j]=M[j*n2+i];
        M[j*n2+i]=c;
    }
    }
    return 0;
}

int transpose_A_to_B(double * Out, double * In, int n1, int n2, int ld1, int ld2){
    
    for(int i=0;i<n1;i++)
    for(int j=0;j<n2;j++)
        Out[i*ld2+j]=In[j*ld1+i];
        
    return 0;
}


int transpose_3d_abc_to_bac(double * Out, double * In, int na, int nb, int nc){
    
    for(int j=0;j<nb;j++)
    for(int i=0;i<na;i++)
    for(int k=0;k<nc;k++)
        Out[(j*na+i)*nc+k]=In[(i*nb+j)*nc+k];
    
    return 0;
}

int transpose_3d_abc_to_bca(double * Out, double * In, int na, int nb, int nc){
    
    for(int j=0;j<nb;j++)
    for(int k=0;k<nc;k++)
    for(int i=0;i<na;i++)
        Out[(j*nc+k)*na+i]=In[(i*nb+j)*nc+k];
    
    return 0;
}

int transpose_3d_abc_to_acb(double * Out, double * In, int na, int nb, int nc){
    
    for(int i=0;i<na;i++)
    for(int k=0;k<nc;k++)
    for(int j=0;j<nb;j++)
        Out[(i*nc+k)*nb+j]=In[(i*nb+j)*nc+k];
    
    return 0;
}

int transpose_3d_abc_to_acb_part_b_sum(double * Out, double * In, int na, int nb0, int nbf, int nb, int nc){
    
#pragma omp parallel for
    for(int i=0  ;i<na ;i++)
    for(int k=0  ;k<nc ;k++)
    for(int j=nb0;j<nbf;j++)
        Out[(i*nc+k)*nb+j]+=In[(i*nb+j)*nc+k];
    
    return 0;
}

int transpose_3d_abc_to_cab_sum(double * Out, double * In, int na, int nb, int nc){
    
    for(int k=0;k<nc;k++)
    for(int i=0;i<na;i++)
    for(int j=0;j<nb;j++)
        Out[(k*na+i)*nb+j]+=In[(i*nb+j)*nc+k];
    
    return 0;
}

int transpose_3d_abc_to_cab_part_b_sum(double * Out, double * In, int na, int nb0, int nbf, int nb, int nc){
    
#pragma omp parallel for
    for(int k=0  ;k<nc ;k++)
    for(int i=0  ;i<na ;i++)
    for(int j=nb0;j<nbf;j++)
        Out[(k*na+i)*nb+j]+=In[(i*nb+j)*nc+k];
    
    return 0;
}


int transpose_3d_abc_to_cab_sum_with_b_mult(double * Out, double * In, int na, int nb, int nc, double m){
    
    for(int k=0;k<nc;k++)
    for(int i=0;i<na;i++)
    for(int j=0;j<nb;j++)
        Out[(k*na+i)*nb+j]+=In[(i*nb+j)*nc+k]*m;
    
    return 0;
}

int transpose_3d_abc_to_abc_part_c_sum(double * Out, double * In, int na, int nb, int nc0, int ncf, int nc){
    
#pragma omp parallel for    
    for(int i=0  ;i<na*nb;i++)
    for(int j=nc0;j<ncf  ;j++)
        Out[i*nc+j]+=In[i*nc+j];
    
    return 0;
}


int transpose_3d_abc_to_bac_int(int * Out, int * In, int na, int nb, int nc){
    
    for(int j=0;j<nb;j++)
    for(int i=0;i<na;i++)
    for(int k=0;k<nc;k++)
        Out[(j*na+i)*nc+k]=In[(i*nb+j)*nc+k];
    
    return 0;
}


int lapack_diag(double * M, double * eval, int n){
    
    if(n==0){
        return 0;
    }
    
	double * h_work;

	lapack_int lwork;
	lapack_int * iwork;
	lapack_int liwork;
	lapack_int info;

	double aux_work [1];
	lapack_int aux_iwork [1];

    //memory query
    lwork = 2*n*n+6*n+1;
    liwork = 5*n+3;
    iwork = (lapack_int *) malloc ( liwork * sizeof ( lapack_int ));
    h_work = (double *) malloc ( lwork * sizeof ( double ));
    
    if(h_work==NULL){
        if(lwork>1024*1024*1024)printf("ERROR: eigensystem() could not allocate %.0f GB of memory\n\n",lwork/1024/1024/1024);
        else                      printf("ERROR: eigensystem() could not allocate %.0f mB of memory\n\n",lwork/1024/1024);
        exit(EXIT_FAILURE);
    }

#ifdef _OPENBLAS
    // printf("ogp = %d\n",openblas_get_parallel());
    // printf("os  = %d\n",OPENBLAS_SEQUENTIAL);
    // printf("oth = %d\n",OPENBLAS_THREAD);
    // printf("ont = %d\n",openblas_get_num_threads());
    // printf("omp = %d\n",OPENBLAS_OPENMP);
    // printf("mth = %d\n",omp_get_max_threads());
    LAPACK_dsyevd("V","L",&n,M,&n,eval,h_work,&lwork,iwork,&liwork,&info);
#endif
#ifdef _MKL
    DSYEVD("V","L",&n,M,&n,eval,h_work,&lwork,iwork,&liwork,&info);
#endif
        
    free(iwork);
    free(h_work);
    
    return 0;
}


double * comp_values;

int compare(const void * i, const void * j){
    
    return (comp_values[*(int*)i]>comp_values[*(int*)j]);

}

int sorter(const void *a, const void *b)
{
    return *(int*)a-*(int*)b;
}


int sort_int_array(int * A, int n){
    
    qsort(A,n,sizeof ( int ),sorter);
    
    return 0;
}


int lapack_right_eig(double * M, double * eval, int n){
    
    
    int info ;
    int ione = 1;
    double * h_work ;
    int lwork ;
//     int nb = magma_get_dgesvd_nb (n,n); // opt . block size
    lwork = 2*n*n+6*n+1;
//     magma_dmalloc_pinned(&h_work,lwork); // host mem . for h_work
    h_work = (double *) malloc ( lwork * sizeof ( double ));
    double * R = (double *) malloc ( n*n * sizeof ( double ));
    double * L = (double *) malloc ( n*n * sizeof ( double ));
    double * ei= (double *) malloc (   n * sizeof ( double ));
    comp_values = (double *) malloc (   n * sizeof ( double ));
    int * int_values = (int *) malloc (   n * sizeof ( int ));
#ifdef _OPENBLAS
    LAPACK_dgeev("V","V",&n,M,&n,eval,ei,L,&n,R,&n,h_work,&lwork,&info );
#endif    
#ifdef _MKL
    DGEEV("V","V",&n,M,&n,eval,ei,L,&n,R,&n,h_work,&lwork,&info );
#endif
//     LAPACK_dgeev(   ,  r, n, ,a,   wr, i,vl,l, , r,     k,* lwor, lapack_int *info );
//     for(int i=0; i<n;i++)printf("ei = %e\n",eval[i]);
    for(int i=0;i<n;i++)comp_values[i]=eval[i];
    for(int i=0;i<n;i++) int_values[i]=i;
    qsort(int_values,n,sizeof ( int ),compare);
    for(int i=0;i<n;i++)eval[i]=comp_values[int_values[i]];
    
    for(int i=0; i<n;i++)if(ei[int_values[i]]>1E-15)printf("WARNING: imaginary eigenvalue ev[%d]=%e+i*%e in dgeev\n",i,eval[i],ei[int_values[i]]);
//     getchar();
    
    for(int i=0;i<n;i++)
    for(int j=0;j<n;j++)M[i*n+j]=R[int_values[i]*n+j];    
          
    free(h_work);
    free(R);
    free(L);
    free(ei);
    free(comp_values);
    free(int_values);

    return 0;
}

int tr_and_diag(double * H, double * M, double * B, double * V, double * E, int dim){
    
    //H = M^t*H*M
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                        dim,dim,dim,1.0,
                        H,dim,
                        M,dim,0.0,
                        V,dim);
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                        dim,dim,dim,1.0,
                        M,dim,
                        V,dim,0.0,
                        B,dim);
    
    //diag(H) 
    lapack_diag(B,E,dim);
    
    //V = M*V
    //V^t=V^t*M - because of row-matrix
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                        dim,dim,dim,1.0,
                        B,dim,
                        M,dim,0.0,
                        V,dim);
    
    
    return 0;
}


int S05_calc(double * S, double * M, double * P, int n){
    
    double * S_eval;
    S_eval = new double[n];
    lapack_diag(S, S_eval, n);

    double * tmp;
    tmp = new double[n*n];
    for(int i = 0;i<n;i++) for(int j=0;j<n;j++) {M[i*n+j]=0.0;P[i*n+j]=0.0;}
    
//     for(int i = 0;i<n;i++)printf("S_eval[%d] = %e\n",i,S_eval[i]);
//     getchar();
    for(int i = 0;i<n;i++) 
        if(S_eval[i]>0.0) {P[i*(n+1)]=sqrt(S_eval[i]);M[i*(n+1)]=1.0/P[i*(n+1)];}
        else fprintf(out_stream, "WARNING: negetive eigenvalue of S matix: %e\n",S_eval[i]);
    
        
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,n,n,n,1.0,P,n,S  ,n,0.0,tmp,n); // P*S^-1
    cblas_dgemm(CblasRowMajor,  CblasTrans,CblasNoTrans,n,n,n,1.0,S,n,tmp,n,0.0,P  ,n); // S*P*S^-1

    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,n,n,n,1.0,M,n,S  ,n,0.0,tmp,n); // M*S^-1
    cblas_dgemm(CblasRowMajor,  CblasTrans,CblasNoTrans,n,n,n,1.0,S,n,tmp,n,0.0,M  ,n); // S*M*S^-1

    
//     mat_mult(gpu_Sm05, gpu_H, gpu_tmp_matrix, n_ao, 0, 1, 1.0, 0.0);
//     mat_mult(gpu_H, gpu_tmp_matrix, gpu_Sm05, n_ao, 0, 0, 1.0, 0.0);
    delete[] S_eval;
    delete []tmp;
    return 0;
}

int S05_calc_p(double * S, double * M, double * P, int n, double eps){
    
    double * S_eval;
    S_eval = new double[n];
    lapack_diag(S, S_eval, n);

    double * tmp;
    tmp = new double[n*n];
    for(int i = 0;i<n;i++) for(int j=0;j<n;j++) {M[i*n+j]=0.0;P[i*n+j]=0.0;}
    
    int non_zero=0;
    for(int i = 0;i<n;i++) {
        if(S_eval[i]>eps){
            P[i*(n+1)]=sqrt(S_eval[i]);
            M[i*(n+1)]=1.0/P[i*(n+1)];
            non_zero++;
        }
    }
        
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,n,n,n,1.0,P,n,S  ,n,0.0,tmp,n); // P*S^-1
    cblas_dgemm(CblasRowMajor,  CblasTrans,CblasNoTrans,n,n,n,1.0,S,n,tmp,n,0.0,P  ,n); // S*P*S^-1

    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,n,n,n,1.0,M,n,S  ,n,0.0,tmp,n); // M*S^-1
    cblas_dgemm(CblasRowMajor,  CblasTrans,CblasNoTrans,n,n,n,1.0,S,n,tmp,n,0.0,M  ,n); // S*M*S^-1

    
//     mat_mult(gpu_Sm05, gpu_H, gpu_tmp_matrix, n_ao, 0, 1, 1.0, 0.0);
//     mat_mult(gpu_H, gpu_tmp_matrix, gpu_Sm05, n_ao, 0, 0, 1.0, 0.0);
    delete[] S_eval;
    delete []tmp;
    return non_zero;
}


int lapack_svd(double * mat, double * L, double * R, double * svd, int n){
    
    int info ;
    int ione = 1;
    double * h_work ;
    int lwork ;
//     int nb = magma_get_dgesvd_nb (n,n); // opt . block size
    lwork = 2*n*n+6*n+1;
//     magma_dmalloc_pinned(&h_work,lwork); // host mem . for h_work
    h_work = (double *) malloc ( lwork * sizeof ( double ));

#ifdef _OPENBLAS
    LAPACK_dgesvd("A","A",&n,&n,mat,&n,svd,L,&n,R,&n,h_work,&lwork,&info );
#endif
#ifdef _MKL
    DGESVD("A","A",&n,&n,mat,&n,svd,L,&n,R,&n,h_work,&lwork,&info );
#endif
            
    free(h_work);


    return 0;
    
    
}

int adj_constr_by_svd(double * O, double * L, double * R, double * S, double * B, int n, double sc){
    
    set_zero_matr(B, n*n);
    
//     PrintMatr(L,n,n,1);
//     
//     PrintMatr(R,n,n,1);
//     
//     PrintMatr(S,n,1,1);
//     
//     PrintMatr(B,n,n,1);
    
    double pl;
    
    for(int i=0; i<n; i++){
        pl=sc;
        for(int j=0; j<n; j++)if(j!=i)pl=pl*S[j];
//         printf("pl= %e\n",pl);
        cblas_daxpy(n, pl, L+i*n,1,B+i*n,1);
//         for(int j=0; j<n; j++)
//                 B[i*n+j]+=L[i*n+j]*S[i];
        
//         for(int j=0; j<n; j++)
//             for(int k=0; k<n; k++)
//                 O[i*n+j]+=L[k*n+j]*S[k]*R[i*n+k];
    }
    
//     printf("new\n");
//     
//     PrintMatr(L,n,n,1);
//     
//     PrintMatr(R,n,n,1);
//     
//     PrintMatr(S,n,1,1);
//     
//     PrintMatr(B,n,n,1);
      
    
    
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,n,n,n,1.0,R,n,B,n,0.0,O,n);
    
    return 0;
    
}


int adj2_constr_by_svd(double * O, const double * __restrict__ L, const double * __restrict__ R, double * S, double *B, int n, double sc){
    
    set_zero_matr(O, n*n*n*n);
    
    for(int a=0; a<n; a++)
    for(int b=0; b<n; b++){
        B[a*n+b]=sc;
        for(int c=0; c<n; c++)if(c!=a)if(c!=b)B[a*n+b]=B[a*n+b]*S[c];
    }
    
    double * T1 = B+n*n*n*n;
    
    for(int k=0; k<n; k++)
    for(int l=0; l<n; l++)
    for(int a=0; a<n; a++)
    for(int b=0; b<n; b++)
    {
        T1[((a*n+b)*n+k)*n+l]=B[a*n+b]*(L[a*n+k]*L[b*n+l]-L[a*n+l]*L[b*n+k]);
    }
    
    double * T2 = T1+n*n*n*n;
    set_zero_matr(T2,n*n*n*n);
    
    for(int i=0; i<n; i++)
    for(int k=0; k<n; k++)
    for(int l=0; l<n; l++)
    for(int a=0; a<n; a++)
    for(int b=0; b<n; b++)
    {
        T2[((i*n+b)*n+k)*n+l]+=R[i*n+a]*T1[((a*n+b)*n+k)*n+l];
    }
    
    
    for(int i=0; i<n; i++)
    for(int j=0; j<n; j++)
    for(int k=0; k<n; k++)
    for(int a=0; a<n; a++)
    for(int l=0; l<n; l++)
    for(int b=0; b<n; b++)
    {
//         O[((i*n+j)*n+k)*n+l]+=R[j*n+b]*T1[((i*n+b)*n+k)*n+l];
        O[((i*n+j)*n+k)*n+l]+=B[a*n+b]*R[i*n+a]*R[j*n+b]*(L[a*n+k]*L[b*n+l]-L[a*n+l]*L[b*n+k]);//correct but slow way
//         O[((i*n+j)*n+k)*n+l]+=pl*R[i*n+a]*R[i*n+a]*(L[k*n+a]*L[l*n+b]-L[l*n+a]*L[l*n+b]);//local test
    }
    
//     delete[] T1;
//     delete[] T2;
    
    return 0;
    
}


int matr_constr_by_svd(double * O, double * L, double * R, double * S, double * B, int n){
    
    set_zero_matr(B, n*n);
    
    for(int i=0; i<n; i++){
        cblas_daxpy(n, S[i], L+i*n,1,B+i*n,1);
    }
    
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,n,n,n,1.0,R,n,B,n,0.0,O,n);
    
    return 0;
    
}

int zero_svd_ortogonalization(double * L, double * R, double * S, double * B, double * B2, int n){
    
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,n,n,n,1.0,R,n,L,n,0.0,B,n);
    
    inv_matr_constr(B,n);
    
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,n,n,n,1.0,B,n,R,n,0.0,B2,n);
    
    memcpy(R,B2,n*n*sizeof(double));
    
//     for(int i=0; i<n; i++)
//         if(fabs(S[i])<eps) for(int j=0; j<n; j++) R[i*n+j]=L[i*n+j];
    return 0;
    
}



double minor_det(double* MATR, int dim, int a, int b){
    
    double * minor = new double[(dim-1)*(dim-1)];
    double result;
    int k,l;
    
    for(k=0;k<a;k++)for(l=0;l<b;l++)minor[k*(dim-1)+l]=MATR[k*dim+l];
    for(k=a+1;k<dim;k++)for(l=0;l<b;l++)minor[(k-1)*(dim-1)+l]=MATR[k*dim+l];
    for(k=0;k<a;k++)for(l=b+1;l<dim;l++)minor[k*(dim-1)+l-1]=MATR[k*dim+l];
    for(k=a+1;k<dim;k++)for(l=b+1;l<dim;l++)minor[(k-1)*(dim-1)+l-1]=MATR[k*dim+l];
//     PrintMatr(minor,dim-1,dim-1,0);			
    result =mat_det_calc_lapack(minor,dim-1);
    
//     for(k=0;k<a;k++)for(l=0;l<b;l++)minor[k*(dim-1)+l]=MATR[k*dim+l];
//     for(k=a+1;k<dim;k++)for(l=0;l<b;l++)minor[(k-1)*(dim-1)+l]=MATR[k*dim+l];
//     for(k=0;k<a;k++)for(l=b+1;l<dim;l++)minor[k*(dim-1)+l-1]=MATR[k*dim+l];
//     for(k=a+1;k<dim;k++)for(l=b+1;l<dim;l++)minor[(k-1)*(dim-1)+l-1]=MATR[k*dim+l];
    
//     printf("%e %e\n",result,mat_det_calc(minor,dim-1));
    
    delete[] minor;
    return result;  
}

double adj_calc(double* MATR, int dim, int a, int b){
    
    if(dim==1)return 1.0;
    if(dim<1)return 0.0;
    
    double * minor = new double[(dim-1)*(dim-1)];
    double result;
    int k,l;
    
    for(k=0;  k<a;  k++)for(l=0;  l<b;  l++)minor[k*    (dim-1)+l] =MATR[k*dim+l];
    for(k=a+1;k<dim;k++)for(l=0;  l<b;  l++)minor[(k-1)*(dim-1)+l] =MATR[k*dim+l];
    for(k=0;  k<a;  k++)for(l=b+1;l<dim;l++)minor[k*    (dim-1)+l-1]=MATR[k*dim+l];
    for(k=a+1;k<dim;k++)for(l=b+1;l<dim;l++)minor[(k-1)*(dim-1)+l-1]=MATR[k*dim+l];
//     PrintMatr(minor,dim-1,dim-1,0);			
    result =mat_det_calc_lapack(minor,dim-1);
    
//     for(k=0;k<a;k++)for(l=0;l<b;l++)minor[k*(dim-1)+l]=MATR[k*dim+l];
//     for(k=a+1;k<dim;k++)for(l=0;l<b;l++)minor[(k-1)*(dim-1)+l]=MATR[k*dim+l];
//     for(k=0;k<a;k++)for(l=b+1;l<dim;l++)minor[k*(dim-1)+l-1]=MATR[k*dim+l];
//     for(k=a+1;k<dim;k++)for(l=b+1;l<dim;l++)minor[(k-1)*(dim-1)+l-1]=MATR[k*dim+l];
    
//     printf("%e %e\n",result,mat_det_calc(minor,dim-1));
    

    if((a+b)%2)result=-result;
    
    delete[] minor;
    return result;  
}


double double_minor_det(double* MATR, int dim, int a, int b, int c, int d){
    
    double * minor = new double[(dim-1)*(dim-1)];
    double result;
    int k,l;
    
    for(k=0;k<a;k++)for(l=0;l<b;l++)minor[k*(dim-1)+l]=MATR[k*dim+l];
    for(k=a+1;k<dim;k++)for(l=0;l<b;l++)minor[(k-1)*(dim-1)+l]=MATR[k*dim+l];
    for(k=0;k<a;k++)for(l=b+1;l<dim;l++)minor[k*(dim-1)+l-1]=MATR[k*dim+l];
    for(k=a+1;k<dim;k++)for(l=b+1;l<dim;l++)minor[(k-1)*(dim-1)+l-1]=MATR[k*dim+l];
//     PrintMatr(minor,dim-1,dim-1,0);
    if(c==a){printf("ERROR in double_minor_det a=c=%d\n",a);exit(1);}
    if(b==d){printf("ERROR in double_minor_det b=d=%d\n",b);exit(1);}
    if(c>a)c--;
    if(d>b)d--;
    
    result = minor_det(minor, dim-1, c, d);
//     getchar();
//     for(k=0;k<a;k++)for(l=0;l<b;l++)minor[k*(dim-1)+l]=MATR[k*dim+l];
//     for(k=a+1;k<dim;k++)for(l=0;l<b;l++)minor[(k-1)*(dim-1)+l]=MATR[k*dim+l];
//     for(k=0;k<a;k++)for(l=b+1;l<dim;l++)minor[k*(dim-1)+l-1]=MATR[k*dim+l];
//     for(k=a+1;k<dim;k++)for(l=b+1;l<dim;l++)minor[(k-1)*(dim-1)+l-1]=MATR[k*dim+l];
    
//     printf("%e %e\n",result,mat_det_calc(minor,dim-1));
    
    delete[] minor;
    return result;  
}

double double_adj_calc(double* MATR, int dim, int a, int b, int c, int d){
    
//     if(dim==2)return 1.0;
    if(dim<2)return 0.0;
        
    
    double * minor = new double[(dim-1)*(dim-1)];
    double result;
    int k,l;
    
    for(k=0;k<a;k++)for(l=0;l<b;l++)minor[k*(dim-1)+l]=MATR[k*dim+l];
    for(k=a+1;k<dim;k++)for(l=0;l<b;l++)minor[(k-1)*(dim-1)+l]=MATR[k*dim+l];
    for(k=0;k<a;k++)for(l=b+1;l<dim;l++)minor[k*(dim-1)+l-1]=MATR[k*dim+l];
    for(k=a+1;k<dim;k++)for(l=b+1;l<dim;l++)minor[(k-1)*(dim-1)+l-1]=MATR[k*dim+l];
//     PrintMatr(minor,dim-1,dim-1,0);
    if(c==a){printf("ERROR in double_adj_calc a=c=%d\n",a);exit(1);}
    if(b==d){printf("ERROR in double_adj_calc b=d=%d\n",b);exit(1);}
    if(c>a)c--;
    if(d>b)d--;
    
    result = adj_calc(minor, dim-1, c, d);
//     getchar();
//     for(k=0;k<a;k++)for(l=0;l<b;l++)minor[k*(dim-1)+l]=MATR[k*dim+l];
//     for(k=a+1;k<dim;k++)for(l=0;l<b;l++)minor[(k-1)*(dim-1)+l]=MATR[k*dim+l];
//     for(k=0;k<a;k++)for(l=b+1;l<dim;l++)minor[k*(dim-1)+l-1]=MATR[k*dim+l];
//     for(k=a+1;k<dim;k++)for(l=b+1;l<dim;l++)minor[(k-1)*(dim-1)+l-1]=MATR[k*dim+l];
    
//     printf("%e %e\n",result,mat_det_calc(minor,dim-1));
    if((a+b)%2)result=-result;
    
    delete[] minor;
    return result;  
}

int double_adj_mat_constr(double* MATR, int dim, int a, int b, double* OUT){
    
//     if(dim==2)return 1.0;
        if(dim==2){
            OUT[0]=1.0;
            return 0;
        }
        
    double result;
    int k,l;
    
//     for(k=0;k<a;k++)    for(l=0;l<b;l++)    OUT[k*(dim-1)+l]=MATR[k*dim+l];
//     for(k=a+1;k<dim;k++)for(l=0;l<b;l++)    OUT[(k-1)*(dim-1)+l]=MATR[k*dim+l];
//     for(k=0;k<a;k++)    for(l=b+1;l<dim;l++)OUT[k*(dim-1)+l-1]=MATR[k*dim+l];
//     for(k=a+1;k<dim;k++)for(l=b+1;l<dim;l++)OUT[(k-1)*(dim-1)+l-1]=MATR[k*dim+l];
//     #pragma omp parallel for                                                                
    for(k=0;k<a;k++)     memcpy(OUT+k*(dim-1),         MATR+k*dim  ,   b*sizeof(double));
//     #pragma omp parallel for
    for(k=0;k<a;k++)     memcpy(OUT+k*(dim-1)+b,     MATR+k*dim+b+1,   (dim-b-1)*sizeof(double));
//     #pragma omp parallel for 
    for(k=a+1;k<dim;k++) memcpy(OUT+(k-1)*(dim-1),     MATR+k*dim  ,   b*sizeof(double));
//     #pragma omp parallel for 
    for(k=a+1;k<dim;k++) memcpy(OUT+(k-1)*(dim-1)+b, MATR+k*dim+b+1,   (dim-b-1)*sizeof(double));
    
    
//     PrintMatr(OUT,dim-1,dim-1,0);
    
    
    
    int info = adj_matr_constr(OUT, dim-1);
    
//     getchar();
//     for(k=0;k<a;k++)for(l=0;l<b;l++)minor[k*(dim-1)+l]=MATR[k*dim+l];
//     for(k=a+1;k<dim;k++)for(l=0;l<b;l++)minor[(k-1)*(dim-1)+l]=MATR[k*dim+l];
//     for(k=0;k<a;k++)for(l=b+1;l<dim;l++)minor[k*(dim-1)+l-1]=MATR[k*dim+l];
//     for(k=a+1;k<dim;k++)for(l=b+1;l<dim;l++)minor[(k-1)*(dim-1)+l-1]=MATR[k*dim+l];
    
//     printf("%e %e\n",result,mat_det_calc(minor,dim-1));
//     if((a+b)%2)result=-result;
    
//     delete[] minor;
    return info;  
}

double adj_matr_constr(double * IN, int N){
    
    lapack_int * piv , info ; // piv - array of indices of inter -
    piv = new lapack_int[ N * sizeof ( lapack_int )];
#ifdef _OPENBLAS
    LAPACK_dgetrf(&N,&N,IN,&N,piv,&info);
#endif
#ifdef _MKL
    DGETRF(&N,&N,IN,&N,piv,&info);
#endif
    
    double det=IN[0];
    for(int i=1;i<N;i++){
    det  = det*IN[(N+1)*i];
//     printf("%e+i*%e\t",eval[0],eval1[0] );
    }
    for(int i=0;i<N-1;i++)if(piv[i]!=(i+1))det=-1.0*det;
    
    lapack_int lwork;
    double * work;
    
    lwork = 2*N*N+6*N+1;
    work = new double[lwork] ;
        
#ifdef _OPENBLAS
    LAPACK_dgetri(&N,IN,&N,piv,work,&lwork,&info);
#endif
#ifdef _MKL
    DGETRI(&N,IN,&N,piv,work,&lwork,&info);
#endif
    
    cblas_dscal(N*N,det,IN,1);
    
    
    
    delete[] work;
    delete[] piv;
    
    return det;
}


int pseudo_inv_matr_constr(double * IN, int n){
    
    double * S_eval;
    S_eval = new double[n];
    lapack_diag(IN, S_eval, n);

    double * tmp;
    tmp = new double[n*n];
    double * M;
    M = new double[n*n];
    for(int i = 0;i<n;i++) for(int j=0;j<n;j++) M[i*n+j]=0.0;
    
    int non_zero=0;
    for(int i = 0;i<n;i++) {
        if(S_eval[i]>1e-10){
            M[i*(n+1)]=1.0/S_eval[i];
            non_zero++;
        }
    }
        
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,n,n,n,1.0,M,n,IN  ,n,0.0,tmp,n); // M*S^-1
    cblas_dgemm(CblasRowMajor,  CblasTrans,CblasNoTrans,n,n,n,1.0,IN,n,tmp,n,0.0,M ,n); // S*M*S^-1
    memcpy(IN,M,sizeof(double)*n*n);

    
//     mat_mult(gpu_Sm05, gpu_H, gpu_tmp_matrix, n_ao, 0, 1, 1.0, 0.0);
//     mat_mult(gpu_H, gpu_tmp_matrix, gpu_Sm05, n_ao, 0, 0, 1.0, 0.0);
    delete[] S_eval;
    delete []tmp;
    delete[] M;
    return non_zero;    
    
    
    
}

int inv_matr_constr(double * IN, int N){
    
//     pseudo_inv_matr_constr(IN, N);
    
    lapack_int * piv , info ; // piv - array of indices of inter -
    piv = new lapack_int[ N * sizeof ( lapack_int )];
#ifdef _OPENBLAS
    LAPACK_dgetrf(&N,&N,IN,&N,piv,&info);
#endif
#ifdef _MKL
    DGETRF(&N,&N,IN,&N,piv,&info);
#endif

    lapack_int lwork;
    double * work;
    
    lwork = 2*N*N+6*N+1;
    work = new double[lwork] ;
        
#ifdef _OPENBLAS
    LAPACK_dgetri(&N,IN,&N,piv,work,&lwork,&info);
#endif
#ifdef _MKL
    DGETRI(&N,IN,&N,piv,work,&lwork,&info);
#endif
    
    
    delete[] work;
    delete[] piv;
    
    return 0;
}


double mat_det_calc_lapack_no_change(double * mat, int N){//LU
    
    if(N==0) return 1.0;
    double * buf = new double[N*N];
    
    for(int i=0;i<N*N;i++)buf[i]=mat[i];
    
    lapack_int * piv , info ; // piv - array of indices of inter -
    piv =( lapack_int *) malloc ( N * sizeof ( lapack_int ));
#ifdef _OPENBLAS
    LAPACK_dgetrf(&N,&N,buf,&N,piv,&info);
#endif
#ifdef _MKL
    DGETRF(&N,&N,buf,&N,piv,&info);
#endif

//     printf("piv:");
//     for(int i=0;i<N;i++)printf(" %d",piv[i]);
//                         printf("\n");
    double result=buf[0];
    for(int i=1;i<N;i++){
    result  = result*buf[(N+1)*i];
//     printf("%e+i*%e\t",eval[0],eval1[0] );
    }
    for(int i=0;i<N-1;i++)if(piv[i]!=(i+1))result=-1.0*result;
    free(piv);
    delete[] buf;
    return result;
}

double mat_det_calc_lapack(double * mat, int N){//LU
    if(N==0) return 1.0;
//     for(int i =0;i<N;i++){
//     for(int j=0;j<N;j++)printf("%e  ",mat[i*N+j]);
//     printf("\n");}

    lapack_int * piv , info ; // piv - array of indices of inter -
    piv =( lapack_int *) malloc ( N * sizeof ( lapack_int ));
#ifdef _OPENBLAS
    LAPACK_dgetrf(&N,&N,mat,&N,piv,&info);
#endif
#ifdef _MKL
    DGETRF(&N,&N,mat,&N,piv,&info);
#endif

//     printf("piv:");
//     for(int i=0;i<N;i++)printf(" %d",piv[i]);
//                         printf("\n");
    double result=mat[0];
    for(int i=1;i<N;i++){
    result  = result*mat[(N+1)*i];
//     printf("%e+i*%e\t",eval[0],eval1[0] );
    }
    for(int i=0;i<N-1;i++)if(piv[i]!=(i+1))result=-1.0*result;
    free(piv);
    return result;
}


double qtr_mat_det_calc(double * mat, int dim, int incr){
    double tmp_mat[dim*dim];
    int i,j;
    double c;
    double result=1.0;
    for(i=0;i<dim;i++)
    for(j=0;j<dim;j++)tmp_mat[i*dim+j]=mat[i*incr+j];
    
//     PrintMatr(tmp_mat,dim,dim,1);
    int f_el_row[dim];
    for(i=0;i<dim;i++)f_el_row[i]=i*(dim+1);
    for(i=0;i<dim-1;i++){
          if(fabs(tmp_mat[f_el_row[i]])<fabs(tmp_mat[f_el_row[i+1]-1])){
             result=-result;
             j=f_el_row[i];
             f_el_row[i]=f_el_row[i+1]-1;
             f_el_row[i+1]=j+1;
             }
          result=result*tmp_mat[f_el_row[i]];
          if(result==0)return 0;
          c=tmp_mat[f_el_row[i+1]-1]/tmp_mat[f_el_row[i]];
          for(j=0;j<dim-i;j++)tmp_mat[f_el_row[i+1]-1+j]-=c*tmp_mat[f_el_row[i]+j];
//           PrintMatr(tmp_mat,dim,dim,1);
    }
    result=result*tmp_mat[f_el_row[dim-1]];
//     printf("det = %e\n",result);
    return result;
}

int triang_minor_matr_constr(double * mat, double * inv_mat,int N){
    
    int i,j,d_i;
    
    for(i=0;i<N;i++){
    for(j=0;j<N;j++)inv_mat[i*N+j]=15.0;
        }
    
    for(i=0;i<N;i++){
    for(j=i+1;j<N;j++)inv_mat[i*N+j]=0.0;
        }
    
    for(i=0;i<N;i++){
        inv_mat[i*(N+1)]=1.0;
        for(j=0;j<i;j++)inv_mat[i*(N+1)]=inv_mat[i*(N+1)]*mat[j*(N+1)];
        for(j=i+1;j<N;j++)inv_mat[i*(N+1)]=inv_mat[i*(N+1)]*mat[j*(N+1)];
    }
    for(i=1;i<N;i++){
        inv_mat[i*(N+1)-1]=-1.0;
        for(j=0;j<i-1;j++)inv_mat[i*(N+1)-1]=inv_mat[i*(N+1)-1]*mat[j*(N+1)];
        inv_mat[i*(N+1)-1]=inv_mat[i*(N+1)-1]*mat[i*(N+1)-N];
        for(j=i+1;j<N;j++)inv_mat[i*(N+1)-1]=inv_mat[i*(N+1)-1]*mat[j*(N+1)];
    }
    for(d_i=2;d_i<N;d_i++){
        for(i=0;i<N-d_i;i++){//printf("\n\n");/*PrintMatr(mat,N,N,0);*/printf("(%d,%d)\n",i,d_i);
            inv_mat[i*(N+1)+d_i*N]=pow(-1.0,d_i);
//             PrintMatr(inv_mat,N,N,1);
            for(j=0;j<i;j++)inv_mat[i*(N+1)+d_i*N]=inv_mat[i*(N+1)+d_i*N]*mat[j*(N+1)];
//             PrintMatr(inv_mat,N,N,1);
            inv_mat[i*(N+1)+d_i*N]=inv_mat[i*(N+1)+d_i*N]*qtr_mat_det_calc(mat+i*(N+1)+1,d_i,N);
//             PrintMatr(inv_mat,N,N,1);
            for(j=i+d_i+1;j<N;j++)inv_mat[i*(N+1)+d_i*N]=inv_mat[i*(N+1)+d_i*N]*mat[j*(N+1)];
//             PrintMatr(inv_mat,N,N,1);
            
            
        }
        
        
    }
    

    return 0;
}


int ortogonalization(double * matr, int n){
    
    double a;
    int i,j;
    for(i=0;i<n;i++){
    a = cblas_ddot(n,matr+i*n,1,matr+i*n,1);
    a=1/sqrt(a);
    cblas_dscal(n,a,matr+i*n,1);
    for(j=i+1;j<n;j++){
    a = cblas_ddot(n,matr+i*n,1,matr+j*n,1);
    a=-a;
    cblas_daxpy(n, a,matr+i*n,1,matr+j*n,1);
    }
    }
    // Destroy
    
    return 0;
}   

int ortogonalization_no_sq(double * matr, int n, int m){
    
    double a;
    int i,j;
    for(i=0;i<n;i++){
        a = cblas_ddot(m,matr+i*m,1,matr+i*m,1);
//         PrintMatr(matr+i*m,1,m,0);
//         printf("a=%e\n",a);
//         getchar();
        a=1/sqrt(a);
//         if(matr[i*(m+1)]<0.0)a=-a;
        cblas_dscal(m,a,matr+i*m,1);
        for(j=i+1;j<n;j++){
            a = cblas_ddot(m,matr+i*m,1,matr+j*m,1);
//             printf("a=%e\n",a);
//             getchar();
            a=-a;
            cblas_daxpy(m, a,matr+i*m,1,matr+j*m,1);
        }
    }
   
    return 0;
}   

int ortogonalization_no_sq_row_with_base(double * matr, int n, int m, int ld){
    
    double a;
    int i,j;
    for(i=0;i<n;i++){
        a = cblas_ddot(m,matr+i,ld,matr+i,ld);
//         printf("%e\n",a);
        a=1/sqrt(a);
//         if(matr[i*(m+1)]<0.0)a=-a;
        cblas_dscal(m,a,matr+i,ld);
        for(j=i+1;j<n;j++){
            a = - cblas_ddot(m,matr+i,ld,matr+j,ld);
            cblas_daxpy(m, a,matr+i,ld,matr+j,ld);
        }
    }
//     getchar();
    return 0;
}   



int ortonormalization_no_sq(double * matr, int n, int m){
    
    double a;
    int i,j;
    for(i=0;i<n;i++){
    a = cblas_ddot(m,matr+i*m,1,matr+i*m,1);
    a=1/sqrt(a);
//     if(matr[i*(m+1)]<0.0)a=-a;
    cblas_dscal(m,a,matr+i*m,1);
//     for(j=i+1;j<n;j++){
//     a = cblas_ddot(m,matr+i*m,1,matr+j*m,1);
//     a=-a;
//     cblas_daxpy(m, a,matr+i*m,1,matr+j*m,1);
//     }
    }
    // Destroy
    
    return 0;
}   

int ortogonalization_no_sq_by_ref(double ** matr, int n, int m){
    
    double a;
    int i,j,ii;
//     printf("in matr for ortogonalization:\n");
//     for(i=0;i<n;i++)
//         PrintMatr(matr[i],1,m,0);
//     getchar();
    for(i=0;i<n;i++){
    a = cblas_ddot(m,matr[i],1,matr[i],1);
    a=1/sqrt(a);
//     if(matr[i*(m+1)]<0.0)a=-a;
    cblas_dscal(m,a,matr[i],1);
//     printf("ortogonalization of %d row:\n",i);
//     for(ii=0;ii<n;ii++)
//         PrintMatr(matr[ii],1,m,0);
    for(j=i+1;j<n;j++){
    a = cblas_ddot(m,matr[i],1,matr[j],1);
    a=-a;
    cblas_daxpy(m, a,matr[i],1,matr[j],1);
//     printf("ortogonalization of %d and %d rows:\n",i,j);
//     for(ii=0;ii<n;ii++)
//         PrintMatr(matr[ii],1,m,0);
    
    }
    }
//     printf("out matr from ortogonalization:\n");
//     for(i=0;i<n;i++)
//         PrintMatr(matr[i],1,m,0);
    // Destroy
    
    return 0;
}

double make_ort(double *a, int start_a, double *b, int start_b,int dim){
    double skal =0;
    double b_sq =0;
    int i;
//     start=start_a*dim;
    for(i=0;i<dim;i++){skal=skal+a[start_a+i]*b[start_b+i];b_sq=b_sq+b[start_b+i]*b[start_b+i];}
    skal=skal/b_sq;
    //if (skal<0.00001) return 0.0;
    for(i=0;i<dim;i++)a[start_a+i]=a[start_a+i]-skal*b[start_b+i];
    return skal;
      
}

double make_norm(double *a,int start, int dim){
    
    double a_sq =0;
    int i;
    start=start*dim;
    for(i=start;i<start+dim;i++)a_sq=a_sq+a[i]*a[i];
    a_sq=sqrt(a_sq);
    for(i=start;i<start+dim;i++)a[i]=a[i]/a_sq;
    return a_sq;
}

int matr_cpy(double * Q, int l_Q, double * O, int l_O, int l_col, int l_row){
    int i,j;
    for(i=0;i<l_col;i++)
    for(j=0;j<l_row;j++)
        O[i*l_O+j]=Q[i*l_Q+j];
    return 0;
}

int matr_cpy_vf(double * Q, int l_Q, double * O, int l_O, int l_col, int l_row, int * VF){
    int i,j;
    int jj=0;
    for(i=0;i<l_col;i++){
        jj=0;
        for(j=0;j<l_row;j++)
        if(VF[j]){O[i*l_O+jj]=Q[i*l_Q+j];jj++;}
    }
    return 0;
}

int matr_cpy_gf(double * Q, int l_Q, double * O, int l_O, int l_col, int l_row, int * GF){
    int i,j;
    int ii=0;
    for(i=0;i<l_col;i++){
        if(GF[i]){
        for(j=0;j<l_row;j++)O[ii*l_O+j]=Q[i*l_Q+j];
        ii++;}
    }
    return 0;
}


int matr_cpy_gvf(double * Q, int l_Q, double * O, int l_O, int l_col, int l_row, int * GF, int * VF){
    int i,j;
    int jj;
    int ii=0;
    for(i=0;i<l_col;i++){
        if(GF[i]){
        jj=0;
        for(j=0;j<l_row;j++)
        if(VF[j]){O[ii*l_O+jj]=Q[i*l_Q+j];jj++;
//             printf("%e ",Q[i*l_Q+j]);
        }
        ii++; /*printf("\n");*/}
    }
//     getchar();
    return 0;
}

int find_max_abs_el(double * A, int n){
    
    double max = fabs(A[0]);
    int i_max =0;
    for(int i=1;i<n;i++)
        if(fabs(A[i])>max){
            i_max=i;
            max=fabs(A[i]);
        }
    
    return i_max;
    
}

double find_max_abs_value(const double * A, int n, int m, int lda){
    
    double max = 0.0;
//     int i_max =0;
    for(int i=0;i<n;i++)
        for(int j=0;j<m;j++)
            if(fabs(A[i*lda+j])>max){
//                 i_max=i;
                max=fabs(A[i*lda+j]);
            }
    
    return max;
    
}

int find_N_min_els(int * nums, int N, double * A, int n){
    
    comp_values= new double[n];
    
    int * tmp_nums = new int[n];
    
    for(int i=0; i<n;i++){
        tmp_nums[i]=i;
        comp_values[i]=A[i];    
    }
    qsort(tmp_nums,n,sizeof ( int ),compare);
    for(int i=0; i<N;i++){
        nums[i]=tmp_nums[i];
//         printf("m[%d]=%e\n",nums[i],A[nums[i]]);
    }
//     for(int i=0; i<N;i++){
//         if(A[i]<
        
//     exit(0);
    
    delete[] comp_values;
    delete[] tmp_nums;
    
    
    return 0;
    
}

int find_max_coef(int * nums, int N, double * A, int n, int a, int ld){
    
    comp_values= new double[n];
    
    int * tmp_nums = new int[n];
    
    for(int i=0; i<n;i++){
        tmp_nums[i]=i;
        comp_values[i]=-fabs(A[i*ld+a]);
    }
    qsort(tmp_nums,n,sizeof ( int ),compare);
    for(int i=0; i<N;i++){
        nums[i]=tmp_nums[i];
//         printf("m[%d]=%e\n",nums[i],A[nums[i]]);
    }
//     for(int i=0; i<N;i++){
//         if(A[i]<
        
//     exit(0);
    
    delete[] comp_values;
    delete[] tmp_nums;
    
    
    return 0;
    
}

int find_max_coef_in_row(int * nums, int N, double * A, int n, int a) // <-- Maksim
{    
    comp_values = new double[n];
    
    int * tmp_nums = new int[n];
    
    for(int i=0; i<n; i++)
    {
        tmp_nums[i]=i;
        comp_values[i]=-fabs(A[a*n + i]);
    }
    qsort(tmp_nums,n,sizeof ( int ),compare);
    for(int i=0; i<N;i++)
    {
        nums[i]=tmp_nums[i];
    }
    
    delete[] comp_values;
    delete[] tmp_nums;    
    
    return 0;
    
}//   <--Maksim



int transform(double * O, double * In, double * T, double * buf, int dim){
    
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        dim,dim,dim,1.0,
                        In,dim,
                        T,dim,0.0,
                        buf,dim);
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                        dim,dim,dim,1.0,
                        T,dim,
                        buf,dim,0.0,
                        O,dim);
    
    return 0;
    
}

int transform_back(double * O, double * In, double * T, double * buf, int dim){
    
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                        dim,dim,dim,1.0,
                        In,dim,
                        T,dim,0.0,
                        buf,dim);
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                        dim,dim,dim,1.0,
                        T,dim,
                        buf,dim,0.0,
                        O,dim);
    
    return 0;
    
}

int transform_back_nm(double * O, double * In, double * T, double * buf, int dim_i, int dim_o){
    
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                        dim_i,dim_o,dim_i,1.0,
                        In,dim_i,
                        T,dim_o,0.0,
                        buf,dim_o);
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                        dim_o,dim_o,dim_i,1.0,
                        T,dim_o,
                        buf,dim_o,0.0,
                        O,dim_o);
    
    return 0;
    
}



int HC_SCE(double * H, double * S, double * M, double * P, double * ev, int n){
    
//        in         out
//     H  matrix    not changed
//     S  overlap   orthonormal e_vec (C), vectors are rows - C*C^(T) = E
//     M  no        e_vec  - C*S^(-0,5), vectors are rows - P*S*P^(T) = E
//     P  no        inverse e_vec P=(M^(-1))^(T) =C*S^(+0,5)
//     ev no        e_val
//     
        
    
    S05_calc(S,M,P,n);
    
//     printf("M:\n");
//     PrintMatr(M,n,n,0);
//     printf("P:\n");
//     PrintMatr(P,n,n,2);
    
    double * Int_BUF = new double[n*n];
    
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                        n,n,n,1.0,
                        H,n,
                        M,n,0.0,
                        Int_BUF,n);
    cblas_dgemm(CblasRowMajor,CblasTrans/*NoTrans*/,CblasNoTrans,//check it!!!!
                        n,n,n,1.0,
                        M,n,
                        Int_BUF,n,0.0,
                        S,n);
        
    lapack_diag(S,ev,n);

//     printf("S:\n");
//     PrintMatr(S,n,n,2);
// 
//     printf("eps:\n");
//     PrintMatr(ev,n,1,2);

    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        n,n,n,1.0,
                        S,n,
                        S,n,0.0,
                        Int_BUF,n);
    
//     printf("E:\n");
//     PrintMatr(Int_BUF,n,n,2);
    
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                        n,n,n,1.0,
                        S,n,
                        P,n,0.0,
                        Int_BUF,n);
    memcpy(P,Int_BUF,n*n*sizeof(double));
    
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                        n,n,n,1.0,
                        S,n,
                        M,n,0.0,
                        Int_BUF,n);
    memcpy(M,Int_BUF,n*n*sizeof(double));

//     printf("S2:\n");
//     PrintMatr(M,n,n,2);
//     printf("S3:\n");
//     PrintMatr(P,n,n,2);
    
    
//     printf("eps:\n");
//     PrintMatr(ev,n,1,2);

//     cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
//                         n,n,n,1.0,
//                         M,n,
//                         P,n,0.0,
//                         Int_BUF,n);
    
    
//     inv_matr_constr(Int_BUF,n);
//     printf("D:\n");
//     PrintMatr(Int_BUF,n,n,2);
    
   delete[] Int_BUF;
    
    return 1;
    
}

int HC_SCE_test(double * H, double * S, double * C, double * B, double * ev, int n){
    
//        in         out
//     H  matrix    not changed
//     S  overlap   orthonormal e_vec
//     C  no        e_vec
//     B  no        no
//     ev no        e_val
//     
    
    memcpy(B,S, n*n*sizeof(double));

    lapack_diag(S,ev,n);
    printf("C\n");
    PrintMatr(S,n,n,1);
    
    printf("n\n");
    PrintMatr(ev,n,1,1);
    
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        n,n,n,1.0,
                        H,n,
                        S,n,0.0,
                        C,n);
    cblas_dgemm(CblasRowMajor,CblasNoTrans/*NoTrans*/,CblasNoTrans,//check it!!!!
                        n,n,n,1.0,
                        S,n,
                        C,n,0.0,
                        B,n);
    
    printf("H\n");
    PrintMatr(B,n,n,1);
    
    
    for(int i=0;i<n;i++)
        ev[i]=B[i*n+i]/ev[i];
    
    printf("n\n");
    PrintMatr(ev,n,1,1);
    
    
    
    
    return 1;
    
}



// int HC_SCE_test(double * H, double * S, double * C, double * B, double * ev, int n){
//     
// //        in         out
// //     H  matrix    ort_matrix
// //     S  overlap   trash
// //     C  no        trash
// //     B  no        trash
// //     ev no        trash
// //     
//     
// //     printf("C\n");
// //     PrintMatr(C,n,n,1);
//     
//     S05_calc(S,B,C,n);
//     
//     cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
//                         n,n,n,1.0,
//                         H,n,
//                         B,n,0.0,
//                         C,n);
//     cblas_dgemm(CblasRowMajor,CblasTrans/*NoTrans*/,CblasNoTrans,//check it!!!!
//                         n,n,n,1.0,
//                         B,n,
//                         C,n,0.0,
//                         H,n);
//     
// //     printf("H_tr:\n");
// //     PrintMatr(H,n,n,1);
//     
//     
//     return 1;
//     
// }


int HPC_SPCE(double * H, double * S, double * P, double * C, double * B, double * ev, int n){
    
//        in         out
//     H  matrix    not changed
//     S  overlap   orthonormal e_vec
//     P  projector not changed
//     C  no        e_vec
//     B  no        no                 (double size)
//     ev no        e_val
//     
    double * B2 = B+n*n; 
    set_zero_matr(B2,n*n);
    for(int i=0;i<n;i++)
        B2[i*(n+1)]=1.0;//B2=1
        
    transform(C,B2,P,B,n);//P*P^T
    int non_zero = 0;
    
    
    double * S_eval;
    S_eval = new double[n];
    lapack_diag(C, S_eval, n);
//     printf("ev:\n");
//     PrintMatr(S_eval,n,1,1);

    set_zero_matr(B,n*n);
    for(int i = 0;i<n;i++){
        if(S_eval[i]>1e-8){
            for(int j = 0;j<n;j++)
                B[non_zero*n+j]=C[i*n+j]*sqrt(S_eval[i]);
            non_zero++;
        }
    }
    if(non_zero){
        double * H_proj = new double[non_zero*non_zero];
        double * S_proj = new double[non_zero*non_zero];
        
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            n,non_zero,n,1.0,
                            H,n,
                            B,n,0.0,
                            C,non_zero);
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                            non_zero,non_zero,n,1.0,
                            B,n,
                            C,non_zero,0.0,
                            H_proj,non_zero);
        
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            n,non_zero,n,1.0,
                            S,n,
                            B,n,0.0,
                            C,non_zero);
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                            non_zero,non_zero,n,1.0,
                            B,n,
                            C,non_zero,0.0,
                            S_proj,non_zero);
        
        
//         non_zero=0;
//         for(int i = 0;i<n;i++){
//             if(S_eval[i]>1e-8){
//                 for(int j = 0;j<n;j++)
//                     B[non_zero*n+j]=B[non_zero*n+j]/S_eval[i];
//                 non_zero++;
//             }
//         }
        
        HC_SCE(H_proj,S_proj,B2,C,ev,non_zero);
        
//         printf("e_vec:\n");
//         PrintMatr(B2,non_zero,non_zero,1);
        
        
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                            non_zero,n,non_zero,1.0,
                            B2,non_zero,
                            B,n,0.0,
                            C,n);
        
//         printf("e_vec:\n");
//         PrintMatr(C,n,non_zero,1);
        
        
        delete[] H_proj;
        delete[] S_proj;
    }
    delete[] S_eval;
    
    
    return non_zero;
    
    
    return 1;
    
}


int HC_SCE_p(double * H, double * S, double * M, double * P, double * ev, int n, double eps){
    
//        in         out
//     H  matrix    not changed
//     S  overlap   orthonormal e_vec (C), vectors are rows - C*C^(T) = E
//     M  no        e_vec  - C*S^(-0,5), vectors are rows - P*S*P^(T) = E
//     P  no        inverse e_vec P=(M^(-1))^(T) =C*S^(+0,5)
//     ev no        e_val
//     
    
    
    S05_calc_p(S,M,P,n,eps);
    
    
    double * Int_BUF = new double[n*n];

    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                        n,n,n,1.0,
                        H,n,
                        M,n,0.0,
                        Int_BUF,n);
    cblas_dgemm(CblasRowMajor,CblasTrans/*NoTrans*/,CblasNoTrans,//check it!!!!
                        n,n,n,1.0,
                        M,n,
                        Int_BUF,n,0.0,
                        S,n);
        
    lapack_diag(S,ev,n);

//     printf("S:\n");
//     PrintMatr(S,n,n,2);
// 
//     printf("eps:\n");
//     PrintMatr(ev,n,1,2);

    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        n,n,n,1.0,
                        S,n,
                        S,n,0.0,
                        Int_BUF,n);
    
//     printf("E:\n");
//     PrintMatr(Int_BUF,n,n,2);
    
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                        n,n,n,1.0,
                        S,n,
                        P,n,0.0,
                        Int_BUF,n);
    memcpy(P,Int_BUF,n*n*sizeof(double));
    
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                        n,n,n,1.0,
                        S,n,
                        M,n,0.0,
                        Int_BUF,n);
    memcpy(M,Int_BUF,n*n*sizeof(double));

//     printf("S2:\n");
//     PrintMatr(M,n,n,2);
//     printf("S3:\n");
//     PrintMatr(P,n,n,2);
    
    
//     printf("eps:\n");
//     PrintMatr(ev,n,1,2);

//     cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
//                         n,n,n,1.0,
//                         M,n,
//                         P,n,0.0,
//                         Int_BUF,n);
    
    
//     inv_matr_constr(Int_BUF,n);
//     printf("D:\n");
//     PrintMatr(Int_BUF,n,n,2);
    
    delete[] Int_BUF;
    
    return 1;
    
}

int HC_SCE_low(double * H, double * S, double * M, double * P, double * ev, int n, double eps){
    
//        in         out
//     H  matrix    not changed
//     S  overlap   orthonormal e_vec (C), vectors are rows - C*C^(T) = E
//     M  no        e_vec  - C*S^(-0,5), vectors are rows - P*S*P^(T) = E
//     P  no        inverse e_vec P=(M^(-1))^(T) =C*S^(+0,5)
//     ev no        e_val
//     
    
//     printf("S:\n");
//     PrintMatr(S,n,n,0);
//     printf("H:\n");
//     PrintMatr(H,n,n,2);
    
    
    
    S05_calc_p(S,M,P,n,eps);
    
    
    double * Int_BUF = new double[n*n];

    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                        n,n,n,1.0,
                        H,n,
                        M,n,0.0,
                        Int_BUF,n);
    cblas_dgemm(CblasRowMajor,CblasTrans/*NoTrans*/,CblasNoTrans,//check it!!!!
                        n,n,n,1.0,
                        M,n,
                        Int_BUF,n,0.0,
                        S,n);
        
    lapack_diag(S,ev,n);

//     fprintf(out_stream,"S:\n");
//     fPrintMatr(out_stream,S,n,n,2);

//     fprintf(out_stream,"eps:\n");
//     fPrintMatr(out_stream,ev,n,1,0);
    
    for(int i=0;i<n;i++)
        if(fabs(ev[i])<1e-5)
        for(int j=0;j<n;j++)S[i*n+j]=0.0;
    

    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        n,n,n,1.0,
                        S,n,
                        S,n,0.0,
                        Int_BUF,n);
    
//     printf("E:\n");
//     PrintMatr(Int_BUF,n,n,2);
    
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                        n,n,n,1.0,
                        S,n,
                        P,n,0.0,
                        Int_BUF,n);
    memcpy(P,Int_BUF,n*n*sizeof(double));
    
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                        n,n,n,1.0,
                        S,n,
                        M,n,0.0,
                        Int_BUF,n);
    memcpy(M,Int_BUF,n*n*sizeof(double));

//     printf("S2:\n");
//     PrintMatr(M,n,n,0);
//     printf("S3:\n");
//     PrintMatr(P,n,n,2);
    
    
//     printf("eps:\n");
//     PrintMatr(ev,n,1,2);

//     cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
//                         n,n,n,1.0,
//                         M,n,
//                         P,n,0.0,
//                         Int_BUF,n);
    
    
//     inv_matr_constr(Int_BUF,n);
//     printf("D:\n");
//     PrintMatr(Int_BUF,n,n,2);
    
    delete[] Int_BUF;
    
    return 1;
    
}



int HSC_CE(double * H, double * S, double * C, double * B, double * ev, int n){
    
//        in         out
//     H  matrix    not changed
//     S  overlap   orthonormal e_vec
//     C  no        e_vec
//     B  no        no
//     ev no        e_val
//     
    double * Sm05 = new double[n*n]; 
//    printf("B\n");
//    PrintMatr(B,n,n,1);
    
    S05_calc(S,Sm05,B,n);
    
//    PrintMatr(B,n,n,1);
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                        n,n,n,1.0,
                        H,n,
                        B,n,0.0,
                        C,n);
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                        n,n,n,1.0,
                        B,n,
                        C,n,0.0,
                        S,n);
    
//     printf("H_tr:\n");
//     PrintMatr(S,n,n,0);
    
    
    lapack_diag(S,ev,n);
    
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                        n,n,n,1.0,
                        S,n,
                        Sm05,n,0.0,
                        C,n);
    
    
    return 1;
    
}


int HPSPC_PCE(double * H, double * S, double * P, double * C, double * B, double * ev, int n){
    
//        in         out           dim
//     H  matrix    not changed    n*n
//     S  overlap   not changed    n*n
//     C  no        e_vec          n*n
//     B  no        no             2*n*n !!!!
//     ev no        e_val          n
//     
    double * B2 = B+n*n; 
    set_zero_matr(B2,n*n);
    for(int i=0;i<n;i++)
        B2[i*(n+1)]=1.0;
    
    transform(C,B2,P,B,n);
    int non_zero = 0;
    
    double * S_eval;
    S_eval = new double[n];
    lapack_diag(C, S_eval, n);

    set_zero_matr(B,n*n);
    for(int i = 0;i<n;i++){
        if(S_eval[i]>1e-8){
//             printf("%e\n",S_eval[i]);
//             getchar();
            for(int j = 0;j<n;j++)
                B[non_zero*n+j]=C[i*n+j]*sqrt(S_eval[i]);
            non_zero++;
        }
    }

    
    if(non_zero){
        double * H_proj = new double[non_zero*non_zero];
        double * S_proj = new double[non_zero*non_zero];
        
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            n,non_zero,n,1.0,
                            H,n,
                            B,n,0.0,
                            C,non_zero);
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                            non_zero,non_zero,n,1.0,
                            B,n,
                            C,non_zero,0.0,
                            H_proj,non_zero);
        
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            n,non_zero,n,1.0,
                            S,n,
                            B,n,0.0,
                            C,non_zero);
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                            non_zero,non_zero,n,1.0,
                            B,n,
                            C,non_zero,0.0,
                            S_proj,non_zero);
        
        non_zero=0;
        //must be tested!!!!
//         for(int i = 0;i<n;i++){
//             if(S_eval[i]>1e-8){
//                 for(int j = 0;j<n;j++)
//                     B[non_zero*n+j]=B[non_zero*n+j]/S_eval[i];
//                 non_zero++;
//             }
//         }
        
        HSC_CE(H_proj,S_proj,B2,C,ev,non_zero);
        
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                            non_zero,n,non_zero,1.0,
                            B2,non_zero,
                            B,n,0.0,
                            C,n);
        
        
        delete[] H_proj;
        delete[] S_proj;
    }
    
    delete[] S_eval;
    
    
    return non_zero;
    
}


int HSC_CE_p(double * H, double * S, double * C, double * B, double * ev, int n, double eps){
    
//        in         out
//     H  matrix    not changed
//     S  overlap   orthonormal e_vec
//     C  no        e_vec
//     B  no        no
//     ev no        e_val
//     
    double * Sm05 = new double[n*n]; 
//    printf("B\n");
//    PrintMatr(B,n,n,1);
    
    S05_calc_p(S,Sm05,B,n,eps);
    
//    PrintMatr(B,n,n,1);
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                        n,n,n,1.0,
                        H,n,
                        B,n,0.0,
                        C,n);
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                        n,n,n,1.0,
                        B,n,
                        C,n,0.0,
                        S,n);
    
//     printf("H_tr:\n");
//     PrintMatr(S,n,n,0);
    
    
    lapack_diag(S,ev,n);
    
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                        n,n,n,1.0,
                        S,n,
                        Sm05,n,0.0,
                        C,n);
    
    
    return 1;
    
}

int symmetrization(double * A, int n){
    
    for(int i=0  ;i<n;i++)
    for(int j=i+1;j<n;j++){
        A[i*n+j]=0.5*(A[i*n+j]+A[j*n+i]);
        A[j*n+i]=A[i*n+j];
    }
    
    return 0;
}

int symmetrization_with_scaling(double * A, int n, double c){
    
    for(int i=0  ;i<n;i++)
    for(int j=i;j<n;j++){
        A[i*n+j]=0.5*(A[i*n+j]+A[j*n+i])*c;
        A[j*n+i]=A[i*n+j];
    }
    
    return 0;
}

int anti_symmetrization(double * A, int n){
    
    for(int i=0  ;i<n;i++)
    for(int j=i;j<n;j++){
        A[i*n+j]=0.5*(A[i*n+j]-A[j*n+i]);
        A[j*n+i]=-A[i*n+j];
    }
    
    return 0;
}


int anti_symmetrization_with_scaling(double * A, int n, double c){
    
    for(int i=0  ;i<n;i++)
    for(int j=i;j<n;j++){
        A[i*n+j]=0.5*(A[i*n+j]-A[j*n+i])*c;
        A[j*n+i]=-A[i*n+j];
    }
    
    return 0;
}

double sym_check_for_tensor(double * A, int N, int n, int ld, double f){
    
    double R=0;
    double M=0;
    
    for(int i=0; i<N; i++)
    for(int j=0; j<N; j++)
    for(int k=0; k<n; k++){
        M = fabs(A[(i*N+j)*ld+k]+f*A[(j*N+i)*ld+k]);
        if(R<M)R=M;
    }
        
    return R;
    
}

int tensor_symmetrization(double * A, double * s, int N, int n, int ld){
    
    
    for(int i=0; i<N; i++)
    for(int j=i; j<N; j++)
    for(int k=0; k<n; k++){
        A[(i*N+j)*ld+k]=0.5*(A[(i*N+j)*ld+k]+s[k]*A[(j*N+i)*ld+k]);
        A[(j*N+i)*ld+k]=A[(i*N+j)*ld+k]*s[k];
    }
    
    return 0;
}


int Cholesky_step(double * L, double * M, int n, double eps){
    
    
    double M_max=0;
    int i_max=0;
    
    
    for(int i=0;i<n;i++)
        if(M[i*n+i]>M_max){
            i_max=i;
            M_max=M[i*n+i];
        }
        
//     printf("i_max=%d %e\n", i_max, M_max);
//     getchar();
    if(M_max<eps)
        return 1;
    
    double sq = sqrt(M_max);
    for(int i=0;i<n;i++)
        L[i]=M[i_max*n+i]/sq;
    
//     for(int i=0;i<n;i++)
//         printf("%e\n",L[i]);
    
    for(int i=0;i<n;i++)
    for(int j=0;j<n;j++)
        M[i*n+j]-=L[i]*L[j];
    
    
    
    return 0;
    
}
    
int Cholesky(double * L, double * M, int n, double eps){
    
    set_zero_matr(L,n*n);
    int i;
    for(i=0;i<n;i++){
//         printf("i=%d  ",i);
        if(Cholesky_step(L+i*n, M, n, eps)) break;
    }
    
    printf("i=%d/%d\n",i,n);
//     printf("Cholesky\n");
//     PrintMatr(L,n,n,1);
    
    return 0;
}
    
    
    

