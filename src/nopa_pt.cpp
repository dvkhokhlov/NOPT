# include <stdio.h>
# include <stdlib.h>
// 
# include "blas_link.h"
# include "matr.h"
# include "inp_out.h"
# include "nopa_pt.h"
# include "grabbers.h"

#define PT_V
//#define DOCI_V


int nopa_v_data::alloc(int n){
    
    n_block=n;
    
    V  = new double [CI_dim*CI_dim];
    set_zero_matr(V,CI_dim*CI_dim);
    
    T2  = new double [CI_dim*CI_dim];
    
    E_PT=new std::vector<double>[2*n];
    V_A  =new double*[n];
    S_A  =new double*[n];
    V_B  =new double*[n];
    S_B  =new double*[n];
    V_A_0=new double*[n];
    V_B_0=new double*[n];
//     double * T2;//bufer
//     v_name = new char[255];
    
    n_a = new int[n];
    n_b = new int[n];
    
    S_ort = new double[CI_dim*CI_dim]; 
    set_zero_matr(S_ort,CI_dim*CI_dim);
    for(int i=0;i<CI_dim;i++)
        S_ort[i*(CI_dim+1)]=1.0;
    
    
    
    return 0;
}

int nopa_v_data::alloc2(int n, int ext_CI_dim){
    
    CI_dim=ext_CI_dim;
    alloc(n);
    
    return 0;
}


int nopa_v_data::V_read(int i, char * name){
    
    int wrong_data=0;
    char * v_name = new char[255];
    
#ifdef PT_V        
    E_PT[2*i  ] = read_E_XMC(name,"_A.mp2");
    E_PT[2*i+1] = read_E_XMC(name,"_B.mp2");
    
    n_a[i]=E_PT[2*i  ].size();
    n_b[i]=E_PT[2*i+1].size();
#endif
#ifdef DOCI_V
    sprintf(v_name, "%sVa\0", name);
    n_a[i]=read_block_size(v_name);
    E_PT[2*i  ].resize(n_a[i]*n_a[i]);
    if(read_HS_data(v_name,v_name,E_PT[2*i  ].data(),T2,T2,T2,T2,n_a[i],n_a[i],n_a[i]))wrong_data=1;
    fprintf(out_stream,"X_a\n");
    fPrintMatr(out_stream,E_PT[2*i  ].data(),n_a[i],n_a[i],0);
    
    
    sprintf(v_name, "%sVb\0", name);
    n_b[i]=read_block_size(v_name);
    E_PT[2*i+1].resize(n_b[i]*n_b[i]);
    if(read_HS_data(v_name,v_name,E_PT[2*i+1].data(),T2,T2,T2,T2,n_b[i],n_b[i],n_b[i]))wrong_data=1;
    fprintf(out_stream,"X_b\n");
    fPrintMatr(out_stream,E_PT[2*i+1].data(),n_b[i],n_b[i],0);
        
#endif
    if(block_dim[i]!=n_a[i]*n_b[i]){
        fprintf(out_stream," ERROR: inconsistent read in block %d:\n",i);
        fprintf(out_stream,"%d != %d * %d\n", block_dim[i], n_a[i],n_b[i]);
        exit(1);
    }
    V_A  [i]=new double[n_a[i]*n_a[i]];
    S_A  [i]=new double[n_a[i]*n_a[i]];
    V_B  [i]=new double[n_b[i]*n_b[i]];
    S_B  [i]=new double[n_b[i]*n_b[i]];
    V_A_0[i]=new double[n_a[i]*n_a[i]];
    V_B_0[i]=new double[n_b[i]*n_b[i]];
    
    if(read_HS_data(name,"_A",V_A_0[i],S_A[i],T2,T2,T2,n_a[i],n_a[i],n_a[i]))wrong_data=1;
    fprintf(out_stream,"D_a\n");
    fPrintMatr(out_stream,V_A_0[i],n_a[i],n_a[i],0);
     
    
    if(read_HS_data(name,"_B",V_B_0[i],S_B[i],T2,T2,T2,n_b[i],n_b[i],n_b[i]))wrong_data=1;
    fprintf(out_stream,"D_b\n");
    fPrintMatr(out_stream,V_B_0[i],n_a[i],n_a[i],0);
    
    return wrong_data;
}

int nopa_v_data::get_V_from_global_memory(int i, int a){
    
    E_PT[2*i+a] = E_grabbed;
    
    if(a==0){
        n_a  [i]=E_grabbed.size();
        V_A  [i]=new double[n_a[i]*n_a[i]];
        S_A  [i]=new double[n_a[i]*n_a[i]];
        V_A_0[i]=new double[n_a[i]*n_a[i]];
        memcpy(V_A_0[i], H_grabbed, n_a[i]*n_a[i]*sizeof(double));
        memcpy(S_A  [i], S_grabbed, n_a[i]*n_a[i]*sizeof(double));
        // fprintf(out_stream,"D_a\n");
        // fPrintMatr(out_stream,V_A_0[i],n_a[i],n_a[i],0);
    }
    if(a==1){
        n_b  [i]=E_grabbed.size();
        V_B [i]=new double[n_b[i]*n_b[i]];
        S_B [i]=new double[n_b[i]*n_b[i]];
        V_B_0[i]=new double[n_b[i]*n_b[i]];
        memcpy(V_B_0[i], H_grabbed, n_b[i]*n_b[i]*sizeof(double));
        memcpy(S_B  [i], S_grabbed, n_b[i]*n_b[i]*sizeof(double));
        // fprintf(out_stream,"D_b\n");
        // fPrintMatr(out_stream,V_B_0[i],n_b[i],n_b[i],0);
    }
    
    
    return 0;
}

int nopa_v_data::copy_V(int i, int j){
        
    E_PT[2*i  ] = E_PT[2*j  ];
    E_PT[2*i+1] = E_PT[2*j+1];
    
    n_a  [i]=n_a[j];
    V_A  [i]=new double[n_a[i]*n_a[i]];
    S_A  [i]=new double[n_a[i]*n_a[i]];
    V_A_0[i]=new double[n_a[i]*n_a[i]];
    // memcpy(V_A  [i], V_A  [j], n_a[i]*n_a[i]*sizeof(double));
    // memcpy(S_A  [i], S_A  [j], n_a[i]*n_a[i]*sizeof(double));
    memcpy(V_A_0[i], V_A_0[j], n_a[i]*n_a[i]*sizeof(double));
    memcpy(S_A  [i], S_A  [j], n_a[i]*n_a[i]*sizeof(double));
    fprintf(out_stream,"D_a\n");
    fPrintMatr(out_stream,V_A_0[i],n_a[i],n_a[i],0);
    
    
    n_b  [i]=n_b[j];
    V_B  [i]=new double[n_b[i]*n_b[i]];
    S_B  [i]=new double[n_b[i]*n_b[i]];
    V_B_0[i]=new double[n_b[i]*n_b[i]];
    // memcpy(V_B  [i], V_B  [j], n_b[i]*n_b[i]*sizeof(double));
    // memcpy(S_B  [i], S_B  [j], n_b[i]*n_b[i]*sizeof(double));
    memcpy(V_B_0[i], V_B_0[j], n_b[i]*n_b[i]*sizeof(double));
    memcpy(S_B  [i], S_B  [j], n_b[i]*n_b[i]*sizeof(double));
    fprintf(out_stream,"D_b\n");
    fPrintMatr(out_stream,V_B_0[i],n_a[i],n_a[i],0);
    
    return 0;
}

int nopa_v_data::V_read_all(char ** name){
    
    int wrong_data=0;
    for(int i=0; i<n_block; i++){
        if(V_read(i, name[i+1]))wrong_data=1;
    }
    
    if(wrong_data){
        fprintf(out_stream,"There were some errors in non-diagonal blocks, program terminated.\n");
        exit(1);
    }
    
    return 0;
}


int nopa_v_data::V_block_recalc(int i){
    
    for(int ii=0;ii<n_a[i];ii++)
    for(int jj=0;jj<n_a[i];jj++)
#ifdef PT_V        
        V_A[i][ii*n_a[i]+jj]=E_PT[2*i  ][jj]*S_A[i][ii*n_a[i]+jj]-V_A_0[i][ii*n_a[i]+jj];
#endif
#ifdef DOCI_V
        V_A[i][ii*n_a[i]+jj]=E_PT[2*i  ][ii*n_a[i]+jj]-V_A_0[i][ii*n_a[i]+jj];
#endif
        
    for(int ii=0;ii<n_b[i];ii++)
    for(int jj=0;jj<n_b[i];jj++)
#ifdef PT_V        
        V_B[i][ii*n_b[i]+jj]=E_PT[2*i+1][jj]*S_B[i][ii*n_b[i]+jj]-V_B_0[i][ii*n_b[i]+jj];
#endif
#ifdef DOCI_V
        V_B[i][ii*n_b[i]+jj]=E_PT[2*i+1][ii*n_b[i]+jj]-V_B_0[i][ii*n_b[i]+jj];
#endif

    fprintf(out_stream,"V_A1\n");fPrintMatr(out_stream,V_A[i],n_a[i],n_a[i],0);
    fprintf(out_stream,"S_A1\n");fPrintMatr(out_stream,S_A[i],n_a[i],n_a[i],0);
    fprintf(out_stream,"V_B1\n");fPrintMatr(out_stream,V_B[i],n_b[i],n_b[i],0);
    fprintf(out_stream,"S_B1\n");fPrintMatr(out_stream,S_B[i],n_b[i],n_b[i],0);
    
    return 0;
}

int nopa_v_data::V_block_recalc_all(){
    
    for(int i=0; i<n_block; i++){
        V_block_recalc(i);
    }
    
    return 0;
}

int nopa_v_data::V_orthogonal_calc(){
    
    double * V_block;
    int ij,kl;
    for(int i_b=0; i_b<n_block; i_b++)
    {
        V_block= V + dim[i_b]*CI_dim+dim[i_b];
//         fprintf(out_stream,"%d, %d, %d\n",block_dim[i_b],n_a[i_b],n_b[i_b]);
        
//         if(block_dim[i_b]==(n_a[i_b]*n_b[i_b])){
//             fprintf(out_stream,"block %d is full\n",i_b);

            for(int i=0;i<n_a[i_b];i++)
            for(int j=0;j<n_b[i_b];j++)
            for(int k=0;k<n_a[i_b];k++)
            for(int l=0;l<n_b[i_b];l++){
                ij=i*n_b[i_b]+j;
                kl=k*n_b[i_b]+l;
                
                V_block[ij*CI_dim+kl]=V_A[i_b][i*n_a[i_b]+k]*S_B[i_b][j*n_b[i_b]+l]+
                                      S_A[i_b][i*n_a[i_b]+k]*V_B[i_b][j*n_b[i_b]+l];
                
                
                }
        
    }
    
    return 0;
}

int nopa_v_data::mulliken_update(double * M){
    
    set_zero_matr(T2,CI_dim*CI_dim);
    
    for(int i=0;i<CI_dim;i++)
    for(int j=0;j<CI_dim;j++)
    for(int k=0;k<CI_dim;k++){
        T2[i*CI_dim+j]+=M[i*CI_dim+k]*V[k*CI_dim+j];
    }
    
    memcpy(V,T2,CI_dim*CI_dim*sizeof(double));
    
    return 0;
}

int nopa_v_data::leowdin_update(double * M, int inv){
    
    memcpy(T2,M,CI_dim*CI_dim*sizeof(double));
    
    double * Mp05 = new double[CI_dim*CI_dim];
    double * Mm05 = new double[CI_dim*CI_dim];
    
    S05_calc(T2,Mm05,Mp05,CI_dim);
    
    double * Mz05;
    if(inv)Mz05=Mm05;
    else   Mz05=Mp05;
    
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                CI_dim,CI_dim,CI_dim,1.0,
                Mz05,CI_dim,
                V,CI_dim,0.0,
                T2,CI_dim);
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                CI_dim,CI_dim,CI_dim,1.0,
                T2,CI_dim,
                Mz05,CI_dim,0.0,
                V,CI_dim);
    delete[] Mp05;
    delete[] Mm05;
//     memcpy(V,T2,CI_dim*CI_dim*sizeof(double));
    
    return 0;
}


int nopa_v_data::read_1el(int a, int b, char ** name){
    
    char * v_nd_name = new char[1024];
    sprintf(v_nd_name, "%s%s.smo\0", name[a+1],name[b+1]);
    double * v_nd = new double[n_a[a]*n_b[a]*n_a[b]*n_b[b]];
    
    if(ReadMatr10(v_nd, n_a[a]*n_a[b],n_b[a]*n_b[b], v_nd_name, "<fi(A)|fi(B)>:"))
        return 0;
    
    fprintf(out_stream,"%s\n",v_nd_name);
    fPrintMatr(out_stream,v_nd,n_a[a]*n_b[a],n_a[b]*n_b[b],0);
        
//     exit(0);    
    
    int shift_g=0;
    int shift_v=0;
    for(int i=0;i<a;i++){
        shift_g+=n_a[i]*n_b[i]*CI_dim;
        shift_v+=n_a[i]*n_b[i];
    }
    for(int i=0;i<b;i++){
        shift_g+=n_a[i]*n_b[i];
        shift_v+=n_a[i]*n_b[i]*CI_dim;
    }
    
    
    
    for(int i=0;i<n_a[b];i++)//A(-)
    for(int j=0;j<n_b[b];j++)//B(+)
    for(int k=0;k<n_a[a];k++)//A(0)
    for(int l=0;l<n_b[a];l++)//B(0)
    {
        int i_state =i*n_b[b]*CI_dim       +j*CI_dim       +k*n_b[a]+l+shift_v;
//         (i(A)*n+j(A+))*N+i(B-)*n+j(B)
        int i_orb   =i*n_a[b]*n_b[a]*n_b[b]+k*n_b[a]*n_b[b]+l*n_b[a]+j;
        
        S_ort[i_state]-=v_nd[i_orb];
    }
    
    for(int i=0;i<n_a[a];i++)//A(0)
    for(int j=0;j<n_b[a];j++)//B(0)
    for(int k=0;k<n_a[b];k++)//A(-)
    for(int l=0;l<n_b[b];l++)//B(+)
    {
        int i_state =i*n_b[a]*CI_dim       +j*CI_dim       +k*n_b[b]+l+shift_g;
//         (i(A)*n+j(A+))*N+i(B-)*n+j(B)
        int i_orb   =k*n_a[b]*n_b[a]*n_b[b]+i*n_b[a]*n_b[b]+j*n_b[a]+l;
        S_ort[i_state]-=v_nd[i_orb];
    }
    
    for(int i=0;i<n_a[b];i++)//A(-)
    for(int j=0;j<n_b[b];j++)//B(+)
    for(int k=0;k<n_a[a];k++)//A(0)
    for(int l=0;l<n_b[a];l++)//B(0)
    {
        int i_state =i*n_b[b]*CI_dim       +j*CI_dim       +k*n_b[a]+l+shift_v;
//         int i_orb   =i*n_a[b]*n_b[a]*n_b[b]+k*n_b[a]*n_b[b]+l*n_b[a]+j;
        
        for(int n=0;n<n_a[b];n++){
            int i_orb   =i*n_a[b]*n_b[a]*n_b[b]+n*n_b[a]*n_b[b]+l*n_b[a]+j;
            V[i_state]-= V_A[a][k*n_a[b]+n]*v_nd[i_orb];
        }
        for(int n=0;n<n_b[a];n++){
            int i_orb   =i*n_a[b]*n_b[a]*n_b[b]+k*n_b[a]*n_b[b]+l*n_b[a]+n;
            V[i_state]-= V_B[b][j*n_b[a]+n]*v_nd[i_orb];
        }
        
    }
//     fprintf(out_stream,"V1:\n");
//     fPrintMatr10(out_stream,V,CI_dim,CI_dim,1);
    ReadMatr10(v_nd, n_a[a]*n_b[a],n_a[b]*n_b[b], v_nd_name, "<fi(A)|V|fi(B)>:");
//     fprintf(out_stream,"%s\n",v_nd_name);
//     fPrintMatr(out_stream,v_nd,n_a[a]*n_b[a],n_a[b]*n_b[b],0);
    
    
    for(int i=0;i<n_a[b];i++)//A(-)
    for(int j=0;j<n_b[b];j++)//B(+)
    for(int k=0;k<n_a[a];k++)//A(0)
    for(int l=0;l<n_b[a];l++)//B(0)
    {
        int i_state  =i*n_b[b]*CI_dim       +j*CI_dim       +k*n_b[a]+l+shift_v;
        int i_orb   =i*n_a[b]*n_b[a]*n_b[b]+k*n_b[a]*n_b[b]+l*n_b[a]+j;
        
        V[i_state]-=v_nd[i_orb];
        
    }
    
    
    
    for(int i=0;i<n_a[b];i++)//A(+)
    for(int j=0;j<n_b[b];j++)//B(-)
    for(int k=0;k<n_a[a];k++)//A(0)
    for(int l=0;l<n_b[a];l++)//B(0)
    {
        int i_state  =i*n_b[b]*CI_dim       +j*CI_dim       +k*n_b[a]+l+shift_v;
        int i_state2 =k*n_b[a]*CI_dim       +l*CI_dim       +i*n_b[b]+j+shift_g;
        
        V[i_state2]=V[i_state];
        
    }
//     fprintf(out_stream,"V2:\n");
//     fPrintMatr10(out_stream,V,CI_dim,CI_dim,0);
    
    delete[] v_nd_name;
    delete[] v_nd;
    
    return 0;
}

int nopa_v_data::read_1el_all(char ** name){
    
    for(int i=0  ;i<n_block;i++)
    for(int j=i+1;j<n_block;j++)
        read_1el(i,j,name);
    
    return 1;
}

int nopa_v_data::H_update(){
    
    for(int i=0;i<CI_dim*CI_dim;i++){
        H[i]+=V[i];
    }
    
    return 0;
}

int nopa_v_data::H_update_v_ort(double * ext_H, double * ext_S, int argc, char ** argv, int ext_CI_dim, int * ext_dim, int * ext_block_dim){
    
    fprintf(out_stream,"\n");
    fprintf(out_stream,"Start PT hamiltonian update.\n");
    fprintf(out_stream,"\n");
    fprintf(out_stream,"Using update without account of overlap\n\n");
    
    H         = ext_H        ;
    S         = ext_S        ;
    dim       = ext_dim      ;
    block_dim = ext_block_dim;
    CI_dim    = ext_CI_dim   ;
    
    
    alloc(argc);
    
    V_read_all(argv);
    
    V_block_recalc_all();
    
    V_orthogonal_calc();
    
    fprintf(out_stream,"V:\n");
    fPrintMatr(out_stream,V,CI_dim,CI_dim,0);
    
    H_update();
    
    
    return 1;
}


int nopa_v_data::H_update_simple(double * ext_H, double * ext_S, int argc, char ** argv, int ext_CI_dim, int * ext_dim, int * ext_block_dim){
    
    fprintf(out_stream,"\n");
    fprintf(out_stream,"Start PT hamiltonian update.\n");
    fprintf(out_stream,"\n");
    fprintf(out_stream,"Using mulliken projector method without perturbation eigenvectors\n\n");
    
    H         = ext_H        ;
    S         = ext_S        ;
    dim       = ext_dim      ;
    block_dim = ext_block_dim;
    CI_dim    = ext_CI_dim   ;
    
    
    alloc(argc);
    
    V_read_all(argv);
    
    V_block_recalc_all();
    
    V_orthogonal_calc();
    
    mulliken_update(S);
    symmetrization(V,CI_dim);
    fprintf(out_stream,"V:\n");
    fPrintMatr(out_stream,V,CI_dim,CI_dim,0);
    
    H_update();
    
    
    return 1;
}

int nopa_v_data::H_update_simple_with_non_diag_V(double * ext_H, double * ext_S, int argc, char ** argv, int ext_CI_dim, int * ext_dim, int * ext_block_dim){
    
    fprintf(out_stream,"\n");
    fprintf(out_stream,"Start PT hamiltonian update.\n");
    fprintf(out_stream,"\n");
    fprintf(out_stream,"Using mulliken projector method with 1-electron extension\n\n");
    
    H         = ext_H        ;
    S         = ext_S        ;
    dim       = ext_dim      ;
    block_dim = ext_block_dim;
    CI_dim    = ext_CI_dim   ;
    
    
    alloc(argc);
    
    V_read_all(argv);
    
    V_block_recalc_all();
    
    V_orthogonal_calc();
    
    read_1el_all(argv);
    
    fprintf(out_stream,"S(orthogonal functions):\n");
    fPrintMatr(out_stream,S_ort,CI_dim,CI_dim,0);
    
    inv_matr_constr(S_ort,CI_dim);
    mulliken_update(S_ort);
    mulliken_update(S);
    symmetrization(V,CI_dim);
    
    fprintf(out_stream,"V:\n");
    fPrintMatr10(out_stream,V,CI_dim,CI_dim,0);
    
    H_update();
    
    return 1;
}

int nopa_v_data::H_update_leow(double * ext_H, double * ext_S, int argc, char ** argv, int ext_CI_dim, int * ext_dim, int * ext_block_dim){
    
    fprintf(out_stream,"\n");
    fprintf(out_stream,"Start PT hamiltonian update.\n");
    fprintf(out_stream,"\n");
    fprintf(out_stream,"Using leowdin projector method without perturbation eigenvectors\n\n");
    
    H         = ext_H        ;
    S         = ext_S        ;
    dim       = ext_dim      ;
    block_dim = ext_block_dim;
    CI_dim    = ext_CI_dim   ;
    
    
    alloc(argc);
    
    V_read_all(argv);
    
    V_block_recalc_all();
    
    V_orthogonal_calc();
    
    leowdin_update(S,0);
//     symmetrization(V,CI_dim);
    fprintf(out_stream,"V:\n");
    fPrintMatr(out_stream,V,CI_dim,CI_dim,0);
    
    H_update();
    
    
    return 1;
}

int nopa_v_data::H_update_leow_with_non_diag_V(double * ext_H, double * ext_S, int argc, char ** argv, int ext_CI_dim, int * ext_dim, int * ext_block_dim){
    
    fprintf(out_stream,"\n");
    fprintf(out_stream,"Start PT hamiltonian update.\n");
    fprintf(out_stream,"\n");
    fprintf(out_stream,"Using leowdin projector method with 1-electron extension\n\n");
    
    H         = ext_H        ;
    S         = ext_S        ;
    dim       = ext_dim      ;
    block_dim = ext_block_dim;
    CI_dim    = ext_CI_dim   ;
    
    
    alloc(argc);
    
    V_read_all(argv);
    
    V_block_recalc_all();
    
    V_orthogonal_calc();
    
    read_1el_all(argv);
    
    fprintf(out_stream,"S(orthogonal functions):\n");
    fPrintMatr(out_stream,S_ort,CI_dim,CI_dim,0);
    
    leowdin_update(S_ort,1);
    leowdin_update(S,0);
//     symmetrization(V,CI_dim);
    fprintf(out_stream,"V:\n");
    fPrintMatr(out_stream,V,CI_dim,CI_dim,0);
    
    H_update();
    
    
    return 1;
}

int nopa_v_data::H_update_leow_G(double * ext_H, double * ext_S, int ext_CI_dim, int * ext_dim, int * ext_block_dim){
    
    fprintf(out_stream,"\n");
    fprintf(out_stream,"Start PT hamiltonian update.\n");
    fprintf(out_stream,"\n");
    fprintf(out_stream,"Using leowdin projector method without perturbation eigenvectors\n\n");
    
    H         = ext_H        ;
    S         = ext_S        ;
    dim       = ext_dim      ;
    block_dim = ext_block_dim;
    CI_dim    = ext_CI_dim   ;
    
        
    V_block_recalc_all();
    
    V_orthogonal_calc();
    
    leowdin_update(S,0);
//     symmetrization(V,CI_dim);
    fprintf(out_stream,"V:\n");
    fPrintMatr(out_stream,V,CI_dim,CI_dim,0);
    
    H_update();
    
    
    return 1;
}
