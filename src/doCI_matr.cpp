//standart
# include <stdio.h>
#include <iostream>
# include <libint2/engine.h>


//user
# include "blas_link.h"
# include "molecule.h"
# include "molecule2.h"
# include "libint_link.h"
# include "matr.h"
# include "timer.h"
# include "CI.h"
# include "etc.h"
# include "doCI_data.h"
# include "doCI_matr.h"
# include "from_hash.h"
# include "inp_par_read.h"


using std::vector;
using libint2::Shell;
using libint2::Engine;
using libint2::Operator;

// #define READ_H1
#define CALC_H1
// #define WRITE_H1
#define CI_TYPE 1 //0 - doCI, 1 - CI


double svd_eps=1e-9;

extern int num_threads;


double E_1el_calc(double * H, double * DM, int n1, int n2){
    
    double E=0;
    
    for(int i=0;i<n1;i++) for(int j=0;j<n2;j++){ E+=DM[i+n2*j]*H[i*n2+j];
//        printf("E=%f\n",E);
    }
    
    return E;
}

double E_1el_calc_2(double * H, double * DM, int n1, int n2){
    
    double E=0;
    
    for(int i=0;i<n1*n2;i++){ E+=DM[i]*H[i];
//        printf("E=%f\n",E);
    }
    
    return E;
}



double E_1el_2MO_calc(double * P, double * B, double * K, double * BUF, int n){
     
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                    n,1,n,1.0,
                    P,n,
                    K,n,0.0,
                    BUF,1);
        
     return cblas_ddot(n,B,1,BUF,1);
     
}


int gen_S_DD(double *D, double *S, int n_1, int d_1, int n_2, int d_2){
    
    for(int i=0; i<d_2; i++)
        for(int j=0; j<d_1; j++)
            D[i*d_1+j]=S[i*n_1+j];
    
    return 0;
}

int act_to_cor_ort(double * O,double * S,int d_o,int d_n, int v, int n, double * T,double * V, double * L, double * D,int trans){
    if(v==0) return 1;
    double * TMP;
    TMP = new double[v*d_n];
    
//     printf("S\n"); 
//     PrintMatr(S,v,n,1);
//     printf("T\n");
//     PrintMatr(T,d_n,d_n,1);
//     printf("V\n");
//     PrintMatr(V,v,n,1);
//     printf("D\n");
//     PrintMatr(D,d_n,n,1);
//     printf("L\n");
//     PrintMatr(L,1,d_n,1);
    
    CBLAS_TRANSPOSE Trans1=CblasTrans;
    CBLAS_TRANSPOSE Trans2=CblasNoTrans;
    if(trans){
        Trans1=CblasNoTrans;
        Trans2=CblasTrans;
            
    }    
    cblas_dgemm(CblasRowMajor,Trans1,Trans2,
                v,d_n,d_o,1.0,
                S,n,
                T,d_o,0.0,
                TMP,d_n);
    
    for(int i=0;i<v;i++)
        for(int j=0;j<d_n;j++)
            TMP[i*d_n+j]=TMP[i*d_n+j]/L[j];

        
    for(int i=0;i<v;i++)
    for(int j=0;j<n;j++)
            O[j*v+i]=V[i*n+j];
        
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasTrans,
                n,v,d_n,-1.0,
                D,n,
                TMP,d_n,1.0,
                O,v);
    
    for(int i=d_n;i<d_o;i++){
        printf("adding SVD orbital %d to active\n",i);
        for(int j=0;j<n;j++)
        O[j*v+i-d_n]=D[i*n+j];
    }
    
    delete[] TMP;
    
    return 0;
    
}
    
int act_to_act_ort(double * O,double * S, int n_b, int n_t, int ld, int trans){
    
    //n_b - dimension of base vectors
    //n_b - dimension of vectors to be transformed
    //ld  - leading dimension of S
    
    if(n_b==0) return 0;
    
    set_zero_matr(O,n_t*(n_b+n_t));
    for (int i=0;i<n_t;i++) O[i*(n_b+n_t+1)+n_b]=1.0;
    
    double * T;
    if(trans){
        T=S+ld*n_b;
        for(int i=0;i<n_t;i++)
        for(int j=0;j<n_b;j++)
            O[i*(n_b+n_t)+j]=-T[i*ld+j];
    }
    else{
        T=S+n_b;
        for(int i=0;i<n_t;i++)
        for(int j=0;j<n_b;j++)
            O[i*(n_b+n_t)+j]=-T[j*ld+i];
    }
    // printf("O1:\n");
    // PrintMatr(O,n_t,n_b+n_t,1);
    return 1;
    
}

int act_to_act_so_matr(double * O, double * S, double * U, int n_b, int n_t, int ld, int trans){
    
    //n_b - dimension of base vectors
    //n_b - dimension of vectors to be transformed
    //ld  - leading dimension of S
    
    set_zero_matr(O,n_t*(n_b+n_t));
    for (int i=0;i<n_t;i++) 
    for (int j=0;j<n_t;j++) 
        O[i*(n_b+n_t)+j+n_b]=U[i*n_t+j];
    
    double * T;
    if(trans){
        T=S+ld*n_b;
        for(int i=0;i<n_t;i++)
        for(int j=0;j<n_b;j++)
            O[i*(n_b+n_t)+j]=T[i*ld+j];
    }
    else{
        T=S+n_b;
//         printf("O:\n");
//         PrintMatr(O,n_t,n_b+n_t,1);
//         printf("S:\n");
//         PrintMatr(T,ld,ld,1);
        for(int i=0;i<n_t;i++)
        for(int j=0;j<n_b;j++)
            O[i*(n_b+n_t)+j]=T[j*ld+i];
    }
    // printf("O:\n");
    // PrintMatr(O,n_t,n_b+n_t,1);
    return 1;
    
}


int d_o_DM_gen(double * DM, double coef,
                double * L, int dim_l,
                double * R, int dim_r,
                double * sigma, int n_el){
    set_zero_matr(DM, dim_l*dim_r);

        
    
    for(int i=0;i<dim_l;i++)
    for(int j=0;j<dim_r;j++){
    for(int k=0;k<n_el ;k++)if(fabs(sigma[k])>svd_eps){
        DM[i*dim_r+j]+=coef*L[k*dim_l+i]*R[k*dim_r+j]/sigma[k];//printf("+ %e * %e / %e\n",L[k*dim_l+i],R[k*dim_r+j],sigma[k]);
    }//printf("+ %e * %e / %e\n",L[(n_el-1)*dim_l+i],R[(n_el-1)*dim_r+j],sigma[n_el-1]);
//     printf("DM(%d,%d)=%e\n",i,j,DM[i*dim_r+j]);
//     PrintMatr(DM,dim_l,dim_r,1);
    }
//     getchar();
    return 0;
}


int d_o_DMV_gen(double * DM,
                double * L, int dim_l,
                double * R, int dim_r,
                double * TH, int n_v , double * B){
    
    
//     double * B= new double[dim_l*n_v];
        
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                dim_l,n_v,n_v,1.0,
                L ,n_v,
                TH,n_v,0.0,
                B ,n_v);
        
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                dim_l,dim_r,n_v,1.0,
                B ,n_v,
                R ,n_v,0.0,
                DM,dim_r);
    
    
    return 0;
}


double E_active_4MO(double * AB, double * T2, double * MI, int n){
    
    double E=0;
    
    for(int i=0;i<n;i++)
    for(int j=0;j<n;j++)
    for(int k=0;k<n;k++)
    for(int l=0;l<n;l++){
        E+= AB[((i*n+j)*n+k)*n+l]
           *MI[((i*n+k)*n+j)*n+l];
               
        E+= T2[((i*n+j)*n+k)*n+l]
           *MI[((i*n+k)*n+j)*n+l]*0.25;//factor 1/4 must be used for full tensor; 1 must be used for symmetrized tensor
               
        E-= T2[((i*n+j)*n+k)*n+l]
           *MI[((i*n+l)*n+j)*n+k]*0.25;
    
    }
    
    return E;
    
}

inline int n_ext_calc(int *o1, int * o2, int n_v/*, int ** ext_ab*/){
    
    int res=0;
    int e1=0;
    int e2=0;
    for(int i=0;i<n_v;i++)
        if(o1[i]!=o2[i]){
//             ext_ab[o2[i]][///working here!!!!!!!
            res++;
        }
    
//     for(int k=0;k<n_v;k++)printf("%d",o1[k]);
//     printf("  ||||| ");
//     for(int k=0;k<n_v;k++)printf("%d",o2[k]);
//     printf("\n\nn = %d\n",res/2);
//     
    return res/2;
}

int make_doCI_matr(molecule * A, molecule * B, molecule * EF, int have_efrag, FILE * out_file, int dec){
    
    if(dec)if(A->n_frag!=2){
        fprintf(stdout,"WARNING: can not use doCI-wf decomposition (dec=1) for molecule with n_frag(%d) != 2\n",A->n_frag);
        fprintf(stdout,"         NOPT forced dec=0\n");
    }
    
    libint2::initialize();
    
    doCI_data D;
    
    D.first_alloc(A,B);
    D.gen_aldet_data();
    
    
    
        
    printf_timer("preparation");
    printf("\n\n");
    
    A->gen_1el_data();
    
    D.calc_S_MO(A);
    A->calc_H_AO();
    A->calc_D_AO();
    
    D.cor_svd(A);
    
    D.second_alloc();
    
    D.act_ort(A);
    
    if(dec==0){//full linking of fragment wf
        D.calc_S_SV(A);
        
        D.ci_ortogonalization(A->S_AO);
    }
    if(dec==1){//separate fragment wf-arrays with first-order decomposition
        D.link_and_transform_for_separate_spaces(A);
    }
//     if(dec!=1)if(dec!=0){
//        printf("ERROR: wf-decomposition for %d-th order is not realized. You can use:\n",dec);
//        printf("       a) dec=0 for full linking of fragment wf (can be expensive for large active spaces).\n");
//        printf("       b) dec=1 for approximate first order wf decomposition (can be used in the low overlap case).\n"); 
//     }
    
    printf_timer("ci_transformation");
       
    fflush(stdout);
    
    //2-electron engine
   
// #define READ_act_ints
// #define NO_2EL_DEBUG

// #define WRITE_act_ints
    
#ifndef NO_2EL_DEBUG
    calc_2el_CI(A,&D);
#endif
#ifdef WRITE_act_ints
    D.write_act_ints();
#endif
#ifdef READ_act_ints
     D.read_act_ints();
#endif
//     double a[2];
//     double b[2];
    
//     FILE * camm_par;
//     camm_par = fopen("par.txt","r");
//     char * tmp = new char[255];
//     fgets(tmp,255,camm_par);
// //     printf("%s\n", tmp);
//     sscanf(tmp,"%lf %lf %lf %lf",a,a+1,b,b+1);
    
//     fprintf(out_file, "reading external cam potential parameters\n");
//     fprintf(out_file, "V(1)=%f\nV(2)=%f\nc(1)=%f\nc(2)=%f\n\n",a[0],a[1],b[0],b[1]);
    
//     exit(1);
    
     printf_timer("before H_S_calc");
    if(dec==0){//full linking of fragment wf
        D.ci_link();
        
        D.S_calc();
        D.calc_1el_DM();
        
        calc_2e_ints(D.h_2e, D.n_act[0],D.n_alp_el[0],D.n_bet_el[0], D.aldet.coef[0], D.aldet.coef[1], D.act_INTS,0,D.n_st_total[0],D.n_st_total[0],D.n_st_total[1],D.n_st_total[1]);
        
        for(int i=0;i<D.n_st_total[0];i++)    
        for(int j=0;j<D.n_st_total[1];j++){
          
            double S = D.ACT_DET_A[i*D.n_st_total[1]+j];
        
            D.transform_1el_DM(i*D.n_st_total[1]+j);
            
            
            double H_1el = E_1el_calc(A->H_AO, D.DM_T, D.n_ao[0], D.n_ao[1])*D.p_SVD/*+camm_v*/;
            double  H_2el = D.H_2el_calc(S) + D.h_2e[i*D.n_st_total[1]+j]*D.p_SVD;
        
//             V_nuc = 0;
//             H_1el = 0;
            S=S*D.p_SVD;
            double d_x= E_1el_calc(A->Dx_AO, D.DM_T, D.n_ao[0], D.n_ao[1])*D.p_SVD+A->Dx_nuc*S;
            double d_y= E_1el_calc(A->Dy_AO, D.DM_T, D.n_ao[0], D.n_ao[1])*D.p_SVD+A->Dy_nuc*S;
            double d_z= E_1el_calc(A->Dz_AO, D.DM_T, D.n_ao[0], D.n_ao[1])*D.p_SVD+A->Dz_nuc*S;
        
            print_CI_results(out_file, i, j, S, A->V_nuc, H_1el, H_2el, d_x,d_y, d_z);
            
        }
    }
    
    if(dec==1){//separate fragment wf-arrays with first-order decomposition
        D.calc_CI_for_dec_wf();
        
        D.AO_to_act(D.H_ACT,A->H_AO);
        D.AO_to_act(D.Dx_ACT,A->Dx_AO);
        D.AO_to_act(D.Dy_ACT,A->Dy_AO);
        D.AO_to_act(D.Dz_ACT,A->Dz_AO);
        
        for(int i_n2=0;i_n2<D.n_ao[0]*D.n_ao[1];i_n2++)
            A->BUF[i_n2]=2*D.J[i_n2]-D.K[i_n2];
        
        double * JK_act = new double[D.n_act[0]*D.n_act[1]];
        D.AO_to_act(JK_act,A->BUF);
        
        double E_2el_core = E_1el_calc(D.J, D.DM_C, D.n_ao[0], D.n_ao[1])*D.p_SVD*2-
                            E_1el_calc(D.K, D.DM_C, D.n_ao[0], D.n_ao[1])*D.p_SVD;
        
        for(int a=0;a<A->n_states[0];a++)
        for(int b=0;b<A->n_states[1];b++)
        for(int c=0;c<B->n_states[0];c++)
        for(int d=0;d<B->n_states[1];d++){
            int abcd = ((a*A->n_states[1]+b)*B->n_states[0]+c)*B->n_states[1]+d;
            for(int i_n2=0;i_n2<D.n_act[0]*D.n_act[1];i_n2++)D.DM_T[i_n2]=(D.ACT_DMT_A[abcd*D.n_act[0]*D.n_act[1]+i_n2]+
                                                                           D.ACT_DMT_B[abcd*D.n_act[0]*D.n_act[1]+i_n2]);
            D.H_2el[abcd]+= E_2el_core*D.S_new[abcd]+E_1el_calc_2(JK_act, D.DM_T, D.n_act[0], D.n_act[1])*D.p_SVD;
            
            double H_1el = E_1el_calc_2(D.H_ACT, D.DM_T, D.n_act[0], D.n_act[1])*D.p_SVD +D.E_core*D.S_new[abcd]*2;
//             printf_timer("E_calc");
           
            double d_x= E_1el_calc_2(D.Dx_ACT, D.DM_T, D.n_act[0], D.n_act[1])*D.p_SVD-A->Dx_nuc*D.S_new[abcd]*D.p_SVD+D.Dx_core*D.S_new[abcd]*2;
            double d_y= E_1el_calc_2(D.Dy_ACT, D.DM_T, D.n_act[0], D.n_act[1])*D.p_SVD-A->Dy_nuc*D.S_new[abcd]*D.p_SVD+D.Dy_core*D.S_new[abcd]*2;
            double d_z= E_1el_calc_2(D.Dz_ACT, D.DM_T, D.n_act[0], D.n_act[1])*D.p_SVD-A->Dz_nuc*D.S_new[abcd]*D.p_SVD+D.Dz_core*D.S_new[abcd]*2;
//             printf_timer("D_calc");
            int i=a*A->n_states[1]+b;
            int j=c*B->n_states[1]+d;
            D.S_new[abcd] = D.S_new[abcd] * D.p_SVD;
            print_CI_results(out_file, i, j, D.S_new[abcd], A->V_nuc, H_1el, D.H_2el[abcd], d_x,d_y, d_z);
            
        }
        delete[] JK_act;
    }
    
    printf_timer("H and S calculation");
    
    D.clear();
    
    
    return 0;
}


int make_doCI_matr_with_1_order_decomp(molecule * A, molecule * B, molecule * EF, int have_efrag, FILE * out_file){
    
    libint2::initialize();
    
    doCI_data D;
    
    D.first_alloc(A,B);
//     D.gen_aldet_data();
    
    
    
    
    A->gen_1el_data();
    
     D.calc_S_MO(A);
    A->calc_H_AO();
    A->calc_D_AO();
    
    D.cor_svd(A);
    
    D.second_alloc();
    
    D.act_ort(A);
    
    D.link_and_transform_for_separate_spaces(A);
    printf_timer("ci_transformation");
    
// #define READ_act_ints
// #define NO_2EL_DEBUG

// #define WRITE_act_ints
    
#ifndef NO_2EL_DEBUG
    calc_2el_CI(A,&D);
#endif
#ifdef WRITE_act_ints
    D.write_act_ints();
#endif
#ifdef READ_act_ints
     D.read_act_ints();
#endif
    
    D.calc_CI_for_dec_wf();
    
    
    for(int a=0;a<A->n_states[0];a++)
    for(int b=0;b<A->n_states[1];b++)
    for(int c=0;c<B->n_states[0];c++)
    for(int d=0;d<B->n_states[1];d++){
        int abcd = ((a*A->n_states[1]+b)*B->n_states[0]+c)*B->n_states[1]+d;
//         ac   =   a*B->n_states[0]+c;
//         bd   =   b*B->n_states[1]+d;        
        if(D.n_act[0])
        d_o_DMV_gen(D.DM_A,
                    D.L_ACT, D.n_ao[0],
                    D.R_ACT, D.n_ao[1],
                    D.ACT_DMT_A+abcd*D.n_act[0]*D.n_act[1], D.n_act[0],D.BUF);
//         printf_timer("DMA2");
        if(D.n_act[0])
        d_o_DMV_gen(D.DM_B,
                    D.L_ACT, D.n_ao[0],
                    D.R_ACT, D.n_ao[1],
                    D.ACT_DMT_B+abcd*D.n_act[0]*D.n_act[1], D.n_act[0],D.BUF);
                
        for(int i_n2=0;i_n2<D.n_ao[0]*D.n_ao[1];i_n2++)D.DM_T[i_n2]=D.DM_C[i_n2]*D.S_new[abcd]*2+(D.DM_B[i_n2]+D.DM_A[i_n2]);
        
        D.H_2el[abcd]+= E_1el_calc(D.J, D.DM_C, D.n_ao[0], D.n_ao[1])*D.S_new[abcd]*D.p_SVD*2+
                        E_1el_calc(D.J, D.DM_A,D.n_ao[0], D.n_ao[1])*D.p_SVD*2+
                        E_1el_calc(D.J, D.DM_B,D.n_ao[0], D.n_ao[1])*D.p_SVD*2-
                        E_1el_calc(D.K, D.DM_C, D.n_ao[0], D.n_ao[1])*D.S_new[abcd]*D.p_SVD-
                        E_1el_calc(D.K, D.DM_A,D.n_ao[0], D.n_ao[1])*D.p_SVD-
                        E_1el_calc(D.K, D.DM_B,D.n_ao[0], D.n_ao[1])*D.p_SVD;
        double H_1el = E_1el_calc(A->H_AO, D.DM_T, D.n_ao[0], D.n_ao[1])*D.p_SVD;
        
        D.S_new[abcd] = D.S_new[abcd] * D.p_SVD;
        double d_x= E_1el_calc(A->Dx_AO, D.DM_T, D.n_ao[0], D.n_ao[1])*D.p_SVD-A->Dx_nuc*D.S_new[abcd];
        double d_y= E_1el_calc(A->Dy_AO, D.DM_T, D.n_ao[0], D.n_ao[1])*D.p_SVD-A->Dy_nuc*D.S_new[abcd];
        double d_z= E_1el_calc(A->Dz_AO, D.DM_T, D.n_ao[0], D.n_ao[1])*D.p_SVD-A->Dz_nuc*D.S_new[abcd];
        
        int i=a*A->n_states[1]+b;
        int j=c*B->n_states[1]+d;
        print_CI_results(out_file, i, j, D.S_new[abcd], A->V_nuc, H_1el, D.H_2el[abcd], d_x,d_y, d_z);
        
    }
        printf_timer("H and S calculation");

    return 0;
}

int Sorb_calc(molecule * A, molecule * Ap){
    
    int calc_type=0;
    
    if(A->n_frag!=1){
       printf("ERROR: n_frag=%d is impossible in S_orb calculation\n",A->n_frag);
       printf("       check that you use n_mol=1 in 1-st molecule defined by $MOL group\n");
       exit(1);
    }
    if(Ap->n_frag!=1){
       printf("ERROR: n_frag=%d is impossible in S_orb calculation\n",Ap->n_frag);
       printf("       check that you use n_mol=1 in 2-nd molecule defined by $MOL group\n");
       exit(1);
    }
    if((A->n_act_el_alp[0]-Ap->n_act_el_alp[0])==1)
    if((A->n_act_el_bet[0]-Ap->n_act_el_bet[0])==0){
        calc_type=1;
        printf("First molecule has 1 extra alpha-electron\n");
        printf("Continue calculations with type A\n\n");
    }
    
    if((A->n_act_el_alp[0]-Ap->n_act_el_alp[0])==0)
    if((A->n_act_el_bet[0]-Ap->n_act_el_bet[0])==1){
        calc_type=2;
        printf("First molecule has 1 extra beta-electron\n");
        printf("Continue calculations with type B\n\n");
    }
    
    
    if(calc_type==0){
       printf("ERROR: wrong number of electrons\n");
       printf("       n_alp_el(A ) = %d\n",A ->n_act_el_alp[0]);
       printf("       n_alp_el(A+) = %d\n",Ap->n_act_el_alp[0]);
       printf("       n_bet_el(A ) = %d\n",A ->n_act_el_bet[0]);
       printf("       n_bet_el(A+) = %d\n",Ap->n_act_el_bet[0]);
       printf("       first molecule must have 1 extra electron\n");
       exit(1);
    }
    
    
    libint2::initialize();
    
    doCI_data D;
    
    
    D.first_alloc(A,Ap);
    D.two_charge_states=1;
    
    D.gen_aldet_data();
    
    
    
        
    printf_timer("preparation");
    printf("\n\n");
    
    A->gen_1el_data();
    D.calc_S_MO(A);
    
    D.cor_svd(A);
    
    D.second_alloc();
    
    D.act_ort(A);

    D.calc_S_SV(A);
    
    D.ci_ortogonalization(A->S_AO);
    
    printf_timer("ci_transformation");
    
    D.ci_link();
    
    int n_smo = A->n_states[0]*Ap->n_states[0];
    double * S_ORB   = new double[A->n_act_orb[0]*n_smo];
    double * S_PRINT = new double[A->n_ao        *n_smo];
    
    if(calc_type==1){
        printf("type A calculation is yet disabled\n");
        exit(1);
    }
    
    if(calc_type==2){
        aldet_gen_1_el_vec_arr(S_ORB ,&(D.aldet  ), 0, A->n_states[0], 
                                      &(D.aldet_2), 0, Ap->n_states[0], 
                                      0, A->n_act_orb[0], A->n_act_orb[0], 1);
        
//         for(int i=0;i<n_smo;i++)
//             PrintMatr(S_ORB+i*A->n_act_orb[0],1,A->n_act_orb[0],0);
        
        
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                     n_smo,A->n_ao,A->n_act_orb[0],1.0,
                     S_ORB,A->n_act_orb[0],
                     D.L_ACT,A->n_act_orb[0],0.0,
                     S_PRINT,A->n_ao);
        
//         for(int i=0;i<n_smo;i++)
//             PrintMatr(S_PRINT+i*A->n_ao,1,A->n_ao,0);
        
        
//         exit(1);
    }
    delete[] A->MO_VEC;
    A->MO_VEC = S_PRINT;
    A->n_mo=n_smo;
    
    A->calc_S_AO();
    
    double * BUF  = new double[A->n_ao*n_smo];
    double * S_MO = new double[n_smo*n_smo];
    
//     gen_matr_from_2shells(S_AO,s[1],s[1],A->n_ao,A->n_ao,&(s_engine),0);
//     printf("S_AO:\n");
//     PrintMatr(A->S_AO,A->basis[0].n_ao,A->basis[0].n_ao,1);
    
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        A->n_ao,n_smo,A->n_ao,1.0,
                        A->S_AO,A->n_ao,
                        A->MO_VEC,A->n_ao,0.0,
                        BUF,n_smo);
    
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                        n_smo,n_smo,A->n_ao,1.0,
                        A->MO_VEC,A->n_ao,
                        BUF,n_smo,0.0,
                        S_MO,n_smo);
    printf("S:\n");
    PrintMatr(S_MO, n_smo, n_smo,0);
    
    delete[] A->orb_energy;
    A->orb_energy = new double[n_smo];
    for(int i=0; i < n_smo; i++)
        A->orb_energy[i]=S_MO[i*(n_smo+1)];
//     printf("S_diag:\n");
//     PrintMatr(A->orb_energy[0],1,n_smo,0);
    
    D.clear();
    
    
    delete[] BUF ;
    delete[] S_MO;

    return 1;
}

double gen_doCI_DM(double ** DM, double * S_MO, molecule * A, molecule * B, int n_ao,double ** C, int index1, int index2, int n_s, int dec){
    
    double S_total=0;
    
    libint2::initialize();
    
    doCI_data D;
    
    D.first_alloc(A,B);
    D.gen_aldet_data();
    
    
    
    printf_timer("preparation");
    printf("\n\n");
    
    A->gen_1el_data();
    A->calc_S_AO();
    memcpy(S_MO, A->S_AO, A->n_ao*A->n_ao*sizeof(double));
    
    D.calc_S_MO(A);
    
    D.cor_svd(A);
    
    D.second_alloc();
    
    D.act_ort(A);
    
//     D.calc_S_SV(A);
//     D.ci_ortogonalization(A->S_AO);
    
    double * C1;
    double * C2;
    double * DM_ij;
    
    if(dec==0){//full linking of fragment wf
        D.calc_S_SV(A);
        
        D.ci_ortogonalization(A->S_AO);
                    
        printf_timer("ci_transformation");
           
        fflush(stdout);
        
        D.ci_link();
        
        D.S_calc();
        
        D.calc_1el_DM();
        
        for(int i=0;i<D.n_st_total[0];i++)    
        for(int j=0;j<D.n_st_total[1];j++){
            
            double S = D.ACT_DET_A[i*D.n_st_total[1]+j];
        
            D.transform_1el_DM(i*D.n_st_total[1]+j);
            for(int i_s=0;i_s<n_s;i_s++)
            for(int j_s=0;j_s<n_s;j_s++){
                C1 =C[i_s]+index1;
                C2 =C[j_s]+index2;
                DM_ij = DM[i_s*n_s+j_s];
                for(int i_n2=0;i_n2<D.n_ao[0]*D.n_ao[1];i_n2++)DM_ij[i_n2]+=D.DM_T[i_n2]*C1[i]*C2[j]*D.p_SVD;
            }
            
        }
        printf("\n");
    
    }
    if(dec==1){//separate fragment wf-arrays with first-order decomposition
        D.link_and_transform_for_separate_spaces(A);
        D.calc_DM_for_dec_wf(/*num_threads*/);
        
        
        for(int a=0;a<A->n_states[0];a++)
        for(int b=0;b<A->n_states[1];b++)
        for(int c=0;c<B->n_states[0];c++)
        for(int d=0;d<B->n_states[1];d++){
            int abcd = ((a*A->n_states[1]+b)*B->n_states[0]+c)*B->n_states[1]+d;
            if(D.n_act[0])
            d_o_DMV_gen(D.DM_A,
                        D.L_ACT, D.n_ao[0],
                        D.R_ACT, D.n_ao[1],
                        D.ACT_DMT_A+abcd*D.n_act[0]*D.n_act[1], D.n_act[0],D.BUF);
            if(D.n_act[0])
            d_o_DMV_gen(D.DM_B,
                        D.L_ACT, D.n_ao[0],
                        D.R_ACT, D.n_ao[1],
                        D.ACT_DMT_B+abcd*D.n_act[0]*D.n_act[1], D.n_act[0],D.BUF);
                    
            for(int i_n2=0;i_n2<D.n_ao[0]*D.n_ao[1];i_n2++)D.DM_T[i_n2]=D.DM_C[i_n2]*D.S_new[abcd]*2+(D.DM_B[i_n2]+D.DM_A[i_n2]);
            
            D.S_new[abcd] = D.S_new[abcd] * D.p_SVD;
            
            int i=a*A->n_states[1]+b;
            int j=c*B->n_states[1]+d;
            for(int i_s=0;i_s<n_s;i_s++)
            for(int j_s=0;j_s<n_s;j_s++){
                C1 =C[i_s]+index1;
                C2 =C[j_s]+index2;
                DM_ij = DM[i_s*n_s+j_s];
                for(int i_n2=0;i_n2<D.n_ao[0]*D.n_ao[1];i_n2++)DM_ij[i_n2]+=D.DM_T[i_n2]*C1[i]*C2[j]*D.p_SVD;
//                 for(int i_n2=0;i_n2<D.n_ao[0]*D.n_ao[1];i_n2++)DM_ij[i_n2]+=D.DM_T[i_n2]*C1[i]*C2[j]*D.p_SVD;
            }
            
            
        }
        
    }
    

//     printf("start cleaning D\n");
    D.clear();
        
    

    
    return S_total;
}

double gen_doCI_DM_1block(double ** DM, double * S_WF, double * S_MO, molecule * A, molecule * B, int dec){
    
    double S_total=0;
    
    libint2::initialize();
    
    doCI_data D;
    
    D.first_alloc(A,B);
    D.gen_aldet_data();
    
    
    
    printf_timer("preparation");
    printf("\n\n");
    
    A->gen_1el_data();
    
    A->calc_S_AO();
    memcpy(S_MO, A->S_AO, A->n_ao*A->n_ao*sizeof(double));
    
    D.calc_S_MO(A);
    
    D.cor_svd(A);
    
    D.second_alloc();
    
    D.act_ort(A);
    
    double * DM_ij;
    
    if(dec==0){//full linking of fragment wf
        D.calc_S_SV(A);
        
        D.ci_ortogonalization(A->S_AO);
                    
        printf_timer("ci_transformation");
           
        fflush(stdout);
        
        D.ci_link();
        
        D.S_calc();
        
        D.calc_1el_DM();
        
        for(int i=0;i<D.n_st_total[0];i++)    
        for(int j=0;j<D.n_st_total[1];j++){
            
            double S = D.ACT_DET_A[i*D.n_st_total[1]+j];
            
            S_WF[i*D.n_st_total[1]+j]=S*D.p_SVD;
        
            D.transform_1el_DM(i*D.n_st_total[1]+j);
            DM_ij = DM[i*D.n_st_total[1]+j];
            for(int i_n2=0;i_n2<D.n_ao[0]*D.n_ao[1];i_n2++)DM_ij[i_n2]+=D.DM_T[i_n2]*D.p_SVD;
            
            
        }
        printf("\n");
    
    }
    if(dec==1){//separate fragment wf-arrays with first-order decomposition
        D.link_and_transform_for_separate_spaces(A);
        D.calc_DM_for_dec_wf(/*num_threads*/);
        
        
        for(int a=0;a<A->n_states[0];a++)
        for(int b=0;b<A->n_states[1];b++)
        for(int c=0;c<B->n_states[0];c++)
        for(int d=0;d<B->n_states[1];d++){
            int abcd = ((a*A->n_states[1]+b)*B->n_states[0]+c)*B->n_states[1]+d;
            if(D.n_act[0])
            d_o_DMV_gen(D.DM_A,
                        D.L_ACT, D.n_ao[0],
                        D.R_ACT, D.n_ao[1],
                        D.ACT_DMT_A+abcd*D.n_act[0]*D.n_act[1], D.n_act[0],D.BUF);
            if(D.n_act[0])
            d_o_DMV_gen(D.DM_B,
                        D.L_ACT, D.n_ao[0],
                        D.R_ACT, D.n_ao[1],
                        D.ACT_DMT_B+abcd*D.n_act[0]*D.n_act[1], D.n_act[0],D.BUF);
                    
            for(int i_n2=0;i_n2<D.n_ao[0]*D.n_ao[1];i_n2++)D.DM_T[i_n2]=D.DM_C[i_n2]*D.S_new[abcd]*2+(D.DM_B[i_n2]+D.DM_A[i_n2]);
            
            D.S_new[abcd] = D.S_new[abcd] * D.p_SVD;
            
            
            int i=a*A->n_states[1]+b;
            int j=c*B->n_states[1]+d;
            
            S_WF[i*D.n_st_total[1]+j]=D.S_new[abcd];
            
            DM_ij = DM[i*D.n_st_total[1]+j];
            
            for(int i_n2=0;i_n2<D.n_ao[0]*D.n_ao[1];i_n2++)DM_ij[i_n2]+=D.DM_T[i_n2]*D.p_SVD;
            
        }
        
    }
    

//     printf("start cleaning D\n");
    D.clear();
        
    

    
    return S_total;
}

double gen_doCI_DM_1block_PT(double ** DM, double * S_WF, double * S_MO, molecule * A, molecule * B, int dec){
    
    double S_total=0;
    
    libint2::initialize();
    
    doCI_data D;
    
    D.first_alloc(A,B);
    D.gen_aldet_data();
    
    
    
    printf_timer("preparation");
    printf("\n\n");
    
    A->gen_1el_data();
    
    A->calc_S_AO();
    memcpy(S_MO, A->S_AO, A->n_ao*A->n_ao*sizeof(double));
    
    D.calc_S_MO(A);
    
    D.cor_svd_PT(A);
    
    D.second_alloc();
    
    
    double * DM_ij;
    
    D.link_and_transform_for_PT(A);
    D.calc_DM_for_PT(/*num_threads*/);
    
    
    double * DM_A_1st_L= new double[D.n_ao[0]*D.n_ao[1]];
    double * DM_B_1st_L= new double[D.n_ao[0]*D.n_ao[1]];
    double * DM_A_1st_R= new double[D.n_ao[0]*D.n_ao[1]];
    double * DM_B_1st_R= new double[D.n_ao[0]*D.n_ao[1]];
        
    for(int a=0;a<A->n_states[0];a++)
    for(int b=0;b<A->n_states[1];b++)
    for(int c=0;c<B->n_states[0];c++)
    for(int d=0;d<B->n_states[1];d++){
        int abcd = ((a*A->n_states[1]+b)*B->n_states[0]+c)*B->n_states[1]+d;
        if(D.n_act[0])
        d_o_DMV_gen(D.DM_A,
                    D.L_ACT, D.n_ao[0],
                    D.R_ACT, D.n_ao[1],
                    D.ACT_DMT_A+abcd*D.n_act[0]*D.n_act[1], D.n_act[0],D.BUF);
        if(D.n_act[0])
        d_o_DMV_gen(D.DM_B,
                    D.L_ACT, D.n_ao[0],
                    D.R_ACT, D.n_ao[1],
                    D.ACT_DMT_B+abcd*D.n_act[0]*D.n_act[1], D.n_act[0],D.BUF);
                
        d_o_DMV_gen(DM_A_1st_L,
                    D.L_ACT_1st_order, D.n_ao[0],
                    D.R_ACT, D.n_ao[1],
                    D.ACT_DMT_A_1st_order+abcd*D.n_act[0]*D.n_act[1], D.n_act[0],D.BUF);
        
        d_o_DMV_gen(DM_B_1st_L,
                    D.L_ACT_1st_order, D.n_ao[0],
                    D.R_ACT, D.n_ao[1],
                    D.ACT_DMT_B_1st_order+abcd*D.n_act[0]*D.n_act[1], D.n_act[0],D.BUF);

        d_o_DMV_gen(DM_A_1st_R,
                    D.L_ACT, D.n_ao[0],
                    D.R_ACT_1st_order, D.n_ao[1],
                    D.ACT_DMT_A_1st_order+abcd*D.n_act[0]*D.n_act[1], D.n_act[0],D.BUF);
        
        d_o_DMV_gen(DM_B_1st_R,
                    D.L_ACT, D.n_ao[0],
                    D.R_ACT_1st_order, D.n_ao[1],
                    D.ACT_DMT_B_1st_order+abcd*D.n_act[0]*D.n_act[1], D.n_act[0],D.BUF);

        
        
        for(int i_n2=0;i_n2<D.n_ao[0]*D.n_ao[1];i_n2++)
            D.DM_T[i_n2]=D.DM_C[i_n2]*D.S_new[abcd]*2+(D.DM_B[i_n2]+D.DM_A[i_n2])
                                                     +(DM_B_1st_L[i_n2]+DM_A_1st_L[i_n2])
                                                     +(DM_B_1st_R[i_n2]+DM_A_1st_R[i_n2]);
        
        
        
        D.S_new[abcd] = D.S_new[abcd] * D.p_SVD;
        
        
        int i=a*A->n_states[1]+b;
        int j=c*B->n_states[1]+d;
        
        S_WF[i*D.n_st_total[1]+j]=D.S_new[abcd];
        
        DM_ij = DM[i*D.n_st_total[1]+j];
        
        for(int i_n2=0;i_n2<D.n_ao[0]*D.n_ao[1];i_n2++)DM_ij[i_n2]+=D.DM_T[i_n2]*D.p_SVD;
        
    }
        
    
    

//     printf("start cleaning D\n");
    D.clear();
        
    delete[] DM_A_1st_L;
    delete[] DM_B_1st_L;
    delete[] DM_A_1st_R;
    delete[] DM_B_1st_R;
    

    
    return S_total;
}


int make_doCI_Smatr(molecule * A, molecule * B, molecule * EF, int have_efrag, FILE * out_file, int dec){
    
    libint2::initialize();
    
    doCI_data D;
    
    D.first_alloc(A,B);
    D.gen_aldet_data();
    
    
    
        
    printf_timer("preparation");
    printf("\n\n");
    
    A->gen_1el_data();
    
    
    D.calc_S_MO(A);
    
    D.cor_svd(A);
    
    D.second_alloc();
    
    D.act_ort(A);
    
    if(dec==0){//full linking of fragment wf
        D.calc_S_SV(A);
        
        D.ci_ortogonalization(A->S_AO);
    }
    if(dec==1){//separate fragment wf-arrays with first-order decomposition
        D.link_and_transform_for_separate_spaces(A);
    }
//     if(dec!=1)if(dec!=0){
//        printf("ERROR: wf-decomposition for %d-th order is not realized. You can use:\n",dec);
//        printf("       a) dec=0 for full linking of fragment wf (can be expensive for large active spaces).\n");
//        printf("       b) dec=1 for approximate first order wf decomposition (can be used in the low overlap case).\n"); 
//     }
    
    printf_timer("ci_transformation");
       
    fflush(stdout);
    
    //2-electron engine
   
// #define READ_act_ints
// #define NO_2EL_DEBUG

// #define WRITE_act_ints
    if(dec==0){//full linking of fragment wf
        D.ci_link();
        printf("not tested!!!!\n");
        D.S_calc();
//         D.calc_1el_DM();
        
//         calc_2e_ints(D.h_2e, D.n_act[0],D.n_alp_el[0],D.n_bet_el[0], D.aldet.coef[0], D.aldet.coef[1], D.act_INTS,0,D.n_st_total[0],D.n_st_total[0],D.n_st_total[1],D.n_st_total[1]);
        
        for(int i=0;i<D.n_st_total[0];i++)    
        for(int j=0;j<D.n_st_total[1];j++){
          
            double S = D.ACT_DET_A[i*D.n_st_total[1]+j];
        
            S=S*D.p_SVD;
            print_CI_results(out_file, i, j, S, 0, 0, 0, 0, 0, 0);
            printf_timer("S calculation");
            printf("\n");
            fflush(stdout);
        
        
        }
    }
    
    if(dec==1){//separate fragment wf-arrays with first-order decomposition
        D.calc_S_for_dec_wf(/*num_threads*/);
        
        
        for(int a=0;a<A->n_states[0];a++)
        for(int b=0;b<A->n_states[1];b++)
        for(int c=0;c<B->n_states[0];c++)
        for(int d=0;d<B->n_states[1];d++){
            int abcd = ((a*A->n_states[1]+b)*B->n_states[0]+c)*B->n_states[1]+d;
//             ac   =   a*B->n_states[0]+c;
//             bd   =   b*B->n_states[1]+d;        
            
            D.S_new[abcd] = D.S_new[abcd] * D.p_SVD;
            int i=a*A->n_states[1]+b;
            int j=c*B->n_states[1]+d;
            print_CI_results(out_file, i, j, D.S_new[abcd], 0, 0, 0, 0, 0, 0);
            
        }
        printf_timer("H and S calculation");
    }
    
    
    D.clear();
    
    
    return 0;
}

int cut_AA_matr(double * Out, double * Inp, int d_a, int d_b, int v_a, int v_b, int n_a, int n_b,int ld){
    
    set_zero_matr(Out, n_a*n_a);
    
    for(int i=0; i<d_a;i++){
        for(int j=0; j<d_a;j++)
            Out[i*n_a+j]=Inp[i*ld+j];
        for(int j=0; j<v_a;j++)
            Out[i*n_a+j+d_a]=Inp[i*ld+j+d_a+d_b];
        for(int j=0; j<n_a-v_a-d_a;j++)
            Out[i*n_a+j+d_a+v_a]=Inp[i*ld+j+d_a+d_b+v_a+v_b];
    }
    for(int i=0; i<v_a;i++){
        for(int j=0; j<d_a;j++)
            Out[(i+d_a)*n_a+j]=Inp[(i+d_a+d_b)*ld+j];
        for(int j=0; j<v_a;j++)
            Out[(i+d_a)*n_a+j+d_a]=Inp[(i+d_a+d_b)*ld+j+d_a+d_b];
        for(int j=0; j<n_a-v_a-d_a;j++)
            Out[(i+d_a)*n_a+j+d_a+v_a]=Inp[(i+d_a+d_b)*ld+j+d_a+d_b+v_a+v_b];
    }
    for(int i=0; i<n_a-v_a-d_a;i++){
        for(int j=0; j<d_a;j++)
            Out[(i+d_a+v_a)*n_a+j]=Inp[(i+d_a+d_b+v_a+v_b)*ld+j];
        for(int j=0; j<v_a;j++)
            Out[(i+d_a+v_a)*n_a+j+d_a]=Inp[(i+d_a+d_b+v_a+v_b)*ld+j+d_a+d_b];
        for(int j=0; j<n_a-v_a-d_a;j++)
            Out[(i+d_a+v_a)*n_a+j+d_a+v_a]=Inp[(i+d_a+d_b+v_a+v_b)*ld+j+d_a+d_b+v_a+v_b];
    }
    
    return 0;
}

int cut_AB_matr(double * Out, double * Inp, int d_a, int d_b, int v_a, int v_b, int n_a, int n_b,int ld){
    
    set_zero_matr(Out, n_a*n_b);
    
    for(int i=0; i<d_a;i++){
        for(int j=0; j<d_b;j++)
            Out[i*n_b+j]=Inp[i*ld+j+d_a];
        for(int j=0; j<v_b;j++)
            Out[i*n_b+j+d_b]=Inp[i*ld+j+d_a+d_b+v_a];
        for(int j=0; j<n_b-v_b-d_b;j++)
            Out[i*n_b+j+d_b+v_b]=Inp[i*ld+j+ld/2+d_b+v_b];
    }
    for(int i=0; i<v_a;i++){
        for(int j=0; j<d_b;j++)
            Out[(i+d_a)*n_b+j]=Inp[(i+d_a+d_b)*ld+j+d_a];
        for(int j=0; j<v_b;j++)
            Out[(i+d_a)*n_b+j+d_b]=Inp[(i+d_a+d_b)*ld+j+d_a+d_b+v_a];
        for(int j=0; j<n_b-v_b-d_b;j++)
            Out[(i+d_a)*n_b+j+d_b+v_b]=Inp[(i+d_a+d_b)*ld+j+ld/2+d_b+v_b];
    }
    for(int i=0; i<n_a-v_a-d_a;i++){
        for(int j=0; j<d_b;j++)
            Out[(i+d_a+v_a)*n_b+j]=Inp[(i+d_a+d_b+v_a+v_b)*ld+j+d_a];
        for(int j=0; j<v_b;j++)
            Out[(i+d_a+v_a)*n_b+j+d_b]=Inp[(i+d_a+d_b+v_a+v_b)*ld+j+d_a+d_b+v_a];
        for(int j=0; j<n_b-v_b-d_b;j++)
            Out[(i+d_a+v_a)*n_b+j+d_b+v_b]=Inp[(i+d_a+d_b+v_a+v_b)*ld+j+ld/2+d_b+v_b];
    }
    
    return 0;
}

int cut_BB_matr(double * Out, double * Inp, int d_a, int d_b, int v_a, int v_b, int n_a, int n_b,int ld){
    
    set_zero_matr(Out, n_b*n_b);
    
    for(int i=0; i<d_b;i++){
        for(int j=0; j<d_b;j++)
            Out[i*n_b+j]=Inp[(i+d_a)*ld+j+d_a];
        for(int j=0; j<v_b;j++)
            Out[i*n_b+j+d_b]=Inp[(i+d_a)*ld+j+d_a+d_b+v_a];
        for(int j=0; j<n_b-v_b-d_b;j++)
            Out[i*n_b+j+d_b+v_b]=Inp[(i+d_a)*ld+j+ld/2+d_b+v_b];
    }
    for(int i=0; i<v_b;i++){
        for(int j=0; j<d_b;j++)
            Out[(i+d_b)*n_b+j]=Inp[(i+d_a+d_b+v_a)*ld+j+d_a];
        for(int j=0; j<v_b;j++)
            Out[(i+d_b)*n_b+j+d_b]=Inp[(i+d_a+d_b+v_a)*ld+j+d_a+d_b+v_a];
        for(int j=0; j<n_b-v_b-d_b;j++)
            Out[(i+d_b)*n_b+j+d_b+v_b]=Inp[(i+d_a+d_b+v_a)*ld+j+ld/2+d_b+v_b];
    }
    for(int i=0; i<n_b-v_b-d_b;i++){
        for(int j=0; j<d_b;j++)
            Out[(i+d_b+v_b)*n_b+j]=Inp[(i+ld/2+d_b+v_b)*ld+j+d_a];
        for(int j=0; j<v_b;j++)
            Out[(i+d_b+v_b)*n_b+j+d_b]=Inp[(i+ld/2+d_b+v_b)*ld+j+d_a+d_b+v_a];
        for(int j=0; j<n_b-v_b-d_b;j++)
            Out[(i+d_b+v_b)*n_b+j+d_b+v_b]=Inp[(i+ld/2+d_b+v_b)*ld+j+ld/2+d_b+v_b];
    }
    
    return 0;
}

int vec_reordr(double * Out, double * Inp, int d_a, int d_b, int v_a, int v_b, int n_a, int n_b,int ld){
    
    set_zero_matr(Out, n_a*n_a);
    
    for(int i=0; i<d_a;i++){
        for(int j=0; j<ld;j++)
            Out[i*ld+j]=Inp[i*ld+j];
    }
    for(int i=0; i<v_a;i++){
        for(int j=0; j<ld;j++)
            Out[(i+d_a)*ld+j]=Inp[(i+d_a+d_b)*ld+j];
    }
    for(int i=0; i<n_a-v_a-d_a;i++){
        for(int j=0; j<ld;j++)
            Out[(i+d_a+v_a)*ld+j]=Inp[(i+d_a+d_b+v_a+v_b)*ld+j];
    }
    for(int i=0; i<d_b;i++){
        for(int j=0; j<ld;j++)
            Out[(i+n_a)*ld+j]=Inp[(i+d_a)*ld+j];
    }
    for(int i=0; i<v_b;i++){
        for(int j=0; j<ld;j++)
            Out[(i+n_a+d_b)*ld+j]=Inp[(i+d_a+d_b+v_a)*ld+j];
    }
    for(int i=0; i<n_b-v_b-d_b;i++){
        for(int j=0; j<ld;j++)
            Out[(i+n_a+d_b+v_b)*ld+j]=Inp[(i+ld/2+d_b+v_b)*ld+j];
    }
    
    return 0;
}



int make_doCI_H1el(molecule * A, molecule * B, molecule * EF, int have_efrag){
    
    libint2::initialize();
    
    doCI_data D;
    
    D.first_alloc(A,B);
    D.gen_aldet_data();
    
    
    
        
    printf_timer("preparation");
    printf("\n\n");
    
    A->gen_1el_data();
    
    
    D.calc_S_MO(A);
    
    A->calc_H_AO();

    A->write_H_AO();
    
    return 0;
}

int make_doCI_PT(molecule * A, molecule * B, molecule * EF, int have_efrag, FILE * out_file, int dec, int nshA, int nshB){
    
    libint2::initialize();
    
    doCI_data D;
    
    D.first_alloc(A,B);
    D.gen_aldet_data();
    
    
    
        
    printf_timer("preparation");
    printf("\n\n");
    
    A->gen_1el_data();
    
     D.calc_S_MO(A);
    A->calc_H_AO();
    A->calc_D_AO();
    
    D.cor_svd_PT(A);
    
    D.second_alloc();
    
//     D.act_ort(A);
    D.link_and_transform_for_PT(A);
    
   printf_timer("ci_transformation");
       
    fflush(stdout);
    
    //2-electron engine
   
// #define READ_act_ints
// #define NO_2EL_DEBUG

// #define WRITE_act_ints
    
// #ifndef NO_2EL_DEBUG
//     calc_2el_CI(A,&D);
    
    calc_2el_CI_1st_order(A,&D, nshA, nshB);
// #endif
// #ifdef WRITE_act_ints
//     D.write_act_ints();
// #endif
// #ifdef READ_act_ints
//      D.read_act_ints();
// #endif
    
    D.calc_CI_for_PT();
        
//     exit(1);
    
    double * DM_A_1st_L= new double[D.n_ao[0]*D.n_ao[1]];
    double * DM_B_1st_L= new double[D.n_ao[0]*D.n_ao[1]];
    double * DM_A_1st_R= new double[D.n_ao[0]*D.n_ao[1]];
    double * DM_B_1st_R= new double[D.n_ao[0]*D.n_ao[1]];
    
    double * K_1st = D.K_L_1st;
    
//     printf("<%d,%d+%d,%d|%d,%d+%d,%d>\n",D.n_alp_el_f[0][0],D.n_bet_el_f[0][0]
//                                         ,D.n_alp_el_f[0][1],D.n_bet_el_f[0][1]
//                                         ,D.n_alp_el_f[1][0],D.n_bet_el_f[1][0]
//                                         ,D.n_alp_el_f[1][1],D.n_bet_el_f[1][1]);
    if((D.n_alp_el_f[0][0]+D.n_bet_el_f[0][0])==(D.n_alp_el_f[1][0]+D.n_bet_el_f[1][0]+1)){
        K_1st = D.K_R_1st;
//         printf("K changed!!\n");
    }
//     exit(0);
    for(int a=0;a<A->n_states[0];a++)
    for(int b=0;b<A->n_states[1];b++)
    for(int c=0;c<B->n_states[0];c++)
    for(int d=0;d<B->n_states[1];d++){
        int abcd = ((a*A->n_states[1]+b)*B->n_states[0]+c)*B->n_states[1]+d;
//         ac   =   a*B->n_states[0]+c;
//         bd   =   b*B->n_states[1]+d;        
        if(D.n_act[0])
        d_o_DMV_gen(D.DM_A,
                    D.L_ACT, D.n_ao[0],
                    D.R_ACT, D.n_ao[1],
                    D.ACT_DMT_A+abcd*D.n_act[0]*D.n_act[1], D.n_act[0],D.BUF);

        if(D.n_act[0])
        d_o_DMV_gen(D.DM_B,
                    D.L_ACT, D.n_ao[0],
                    D.R_ACT, D.n_ao[1],
                    D.ACT_DMT_B+abcd*D.n_act[0]*D.n_act[1], D.n_act[0],D.BUF);

        d_o_DMV_gen(DM_A_1st_L,
                    D.L_ACT_1st_order, D.n_ao[0],
                    D.R_ACT, D.n_ao[1],
                    D.ACT_DMT_A_1st_order+abcd*D.n_act[0]*D.n_act[1], D.n_act[0],D.BUF);
        
        d_o_DMV_gen(DM_B_1st_L,
                    D.L_ACT_1st_order, D.n_ao[0],
                    D.R_ACT, D.n_ao[1],
                    D.ACT_DMT_B_1st_order+abcd*D.n_act[0]*D.n_act[1], D.n_act[0],D.BUF);

        d_o_DMV_gen(DM_A_1st_R,
                    D.L_ACT, D.n_ao[0],
                    D.R_ACT_1st_order, D.n_ao[1],
                    D.ACT_DMT_A_1st_order+abcd*D.n_act[0]*D.n_act[1], D.n_act[0],D.BUF);
        
        d_o_DMV_gen(DM_B_1st_R,
                    D.L_ACT, D.n_ao[0],
                    D.R_ACT_1st_order, D.n_ao[1],
                    D.ACT_DMT_B_1st_order+abcd*D.n_act[0]*D.n_act[1], D.n_act[0],D.BUF);

        
        
        for(int i_n2=0;i_n2<D.n_ao[0]*D.n_ao[1];i_n2++)
            D.DM_T[i_n2]=D.DM_C[i_n2]*D.S_new[abcd]*2+(D.DM_B[i_n2]+D.DM_A[i_n2])
                                                     +(DM_B_1st_L[i_n2]+DM_A_1st_L[i_n2])
                                                     +(DM_B_1st_R[i_n2]+DM_A_1st_R[i_n2]);
        
        D.H_2el[abcd]+= E_1el_calc(D.J, D.DM_C, D.n_ao[0], D.n_ao[1])*D.S_new[abcd]*D.p_SVD*2+
                        E_1el_calc(D.J, D.DM_A,D.n_ao[0], D.n_ao[1])*D.p_SVD*2+
                        E_1el_calc(D.J, D.DM_B,D.n_ao[0], D.n_ao[1])*D.p_SVD*2+
                        E_1el_calc(D.J, DM_A_1st_L,D.n_ao[0], D.n_ao[1])*D.p_SVD*2+
                        E_1el_calc(D.J, DM_A_1st_R,D.n_ao[0], D.n_ao[1])*D.p_SVD*2+
                        E_1el_calc(D.J, DM_B_1st_L,D.n_ao[0], D.n_ao[1])*D.p_SVD*2+
                        E_1el_calc(D.J, DM_B_1st_R,D.n_ao[0], D.n_ao[1])*D.p_SVD*2-
                        E_1el_calc(D.K, D.DM_C, D.n_ao[0], D.n_ao[1])*D.S_new[abcd]*D.p_SVD-
                        E_1el_calc(D.K, D.DM_A,D.n_ao[0], D.n_ao[1])*D.p_SVD-
                        E_1el_calc(D.K, D.DM_B,D.n_ao[0], D.n_ao[1])*D.p_SVD-
                        E_1el_calc(D.K, DM_A_1st_L,D.n_ao[0], D.n_ao[1])*D.p_SVD-
                        E_1el_calc(D.K, DM_A_1st_R,D.n_ao[0], D.n_ao[1])*D.p_SVD-
                        E_1el_calc(D.K, DM_B_1st_L,D.n_ao[0], D.n_ao[1])*D.p_SVD-
                        E_1el_calc(D.K, DM_B_1st_R,D.n_ao[0], D.n_ao[1])*D.p_SVD;
        double H_1el = E_1el_calc(A->H_AO, D.DM_T, D.n_ao[0], D.n_ao[1])*D.p_SVD;
        
        
        d_o_DMV_gen(DM_A_1st_L,
                    D.L_ACT, D.n_ao[0],
                    D.R_ACT, D.n_ao[1],
                    D.ACT_DMT_A_1st_order+abcd*D.n_act[0]*D.n_act[1], D.n_act[0],D.BUF);
        
        d_o_DMV_gen(DM_B_1st_L,
                    D.L_ACT, D.n_ao[0],
                    D.R_ACT, D.n_ao[1],
                    D.ACT_DMT_B_1st_order+abcd*D.n_act[0]*D.n_act[1], D.n_act[0],D.BUF);     
        
        D.H_2el[abcd]+=-E_1el_calc(K_1st, DM_A_1st_L,D.n_ao[0], D.n_ao[1])*D.p_SVD
                       -E_1el_calc(K_1st, DM_B_1st_L,D.n_ao[0], D.n_ao[1])*D.p_SVD;
        
        
        
        
        
        D.S_new[abcd] = D.S_new[abcd] * D.p_SVD;
        double d_x= E_1el_calc(A->Dx_AO, D.DM_T, D.n_ao[0], D.n_ao[1])*D.p_SVD-A->Dx_nuc*D.S_new[abcd];
        double d_y= E_1el_calc(A->Dy_AO, D.DM_T, D.n_ao[0], D.n_ao[1])*D.p_SVD-A->Dy_nuc*D.S_new[abcd];
        double d_z= E_1el_calc(A->Dz_AO, D.DM_T, D.n_ao[0], D.n_ao[1])*D.p_SVD-A->Dz_nuc*D.S_new[abcd];
        
        int i=a*A->n_states[1]+b;
        int j=c*B->n_states[1]+d;
        print_CI_results(out_file, i, j, D.S_new[abcd], A->V_nuc, H_1el, D.H_2el[abcd], d_x,d_y, d_z);
        
    }
    printf_timer("H and S calculation");
    
    
    
    D.clear();
    
    delete[] DM_A_1st_L;
    delete[] DM_B_1st_L;
    delete[] DM_A_1st_R;
    delete[] DM_B_1st_R;

    
    return 0;
}


