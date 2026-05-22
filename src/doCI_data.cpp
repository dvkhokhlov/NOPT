//standart
# include <stdio.h>

//user
# include "blas_link.h"
# include "molecule.h"
# include "doCI_data.h"
# include "doCI_matr.h"
# include "libint_link.h"
# include "matr.h"
# include "timer.h"
# include "etc.h"
# include "from_hash.h"
# include "defaults.h"

int doCI_data::first_alloc(molecule * ext_A, molecule * ext_B){
    
    two_charge_states=0;
    
    A = ext_A;
    B = ext_B;
    
    n_ao     = new int[2];
    n_alp_el = new int[2];
    n_bet_el = new int[2];
    n_act    = new int[2];
    VEC = new double *[2];
    
    
    n_alp_el_f = new  int*[2];
    n_bet_el_f = new  int*[2];
    n_act_f    = new  int*[2];
    n_cor      = new  int [2];
    n_alp_el_f[0] = A->n_act_el_alp;
    n_alp_el_f[1] = B->n_act_el_alp;
    n_bet_el_f[0] = A->n_act_el_bet;
    n_bet_el_f[1] = B->n_act_el_bet;
    n_act_f   [0] = A->n_act_orb   ;
    n_act_f   [1] = B->n_act_orb   ;
    n_cor     [0] = A->n_cor_orb   ;
    n_cor     [1] = B->n_cor_orb   ;
    
    
    n_alp_el[0] = 0; for(int i=0; i<A->n_frag;i++) n_alp_el[0]+=A->n_act_el_alp[i];
    n_alp_el[1] = 0; for(int i=0; i<B->n_frag;i++) n_alp_el[1]+=B->n_act_el_alp[i];
    
    n_bet_el[0] = 0; for(int i=0; i<A->n_frag;i++) n_bet_el[0]+=A->n_act_el_bet[i];
    n_bet_el[1] = 0; for(int i=0; i<B->n_frag;i++) n_bet_el[1]+=B->n_act_el_bet[i];
    
    n_act   [0] = 0; for(int i=0; i<A->n_frag;i++) n_act   [0]+=A->n_act_orb   [i];
    n_act   [1] = 0; for(int i=0; i<B->n_frag;i++) n_act   [1]+=B->n_act_orb   [i];    
    
//     for(int i=0; i<2; i++) n_ao    [i] = A->basis[i].n_ao;
    n_ao[0] = A->n_ao;
    n_ao[1] = B->n_ao;
    
    VEC[0]=A->MO_VEC;
    VEC[1]=B->MO_VEC;

//     fprintf(out_stream,"d1 = %d\nv1 = %d\na1 = %d\nb1 = %d\n\n",n_cor[0],n_act[0],n_alp_el[0],n_bet_el[0]);
//     fprintf(out_stream,"d2 = %d\nv2 = %d\na2 = %d\nb2 = %d\n\n",n_cor[1],n_act[1],n_alp_el[1],n_bet_el[1]);
//     getchar();
//     fprintf(out_stream,"\n\n");

    
    
    n_svd=n_cor[0];
    n_zero_svd=0;
    n_frag = A->n_frag;
    
//     if(n_cor[0]>n_cor[1]){
//         fprintf(out_stream,"WARNING: different n_cor\n");
//         A->add_cor_2_AS(0,0,n_cor[0]-n_cor[1]);
//     }
//     if(n_cor[0]<n_cor[1]){
//         fprintf(out_stream,"WARNING: different n_cor\n");
//         A->add_cor_2_AS(1,0,n_cor[1]-n_cor[0]);
//     }
//     
//     if(n_act[0]>n_act[1]){
//         fprintf(out_stream,"WARNING: different n_act\n");
//         A->add_vac_2_AS(1,0,n_act[0]-n_act[1]);
//     }
//     if(n_act[0]<n_act[1]){
//         fprintf(out_stream,"WARNING: different n_act\n");
//         A->add_vac_2_AS(0,0,n_act[1]-n_act[0]);
//     }
    
        
    //fisrt allocation
    J    = new double[n_ao[0]*n_ao[1]];
    K    = new double[n_ao[0]*n_ao[1]];
    BUF  = new double[n_ao[0]*n_ao[0] ];
    S_DD = new double[n_cor[0]*n_cor[1]];
    L_MO = new double[n_cor[0]*n_svd];
    R_MO = new double[n_cor[1]*n_svd];
    SVD  = new double[n_svd];
    L_COR= new double[(n_cor[0]/*+n_act[0]*/)*n_ao[0]];
    R_COR= new double[(n_cor[1]/*+n_act[1]*/)*n_ao[1]];
    
    DM_C = new double[n_ao[0]*n_ao[1]];
    DM_A= new double[n_ao[0]*n_ao[1]];
    DM_B= new double[n_ao[0]*n_ao[1]];
    DM_T = new double[n_ao[0]*n_ao[1]];
    DM_TA= new double[n_ao[0]*n_ao[1]];
    DM_TB= new double[n_ao[0]*n_ao[1]];
    
    ci1_f = new ci_map*[A->n_frag];
    ci2_f = new ci_map*[B->n_frag];
    for(int i_f=0; i_f<A->n_frag; i_f++)ci1_f[i_f] = new ci_map[A->n_states[i_f]];
    for(int i_f=0; i_f<B->n_frag; i_f++)ci2_f[i_f] = new ci_map[B->n_states[i_f]];
//     for(int i_f=0; i_f<B->n_frag; i_f++){fprintf(out_stream,"A_nstate[%d]= %d\n",i_f,A->n_states[i_f]);getchar();}
//     for(int i_f=0; i_f<B->n_frag; i_f++){fprintf(out_stream,"B_nstate[%d]= %d\n",i_f,B->n_states[i_f]);getchar();}
    
    
    n_st_total = new int[2];
    
    n_st_total[0]=1;
    n_st_total[1]=1;
    
    n_states = new int*[2];
    n_states[0] = A->n_states;
    n_states[1] = B->n_states;
    
   for(int i=0,n_vt=0; i<n_frag;i++){//delete n_vt
        n_st_total[0]=n_st_total[0]*n_states[0][i];
        n_st_total[1]=n_st_total[1]*n_states[1][i];
    }
    //write proper size !!!
    h_2e      = new double[                  n_st_total[0]*n_st_total[1]*10];
    ACT_DET_A = new double[                  n_st_total[0]*n_st_total[1]*10];
    ACT_ADJ_A = new double[n_act[0]*n_act[1]*n_st_total[0]*n_st_total[1]*10];
    ACT_ADJ_B = new double[n_act[0]*n_act[1]*n_st_total[0]*n_st_total[1]*10];
    
    set_zero_matr(ACT_DET_A,                  n_st_total[0]*n_st_total[1]*10);
    set_zero_matr(ACT_ADJ_A,n_act[0]*n_act[1]*n_st_total[0]*n_st_total[1]*10);
    set_zero_matr(ACT_ADJ_B,n_act[0]*n_act[1]*n_st_total[0]*n_st_total[1]*10);
    
    
    
    S_new=NULL;
    H_2el=NULL;
    
    ACT_DMT_A=NULL;
    ACT_DMT_B=NULL;
    
    return 0;
}

int doCI_data::gen_aldet_data(){
    
    fprintf(out_stream,"\nWARNING: doCI_data uses default value CAS_PRINT_NUMBER_DEFAULT\n\n");
    if(two_charge_states==0)
        aldet.get_dim(n_act[0], n_alp_el[0], n_bet_el[0], 2, 0, CAS_PRINT_NUMBER_DEFAULT);
    if(two_charge_states==1){
        aldet.get_dim(n_act[0], n_alp_el[0], n_bet_el[0], 1, 0, CAS_PRINT_NUMBER_DEFAULT);
        aldet_2.get_dim(n_act[1], n_alp_el[1], n_bet_el[1], 1, 0, CAS_PRINT_NUMBER_DEFAULT);
    }
    
    
    if(n_frag!=1) return 0;
    
    ///READING FROM MOLECULE IS TURNED OFF!!!!!!!!!
//     aldet.read_set_from_mol(A, 0, 0, n_st_total[0], 0);
    aldet.copy_coef(0, A->CI+0, A->CI[0].n_states[0], 0, 1);
//     printf_timer("read A");
    if(two_charge_states==0)
        aldet.copy_coef(1, B->CI+0, A->CI[0].n_states[0], 0, 1);
//         aldet.read_set_from_mol(B, 0, 0, n_st_total[1], 1);
    if(two_charge_states==1){
        aldet_2.copy_coef(0, B->CI+0, A->CI[0].n_states[0], 0, 1);///????
//         aldet_2.read_set_from_mol(B, 0, 0, n_st_total[1], 0);
    }
    
    
//     printf_timer("red B");
    
    
    return 0;
}


int doCI_data::calc_S_MO(molecule * M){///to be deleted 
    
    M->calc_S_AO();
    
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                         n_ao[0],n_ao[0],n_ao[0],1.0,
                         M->S_AO,n_ao[0],
                         VEC[1],n_ao[0],0.0,
                         BUF,n_ao[0]);
     
     cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                         n_ao[0],n_ao[0],n_ao[0],1.0,
                         //A->MO_VEC[0],A->basis[0].n_ao,    ///TRY TO GET DIFFERENT VEC!!!!!!!!!!!!!!
                         VEC[0],n_ao[0],  ///CHECK HERE (see 1 line up) !!!!!!!!!!
                         BUF,n_ao[0],0.0,
                         M->S_MO,n_ao[0]);
     
     return 0;
     
}

int doCI_data::AO_to_MO(double * M){
    
    
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                         n_ao[0],n_ao[0],n_ao[0],1.0,
                         M,n_ao[0],
                         VEC[1],n_ao[0],0.0,
                         BUF,n_ao[0]);
     
     cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                         n_ao[0],n_ao[0],n_ao[0],1.0,
                         //A->MO_VEC[0],A->basis[0].n_ao,    ///TRY TO GET DIFFERENT VEC!!!!!!!!!!!!!!
                         VEC[0],n_ao[0],  ///CHECK HERE (see 1 line up) !!!!!!!!!!
                         BUF,n_ao[0],0.0,
                         M,n_ao[0]);
//      for(int i=0;i<216;i++){
//          for(int j=0;j<216;j++)
//              fprintf(out_stream," % .3e  ",M[i*n_ao[0]+j]);
//          fprintf(out_stream,"\n");
//      }
//      fprintf(out_stream,"stopped at SMO\n");
//      exit(1);
     
     
     return 0;
     
}

int doCI_data::cor_svd(molecule * M){
    
    gen_S_DD(S_DD,M->S_MO,n_ao[0],n_cor[0],n_ao[1],n_cor[1]);
//     fprintf(out_stream,"S_DD:\n");
//     PrintMatr(S_DD,n_cor[0],n_cor[0],1);
    
    if(n_cor[0]!=0){
        lapack_svd(S_DD,R_MO,L_MO,SVD,n_cor[0]);////i!=j - ????????????????????
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                    n_cor[0],n_ao[1],n_cor[0],1.0,
                    R_MO,n_cor[0],    
                    VEC[1],n_ao[1],0.0,
                    R_COR,n_ao[1]);
        cblas_dgemm(CblasRowMajor,CblasTrans ,CblasNoTrans,
                    n_cor[0],n_ao[1],n_cor[0],1.0,
                    L_MO,n_cor[0],    
                    VEC[0],n_ao[1],0.0,
                    L_COR,n_ao[1]);
    }
    for(int i_svd=0;i_svd<n_cor[0];i_svd++)
        if(SVD[i_svd]<svd_eps){
            fprintf(out_stream,"warning: SVD[%d]= %e < cut-off value(%e), check your input and further calculations\n",i_svd,SVD[i_svd],svd_eps);
            n_zero_svd++;
        }
    if(n_zero_svd){
        fprintf(out_stream,"n_zero_svd !=0. block not written\n");
        exit(1);
//         A->add_cor_2_AS(0,0,n_zero_svd);
//         A->add_cor_2_AS(1,0,n_zero_svd);
    }
    
    p_SVD=1.0;
    for(int i_d=0;i_d<n_cor[0];i_d++)p_SVD=p_SVD*SVD[i_d]*SVD[i_d];
    d_o_DM_gen(DM_C, 1, L_COR, n_ao[0], R_COR, n_ao[1], SVD, n_cor[0]);//!!!!!!
    
     E_core = E_1el_calc(M->H_AO , DM_C, n_ao[0], n_ao[1])*p_SVD;
    Dx_core = E_1el_calc(M->Dx_AO, DM_C, n_ao[0], n_ao[1])*p_SVD;
    Dy_core = E_1el_calc(M->Dy_AO, DM_C, n_ao[0], n_ao[1])*p_SVD;
    Dz_core = E_1el_calc(M->Dz_AO, DM_C, n_ao[0], n_ao[1])*p_SVD;
    
    
    return 0;
}

int doCI_data::cor_svd_PT(molecule * M){
    
//     fprintf(out_stream,"S_DVDV:\n");
//     for(int i=0;i<n_cor[0]+n_act[0];i++){
//         for(int j=0;j<n_cor[0]+n_act[0];j++)
//             fprintf(out_stream,"%.4e ",M->S_MO[i*n_ao[0]+j]);
//         fprintf(out_stream,"\n");
//     }
    
    
    gen_S_DD(S_DD,M->S_MO,n_ao[0],n_cor[0],n_ao[1],n_cor[1]);
//     fprintf(out_stream,"S_DD:\n");
//     PrintMatr(S_DD,n_cor[0],n_cor[0],1);
    
    int n_cor_A = A->n_cor_orb_f[0];
    int n_cor_B = A->n_cor_orb_f[1];
    
    double * S_AA = new double[n_cor_A*n_cor_A];
    double * S_BB = new double[n_cor_B*n_cor_B];
    
    for(int i=0;i<n_cor_A;i++)
    for(int j=0;j<n_cor_A;j++)
        S_AA[i*n_cor_A+j]=S_DD[i*n_cor[0]+j];
    
    for(int i=0;i<n_cor_B;i++)
    for(int j=0;j<n_cor_B;j++)
        S_BB[i*n_cor_B+j]=S_DD[(i+n_cor_A)*n_cor[0]+j+n_cor_A];
//     
//     fprintf(out_stream,"S_AA:\n");
//     PrintMatr(S_AA,n_cor_A,n_cor_A,1);
    
//     fprintf(out_stream,"S_BB:\n");
//     PrintMatr(S_BB,n_cor_B,n_cor_B,1);
    
    double * L_MO_A= new double[n_cor_A*n_cor_A];
    double * R_MO_A= new double[n_cor_A*n_cor_A];
    
    double * L_MO_B= new double[n_cor_B*n_cor_B];
    double * R_MO_B= new double[n_cor_B*n_cor_B];
    
    lapack_svd(S_AA,R_MO_A,L_MO_A,SVD,n_cor_A);////i!=j - ????????????????????
    
//     fprintf(out_stream,"L_MO_A:\n");
//     PrintMatr(L_MO_A,n_cor_A,n_cor_A,1);
    
//     fprintf(out_stream,"R_MO_A:\n");
//     PrintMatr(R_MO_A,n_cor_A,n_cor_A,1);
    
//     fprintf(out_stream,"SVD_A:\n");
//     PrintMatr(SVD,n_cor_A,1,0);
    
    lapack_svd(S_BB,R_MO_B,L_MO_B,SVD+n_cor_A,n_cor_B);////i!=j - ????????????????????
    
//     fprintf(out_stream,"L_MO_B:\n");
//     PrintMatr(L_MO_B,n_cor_B,n_cor_B,1);
    
//     fprintf(out_stream,"R_MO_B:\n");
//     PrintMatr(R_MO_B,n_cor_B,n_cor_B,1);
    
//     fprintf(out_stream,"SVD_B:\n");
//     PrintMatr(SVD,n_cor[0],1,0);
    
    set_zero_matr(R_COR,n_cor[0]*n_ao[1]);
    set_zero_matr(L_COR,n_cor[0]*n_ao[1]);
//     fprintf(out_stream,"R_COR:\n");
//     PrintMatr(R_COR,n_cor_A,n_ao[1],1);
    
    
    
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
            n_cor_A,n_ao[1],n_cor_A,1.0,
            R_MO_A,n_cor_A,    
            VEC[1],n_ao[1],0.0,
            R_COR,n_ao[1]);
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
            n_cor_B,n_ao[1],n_cor_B,1.0,
            R_MO_B,n_cor_B,    
            VEC[1]+n_ao[1]*n_cor_A,n_ao[1],0.0,
            R_COR+n_ao[1]*n_cor_A,n_ao[1]);
    
    
    cblas_dgemm(CblasRowMajor,CblasTrans ,CblasNoTrans,
            n_cor_A,n_ao[1],n_cor_A,1.0,
            L_MO_A,n_cor_A,    
            VEC[0],n_ao[1],0.0,
            L_COR,n_ao[1]);
    
    cblas_dgemm(CblasRowMajor,CblasTrans ,CblasNoTrans,
            n_cor_B,n_ao[1],n_cor_B,1.0,
            L_MO_B,n_cor_B,    
            VEC[0]+n_ao[1]*n_cor_A,n_ao[1],0.0,
            L_COR+n_ao[1]*n_cor_A,n_ao[1]);
    
    
//     fprintf(out_stream,"L_COR:\n");
//     PrintMatr(L_COR,n_cor[0],n_ao[1],1);
    
    
//     exit(0);
    
    
    for(int i_svd=0;i_svd<n_cor[0];i_svd++)
        if(SVD[i_svd]<svd_eps){
            fprintf(out_stream,"warning: SVD[%d]= %e < cut-off value(%e), check your input and further calculations\n",i_svd,SVD[i_svd],svd_eps);
            n_zero_svd++;
        }
    if(n_zero_svd){
        fprintf(out_stream,"n_zero_svd !=0. block not written\n");
        exit(1);
//         A->add_cor_2_AS(0,0,n_zero_svd);
//         A->add_cor_2_AS(1,0,n_zero_svd);
    }
    
    p_SVD=SVD[0]*SVD[0];
        for(int i_d=1;i_d<n_cor[0];i_d++)p_SVD=p_SVD*SVD[i_d]*SVD[i_d];
    
    DM_C_F = new double *[2];
    
    DM_C_F[0] = new double [n_ao[0]*n_ao[1]];
    DM_C_F[1] = new double [n_ao[0]*n_ao[1]];
    
    d_o_DM_gen(DM_C_F[0], 1, L_COR                , n_ao[0], R_COR                , n_ao[1], SVD        , n_cor_A);//!!!!!!
    d_o_DM_gen(DM_C_F[1], 1, L_COR+n_cor_A*n_ao[0], n_ao[0], R_COR+n_cor_A*n_ao[1], n_ao[1], SVD+n_cor_A, n_cor_B);//!!!!!!
    
    for(int i=0; i<n_ao[0]*n_ao[1];i++)
        DM_C[i]=DM_C_F[0][i]+DM_C_F[1][i];
    
    delete[] L_MO_A;
    delete[] R_MO_A;
    delete[] L_MO_B;
    delete[] R_MO_B;
    
    
    delete[] S_AA;
    delete[] S_BB;
    
    M->calc_S_AO();
    
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                        n_cor[0],n_ao[0],n_ao[0],1.0,
                        L_COR,n_ao[0],//no change in lda!!
                        M->S_MO,n_ao[0],0.0,
                        BUF,n_ao[0]);
    
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        n_cor[0],n_cor[0],n_ao[0],1.0,
                        BUF,n_ao[0],
                        R_COR,n_ao[0],0.0,//no change in lda!!
                        S_DD,n_cor[0]);
//     fprintf(out_stream,"S_DD_final:\n");
//     PrintMatr(S_DD,n_cor[0],n_cor[0],1);
    
    L_COR_1st_order = new double[n_cor[0]*n_ao[0]];
    R_COR_1st_order = new double[n_cor[0]*n_ao[0]];
    
    set_zero_matr(L_COR_1st_order, n_cor[0]*n_ao[0]);
    set_zero_matr(R_COR_1st_order, n_cor[0]*n_ao[0]);
    
    for(int i=0; i<n_cor_B;i++)
    for(int k=0; k<n_ao[0];k++)
    for(int j=0; j<n_cor_A;j++)
        L_COR_1st_order[(i+n_cor_B)*n_ao[0]+k]-=L_COR[j*n_ao[0]+k]*S_DD[(i+n_cor_B)*(n_cor_A+n_cor_B)+j]/S_DD[j*(n_cor_A+n_cor_B)+j];
    
    for(int i=0; i<n_cor_B;i++)
    for(int k=0; k<n_ao[0];k++)
    for(int j=0; j<n_cor_A;j++)
        R_COR_1st_order[(i+n_cor_B)*n_ao[0]+k]-=R_COR[j*n_ao[0]+k]*S_DD[j*(n_cor_A+n_cor_B)+i+n_cor_B]/S_DD[j*(n_cor_A+n_cor_B)+j];
    
    
    DM_R_1st = new double[n_ao[0]*n_ao[1]];
    DM_L_1st = new double[n_ao[0]*n_ao[1]];
    
    d_o_DM_gen(DM_R_1st, 1, L_COR, n_ao[0], R_COR_1st_order, n_ao[1], SVD, n_cor[0]);//!!!!!!
    d_o_DM_gen(DM_L_1st, 1, L_COR_1st_order, n_ao[0], R_COR, n_ao[1], SVD, n_cor[0]);//!!!!!!
    
    
//     cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
//                         n_cor[0],n_ao[0],n_ao[0],1.0,
//                         L_COR_1st_order,n_ao[0],//no change in lda!!
//                         M->S_MO,n_ao[0],0.0,
//                         BUF,n_ao[0]);
//     
//     cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
//                         n_cor[0],n_cor[0],n_ao[0],1.0,
//                         BUF,n_ao[0],
//                         R_COR,n_ao[0],0.0,//no change in lda!!
//                         S_DD,n_cor[0]);
//     fprintf(out_stream,"S_L(1)R(0):\n");
//     PrintMatr(S_DD,n_cor[0],n_cor[0],1);
//     
//     cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
//                         n_cor[0],n_ao[0],n_ao[0],1.0,
//                         L_COR,n_ao[0],//no change in lda!!
//                         M->S_MO,n_ao[0],0.0,
//                         BUF,n_ao[0]);
//     
//     cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
//                         n_cor[0],n_cor[0],n_ao[0],1.0,
//                         BUF,n_ao[0],
//                         R_COR_1st_order,n_ao[0],0.0,//no change in lda!!
//                         S_DD,n_cor[0]);
//     fprintf(out_stream,"S_L(0)R(1):\n");
//     PrintMatr(S_DD,n_cor[0],n_cor[0],1);



    
    
    return 0;
}



// int doCI_data::calc_DM_C_MO(double * C){
//     
//     
//     d_o_DM_gen(C, 1, L_MO, n_cor[0], R_MO, n_cor[1], SVD, n_cor[0]);
//     fprintf(out_stream,"L_MO:\n");
//     PrintMatr(L_MO,n_cor[0],n_cor[1],0);
//     
//     
//     fprintf(out_stream,"DM_C_MO:\n");
//     PrintMatr(C,n_cor[0],n_cor[1],0);
//     exit(1);
//     
//     return 0;
//     
// }
// 

int doCI_data::DM_gen_ortogonal(int s, double coef){
    
    set_zero_matr(DM_C, n_ao[s]*n_ao[s]);        
    
    for(int i=0;i<n_ao[s];i++)
    for(int j=0;j<n_ao[s];j++){
    for(int k=0;k<n_cor[s] ;k++){
        DM_C[i*n_ao[s]+j]+=coef*VEC[s][k*n_ao[s]+i]*VEC[s][k*n_ao[s]+j];
    }//fprintf(out_stream,"+ %e * %e / %e\n",L[(n_cor[s]-1)*dim_l+i],R[(n_cor[s]-1)*dim_r+j],sigma[n_cor[s]-1]);
//     fprintf(out_stream,"DM(%d,%d)=%e\n",i,j,DM[i*dim_r+j]);
//     PrintMatr(DM,n,n,1);
    }
//     PrintMatr(DM_C,n_ao[s],n_ao[s],1);
//     getchar();
    return 0;
}


int doCI_data::second_alloc(){
    
    L_ACT= new double[n_act[0]*n_ao[0]];
    R_ACT= new double[n_act[1]*n_ao[1]];
    U_ACT= new double[n_act[0]*n_act[0]];
    V_ACT= new double[n_act[1]*n_act[1]];
    S_SV = new double[n_act[0]*n_act[0]];
    S_FV = new double[n_act[0]*n_act[0]];
    S_SA = new double[n_alp_el[0]*n_alp_el[0]];
    S_SB = new double[n_bet_el[0]*n_bet_el[0]];
    
     H_ACT = new double[n_act[0]*n_act[0]];
    Dx_ACT = new double[n_act[0]*n_act[0]];
    Dy_ACT = new double[n_act[0]*n_act[0]];
    Dz_ACT = new double[n_act[0]*n_act[0]];
    
    
    
    
    L_ACT_F = new double *[n_frag];
    R_ACT_F = new double *[n_frag];

    act_INTS = new double[ n_act[0]*n_act[1]*n_act[0]*n_act[1]];;
    
    return 0;
    
}

int doCI_data::act_ort(molecule * M){
    
    act_to_cor_ort(L_ACT,M->S_MO+n_cor[0]*n_ao[0],n_svd,n_cor[0],n_act[0],n_ao[0],R_MO, VEC[0]+n_cor[0]*n_ao[0],SVD,L_COR,1);
    act_to_cor_ort(R_ACT,M->S_MO+n_cor[1]        ,n_svd,n_cor[1],n_act[1],n_ao[1],L_MO, VEC[1]+n_cor[1]*n_ao[1],SVD,R_COR,0);
    
//     if(n_act[0])
//     cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
//                         n_act[0],n_ao[0],n_ao[0],1.0,
//                         L_ACT,n_act[0],//no change in lda!!
//                         M->H_AO,n_ao[0],0.0,
//                         BUF,n_ao[0]);
//     if(n_act[0])
//     cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
//                         n_act[0],n_act[0],n_ao[0],1.0,
//                         BUF,n_ao[0],
//                         R_ACT,n_act[0],0.0,//no change in lda!!
//                         H_ACT,n_act[0]);
//     
//     if(n_act[0])
//     cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
//                         n_act[0],n_ao[0],n_ao[0],1.0,
//                         L_ACT,n_act[0],//no change in lda!!
//                         M->Dx_AO,n_ao[0],0.0,
//                         BUF,n_ao[0]);
//     if(n_act[0])
//     cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
//                         n_act[0],n_act[0],n_ao[0],1.0,
//                         BUF,n_ao[0],
//                         R_ACT,n_act[0],0.0,//no change in lda!!
//                         Dx_ACT,n_act[0]);
//     
//     if(n_act[0])
//     cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
//                         n_act[0],n_ao[0],n_ao[0],1.0,
//                         L_ACT,n_act[0],//no change in lda!!
//                         M->Dy_AO,n_ao[0],0.0,
//                         BUF,n_ao[0]);
//     if(n_act[0])
//     cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
//                         n_act[0],n_act[0],n_ao[0],1.0,
//                         BUF,n_ao[0],
//                         R_ACT,n_act[0],0.0,//no change in lda!!
//                         Dy_ACT,n_act[0]);
//     
//     if(n_act[0])
//     cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
//                         n_act[0],n_ao[0],n_ao[0],1.0,
//                         L_ACT,n_act[0],//no change in lda!!
//                         M->Dz_AO,n_ao[0],0.0,
//                         BUF,n_ao[0]);
//     if(n_act[0])
//     cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
//                         n_act[0],n_act[0],n_ao[0],1.0,
//                         BUF,n_ao[0],
//                         R_ACT,n_act[0],0.0,//no change in lda!!
//                         Dz_ACT,n_act[0]);
//     
    
    
    return 0;
}


int doCI_data::AO_to_act(double * M_act, double *M_AO){
    
    if(n_act[0])
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                        n_act[0],n_ao[0],n_ao[0],1.0,
                        L_ACT,n_act[0],//no change in lda!!
                        M_AO,n_ao[0],0.0,
                        BUF,n_ao[0]);
    if(n_act[0])
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                        n_act[0],n_act[0],n_ao[0],1.0,
                        BUF,n_ao[0],
                        R_ACT,n_act[0],0.0,//no change in lda!!
                        M_act,n_act[0]);
    
    
    return 0;
}




int doCI_data::cpy_VEC_to_ACT(){
    //copy with transposition
    for(int i=0; i<n_ao [0];i++)
    for(int j=0; j<n_act[0];j++)
        L_ACT[i*n_act[0]+j]=VEC[0][(j+n_cor[0])*n_ao[0]+i];
    
    for(int i=0; i<n_ao[0]*n_act[0];i++)
        R_ACT[i]=L_ACT[i];
    
    return 0;
}

int doCI_data::cpy_ACT_VEC_to_ACT(double * ACT_VEC){
    //copy with transposition
    for(int i=0; i<n_ao[0]*n_act[0];i++)
        L_ACT[i]=ACT_VEC[i];
    for(int i=0; i<n_ao[0]*n_act[0];i++)
        R_ACT[i]=ACT_VEC[i];
    
    return 0;
}


int doCI_data::calc_S_SV(molecule * M){
    
    M->calc_S_AO();
    if(n_act[0])
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                        n_act[0],n_ao[0],n_ao[0],1.0,
                        L_ACT,n_act[0],//no change in lda!!
                        M->S_AO,n_ao[0],0.0,
                        BUF,n_ao[0]);
    if(n_act[0])
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                        n_act[0],n_act[0],n_ao[0],1.0,
                        BUF,n_ao[0],
                        R_ACT,n_act[0],0.0,//no change in lda!!
                        S_SV,n_act[0]);
    
    // fprintf(out_stream,"S_SV_full:\n");
    // PrintMatr(S_SV,n_act[0],n_act[0],1);
    
    
    return 0;
}

int doCI_data::ci_ortogonalization(double * S_MO){
    

    
    for(int i=0,n_vt=0; i<n_frag;i++){
        // calculation of activent orbital overlap
        L_ACT_F[i] = L_ACT + n_vt /** n_act[0]*/;
        R_ACT_F[i] = R_ACT + n_vt /** n_act[0]*/;
        int ld = n_vt+n_act_f[0][i];
        double * L_O_V = new double[n_act_f[0][i]*ld];
        double * R_O_V = new double[n_act_f[0][i]*ld];
        double * U_malm = new double[ld*ld];
        double * V_malm = new double[ld*ld];
        set_zero_matr(V_malm,ld*ld);for(int ii=0;ii<ld;ii++)V_malm[ii*(ld+1)]=1.0;
        set_zero_matr(U_malm,ld*ld);for(int ii=0;ii<ld;ii++)U_malm[ii*(ld+1)]=1.0;
        
        if(act_to_act_ort(L_O_V,S_SV,n_vt,n_act_f[0][i],n_act[0],1)){
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                                n_ao[0],n_act_f[0][i],n_act_f[0][i]+n_vt,1.0,
                                L_ACT,n_act[0],
                                L_O_V,n_act_f[0][i]+n_vt,0.0,
                                BUF,n_act_f[0][i]);
            for(int j=0;j<n_ao[0];j++)memcpy(L_ACT_F[i]+j*n_act[0],BUF+j*n_act_f[0][i],n_act_f[0][i]*sizeof(double));
            for(int jj=   0;jj<ld;jj++)
            for(int ii=n_vt;ii<ld;ii++)
                U_malm[jj*ld+ii]=L_O_V[(ii-n_vt)*ld+jj];
        }
        if(act_to_act_ort(R_O_V,S_SV,n_vt,n_act_f[0][i],n_act[0],0)){
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                                n_ao[0],n_act_f[0][i],n_act_f[0][i]+n_vt,1.0,
                                R_ACT,n_act[0],
                                R_O_V,n_act_f[0][i]+n_vt,0.0,
                                BUF,n_act_f[0][i]);
            for(int j=0;j<n_ao[0];j++)memcpy(R_ACT_F[i]+j*n_act[0],BUF+j*n_act_f[0][i],n_act_f[0][i]*sizeof(double));
            for(int jj=   0;jj<ld;jj++)
            for(int ii=n_vt;ii<ld;ii++)
                V_malm[jj*ld+ii]=R_O_V[(ii-n_vt)*ld+jj];
        }
        
        if(n_act_f[0][i])
        cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                            n_act_f[0][i],n_ao[0],n_ao[0],1.0,
                            L_ACT_F[i],n_act[0],//no change in lda!!
                            S_MO,n_ao[0],0.0,
                            BUF,n_ao[0]);
        if(n_act_f[0][i])
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                            n_act_f[0][i],n_act_f[0][i],n_ao[0],1.0,
                            BUF,n_ao[0],
                            R_ACT_F[i],n_act[0],0.0,//no change in lda!!
                            S_FV,n_act_f[0][i]);
        
        
        lapack_svd(S_FV,V_ACT,U_ACT,SVD,n_act_f[0][i]);
        
        //transformation of orthogonalized active orbitals to singular active orbitals
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                            n_ao[0],n_act_f[0][i],n_act_f[0][i],1.0,
                            L_ACT_F[i],n_act[0],
                            U_ACT,n_act_f[0][i],0.0,
                            BUF,n_act_f[0][i]);
        for(int j=0;j<n_ao[0];j++)memcpy(L_ACT_F[i]+j*n_act[0],BUF+j*n_act_f[0][i],n_act_f[0][i]*sizeof(double));
        
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            n_ao[0],n_act_f[0][i],n_act_f[0][i],1.0,
                            R_ACT_F[i],n_act[0],
                            V_ACT,n_act_f[0][i],0.0,
                            BUF,n_act_f[0][i]);
        
        
        for(int j=0;j<n_ao [0]     ;j++)
        for(int k=0;k<n_act_f[0][i];k++){
            R_ACT_F[i][j*n_act[0]+k]=BUF[j*n_act_f[0][i]+k]/SVD[k];
        }
        
//         V_ACT -> VT
        for(int k=0;k<n_act_f[0][i];k++)
        for(int j=0;j<n_act_f[0][i];j++){
            V_malm[(k+n_vt)*ld+j+n_vt]=V_ACT[j*n_act_f[0][i]+k]/SVD[j];
        }
        for(int j=0;j<n_act_f[0][i];j++)
        for(int k=0;k<n_act_f[0][i];k++){
            U_malm[(j+n_vt)*ld+k+n_vt]=U_ACT[j*n_act_f[0][i]+k];
        }
        
        
        for(int j=0;j<n_act_f[0][i];j++)
        for(int k=j;k<n_act_f[0][i];k++){
            BUF[0]=V_ACT[j*n_act_f[0][i]+k]*SVD[j];
            V_ACT[j*n_act_f[0][i]+k]=V_ACT[k*n_act_f[0][i]+j]*SVD[k];
            V_ACT[k*n_act_f[0][i]+j]=BUF[0];
        }
        
        for(int i_S=0; i_S<n_states[0][i]; i_S++) {
            double * ci_new_1;
            double * ci_new_2;
//             fprintf(out_stream,"i = %d\n",i_S);
            init_ci_map(&(ci1_f[i][i_S]));
            init_ci_map(&(ci2_f[i][i_S]));
            
            
            ci_from_ci(n_act_f[0][i]+n_vt,n_alp_el_f[0][i], n_bet_el_f[0][i], A->CI+i, 0, i, i_S, ci_new_1);
            ci_from_ci(n_act_f[0][i]+n_vt,n_alp_el_f[0][i], n_bet_el_f[0][i], B->CI+i, 0, i, i_S, ci_new_2);
            
            malmqvist   (n_act_f[0][i]+n_vt,n_alp_el_f[0][i], n_bet_el_f[0][i],ci_new_1, ci_new_2, U_malm, V_malm);
            
            
            hash_from_ci(n_act_f[0][i]+n_vt, 0, n_alp_el_f[0][i], n_bet_el_f[0][i], ci1_f[i][i_S], i_S, ci_new_1);
            hash_from_ci(n_act_f[0][i]+n_vt, 0, n_alp_el_f[0][i], n_bet_el_f[0][i], ci2_f[i][i_S], i_S, ci_new_2);

            delete[] ci_new_1;
            delete[] ci_new_2;
            
        }

        
        
        
        if(n_act_f[0][i])
        cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                            n_act[0],n_ao[0],n_ao[0],1.0,
                            L_ACT,n_act[0],
                            S_MO,n_ao[0],0.0,
                            BUF,n_ao[0]);
        if(n_act[0])
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                            n_act[0],n_act[0],n_ao[0],1.0,
                            BUF,n_ao[0],
                            R_ACT,n_act[0],0.0,
                            S_SV,n_act[0]);
        
        
        n_vt+=n_act_f[0][i];
        
        delete [] L_O_V;
        delete [] R_O_V;
        delete[] U_malm;
        delete[] V_malm;
    }
    
    
    return 0;
}

int doCI_data::read_act_ints(char * prefix){
    
    printf("read_act_ints is possibly unneeded\n");
    exit(0);
    
//     char * full_name;
//     get_full_name(&full_name, prefix, "act_ints.dat");
//     
//     FILE * vo;
//     fprintf(stderr,"reading %s\n", full_name);
//     vo = fopen(full_name,"r");
//     if(vo==NULL){
// 		fprintf(out_stream,"ERROR: could not open input file %s\n\n", full_name);
// 	   	exit(EXIT_FAILURE);
//     };
//     set_zero_matr(act_INTS,n_act[0]*n_act[0]*n_act[1]*n_act[1]);
//     XDR out_xdr2;
//     xdrstdio_create (&out_xdr2, vo, XDR_DECODE);
//     
//     for(int ii=0; ii<n_act[0]*n_act[0]*n_act[1]*n_act[1];ii++) xdr_double(&out_xdr2,&act_INTS[ii]);
//     for(int ii=0; ii<n_ao[0]*n_ao[1];ii++) xdr_double(&out_xdr2,&J[ii]);
//     for(int ii=0; ii<n_ao[0]*n_ao[1];ii++) xdr_double(&out_xdr2,&K[ii]);
//     
//     xdr_destroy(&out_xdr2);
//     fclose(vo);
    
    return 0;
}

int doCI_data::write_act_ints(char * prefix){
    
    printf("write_act_ints is possibly unneeded\n");
    exit(0);
    
    
// //     fprintf(stderr,"write act_ints?\n");
// //     getchar();
//     
//     char * full_name;
//     get_full_name(&full_name, prefix, "act_ints.dat");
//     
//     FILE * vo;
//     vo = fopen(full_name,"w");
//     if(vo==NULL){
// 		fprintf(out_stream,"ERROR: could not open output file %s\n\n", full_name);
// 	   	exit(EXIT_FAILURE);
//     };
//     
//     XDR out_xdr2;
//     xdrstdio_create (&out_xdr2, vo, XDR_ENCODE);
//     
//     for(int ii=0; ii<n_act[0]*n_act[0]*n_act[1]*n_act[1];ii++) xdr_double(&out_xdr2,&act_INTS[ii]);
//     for(int ii=0; ii<n_ao[0]*n_ao[1];ii++) xdr_double(&out_xdr2,&J[ii]);
//     for(int ii=0; ii<n_ao[0]*n_ao[1];ii++) xdr_double(&out_xdr2,&K[ii]);
//     
//     xdr_destroy(&out_xdr2);
//     fclose(vo);
    
    return 0;
}

int doCI_data::ci_link(){
    
    init_ci_map_arr(&ci1);
    init_ci_map_arr(&ci2);
    
    multi_iterator i_S;
    multi_iterator j_S;
    i_S.set(n_states[0],n_frag);
    j_S.set(n_states[1],n_frag);
    
    int i,j;
    
    for(i_S.zero(),i=0; i_S.not_ended();i_S.next(),i++){ 
        // fprintf(stderr,"linking for %d bra state\n",i);
        ci_array_link_to_array(&ci1,ci1_f,n_frag,i_S.number,n_act_f[0],i,n_st_total[0]);
        
    }
    
    aldet.read_set_from_ci_map_arr(&ci1, n_st_total[0], 0);
    
    for(j_S.zero(),j=0; j_S.not_ended();j_S.next(),j++){ 
        // fprintf(stderr,"linking for %d ket state\n",j);
        ci_array_link_to_array(&ci2,ci2_f,n_frag,j_S.number,n_act_f[0],j,n_st_total[1]);
        
    }
    if(two_charge_states==0)
        aldet.read_set_from_ci_map_arr(&ci2, n_st_total[1], 1);
    if(two_charge_states==1)
        aldet_2.read_set_from_ci_map_arr(&ci2, n_st_total[1], 0);
    
    return 0;
}

int doCI_data::link_and_transform_for_separate_spaces(molecule * M){
    
    M->calc_S_AO();
    
    if(n_act[0])
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                        n_act[0],n_ao[0],n_ao[0],1.0,
                        L_ACT,n_act[0],//no change in lda!!
                        M->S_MO,n_ao[0],0.0,
                        BUF,n_ao[0]);
    if(n_act[0])
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                        n_act[0],n_act[0],n_ao[0],1.0,
                        BUF,n_ao[0],
                        R_ACT,n_act[0],0.0,//no change in lda!!
                        S_SV,n_act[0]);
    

    if(n_act[0])
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                        n_act[0],n_ao[0],n_ao[0],1.0,
                        L_ACT,n_act[0],//no change in lda!!
                        M->S_MO,n_ao[0],0.0,
                        BUF,n_ao[0]);
    if(n_act[0])
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                        n_act[0],n_act[0],n_ao[0],1.0,
                        BUF,n_ao[0],
                        R_ACT,n_act[0],0.0,//no change in lda!!
                        S_FV,n_act[0]);
    if(A->n_frag>2){
        fprintf(stdout,"doCI_data::link_and_transform_for_separate_space is unable to use more then 2 fragments: n_frag = %d\n",A->n_frag);
        exit(1);
    }
    
    if(A->n_frag<2){
        fprintf(stdout,"doCI_data::link_and_transform_for_separate_space is unable to use less then 2 fragments: n_frag = %d\n",A->n_frag);
        exit(1);
    }
    
    double * L_O = new double[n_act_f[0][1]*(n_act_f[0][0]+n_act_f[0][1])];
    double * R_O = new double[n_act_f[0][1]*(n_act_f[0][0]+n_act_f[0][1])];
            
    
    aldet_data * aldet_bra = new aldet_data[A->n_frag];
    aldet_data * aldet_ket = new aldet_data[A->n_frag];
    
    for(int i=0,n_vt=0; i<A->n_frag;i++){
        // calculation of active orbital overlap
//         if(i)fprintf(out_stream,"n_vt = %d\n",n_vt);
        L_ACT_F[i] = L_ACT + n_vt /** n_act[0]*/;
        R_ACT_F[i] = R_ACT + n_vt /** n_act[0]*/;
        double * L_O_V = new double[n_act_f[0][i]*(n_vt+n_act_f[0][i])];
        double * R_O_V = new double[n_act_f[0][i]*(n_vt+n_act_f[0][i])];
        
        
//         double * L_FV = new double[n_act_f[0][i]*n_act_f[0][i]];
//         double * R_FV = new double[n_act_f[0][i]*n_act_f[0][i]];
        
        if(act_to_act_ort(L_O_V,S_SV,n_vt,n_act_f[0][i],n_act[0],1)){
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                                n_ao[0],n_act_f[0][i],n_act_f[0][i]+n_vt,1.0,
                                L_ACT,n_act[0],
                                L_O_V,n_act_f[0][i]+n_vt,0.0,
                                BUF,n_act_f[0][i]);
            for(int j=0;j<n_ao[0];j++)memcpy(L_ACT_F[i]+j*n_act[0],BUF+j*n_act_f[0][i],n_act_f[0][i]*sizeof(double));
            
        }
        if(act_to_act_ort(R_O_V,S_SV,n_vt,n_act_f[0][i],n_act[0],0)){
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                                n_ao[0],n_act_f[0][i],n_act_f[0][i]+n_vt,1.0,
                                R_ACT,n_act[0],
                                R_O_V,n_act_f[0][i]+n_vt,0.0,
                                BUF,n_act_f[0][i]);
            for(int j=0;j<n_ao[0];j++)memcpy(R_ACT_F[i]+j*n_act[0],BUF+j*n_act_f[0][i],n_act_f[0][i]*sizeof(double));
            
        }
        
        if(n_act_f[0][i])
        cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                            n_act_f[0][i],n_ao[0],n_ao[0],1.0,
                            L_ACT_F[i],n_act[0],//no change in lda!!
                            M->S_MO,n_ao[0],0.0,
                            BUF,n_ao[0]);
        if(n_act_f[0][i])
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                            n_act_f[0][i],n_act_f[0][i],n_ao[0],1.0,
                            BUF,n_ao[0],
                            R_ACT_F[i],n_act[0],0.0,//no change in lda!!
                            S_FV,n_act_f[0][i]);
        
//         fprintf(out_stream,"S_FV(%d)\n",i);
//         PrintMatr(S_FV,n_act_f[0][i],n_act_f[0][i],1);
        lapack_svd(S_FV,V_ACT,U_ACT,SVD,n_act_f[0][i]);

                
        //transformation of orthogonalized active orbitals to singular active orbitals
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                            n_ao[0],n_act_f[0][i],n_act_f[0][i],1.0,
                            L_ACT_F[i],n_act[0],
                            U_ACT,n_act_f[0][i],0.0,
                            BUF,n_act_f[0][i]);
        for(int j=0;j<n_ao [0]     ;j++)
        for(int k=0;k<n_act_f[0][i];k++){
            L_ACT_F[i][j*n_act[0]+k]=BUF[j*n_act_f[0][i]+k];///sqrt(SVD[k]);
        }
        
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            n_ao[0],n_act_f[0][i],n_act_f[0][i],1.0,
                            R_ACT_F[i],n_act[0],
                            V_ACT,n_act_f[0][i],0.0,
                            BUF,n_act_f[0][i]);
    //     memcpy(R_ACT_F[i],BUF,n_act[0]*n_ao[0]*sizeof(double));
        for(int j=0;j<n_ao [0]     ;j++)
        for(int k=0;k<n_act_f[0][i];k++){
            R_ACT_F[i][j*n_act[0]+k]=BUF[j*n_act_f[0][i]+k]/SVD[k];/*/sqrt(SVD[k])*/
        }
        
        double * V_ACT2 = new double[n_act_f[0][i]*n_act_f[0][i]];// transformation Vmatrix for malmqvist
        
        for(int j=0;j<n_act_f[0][i];j++)
        for(int k=j;k<n_act_f[0][i];k++){
            BUF[0]=V_ACT[j*n_act_f[0][i]+k]/SVD[j];
            V_ACT2[j*n_act_f[0][i]+k]=V_ACT[k*n_act_f[0][i]+j]/SVD[k];
            V_ACT2[k*n_act_f[0][i]+j]=BUF[0];
        }
        
        for(int j=0;j<n_act_f[0][i];j++)
        for(int k=j;k<n_act_f[0][i];k++){
            BUF[0]=V_ACT[j*n_act_f[0][i]+k]*SVD[j];
            V_ACT[j*n_act_f[0][i]+k]=V_ACT[k*n_act_f[0][i]+j]*SVD[k];////sqrt!!!!!!!
            V_ACT[k*n_act_f[0][i]+j]=BUF[0];
        }

        act_to_act_so_matr(L_O_V,S_SV,U_ACT,n_vt,n_act_f[0][i],n_act[0],1);
        act_to_act_so_matr(R_O_V,S_SV,V_ACT,n_vt,n_act_f[0][i],n_act[0],0);
        
        if(i>0){
            cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                            n_act_f[0][i],n_act_f[0][i]+n_vt,n_act_f[0][i],1.0,
                            U_ACT,     n_act_f[0][i],
                            L_O_V,n_vt+n_act_f[0][i],0.0,
                            L_O,  n_vt+n_act_f[0][i]);
            
            cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                            n_act_f[0][i],n_act_f[0][i]+n_vt,n_act_f[0][i],1.0,
                            V_ACT,     n_act_f[0][i], //// V_ACT inv!!!!!!!!!!!
                            R_O_V,n_vt+n_act_f[0][i],0.0,
                            R_O,  n_vt+n_act_f[0][i]);
            
            for(int j=0; j<n_act_f[0][i];j++)
                cblas_dscal(n_vt+n_act_f[0][i],1/(SVD[j]*SVD[j]),R_O+j*(n_vt+n_act_f[0][i]),1);
            
        }
        
        if(n_act_f[0][i]!=n_act_f[1][i]){
            fprintf(out_stream,"ERROR:different size of active space in fragment %d\n",i);
            exit(1);
        }
        
        aldet_bra[i].get_dim(n_act_f[0][i],n_alp_el_f[0][i],n_bet_el_f[0][i],1, 0, CAS_PRINT_NUMBER_DEFAULT);
        aldet_ket[i].get_dim(n_act_f[1][i],n_alp_el_f[1][i],n_bet_el_f[1][i],1, 0, CAS_PRINT_NUMBER_DEFAULT);

        aldet_bra[i].copy_coef(0, A->CI+i, A->n_states[i], 0, 1);
        aldet_ket[i].copy_coef(0, B->CI+i, B->n_states[i], 0, 1);
        
        aldet_bra[i].malmqvist(0,U_ACT);
        
        aldet_ket[i].malmqvist(0,V_ACT2);
        delete[] V_ACT2;
        
        if(n_act_f[0][i])
        cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                            n_act[0],n_ao[0],n_ao[0],1.0,
                            L_ACT,n_act[0],
                            M->S_MO,n_ao[0],0.0,
                            BUF,n_ao[0]);
        if(n_act[0])
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                            n_act[0],n_act[0],n_ao[0],1.0,
                            BUF,n_ao[0],
                            R_ACT,n_act[0],0.0,
                            S_SV,n_act[0]);

        
        
        n_vt+=n_act_f[0][i];
        
        delete [] L_O_V;
        delete [] R_O_V;
    }
    
    for(int i=0;i<n_act_f[0][1];i++)
    for(int j=0;j<n_act_f[0][0];j++)
        BUF[i*n_act_f[0][0]+j]=L_O[i*(n_act_f[0][1]+n_act_f[0][0])+j];
    
    ci_pr_aldet_first_order(&ci1_lt_aldet, aldet_bra+1, n_act_f[0][1], &BUF, n_alp_el_f[0][0], n_bet_el_f[0][0],A->n_states[1]);
    
    for(int i=0;i<n_act_f[0][1];i++)
    for(int j=0;j<n_act_f[0][1];j++)
        BUF[i*n_act_f[0][1]+j]=R_O[i*(n_act_f[0][1]+n_act_f[0][0])+j];

    //check alpha and beta
    ci_pr_aldet_first_order(&ci2_lt_aldet, aldet_ket+1, n_act_f[0][1], &BUF, n_alp_el_f[1][0], n_bet_el_f[1][0],B->n_states[1]);
    

    fprintf(out_stream,"making first product list\n");
    
    
    
    for(int i=0;i<ci1_lt_aldet.size();i++){
        ci1_lt_aldet[i].add_A(aldet_bra+0,0,0);
    }
    printf_timer("add_A1");
    
    fprintf(out_stream,"making second product list\n");
    
    
    for(int i=0;i<ci2_lt_aldet.size();i++){
        ci2_lt_aldet[i].add_A(aldet_ket+0,0,0);
    }


    /*for(int i=0;i<A->n_frag;i++)aldet_bra[i].clear();*/delete[]aldet_bra;
    /*for(int i=0;i<A->n_frag;i++)aldet_ket[i].clear();*/delete[]aldet_ket;
    
    delete[] L_O;
    delete[] R_O;
    
    
    return 0;
}

int doCI_data::link_and_transform_for_PT(molecule * M){
    
    
    int n_cor_A = A->n_cor_orb_f[0];
    int n_cor_B = A->n_cor_orb_f[1];
    
    
    for(int i=0; i<n_ao [0];i++)
    for(int j=0; j<n_act[0];j++)
        L_ACT[i*n_act[0]+j]=VEC[0][(j+n_cor[0])*n_ao[0]+i];
    
    for(int i=0; i<n_ao [0];i++)
    for(int j=0; j<n_act[0];j++)
        R_ACT[i*n_act[0]+j]=VEC[1][(j+n_cor[0])*n_ao[0]+i];
    

    double * S_DV = new double[n_cor[0]*(n_act_f[0][0]+n_act_f[0][1])];
    M->calc_S_AO();
    
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                        n_act[0],n_ao[0],n_ao[0],1.0,
                        L_ACT,n_act[0],//no change in lda!!
                        M->S_MO,n_ao[0],0.0,
                        BUF,n_ao[0]);
    
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        n_act[0],n_cor[0],n_ao[0],1.0,
                        BUF,n_ao[0],
                        R_COR,n_ao[0],0.0,//no change in lda!!
                        S_DV,n_cor[0]);
//     fprintf(out_stream,"S_DV:\n");
//     PrintMatr(S_DV,n_act[0],n_cor[0],1);
    
//     L_COR_1st_order = new double[n_cor[0]*n_ao[0]];
//     R_COR_1st_order = new double[n_cor[0]*n_ao[0]];
    
//     set_zero_matr(L_COR_1st_order, n_cor[0]*n_ao[0]);
//     set_zero_matr(R_COR_1st_order, n_cor[0]*n_ao[0]);
    
    for(int i=0; i<n_act_f[0][0];i++)
    for(int k=0; k<n_ao[0];k++)
    for(int j=0; j<n_cor_A;j++)
        L_ACT[k*n_act[0]+i]-=L_COR[j*n_ao[0]+k]*S_DV[i*n_cor[0]+j]/SVD[j];
    
    for(int i=n_act_f[0][0]; i<n_act_f[0][0]+n_act_f[0][1];i++)
    for(int k=0; k<n_ao[0];k++)
    for(int j=n_cor_A; j<n_cor_A+n_cor_B;j++)
        L_ACT[k*n_act[0]+i]-=L_COR[j*n_ao[0]+k]*S_DV[i*n_cor[0]+j]/SVD[j];
    
    
//     cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
//                         n_act[0],n_ao[0],n_ao[0],1.0,
//                         L_ACT,n_act[0],//no change in lda!!
//                         M->S_MO,n_ao[0],0.0,
//                         BUF,n_ao[0]);
//     
//     cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
//                         n_act[0],n_cor[0],n_ao[0],1.0,
//                         BUF,n_ao[0],
//                         R_COR,n_ao[0],0.0,//no change in lda!!
//                         S_DV,n_cor[0]);
//     fprintf(out_stream,"S_DV:\n");
//     PrintMatr(S_DV,n_act[0],n_cor[0],1);
    
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                        n_cor[0],n_ao[0],n_ao[0],1.0,
                        L_COR,n_ao[0],//no change in lda!!
                        M->S_MO,n_ao[0],0.0,
                        BUF,n_ao[0]);
    
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                        n_cor[0],n_act[0],n_ao[0],1.0,
                        BUF,n_ao[0],
                        R_ACT,n_act[0],0.0,//no change in lda!!
                        S_DV,n_act[0]);
//     fprintf(out_stream,"S_DV:\n");
//     PrintMatr(S_DV,n_cor[0],n_act[0],1);
    

    for(int i=0; i<n_act_f[0][0];i++)
    for(int k=0; k<n_ao[0];k++)
    for(int j=0; j<n_cor_A;j++)
        R_ACT[k*n_act[0]+i]-=R_COR[j*n_ao[0]+k]*S_DV[j*n_act[0]+i]/SVD[j];
    
    for(int i=n_act_f[0][0]; i<n_act_f[0][0]+n_act_f[0][1];i++)
    for(int k=0; k<n_ao[0];k++)
    for(int j=n_cor_A; j<n_cor_A+n_cor_B;j++)
        R_ACT[k*n_act[0]+i]-=R_COR[j*n_ao[0]+k]*S_DV[j*n_act[0]+i]/SVD[j];
    
    
//     cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
//                         n_cor[0],n_ao[0],n_ao[0],1.0,
//                         L_COR,n_ao[0],//no change in lda!!
//                         M->S_MO,n_ao[0],0.0,
//                         BUF,n_ao[0]);
//     
//     cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
//                         n_cor[0],n_act[0],n_ao[0],1.0,
//                         BUF,n_ao[0],
//                         R_ACT,n_act[0],0.0,//no change in lda!!
//                         S_DV,n_act[0]);
//     fprintf(out_stream,"S_DV:\n");
//     PrintMatr(S_DV,n_cor[0],n_act[0],1);
//     
// 

    
    M->calc_S_AO();
    
    if(n_act[0])
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                        n_act[0],n_ao[0],n_ao[0],1.0,
                        L_ACT,n_act[0],//no change in lda!!
                        M->S_MO,n_ao[0],0.0,
                        BUF,n_ao[0]);
    if(n_act[0])
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                        n_act[0],n_act[0],n_ao[0],1.0,
                        BUF,n_ao[0],
                        R_ACT,n_act[0],0.0,//no change in lda!!
                        S_SV,n_act[0]);
//     fprintf(out_stream,"S_SV:\n");
//     PrintMatr(S_SV,n_act[0],n_act[0],1);
    
    
    
    if(n_act[0])
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                        n_act[0],n_ao[0],n_ao[0],1.0,
                        L_ACT,n_act[0],//no change in lda!!
                        M->S_MO,n_ao[0],0.0,
                        BUF,n_ao[0]);
    if(n_act[0])
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                        n_act[0],n_act[0],n_ao[0],1.0,
                        BUF,n_ao[0],
                        R_ACT,n_act[0],0.0,//no change in lda!!
                        S_FV,n_act[0]);
//     fprintf(out_stream,"S_FV:\n");
//     PrintMatr(S_FV,n_act[0],n_act[0],1);
    
    
    if(A->n_frag>2){
        fprintf(out_stream,"Unable to use more then 2 fragments: n_frag = %d\n",A->n_frag);
        exit(1);
    }
    
    double * L_O = new double[n_act_f[0][1]*(n_act_f[0][0]+n_act_f[0][1])];
    double * R_O = new double[n_act_f[0][1]*(n_act_f[0][0]+n_act_f[0][1])];
                
    
    aldet_data * aldet_bra = new aldet_data[A->n_frag];
    aldet_data * aldet_ket = new aldet_data[A->n_frag];
    
    double * SVD_ACT=new double[n_act[0]];
    
    for(int i=0,n_vt=0; i<A->n_frag;i++){
        // calculation of active orbital overlap
//         if(i)fprintf(out_stream,"n_vt = %d\n",n_vt);
        L_ACT_F[i] = L_ACT + n_vt /** n_act[0]*/;
        R_ACT_F[i] = R_ACT + n_vt /** n_act[0]*/;
        double * L_O_V = new double[n_act_f[0][i]*(n_vt+n_act_f[0][i])];
        double * R_O_V = new double[n_act_f[0][i]*(n_vt+n_act_f[0][i])];
        
        
//         double * L_FV = new double[n_act_f[0][i]*n_act_f[0][i]];
//         double * R_FV = new double[n_act_f[0][i]*n_act_f[0][i]];
        
//         if(act_to_act_ort(L_O_V,S_SV,n_vt,n_act_f[0][i],n_act[0],1)){
//             cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
//                                 n_ao[0],n_act_f[0][i],n_act_f[0][i]+n_vt,1.0,
//                                 L_ACT,n_act[0],
//                                 L_O_V,n_act_f[0][i]+n_vt,0.0,
//                                 BUF,n_act_f[0][i]);
//             for(int j=0;j<n_ao[0];j++)memcpy(L_ACT_F[i]+j*n_act[0],BUF+j*n_act_f[0][i],n_act_f[0][i]*sizeof(double));
//             
//         }
//         if(act_to_act_ort(R_O_V,S_SV,n_vt,n_act_f[0][i],n_act[0],0)){
//             cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
//                                 n_ao[0],n_act_f[0][i],n_act_f[0][i]+n_vt,1.0,
//                                 R_ACT,n_act[0],
//                                 R_O_V,n_act_f[0][i]+n_vt,0.0,
//                                 BUF,n_act_f[0][i]);
//             for(int j=0;j<n_ao[0];j++)memcpy(R_ACT_F[i]+j*n_act[0],BUF+j*n_act_f[0][i],n_act_f[0][i]*sizeof(double));
//             
//         }
        
        if(n_act_f[0][i])
        cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                            n_act_f[0][i],n_ao[0],n_ao[0],1.0,
                            L_ACT_F[i],n_act[0],//no change in lda!!
                            M->S_MO,n_ao[0],0.0,
                            BUF,n_ao[0]);
        if(n_act_f[0][i])
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                            n_act_f[0][i],n_act_f[0][i],n_ao[0],1.0,
                            BUF,n_ao[0],
                            R_ACT_F[i],n_act[0],0.0,//no change in lda!!
                            S_FV,n_act_f[0][i]);
        
//         fprintf(out_stream,"S_FV(%d)\n",i);
//         PrintMatr(S_FV,n_act_f[0][i],n_act_f[0][i],1);
//         exit(0);
        lapack_svd(S_FV,V_ACT,U_ACT,SVD_ACT,n_act_f[0][i]);

                
        //transformation of orthogonalized active orbitals to singular active orbitals
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                            n_ao[0],n_act_f[0][i],n_act_f[0][i],1.0,
                            L_ACT_F[i],n_act[0],
                            U_ACT,n_act_f[0][i],0.0,
                            BUF,n_act_f[0][i]);
        for(int j=0;j<n_ao [0]     ;j++)
        for(int k=0;k<n_act_f[0][i];k++){
            L_ACT_F[i][j*n_act[0]+k]=BUF[j*n_act_f[0][i]+k];///sqrt(SVD_ACT[k]);
        }
        
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            n_ao[0],n_act_f[0][i],n_act_f[0][i],1.0,
                            R_ACT_F[i],n_act[0],
                            V_ACT,n_act_f[0][i],0.0,
                            BUF,n_act_f[0][i]);
    //     memcpy(R_ACT_F[i],BUF,n_act[0]*n_ao[0]*sizeof(double));
        for(int j=0;j<n_ao [0]     ;j++)
        for(int k=0;k<n_act_f[0][i];k++){
            R_ACT_F[i][j*n_act[0]+k]=BUF[j*n_act_f[0][i]+k]/SVD_ACT[k];/*/sqrt(SVD_ACT[k])*/
        }
        
        double * V_ACT2 = new double[n_act_f[0][i]*n_act_f[0][i]];// transformation Vmatrix for malmqvist
        
        for(int j=0;j<n_act_f[0][i];j++)
        for(int k=j;k<n_act_f[0][i];k++){
            BUF[0]=V_ACT[j*n_act_f[0][i]+k]/SVD_ACT[j];
            V_ACT2[j*n_act_f[0][i]+k]=V_ACT[k*n_act_f[0][i]+j]/SVD_ACT[k];
            V_ACT2[k*n_act_f[0][i]+j]=BUF[0];
        }
        
        for(int j=0;j<n_act_f[0][i];j++)
        for(int k=j;k<n_act_f[0][i];k++){
            BUF[0]=V_ACT[j*n_act_f[0][i]+k]*SVD_ACT[j];
            V_ACT[j*n_act_f[0][i]+k]=V_ACT[k*n_act_f[0][i]+j]*SVD_ACT[k];////sqrt!!!!!!!
            V_ACT[k*n_act_f[0][i]+j]=BUF[0];
        }

        act_to_act_so_matr(L_O_V,S_SV,U_ACT,n_vt,n_act_f[0][i],n_act[0],1);
        act_to_act_so_matr(R_O_V,S_SV,V_ACT,n_vt,n_act_f[0][i],n_act[0],0);
        
        if(i>0){
            cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                            n_act_f[0][i],n_act_f[0][i]+n_vt,n_act_f[0][i],1.0,
                            U_ACT,     n_act_f[0][i],
                            L_O_V,n_vt+n_act_f[0][i],0.0,
                            L_O,  n_vt+n_act_f[0][i]);
            
            cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                            n_act_f[0][i],n_act_f[0][i]+n_vt,n_act_f[0][i],1.0,
                            V_ACT,     n_act_f[0][i], //// V_ACT inv!!!!!!!!!!!
                            R_O_V,n_vt+n_act_f[0][i],0.0,
                            R_O,  n_vt+n_act_f[0][i]);
            
            for(int j=0; j<n_act_f[0][i];j++)
                cblas_dscal(n_vt+n_act_f[0][i],1/(SVD_ACT[j]*SVD_ACT[j]),R_O+j*(n_vt+n_act_f[0][i]),1);
            
        }
        
        if(n_act_f[0][i]!=n_act_f[1][i]){
            fprintf(out_stream,"ERROR:different size of active space in fragment %d\n",i);
            exit(1);
        }
        
        aldet_bra[i].get_dim(n_act_f[0][i],n_alp_el_f[0][i],n_bet_el_f[0][i],1, 0, CAS_PRINT_NUMBER_DEFAULT);
        aldet_ket[i].get_dim(n_act_f[1][i],n_alp_el_f[1][i],n_bet_el_f[1][i],1, 0, CAS_PRINT_NUMBER_DEFAULT);

        fprintf(out_stream,"READING FROM MOLECULE IS TURNED OFF in NOPT 1.2\n");
        fprintf(out_stream,"CONTINUE writting code!.... pleeeeease :(\n");
        exit(0);
//         aldet_bra[i].read_set_from_mol(A, 0, i, A->n_states[i], 0);
//         aldet_ket[i].read_set_from_mol(B, 0, i, B->n_states[i], 0);
        
        aldet_bra[i].malmqvist(0,U_ACT);
        
        aldet_ket[i].malmqvist(0,V_ACT2);
        delete[] V_ACT2;
        
        if(n_act_f[0][i])
        cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                            n_act[0],n_ao[0],n_ao[0],1.0,
                            L_ACT,n_act[0],
                            M->S_MO,n_ao[0],0.0,
                            BUF,n_ao[0]);
        if(n_act[0])
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                            n_act[0],n_act[0],n_ao[0],1.0,
                            BUF,n_ao[0],
                            R_ACT,n_act[0],0.0,
                            S_SV,n_act[0]);

//         fprintf(out_stream,"S_SV final:\n");
//         PrintMatr(S_SV,n_act[0],n_act[0],0);
    
        
        n_vt+=n_act_f[0][i];
        
        delete [] L_O_V;
        delete [] R_O_V;
    }
    
    /// make ci_pr_aldet_first_order
    
    for(int i=0;i<n_act_f[0][1];i++)
    for(int j=0;j<n_act_f[0][0];j++)
        BUF[i*n_act_f[0][0]+j]=L_O[i*(n_act_f[0][1]+n_act_f[0][0])+j];
    
    ci_pr_aldet_first_order(&ci1_lt_aldet, aldet_bra+1, n_act_f[0][1], &BUF, n_alp_el_f[0][0], n_bet_el_f[0][0],A->n_states[1]);
    
    for(int i=0;i<n_act_f[0][1];i++)
    for(int j=0;j<n_act_f[0][1];j++)
        BUF[i*n_act_f[0][1]+j]=R_O[i*(n_act_f[0][1]+n_act_f[0][0])+j];

    //check alpha and beta
    ci_pr_aldet_first_order(&ci2_lt_aldet, aldet_ket+1, n_act_f[0][1], &BUF, n_alp_el_f[1][0], n_bet_el_f[1][0],B->n_states[1]);
    

    fprintf(out_stream,"making first product list\n");
    
    
    
    for(int i=0;i<ci1_lt_aldet.size();i++){
        ci1_lt_aldet[i].add_A(aldet_bra+0,0,0);
    }
    printf_timer("add_A1");
    
    fprintf(out_stream,"making second product list\n");
    
    
    for(int i=0;i<ci2_lt_aldet.size();i++){
        ci2_lt_aldet[i].add_A(aldet_ket+0,0,0);
    }

    
    
    
    
    L_ACT_1st_order = new double[n_ao[0]*n_act[0]];
    R_ACT_1st_order = new double[n_ao[0]*n_act[0]];
    
    set_zero_matr(L_ACT_1st_order, n_ao[0]*n_act[0]);
    set_zero_matr(R_ACT_1st_order, n_ao[0]*n_act[0]);
    
//     for(int i=0; i<n_act_f[0][0];i++)
//     for(int k=0; k<n_ao[0];k++)
//     for(int j=0; j<n_cor_A;j++)
//         L_ACT[k*n_act[0]+i]-=L_COR[j*n_ao[0]+k]*S_DV[i*n_cor[0]+j]/SVD[j];
//     
    for(int i=n_act_f[0][0]; i<n_act_f[0][0]+n_act_f[0][1];i++)
    for(int k=0; k<n_ao[0];k++)
    for(int j=0; j<n_act_f[0][0];j++)
        L_ACT_1st_order[k*n_act[0]+i]-=L_ACT[k*n_act[0]+j]*S_SV[i*n_act[0]+j];
    
    for(int i=n_act_f[0][0]; i<n_act_f[0][0]+n_act_f[0][1];i++)
    for(int k=0; k<n_ao[0];k++)
    for(int j=0; j<n_act_f[0][0];j++)
        R_ACT_1st_order[k*n_act[0]+i]-=R_ACT[k*n_act[0]+j]*S_SV[j*n_act[0]+i];
    
    
//     cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
//                         n_act[0],n_ao[0],n_ao[0],1.0,
//                         L_ACT,n_act[0],
//                         M->S_MO,n_ao[0],0.0,
//                         BUF,n_ao[0]);
//     if(n_act[0])
//     cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
//                         n_act[0],n_act[0],n_ao[0],1.0,
//                         BUF,n_ao[0],
//                         R_ACT_1st_order,n_act[0],0.0,
//                         S_SV,n_act[0]);
// 
//     fprintf(out_stream,"S_V_L(0)R(1):\n");
//     PrintMatr(S_SV,n_act[0],n_act[0],1);

    M->calc_S_AO();
    
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                        n_act[0],n_ao[0],n_ao[0],1.0,
                        L_ACT,n_act[0],//no change in lda!!
                        M->S_MO,n_ao[0],0.0,
                        BUF,n_ao[0]);
    
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        n_act[0],n_cor[0],n_ao[0],1.0,
                        BUF,n_ao[0],
                        R_COR,n_ao[0],0.0,//no change in lda!!
                        S_DV,n_cor[0]);

//     fprintf(out_stream,"S_DV:\n");
//     PrintMatr(S_DV,n_act[0],n_cor[0],1);

    for(int i=0; i<n_act_f[0][0];i++)
    for(int k=0; k<n_ao[0];k++)
    for(int j=n_cor_A; j<n_cor_A+n_cor_B;j++)
        L_ACT_1st_order[k*n_act[0]+i]-=L_COR[j*n_ao[0]+k]*S_DV[i*n_cor[0]+j]/SVD[j];
    
    for(int i=n_act_f[0][0]; i<n_act_f[0][0]+n_act_f[0][1];i++)
    for(int k=0; k<n_ao[0];k++)
    for(int j=0; j<n_cor_A;j++)
        L_ACT_1st_order[k*n_act[0]+i]-=L_COR[j*n_ao[0]+k]*S_DV[i*n_cor[0]+j]/SVD[j];
    
//     cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
//                         n_act[0],n_ao[0],n_ao[0],1.0,
//                         L_ACT_1st_order,n_act[0],//no change in lda!!
//                         M->S_MO,n_ao[0],0.0,
//                         BUF,n_ao[0]);
//     
//     cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
//                         n_act[0],n_cor[0],n_ao[0],1.0,
//                         BUF,n_ao[0],
//                         R_COR,n_ao[0],0.0,//no change in lda!!
//                         S_DV,n_cor[0]);
//     fprintf(out_stream,"S_DV:\n");
//     PrintMatr(S_DV,n_act[0],n_cor[0],1);


    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                        n_cor[0],n_ao[0],n_ao[0],1.0,
                        L_COR,n_ao[0],//no change in lda!!
                        M->S_MO,n_ao[0],0.0,
                        BUF,n_ao[0]);
    
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                        n_cor[0],n_act[0],n_ao[0],1.0,
                        BUF,n_ao[0],
                        R_ACT,n_act[0],0.0,//no change in lda!!
                        S_DV,n_act[0]);
//     fprintf(out_stream,"S_DV:\n");
//     PrintMatr(S_DV,n_cor[0],n_act[0],1);
//     exit(0);

    for(int i=0; i<n_act_f[0][0];i++)
    for(int k=0; k<n_ao[0];k++)
    for(int j=n_cor_A; j<n_cor_A+n_cor_B;j++)
        R_ACT_1st_order[k*n_act[0]+i]-=R_COR[j*n_ao[0]+k]*S_DV[j*n_act[0]+i]/SVD[j];
    
    for(int i=n_act_f[0][0]; i<n_act_f[0][0]+n_act_f[0][1];i++)
    for(int k=0; k<n_ao[0];k++)
    for(int j=0; j<n_cor_A;j++)
        R_ACT_1st_order[k*n_act[0]+i]-=R_COR[j*n_ao[0]+k]*S_DV[j*n_act[0]+i]/SVD[j];
    
    
//     cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
//                         n_cor[0],n_ao[0],n_ao[0],1.0,
//                         L_COR,n_ao[0],//no change in lda!!
//                         M->S_MO,n_ao[0],0.0,
//                         BUF,n_ao[0]);
    
//     cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
//                         n_cor[0],n_act[0],n_ao[0],1.0,
//                         BUF,n_ao[0],
//                         R_ACT_1st_order,n_act[0],0.0,//no change in lda!!
//                         S_DV,n_act[0]);
//     fprintf(out_stream,"S_DV:\n");
//     PrintMatr(S_DV,n_cor[0],n_act[0],1);








    /*for(int i=0;i<A->n_frag;i++)aldet_bra[i].clear();*/delete[]aldet_bra;
    /*for(int i=0;i<A->n_frag;i++)aldet_ket[i].clear();*/delete[]aldet_ket;
    
    delete[] L_O;
    delete[] R_O;
    
    
    
    
    
    return 0;
}

int doCI_data::calc_CI_for_dec_wf(){
    
    double * act_c_INTS_A=new double[n_act_f[0][0]*n_act_f[0][0]*n_act_f[0][0]*n_act_f[0][0]];
    double * act_c_INTS_B=new double[n_act_f[0][1]*n_act_f[0][1]*n_act_f[0][1]*n_act_f[0][1]];
    
    for(int i=0; i<n_act_f[0][0];i++)
    for(int j=0; j<n_act_f[0][0];j++)
    for(int k=0; k<n_act_f[0][0];k++)
    for(int l=0; l<n_act_f[0][0];l++)
        act_c_INTS_A[((i*n_act_f[0][0]+j)*n_act_f[0][0]+k)*n_act_f[0][0]+l]=
        act_INTS    [((i*n_act  [0]   +j)*n_act  [0]   +k)*n_act  [0]   +l];
    
    for(int i=0,ii=n_act_f[0][0]; i<n_act_f[0][1];i++,ii++)
    for(int j=0,jj=n_act_f[0][0]; j<n_act_f[0][1];j++,jj++)
    for(int k=0,kk=n_act_f[0][0]; k<n_act_f[0][1];k++,kk++)
    for(int l=0,ll=n_act_f[0][0]; l<n_act_f[0][1];l++,ll++)
        act_c_INTS_B[((i *n_act_f[0][0]+j )*n_act_f[0][0]+k )*n_act_f[0][0]+l ]=
        act_INTS    [((ii*n_act  [0]   +jj)*n_act  [0]   +kk)*n_act  [0]   +ll];
    
    double * s1_ext_vec;
    double * s2_ext_vec;
    double * v1_ext_vec;
    double * v2_ext_vec;
    double * v1b_ext_vec;
    double * v2b_ext_vec;
    double * s1_ext_mat;
    double * s2_ext_mat;
    
    double sign = 1.0;
    
    ACT_DMT_A = new double[n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1])];
    ACT_DMT_B = new double[n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1])];
    s1_ext_vec= new double[n_act[0]         *(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1])];
    s2_ext_vec= new double[         n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1])];
    v1_ext_vec= new double[n_act[0]         *(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1])];
    v2_ext_vec= new double[         n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1])];
    v1b_ext_vec= new double[n_act[0]         *(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1])];
    v2b_ext_vec= new double[         n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1])];
    s1_ext_mat= new double[n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1])];
    s2_ext_mat= new double[n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1])];
    
    double * ACT_DET_A2= ACT_DET_A+  A->n_states[0]*B->n_states[0];
    
    
    set_zero_matr(h_2e     ,                  A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]);
    set_zero_matr(ACT_DET_A,                  A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]);
    set_zero_matr(ACT_ADJ_A,n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
    set_zero_matr(ACT_ADJ_B,n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
    
    set_zero_matr(ACT_DMT_A,n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1]));
    set_zero_matr(ACT_DMT_B,n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1]));
    set_zero_matr(s1_ext_vec,n_act[0]         *(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
    set_zero_matr(s2_ext_vec,         n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
    set_zero_matr(v1_ext_vec,n_act[0]         *(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
    set_zero_matr(v2_ext_vec,         n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
    set_zero_matr(v1b_ext_vec,n_act[0]         *(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
    set_zero_matr(v2b_ext_vec,         n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
    
    if(S_new==NULL)S_new= new double[A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1]];
    if(H_2el==NULL)H_2el= new double[A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1]];
    
    set_zero_matr(S_new,A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1]);
    set_zero_matr(H_2el,A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1]);
    
    
//     double H_2el_new = 0;
    
    
    double * TH_1A = ACT_ADJ_A+(/*i*n_st_total[1]+j */0                          )*n_act[0]*n_act[0];
    double * TH_1B = ACT_ADJ_B+(/*i*n_st_total[1]+j */0                          )*n_act[0]*n_act[0];
    double * TH_2A = ACT_ADJ_A+(/*i*n_st_total[1]+j*/+A->n_states[0]*B->n_states[0])*n_act[0]*n_act[1];
    double * TH_2B = ACT_ADJ_B+(/*i*n_st_total[1]+j*/+A->n_states[0]*B->n_states[0])*n_act[0]*n_act[1];
    
    int abcd;
    int ac;
    int bd;
    int mc_2 = n_act_f[0][0]*n_act[1]+n_act_f[0][0]; // beginning of 2-2 block
    
//     getchar();
    double * act_INTS_B=act_INTS+n_act[0]*n_act[0]*n_act[0]*n_act_f[0][1]+n_act[0]*n_act[0]*n_act_f[0][1]+n_act[0]*n_act_f[0][1]+n_act_f[0][1];
    for(int i=0;i<ci1_lt_aldet.size();i++){
    for(int j=0;j<ci2_lt_aldet.size();j++){

//         if(false)//0
        if(ci1_lt_aldet[i].n_A==ci2_lt_aldet[j].n_A)
        if(ci1_lt_aldet[i].n_B==ci2_lt_aldet[j].n_B){
            aldet_S_calc(ACT_DET_A ,&(ci1_lt_aldet[i].A), 0, &(ci2_lt_aldet[j].A), 0);
            aldet_S_calc(ACT_DET_A2 ,&(ci1_lt_aldet[i].B), 0, &(ci2_lt_aldet[j].B), 0);
            
            calc_2e_ints(h_2e, ci1_lt_aldet[i].A.n_act,ci1_lt_aldet[i].A.na,ci1_lt_aldet[i].A.nb, ci1_lt_aldet[i].A.coef[0], ci2_lt_aldet[j].A.coef[0], act_c_INTS_A,0,ci1_lt_aldet[i].A.n_states[0],ci1_lt_aldet[i].A.n_states[0],ci2_lt_aldet[j].A.n_states[0],ci2_lt_aldet[j].A.n_states[0]);
            
            calc_2e_ints(h_2e+A->n_states[0]*B->n_states[0], ci1_lt_aldet[i].B.n_act,ci1_lt_aldet[i].B.na,ci1_lt_aldet[i].B.nb, ci1_lt_aldet[i].B.coef[0], ci2_lt_aldet[j].B.coef[0], act_c_INTS_B,0,ci1_lt_aldet[i].B.n_states[0],ci1_lt_aldet[i].B.n_states[0],ci2_lt_aldet[j].B.n_states[0],ci2_lt_aldet[j].B.n_states[0]);
            
            aldet_DMA_calc(ACT_ADJ_A ,&(ci1_lt_aldet[i].A), 0, &(ci2_lt_aldet[j].A), 0, n_act[0]);
            aldet_DMB_calc(ACT_ADJ_B ,&(ci1_lt_aldet[i].A), 0, &(ci2_lt_aldet[j].A), 0, n_act[0]);

            aldet_DMA_calc(ACT_ADJ_A+A->n_states[0]*B->n_states[0]*n_act[0]*n_act[1] ,&(ci1_lt_aldet[i].B), 0, &(ci2_lt_aldet[j].B), 0, n_act[0]);
            aldet_DMB_calc(ACT_ADJ_B+A->n_states[0]*B->n_states[0]*n_act[0]*n_act[1] ,&(ci1_lt_aldet[i].B), 0, &(ci2_lt_aldet[j].B), 0, n_act[0]);
            
            for(int a=0;a<A->n_states[0];a++)
            for(int b=0;b<A->n_states[1];b++)
            for(int c=0;c<B->n_states[0];c++)
            for(int d=0;d<B->n_states[1];d++){
                abcd = ((a*A->n_states[1]+b)*B->n_states[0]+c)*B->n_states[1]+d;
                ac   =   a*B->n_states[0]+c;
                bd   =   b*B->n_states[1]+d;
                
                for(int ii=0;ii<n_act[0]*n_act[1];ii++)
                    ACT_DMT_A[ii+     abcd*n_act[0]*n_act[1]]+=TH_1A[ii+ac*n_act[0]*n_act[1]]*ACT_DET_A2[bd];
                                   
                for(int ii=0;ii<n_act[0]*n_act[1]-mc_2;ii++)
                    ACT_DMT_A[ii+mc_2+abcd*n_act[0]*n_act[1]]+=TH_2A[ii+bd*n_act[0]*n_act[1]]*ACT_DET_A[ac];
                
                for(int ii=0;ii<n_act[0]*n_act[1];ii++)
                    ACT_DMT_B[ii+     abcd*n_act[0]*n_act[1]]+=TH_1B[ii+ac*n_act[0]*n_act[1]]*ACT_DET_A2[bd];
                                   
                for(int ii=0;ii<n_act[0]*n_act[1]-mc_2;ii++)
                    ACT_DMT_B[ii+mc_2+abcd*n_act[0]*n_act[1]]+=TH_2B[ii+bd*n_act[0]*n_act[1]]*ACT_DET_A[ac];
                
                H_2el[abcd]+= V_AB_22_J_calc(act_INTS,TH_1A+ac*n_act[0]*n_act[1], TH_2A+bd*n_act[0]*n_act[1], 0, n_act_f[0][0], n_act_f[0][0], n_act_f[0][0]+n_act_f[0][1], n_act[0])*p_SVD;
                H_2el[abcd]+= V_AB_22_J_calc(act_INTS,TH_1A+ac*n_act[0]*n_act[1], TH_2B+bd*n_act[0]*n_act[1], 0, n_act_f[0][0], n_act_f[0][0], n_act_f[0][0]+n_act_f[0][1], n_act[0])*p_SVD;
                H_2el[abcd]+= V_AB_22_J_calc(act_INTS,TH_1B+ac*n_act[0]*n_act[1], TH_2A+bd*n_act[0]*n_act[1], 0, n_act_f[0][0], n_act_f[0][0], n_act_f[0][0]+n_act_f[0][1], n_act[0])*p_SVD;
                H_2el[abcd]+= V_AB_22_J_calc(act_INTS,TH_1B+ac*n_act[0]*n_act[1], TH_2B+bd*n_act[0]*n_act[1], 0, n_act_f[0][0], n_act_f[0][0], n_act_f[0][0]+n_act_f[0][1], n_act[0])*p_SVD;
                H_2el[abcd]-= V_AB_22_K_calc(act_INTS,TH_1A+ac*n_act[0]*n_act[1], TH_2A+bd*n_act[0]*n_act[1], 0, n_act_f[0][0], n_act_f[0][0], n_act_f[0][0]+n_act_f[0][1], n_act[0])*p_SVD;
                H_2el[abcd]-= V_AB_22_K_calc(act_INTS,TH_1B+ac*n_act[0]*n_act[1], TH_2B+bd*n_act[0]*n_act[1], 0, n_act_f[0][0], n_act_f[0][0], n_act_f[0][0]+n_act_f[0][1], n_act[0])*p_SVD;
                
                
                S_new  [abcd] += ACT_DET_A [ac]*ACT_DET_A2[bd];
                
                H_2el  [abcd] += ACT_DET_A [ac]*h_2e      [bd+A->n_states[0]*B->n_states[0]]*p_SVD;
                
                H_2el  [abcd] += h_2e      [ac]*ACT_DET_A2[bd]*p_SVD;
            }
        }
//         if(false)//1
        if(ci1_lt_aldet[i].n_A==(ci2_lt_aldet[j].n_A+1))
        if(ci1_lt_aldet[i].n_B== ci2_lt_aldet[j].n_B  ){
//             //fprintf(out_stream,"into 1\n");
            
            sign = pow(-1,ci1_lt_aldet[i].n_A-1);

            set_zero_matr(v1_ext_vec, n_act[0]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            set_zero_matr(v2_ext_vec, n_act[0]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            set_zero_matr(v1b_ext_vec,n_act[0]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            set_zero_matr(v2b_ext_vec,n_act[0]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            
            aldet_calc_CI_2el_AAA0(v1_ext_vec , &(ci1_lt_aldet[i].A), 0, &(ci2_lt_aldet[j].A), 0, act_INTS, 0, n_act_f[0][0], n_act[0],        A->n_states[0],B->n_states[0],0);
            aldet_calc_CI_2el_AB0B(v1b_ext_vec, &(ci1_lt_aldet[i].A), 0, &(ci2_lt_aldet[j].A), 0, act_INTS, 0, n_act_f[0][0], n_act[0],        A->n_states[0],B->n_states[0],0);
            aldet_calc_CI_2el_A0AA(v2_ext_vec , &(ci1_lt_aldet[i].B), 0, &(ci2_lt_aldet[j].B), 0, act_INTS, n_act_f[0][0], n_act[0], n_act[0], A->n_states[1],B->n_states[1],0);
            aldet_calc_CI_2el_0BAB(v2b_ext_vec, &(ci1_lt_aldet[i].B), 0, &(ci2_lt_aldet[j].B), 0, act_INTS, n_act_f[0][0], n_act[0], n_act[0], A->n_states[1],B->n_states[1],0);
            
            aldet_gen_1_el_vec_arr  (s1_ext_vec,&(ci1_lt_aldet[i].A),0,A->n_states[0],&(ci2_lt_aldet[j].A),0,B->n_states[0], 0, n_act_f[0][0], n_act[0], 0);
            aldet_gen_1_el_vec_m_arr(s2_ext_vec,&(ci1_lt_aldet[i].B),0,A->n_states[1],&(ci2_lt_aldet[j].B),0,B->n_states[1], n_act_f[0][0], n_act[0], n_act[0], 0);
            
            for(int a=0;a<A->n_states[0];a++)
            for(int b=0;b<A->n_states[1];b++)
            for(int c=0;c<B->n_states[0];c++)
            for(int d=0;d<B->n_states[1];d++){
                abcd = ((a*A->n_states[1]+b)*B->n_states[0]+c)*B->n_states[1]+d;
                ac   =   a*B->n_states[0]+c;
                bd   =   b*B->n_states[1]+d;
            
               for(int ii=0;ii<n_act[0];ii++)
               for(int jj=0;jj<n_act[0];jj++)
                   ACT_DMT_A[ii*n_act[0]+jj+abcd*n_act[0]*n_act[1]]+=s1_ext_vec[ii+ac*n_act[0]]*
                                                                     s2_ext_vec[jj+bd*n_act[0]]*sign;
               
               for(int ii=0;ii<n_act[0];ii++)H_2el[abcd] -=  v1_ext_vec[ii+ac*n_act[0]]*s2_ext_vec[ii+bd*n_act[0]]*p_SVD*sign;
               for(int ii=0;ii<n_act[0];ii++)H_2el[abcd] += v1b_ext_vec[ii+ac*n_act[0]]*s2_ext_vec[ii+bd*n_act[0]]*p_SVD*sign;
               for(int ii=0;ii<n_act[0];ii++)H_2el[abcd] -=  v2_ext_vec[ii+bd*n_act[0]]*s1_ext_vec[ii+ac*n_act[0]]*p_SVD*sign;
               for(int ii=0;ii<n_act[0];ii++)H_2el[abcd] += v2b_ext_vec[ii+bd*n_act[0]]*s1_ext_vec[ii+ac*n_act[0]]*p_SVD*sign;
            }    
           
        }
//         if(false)//2
        if(ci1_lt_aldet[i].n_A== ci2_lt_aldet[j].n_A   )
        if(ci1_lt_aldet[i].n_B==(ci2_lt_aldet[j].n_B+1)){
            //fprintf(out_stream,"into 2\n");
            sign = pow(-1,ci1_lt_aldet[i].n_B-1);
            set_zero_matr(v1_ext_vec, n_act[0]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            set_zero_matr(v2_ext_vec, n_act[0]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            set_zero_matr(v1b_ext_vec,n_act[0]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            set_zero_matr(v2b_ext_vec,n_act[0]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            
            aldet_calc_CI_2el_AAA0(v1_ext_vec, &(ci1_lt_aldet[i].A),0, &(ci2_lt_aldet[j].A),0, act_INTS, 0, n_act_f[0][0], n_act[0],        A->n_states[0],B->n_states[0],1);
            aldet_calc_CI_2el_AB0B(v1b_ext_vec,&(ci1_lt_aldet[i].A),0, &(ci2_lt_aldet[j].A),0, act_INTS, 0, n_act_f[0][0], n_act[0],        A->n_states[0],B->n_states[0],1);
            aldet_calc_CI_2el_A0AA(v2_ext_vec, &(ci1_lt_aldet[i].B),0, &(ci2_lt_aldet[j].B),0, act_INTS, n_act_f[0][0], n_act[0], n_act[0], A->n_states[1],B->n_states[1],1);
            aldet_calc_CI_2el_0BAB(v2b_ext_vec,&(ci1_lt_aldet[i].B),0, &(ci2_lt_aldet[j].B),0, act_INTS, n_act_f[0][0], n_act[0], n_act[0], A->n_states[1],B->n_states[1],1);
            
            aldet_gen_1_el_vec_arr  (s1_ext_vec,&(ci1_lt_aldet[i].A),0,A->n_states[0],&(ci2_lt_aldet[j].A),0,B->n_states[0], 0, n_act_f[0][0], n_act[0],1);
            aldet_gen_1_el_vec_m_arr(s2_ext_vec,&(ci1_lt_aldet[i].B),0,A->n_states[1],&(ci2_lt_aldet[j].B),0,B->n_states[1], n_act_f[0][0], n_act[0], n_act[0],1);
                        
            for(int a=0;a<A->n_states[0];a++)
            for(int b=0;b<A->n_states[1];b++)
            for(int c=0;c<B->n_states[0];c++)
            for(int d=0;d<B->n_states[1];d++){
                abcd = ((a*A->n_states[1]+b)*B->n_states[0]+c)*B->n_states[1]+d;
                ac   =   a*B->n_states[0]+c;
                bd   =   b*B->n_states[1]+d;

                for(int ii=0;ii<n_act[0];ii++)
                for(int jj=0;jj<n_act[0];jj++)
                    ACT_DMT_B[ii*n_act[0]+jj+abcd*n_act[0]*n_act[1]]+=s1_ext_vec[ii+ac*n_act[0]]*
                                                                      s2_ext_vec[jj+bd*n_act[0]]*sign;
                
                for(int ii=0;ii<n_act[0];ii++)H_2el[abcd] -=  v1_ext_vec[ii+ac*n_act[0]]*s2_ext_vec[ii+bd*n_act[0]]*p_SVD*sign;
                for(int ii=0;ii<n_act[0];ii++)H_2el[abcd] += v1b_ext_vec[ii+ac*n_act[0]]*s2_ext_vec[ii+bd*n_act[0]]*p_SVD*sign;
                for(int ii=0;ii<n_act[0];ii++)H_2el[abcd] -=  v2_ext_vec[ii+bd*n_act[0]]*s1_ext_vec[ii+ac*n_act[0]]*p_SVD*sign;
                for(int ii=0;ii<n_act[0];ii++)H_2el[abcd] += v2b_ext_vec[ii+bd*n_act[0]]*s1_ext_vec[ii+ac*n_act[0]]*p_SVD*sign;
            }
            
        }
//         if(false)//3
        if(ci1_lt_aldet[i].n_A==(ci2_lt_aldet[j].n_A-1))
        if(ci1_lt_aldet[i].n_B== ci2_lt_aldet[j].n_B  ){
//             fprintf(out_stream,"into 3 %d %d\n",i,j);
            sign = pow(-1,ci1_lt_aldet[i].n_A);
            
            set_zero_matr(v1_ext_vec, n_act[0]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            set_zero_matr(v2_ext_vec, n_act[0]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            set_zero_matr(v1b_ext_vec,n_act[0]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            set_zero_matr(v2b_ext_vec,n_act[0]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            
            aldet_calc_CI_2el_A0AA(v1_ext_vec, &(ci1_lt_aldet[i].A),0, &(ci2_lt_aldet[j].A),0, act_INTS, 0, n_act_f[0][0], n_act[0],        A->n_states[0],B->n_states[0],0);
            aldet_calc_CI_2el_0BAB(v1b_ext_vec,&(ci1_lt_aldet[i].A),0, &(ci2_lt_aldet[j].A),0, act_INTS, 0, n_act_f[0][0], n_act[0],        A->n_states[0],B->n_states[0],0);
            aldet_calc_CI_2el_AAA0(v2_ext_vec, &(ci1_lt_aldet[i].B),0, &(ci2_lt_aldet[j].B),0, act_INTS, n_act_f[0][0], n_act[0], n_act[0], A->n_states[1],B->n_states[1],0);
            aldet_calc_CI_2el_AB0B(v2b_ext_vec,&(ci1_lt_aldet[i].B),0, &(ci2_lt_aldet[j].B),0, act_INTS, n_act_f[0][0], n_act[0], n_act[0], A->n_states[1],B->n_states[1],0);
            
            aldet_gen_1_el_vec_m_arr(s1_ext_vec,&(ci1_lt_aldet[i].A),0,A->n_states[0],&(ci2_lt_aldet[j].A),0,B->n_states[0], 0, n_act_f[0][0], n_act[0], 0);
            aldet_gen_1_el_vec_arr  (s2_ext_vec,&(ci1_lt_aldet[i].B),0,A->n_states[1],&(ci2_lt_aldet[j].B),0,B->n_states[1], n_act_f[0][0], n_act[0], n_act[0], 0);
            
            
            for(int a=0;a<A->n_states[0];a++)
            for(int b=0;b<A->n_states[1];b++)
            for(int c=0;c<B->n_states[0];c++)
            for(int d=0;d<B->n_states[1];d++){
                abcd = ((a*A->n_states[1]+b)*B->n_states[0]+c)*B->n_states[1]+d;
                ac   =   a*B->n_states[0]+c;
                bd   =   b*B->n_states[1]+d;

                for(int ii=0;ii<n_act[0];ii++)
                for(int jj=0;jj<n_act[0];jj++)
                    ACT_DMT_A[ii*n_act[0]+jj+abcd*n_act[0]*n_act[1]]+=s1_ext_vec[jj+ac*n_act[0]]*
                                                                      s2_ext_vec[ii+bd*n_act[0]]*sign;
                
                for(int ii=0;ii<n_act[0];ii++)H_2el[abcd] -=  v1_ext_vec[ii+ac*n_act[0]]*s2_ext_vec[ii+bd*n_act[0]]*p_SVD*sign;
                for(int ii=0;ii<n_act[0];ii++)H_2el[abcd] += v1b_ext_vec[ii+ac*n_act[0]]*s2_ext_vec[ii+bd*n_act[0]]*p_SVD*sign;
                for(int ii=0;ii<n_act[0];ii++)H_2el[abcd] -=  v2_ext_vec[ii+bd*n_act[0]]*s1_ext_vec[ii+ac*n_act[0]]*p_SVD*sign;
                for(int ii=0;ii<n_act[0];ii++)H_2el[abcd] += v2b_ext_vec[ii+bd*n_act[0]]*s1_ext_vec[ii+ac*n_act[0]]*p_SVD*sign;
            }
            
        }
//         if(false)//4
        if(ci1_lt_aldet[i].n_A== ci2_lt_aldet[j].n_A   )
        if(ci1_lt_aldet[i].n_B==(ci2_lt_aldet[j].n_B-1)){
            //fprintf(out_stream,"into 4\n");
            sign = pow(-1,ci1_lt_aldet[i].n_B);
            set_zero_matr(v1_ext_vec, n_act[0]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            set_zero_matr(v2_ext_vec, n_act[0]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            set_zero_matr(v1b_ext_vec,n_act[0]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            set_zero_matr(v2b_ext_vec,n_act[0]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            
            aldet_calc_CI_2el_A0AA(v1_ext_vec, &(ci1_lt_aldet[i].A),0, &(ci2_lt_aldet[j].A),0, act_INTS, 0, n_act_f[0][0], n_act[0],        A->n_states[0],B->n_states[0],1);
            aldet_calc_CI_2el_0BAB(v1b_ext_vec,&(ci1_lt_aldet[i].A),0, &(ci2_lt_aldet[j].A),0, act_INTS, 0, n_act_f[0][0], n_act[0],        A->n_states[0],B->n_states[0],1);
            aldet_calc_CI_2el_AAA0(v2_ext_vec, &(ci1_lt_aldet[i].B),0, &(ci2_lt_aldet[j].B),0, act_INTS, n_act_f[0][0], n_act[0], n_act[0], A->n_states[1],B->n_states[1],1);
            aldet_calc_CI_2el_AB0B(v2b_ext_vec,&(ci1_lt_aldet[i].B),0, &(ci2_lt_aldet[j].B),0, act_INTS, n_act_f[0][0], n_act[0], n_act[0], A->n_states[1],B->n_states[1],1);

            aldet_gen_1_el_vec_m_arr(s1_ext_vec,&(ci1_lt_aldet[i].A),0,A->n_states[0],&(ci2_lt_aldet[j].A),0,B->n_states[0], 0, n_act_f[0][0], n_act[0], 1);
            aldet_gen_1_el_vec_arr  (s2_ext_vec,&(ci1_lt_aldet[i].B),0,A->n_states[1],&(ci2_lt_aldet[j].B),0,B->n_states[1], n_act_f[0][0], n_act[0], n_act[0], 1);
            
            for(int a=0;a<A->n_states[0];a++)
            for(int b=0;b<A->n_states[1];b++)
            for(int c=0;c<B->n_states[0];c++)
            for(int d=0;d<B->n_states[1];d++){
                abcd = ((a*A->n_states[1]+b)*B->n_states[0]+c)*B->n_states[1]+d;
                ac   =   a*B->n_states[0]+c;
                bd   =   b*B->n_states[1]+d;

                for(int ii=0;ii<n_act[0];ii++)
                for(int jj=0;jj<n_act[0];jj++)
                    ACT_DMT_B[ii*n_act[0]+jj+abcd*n_act[0]*n_act[1]]+=s1_ext_vec[jj+ac*n_act[0]]*
                                                                      s2_ext_vec[ii+bd*n_act[0]]*sign;
                
                for(int ii=0;ii<n_act[0];ii++)H_2el[abcd] -=  v1_ext_vec[ii+ac*n_act[0]]*s2_ext_vec[ii+bd*n_act[0]]*p_SVD*sign;
                for(int ii=0;ii<n_act[0];ii++)H_2el[abcd] += v1b_ext_vec[ii+ac*n_act[0]]*s2_ext_vec[ii+bd*n_act[0]]*p_SVD*sign;
                for(int ii=0;ii<n_act[0];ii++)H_2el[abcd] -=  v2_ext_vec[ii+bd*n_act[0]]*s1_ext_vec[ii+ac*n_act[0]]*p_SVD*sign;
                for(int ii=0;ii<n_act[0];ii++)H_2el[abcd] += v2b_ext_vec[ii+bd*n_act[0]]*s1_ext_vec[ii+ac*n_act[0]]*p_SVD*sign;
            }
            
        }
//         if(false)//5
        if(ci1_lt_aldet[i].n_A==(ci2_lt_aldet[j].n_A+2))
        if(ci1_lt_aldet[i].n_B== ci2_lt_aldet[j].n_B  ){
//             fprintf(out_stream,"into 5\n");
            sign = pow(-1,2*ci1_lt_aldet[i].n_A-2);
            
            set_zero_matr(s1_ext_mat,n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            set_zero_matr(s2_ext_mat,n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            
            aldet_calc_CI_2el_00AA(s2_ext_mat, &(ci1_lt_aldet[i].B),0, &(ci2_lt_aldet[j].B),0, act_INTS, 0, n_act_f[0][1], n_act[0], A->n_states[1], B->n_states[1],0,0);
            aldet_calc_CI_2el_AA00(s1_ext_mat, &(ci1_lt_aldet[i].A),0, &(ci2_lt_aldet[j].A),0, act_INTS, 0, n_act_f[0][0], n_act[0], A->n_states[0], B->n_states[0],0,0);
            
            for(int a=0;a<A->n_states[0];a++)
            for(int b=0;b<A->n_states[1];b++)
            for(int c=0;c<B->n_states[0];c++)
            for(int d=0;d<B->n_states[1];d++){
                abcd = ((a*A->n_states[1]+b)*B->n_states[0]+c)*B->n_states[1]+d;
                ac   =   a*B->n_states[0]+c;
                bd   =   b*B->n_states[1]+d;

                for(int i=0;i<n_act_f[0][0];i++)
                for(int j=0;j<n_act_f[0][0];j++)
                for(int k=0;k<n_act_f[0][1];k++)
                for(int l=0;l<n_act_f[0][1];l++)
                    H_2el[abcd]+=(act_INTS[((i*n_act[0]+k+n_act_f[0][0])*n_act[0]+j)*n_act[0]+l+n_act_f[0][0]]-
                                  act_INTS[((i*n_act[0]+l+n_act_f[0][0])*n_act[0]+j)*n_act[0]+k+n_act_f[0][0]])
                                *s1_ext_mat[i*n_act[0]+j+ac*n_act[0]*n_act[1]]*
                                 s2_ext_mat[k*n_act[0]+l+bd*n_act[0]*n_act[1]]*p_SVD*sign;
            }
        }
//         if(false)//6
        if(ci1_lt_aldet[i].n_A== ci2_lt_aldet[j].n_A   )
        if(ci1_lt_aldet[i].n_B==(ci2_lt_aldet[j].n_B+2)){
//             fprintf(out_stream,"into 6\n");
            sign = pow(-1,2*ci1_lt_aldet[i].n_B-2);
            
            set_zero_matr(s1_ext_mat,n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            set_zero_matr(s2_ext_mat,n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            
            aldet_calc_CI_2el_00AA(s2_ext_mat, &(ci1_lt_aldet[i].B),0, &(ci2_lt_aldet[j].B),0, act_INTS, 0, n_act_f[0][1], n_act[0], A->n_states[1], B->n_states[1],1,1);
            aldet_calc_CI_2el_AA00(s1_ext_mat, &(ci1_lt_aldet[i].A),0, &(ci2_lt_aldet[j].A),0, act_INTS, 0, n_act_f[0][0], n_act[0], A->n_states[0], B->n_states[0],1,1);
            
            for(int a=0;a<A->n_states[0];a++)
            for(int b=0;b<A->n_states[1];b++)
            for(int c=0;c<B->n_states[0];c++)
            for(int d=0;d<B->n_states[1];d++){
                abcd = ((a*A->n_states[1]+b)*B->n_states[0]+c)*B->n_states[1]+d;
                ac   =   a*B->n_states[0]+c;
                bd   =   b*B->n_states[1]+d;

                for(int i=0;i<n_act_f[0][0];i++)
                for(int j=0;j<n_act_f[0][0];j++)
                for(int k=0;k<n_act_f[1][0];k++)
                for(int l=0;l<n_act_f[1][0];l++)
                    H_2el[abcd]+=(act_INTS[((i*n_act[0]+k+n_act_f[0][0])*n_act[0]+j)*n_act[0]+l+n_act_f[0][0]]-
                                  act_INTS[((i*n_act[0]+l+n_act_f[0][0])*n_act[0]+j)*n_act[0]+k+n_act_f[0][0]])
                                *s1_ext_mat[i*n_act[0]+j+ac*n_act[0]*n_act[1]]*
                                 s2_ext_mat[k*n_act[0]+l+bd*n_act[0]*n_act[1]]*p_SVD*sign;
            }
            
        }
//         if(false)//7
        if(ci1_lt_aldet[i].n_A==(ci2_lt_aldet[j].n_A-2))
        if(ci1_lt_aldet[i].n_B== ci2_lt_aldet[j].n_B  ){
//             fprintf(out_stream,"into 7\n");
            sign = pow(-1,2*ci1_lt_aldet[i].n_A);
            
            set_zero_matr(s1_ext_mat,n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            set_zero_matr(s2_ext_mat,n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            
            aldet_calc_CI_2el_00AA(s1_ext_mat, &(ci1_lt_aldet[i].A),0, &(ci2_lt_aldet[j].A),0, act_INTS, 0, n_act_f[0][0], n_act[0], A->n_states[0], B->n_states[0],0,0);
            aldet_calc_CI_2el_AA00(s2_ext_mat, &(ci1_lt_aldet[i].B),0, &(ci2_lt_aldet[j].B),0, act_INTS, 0, n_act_f[0][1], n_act[0], A->n_states[1], B->n_states[1],0,0);
            
            for(int a=0;a<A->n_states[0];a++)
            for(int b=0;b<A->n_states[1];b++)
            for(int c=0;c<B->n_states[0];c++)
            for(int d=0;d<B->n_states[1];d++){
                abcd = ((a*A->n_states[1]+b)*B->n_states[0]+c)*B->n_states[1]+d;
                ac   =   a*B->n_states[0]+c;
                bd   =   b*B->n_states[1]+d;

                for(int i=0;i<n_act_f[0][1];i++)
                for(int j=0;j<n_act_f[0][1];j++)
                for(int k=0;k<n_act_f[0][0];k++)
                for(int l=0;l<n_act_f[0][1];l++)
                    H_2el[abcd]+=(act_INTS[(((i+n_act_f[0][0])*n_act[0]+k)*n_act[0]+j+n_act_f[0][0])*n_act[0]+l]-
                                  act_INTS[(((i+n_act_f[0][0])*n_act[0]+l)*n_act[0]+j+n_act_f[0][0])*n_act[0]+k])
                                *s1_ext_mat[k*n_act[0]+l+ac*n_act[0]*n_act[1]]*
                                 s2_ext_mat[i*n_act[0]+j+bd*n_act[0]*n_act[1]]*p_SVD*sign;
            }
        }
//         if(false)//8
        if(ci1_lt_aldet[i].n_A== ci2_lt_aldet[j].n_A   )
        if(ci1_lt_aldet[i].n_B==(ci2_lt_aldet[j].n_B-2)){
//             fprintf(out_stream,"into 8\n");
            sign = pow(-1,2*ci1_lt_aldet[i].n_B);
            
            set_zero_matr(s1_ext_mat,n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            set_zero_matr(s2_ext_mat,n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            
            aldet_calc_CI_2el_00AA(s1_ext_mat, &(ci1_lt_aldet[i].A),0, &(ci2_lt_aldet[j].A),0, act_INTS, 0, n_act_f[0][0], n_act[0], A->n_states[0], B->n_states[0],1,1);
            aldet_calc_CI_2el_AA00(s2_ext_mat, &(ci1_lt_aldet[i].B),0, &(ci2_lt_aldet[j].B),0, act_INTS, 0, n_act_f[0][1], n_act[0], A->n_states[1], B->n_states[1],1,1);
            
            for(int a=0;a<A->n_states[0];a++)
            for(int b=0;b<A->n_states[1];b++)
            for(int c=0;c<B->n_states[0];c++)
            for(int d=0;d<B->n_states[1];d++){
                abcd = ((a*A->n_states[1]+b)*B->n_states[0]+c)*B->n_states[1]+d;
                ac   =   a*B->n_states[0]+c;
                bd   =   b*B->n_states[1]+d;

                for(int i=0;i<n_act_f[1][0];i++)
                for(int j=0;j<n_act_f[1][0];j++)
                for(int k=0;k<n_act_f[0][0];k++)
                for(int l=0;l<n_act_f[0][0];l++)
                    H_2el[abcd]+=(act_INTS[(((i+n_act_f[0][0])*n_act[0]+k)*n_act[0]+j+n_act_f[0][0])*n_act[0]+l]-
                                  act_INTS[(((i+n_act_f[0][0])*n_act[0]+l)*n_act[0]+j+n_act_f[0][0])*n_act[0]+k])
                                *s1_ext_mat[k*n_act[0]+l+ac*n_act[0]*n_act[1]]*
                                 s2_ext_mat[i*n_act[0]+j+bd*n_act[0]*n_act[1]]*p_SVD*sign;
            }
        }
//         if(false)//9
        if(ci1_lt_aldet[i].n_A==(ci2_lt_aldet[j].n_A+1))
        if(ci1_lt_aldet[i].n_B==(ci2_lt_aldet[j].n_B+1)){
//             fprintf(out_stream,"into 9\n");
            sign = pow(-1,ci1_lt_aldet[i].n_A+ci1_lt_aldet[i].n_B-2);
            
            set_zero_matr(s1_ext_mat,n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            set_zero_matr(s2_ext_mat,n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            
            aldet_calc_CI_2el_00AA(s2_ext_mat, &(ci1_lt_aldet[i].B),0, &(ci2_lt_aldet[j].B),0, act_INTS, 0, n_act_f[0][1], n_act[0], A->n_states[1], B->n_states[1],0,1);
            aldet_calc_CI_2el_AA00(s1_ext_mat, &(ci1_lt_aldet[i].A),0, &(ci2_lt_aldet[j].A),0, act_INTS, 0, n_act_f[0][0], n_act[0], A->n_states[0], B->n_states[0],0,1);
            
            for(int a=0;a<A->n_states[0];a++)
            for(int b=0;b<A->n_states[1];b++)
            for(int c=0;c<B->n_states[0];c++)
            for(int d=0;d<B->n_states[1];d++){
                abcd = ((a*A->n_states[1]+b)*B->n_states[0]+c)*B->n_states[1]+d;
                ac   =   a*B->n_states[0]+c;
                bd   =   b*B->n_states[1]+d;

                for(int i=0;i<n_act_f[0][0];i++)
                for(int j=0;j<n_act_f[0][0];j++)
                for(int k=0;k<n_act_f[0][1];k++)
                for(int l=0;l<n_act_f[0][1];l++)
                    H_2el[abcd]+=act_INTS[((i*n_act[0]+k+n_act_f[0][0])*n_act[0]+j)*n_act[0]+l+n_act_f[0][0]]
                                *s1_ext_mat[i*n_act[0]+j+ac*n_act[0]*n_act[1]]*
                                 s2_ext_mat[k*n_act[0]+l+bd*n_act[0]*n_act[1]]*p_SVD*sign;
            }
        }
//         if(false)//10
        if(ci1_lt_aldet[i].n_A==(ci2_lt_aldet[j].n_A-1))
        if(ci1_lt_aldet[i].n_B==(ci2_lt_aldet[j].n_B-1)){
//             fprintf(out_stream,"into 10\n");
            sign = pow(-1,ci1_lt_aldet[i].n_A+ci1_lt_aldet[i].n_B);
            
            set_zero_matr(s1_ext_mat,n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            set_zero_matr(s2_ext_mat,n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            
            aldet_calc_CI_2el_AA00(s2_ext_mat, &(ci1_lt_aldet[i].B),0, &(ci2_lt_aldet[j].B),0, act_INTS, 0, n_act_f[0][1], n_act[0], A->n_states[1], B->n_states[1],0,1);
            aldet_calc_CI_2el_00AA(s1_ext_mat, &(ci1_lt_aldet[i].A),0, &(ci2_lt_aldet[j].A),0, act_INTS, 0, n_act_f[0][0], n_act[0], A->n_states[0], B->n_states[0],0,1);
            
            for(int a=0;a<A->n_states[0];a++)
            for(int b=0;b<A->n_states[1];b++)
            for(int c=0;c<B->n_states[0];c++)
            for(int d=0;d<B->n_states[1];d++){
                abcd = ((a*A->n_states[1]+b)*B->n_states[0]+c)*B->n_states[1]+d;
                ac   =   a*B->n_states[0]+c;
                bd   =   b*B->n_states[1]+d;

                for(int i=0;i<n_act_f[0][1];i++)
                for(int j=0;j<n_act_f[0][1];j++)
                for(int k=0;k<n_act_f[0][0];k++)
                for(int l=0;l<n_act_f[0][0];l++)
                    H_2el[abcd]+=act_INTS[(((i+n_act_f[0][0])*n_act[0]+k)*n_act[0]+j+n_act_f[0][0])*n_act[0]+l]
                                *s1_ext_mat[k*n_act[0]+l+ac*n_act[0]*n_act[1]]
                                *s2_ext_mat[i*n_act[0]+j+bd*n_act[0]*n_act[1]]*p_SVD*sign;
            }
        }
//         if(false)//11
        if(ci1_lt_aldet[i].n_A==(ci2_lt_aldet[j].n_A+1))
        if(ci1_lt_aldet[i].n_B==(ci2_lt_aldet[j].n_B-1)){
//             fprintf(out_stream,"into 11\n");
            sign = pow(-1,ci1_lt_aldet[i].n_A+ci1_lt_aldet[i].n_B-1);
            
            set_zero_matr(s1_ext_mat,n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            set_zero_matr(s2_ext_mat,n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            
            aldet_calc_CI_2el_0A0B(s1_ext_mat, &(ci1_lt_aldet[i].A),0, &(ci2_lt_aldet[j].A),0, act_INTS, 0, n_act_f[0][0], n_act[0], A->n_states[0], B->n_states[0],0);
            aldet_calc_CI_2el_0A0B(s2_ext_mat, &(ci1_lt_aldet[i].B),0, &(ci2_lt_aldet[j].B),0, act_INTS, 0, n_act_f[0][1], n_act[0], A->n_states[1], B->n_states[1],1);
            
            for(int a=0;a<A->n_states[0];a++)
            for(int b=0;b<A->n_states[1];b++)
            for(int c=0;c<B->n_states[0];c++)
            for(int d=0;d<B->n_states[1];d++){
                abcd = ((a*A->n_states[1]+b)*B->n_states[0]+c)*B->n_states[1]+d;
                ac   =   a*B->n_states[0]+c;
                bd   =   b*B->n_states[1]+d;

                for(int i=0;i<n_act_f[0][1];i++)
                for(int j=0;j<n_act_f[0][0];j++)
                for(int k=0;k<n_act_f[0][0];k++)
                for(int l=0;l<n_act_f[0][1];l++)
                    H_2el[abcd]+=act_INTS[(((i+n_act_f[0][0])*n_act[0]+k)*n_act[0]+j)*n_act[0]+l+n_act_f[0][0]]
                                *s1_ext_mat[j*n_act[0]+k+ac*n_act[0]*n_act[1]]
                                *s2_ext_mat[i*n_act[0]+l+bd*n_act[0]*n_act[1]]*p_SVD*sign;
            }
        }
//         if(false)//12
        if(ci1_lt_aldet[i].n_A==(ci2_lt_aldet[j].n_A-1))
        if(ci1_lt_aldet[i].n_B==(ci2_lt_aldet[j].n_B+1)){
//             fprintf(out_stream,"into 12\n");
            sign = pow(-1,ci1_lt_aldet[i].n_A+ci1_lt_aldet[i].n_B-1);
            
            set_zero_matr(s1_ext_mat,n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            set_zero_matr(s2_ext_mat,n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            
            aldet_calc_CI_2el_0A0B(s1_ext_mat, &(ci1_lt_aldet[i].A),0, &(ci2_lt_aldet[j].A),0, act_INTS, 0, n_act_f[0][0], n_act[0], A->n_states[0], B->n_states[0],1);
            aldet_calc_CI_2el_0A0B(s2_ext_mat, &(ci1_lt_aldet[i].B),0, &(ci2_lt_aldet[j].B),0, act_INTS, 0, n_act_f[0][1], n_act[0], A->n_states[1], B->n_states[1],0);
            
            for(int a=0;a<A->n_states[0];a++)
            for(int b=0;b<A->n_states[1];b++)
            for(int c=0;c<B->n_states[0];c++)
            for(int d=0;d<B->n_states[1];d++){
                abcd = ((a*A->n_states[1]+b)*B->n_states[0]+c)*B->n_states[1]+d;
                ac   =   a*B->n_states[0]+c;
                bd   =   b*B->n_states[1]+d;

                for(int i=0;i<n_act_f[0][1];i++)
                for(int j=0;j<n_act_f[0][0];j++)
                for(int k=0;k<n_act_f[0][0];k++)
                for(int l=0;l<n_act_f[0][1];l++)
                    H_2el[abcd]+=act_INTS[(((i+n_act_f[0][0])*n_act[0]+k)*n_act[0]+j)*n_act[0]+l+n_act_f[0][0]]
                                *s1_ext_mat[j*n_act[0]+k+ac*n_act[0]*n_act[1]]
                                *s2_ext_mat[i*n_act[0]+l+bd*n_act[0]*n_act[1]]*p_SVD*sign;
            }
            
        }

        
    }
            // fprintf(stderr,"\r%d/%d",i , ci1_lt_aldet.size());
    }
    // fprintf(stderr,"\n");
    
    
    return 0;
}

int doCI_data::calc_CI_for_PT(){
    
    double * act_c_INTS_A=new double[n_act_f[0][0]*n_act_f[0][0]*n_act_f[0][0]*n_act_f[0][0]];
    double * act_c_INTS_B=new double[n_act_f[0][1]*n_act_f[0][1]*n_act_f[0][1]*n_act_f[0][1]];
    
    for(int i=0; i<n_act_f[0][0];i++)
    for(int j=0; j<n_act_f[0][0];j++)
    for(int k=0; k<n_act_f[0][0];k++)
    for(int l=0; l<n_act_f[0][0];l++)
        act_c_INTS_A[((i*n_act_f[0][0]+j)*n_act_f[0][0]+k)*n_act_f[0][0]+l]=
        act_INTS    [((i*n_act  [0]   +j)*n_act  [0]   +k)*n_act  [0]   +l];
    
    for(int i=0,ii=n_act_f[0][0]; i<n_act_f[0][1];i++,ii++)
    for(int j=0,jj=n_act_f[0][0]; j<n_act_f[0][1];j++,jj++)
    for(int k=0,kk=n_act_f[0][0]; k<n_act_f[0][1];k++,kk++)
    for(int l=0,ll=n_act_f[0][0]; l<n_act_f[0][1];l++,ll++)
        act_c_INTS_B[((i *n_act_f[0][0]+j )*n_act_f[0][0]+k )*n_act_f[0][0]+l ]=
        act_INTS    [((ii*n_act  [0]   +jj)*n_act  [0]   +kk)*n_act  [0]   +ll];
    
    double * act_INTS_act=act_INTS_1st_order;
    
    double * s1_ext_vec;
    double * s2_ext_vec;
    double * v1_ext_vec;
    double * v2_ext_vec;
    double * v1b_ext_vec;
    double * v2b_ext_vec;
    double * s1_ext_mat;
    double * s2_ext_mat;
    
    double sign = 1.0;
    
    ACT_DMT_A = new double[n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1])];
    ACT_DMT_B = new double[n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1])];
    s1_ext_vec= new double[n_act[0]         *(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1])];
    s2_ext_vec= new double[         n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1])];
    v1_ext_vec= new double[n_act[0]         *(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1])];
    v2_ext_vec= new double[         n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1])];
    v1b_ext_vec= new double[n_act[0]         *(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1])];
    v2b_ext_vec= new double[         n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1])];
    s1_ext_mat= new double[n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1])];
    s2_ext_mat= new double[n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1])];
    
    
    ACT_DMT_A_1st_order = new double[n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1])];
    ACT_DMT_B_1st_order = new double[n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1])];
    
    
    
    
    
    double * ACT_DET_A2= ACT_DET_A+  A->n_states[0]*B->n_states[0];
    
    
    set_zero_matr(h_2e     ,                  A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]);
    set_zero_matr(ACT_DET_A,                  A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]);
    set_zero_matr(ACT_ADJ_A,n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
    set_zero_matr(ACT_ADJ_B,n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
    
    set_zero_matr(ACT_DMT_A,n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1]));
    set_zero_matr(ACT_DMT_B,n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1]));
    set_zero_matr(s1_ext_vec,n_act[0]         *(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
    set_zero_matr(s2_ext_vec,         n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
    set_zero_matr(v1_ext_vec,n_act[0]         *(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
    set_zero_matr(v2_ext_vec,         n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
    set_zero_matr(v1b_ext_vec,n_act[0]         *(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
    set_zero_matr(v2b_ext_vec,         n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
    
    set_zero_matr(ACT_DMT_A_1st_order, n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1]));
    set_zero_matr(ACT_DMT_B_1st_order, n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1]));
    
    
    
    
    if(S_new==NULL)S_new= new double[A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1]];
    if(H_2el==NULL)H_2el= new double[A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1]];
    
    set_zero_matr(S_new,A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1]);
    set_zero_matr(H_2el,A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1]);
    
    
//     double H_2el_new = 0;
    
    
    double * TH_1A = ACT_ADJ_A+(/*i*n_st_total[1]+j */0                          )*n_act[0]*n_act[0];
    double * TH_1B = ACT_ADJ_B+(/*i*n_st_total[1]+j */0                          )*n_act[0]*n_act[0];
    double * TH_2A = ACT_ADJ_A+(/*i*n_st_total[1]+j*/+A->n_states[0]*B->n_states[0])*n_act[0]*n_act[1];
    double * TH_2B = ACT_ADJ_B+(/*i*n_st_total[1]+j*/+A->n_states[0]*B->n_states[0])*n_act[0]*n_act[1];
    
    int abcd;
    int ac;
    int bd;
    int mc_2 = n_act_f[0][0]*n_act[1]+n_act_f[0][0]; // beginning of 2-2 block
    
//     getchar();
    double * act_INTS_B=act_INTS+n_act[0]*n_act[0]*n_act[0]*n_act_f[0][1]+n_act[0]*n_act[0]*n_act_f[0][1]+n_act[0]*n_act_f[0][1]+n_act_f[0][1];
    for(int i=0;i<ci1_lt_aldet.size();i++){
    for(int j=0;j<ci2_lt_aldet.size();j++){

//         if(false)//0
        if(ci1_lt_aldet[i].n_A==ci2_lt_aldet[j].n_A)
        if(ci1_lt_aldet[i].n_B==ci2_lt_aldet[j].n_B)
        if((ci1_lt_aldet[i].order+ci2_lt_aldet[j].order)<2){
            
            fprintf(out_stream,"PT(%d+%d)\n",ci1_lt_aldet[i].order, ci2_lt_aldet[j].order);
            fprintf(out_stream,"n(0) = %d %d %d %d \n", ci1_lt_aldet[i].n_A, ci2_lt_aldet[j].n_A,ci1_lt_aldet[i].n_B,ci2_lt_aldet[j].n_B);
//             getchar();
            
            aldet_S_calc(ACT_DET_A ,&(ci1_lt_aldet[i].A), 0, &(ci2_lt_aldet[j].A), 0);
            aldet_S_calc(ACT_DET_A2 ,&(ci1_lt_aldet[i].B), 0, &(ci2_lt_aldet[j].B), 0);
            
            calc_2e_ints(h_2e, ci1_lt_aldet[i].A.n_act,ci1_lt_aldet[i].A.na,ci1_lt_aldet[i].A.nb, ci1_lt_aldet[i].A.coef[0], ci2_lt_aldet[j].A.coef[0], act_c_INTS_A,0,ci1_lt_aldet[i].A.n_states[0],ci1_lt_aldet[i].A.n_states[0],ci2_lt_aldet[j].A.n_states[0],ci2_lt_aldet[j].A.n_states[0]);
            
            calc_2e_ints(h_2e+A->n_states[0]*B->n_states[0], ci1_lt_aldet[i].B.n_act,ci1_lt_aldet[i].B.na,ci1_lt_aldet[i].B.nb, ci1_lt_aldet[i].B.coef[0], ci2_lt_aldet[j].B.coef[0], act_c_INTS_B,0,ci1_lt_aldet[i].B.n_states[0],ci1_lt_aldet[i].B.n_states[0],ci2_lt_aldet[j].B.n_states[0],ci2_lt_aldet[j].B.n_states[0]);
            
            aldet_DMA_calc(ACT_ADJ_A ,&(ci1_lt_aldet[i].A), 0, &(ci2_lt_aldet[j].A), 0, n_act[0]);
            aldet_DMB_calc(ACT_ADJ_B ,&(ci1_lt_aldet[i].A), 0, &(ci2_lt_aldet[j].A), 0, n_act[0]);

            aldet_DMA_calc(ACT_ADJ_A+A->n_states[0]*B->n_states[0]*n_act[0]*n_act[1] ,&(ci1_lt_aldet[i].B), 0, &(ci2_lt_aldet[j].B), 0, n_act[0]);
            aldet_DMB_calc(ACT_ADJ_B+A->n_states[0]*B->n_states[0]*n_act[0]*n_act[1] ,&(ci1_lt_aldet[i].B), 0, &(ci2_lt_aldet[j].B), 0, n_act[0]);
            
            for(int a=0;a<A->n_states[0];a++)
            for(int b=0;b<A->n_states[1];b++)
            for(int c=0;c<B->n_states[0];c++)
            for(int d=0;d<B->n_states[1];d++){
                abcd = ((a*A->n_states[1]+b)*B->n_states[0]+c)*B->n_states[1]+d;
                ac   =   a*B->n_states[0]+c;
                bd   =   b*B->n_states[1]+d;
                
                for(int ii=0;ii<n_act[0]*n_act[1];ii++)
                    ACT_DMT_A[ii+     abcd*n_act[0]*n_act[1]]+=TH_1A[ii+ac*n_act[0]*n_act[1]]*ACT_DET_A2[bd];
                                   
                for(int ii=0;ii<n_act[0]*n_act[1]-mc_2;ii++)
                    ACT_DMT_A[ii+mc_2+abcd*n_act[0]*n_act[1]]+=TH_2A[ii+bd*n_act[0]*n_act[1]]*ACT_DET_A[ac];
                
                for(int ii=0;ii<n_act[0]*n_act[1];ii++)
                    ACT_DMT_B[ii+     abcd*n_act[0]*n_act[1]]+=TH_1B[ii+ac*n_act[0]*n_act[1]]*ACT_DET_A2[bd];
                                   
                for(int ii=0;ii<n_act[0]*n_act[1]-mc_2;ii++)
                    ACT_DMT_B[ii+mc_2+abcd*n_act[0]*n_act[1]]+=TH_2B[ii+bd*n_act[0]*n_act[1]]*ACT_DET_A[ac];
                
                H_2el[abcd]+= V_AB_22_J_calc(act_INTS,TH_1A+ac*n_act[0]*n_act[1], TH_2A+bd*n_act[0]*n_act[1], 0, n_act_f[0][0], n_act_f[0][0], n_act_f[0][0]+n_act_f[0][1], n_act[0])*p_SVD;
                H_2el[abcd]+= V_AB_22_J_calc(act_INTS,TH_1A+ac*n_act[0]*n_act[1], TH_2B+bd*n_act[0]*n_act[1], 0, n_act_f[0][0], n_act_f[0][0], n_act_f[0][0]+n_act_f[0][1], n_act[0])*p_SVD;
                H_2el[abcd]+= V_AB_22_J_calc(act_INTS,TH_1B+ac*n_act[0]*n_act[1], TH_2A+bd*n_act[0]*n_act[1], 0, n_act_f[0][0], n_act_f[0][0], n_act_f[0][0]+n_act_f[0][1], n_act[0])*p_SVD;
                H_2el[abcd]+= V_AB_22_J_calc(act_INTS,TH_1B+ac*n_act[0]*n_act[1], TH_2B+bd*n_act[0]*n_act[1], 0, n_act_f[0][0], n_act_f[0][0], n_act_f[0][0]+n_act_f[0][1], n_act[0])*p_SVD;
//                 H_2el[abcd]-= V_AB_22_K_calc(act_INTS,TH_1A+ac*n_act[0]*n_act[1], TH_2A+bd*n_act[0]*n_act[1], 0, n_act_f[0][0], n_act_f[0][0], n_act_f[0][0]+n_act_f[0][1], n_act[0])*p_SVD;
//                 H_2el[abcd]-= V_AB_22_K_calc(act_INTS,TH_1B+ac*n_act[0]*n_act[1], TH_2B+bd*n_act[0]*n_act[1], 0, n_act_f[0][0], n_act_f[0][0], n_act_f[0][0]+n_act_f[0][1], n_act[0])*p_SVD;
                
                
                S_new  [abcd] += ACT_DET_A [ac]*ACT_DET_A2[bd];
                
                H_2el  [abcd] += ACT_DET_A [ac]*h_2e      [bd+A->n_states[0]*B->n_states[0]]*p_SVD;
                
                H_2el  [abcd] += h_2e      [ac]*ACT_DET_A2[bd]*p_SVD;
            }
        }
//         if(false)//1
        if(ci1_lt_aldet[i].n_A==(ci2_lt_aldet[j].n_A+1))
        if(ci1_lt_aldet[i].n_B== ci2_lt_aldet[j].n_B  )
        if((ci1_lt_aldet[i].order+ci2_lt_aldet[j].order)<1){
            fprintf(out_stream,"into 1\n");
            fprintf(out_stream,"PT(%d+%d)\n",ci1_lt_aldet[i].order, ci2_lt_aldet[j].order);
            fprintf(out_stream,"n = %d %d %d %d \n", ci1_lt_aldet[i].n_A, ci2_lt_aldet[j].n_A,ci1_lt_aldet[i].n_B,ci2_lt_aldet[j].n_B);
//             getchar();
//             act_INTS_act=act_INTS;
//             if((ci1_lt_aldet[i].order+ci2_lt_aldet[j].order)==0)
            
            sign = pow(-1,ci1_lt_aldet[i].n_A-1);

            set_zero_matr(v1_ext_vec, n_act[0]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            set_zero_matr(v2_ext_vec, n_act[0]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            set_zero_matr(v1b_ext_vec,n_act[0]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            set_zero_matr(v2b_ext_vec,n_act[0]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            
            aldet_calc_CI_2el_AAA0(v1_ext_vec , &(ci1_lt_aldet[i].A), 0, &(ci2_lt_aldet[j].A), 0, act_INTS_act, 0, n_act_f[0][0], n_act[0],        A->n_states[0],B->n_states[0],0);
            aldet_calc_CI_2el_AB0B(v1b_ext_vec, &(ci1_lt_aldet[i].A), 0, &(ci2_lt_aldet[j].A), 0, act_INTS_act, 0, n_act_f[0][0], n_act[0],        A->n_states[0],B->n_states[0],0);
            aldet_calc_CI_2el_A0AA(v2_ext_vec , &(ci1_lt_aldet[i].B), 0, &(ci2_lt_aldet[j].B), 0, act_INTS_act, n_act_f[0][0], n_act[0], n_act[0], A->n_states[1],B->n_states[1],0);
            aldet_calc_CI_2el_0BAB(v2b_ext_vec, &(ci1_lt_aldet[i].B), 0, &(ci2_lt_aldet[j].B), 0, act_INTS_act, n_act_f[0][0], n_act[0], n_act[0], A->n_states[1],B->n_states[1],0);
            
            aldet_gen_1_el_vec_arr  (s1_ext_vec,&(ci1_lt_aldet[i].A),0,A->n_states[0],&(ci2_lt_aldet[j].A),0,B->n_states[0], 0, n_act_f[0][0], n_act[0], 0);
            aldet_gen_1_el_vec_m_arr(s2_ext_vec,&(ci1_lt_aldet[i].B),0,A->n_states[1],&(ci2_lt_aldet[j].B),0,B->n_states[1], n_act_f[0][0], n_act[0], n_act[0], 0);
            
            for(int a=0;a<A->n_states[0];a++)
            for(int b=0;b<A->n_states[1];b++)
            for(int c=0;c<B->n_states[0];c++)
            for(int d=0;d<B->n_states[1];d++){
                abcd = ((a*A->n_states[1]+b)*B->n_states[0]+c)*B->n_states[1]+d;
                ac   =   a*B->n_states[0]+c;
                bd   =   b*B->n_states[1]+d;
            
               for(int ii=0;ii<n_act[0];ii++)
               for(int jj=0;jj<n_act[0];jj++)
                   ACT_DMT_A[ii*n_act[0]+jj+abcd*n_act[0]*n_act[1]]+=s1_ext_vec[ii+ac*n_act[0]]*
                                                                     s2_ext_vec[jj+bd*n_act[0]]*sign;
               
               for(int ii=0;ii<n_act[0];ii++)
               for(int jj=0;jj<n_act[0];jj++)
                   ACT_DMT_A_1st_order[ii*n_act[0]+jj+abcd*n_act[0]*n_act[1]]+=s1_ext_vec[ii+ac*n_act[0]]*
                                                                               s2_ext_vec[jj+bd*n_act[0]]*sign;
                                                                     
                                                                     
               for(int ii=0;ii<n_act[0];ii++)H_2el[abcd] -=  v1_ext_vec[ii+ac*n_act[0]]*s2_ext_vec[ii+bd*n_act[0]]*p_SVD*sign;
               for(int ii=0;ii<n_act[0];ii++)H_2el[abcd] += v1b_ext_vec[ii+ac*n_act[0]]*s2_ext_vec[ii+bd*n_act[0]]*p_SVD*sign;
               for(int ii=0;ii<n_act[0];ii++)H_2el[abcd] -=  v2_ext_vec[ii+bd*n_act[0]]*s1_ext_vec[ii+ac*n_act[0]]*p_SVD*sign;
               for(int ii=0;ii<n_act[0];ii++)H_2el[abcd] += v2b_ext_vec[ii+bd*n_act[0]]*s1_ext_vec[ii+ac*n_act[0]]*p_SVD*sign;
            }    
           
        }
//         if(false)//2
        if(ci1_lt_aldet[i].n_A== ci2_lt_aldet[j].n_A   )
        if(ci1_lt_aldet[i].n_B==(ci2_lt_aldet[j].n_B+1))
        if((ci1_lt_aldet[i].order+ci2_lt_aldet[j].order)<1){
            fprintf(out_stream,"into 2\n");
            fprintf(out_stream,"PT(%d+%d)\n",ci1_lt_aldet[i].order, ci2_lt_aldet[j].order);
            fprintf(out_stream,"n = %d %d %d %d \n", ci1_lt_aldet[i].n_A, ci2_lt_aldet[j].n_A,ci1_lt_aldet[i].n_B,ci2_lt_aldet[j].n_B);
//             getchar();
            
            
            sign = pow(-1,ci1_lt_aldet[i].n_B-1);
            set_zero_matr(v1_ext_vec, n_act[0]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            set_zero_matr(v2_ext_vec, n_act[0]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            set_zero_matr(v1b_ext_vec,n_act[0]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            set_zero_matr(v2b_ext_vec,n_act[0]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            
            aldet_calc_CI_2el_AAA0(v1_ext_vec, &(ci1_lt_aldet[i].A),0, &(ci2_lt_aldet[j].A),0, act_INTS_act, 0, n_act_f[0][0], n_act[0],        A->n_states[0],B->n_states[0],1);
            aldet_calc_CI_2el_AB0B(v1b_ext_vec,&(ci1_lt_aldet[i].A),0, &(ci2_lt_aldet[j].A),0, act_INTS_act, 0, n_act_f[0][0], n_act[0],        A->n_states[0],B->n_states[0],1);
            aldet_calc_CI_2el_A0AA(v2_ext_vec, &(ci1_lt_aldet[i].B),0, &(ci2_lt_aldet[j].B),0, act_INTS_act, n_act_f[0][0], n_act[0], n_act[0], A->n_states[1],B->n_states[1],1);
            aldet_calc_CI_2el_0BAB(v2b_ext_vec,&(ci1_lt_aldet[i].B),0, &(ci2_lt_aldet[j].B),0, act_INTS_act, n_act_f[0][0], n_act[0], n_act[0], A->n_states[1],B->n_states[1],1);
            
            aldet_gen_1_el_vec_arr  (s1_ext_vec,&(ci1_lt_aldet[i].A),0,A->n_states[0],&(ci2_lt_aldet[j].A),0,B->n_states[0], 0, n_act_f[0][0], n_act[0],1);
            aldet_gen_1_el_vec_m_arr(s2_ext_vec,&(ci1_lt_aldet[i].B),0,A->n_states[1],&(ci2_lt_aldet[j].B),0,B->n_states[1], n_act_f[0][0], n_act[0], n_act[0],1);
                        
            for(int a=0;a<A->n_states[0];a++)
            for(int b=0;b<A->n_states[1];b++)
            for(int c=0;c<B->n_states[0];c++)
            for(int d=0;d<B->n_states[1];d++){
                abcd = ((a*A->n_states[1]+b)*B->n_states[0]+c)*B->n_states[1]+d;
                ac   =   a*B->n_states[0]+c;
                bd   =   b*B->n_states[1]+d;

                for(int ii=0;ii<n_act[0];ii++)
                for(int jj=0;jj<n_act[0];jj++)
                    ACT_DMT_B[ii*n_act[0]+jj+abcd*n_act[0]*n_act[1]]+=s1_ext_vec[ii+ac*n_act[0]]*
                                                                      s2_ext_vec[jj+bd*n_act[0]]*sign;
                for(int ii=0;ii<n_act[0];ii++)
                for(int jj=0;jj<n_act[0];jj++)
                    ACT_DMT_B_1st_order[ii*n_act[0]+jj+abcd*n_act[0]*n_act[1]]+=s1_ext_vec[ii+ac*n_act[0]]*
                                                                                s2_ext_vec[jj+bd*n_act[0]]*sign;
               
                for(int ii=0;ii<n_act[0];ii++)H_2el[abcd] -=  v1_ext_vec[ii+ac*n_act[0]]*s2_ext_vec[ii+bd*n_act[0]]*p_SVD*sign;
                for(int ii=0;ii<n_act[0];ii++)H_2el[abcd] += v1b_ext_vec[ii+ac*n_act[0]]*s2_ext_vec[ii+bd*n_act[0]]*p_SVD*sign;
                for(int ii=0;ii<n_act[0];ii++)H_2el[abcd] -=  v2_ext_vec[ii+bd*n_act[0]]*s1_ext_vec[ii+ac*n_act[0]]*p_SVD*sign;
                for(int ii=0;ii<n_act[0];ii++)H_2el[abcd] += v2b_ext_vec[ii+bd*n_act[0]]*s1_ext_vec[ii+ac*n_act[0]]*p_SVD*sign;
            }
            
        }
//         if(false)//3
        if(ci1_lt_aldet[i].n_A==(ci2_lt_aldet[j].n_A-1))
        if(ci1_lt_aldet[i].n_B== ci2_lt_aldet[j].n_B  )
        if((ci1_lt_aldet[i].order+ci2_lt_aldet[j].order)<1){
            fprintf(out_stream,"into 3 %d %d\n",i,j);
            fprintf(out_stream,"PT(%d+%d)\n",ci1_lt_aldet[i].order, ci2_lt_aldet[j].order);
            fprintf(out_stream,"n = %d %d %d %d \n", ci1_lt_aldet[i].n_A, ci2_lt_aldet[j].n_A,ci1_lt_aldet[i].n_B,ci2_lt_aldet[j].n_B);
//             getchar();
            
            
            sign = pow(-1,ci1_lt_aldet[i].n_A);
            
            set_zero_matr(v1_ext_vec, n_act[0]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            set_zero_matr(v2_ext_vec, n_act[0]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            set_zero_matr(v1b_ext_vec,n_act[0]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            set_zero_matr(v2b_ext_vec,n_act[0]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            
            aldet_calc_CI_2el_A0AA(v1_ext_vec, &(ci1_lt_aldet[i].A),0, &(ci2_lt_aldet[j].A),0, act_INTS_act, 0, n_act_f[0][0], n_act[0],        A->n_states[0],B->n_states[0],0);
            aldet_calc_CI_2el_0BAB(v1b_ext_vec,&(ci1_lt_aldet[i].A),0, &(ci2_lt_aldet[j].A),0, act_INTS_act, 0, n_act_f[0][0], n_act[0],        A->n_states[0],B->n_states[0],0);
            aldet_calc_CI_2el_AAA0(v2_ext_vec, &(ci1_lt_aldet[i].B),0, &(ci2_lt_aldet[j].B),0, act_INTS_act, n_act_f[0][0], n_act[0], n_act[0], A->n_states[1],B->n_states[1],0);
            aldet_calc_CI_2el_AB0B(v2b_ext_vec,&(ci1_lt_aldet[i].B),0, &(ci2_lt_aldet[j].B),0, act_INTS_act, n_act_f[0][0], n_act[0], n_act[0], A->n_states[1],B->n_states[1],0);
            
            aldet_gen_1_el_vec_m_arr(s1_ext_vec,&(ci1_lt_aldet[i].A),0,A->n_states[0],&(ci2_lt_aldet[j].A),0,B->n_states[0], 0, n_act_f[0][0], n_act[0], 0);
            aldet_gen_1_el_vec_arr  (s2_ext_vec,&(ci1_lt_aldet[i].B),0,A->n_states[1],&(ci2_lt_aldet[j].B),0,B->n_states[1], n_act_f[0][0], n_act[0], n_act[0], 0);
            
            
            for(int a=0;a<A->n_states[0];a++)
            for(int b=0;b<A->n_states[1];b++)
            for(int c=0;c<B->n_states[0];c++)
            for(int d=0;d<B->n_states[1];d++){
                abcd = ((a*A->n_states[1]+b)*B->n_states[0]+c)*B->n_states[1]+d;
                ac   =   a*B->n_states[0]+c;
                bd   =   b*B->n_states[1]+d;

                for(int ii=0;ii<n_act[0];ii++)
                for(int jj=0;jj<n_act[0];jj++)
                    ACT_DMT_A[ii*n_act[0]+jj+abcd*n_act[0]*n_act[1]]+=s1_ext_vec[jj+ac*n_act[0]]*
                                                                      s2_ext_vec[ii+bd*n_act[0]]*sign;
                for(int ii=0;ii<n_act[0];ii++)
                for(int jj=0;jj<n_act[0];jj++)
                    ACT_DMT_A_1st_order[ii*n_act[0]+jj+abcd*n_act[0]*n_act[1]]+=s1_ext_vec[jj+ac*n_act[0]]*
                                                                                s2_ext_vec[ii+bd*n_act[0]]*sign;
                
                for(int ii=0;ii<n_act[0];ii++)H_2el[abcd] -=  v1_ext_vec[ii+ac*n_act[0]]*s2_ext_vec[ii+bd*n_act[0]]*p_SVD*sign;
                for(int ii=0;ii<n_act[0];ii++)H_2el[abcd] += v1b_ext_vec[ii+ac*n_act[0]]*s2_ext_vec[ii+bd*n_act[0]]*p_SVD*sign;
                for(int ii=0;ii<n_act[0];ii++)H_2el[abcd] -=  v2_ext_vec[ii+bd*n_act[0]]*s1_ext_vec[ii+ac*n_act[0]]*p_SVD*sign;
                for(int ii=0;ii<n_act[0];ii++)H_2el[abcd] += v2b_ext_vec[ii+bd*n_act[0]]*s1_ext_vec[ii+ac*n_act[0]]*p_SVD*sign;
            }
            
        }
//         if(false)//4
        if(ci1_lt_aldet[i].n_A== ci2_lt_aldet[j].n_A   )
        if(ci1_lt_aldet[i].n_B==(ci2_lt_aldet[j].n_B-1))
        if((ci1_lt_aldet[i].order+ci2_lt_aldet[j].order)<1){
            fprintf(out_stream,"into 4\n");
            fprintf(out_stream,"PT(%d+%d)\n",ci1_lt_aldet[i].order, ci2_lt_aldet[j].order);
            fprintf(out_stream,"n = %d %d %d %d \n", ci1_lt_aldet[i].n_A, ci2_lt_aldet[j].n_A,ci1_lt_aldet[i].n_B,ci2_lt_aldet[j].n_B);
//             getchar();
            
            
            sign = pow(-1,ci1_lt_aldet[i].n_B);
            set_zero_matr(v1_ext_vec, n_act[0]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            set_zero_matr(v2_ext_vec, n_act[0]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            set_zero_matr(v1b_ext_vec,n_act[0]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            set_zero_matr(v2b_ext_vec,n_act[0]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
            
            aldet_calc_CI_2el_A0AA(v1_ext_vec, &(ci1_lt_aldet[i].A),0, &(ci2_lt_aldet[j].A),0, act_INTS_act, 0, n_act_f[0][0], n_act[0],        A->n_states[0],B->n_states[0],1);
            aldet_calc_CI_2el_0BAB(v1b_ext_vec,&(ci1_lt_aldet[i].A),0, &(ci2_lt_aldet[j].A),0, act_INTS_act, 0, n_act_f[0][0], n_act[0],        A->n_states[0],B->n_states[0],1);
            aldet_calc_CI_2el_AAA0(v2_ext_vec, &(ci1_lt_aldet[i].B),0, &(ci2_lt_aldet[j].B),0, act_INTS_act, n_act_f[0][0], n_act[0], n_act[0], A->n_states[1],B->n_states[1],1);
            aldet_calc_CI_2el_AB0B(v2b_ext_vec,&(ci1_lt_aldet[i].B),0, &(ci2_lt_aldet[j].B),0, act_INTS_act, n_act_f[0][0], n_act[0], n_act[0], A->n_states[1],B->n_states[1],1);

            aldet_gen_1_el_vec_m_arr(s1_ext_vec,&(ci1_lt_aldet[i].A),0,A->n_states[0],&(ci2_lt_aldet[j].A),0,B->n_states[0], 0, n_act_f[0][0], n_act[0], 1);
            aldet_gen_1_el_vec_arr  (s2_ext_vec,&(ci1_lt_aldet[i].B),0,A->n_states[1],&(ci2_lt_aldet[j].B),0,B->n_states[1], n_act_f[0][0], n_act[0], n_act[0], 1);
            
            for(int a=0;a<A->n_states[0];a++)
            for(int b=0;b<A->n_states[1];b++)
            for(int c=0;c<B->n_states[0];c++)
            for(int d=0;d<B->n_states[1];d++){
                abcd = ((a*A->n_states[1]+b)*B->n_states[0]+c)*B->n_states[1]+d;
                ac   =   a*B->n_states[0]+c;
                bd   =   b*B->n_states[1]+d;

                for(int ii=0;ii<n_act[0];ii++)
                for(int jj=0;jj<n_act[0];jj++)
                    ACT_DMT_B[ii*n_act[0]+jj+abcd*n_act[0]*n_act[1]]+=s1_ext_vec[jj+ac*n_act[0]]*
                                                                      s2_ext_vec[ii+bd*n_act[0]]*sign;
                for(int ii=0;ii<n_act[0];ii++)
                for(int jj=0;jj<n_act[0];jj++)
                    ACT_DMT_B_1st_order[ii*n_act[0]+jj+abcd*n_act[0]*n_act[1]]+=s1_ext_vec[jj+ac*n_act[0]]*
                                                                                s2_ext_vec[ii+bd*n_act[0]]*sign;
                
                for(int ii=0;ii<n_act[0];ii++)H_2el[abcd] -=  v1_ext_vec[ii+ac*n_act[0]]*s2_ext_vec[ii+bd*n_act[0]]*p_SVD*sign;
                for(int ii=0;ii<n_act[0];ii++)H_2el[abcd] += v1b_ext_vec[ii+ac*n_act[0]]*s2_ext_vec[ii+bd*n_act[0]]*p_SVD*sign;
                for(int ii=0;ii<n_act[0];ii++)H_2el[abcd] -=  v2_ext_vec[ii+bd*n_act[0]]*s1_ext_vec[ii+ac*n_act[0]]*p_SVD*sign;
                for(int ii=0;ii<n_act[0];ii++)H_2el[abcd] += v2b_ext_vec[ii+bd*n_act[0]]*s1_ext_vec[ii+ac*n_act[0]]*p_SVD*sign;
            }
            
        }

        
    }
            // fprintf(stderr,"\r%d/%d",i , ci1_lt_aldet.size());
    }
    // fprintf(stderr,"\n");
    
    
    return 0;
}

int doCI_data::calc_DM_for_dec_wf(){
    
    
    double * s1_ext_vec;
    double * s2_ext_vec;
    
    double sign = 1.0;
    
    ACT_DMT_A = new double[n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1])];
    ACT_DMT_B = new double[n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1])];
    s1_ext_vec= new double[n_act[0]         *(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1])];
    s2_ext_vec= new double[         n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1])];
    
    double * ACT_DET_A2= ACT_DET_A+  A->n_states[0]*B->n_states[0];
    
    
    set_zero_matr(h_2e     ,                  A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]);
    set_zero_matr(ACT_DET_A,                  A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]);
    set_zero_matr(ACT_ADJ_A,n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
    set_zero_matr(ACT_ADJ_B,n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
    
    set_zero_matr(ACT_DMT_A,n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1]));
    set_zero_matr(ACT_DMT_B,n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1]));
    set_zero_matr(s1_ext_vec,n_act[0]         *(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
    set_zero_matr(s2_ext_vec,         n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
    
    
    if(S_new==NULL)S_new= new double[A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1]];
    if(H_2el==NULL)H_2el= new double[A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1]];
    
    set_zero_matr(S_new,A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1]);
    set_zero_matr(H_2el,A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1]);
    
    
//     double H_2el_new = 0;
    
    
    double * TH_1A = ACT_ADJ_A+(/*i*n_st_total[1]+j */0                          )*n_act[0]*n_act[0];
    double * TH_1B = ACT_ADJ_B+(/*i*n_st_total[1]+j */0                          )*n_act[0]*n_act[0];
    double * TH_2A = ACT_ADJ_A+(/*i*n_st_total[1]+j*/+A->n_states[0]*B->n_states[0])*n_act[0]*n_act[1];
    double * TH_2B = ACT_ADJ_B+(/*i*n_st_total[1]+j*/+A->n_states[0]*B->n_states[0])*n_act[0]*n_act[1];
    
    int abcd;
    int ac;
    int bd;
    int mc_2 = n_act_f[0][0]*n_act[1]+n_act_f[0][0]; // beginning of 2-2 block
    
    for(int i=0;i<ci1_lt_aldet.size();i++){
    for(int j=0;j<ci2_lt_aldet.size();j++){

//         if(false)//0
        if(ci1_lt_aldet[i].n_A==ci2_lt_aldet[j].n_A)
        if(ci1_lt_aldet[i].n_B==ci2_lt_aldet[j].n_B){
            aldet_S_calc(ACT_DET_A ,&(ci1_lt_aldet[i].A), 0, &(ci2_lt_aldet[j].A), 0);
            aldet_S_calc(ACT_DET_A2 ,&(ci1_lt_aldet[i].B), 0, &(ci2_lt_aldet[j].B), 0);
            
            aldet_DMA_calc(ACT_ADJ_A ,&(ci1_lt_aldet[i].A), 0, &(ci2_lt_aldet[j].A), 0, n_act[0]);
            aldet_DMB_calc(ACT_ADJ_B ,&(ci1_lt_aldet[i].A), 0, &(ci2_lt_aldet[j].A), 0, n_act[0]);

            aldet_DMA_calc(ACT_ADJ_A+A->n_states[0]*B->n_states[0]*n_act[0]*n_act[1] ,&(ci1_lt_aldet[i].B), 0, &(ci2_lt_aldet[j].B), 0, n_act[0]);
            aldet_DMB_calc(ACT_ADJ_B+A->n_states[0]*B->n_states[0]*n_act[0]*n_act[1] ,&(ci1_lt_aldet[i].B), 0, &(ci2_lt_aldet[j].B), 0, n_act[0]);
            
            for(int a=0;a<A->n_states[0];a++)
            for(int b=0;b<A->n_states[1];b++)
            for(int c=0;c<B->n_states[0];c++)
            for(int d=0;d<B->n_states[1];d++){
                abcd = ((a*A->n_states[1]+b)*B->n_states[0]+c)*B->n_states[1]+d;
                ac   =   a*B->n_states[0]+c;
                bd   =   b*B->n_states[1]+d;
                
                for(int ii=0;ii<n_act[0]*n_act[1];ii++)
                    ACT_DMT_A[ii+     abcd*n_act[0]*n_act[1]]+=TH_1A[ii+ac*n_act[0]*n_act[1]]*ACT_DET_A2[bd];
                                   
                for(int ii=0;ii<n_act[0]*n_act[1]-mc_2;ii++)
                    ACT_DMT_A[ii+mc_2+abcd*n_act[0]*n_act[1]]+=TH_2A[ii+bd*n_act[0]*n_act[1]]*ACT_DET_A[ac];
                
                for(int ii=0;ii<n_act[0]*n_act[1];ii++)
                    ACT_DMT_B[ii+     abcd*n_act[0]*n_act[1]]+=TH_1B[ii+ac*n_act[0]*n_act[1]]*ACT_DET_A2[bd];
                                   
                for(int ii=0;ii<n_act[0]*n_act[1]-mc_2;ii++)
                    ACT_DMT_B[ii+mc_2+abcd*n_act[0]*n_act[1]]+=TH_2B[ii+bd*n_act[0]*n_act[1]]*ACT_DET_A[ac];
                
                
                S_new  [abcd] += ACT_DET_A [ac]*ACT_DET_A2[bd];
                
            }
        }
//         if(false)//1
        if(ci1_lt_aldet[i].n_A==(ci2_lt_aldet[j].n_A+1))
        if(ci1_lt_aldet[i].n_B== ci2_lt_aldet[j].n_B  ){
//             //fprintf(out_stream,"into 1\n");
            
            sign = pow(-1,ci1_lt_aldet[i].n_A-1);
            
            aldet_gen_1_el_vec_arr  (s1_ext_vec,&(ci1_lt_aldet[i].A),0,A->n_states[0],&(ci2_lt_aldet[j].A),0,B->n_states[0], 0, n_act_f[0][0], n_act[0], 0);
            aldet_gen_1_el_vec_m_arr(s2_ext_vec,&(ci1_lt_aldet[i].B),0,A->n_states[1],&(ci2_lt_aldet[j].B),0,B->n_states[1], n_act_f[0][0], n_act[0], n_act[0], 0);
            
            for(int a=0;a<A->n_states[0];a++)
            for(int b=0;b<A->n_states[1];b++)
            for(int c=0;c<B->n_states[0];c++)
            for(int d=0;d<B->n_states[1];d++){
                abcd = ((a*A->n_states[1]+b)*B->n_states[0]+c)*B->n_states[1]+d;
                ac   =   a*B->n_states[0]+c;
                bd   =   b*B->n_states[1]+d;
            
               for(int ii=0;ii<n_act[0];ii++)
               for(int jj=0;jj<n_act[0];jj++)
                   ACT_DMT_A[ii*n_act[0]+jj+abcd*n_act[0]*n_act[1]]+=s1_ext_vec[ii+ac*n_act[0]]*
                                                                     s2_ext_vec[jj+bd*n_act[0]]*sign;
               
            }    
           
        }
//         if(false)//2
        if(ci1_lt_aldet[i].n_A== ci2_lt_aldet[j].n_A   )
        if(ci1_lt_aldet[i].n_B==(ci2_lt_aldet[j].n_B+1)){
            //fprintf(out_stream,"into 2\n");
            sign = pow(-1,ci1_lt_aldet[i].n_B-1);
            
            aldet_gen_1_el_vec_arr  (s1_ext_vec,&(ci1_lt_aldet[i].A),0,A->n_states[0],&(ci2_lt_aldet[j].A),0,B->n_states[0], 0, n_act_f[0][0], n_act[0],1);
            aldet_gen_1_el_vec_m_arr(s2_ext_vec,&(ci1_lt_aldet[i].B),0,A->n_states[1],&(ci2_lt_aldet[j].B),0,B->n_states[1], n_act_f[0][0], n_act[0], n_act[0],1);
                        
            for(int a=0;a<A->n_states[0];a++)
            for(int b=0;b<A->n_states[1];b++)
            for(int c=0;c<B->n_states[0];c++)
            for(int d=0;d<B->n_states[1];d++){
                abcd = ((a*A->n_states[1]+b)*B->n_states[0]+c)*B->n_states[1]+d;
                ac   =   a*B->n_states[0]+c;
                bd   =   b*B->n_states[1]+d;

                for(int ii=0;ii<n_act[0];ii++)
                for(int jj=0;jj<n_act[0];jj++)
                    ACT_DMT_B[ii*n_act[0]+jj+abcd*n_act[0]*n_act[1]]+=s1_ext_vec[ii+ac*n_act[0]]*
                                                                      s2_ext_vec[jj+bd*n_act[0]]*sign;
            }
            
        }
//         if(false)//3
        if(ci1_lt_aldet[i].n_A==(ci2_lt_aldet[j].n_A-1))
        if(ci1_lt_aldet[i].n_B== ci2_lt_aldet[j].n_B  ){
//             fprintf(out_stream,"into 3 %d %d\n",i,j);
            sign = pow(-1,ci1_lt_aldet[i].n_A);
            
            aldet_gen_1_el_vec_m_arr(s1_ext_vec,&(ci1_lt_aldet[i].A),0,A->n_states[0],&(ci2_lt_aldet[j].A),0,B->n_states[0], 0, n_act_f[0][0], n_act[0], 0);
            aldet_gen_1_el_vec_arr  (s2_ext_vec,&(ci1_lt_aldet[i].B),0,A->n_states[1],&(ci2_lt_aldet[j].B),0,B->n_states[1], n_act_f[0][0], n_act[0], n_act[0], 0);
            
            
            for(int a=0;a<A->n_states[0];a++)
            for(int b=0;b<A->n_states[1];b++)
            for(int c=0;c<B->n_states[0];c++)
            for(int d=0;d<B->n_states[1];d++){
                abcd = ((a*A->n_states[1]+b)*B->n_states[0]+c)*B->n_states[1]+d;
                ac   =   a*B->n_states[0]+c;
                bd   =   b*B->n_states[1]+d;

                for(int ii=0;ii<n_act[0];ii++)
                for(int jj=0;jj<n_act[0];jj++)
                    ACT_DMT_A[ii*n_act[0]+jj+abcd*n_act[0]*n_act[1]]+=s1_ext_vec[jj+ac*n_act[0]]*
                                                                      s2_ext_vec[ii+bd*n_act[0]]*sign;
                
            }
            
        }
//         if(false)//4
        if(ci1_lt_aldet[i].n_A== ci2_lt_aldet[j].n_A   )
        if(ci1_lt_aldet[i].n_B==(ci2_lt_aldet[j].n_B-1)){
            //fprintf(out_stream,"into 4\n");
            sign = pow(-1,ci1_lt_aldet[i].n_B);
            
            aldet_gen_1_el_vec_m_arr(s1_ext_vec,&(ci1_lt_aldet[i].A),0,A->n_states[0],&(ci2_lt_aldet[j].A),0,B->n_states[0], 0, n_act_f[0][0], n_act[0], 1);
            aldet_gen_1_el_vec_arr  (s2_ext_vec,&(ci1_lt_aldet[i].B),0,A->n_states[1],&(ci2_lt_aldet[j].B),0,B->n_states[1], n_act_f[0][0], n_act[0], n_act[0], 1);
            
            for(int a=0;a<A->n_states[0];a++)
            for(int b=0;b<A->n_states[1];b++)
            for(int c=0;c<B->n_states[0];c++)
            for(int d=0;d<B->n_states[1];d++){
                abcd = ((a*A->n_states[1]+b)*B->n_states[0]+c)*B->n_states[1]+d;
                ac   =   a*B->n_states[0]+c;
                bd   =   b*B->n_states[1]+d;

                for(int ii=0;ii<n_act[0];ii++)
                for(int jj=0;jj<n_act[0];jj++)
                    ACT_DMT_B[ii*n_act[0]+jj+abcd*n_act[0]*n_act[1]]+=s1_ext_vec[jj+ac*n_act[0]]*
                                                                      s2_ext_vec[ii+bd*n_act[0]]*sign;
                
            }
            
        }

    }
            // fprintf(stderr,"\r%d/%d",i , ci1_lt_aldet.size());
    }
    // fprintf(stderr,"\n");
    
    
    delete[] s1_ext_vec;
    delete[] s2_ext_vec;
    
    return 0;
}

int doCI_data::calc_DM_for_PT(){
    
    
    double * s1_ext_vec;
    double * s2_ext_vec;
    
    double sign = 1.0;
    
    ACT_DMT_A = new double[n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1])];
    ACT_DMT_B = new double[n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1])];
    s1_ext_vec= new double[n_act[0]         *(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1])];
    s2_ext_vec= new double[         n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1])];
    
    ACT_DMT_A_1st_order = new double[n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1])];
    ACT_DMT_B_1st_order = new double[n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1])];

    
    double * ACT_DET_A2= ACT_DET_A+  A->n_states[0]*B->n_states[0];
    
    
    set_zero_matr(h_2e     ,                  A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]);
    set_zero_matr(ACT_DET_A,                  A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]);
    set_zero_matr(ACT_ADJ_A,n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
    set_zero_matr(ACT_ADJ_B,n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
    
    set_zero_matr(ACT_DMT_A,n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1]));
    set_zero_matr(ACT_DMT_B,n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1]));
    set_zero_matr(s1_ext_vec,n_act[0]         *(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
    set_zero_matr(s2_ext_vec,         n_act[1]*(A->n_states[0]*B->n_states[0]+A->n_states[1]*B->n_states[1]));
    
    set_zero_matr(ACT_DMT_A_1st_order, n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1]));
    set_zero_matr(ACT_DMT_B_1st_order, n_act[0]*n_act[1]*(A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1]));
    
    if(S_new==NULL)S_new= new double[A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1]];
    if(H_2el==NULL)H_2el= new double[A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1]];
    
    set_zero_matr(S_new,A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1]);
    set_zero_matr(H_2el,A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1]);
    
    
//     double H_2el_new = 0;
    
    
    double * TH_1A = ACT_ADJ_A+(/*i*n_st_total[1]+j */0                          )*n_act[0]*n_act[0];
    double * TH_1B = ACT_ADJ_B+(/*i*n_st_total[1]+j */0                          )*n_act[0]*n_act[0];
    double * TH_2A = ACT_ADJ_A+(/*i*n_st_total[1]+j*/+A->n_states[0]*B->n_states[0])*n_act[0]*n_act[1];
    double * TH_2B = ACT_ADJ_B+(/*i*n_st_total[1]+j*/+A->n_states[0]*B->n_states[0])*n_act[0]*n_act[1];
    
    int abcd;
    int ac;
    int bd;
    int mc_2 = n_act_f[0][0]*n_act[1]+n_act_f[0][0]; // beginning of 2-2 block
    
    for(int i=0;i<ci1_lt_aldet.size();i++){
    for(int j=0;j<ci2_lt_aldet.size();j++){

//         if(false)//0
        if(ci1_lt_aldet[i].n_A==ci2_lt_aldet[j].n_A)
        if(ci1_lt_aldet[i].n_B==ci2_lt_aldet[j].n_B)
        if((ci1_lt_aldet[i].order+ci2_lt_aldet[j].order)<2){
            aldet_S_calc(ACT_DET_A ,&(ci1_lt_aldet[i].A), 0, &(ci2_lt_aldet[j].A), 0);
            aldet_S_calc(ACT_DET_A2 ,&(ci1_lt_aldet[i].B), 0, &(ci2_lt_aldet[j].B), 0);
            
            aldet_DMA_calc(ACT_ADJ_A ,&(ci1_lt_aldet[i].A), 0, &(ci2_lt_aldet[j].A), 0, n_act[0]);
            aldet_DMB_calc(ACT_ADJ_B ,&(ci1_lt_aldet[i].A), 0, &(ci2_lt_aldet[j].A), 0, n_act[0]);

            aldet_DMA_calc(ACT_ADJ_A+A->n_states[0]*B->n_states[0]*n_act[0]*n_act[1] ,&(ci1_lt_aldet[i].B), 0, &(ci2_lt_aldet[j].B), 0, n_act[0]);
            aldet_DMB_calc(ACT_ADJ_B+A->n_states[0]*B->n_states[0]*n_act[0]*n_act[1] ,&(ci1_lt_aldet[i].B), 0, &(ci2_lt_aldet[j].B), 0, n_act[0]);
            
            for(int a=0;a<A->n_states[0];a++)
            for(int b=0;b<A->n_states[1];b++)
            for(int c=0;c<B->n_states[0];c++)
            for(int d=0;d<B->n_states[1];d++){
                abcd = ((a*A->n_states[1]+b)*B->n_states[0]+c)*B->n_states[1]+d;
                ac   =   a*B->n_states[0]+c;
                bd   =   b*B->n_states[1]+d;
                
                for(int ii=0;ii<n_act[0]*n_act[1];ii++)
                    ACT_DMT_A[ii+     abcd*n_act[0]*n_act[1]]+=TH_1A[ii+ac*n_act[0]*n_act[1]]*ACT_DET_A2[bd];
                                   
                for(int ii=0;ii<n_act[0]*n_act[1]-mc_2;ii++)
                    ACT_DMT_A[ii+mc_2+abcd*n_act[0]*n_act[1]]+=TH_2A[ii+bd*n_act[0]*n_act[1]]*ACT_DET_A[ac];
                
                for(int ii=0;ii<n_act[0]*n_act[1];ii++)
                    ACT_DMT_B[ii+     abcd*n_act[0]*n_act[1]]+=TH_1B[ii+ac*n_act[0]*n_act[1]]*ACT_DET_A2[bd];
                                   
                for(int ii=0;ii<n_act[0]*n_act[1]-mc_2;ii++)
                    ACT_DMT_B[ii+mc_2+abcd*n_act[0]*n_act[1]]+=TH_2B[ii+bd*n_act[0]*n_act[1]]*ACT_DET_A[ac];
                
                
                S_new  [abcd] += ACT_DET_A [ac]*ACT_DET_A2[bd];
                
            }
        }
//         if(false)//1
        if(ci1_lt_aldet[i].n_A==(ci2_lt_aldet[j].n_A+1))
        if(ci1_lt_aldet[i].n_B== ci2_lt_aldet[j].n_B  )
        if((ci1_lt_aldet[i].order+ci2_lt_aldet[j].order)<1){
//             //fprintf(out_stream,"into 1\n");
            
            sign = pow(-1,ci1_lt_aldet[i].n_A-1);
            
            aldet_gen_1_el_vec_arr  (s1_ext_vec,&(ci1_lt_aldet[i].A),0,A->n_states[0],&(ci2_lt_aldet[j].A),0,B->n_states[0], 0, n_act_f[0][0], n_act[0], 0);
            aldet_gen_1_el_vec_m_arr(s2_ext_vec,&(ci1_lt_aldet[i].B),0,A->n_states[1],&(ci2_lt_aldet[j].B),0,B->n_states[1], n_act_f[0][0], n_act[0], n_act[0], 0);
            
            for(int a=0;a<A->n_states[0];a++)
            for(int b=0;b<A->n_states[1];b++)
            for(int c=0;c<B->n_states[0];c++)
            for(int d=0;d<B->n_states[1];d++){
                abcd = ((a*A->n_states[1]+b)*B->n_states[0]+c)*B->n_states[1]+d;
                ac   =   a*B->n_states[0]+c;
                bd   =   b*B->n_states[1]+d;
            
               for(int ii=0;ii<n_act[0];ii++)
               for(int jj=0;jj<n_act[0];jj++)
                   ACT_DMT_A[ii*n_act[0]+jj+abcd*n_act[0]*n_act[1]]+=s1_ext_vec[ii+ac*n_act[0]]*
                                                                     s2_ext_vec[jj+bd*n_act[0]]*sign;
               for(int ii=0;ii<n_act[0];ii++)
               for(int jj=0;jj<n_act[0];jj++)
                   ACT_DMT_A_1st_order[ii*n_act[0]+jj+abcd*n_act[0]*n_act[1]]+=s1_ext_vec[ii+ac*n_act[0]]*
                                                                               s2_ext_vec[jj+bd*n_act[0]]*sign;
               
            }    
           
        }
//         if(false)//2
        if(ci1_lt_aldet[i].n_A== ci2_lt_aldet[j].n_A   )
        if(ci1_lt_aldet[i].n_B==(ci2_lt_aldet[j].n_B+1))
        if((ci1_lt_aldet[i].order+ci2_lt_aldet[j].order)<1){
            //fprintf(out_stream,"into 2\n");
            sign = pow(-1,ci1_lt_aldet[i].n_B-1);
            
            aldet_gen_1_el_vec_arr  (s1_ext_vec,&(ci1_lt_aldet[i].A),0,A->n_states[0],&(ci2_lt_aldet[j].A),0,B->n_states[0], 0, n_act_f[0][0], n_act[0],1);
            aldet_gen_1_el_vec_m_arr(s2_ext_vec,&(ci1_lt_aldet[i].B),0,A->n_states[1],&(ci2_lt_aldet[j].B),0,B->n_states[1], n_act_f[0][0], n_act[0], n_act[0],1);
                        
            for(int a=0;a<A->n_states[0];a++)
            for(int b=0;b<A->n_states[1];b++)
            for(int c=0;c<B->n_states[0];c++)
            for(int d=0;d<B->n_states[1];d++){
                abcd = ((a*A->n_states[1]+b)*B->n_states[0]+c)*B->n_states[1]+d;
                ac   =   a*B->n_states[0]+c;
                bd   =   b*B->n_states[1]+d;

                for(int ii=0;ii<n_act[0];ii++)
                for(int jj=0;jj<n_act[0];jj++)
                    ACT_DMT_B[ii*n_act[0]+jj+abcd*n_act[0]*n_act[1]]+=s1_ext_vec[ii+ac*n_act[0]]*
                                                                      s2_ext_vec[jj+bd*n_act[0]]*sign;
                for(int ii=0;ii<n_act[0];ii++)
                for(int jj=0;jj<n_act[0];jj++)
                    ACT_DMT_B_1st_order[ii*n_act[0]+jj+abcd*n_act[0]*n_act[1]]+=s1_ext_vec[ii+ac*n_act[0]]*
                                                                                s2_ext_vec[jj+bd*n_act[0]]*sign;
            }
            
        }
//         if(false)//3
        if(ci1_lt_aldet[i].n_A==(ci2_lt_aldet[j].n_A-1))
        if(ci1_lt_aldet[i].n_B== ci2_lt_aldet[j].n_B  )
        if((ci1_lt_aldet[i].order+ci2_lt_aldet[j].order)<1){
//             fprintf(out_stream,"into 3 %d %d\n",i,j);
            sign = pow(-1,ci1_lt_aldet[i].n_A);
            
            aldet_gen_1_el_vec_m_arr(s1_ext_vec,&(ci1_lt_aldet[i].A),0,A->n_states[0],&(ci2_lt_aldet[j].A),0,B->n_states[0], 0, n_act_f[0][0], n_act[0], 0);
            aldet_gen_1_el_vec_arr  (s2_ext_vec,&(ci1_lt_aldet[i].B),0,A->n_states[1],&(ci2_lt_aldet[j].B),0,B->n_states[1], n_act_f[0][0], n_act[0], n_act[0], 0);
            
            
            for(int a=0;a<A->n_states[0];a++)
            for(int b=0;b<A->n_states[1];b++)
            for(int c=0;c<B->n_states[0];c++)
            for(int d=0;d<B->n_states[1];d++){
                abcd = ((a*A->n_states[1]+b)*B->n_states[0]+c)*B->n_states[1]+d;
                ac   =   a*B->n_states[0]+c;
                bd   =   b*B->n_states[1]+d;

                for(int ii=0;ii<n_act[0];ii++)
                for(int jj=0;jj<n_act[0];jj++)
                    ACT_DMT_A[ii*n_act[0]+jj+abcd*n_act[0]*n_act[1]]+=s1_ext_vec[jj+ac*n_act[0]]*
                                                                      s2_ext_vec[ii+bd*n_act[0]]*sign;
                for(int ii=0;ii<n_act[0];ii++)
                for(int jj=0;jj<n_act[0];jj++)
                    ACT_DMT_A_1st_order[ii*n_act[0]+jj+abcd*n_act[0]*n_act[1]]+=s1_ext_vec[jj+ac*n_act[0]]*
                                                                                s2_ext_vec[ii+bd*n_act[0]]*sign;
                
            }
            
        }
//         if(false)//4
        if(ci1_lt_aldet[i].n_A== ci2_lt_aldet[j].n_A   )
        if(ci1_lt_aldet[i].n_B==(ci2_lt_aldet[j].n_B-1)){
            //fprintf(out_stream,"into 4\n");
            sign = pow(-1,ci1_lt_aldet[i].n_B);
            
            aldet_gen_1_el_vec_m_arr(s1_ext_vec,&(ci1_lt_aldet[i].A),0,A->n_states[0],&(ci2_lt_aldet[j].A),0,B->n_states[0], 0, n_act_f[0][0], n_act[0], 1);
            aldet_gen_1_el_vec_arr  (s2_ext_vec,&(ci1_lt_aldet[i].B),0,A->n_states[1],&(ci2_lt_aldet[j].B),0,B->n_states[1], n_act_f[0][0], n_act[0], n_act[0], 1);
            
            for(int a=0;a<A->n_states[0];a++)
            for(int b=0;b<A->n_states[1];b++)
            for(int c=0;c<B->n_states[0];c++)
            for(int d=0;d<B->n_states[1];d++){
                abcd = ((a*A->n_states[1]+b)*B->n_states[0]+c)*B->n_states[1]+d;
                ac   =   a*B->n_states[0]+c;
                bd   =   b*B->n_states[1]+d;

                for(int ii=0;ii<n_act[0];ii++)
                for(int jj=0;jj<n_act[0];jj++)
                    ACT_DMT_B[ii*n_act[0]+jj+abcd*n_act[0]*n_act[1]]+=s1_ext_vec[jj+ac*n_act[0]]*
                                                                      s2_ext_vec[ii+bd*n_act[0]]*sign;
                for(int ii=0;ii<n_act[0];ii++)
                for(int jj=0;jj<n_act[0];jj++)
                    ACT_DMT_B_1st_order[ii*n_act[0]+jj+abcd*n_act[0]*n_act[1]]+=s1_ext_vec[jj+ac*n_act[0]]*
                                                                                s2_ext_vec[ii+bd*n_act[0]]*sign;
                
            }
            
        }

    }
            // fprintf(stderr,"\r%d/%d",i , ci1_lt_aldet.size());
    }
    // fprintf(stderr,"\n");
    
    
    delete[] s1_ext_vec;
    delete[] s2_ext_vec;
    
    return 0;
}


int doCI_data::calc_S_for_dec_wf(){
        
    double sign = 1.0;
    
    
    double * ACT_DET_A2= ACT_DET_A+  A->n_states[0]*B->n_states[0];
    
    
    if(S_new==NULL)S_new= new double[A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1]];
    
    set_zero_matr(S_new,A->n_states[0]*B->n_states[0]*A->n_states[1]*B->n_states[1]);
    
    
//     double H_2el_new = 0;
    
    
    int abcd;
    int ac;
    int bd;
    int mc_2 = n_act_f[0][0]*n_act[1]+n_act_f[0][0]; // beginning of 2-2 block
    
//     getchar();
    double * act_INTS_B=act_INTS+n_act[0]*n_act[0]*n_act[0]*n_act_f[0][1]+n_act[0]*n_act[0]*n_act_f[0][1]+n_act[0]*n_act_f[0][1]+n_act_f[0][1];
    for(int i=0;i<ci1_lt_aldet.size();i++){
    for(int j=0;j<ci2_lt_aldet.size();j++){

//         if(false)//0
        if(ci1_lt_aldet[i].n_A==ci2_lt_aldet[j].n_A)
        if(ci1_lt_aldet[i].n_B==ci2_lt_aldet[j].n_B){
            aldet_S_calc(ACT_DET_A ,&(ci1_lt_aldet[i].A), 0, &(ci2_lt_aldet[j].A), 0);
            aldet_S_calc(ACT_DET_A2 ,&(ci1_lt_aldet[i].B), 0, &(ci2_lt_aldet[j].B), 0);
            
            for(int a=0;a<A->n_states[0];a++)
            for(int b=0;b<A->n_states[1];b++)
            for(int c=0;c<B->n_states[0];c++)
            for(int d=0;d<B->n_states[1];d++){
                abcd = ((a*A->n_states[1]+b)*B->n_states[0]+c)*B->n_states[1]+d;
                ac   =   a*B->n_states[0]+c;
                bd   =   b*B->n_states[1]+d;
                
                S_new  [abcd] += ACT_DET_A [ac]*ACT_DET_A2[bd];
                
            }
        }
        
    }
            // fprintf(stderr,"\r%d/%d",i , ci1_lt_aldet.size());
    }
    // fprintf(stderr,"\n");
    
    
    return 0;
}


int doCI_data::S_calc(){
    
    aldet.calc_S(ACT_DET_A,0,1);
        
    return 0;
}

int doCI_data::calc_1el_DM(){
    
    set_zero_matr(ACT_ADJ_A,n_st_total[0]*n_act[0]*n_st_total[1]*n_act[1]);
    aldet.calc_DMA(ACT_ADJ_A,0,1);
//     printf_timer("1A matr calculation by aldet");

    set_zero_matr(ACT_ADJ_B,n_st_total[0]*n_act[0]*n_st_total[1]*n_act[1]);
    aldet.calc_DMB(ACT_ADJ_B,0,1);
//     printf_timer("1B matr calculation by aldet");
    
    
    return 0;
    
}

int doCI_data::average_DM(std::vector<double> avecoe){
    
    
    double norm=avecoe[0];
    for(int j=0;j<n_act[0]*n_act[0];j++){
        ACT_ADJ_A[j]=ACT_ADJ_A[j]*avecoe[0];
        ACT_ADJ_B[j]=ACT_ADJ_B[j]*avecoe[0];
    }
    for(int i=1;i<avecoe.size();i++){
        norm+=avecoe[i];
        for(int j=0;j<n_act[0]*n_act[0];j++){
            ACT_ADJ_A[j]+=ACT_ADJ_A[(i*n_st_total[1]+i)*n_act[0]*n_act[0]+j]*avecoe[i];
            ACT_ADJ_B[j]+=ACT_ADJ_B[(i*n_st_total[1]+i)*n_act[0]*n_act[0]+j]*avecoe[i];
        }
    }
    for(int j=0;j<n_act[0]*n_act[0];j++)
        ACT_ADJ_A[j]=ACT_ADJ_A[j]/norm;
    for(int j=0;j<n_act[0]*n_act[0];j++)
        ACT_ADJ_B[j]=ACT_ADJ_B[j]/norm;
        
//     fprintf(out_stream,"alp:\n");
//     PrintMatr(ACT_ADJ_A,n_act[0],n_act[0],1);
//     fprintf(out_stream,"bet:\n");
//     PrintMatr(ACT_ADJ_B,n_act[0],n_act[0],1);
    
//     getchar();
//     fprintf(out_stream,"norm = %e\n",norm);
    return 0;
}

int doCI_data::nat_orb_calc_1mol(int i, int renorm){
    
    for(int j=0;j<n_act[0]*n_act[0];j++)
        BUF[j]=ACT_ADJ_A[(i*n_st_total[1]+i)*n_act[0]*n_act[0]+j]+
               ACT_ADJ_B[(i*n_st_total[1]+i)*n_act[0]*n_act[0]+j];
    
    A->nat_orb_calc(BUF, renorm,'T');
    
    return 0;
}

int doCI_data::transform_1el_DM(int N){
    
//     if(n_act[0])
    d_o_DMV_gen(DM_A,
                L_ACT, n_ao[0],
                R_ACT, n_ao[1],
                ACT_ADJ_A+N*n_act[0]*n_act[0], n_act[0],BUF);
        
//     if(n_act[0])
    d_o_DMV_gen(DM_B,
                L_ACT, n_ao[0],
                R_ACT, n_ao[1],
                ACT_ADJ_B+N*n_act[0]*n_act[0], n_act[0],BUF);
    
    
    for(int i_n2=0;i_n2<n_ao[0]*n_ao[1];i_n2++)DM_T [i_n2]=DM_C[i_n2]*ACT_DET_A[N]*2+DM_B[i_n2]+DM_A[i_n2];
    for(int i_n2=0;i_n2<n_ao[0]*n_ao[1];i_n2++)DM_TA[i_n2]=DM_C[i_n2]*ACT_DET_A[N]             +DM_A[i_n2];
    for(int i_n2=0;i_n2<n_ao[0]*n_ao[1];i_n2++)DM_TB[i_n2]=DM_C[i_n2]*ACT_DET_A[N]  +DM_B[i_n2]            ;
     
    return 0;
    
}   

int doCI_data::V_2el_calc(){
    
//     double * ci_new_1;
//     double * ci_new_2;
    
//     fprintf(out_stream,"1\n");
    
//     ci_from_hash(n_act[0],n_alp_el[0],n_bet_el[0], ci1, 0, ci_new_1);
//     fprintf(out_stream,"2\n");
//     ci_from_hash(n_act[0],n_alp_el[0],n_bet_el[0], ci2, 0, ci_new_2);
//     fprintf(out_stream,"3\n");
//     h_2e[0]=calc_2e_ints(12, 6, 6, ci_new_1, ci_new_2, act_INTS,0,4);
//     fprintf(out_stream,"V_new = %.16e\n",calc_2e_ints(n_act[0],n_alp_el[0],n_bet_el[0], ci_new_1, ci_new_2, act_INTS,0,4));
//     
//     delete[]ci_new_1;
//     delete[]ci_new_2;
    
    
    calc_CI_2el_all(h_2e, &ci1, &ci2, act_INTS, n_act[0], n_act[0], n_st_total[0], n_st_total[1]);
//     fprintf(out_stream,"V_old = %.16e\n",h_2e[0]);
    
    
    printf_timer("V calculation");
//     exit(1);

    return 0;
}

double doCI_data::H_2el_calc(double S){

    return E_1el_calc(J, DM_C, n_ao[0], n_ao[1])*S*p_SVD*2+
           E_1el_calc(J, DM_A,n_ao[0], n_ao[1])*  p_SVD*2+
           E_1el_calc(J, DM_B,n_ao[0], n_ao[1])*  p_SVD*2-
           E_1el_calc(K, DM_C, n_ao[0], n_ao[1])*S*p_SVD-
           E_1el_calc(K, DM_A,n_ao[0], n_ao[1])*  p_SVD-
           E_1el_calc(K, DM_B,n_ao[0], n_ao[1])*  p_SVD;
}

int doCI_data::clear(){// to be changed to the destructor!
    
    delete[] n_alp_el_f;
    delete[] n_bet_el_f;
    delete[] n_act_f   ;
    delete[] n_cor     ;
    delete[]n_ao      ;
    delete[] n_alp_el  ;
    delete[] n_bet_el  ;
    delete[] n_act     ;
    
    delete[] S_DD ;
    delete[] L_MO ;
    delete[] R_MO ;
    delete[] L_COR;
    delete[] R_COR;
    delete[] L_ACT;
    delete[] R_ACT;
    delete[] SVD  ;
    delete[] BUF  ;
    delete[] S_SV ;
    delete[] L_ACT_F;
    delete[] R_ACT_F;
    delete[] S_FV ;
    delete[] U_ACT;
    delete[] V_ACT;
    delete[] S_SA ;
    delete[] S_SB ;
    delete[]  H_ACT;
    delete[] Dx_ACT;
    delete[] Dy_ACT;
    delete[] Dz_ACT;        
        
    
    
    delete[] DM_C ;//double occ
    delete[] DM_A;//alp
    delete[] DM_B;//bet
    delete[] DM_T ;//total for 1-el ints
    delete[] DM_TA;//total for 1-el ints alp
    delete[] DM_TB;//total for 1-el ints bet
    delete[] n_st_total;
    
    for(int i=0; i<n_frag;i++){
        for(int i_S=0; i_S<n_states[0][i]; i_S++) {
            ci1_f[i][i_S].clear();
        }
        delete[] ci1_f[i];
    }
    delete[] ci1_f;
    
    for(int i=0; i<n_frag;i++){
        for(int i_S=0; i_S<n_states[1][i]; i_S++) {
            ci2_f[i][i_S].clear();
        }
        delete[] ci2_f[i];
    }
    delete[] ci2_f;
    
    delete[] n_states;
    ci1.clear();
    ci2.clear();
    delete[] h_2e;
    delete[] ACT_DET_A;
    delete[] ACT_ADJ_A;
    delete[] ACT_ADJ_B;
    delete[] J;
    delete[] K;
    delete[] act_INTS;
    
    delete[] VEC;
    
    
    for(int i=0; i<ci1_lt_aldet.size();i++)ci1_lt_aldet[i].clear();
    for(int i=0; i<ci2_lt_aldet.size();i++)ci2_lt_aldet[i].clear();
    
    ci1_lt_aldet.clear();
    ci2_lt_aldet.clear();
    
    if(ACT_DMT_A!=NULL)delete[] ACT_DMT_A;
    if(ACT_DMT_B!=NULL)delete[] ACT_DMT_B;
    if(S_new    !=NULL)delete[] S_new;
    if(H_2el    !=NULL)delete[] H_2el;
    
//     aldet.~aldet_data();
    
    return 0;
}
    
