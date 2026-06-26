# include "blas_link.h"
# include "libint_link.h"
# include "matr.h"
# include "timer.h"
# include "defaults.h"
# include "doCI_data.h"
# include "UPT_tensors.h"
# include "RI.h"
# include "inp_out.h"
# include "XMCQDPT.h"
# include "CAS.h"

# define max(a,b)  (((a)<(b))?(b):(a))

int three_block_diagonalization_no_change_A
                               (double * M, double * V, double * e1,
                                int d1, int d2, int d3, int d_all,
                                double * O1, double * O2, double* O3, double * O4){
    
//     int d3 = d_all-d1-d2;
    int dm = max(d3,max(d1,d2));
    double * e2 = e1+d1;
    double * e3 = e2+d2;
//     fprintf(out_stream,"M:\n");
//     fPrintMatr(out_stream,M,d_all,d_all,1);
    
    double * M1,* M2,* M3;
    M1 = new double[d1*d1];
    M2 = new double[d2*d2];
    M3 = new double[d3*d3];
    
    for(int i=0;i<d1;i++)
    for(int j=0;j<d1;j++)
        M1[i*d1+j]=M[i*d_all+j];
    
    for(int i=0;i<d2;i++)
    for(int j=0;j<d2;j++)
        M2[i*d2+j]=M[(i+d1)*d_all+j+d1];
        
    for(int i=0;i<d3;i++)
    for(int j=0;j<d3;j++)
        M3[i*d3+j]=M[(i+d1+d2)*d_all+j+d1+d2];
    
    if(d1!=0)lapack_diag(M1,e1,d1);
    if(d2!=0)lapack_diag(M2,e2,d2);
    if(d3!=0)lapack_diag(M3,e3,d3);
    
    set_zero_matr(O4,d2*d2);
    for(int i=0;i<d2;i++)
        O4[i*d2+i]=1.0;

    if(d1!=0)
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasTrans,
                     d_all,d1,d1,1.0,
                     V,d_all,
                     M1,d1,0.0,
                     O1,d1);
//     if(d2!=0)
//     cblas_dgemm(CblasRowMajor,CblasTrans,CblasTrans,
//                      d_all,d2,d2,1.0,
//                      V+d_all*d1,d_all,
//                      M2,d2,0.0,
//                      O2,d2);
    if(d3!=0)
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasTrans,
                     d_all,d3,d3,1.0,
                     V+d_all*(d1+d2),d_all,
                     M3,d3,0.0,
                     O3,d3);
    
    for(int i=0;i<d1   ;i++)
    for(int j=0;j<d_all;j++)
        V[i*d_all+j]=O1[j*d1+i];
    
//     for(int i=0;i<d2   ;i++)
//     for(int j=0;j<d_all;j++)
//         V[(i+d1)*d_all+j]=O2[j*d2+i];
    
    for(int i=0;i<d2   ;i++)
    for(int j=0;j<d_all;j++)
        O2[j*d2+i]=V[(i+d1)*d_all+j];

    for(int i=0;i<d3   ;i++)
    for(int j=0;j<d_all;j++)
        V[(i+d1+d2)*d_all+j]=O3[j*d3+i];
    
    
    delete[] M1;
    delete[] M2;
    delete[] M3;
//     delete[] BUF;
    
    return 1;
    
}



# define max(a,b)  (((a)<(b))?(b):(a))

extern int num_threads;

// extern int testing;

inline double ED_with_shift(double E, double edshift){
    
    return E/(E*E+edshift);
}

int print_orb_energies(const char *name, double * eps, int n_cor, int n_act, int n_virt, cdas_par * cdas, int n_HOMO){
    
    double * eps_a;
    double * eps_e;
    
    eps_a = eps   + n_cor;
    eps_e = eps_a + n_act;
    
    printf("\n");
    printf("Orbital energies (3 blocks) %s:\n", name);
    printf("core     :");PrintMatr(eps  ,1,n_cor,0);
    printf("active F :");PrintMatr(eps_a,1,n_act,0);
    if(cdas->have_eps){
        for(int i=0; i<n_act;i++)eps_a[i]=cdas->eps[i];
    }
    if(cdas->HOMO){
        for(int i=0; i<n_act;i++)eps_a[i]=eps_a[n_HOMO];//eps_a[M->n_act_el_alp[0]-1];
    }
    if(cdas->sing_e){
        for(int i=0; i<n_act;i++)eps_a[i]=cdas->eps[0];
    }
    if(cdas->orb_e){
        for(int i=0; i<n_act;i++)eps_a[i]=eps_a[cdas->n_orb-1];
    }
    
    if(cdas->have_eps+cdas->HOMO+cdas->sing_e+cdas->orb_e){
        printf("active PT:");PrintMatr(eps_a,1,n_act,0);
    }
    printf("vacant   :");PrintMatr(eps_e,1,n_virt,0);
    printf("\n\n");
    
    
    return 0;
    
}

int U_three_block_diagonalization(double * M, double * V, double * e1,
                                  int d1, int d2, int d3, int d_all,
                                  double * O1, double * O2, double* O3, double * O4){
    
//     int d3 = d_all-d1-d2;
    int dm = max(d3,max(d1,d2));
    double * e2 = e1+d1;
    double * e3 = e2+d2;
//     printf("M:\n");
//     PrintMatr(M,d_all,d_all,1);
    
    double * M1,* M2,* M3;
    M1 = new double[d1*d1];
    M2 = new double[d2*d2];
    M3 = new double[d3*d3];
    
    for(int i=0;i<d1;i++)
    for(int j=0;j<d1;j++)
        M1[i*d1+j]=M[i*d_all+j];
    
    for(int i=0;i<d2;i++)
    for(int j=0;j<d2;j++)
        M2[i*d2+j]=M[(i+d1)*d_all+j+d1];
        
    for(int i=0;i<d3;i++)
    for(int j=0;j<d3;j++)
        M3[i*d3+j]=M[(i+d1+d2)*d_all+j+d1+d2];
    
    if(d1!=0)lapack_diag(M1,e1,d1);
    if(d2!=0)lapack_diag(M2,e2,d2);
    if(d3!=0)lapack_diag(M3,e3,d3);
    
    for(int i=0;i<d2;i++)
    for(int j=0;j<d2;j++)
        O4[i*d2+j]=M2[j*d2+i];

    if(d1!=0)
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasTrans,
                     d_all,d1,d1,1.0,
                     V,d_all,
                     M1,d1,0.0,
                     O1,d1);
//     if(d2!=0)
//     cblas_dgemm(CblasRowMajor,CblasTrans,CblasTrans,
//                      d_all,d2,d2,1.0,
//                      V+d_all*d1,d_all,
//                      M2,d2,0.0,
//                      O2,d2);
    if(d3!=0)
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasTrans,
                     d_all,d3,d3,1.0,
                     V+d_all*(d1+d2),d_all,
                     M3,d3,0.0,
                     O3,d3);
    
    for(int i=0;i<d1   ;i++)
    for(int j=0;j<d_all;j++)
        V[i*d_all+j]=O1[j*d1+i];
    
    for(int i=0;i<d2   ;i++)
    for(int j=0;j<d_all;j++)
        O2[j*d2+i]=V[(i+d1)*d_all+j];
    
    for(int i=0;i<d3   ;i++)
    for(int j=0;j<d_all;j++)
        V[(i+d1+d2)*d_all+j]=O3[j*d3+i];
    
    
    delete[] M1;
    delete[] M2;
    delete[] M3;
//     delete[] BUF;
    
    return 1;
    
}

int U_CDAS_PT2(molecule * M, cdas_par * cdas, char * job_name){
    
    printf("U_CDAS_PT2 is under refactoring\n");
    
    
   
    if(RI){
        gen_RI_AA(M);
//         printf("\n\n");
//         printf_timer("RI integrals calculation");
    }
    else{
        fprintf(out_stream,"ERROR: RI=0 is not supported for any U_PT\n");
        exit(0);
    }
    
    
    cdas->write_info(M->n_act_el_alp[0],
                   M->n_act_el_bet[0],
                   M->n_act_orb   [0],
                   M->CI[0].mult);
    if(fabs(cdas->edshift)>1E-8)
        fprintf(out_stream,"\n\nWARNING: nonzero ISA denominator shift is not typical for CDAS-PT2\n\n");
    
//     libint2::initialize();
    
    
    char * name = new char[BUF_LINE_LENGTH];
    doCI_data D;
    
    int n_s=cdas->cas->n_s;
    M->n_states[0]=n_s;
    D.first_alloc(M,M);
    
    D.gen_aldet_data();
    
    
        
    M->mc=1;//set multi-configuration
    
    int n_ao  = D.n_ao[0];
    int n_act = M->n_act_orb[0];
    int n_cor = M->n_cor_orb;
    int n_mo  = n_ao;
    if(D5){
        M->calc_d5_n_ao();
        n_mo = M->n_ao_d5;
    }
    
    RI_core_realloc(n_cor+n_act, n_ao);
        
    
    int n_virt = n_mo-n_cor-n_act;
    
    if(cdas->have_eps){
        if(n_act!=cdas->eps.size()){
            fprintf(out_stream,"ERROR:wrong dim of eps(%d) != n_act(%d)\n",cdas->eps.size(),n_act);
            exit(1);
        }
    }
    if(cdas->orb_e){
        if((cdas->n_orb<1)||(cdas->n_orb>(n_act+1))){
            fprintf(out_stream,"ERROR: could not use orbital %d for active space of size %d\n",cdas->n_orb,n_act);
            exit(0);
        }
    }
    
    D.DM_gen_ortogonal(0, 1.0000);    
    
    D.second_alloc();
    D.cpy_VEC_to_ACT();
    
    D.S_calc();
//     fprintf(out_stream,"Initial CI vectors overlap matrix:\n");
//     fPrintMatr(out_stream,D.ACT_DET_A,D.n_st_total[0],D.n_st_total[1],0);
    
    

    
    D.calc_1el_DM();
    D.average_DM(cdas->cas->w_state);
    M->nat_orb_calc(D.ACT_ADJ_A, 1,'A');
    M->nat_orb_calc(D.ACT_ADJ_B, 1,'B');
    D.transform_1el_DM(0);
    
    fprintf(out_stream,"\n\n");
    fprintf(out_stream,"Calculation Fock-like operator...\n");
    fflush(out_stream);
    
    double * F_A;
    F_A = new double[n_ao*n_ao];
    double * F_B;
    F_B = new double[n_ao*n_ao];
    
    M->calc_U_F_AO(F_A, F_B, D.DM_TA, D.DM_TB, 1);//DM is normalized by 1 in J and by 1 in K
    D.AO_to_MO(F_A);
    D.AO_to_MO(F_B);
//     for(int i=0;i<n_ao*n_ao;i++)F_A[i]=0.5*(F_A[i]+F_B[i]);////!!!!!!!!!!!!!!!!!!!
//     for(int i=0;i<n_ao*n_ao;i++)F_B[i]=     F_A[i]        ;
       
    double * eps_A;
    double * eps_B;
    eps_A = new double[n_ao];
    eps_B = new double[n_ao];
    
//     //transposed vectors for calculation of 2e integrals
    double * COR_VEC_A;//core
    double * ACT_VEC_A;//active
    double * VIRT_VEC_A;//virtual
    
    double * COR_VEC_B;//core
    double * ACT_VEC_B;//active
    double * VIRT_VEC_B;//virtual
    
    int n_int=n_cor+n_act;
    
    COR_VEC_A = new double[n_ao * n_cor];
    ACT_VEC_A = new double[n_ao * n_act];
    VIRT_VEC_A = new double[n_ao * n_virt];

    COR_VEC_B = new double[n_ao * n_cor];
    ACT_VEC_B = new double[n_ao * n_act];
    VIRT_VEC_B = new double[n_ao * n_virt];

    
    double * U_A;
    double * U_B;
    U_A=new double[n_act*n_act];
    U_B=new double[n_act*n_act];
    
    printf("Diagonalization of Fock-like operator in 3 subspaces...\n");
    fflush(stdout);
    
    if(M->MO_VEC==M->MO_VEC_B){
        M->MO_VEC_B=new double[n_ao*n_ao];
        memcpy(M->MO_VEC_B,M->MO_VEC,n_ao*n_ao*sizeof(double));
    }
    printf("active orbitals are not rotated\n\n");
    
    three_block_diagonalization_no_change_A(F_A,M->MO_VEC  ,eps_A,n_cor,n_act, n_virt,n_ao,COR_VEC_A,ACT_VEC_A,VIRT_VEC_A, U_A);
    three_block_diagonalization_no_change_A(F_B,M->MO_VEC_B,eps_B,n_cor,n_act, n_virt,n_ao,COR_VEC_B,ACT_VEC_B,VIRT_VEC_B, U_B);
    
    print_orb_energies( "for alpha subset", eps_A, n_cor, n_act, n_virt, cdas, M->n_act_el_alp[0]-1);//getchar();
//     print_orb_energies( "for beta  subset", eps_B, n_cor, n_act, n_virt, cdas, M->n_act_el_bet[0]-1);//getchar();
    //testing!!!!!!!!!!!!!!
    print_orb_energies( "for beta  subset", eps_B, n_cor, n_act, n_virt, cdas, M->n_act_el_alp[0]-1);//getchar();
    
// //     printf("Writing orbitals:\n");
// //     VEC_energy_cpy(M->MO_VEC,D.VEC[0], M->orb_energy, eps, n_ao);
// //     M->check_orb_symmetry();
// //     M->MO_libint_back_reordr();
// //     sprintf(name,"%s_CDAS_F.out\0",job_name);
// //     M->GAMESS_type_out_print(name,-1);
// //     printf("visualization file: %s\n",name);
//     
// //     sprintf(name,"%s_CDAS_F.orb\0",job_name);
// //     M->MO_print(name);
// //     printf("data file         : %s\n",name);
// //     M->MO_libint_reordr();
// 
//     
    printf_timer("Fock matrix calculation");
    fprintf(out_stream,"_______________________________________________________________________\n\n\n");
    D.n_st_total[0]=n_s;
    D.n_st_total[1]=n_s;
    
    //CI coef transformation
    D.aldet.malmqvist(0,U_A);
    D.aldet.malmqvist(1,U_A);
//     PrintMatr(U,n_act,n_act,1);
    D.S_calc();

    D.cpy_ACT_VEC_to_ACT(ACT_VEC_A); // copies ACT to L_ACT and R_ACT

    D.calc_1el_DM();// possibly unneeded?????
    
    double * H_AO;
    H_AO=F_A;
    
    double * J = new double[n_ao*n_ao];
    double * K = new double[n_ao*n_ao];
    double * act_INTS_A = new double[n_act*n_act*n_act*n_act];
    double * act_INTS_B = new double[n_act*n_act*n_act*n_act];
    double * act_INTS_AB= new double[n_act*n_act*n_act*n_act];
    set_zero_matr(act_INTS_A ,n_act*n_act*n_act*n_act);
    set_zero_matr(act_INTS_B ,n_act*n_act*n_act*n_act);
    set_zero_matr(act_INTS_AB,n_act*n_act*n_act*n_act);
    set_zero_matr(J ,n_ao*n_ao);
    set_zero_matr(K ,n_ao*n_ao);
    
    
    calc_2el_MO_INTS_RI(   n_ao, D.DM_C, J, K,  act_INTS_A, M->MO_VEC  , M->MO_VEC  +n_ao*n_cor, M->MO_VEC  +n_ao*n_cor, n_cor, n_act,0);
    calc_2el_MO_INTS_RI(   n_ao, NULL,NULL,NULL,act_INTS_B, M->MO_VEC_B, M->MO_VEC_B+n_ao*n_cor, M->MO_VEC_B+n_ao*n_cor, n_cor, n_act,0);
    
    for(int i=0;i<n_ao*n_ao;i++)H_AO[i]=M->H_AO[i]+2*J[i]-K[i];

//     //H1 
    double * H_AV_A;
    double * H_AA_A;
    double * H_CA_A;
    double * H_CV_A;
    H_AV_A = new double[n_act*n_virt];
    H_AA_A = new double[n_act*n_act];
    H_CA_A = new double[n_cor*n_act];
    H_CV_A = new double[n_cor*n_virt];

    transform_from_col_MO(H_AV_A, H_AO, n_ao, ACT_VEC_A, n_act, VIRT_VEC_A, n_virt);
    transform_from_col_MO(H_AA_A, H_AO, n_ao, ACT_VEC_A, n_act, ACT_VEC_A, n_act);
    transform_from_col_MO(H_CA_A, H_AO, n_ao, COR_VEC_A, n_cor, ACT_VEC_A, n_act);
    transform_from_col_MO(H_CV_A, H_AO, n_ao, COR_VEC_A, n_cor, VIRT_VEC_A, n_virt);

    double * H_AV_B;
    double * H_AA_B;
    double * H_CA_B;
    double * H_CV_B;
    H_AV_B = new double[n_act*n_virt];
    H_AA_B = new double[n_act*n_act];
    H_CA_B = new double[n_cor*n_act];
    H_CV_B = new double[n_cor*n_virt];

    transform_from_col_MO(H_AV_B, H_AO, n_ao, ACT_VEC_B, n_act, VIRT_VEC_B, n_virt);
    transform_from_col_MO(H_AA_B, H_AO, n_ao, ACT_VEC_B, n_act, ACT_VEC_B, n_act);
    transform_from_col_MO(H_CA_B, H_AO, n_ao, COR_VEC_B, n_cor, ACT_VEC_B, n_act);
    transform_from_col_MO(H_CV_B, H_AO, n_ao, COR_VEC_B, n_cor, VIRT_VEC_B, n_virt);

    
        
    printf_timer("before PT");



    printf("Start calculation of PT2 correction (unrestricted)\n");
    
    RI_data R_A;
    RI_data R_B;
    R_A.set_par(M,n_cor, n_act, n_virt);
    R_B.set_par(M,n_cor, n_act, n_virt);
    R_A.MO_calc(M->MO_VEC  , n_ao);
    R_B.MO_calc(M->MO_VEC_B, n_ao);
    
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                   n_act*n_act,n_act*n_act,R_A.aux_n_ao,1.0,
                   R_A.AA_RI_M,R_A.aux_n_ao,
                   R_B.AA_RI_M,R_A.aux_n_ao,0.0,
                   act_INTS_AB,n_act*n_act);//check order!!!!    
    
    
    printf("\n");
    printf("\n");
    printf_timer("RI orbital transformation");
    printf("\n");
    fflush(stdout);
    
    
    //PT
    UPT_tensors T;
    
    T.set_par(&R_A, &R_B, 
              eps_A, eps_B, 
              n_cor, n_act, n_virt, 
              H_AV_A, H_CA_A, H_CV_A, 
              H_AV_B, H_CA_B, H_CV_B, 
              cdas->edshift);
    
    //calculation of IP/EA Fockian matrix
//     double zero=0;
    D.aldet.U_simple_import_data(act_INTS_A, act_INTS_AB, act_INTS_B, 
                                 H_AA_A,H_AA_B,
                                 0);///change it !!!!
    
//     if(cdas->IPEA){
//         T.IPEA(&(D.aldet), 0,cdas->cas->w_state);
//         T.E2_calc_IPEA();
//     }
//     else 
    if(cdas->MPPT){
        T.MPPT(&(D.aldet), 0,cdas->cas->w_state);
        T.E2_calc_EE();
    }
    else{
        T.E2_calc_EE();
    }
        
    printf_timer("PT tensors calculation");
    printf("_______________________________________________________________________\n\n\n");
// 
//     
    CAS_engine CAS;
    CAS.init(cdas->cas ,M);
//     
//     CAS.E_core =  E_1el_calc(M->H_AO, D.DM_C, n_ao, n_ao)*2
//                  +E_1el_calc(    J  , D.DM_C, n_ao, n_ao)*2
//                  -E_1el_calc(    K  , D.DM_C, n_ao, n_ao)  
//                  +M->V_nuc ;
// 
//     CAS.CI->as_aldet()-> (act_INTS_A, act_INTS_AB, act_INTS_B, 
//                                  H_AA_A,H_AA_B,
//                                  CAS.E_core);
//     
    CAS.CI->as_aldet()->UPT2_import_data(T.RF_P3_JK_A,   
                             T.RF_P3_AB,
                             T.RF_P3_BA,
                             T.RF_P3_JK_B,   
                             T.RF_PV_JK_A,
                             T.RF_PV_AB,
                             T.RF_PV_JK_B,
                             T.RF_PH_A,
                             T.RF_PH_B,
                             T.RF_PS);
    
    CAS.CI_calc(1,0,1);
    
    fprintf(out_stream,"\n\nCDAS-PT2 Energy summary:\n");
    PrintEnergy(CAS.CI->as_aldet()->E_states[0],CAS.n_s,1);

                 
//     davidson_solver DAV2;
//     
//     DAV2.set_par(CAS.CI->as_aldet(),CAS.dav);
//     DAV2.V.H_diag_calc();
    
// //     PrintMatr(DAV2.V.H_diag,DAV2.V.Nd,1,0);
// //     exit(0);
// 
//     int n;
//     n = DAV2.run(1,0);
//     
//     printf("\n\nU-CDAS-PT2 Energy summary:\n");
//     PrintEnergy(CAS.CI->as_aldet()->E_states[0],CAS.n_s,1);
//         
//     CAS.Prop_calc();
// //     getchar();
//     printf("\n");
//     printf("Writing CDAS-PT2 WaveFunctions:\n");
//     sprintf(name,"%s_CDAS.ci\0",job_name);
//     CAS.CI->as_aldet()->write_civec(0, name);
//     printf("data file         : %s\n",name);
// 
//     
// //     if(cdas->IPEA==0){
// //         printf("\n\nCalculation of bias correction (derivative by orbital energy):\n\n");
// //         double * F_act_bc      = new double[n_act*n_act];
// //         double * act_INTS_bc   = new double[n_act*n_act*n_act*n_act];
// //         double * act_INTS_2_bc = new double[n_act*n_act*n_act*n_act];
// // //         double * bc_orb_der   = new double[n_act*n_s];
// //         aldet_data CI_bc;
// //         CI_bc.get_dim(CAS.CI->as_aldet()->n_act, CAS.CI->as_aldet()->na, CAS.CI->as_aldet()->nb, 1, CAS.CI->as_aldet()->mult, CAS.CI->as_aldet()->print_number);
// //         CI_bc.init_zero_vec(n_s,0);
// //         
// //         double E_core_bc=0;
// //         
// //         set_zero_matr(F_act_bc ,n_act*n_act);
// //         set_zero_matr(act_INTS_bc ,n_act*n_act*n_act*n_act);
// //         set_zero_matr(act_INTS_2_bc ,n_act*n_act*n_act*n_act);
// //         CI_bc.E_core  = &E_core_bc ;
// //         CI_bc.F_act     =F_act_bc     ;
// //         CI_bc.act_INTS  =act_INTS_bc  ;
// //         CI_bc.act_INTS_2=act_INTS_2_bc;
// //         
// //         PT_tensors T2;
// //         T2.set_par(&R, eps, n_cor, n_act, n_virt, H_AV, H_CA, H_CV, cdas->edshift);
// //         T2.E2_calc_bc();
// //         
// //         CI_bc.PT2_import_data(T2.RF_P3_JK,
// //                               T2.RF_P3_AB,
// //                               T2.RF_PV_JK,
// //                               T2.RF_PV_AB,
// //                               T2.RF_PH,
// //                               T2.RF_PS);
// //         
// //         CI_bc.gen_ext_ind();
// //         CI_bc.H_diag_calc();
// //         CI_bc.do_PT=1;
// //         CI_bc.PT_update();
// //         memcpy(CI_bc.coef[0],CAS.CI->as_aldet()->coef[0],n_s*CI_bc.Nd*sizeof(double));
// //         CI_bc.transpose_ci(0);
// //         CI_bc.H_mult(0,CAS.n_s);
// //         cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
// //                     n_s,n_s,CI_bc.Nd,1.0,
// //                     CI_bc.coef[0],n_s,
// //                     CI_bc.Hv,n_s,0.0,
// //                     H_bc,n_s);
// //         PrintEnergy_derivative_corrected(CAS.CI->as_aldet()->E_states[0],H_bc,CAS.n_s,0);
// //         
// //         delete[] F_act_bc     ;
// //         delete[] act_INTS_bc  ;
// //         delete[] act_INTS_2_bc;
// // //         delete[] bc_orb_der  ;
// //         
// //     }
//     
//     CAS.print_properties("U-CDAS-PT(0)");    
//         
//     printf_timer("U-CDAS-PT2");
// //     getchar();
//     
//     
//     D.clear();
//     
//     delete[] H_CV_A;
//     delete[] H_CA_A;
//     delete[] H_AA_A;
//     delete[] H_AV_A;
//     delete[] H_CV_B;
//     delete[] H_CA_B;
//     delete[] H_AA_B;
//     delete[] H_AV_B;
//     delete[] F_A;
//     delete[] F_B;
//     delete[] COR_VEC_A;
//     delete[] ACT_VEC_A;
//     delete[] VIRT_VEC_A;
//     delete[] COR_VEC_B;
//     delete[] ACT_VEC_B;
//     delete[] VIRT_VEC_B;
//     
//     delete[] H_bc;
//     delete[] eps_A;
//     delete[] eps_B;
// 
// //     delete[] eps_EA;
//     
//     delete[] U_A;
//     delete[] U_B;
//     delete[] name;
//     
//     delete[] J         ;
//     delete[] K         ;
//     delete[] act_INTS_A ;
//     delete[] act_INTS_B ;
//     delete[] act_INTS_AB;
//     
    return 0;
 
}
