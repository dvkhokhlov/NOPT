# include "blas_link.h"
# include "libint_link.h"
# include "matr.h"
# include "timer.h"
#include "defaults.h"
# include "doCI_data.h"
# include "PT_tensors.h"
# include "RI.h"
# include "inp_out.h"
# include "XMCQDPT.h"
# include "CAS.h"
# include "ecp.h"
# include "aldet_rel.h"

# define max(a,b)  (((a)<(b))?(b):(a))

extern int num_threads;

double * H_ci;
void print_matrix_lower_triangle(char *path, int dim, double *matrix)
{
    double const PRINT_THRESH = 1e-14;

    FILE *out_file = fopen(path, "w");

    fprintf(out_file, "%6d\n", dim);

    for (int i = 0; i < dim; i++) {
        for (int j = 0; j <= i; j++) {
            if (fabs(matrix[i * dim + j]) > PRINT_THRESH) {
                fprintf(out_file, "%6d%6d%30.8E\n", i + 1, j + 1, matrix[i * dim + j]);
            }
        }
    }

    fclose(out_file);
}

// extern int testing;

inline double ED_with_shift(double E, double edshift){
    
    return E/(E*E+edshift);
}

int H_add_SO(double * Hr, double * Hi, double * H_sc, double * Sx, double * Sy, double * Sz, int n){
    
    
    set_zero_matr(Hr, 4*n*n);
    set_zero_matr(Hi, 4*n*n);
      
    if(H_sc!=nullptr)
    for(int i=0;i<n;i++)
    for(int j=0;j<n;j++)
        Hr[i*2*n+j]=H_sc[i*n+j];
    if(H_sc!=nullptr)
    for(int i=0;i<n;i++)
    for(int j=0;j<n;j++)
        Hr[(i+n)*2*n+j+n]=H_sc[i*n+j];
    
    
    for(int i=0;i<n;i++)
    for(int j=0;j<n;j++)
        Hi[i*2*n+j] += Sz[i*n+j]*0.5;
    
    for(int i=0;i<n;i++)
    for(int j=0;j<n;j++)
        Hi[(i+n)*2*n+j+n]-=Sz[i*n+j]*0.5;
    
    
    
    for(int i=0;i<n;i++)
    for(int j=0;j<n;j++)
        Hi[i*2*n+j+n]+= Sx[i*n+j]*0.5;
    
    for(int i=0;i<n;i++)
    for(int j=0;j<n;j++)
        Hi[(i+n)*2*n+j]+=Sx[i*n+j]*0.5;
    
    
    for(int i=0;i<n;i++)
    for(int j=0;j<n;j++)
        Hr[i*2*n+j+n]+= Sy[i*n+j]*0.5;
    
    for(int i=0;i<n;i++)
    for(int j=0;j<n;j++)
        Hr[(i+n)*2*n+j]-=Sy[i*n+j]*0.5;
    
    
    
    
    
    return 0;
}


int CDAS_PT2_rel(molecule * M, cdas_par * cdas, char * job_name){
#ifdef _USE_GRPP    
    if(RI){
        gen_RI_AA(M);
    }
    else{
        fprintf(out_stream,"ERROR: RI=0 is not supported for any PT\n");
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
    
    
    double * F   = new double[n_ao*n_ao];
    double * eps = new double[n_ao];
    
    //transposed vectors for calculation of 2e integrals
    double * COR_VEC;//core
    double * ACT_VEC;//active
    double * VIRT_VEC;//virtual
    double * OCC_VEC;//occupied=cor+act
    
    int n_occ=n_cor+n_act;
    
    COR_VEC = new double[n_ao * n_cor];
    ACT_VEC = new double[n_ao * n_act];
    OCC_VEC = new double[n_ao * n_occ];
    VIRT_VEC = new double[n_ao * n_virt];
    
    double * U = new double[n_act*n_act];
    
    double * J = new double[n_ao*n_ao];
    double * K = new double[n_ao*n_ao];
    
    double * act_INTS   = new double[n_act*n_act*n_act*n_act];
    double * AAAG_INTS = new double[n_ao *n_act*n_act*n_act];
    
    
    double * d_x1= new double[n_s*n_s];
    double * d_y1= new double[n_s*n_s];
    double * d_z1= new double[n_s*n_s];
    double * d_x_CV = new double[n_cor*n_virt];
    double * d_y_CV = new double[n_cor*n_virt];
    double * d_z_CV = new double[n_cor*n_virt];
    double * d_x_AV = new double[n_act*n_virt];
    double * d_y_AV = new double[n_act*n_virt];
    double * d_z_AV = new double[n_act*n_virt];
    double * d_x_CA = new double[n_cor*n_act];
    double * d_y_CA = new double[n_cor*n_act];
    double * d_z_CA = new double[n_cor*n_act];
        
    //H1 
    double * H_AV;
    double * H_AA;
    double * H_CA;
    double * H_CV;
    H_AV = new double[n_act*n_virt];
    H_AA = new double[n_act*n_act];
    H_CA = new double[n_cor*n_act];
    H_CV = new double[n_cor*n_virt];
    
    
    CAS_engine CAS;
    CAS.init(cdas->cas,M);
    
//     for(int i_s=0;i_s<3;i_s++)
//     {
        
//         for(auto &w: cdas->cas->w_state)w=0;
//         if(i_s==0)cdas->cas->w_state[0]=1.0;
//         if(i_s==1)cdas->cas->w_state[2]=1.0;
//         if(i_s==2)cdas->cas->w_state[1]=1.0;
        fprintf(out_stream,"Weights:\t");
        for(auto &w: cdas->cas->w_state) fprintf(out_stream,"% .2e ",w);
        fprintf(out_stream,"\n");
        
        
        D.calc_1el_DM();
        D.average_DM(cdas->cas->w_state);
        D.nat_orb_calc_1mol(0,1);
        D.transform_1el_DM(0);
        
        fprintf(out_stream,"\n\n");
        fprintf(out_stream,"Calculation of Fock-like operator...\n");
        fflush(out_stream);
        
        
        M->calc_F_AO(F, D.DM_T, 0.5);
        double * F_backup = new double[n_ao*n_ao];
        memcpy(F_backup,F, n_ao*n_ao*sizeof(double));
        
        D.AO_to_MO(F);
           
        
        
        fprintf(out_stream,"Diagonalization of Fock-like operator in 3 subspaces...\n");
        fflush(out_stream);
//         if(cdas->have_eps)
            three_block_diagonalization(F,D.VEC[0],eps,n_cor,n_act, n_virt,n_ao,COR_VEC,ACT_VEC,VIRT_VEC, U);
//         else
//             three_block_diagonalization_no_change_A(F,D.VEC[0],eps,n_cor,n_act, n_virt,n_ao,COR_VEC,ACT_VEC,VIRT_VEC, U);
//         for(int i = 0; i<2   ;i++)
//         for(int j = 0; j<n_ao;j++)
//             D.VEC[0][(i+n_cor+n_act)*n_ao+j]=0;
            
            
        double * eps_a = eps   + n_cor;
        double * eps_e = eps_a + n_act;

        fprintf(out_stream,"\n");
        fprintf(out_stream,"Orbital energies (3 blocks):\n");
        fprintf(out_stream,"core     :");fPrintMatr(out_stream,eps  ,1,n_cor,0);
        fprintf(out_stream,"active F :");fPrintMatr(out_stream,eps_a,1,n_act,0);
        if(cdas->have_eps){
            for(int i=0; i<n_act;i++)eps_a[i]=cdas->eps[i];
        }
        if(cdas->HOMO){
            for(int i=0; i<n_act;i++)eps_a[i]=eps_a[M->n_act_el_alp[0]-1];
        }
        if(cdas->sing_e){
            for(int i=0; i<n_act;i++)eps_a[i]=cdas->eps[0];
        }
        if(cdas->orb_e){
            for(int i=0; i<n_act;i++)eps_a[i]=eps_a[cdas->n_orb-1];
        }
        
        if(cdas->have_eps+cdas->HOMO+cdas->sing_e+cdas->orb_e){
            fprintf(out_stream,"active PT:");fPrintMatr(out_stream,eps_a,1,n_act,0);
        }
        
        
        
        fprintf(out_stream,"vacant   :");fPrintMatr(out_stream,eps_e,1,n_virt,0);
        fprintf(out_stream,"\n\n");
        
        VEC_energy_cpy(M->MO_VEC,D.VEC[0], M->orb_energy, eps, n_ao);
        if(write_orbs){
            fprintf(out_stream,"Writing orbitals:\n");
            if(IS_SYM==1)M->check_orb_symmetry();
            M->MO_gamess_format();
            sprintf(name,"%s_CDAS_F.out\0",job_name);
            M->GAMESS_type_out_print(name,-1);
            fprintf(out_stream,"visualization file: %s\n",name);
            
            sprintf(name,"%s_CDAS_F.orb\0",job_name);
            M->MO_print(name);
            fprintf(out_stream,"data file         : %s\n",name);
//             M->MO_libint_reordr();
        }
        
        printf_timer("Fock matrix calculation");
        fprintf(out_stream,"_______________________________________________________________________\n\n\n");
        D.n_st_total[0]=n_s;
        D.n_st_total[1]=n_s;
        
        //CI coef transformation
        D.aldet.malmqvist(0,U);
        D.aldet.malmqvist(1,U);
        
        D.S_calc();
        
        D.cpy_ACT_VEC_to_ACT(ACT_VEC); // copies ACT to L_ACT and R_ACT
        
        D.calc_1el_DM();
        
        double * H_core;
        H_core=F;
        
        set_zero_matr(act_INTS,n_act*n_act*n_act*n_act);
        set_zero_matr(J ,n_ao*n_ao);
        set_zero_matr(K ,n_ao*n_ao);
        
        calc_2el_MO_INTS_RI(   n_ao, D.DM_C, J, K, act_INTS, M->MO_VEC, M->MO_VEC+n_ao*n_cor, M->MO_VEC+n_ao*n_cor, n_cor, n_act,0);
        
        for(int i=0;i<n_ao*n_ao;i++)H_core[i]=M->H_AO[i]+2*J[i]-K[i];//-F_backup[i];
        
        
        
        
        set_zero_matr(d_x1, n_s*n_s);
        set_zero_matr(d_y1, n_s*n_s);
        set_zero_matr(d_z1, n_s*n_s);    
        
        transform_from_col_MO(d_x_CV, M->Dx_AO, n_ao, COR_VEC, n_cor, VIRT_VEC, n_virt);
        transform_from_col_MO(d_y_CV, M->Dy_AO, n_ao, COR_VEC, n_cor, VIRT_VEC, n_virt);
        transform_from_col_MO(d_z_CV, M->Dz_AO, n_ao, COR_VEC, n_cor, VIRT_VEC, n_virt);
        
        transform_from_col_MO(d_x_AV, M->Dx_AO, n_ao, ACT_VEC, n_act, VIRT_VEC, n_virt);
        transform_from_col_MO(d_y_AV, M->Dy_AO, n_ao, ACT_VEC, n_act, VIRT_VEC, n_virt);
        transform_from_col_MO(d_z_AV, M->Dz_AO, n_ao, ACT_VEC, n_act, VIRT_VEC, n_virt);
        
        transform_from_col_MO(d_x_CA, M->Dx_AO, n_ao, COR_VEC, n_cor, ACT_VEC, n_act);
        transform_from_col_MO(d_y_CA, M->Dy_AO, n_ao, COR_VEC, n_cor, ACT_VEC, n_act);
        transform_from_col_MO(d_z_CA, M->Dz_AO, n_ao, COR_VEC, n_cor, ACT_VEC, n_act);
        
        
        
        //H1
        transform_from_col_MO(H_AV, H_core, n_ao, ACT_VEC, n_act, VIRT_VEC, n_virt);
//         transform_from_col_MO(H_AA, H_core, n_ao, ACT_VEC, n_act, ACT_VEC, n_act);
        transform_from_col_MO(H_CA, H_core, n_ao, COR_VEC, n_cor, ACT_VEC, n_act);
        transform_from_col_MO(H_CV, H_core, n_ao, COR_VEC, n_cor, VIRT_VEC, n_virt);
        
        for(int i=0;i<n_ao*n_ao;i++)H_core[i]=M->H_AO[i]+2*J[i]-K[i];
        transform_from_col_MO(H_AA, H_core, n_ao, ACT_VEC, n_act, ACT_VEC, n_act);
        
// //         D.AO_to_MO(H_core);
//         calc_2el_MO_INTS_AAAG_RI(AAAG_INTS, M->MO_VEC, n_cor, n_act, n_ao);
        
        
        
        
        printf_timer("before PT");
        
        
        
        fprintf(out_stream,"Start calculation of PT2 correction\n");
        
        RI_data R;
        R.set_par(M,n_cor, n_act, n_virt);
        R.MO_calc(D.VEC[0], n_ao);
        fprintf(out_stream,"\n");
        fprintf(out_stream,"\n");
        printf_timer("RI orbital transformation");
        fprintf(out_stream,"\n");
        fflush(out_stream);
        
        
        //PT
        PT_tensors T;
        T.set_par(&R, eps, n_cor, n_act, n_virt, H_AV, H_CA, H_CV, cdas->edshift);
        
        //calculation of IP/EA Fockian matrix
//         double zero=0;
        D.aldet.simple_import_data(act_INTS, act_INTS, H_AA, 0);
        double * H_IP = new double[(n_cor+n_act)*(n_cor+n_act)];
        double * U_IP = new double[(n_cor+n_act)*(n_cor+n_act)];
        double * H_EA = new double[(n_virt+n_act)*(n_virt+n_act)];
        double * U_EA = new double[(n_virt+n_act)*(n_virt+n_act)];
        
        
//         D.aldet.calc_IPEA_single_CA(U_IP,H_IP,U_EA,H_EA,n_cor,n_mo,n_ao,0,cdas->cas->w_state,H_core,AAAG_INTS);
        
        
//         F_IPEA(T.e_IP,T.IP_U,T.IP_Um,T.IP_Up,-1, &(D.aldet),0,cdas->cas->w_state);//IP
        
        if(cdas->IPEA){
            T.IPEA(&(D.aldet), 0,cdas->cas->w_state);
            T.E2_calc_IPEA();
        }
        else if(cdas->MPPT){
            T.MPPT(&(D.aldet), 0,cdas->cas->w_state);
            T.E2_calc_EE();
        }
        else{
            T.E2_calc_EE();
        }
        
        printf_timer("PT tensors calculation");
        fprintf(out_stream,"_______________________________________________________________________\n\n\n");
        
//                 fPrintMatr(out_stream,T.RF_PH,n_act,n_act,1);
        if(CAS.CI->as_aldet() == nullptr){
            fprintf(out_stream,"ERROR: CDAS-PT (rel) requires the determinant CI backend (cisolver=aldet)\n");
            exit(EXIT_FAILURE);
        }
        CAS.CI->as_aldet()->PT2_import_data(T.RF_P3_JK,
                                T.RF_P3_AB,
                                T.RF_PV_JK,
                                T.RF_PV_AB,
                                T.RF_PH,
                                T.RF_PS);
         CAS.tensors_recalc(0);
//         CAS.CI_calc(1,1);
        
        double E0 = CAS.CI->as_aldet()->E_core;//+CAS.CI->as_aldet()->T0;
//         for(int i=0;i<n_act*n_act;i++)H_AA[i]+=T.RF_PH[i];
        
        double * so_x_vv = new double[n_act*n_act];
        double * so_y_vv = new double[n_act*n_act];
        double * so_z_vv = new double[n_act*n_act];
        double * so_x = new double[n_ao*n_ao];set_zero_matr(so_x, n_ao*n_ao);
        double * so_y = new double[n_ao*n_ao];set_zero_matr(so_y, n_ao*n_ao);
        double * so_z = new double[n_ao*n_ao];set_zero_matr(so_z, n_ao*n_ao);

        set_zero_matr(so_x_vv,n_act*n_act);
        set_zero_matr(so_y_vv,n_act*n_act);
        set_zero_matr(so_z_vv,n_act*n_act);
        grpp_engine.calc_SO_AO(so_x,so_y,so_z);
//         print_matrix_lower_triangle("NOPT_so_x.txt", n_ao, so_x);
//         print_matrix_lower_triangle("NOPT_so_y.txt", n_ao, so_y);
//         print_matrix_lower_triangle("NOPT_so_z.txt", n_ao, so_z);
    
        transform_from_col_MO(so_x_vv,so_x, n_ao, ACT_VEC, n_act, ACT_VEC, n_act);
        transform_from_col_MO(so_y_vv,so_y, n_ao, ACT_VEC, n_act, ACT_VEC, n_act);
        transform_from_col_MO(so_z_vv,so_z, n_ao, ACT_VEC, n_act, ACT_VEC, n_act);
        
//         printf("so_x:\n");PrintMatr(so_x_vv,n_act,n_act,0);
//         printf("so_y:\n");PrintMatr(so_y_vv,n_act,n_act,0);
//         printf("so_z:\n");PrintMatr(so_z_vv,n_act,n_act,0);
        
        
        
//         int dim;
//         
//         dim=2*n_act;
//         
//         double* H_ci2   = new double [dim*dim];
//         double* H_ci2_i = new double [dim*dim];
//         double * E_ci = new double[dim];
//         
//         for(int i=0; i<dim;i++)
//             H_ci2[i*(dim+1)]+=E0;
//         
//         
// 
//         H_add_SO(H_ci2, H_ci2_i, H_AA, so_x_vv, so_y_vv, so_z_vv, n_act);
//         
        
        
        aldet_rel_data CI_rel;
        
        CI_rel.get_dim(CAS.CI->as_aldet()->n_act,CAS.CI->as_aldet()->na, CAS.CI->as_aldet()->nb, 2, CAS.CI->as_aldet()->print_number);
        CI_rel.init_zero_vec(-1,0);
        
        double * AA = new double[n_act*n_act*n_act*n_act];
        double * AB = new double[n_act*n_act*n_act*n_act];
        
        for(int i=0; i<n_act*n_act*n_act*n_act;i++)AA[i]=CAS.CI->as_aldet()->act_INTS_AA[i];//+T.RF_PV_JK[i];
        for(int i=0; i<n_act*n_act*n_act*n_act;i++)AB[i]=CAS.CI->as_aldet()->act_INTS_AB[i];//+T.RF_PV_AB[i];
        
        CI_rel.simple_import_data(AA,
                                  AB,
                                  CAS.CI->as_aldet()->F_act_A,
                                  E0);
        CI_rel.SO_update(so_x_vv,so_y_vv,so_z_vv);
        
        
        CI_rel.gen_ext_ind();
        CI_rel.PT2_import_data(T.RF_P3_JK,
                               T.RF_P3_AB,
                               T.RF_PV_JK,
                               T.RF_PV_AB,
                               T.RF_PH,
                               T.RF_PS);
        CI_rel.PT_update();
        
//         double * H_ci3   = new double[CI_rel.Nd*CI_rel.Nd];
//         double * H_ci3_i = new double[CI_rel.Nd*CI_rel.Nd];
//         set_zero_matr(H_ci3_i, CI_rel.Nd*CI_rel.Nd);
//         double * E_ci3 = new double[CI_rel.Nd];
//         CI_rel.H_full_calc  (H_ci3);
//         CI_rel.H_full_calc_i(H_ci3_i);

//         lapack_herm_diag(H_ci3,H_ci3_i,E_ci3, CI_rel.Nd);
        
        CI_rel.H_full_calc_and_diag(0);
        
        CI_rel.print_states(0,CI_rel.Nd,1);
        
//         for(int i=0;i<CI_rel.Nd;i++)
//         for(int j=0;j<n_s      ;j++)
//             CI_rel.coef_r[0][i*n_s+j]=H_ci3[j*CI_rel.Nd+i];
        
//         CI_rel.print_states(0,n_s,1);
//     }
    fprintf(out_stream,"\n\nCDAS-PT2 Energy summary:\n");
    PrintEnergy(CI_rel.E_states[0],CI_rel.Nd,1);
    // exit(0);
    // CAS.Prop_calc();
//     getchar();
    fprintf(out_stream,"\n");
    // if(write_ci){
        // fprintf(out_stream,"Writing CDAS-PT2 WaveFunctions:\n");
        // sprintf(name,"%s_CDAS.ci\0",job_name);
        // CAS.CI->as_aldet()->write_civec(0, name);
        // fprintf(out_stream,"data file         : %s\n",name);
    // }
    

    // CAS.print_properties("CDAS-PT(0)");    
    
//     double * d_x = CAS.Prop_value                  ;
//     double * d_y = CAS.Prop_value+CAS.n_s*CAS.n_s  ;
//     double * d_z = CAS.Prop_value+CAS.n_s*CAS.n_s*2;
//     
//     
//     set_zero_matr(d_x1,n_s*n_s);
//     set_zero_matr(d_y1,n_s*n_s);
//     set_zero_matr(d_z1,n_s*n_s);
// #ifdef _CDAS_FO
//     {
//         fprintf(out_stream,"PT1 dipole moment - d(1):\n\n");
//         fprintf(out_stream,"d_x(1)\n");
//         T.P1_calc(d_x1, CAS.CI->as_aldet(), d_x_AV, d_x_CA, d_x_CV, 1, 1, 1);
//         fprintf(out_stream,"\n");
//         fprintf(out_stream,"d_y(1)\n");
//         T.P1_calc(d_y1, CAS.CI->as_aldet(), d_y_AV, d_y_CA, d_y_CV, 1, 1, 1);
//         fprintf(out_stream,"\n");
//         fprintf(out_stream,"d_z(1)\n");
//         T.P1_calc(d_z1, CAS.CI->as_aldet(), d_z_AV, d_z_CA, d_z_CV, 1, 1, 1);
//         
//         symmetrization_with_scaling(d_x1,n_s,2.0);
//         symmetrization_with_scaling(d_y1,n_s,2.0);
//         symmetrization_with_scaling(d_z1,n_s,2.0);
//         
//         for(int i_s=0;i_s<n_s*n_s;i_s++){
//             d_x[i_s] += d_x1[i_s];
//             d_y[i_s] += d_y1[i_s];
//             d_z[i_s] += d_z1[i_s];
//         }
//         
//         fprintf(out_stream,"\n");
//         fprintf(out_stream,"\n\nDipole CDAS(0+1):\n");
//         PrintDipole(d_x,d_y,d_z,n_s);
//     }
// 
// #endif
//         
//     fprintf(out_stream,"_______________________________________________________________________\n\n\n");
// //         fprintf(out_stream,"%e %e %e\n",M->Dx_nuc,M->Dy_nuc,M->Dz_nuc);
// //         exit(0);
// #ifdef _TRrCAMM
//     double *A_DM = new double[n_ao*n_ao];
//     double q;
//     for(int i_a=0;i_a<M->n_atoms;i_a++){
//         
//         gen_atomic_DM(A_DM, M->S_MO, i_a, M->basis[0].ao_center, n_ao);
//         for(int i=0;i<n_ao*n_ao;i++)A_DM[i]=-A_DM[i];
//         fprintf(out_stream,"TrC:\n");
//         for(int i_s=0;i_s<n_s;i_s++){
//             for(int j_s=0;j_s<n_s;j_s++){
//                 D.transform_1el_DM(i_s*n_s+j_s);
//                 q=E_1el_calc(A_DM, D.DM_T, D.n_ao[0], D.n_ao[1])*D.p_SVD+D.ACT_DET_A[i_s*n_s+j_s]*M->nucl_charges[i_a];
//                 
//                 fprintf(out_stream,"% .5e ", q);
//             }
//             fprintf(out_stream,"\n");
//         }
//         getchar();
//         
//         double * q1 = new double[n_s*n_s];
//         set_zero_matr(q1, n_s*n_s);
//         
//         double * q_AV;
//         double * q_CA;
//         double * q_CV;
//         q_AV = new double[n_act*n_virt];
//         q_CA = new double[n_cor*n_act];
//         q_CV = new double[n_cor*n_virt];
//         transform_from_col_MO(q_AV, A_DM, n_ao, ACT_VEC, n_act, VIRT_VEC, n_virt);
//         transform_from_col_MO(q_CA, A_DM, n_ao, COR_VEC, n_cor, ACT_VEC, n_act);
//         transform_from_col_MO(q_CV, A_DM, n_ao, COR_VEC, n_cor, VIRT_VEC, n_virt);
//         
//     
//         T.P1_calc(q1, q_AV, q_CA, q_CV, 1, 1, 1);
//         
//         
//         symmetrization_with_scaling(q1,n_s,2.0);
//         fPrintMatr(out_stream,q1,n_s,n_s,1);
//         
//         
//     }
// #endif
        
    printf_timer("CDAS-PT2");
    
    
    D.clear();
    
    delete[] H_CV;
    delete[] H_CA;
    delete[] H_AA;
    delete[] H_AV;
    delete[] F;
    delete[] COR_VEC;
    delete[] ACT_VEC;
    delete[] OCC_VEC;
    delete[] VIRT_VEC;
    
    delete[] d_x1  ;
    delete[] d_y1  ;
    delete[] d_z1  ;
    delete[] d_x_CV;
    delete[] d_y_CV;
    delete[] d_z_CV;
    delete[] d_x_AV;
    delete[] d_y_AV;
    delete[] d_z_AV;
    delete[] d_x_CA;
    delete[] d_y_CA;
    delete[] d_z_CA;
    delete[] eps;

//     delete[] eps_EA;
    
    delete[] U;
    delete[] name;
    
    delete[] J       ;
    delete[] K       ;
    delete[] act_INTS ;
#endif
#ifndef _USE_GRPP
    printf("CDAS_PT2_rel can not be used without GRPP\n");
    exit(0);
#endif
    return 0;
 
}
