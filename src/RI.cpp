
# include "blas_link.h"
# include "common_vars.h"
# include "RI.h"
# include "libint_link.h"
# include "matr.h"
# include "timer.h"

using libint2::Engine;
using libint2::Operator;

//A - atomic
//a - active
//C - core

double * __restrict__ RI_AA  =nullptr;
double * __restrict__ RI_CA =nullptr;
double * __restrict__ RI_CA2=nullptr;
double * __restrict__ AUX_Vm05=nullptr;
double * __restrict__ RI_aA  =nullptr;
double * __restrict__ RI_aA2  =nullptr;
double * __restrict__ RI_aa  =nullptr;
double * __restrict__ RI_pp  =nullptr;
double * __restrict__ TMP  =nullptr; /// change it

double * RI_J  = nullptr;
double * RI_J2 = nullptr;
double * aaag2 = nullptr;



long aux_n_ao;

int gen_RI_AA(molecule * A){
    
    if(RI==0){
        fprintf(stdout,"ERROR: NOPT tried to run gen_RI_AA() with RI=0\n");
        exit(1);
    }
    
    if(RI_AA!=nullptr)
        return 0;
    
    long n_ao = A->n_ao;
    aux_n_ao = 0;
    for(auto &s: A->ri_s)
        aux_n_ao+=s.size();
    
    double * AUX_V;
    double * AUX_Vp05;
    AUX_V   =new double[aux_n_ao*aux_n_ao];
    AUX_Vp05=new double[aux_n_ao*aux_n_ao];
    if(AUX_Vm05==nullptr)
        AUX_Vm05=new double[aux_n_ao*aux_n_ao];
    
    set_zero_matr(AUX_V,aux_n_ao*aux_n_ao);
    RI_V_calc(A->ri_s, aux_n_ao, AUX_V);
    
    S05_calc(AUX_V, AUX_Vm05, AUX_Vp05, aux_n_ao);

    delete[]AUX_V;
    delete[]AUX_Vp05;
//     printf_timer("dealloc");
    
//     printf_timer("come to RI");
    if(RI_AA==nullptr){
        int n_c  = A->n_el_calc/2;
//         fprintf(out_stream,"n_c=%d\nn_a=%d\n",n_c,n_ao);
        RI_CA  = new double[long(n_c*n_ao*aux_n_ao)];
        RI_CA2 = new double[long(n_c*n_ao*aux_n_ao)];
        RI_AA =new double[long(n_ao*n_ao*aux_n_ao)];
    }
    calc_2el_AO_INTS_RI(A->s, n_ao, A->ri_s, aux_n_ao, RI_AA);
    
    fprintf(out_stream,"\n\n");
    printf_timer("RI integrals calculation");
    
    return 0;
}

int regen_RI_AA(molecule * A){
    
    if(RI_AA!=nullptr){
        delete[] RI_AA ;RI_AA =nullptr;
        delete[] RI_CA ;RI_CA =nullptr;
        delete[] RI_CA2;RI_CA2=nullptr;
    }

    if(AUX_Vm05!=nullptr){delete[] AUX_Vm05; AUX_Vm05=nullptr;}
    if(RI_aA   !=nullptr){delete[] RI_aA   ; RI_aA   =nullptr;}
    if(RI_aA2  !=nullptr){delete[] RI_aA2  ; RI_aA2  =nullptr;}
    if(RI_aa   !=nullptr){delete[] RI_aa   ; RI_aa   =nullptr;}
    if(RI_J    !=nullptr){delete[] RI_J    ; RI_J    =nullptr;}
    if(RI_J2   !=nullptr){delete[] RI_J2   ; RI_J2   =nullptr;}
    if(aaag2   !=nullptr){delete[] aaag2   ; aaag2   =nullptr;}
    if(RI_pp   !=nullptr){delete[] RI_pp   ; RI_pp   =nullptr;}
    
    gen_RI_AA(A);
    
    return 0;
}


int RI_core_realloc(long n_c, long n_ao){
    
    if(RI_CA !=nullptr) delete[] RI_CA ;
    if(RI_CA2!=nullptr) delete[] RI_CA2;
    
    RI_CA  = new double[long(n_c*n_ao*aux_n_ao)];
    RI_CA2 = new double[long(n_c*n_ao*aux_n_ao)];
    
    return 0;
}

int DM_to_F_transform_RI(double * J, double * K, double * DM, double * MO, long n_c, long n_ao){
    
    if(RI_J  == nullptr) RI_J  = new double[aux_n_ao];
    if(RI_J2 == nullptr) RI_J2 = new double[aux_n_ao];
    
    set_zero_matr(RI_J, aux_n_ao);
    
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                1,aux_n_ao,n_ao*n_ao,1.0,
                DM ,n_ao*n_ao,
                RI_AA,aux_n_ao,0.0,
                RI_J,aux_n_ao);
    
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   aux_n_ao,1,aux_n_ao,1.0,
                   AUX_Vm05 ,aux_n_ao,
                   RI_J,1,0.0,
                   RI_J2,1);

    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   aux_n_ao,1,aux_n_ao,1.0,
                   AUX_Vm05 ,aux_n_ao,
                   RI_J2,1,0.0,
                   RI_J,1);

    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   n_ao*n_ao,1,aux_n_ao,1.0,
                   RI_AA ,aux_n_ao,
                   RI_J,1,0.0,
                   J,1);
    
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   n_c,aux_n_ao*n_ao, n_ao,1.0,
                   MO ,n_ao,
                   RI_AA,aux_n_ao*n_ao,0.0,
                   RI_CA,aux_n_ao*n_ao);
    
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   n_c*n_ao,aux_n_ao, aux_n_ao,1.0,
                   RI_CA ,aux_n_ao,
                   AUX_Vm05 ,aux_n_ao,0.0,
                   RI_CA2,aux_n_ao);
    
#pragma omp parallel for
    for(int b=0; b<    n_ao; b++)
    for(int a=0; a<    n_c ; a++)
    for(int k=0; k<aux_n_ao; k++)
    {
        RI_CA[(b*n_c+a)*aux_n_ao +k]=RI_CA2[(a*n_ao+b)*aux_n_ao +k];
    }

    
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                   n_ao,n_ao,aux_n_ao*n_c,1.0,
                   RI_CA,aux_n_ao*n_c,
                   RI_CA,aux_n_ao*n_c,0.0,
                   K,n_ao);    //K[a,b]=sum(k,i)R2[a,k,i]*T2[b,k,i]

    
    
    return 1;
}

int calc_2el_MO_INTS_RI(long n_ao, double * DM_D, double * J, double * K, double * MO_INTS, 
                        double * MO_VEC, double * MO_VEC1, double * MO_VEC2, 
                        long n_c, long n_a, long t){
    
    //V
    if(RI_aa==nullptr){
        RI_aA = new double[long(n_a*n_ao*aux_n_ao)];
        RI_aA2= new double[long(n_a*n_ao*aux_n_ao)];
        RI_aa = new double[long(n_a*n_a *aux_n_ao)];
    }
    
    CBLAS_TRANSPOSE T = CblasNoTrans;
    int ld_MO = n_ao;
    if(t){
        T = CblasTrans;
        ld_MO = n_a;
    }
    
    nopt_par_dgemm(CblasRowMajor,T,CblasNoTrans,
                   n_a,aux_n_ao*n_ao, n_ao,1.0,
                   MO_VEC1, ld_MO,
                   RI_AA,aux_n_ao*n_ao,0.0,
                   RI_aA,aux_n_ao*n_ao);
    
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   n_a*n_ao,aux_n_ao, aux_n_ao,1.0,
                   RI_aA ,aux_n_ao,
                   AUX_Vm05 ,aux_n_ao,0.0,
                   RI_aA2,aux_n_ao);
    
#pragma omp parallel for
    for(int b=0; b<    n_ao; b++)
    for(int a=0; a<    n_a ; a++)
    for(int k=0; k<aux_n_ao; k++)
    {
        RI_aA[(b*n_a+a)*aux_n_ao +k]=RI_aA2[(a*n_ao+b)*aux_n_ao +k];
    }
    
    nopt_par_dgemm(CblasRowMajor,T,CblasNoTrans,
                   n_a,aux_n_ao*n_a, n_ao,1.0,
                   MO_VEC2, ld_MO,
                   RI_aA,aux_n_ao*n_a,0.0,
                   RI_aa,aux_n_ao*n_a);
    
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                   n_a*n_a,n_a*n_a,aux_n_ao,1.0,
                   RI_aa,aux_n_ao,
                   RI_aa,aux_n_ao,0.0,
                   MO_INTS,n_a*n_a);    
    
//     PrintMatr(MO_INTS, n_a*n_a,n_a*n_a,1);
    // J
    if(J==NULL)
        return 0;
    if(K==NULL)
        return 0;
    
    if(RI_J  == nullptr) RI_J  = new double[aux_n_ao];
    if(RI_J2 == nullptr) RI_J2 = new double[aux_n_ao];
    
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                1,aux_n_ao,n_ao*n_ao,1.0,
                DM_D ,n_ao*n_ao,
                RI_AA,aux_n_ao,0.0,
                RI_J,aux_n_ao);
    
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   aux_n_ao,1,aux_n_ao,1.0,
                   AUX_Vm05 ,aux_n_ao,
                   RI_J,1,0.0,
                   RI_J2,1);

    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   aux_n_ao,1,aux_n_ao,1.0,
                   AUX_Vm05 ,aux_n_ao,
                   RI_J2,1,0.0,
                   RI_J,1);

    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   n_ao*n_ao,1,aux_n_ao,1.0,
                   RI_AA ,aux_n_ao,
                   RI_J,1,0.0,
                   J,1);
    
    // K
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   n_c,aux_n_ao*n_ao, n_ao,1.0,
                   MO_VEC ,n_ao,
                   RI_AA,aux_n_ao*n_ao,0.0,
                   RI_CA,aux_n_ao*n_ao);
    
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   n_c*n_ao,aux_n_ao, aux_n_ao,1.0,
                   RI_CA ,aux_n_ao,
                   AUX_Vm05 ,aux_n_ao,0.0,
                   RI_CA2,aux_n_ao);
    
#pragma omp parallel for
    for(int b=0; b<    n_ao; b++)
    for(int a=0; a<    n_c ; a++)
    for(int k=0; k<aux_n_ao; k++)
    {
        RI_CA[(b*n_c+a)*aux_n_ao +k]=RI_CA2[(a*n_ao+b)*aux_n_ao +k];
    }

    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                   n_ao,n_ao,aux_n_ao*n_c,1.0,
                   RI_CA,aux_n_ao*n_c,
                   RI_CA,aux_n_ao*n_c,0.0,
                   K,n_ao);    //K[a,b]=sum(k,i)R2[a,k,i]*T2[b,k,i]
    
    
    return 0;
    
}

int calc_2el_MO_INTS_AAAG_RI(double * MO_INTS, double * MO_VEC, long n_c, long n_a, long n_ao){
    
    if(RI_aa==nullptr){
        RI_aA = new double[long(n_a*n_ao*aux_n_ao)];
        RI_aA2= new double[long(n_a*n_ao*aux_n_ao)];
        RI_aa = new double[long(n_a*n_a *aux_n_ao)];
    }
    if(aaag2 == nullptr)
        aaag2 = new double[long(n_a*n_a*n_a*n_ao)];
    
    //(AB|K)->(tB|K)
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   n_a,aux_n_ao*n_ao, n_ao,1.0,
                   MO_VEC + n_c*n_ao, n_ao,
                   RI_AA,aux_n_ao*n_ao,0.0,
                   RI_aA,aux_n_ao*n_ao);
    
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   n_a*n_ao,aux_n_ao, aux_n_ao,1.0,
                   RI_aA ,aux_n_ao,
                   AUX_Vm05 ,aux_n_ao,0.0,
                   RI_aA2,aux_n_ao);
    //(tA|K)->(At|K)
#pragma omp parallel for
    for(int b=0; b<    n_ao; b++)
    for(int a=0; a<    n_a ; a++)
    for(int k=0; k<aux_n_ao; k++)
    {
        RI_aA[(b*n_a+a)*aux_n_ao +k]=RI_aA2[(a*n_ao+b)*aux_n_ao +k];
    }
    //(Au|K)->(tu|K)
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   n_a,aux_n_ao*n_a, n_ao,1.0,
                   MO_VEC + n_c*n_ao, n_ao,
                   RI_aA,aux_n_ao*n_a,0.0,
                   RI_aa,aux_n_ao*n_a);
    
    //(tu|K)*(vA|K)->(tu|vA)
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                   n_a*n_a,n_a*n_ao,aux_n_ao,1.0,
                   RI_aa,aux_n_ao,
                   RI_aA2,aux_n_ao,0.0,
                   aaag2,n_a*n_ao);    
    //(tu|vA)->(tu|vx)
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                   n_a*n_a*n_a,n_ao,n_ao,1.0,
                   aaag2,n_ao,
                   MO_VEC,n_ao,0.0,
                   MO_INTS,n_ao);    
    
    
    return 0;
}

int calc_SOSCF_hess_INTS_RI (double * aa_J_ints, 
                             double * aa_K_ints,
                             double * cv_J_ints,
                             double * cv_K_ints, double * MO_VEC, long n_c, long n_a, long n_ao){
    if(RI_aa==nullptr){
        RI_aA = new double[long(n_a*n_ao*aux_n_ao)];
        RI_aA2= new double[long(n_a*n_ao*aux_n_ao)];
        RI_aa = new double[long(n_a*n_a *aux_n_ao)];
    }
    if(RI_pp==nullptr)RI_pp = new double[long(n_ao *aux_n_ao*2)];
    double * RI_pp_2 = RI_pp+long(n_ao *aux_n_ao);
        
    //(AB|K)->(tB|K)
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   n_a,aux_n_ao*n_ao, n_ao,1.0,
                   MO_VEC + n_c*n_ao, n_ao,
                   RI_AA,aux_n_ao*n_ao,0.0,
                   RI_aA,aux_n_ao*n_ao);
    
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   n_a*n_ao,aux_n_ao, aux_n_ao,1.0,
                   RI_aA ,aux_n_ao,
                   AUX_Vm05 ,aux_n_ao,0.0,
                   RI_aA2,aux_n_ao);
    //(tA|K)->(At|K)
#pragma omp parallel for
    for(int b=0; b<    n_ao; b++)
    for(int a=0; a<    n_a ; a++)
    for(int k=0; k<aux_n_ao; k++)
    {
        RI_aA[(b*n_a+a)*aux_n_ao +k]=RI_aA2[(a*n_ao+b)*aux_n_ao +k];
    }
    //(Au|K)->(tu|K)
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   n_a,aux_n_ao*n_a, n_ao,1.0,
                   MO_VEC + n_c*n_ao, n_ao,
                   RI_aA,aux_n_ao*n_a,0.0,
                   RI_aa,aux_n_ao*n_a);
    
    set_zero_matr(RI_pp,n_ao *aux_n_ao*2);
    if(TMP==nullptr) TMP = new double[long(n_ao*n_ao*n_ao)];
    for(int p=0; p<n_ao    ; p++){//fprintf(stderr,"RI_pp %d/%d\r",p,n_ao); 
        for(int i=0; i<n_ao    ; i++)
        for(int j=0; j<n_ao    ; j++)
            TMP[(p*n_ao+i)*n_ao+j]=MO_VEC[p*n_ao+i]*MO_VEC[p*n_ao+j];
    }
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   n_ao,aux_n_ao, n_ao*n_ao,1.0,
                   TMP, n_ao*n_ao,
                   RI_AA,aux_n_ao,0.0,
                   RI_pp,aux_n_ao);
//     for(int i=0; i<n_ao    ; i++)
//     for(int j=0; j<n_ao    ; j++)
//     for(int k=0; k<aux_n_ao; k++)
//             RI_pp[p*aux_n_ao+k]+=TMP[i*n_ao+j]*RI_AA[(i*n_ao+j)*aux_n_ao+k];
//     }
    for(int p=0; p<n_ao    ; p++)
    for(int k=0; k<aux_n_ao; k++)
    for(int l=0; l<aux_n_ao; l++)
        RI_pp_2[p*aux_n_ao+k]+=RI_pp[p*aux_n_ao+l]*AUX_Vm05[k*aux_n_ao+l];
    
    set_zero_matr(aa_J_ints, n_ao*n_a*n_a);
    for(int p=0; p<n_ao    ; p++)
    for(int i=0; i<n_a*n_a ; i++)
    for(int k=0; k<aux_n_ao; k++)
        aa_J_ints[p*n_a*n_a+i]+=RI_pp_2[p*aux_n_ao+k]*RI_aa[i*aux_n_ao+k];
    
    
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   n_ao,aux_n_ao*n_a, n_ao,1.0,
                   MO_VEC, n_ao,
                   RI_aA ,aux_n_ao*n_a,0.0,
                   RI_aA2,aux_n_ao*n_a);
    
    set_zero_matr(aa_K_ints, n_ao*n_a*n_a);
    for(int p=0; p<n_ao    ; p++){//fprintf(stderr,"aa_K %d/%d\r",p,n_ao); 
    for(int i=0; i<n_a     ; i++)
    for(int j=0; j<n_a     ; j++)
    for(int k=0; k<aux_n_ao; k++)
        aa_K_ints[p*n_a*n_a+i*n_a+j]+=RI_aA2[(p*n_a+i)*aux_n_ao+k]*RI_aA2[(p*n_a+j)*aux_n_ao+k];
    }
    
    
    
    return 0;
}



RI_data::RI_data(){
    
    VC_RI_M = nullptr;
    CA_RI_M = nullptr;
    CC_RI_M = nullptr;
    VA_RI_M = nullptr;
    AA_RI_M = nullptr;
    
    COR_CVEC = nullptr;
    ACT_CVEC = nullptr;
    VAC_CVEC = nullptr;
    
}

int RI_data::set_par(molecule * A, long ext_n_cor, long ext_n_act, long ext_n_vac){
    
    n_cor = ext_n_cor;
    n_act = ext_n_act;
    n_vac = ext_n_vac;
    n_ao  = n_cor + n_act + n_vac;
    
    s =A->s;
    
    aux_shells = A->ri_s;
    aux_n_ao = 0;
    for(auto &s: A->ri_s)
        aux_n_ao+=s.size();
        
        
    return 0;
}

int RI_data::MO_calc(double * VEC, long ext_n_ao){
    
    
    n_ao=ext_n_ao;
    if(RI){
        double * COR = VEC;
        double * ACT = VEC+n_cor*n_ao;
        double * VAC = VEC+(n_cor+n_act)*n_ao;
        double * BUF = new double[long(n_cor*n_ao*aux_n_ao)];
        double * BUF2= new double[long(n_cor*n_ao*aux_n_ao)];
        
        nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                       n_cor,aux_n_ao*n_ao, n_ao,1.0,
                       COR, n_ao,
                       RI_AA,aux_n_ao*n_ao,0.0,
                       BUF,aux_n_ao*n_ao);
        
        nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                       n_cor*n_ao,aux_n_ao, aux_n_ao,1.0,
                       BUF ,aux_n_ao,
                       AUX_Vm05 ,aux_n_ao,0.0,
                       BUF2,aux_n_ao);
        
#pragma omp parallel for
        for(int b=0; b<    n_ao ; b++)
        for(int a=0; a<    n_cor; a++)
        for(int k=0; k<aux_n_ao ; k++)
        {
            BUF[(b*n_cor+a)*aux_n_ao +k]=BUF2[(a*n_ao+b)*aux_n_ao +k];
        }
        
        CC_RI_M=new double[long(n_cor*n_cor*aux_n_ao)];
        nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                       n_cor,aux_n_ao*n_cor, n_ao,1.0,
                       COR, n_ao,
                       BUF,aux_n_ao*n_cor,0.0,
                       CC_RI_M,aux_n_ao*n_cor);
        
        
        VC_RI_M=new double[long(n_cor*n_vac*aux_n_ao)];
        nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                       n_vac,aux_n_ao*n_cor, n_ao,1.0,
                       VAC, n_ao,
                       BUF,aux_n_ao*n_cor,0.0,
                       VC_RI_M,aux_n_ao*n_cor);
        
        delete[] BUF ;
        delete[] BUF2;
        
        BUF = new double[long(n_act*n_ao*aux_n_ao)];
        BUF2= new double[long(n_act*n_ao*aux_n_ao)];
        
        nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                       n_act,aux_n_ao*n_ao, n_ao,1.0,
                       ACT, n_ao,
                       RI_AA,aux_n_ao*n_ao,0.0,
                       BUF,aux_n_ao*n_ao);
        
        nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                       n_act*n_ao,aux_n_ao, aux_n_ao,1.0,
                       BUF ,aux_n_ao,
                       AUX_Vm05 ,aux_n_ao,0.0,
                       BUF2,aux_n_ao);
        
#pragma omp parallel for
        for(int b=0; b<    n_ao ; b++)
        for(int a=0; a<    n_act; a++)
        for(int k=0; k<aux_n_ao ; k++)
        {
            BUF[(b*n_act+a)*aux_n_ao +k]=BUF2[(a*n_ao+b)*aux_n_ao +k];
        }
        
        AA_RI_M=new double[long(n_act*n_act*aux_n_ao)];
        nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                       n_act,aux_n_ao*n_act, n_ao,1.0,
                       ACT, n_ao,
                       BUF,aux_n_ao*n_act,0.0,
                       AA_RI_M,aux_n_ao*n_act);
        
        CA_RI_M=new double[long(n_act*n_cor*aux_n_ao)];
        nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                       n_cor,aux_n_ao*n_act, n_ao,1.0,
                       COR, n_ao,
                       BUF,aux_n_ao*n_act,0.0,
                       CA_RI_M,aux_n_ao*n_act);
        
        VA_RI_M=new double[long(n_act*n_vac*aux_n_ao)];
        nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                       n_vac,aux_n_ao*n_act, n_ao,1.0,
                       VAC, n_ao,
                       BUF,aux_n_ao*n_act,0.0,
                       VA_RI_M,aux_n_ao*n_act);
        
        delete[] BUF ;
        delete[] BUF2;
    }
    else{
        COR_CVEC = new double[n_cor*n_ao];
        ACT_CVEC = new double[n_act*n_ao];
        VAC_CVEC = new double[n_vac*n_ao];
        
        for(int i=0; i<n_ao ;i++)
        for(int j=0; j<n_cor;j++)
            COR_CVEC[i*n_cor+j]=VEC[j*n_ao+i];
        
        
        for(int i=0; i<n_ao ;i++)
        for(int j=0; j<n_act;j++)
            ACT_CVEC[i*n_act+j]=VEC[(j+n_cor)*n_ao+i];
        
        for(int i=0; i<n_ao ;i++)
        for(int j=0; j<n_vac;j++)
            VAC_CVEC[i*n_vac+j]=VEC[(j+n_cor+n_act)*n_ao+i];
        
        
    }
    
    return 0;
}


int RI_data::VAAA_calc(double * VAAA){
    
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                n_vac*n_act,n_act*n_act,aux_n_ao,1.0,
                VA_RI_M,aux_n_ao,
                AA_RI_M,aux_n_ao,0.0,
                VAAA,n_act*n_act);
    
    return 0;
}

int RI_data::CAAA_calc(double * CAAA){
    
    
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                n_cor*n_act,n_act*n_act,aux_n_ao,1.0,
                CA_RI_M,aux_n_ao,
                AA_RI_M,aux_n_ao,0.0,
                CAAA,n_act*n_act);
        
    return 0;
}

int RI_data::VACC_calc(double * VACC_J, double * VACC_K){
    
#pragma omp parallel for
    for(int a=0; a<n_vac; a++)
    for(int t=0; t<n_act; t++)
    for(int j=0; j<n_cor; j++){
        VACC_J[(a*n_act+t)*n_cor+j] = cblas_ddot(aux_n_ao, VA_RI_M+(a*n_act+t)*aux_n_ao, 1, CC_RI_M+(j*n_cor+j)*aux_n_ao, 1);
        VACC_K[(a*n_act+t)*n_cor+j] = cblas_ddot(aux_n_ao, VC_RI_M+(a*n_cor+j)*aux_n_ao, 1, CA_RI_M+(j*n_act+t)*aux_n_ao, 1);
    }
    
    return 0;
}

int RI_data::CACC_calc(double * CACC_J, double * CACC_K){
    
#pragma omp parallel for
    for(int i=0; i<n_cor; i++)
    for(int t=0; t<n_act; t++)
    for(int j=0; j<n_cor; j++){
        CACC_J[(i*n_act+t)*n_cor+j] = cblas_ddot(aux_n_ao, CA_RI_M+(i*n_act+t)*aux_n_ao, 1, CC_RI_M+(j*n_cor+j)*aux_n_ao, 1);
        CACC_K[(i*n_act+t)*n_cor+j] = cblas_ddot(aux_n_ao, CC_RI_M+(i*n_cor+j)*aux_n_ao, 1, CA_RI_M+(j*n_act+t)*aux_n_ao, 1);
    }
    
    return 0;
}

RI_data::~RI_data(){
    
    if(VC_RI_M != nullptr) delete[] VC_RI_M;
    if(CA_RI_M != nullptr) delete[] CA_RI_M;
    if(CC_RI_M != nullptr) delete[] CC_RI_M;
    if(VA_RI_M != nullptr) delete[] VA_RI_M;
    if(AA_RI_M != nullptr) delete[] AA_RI_M;
    
    if(COR_CVEC != nullptr) delete[] COR_CVEC;
    if(ACT_CVEC != nullptr) delete[] ACT_CVEC;
    if(VAC_CVEC != nullptr) delete[] VAC_CVEC;
    
    
}

RI_complex_data::RI_complex_data(){
    
    VC_RI_M_r = nullptr;
    CA_RI_M_r = nullptr;
    CC_RI_M_r = nullptr;
    VA_RI_M_r = nullptr;
    AA_RI_M_r = nullptr;

    VC_RI_M_i = nullptr;
    CA_RI_M_i = nullptr;
    CC_RI_M_i = nullptr;
    VA_RI_M_i = nullptr;
    AA_RI_M_i = nullptr;


    
}

int RI_complex_data::set_par(molecule * A, long ext_n_cor, long ext_n_act, long ext_n_vac){
    
    n_cor = ext_n_cor;
    n_act = ext_n_act;
    n_vac = ext_n_vac;
    n_ao  = n_cor + n_act + n_vac;
    
//     aux_shells = A->ri_s;
    aux_n_ao = 0;
    for(auto &s: A->ri_s)
        aux_n_ao+=s.size();
        
        
    return 0;
}

int RI_complex_data::MO_calc(double * VEC_r, double * VEC_i, long ext_n_ao){
    
    
    n_ao=ext_n_ao;
    double * COR_r = VEC_r;
    double * ACT_r = VEC_r+2*n_cor*2*n_ao;
    double * VAC_r = VEC_r+(2*n_cor+2*n_act)*2*n_ao;
    double * COR_i = VEC_i;
    double * ACT_i = VEC_i+2*n_cor*2*n_ao;
    double * VAC_i = VEC_i+(2*n_cor+2*n_act)*2*n_ao;
    
    double * BUF_r = new double[long(4*n_cor*n_ao*aux_n_ao)];
    double * BUF_i = new double[long(4*n_cor*n_ao*aux_n_ao)];
    double * BUF2_r= new double[long(4*n_cor*n_ao*aux_n_ao)];
    double * BUF2_i= new double[long(4*n_cor*n_ao*aux_n_ao)];
    
    double * RI_AA_tmp =new double[long(4*n_ao*n_ao*aux_n_ao)];
    
    for(int i=0;i<    n_ao;i++)
    for(int j=0;j<    n_ao;j++)
    for(int k=0;k<aux_n_ao;k++)
        RI_AA_tmp[(i*2*n_ao+j)*aux_n_ao+k]=RI_AA[(i*n_ao+j)*aux_n_ao+k];
    for(int i=0;i<    n_ao;i++)
    for(int j=0;j<    n_ao;j++)
    for(int k=0;k<aux_n_ao;k++)
        RI_AA_tmp[((i+n_ao)*2*n_ao+j+n_ao)*aux_n_ao+k]=RI_AA[(i*n_ao+j)*aux_n_ao+k];
    
    
    
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   2*n_cor,aux_n_ao*2*n_ao, 2*n_ao,1.0,
                   COR_r, 2*n_ao,
                   RI_AA_tmp,aux_n_ao*2*n_ao,0.0,
                   BUF_r,aux_n_ao*2*n_ao);
    
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   2*n_cor,aux_n_ao*2*n_ao, 2*n_ao,1.0,
                   COR_i, 2*n_ao,
                   RI_AA_tmp,aux_n_ao*2*n_ao,0.0,
                   BUF_i,aux_n_ao*2*n_ao);
    
    
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   2*n_cor*2*n_ao,aux_n_ao, aux_n_ao,1.0,
                   BUF_r ,aux_n_ao,
                   AUX_Vm05 ,aux_n_ao,0.0,
                   BUF2_r,aux_n_ao);
    
    
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   2*n_cor*2*n_ao,aux_n_ao, aux_n_ao,1.0,
                   BUF_i ,aux_n_ao,
                   AUX_Vm05 ,aux_n_ao,0.0,
                   BUF2_i,aux_n_ao);
    
#pragma omp parallel for
    for(int b=0; b<  2*n_ao ; b++)
    for(int a=0; a<  2*n_cor; a++)
    for(int k=0; k<aux_n_ao ; k++)
    {
        BUF_r[(b*2*n_cor+a)*aux_n_ao +k]=BUF2_r[(a*2*n_ao+b)*aux_n_ao +k];
        BUF_i[(b*2*n_cor+a)*aux_n_ao +k]=BUF2_i[(a*2*n_ao+b)*aux_n_ao +k];
    }
    
    CC_RI_M_r=new double[long(2*n_cor*2*n_cor*aux_n_ao)];
    CC_RI_M_i=new double[long(2*n_cor*2*n_cor*aux_n_ao)];
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   2*n_cor,aux_n_ao*2*n_cor, 2*n_ao,1.0,
                   COR_r, 2*n_ao,
                   BUF_r,aux_n_ao*2*n_cor,0.0,
                   CC_RI_M_r,aux_n_ao*2*n_cor);
    
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   2*n_cor,aux_n_ao*2*n_cor, 2*n_ao,1.0,///complex conjugate
                   COR_i, 2*n_ao,
                   BUF_i,aux_n_ao*2*n_cor,1.0,
                   CC_RI_M_r,aux_n_ao*2*n_cor);
    
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   2*n_cor,aux_n_ao*2*n_cor, 2*n_ao,1.0,
                   COR_i, 2*n_ao,
                   BUF_r,aux_n_ao*2*n_cor,0.0,
                   CC_RI_M_i,aux_n_ao*2*n_cor);
    
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   2*n_cor,aux_n_ao*2*n_cor, 2*n_ao,-1.0,///complex conjugate
                   COR_i, 2*n_ao,
                   BUF_r,aux_n_ao*2*n_cor,1.0,
                   CC_RI_M_i,aux_n_ao*2*n_cor);
    
    
    
    VC_RI_M_r=new double[long(2*n_cor*2*n_vac*aux_n_ao)];
    VC_RI_M_i=new double[long(2*n_cor*2*n_vac*aux_n_ao)];
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   2*n_vac,aux_n_ao*2*n_cor, 2*n_ao,1.0,
                   VAC_r, 2*n_ao,
                   BUF_r,aux_n_ao*2*n_cor,0.0,
                   VC_RI_M_r,aux_n_ao*2*n_cor);
    
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   2*n_vac,aux_n_ao*2*n_cor, 2*n_ao,1.0,///complex conjugate
                   VAC_i, 2*n_ao,
                   BUF_i,aux_n_ao*2*n_cor,1.0,
                   VC_RI_M_r,aux_n_ao*2*n_cor);
    
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   2*n_vac,aux_n_ao*2*n_cor, 2*n_ao,1.0,
                   VAC_r, 2*n_ao,
                   BUF_i,aux_n_ao*2*n_cor,0.0,
                   VC_RI_M_i,aux_n_ao*2*n_cor);
    
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   2*n_vac,aux_n_ao*2*n_cor, 2*n_ao,-1.0,///complex conjugate
                   VAC_i, 2*n_ao,
                   BUF_r,aux_n_ao*2*n_cor,1.0,
                   VC_RI_M_i,aux_n_ao*2*n_cor);
    
    delete[] BUF_r ;
    delete[] BUF_i ;
    delete[] BUF2_r;
    delete[] BUF2_i;
    
    BUF_r = new double[long(4*n_act*n_ao*aux_n_ao)];
    BUF_i = new double[long(4*n_act*n_ao*aux_n_ao)];
    BUF2_r= new double[long(4*n_act*n_ao*aux_n_ao)];
    BUF2_i= new double[long(4*n_act*n_ao*aux_n_ao)];
    
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   2*n_act,aux_n_ao*2*n_ao, 2*n_ao,1.0,
                   ACT_r, 2*n_ao,
                   RI_AA_tmp,aux_n_ao*2*n_ao,0.0,
                   BUF_r,aux_n_ao*2*n_ao);
    
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   2*n_act,aux_n_ao*2*n_ao, 2*n_ao,1.0,
                   ACT_i, 2*n_ao,
                   RI_AA_tmp,aux_n_ao*2*n_ao,0.0,
                   BUF_i,aux_n_ao*2*n_ao);
    
    delete[] RI_AA_tmp;
    
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   2*n_act*2*n_ao,aux_n_ao, aux_n_ao,1.0,
                   BUF_r ,aux_n_ao,
                   AUX_Vm05 ,aux_n_ao,0.0,
                   BUF2_r,aux_n_ao);
    
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   2*n_act*2*n_ao,aux_n_ao, aux_n_ao,1.0,
                   BUF_i ,aux_n_ao,
                   AUX_Vm05 ,aux_n_ao,0.0,
                   BUF2_i,aux_n_ao);
    
#pragma omp parallel for
    for(int b=0; b<    2*n_ao ; b++)
    for(int a=0; a<    2*n_act; a++)
    for(int k=0; k<aux_n_ao ; k++)
    {
        BUF_r[(b*2*n_act+a)*aux_n_ao +k]=BUF2_r[(a*2*n_ao+b)*aux_n_ao +k];
        BUF_i[(b*2*n_act+a)*aux_n_ao +k]=BUF2_i[(a*2*n_ao+b)*aux_n_ao +k];
    }
    
    AA_RI_M_r=new double[long(2*n_act*2*n_act*aux_n_ao)];
    AA_RI_M_i=new double[long(2*n_act*2*n_act*aux_n_ao)];
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   2*n_act,aux_n_ao*2*n_act, 2*n_ao,1.0,
                   ACT_r, 2*n_ao,
                   BUF_r,aux_n_ao*2*n_act,0.0,
                   AA_RI_M_r,aux_n_ao*2*n_act);
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   2*n_act,aux_n_ao*2*n_act, 2*n_ao,1.0,///complex conjugate
                   ACT_i, 2*n_ao,
                   BUF_i,aux_n_ao*2*n_act,1.0,
                   AA_RI_M_r,aux_n_ao*2*n_act);
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   2*n_act,aux_n_ao*2*n_act, 2*n_ao,1.0,
                   ACT_r, 2*n_ao,
                   BUF_i,aux_n_ao*2*n_act,0.0,
                   AA_RI_M_i,aux_n_ao*2*n_act);
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   2*n_act,aux_n_ao*2*n_act, 2*n_ao,-1.0,///complex conjugate
                   ACT_i, 2*n_ao,
                   BUF_r,aux_n_ao*2*n_act,1.0,
                   AA_RI_M_i,aux_n_ao*2*n_act);
    
    CA_RI_M_r=new double[long(2*n_act*2*n_cor*aux_n_ao)];
    CA_RI_M_i=new double[long(2*n_act*2*n_cor*aux_n_ao)];
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   2*n_cor,aux_n_ao*2*n_act, 2*n_ao,1.0,
                   COR_r, 2*n_ao,
                   BUF_r,aux_n_ao*2*n_act,0.0,
                   CA_RI_M_r,aux_n_ao*2*n_act);
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   2*n_cor,aux_n_ao*2*n_act, 2*n_ao,1.0,///complex conjugate
                   COR_i, 2*n_ao,
                   BUF_i,aux_n_ao*2*n_act,1.0,
                   CA_RI_M_r,aux_n_ao*2*n_act);
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   2*n_cor,aux_n_ao*2*n_act, 2*n_ao,1.0,
                   COR_r, 2*n_ao,
                   BUF_i,aux_n_ao*2*n_act,0.0,
                   CA_RI_M_i,aux_n_ao*2*n_act);
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   2*n_cor,aux_n_ao*2*n_act, 2*n_ao,-1.0,///complex conjugate
                   COR_i, 2*n_ao,
                   BUF_r,aux_n_ao*2*n_act,1.0,
                   CA_RI_M_i,aux_n_ao*2*n_act);
    
    VA_RI_M_r=new double[long(2*n_act*2*n_vac*aux_n_ao)];
    VA_RI_M_i=new double[long(2*n_act*2*n_vac*aux_n_ao)];
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   2*n_vac,aux_n_ao*2*n_act, 2*n_ao,1.0,
                   VAC_r, 2*n_ao,
                   BUF_r,aux_n_ao*2*n_act,0.0,
                   VA_RI_M_r,aux_n_ao*2*n_act);
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   2*n_vac,aux_n_ao*2*n_act, 2*n_ao,1.0,///complex conjugate
                   VAC_i, 2*n_ao,
                   BUF_i,aux_n_ao*2*n_act,1.0,
                   VA_RI_M_r,aux_n_ao*2*n_act);
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   2*n_vac,aux_n_ao*2*n_act, 2*n_ao,1.0,
                   VAC_r, 2*n_ao,
                   BUF_i,aux_n_ao*2*n_act,0.0,
                   VA_RI_M_i,aux_n_ao*2*n_act);
    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   2*n_vac,aux_n_ao*2*n_act, 2*n_ao,-1.0,///complex conjugate
                   VAC_i, 2*n_ao,
                   BUF_r,aux_n_ao*2*n_act,1.0,
                   VA_RI_M_i,aux_n_ao*2*n_act);
    
    delete[] BUF_r ;
    delete[] BUF_i ;
    delete[] BUF2_r;
    delete[] BUF2_i;
    
    return 0;
}

RI_complex_data::~RI_complex_data(){
    
    if(VC_RI_M_r != nullptr) delete[] VC_RI_M_r;
    if(CA_RI_M_r != nullptr) delete[] CA_RI_M_r;
    if(CC_RI_M_r != nullptr) delete[] CC_RI_M_r;
    if(VA_RI_M_r != nullptr) delete[] VA_RI_M_r;
    if(AA_RI_M_r != nullptr) delete[] AA_RI_M_r;
    
    if(VC_RI_M_i != nullptr) delete[] VC_RI_M_i;
    if(CA_RI_M_i != nullptr) delete[] CA_RI_M_i;
    if(CC_RI_M_i != nullptr) delete[] CC_RI_M_i;
    if(VA_RI_M_i != nullptr) delete[] VA_RI_M_i;
    if(AA_RI_M_i != nullptr) delete[] AA_RI_M_i;
    
}
