# include <stdio.h>
# include "blas_link.h"
# include <omp.h>


# include "CI.h"
# include "matr.h"
# include "timer.h"

# include "from_hash.h"
# include "aldet.h"
# include "RI.h"
# include "UPT_tensors.h"


# define max(a,b)  (((a)<(b))?(b):(a))
# define min(a,b)  (((a)>(b))?(b):(a))

extern int num_threads;

int M_average(double * A, double * M, int n, int m, std::vector<double> w);


inline double ED_with_shift(double E, double edshift){
    
    return E/(E*E+edshift);
}


UPT_tensors::UPT_tensors(){
            
//     RF_PV    = nullptr;
    RF_PH_A    = nullptr;
    RF_PH_B    = nullptr;
    RF_PV_JK_A = nullptr;
    RF_PV_JK_B = nullptr;
    RF_PV_AB   = nullptr;
    RF_PV_BA   = nullptr;
    RF_P3_JK_A = nullptr;
    RF_P3_JK_B = nullptr;
    RF_P3_AB   = nullptr;
    RF_P3_BA   = nullptr;
    
    e_c_A      = nullptr;
    e_a_A      = nullptr;
    e_IP_A     = nullptr;
    e_EA_A     = nullptr;
    e_v_A      = nullptr;
               
    e_c_B      = nullptr;
    e_a_B      = nullptr;
    e_IP_B     = nullptr;
    e_EA_B     = nullptr;
    e_v_B      = nullptr;
               
    IP_U_A     = nullptr;
    EA_U_A     = nullptr;
    IP_U_B     = nullptr;
    EA_U_B     = nullptr;
    IP_A_Um    = nullptr;
    EA_A_Um    = nullptr;
    IP_A_Up    = nullptr;
    EA_A_Up    = nullptr;
               
    e_IP_AA    = nullptr;
    e_EA_AA    = nullptr;
    IP_AA_U    = nullptr;
    EA_AA_U    = nullptr;
    IP_AA_Um   = nullptr;
    EA_AA_Um   = nullptr;
    IP_AA_Up   = nullptr;
    EA_AA_Up   = nullptr;
               
    e_IP_BB    = nullptr;
    e_EA_BB    = nullptr;
    IP_BB_U    = nullptr;
    EA_BB_U    = nullptr;
    IP_BB_Um   = nullptr;
    EA_BB_Um   = nullptr;
    IP_BB_Up   = nullptr;
    EA_BB_Up   = nullptr;
               
    e_IP_AB    = nullptr;
    e_EA_AB    = nullptr;
    IP_AB_U    = nullptr;
    EA_AB_U    = nullptr;
    IP_AB_Um   = nullptr;
    EA_AB_Um   = nullptr;
    IP_AB_Up   = nullptr;
    EA_AB_Up   = nullptr;
    
}

UPT_tensors::~UPT_tensors(){
    
    if(RF_PH_A    != nullptr) delete[] RF_PH_A    ;
    if(RF_PH_B    != nullptr) delete[] RF_PH_B    ;
    if(RF_PV_JK_A != nullptr) delete[] RF_PV_JK_A ;
    if(RF_PV_JK_B != nullptr) delete[] RF_PV_JK_B ;
    if(RF_PV_AB   != nullptr) delete[] RF_PV_AB   ;
    if(RF_PV_BA   != nullptr) delete[] RF_PV_BA   ;
    if(RF_P3_JK_A != nullptr) delete[] RF_P3_JK_A ;
    if(RF_P3_JK_B != nullptr) delete[] RF_P3_JK_B ;
    if(RF_P3_AB   != nullptr) delete[] RF_P3_AB   ;
    if(RF_P3_BA   != nullptr) delete[] RF_P3_BA   ;
    
    if(e_IP_A     != nullptr) delete[] e_IP_A     ;
    if(e_EA_A     != nullptr) delete[] e_EA_A     ;
    
    if(e_IP_B     != nullptr) delete[] e_IP_B     ;
    if(e_EA_B     != nullptr) delete[] e_EA_B     ;
    
    if(IP_U_A     != nullptr) delete[] IP_U_A     ;
    if(EA_U_A     != nullptr) delete[] EA_U_A     ;
    if(IP_U_B     != nullptr) delete[] IP_U_B     ;
    if(EA_U_B     != nullptr) delete[] EA_U_B     ;
    if(IP_A_Um    != nullptr) delete[] IP_A_Um    ;
    if(EA_A_Um    != nullptr) delete[] EA_A_Um    ;
    if(IP_A_Up    != nullptr) delete[] IP_A_Up    ;
    if(EA_A_Up    != nullptr) delete[] EA_A_Up    ;
    
    if(e_IP_AA    != nullptr) delete[] e_IP_AA    ;
    if(e_EA_AA    != nullptr) delete[] e_EA_AA    ;
    if(IP_AA_U    != nullptr) delete[] IP_AA_U    ;
    if(EA_AA_U    != nullptr) delete[] EA_AA_U    ;
    if(IP_AA_Um   != nullptr) delete[] IP_AA_Um   ;
    if(EA_AA_Um   != nullptr) delete[] EA_AA_Um   ;
    if(IP_AA_Up   != nullptr) delete[] IP_AA_Up   ;
    if(EA_AA_Up   != nullptr) delete[] EA_AA_Up   ;
    
    if(e_IP_BB    != nullptr) delete[] e_IP_BB    ;
    if(e_EA_BB    != nullptr) delete[] e_EA_BB    ;
    if(IP_BB_U    != nullptr) delete[] IP_BB_U    ;
    if(EA_BB_U    != nullptr) delete[] EA_BB_U    ;
    if(IP_BB_Um   != nullptr) delete[] IP_BB_Um   ;
    if(EA_BB_Um   != nullptr) delete[] EA_BB_Um   ;
    if(IP_BB_Up   != nullptr) delete[] IP_BB_Up   ;
    if(EA_BB_Up   != nullptr) delete[] EA_BB_Up   ;
    
    if(e_IP_AB    != nullptr) delete[] e_IP_AB    ;
    if(e_EA_AB    != nullptr) delete[] e_EA_AB    ;
    if(IP_AB_U    != nullptr) delete[] IP_AB_U    ;
    if(EA_AB_U    != nullptr) delete[] EA_AB_U    ;
    if(IP_AB_Um   != nullptr) delete[] IP_AB_Um   ;
    if(EA_AB_Um   != nullptr) delete[] EA_AB_Um   ;
    if(IP_AB_Up   != nullptr) delete[] IP_AB_Up   ;
    if(EA_AB_Up   != nullptr) delete[] EA_AB_Up   ;

}


int UPT_tensors::set_par(RI_data * ext_RI_A, RI_data * ext_RI_B,
                         double * ext_eps, double * ext_eps_B,
                         int ext_n_c, int ext_n_a, int ext_n_v, 
                         double * ext_H_AV_A, double * ext_H_CA_A, double * ext_H_CV_A,
                         double * ext_H_AV_B, double * ext_H_CA_B, double * ext_H_CV_B,
                         double ext_edshift){
//            
    RI_A     = ext_RI_A      ;
    RI_B     = ext_RI_B      ;
    aux_n_ao = RI_A->aux_n_ao;
    n_c      = ext_n_c       ; 
    n_a      = ext_n_a       ;
    n_v      = ext_n_v       ;
    H_AV_A   = ext_H_AV_A    ;
    H_CA_A   = ext_H_CA_A    ;
    H_CV_A   = ext_H_CV_A    ;
    H_AV_B   = ext_H_AV_B    ;
    H_CA_B   = ext_H_CA_B    ;
    H_CV_B   = ext_H_CV_B    ;
    edshift  = ext_edshift   ;
        
    e_c_A    = ext_eps        ; 
    e_a_A    = ext_eps+n_c    ; 
    e_v_A    = ext_eps+n_c+n_a;
    
    e_c_B   = ext_eps_B        ; 
    e_a_B   = ext_eps_B+n_c    ; 
    e_v_B   = ext_eps_B+n_c+n_a;
    
    e_IP_A= new double[n_a];for(int i=0; i<n_a; i++) e_IP_A[i]=e_a_A[i];
    e_EA_A= new double[n_a];for(int i=0; i<n_a; i++) e_EA_A[i]=e_a_A[i];
    
    e_IP_B= new double[n_a];for(int i=0; i<n_a; i++) e_IP_B[i]=e_a_B[i];
    e_EA_B= new double[n_a];for(int i=0; i<n_a; i++) e_EA_B[i]=e_a_B[i];
    
    IP_U_A= new double[n_a*n_a];
    EA_U_A= new double[n_a*n_a];
    IP_U_B= new double[n_a*n_a];
    EA_U_B= new double[n_a*n_a];
    IP_A_Um= new double[n_a*n_a];
    EA_A_Um= new double[n_a*n_a];
    IP_A_Up= new double[n_a*n_a];
    EA_A_Up= new double[n_a*n_a];

    e_IP_AA= new double[n_a*n_a];
    e_EA_AA= new double[n_a*n_a];
    IP_AA_U= new double[n_a*n_a*n_a*n_a];
    EA_AA_U= new double[n_a*n_a*n_a*n_a];
    IP_AA_Um= new double[n_a*n_a*n_a*n_a];
    EA_AA_Um= new double[n_a*n_a*n_a*n_a];
    IP_AA_Up= new double[n_a*n_a*n_a*n_a];
    EA_AA_Up= new double[n_a*n_a*n_a*n_a];
    
    e_IP_BB= new double[n_a*n_a];
    e_EA_BB= new double[n_a*n_a];
    IP_BB_U= new double[n_a*n_a*n_a*n_a];
    EA_BB_U= new double[n_a*n_a*n_a*n_a];
    IP_BB_Um= new double[n_a*n_a*n_a*n_a];
    EA_BB_Um= new double[n_a*n_a*n_a*n_a];
    IP_BB_Up= new double[n_a*n_a*n_a*n_a];
    EA_BB_Up= new double[n_a*n_a*n_a*n_a];
    
    e_IP_AB= new double[n_a*n_a];
    e_EA_AB= new double[n_a*n_a];
    IP_AB_U= new double[n_a*n_a*n_a*n_a];
    EA_AB_U= new double[n_a*n_a*n_a*n_a];
    IP_AB_Um= new double[n_a*n_a*n_a*n_a];
    EA_AB_Um= new double[n_a*n_a*n_a*n_a];
    IP_AB_Up= new double[n_a*n_a*n_a*n_a];
    EA_AB_Up= new double[n_a*n_a*n_a*n_a];
        

    
//     RF_PV    = new double  [n_a*n_a*n_a*n_a        ];
    RF_PH_A    = new double  [n_a*n_a                ];
    RF_PH_B    = new double  [n_a*n_a                ];
    RF_PV_JK_A = new double  [n_a*n_a*n_a*n_a        ];
    RF_PV_JK_B = new double  [n_a*n_a*n_a*n_a        ];
    RF_PV_AB   = new double  [n_a*n_a*n_a*n_a        ];
    RF_PV_BA   = new double  [n_a*n_a*n_a*n_a        ];
    RF_P3_JK_A = new double  [n_a*n_a*n_a*n_a*n_a*n_a];
    RF_P3_JK_B = new double  [n_a*n_a*n_a*n_a*n_a*n_a];
    RF_P3_AB = new double  [n_a*n_a*n_a*n_a*n_a*n_a];
    RF_P3_BA = new double  [n_a*n_a*n_a*n_a*n_a*n_a];
    
    return 0;
}


int UPT_tensors::MPPT(aldet_data * I, int i_set, std::vector<double> avecoe){
        
    if(n_a!= I->n_act){
        fprintf(out_stream,"ERROR in MPPT\n");
        exit(0);
    }
    
    int n_s =I->n_states[i_set];
    
    double * IP_H_A   = new double[n_a*n_a];
    double * IP_H_B   = new double[n_a*n_a];
    double * EA_H_A   = new double[n_a*n_a];
    double * EA_H_B   = new double[n_a*n_a];
    
    
    I->calc_IPEA_single_AB(IP_U_A, IP_H_A, EA_U_A, EA_H_A, 
                           IP_U_B, IP_H_B, EA_U_B, EA_H_B, 
                           0 ,avecoe);
    
    
    
    
    for(int i=0;i<n_a;i++)e_IP_A[i]=IP_H_A[i*n_a+i]/IP_U_A[i*n_a+i];
    for(int i=0;i<n_a;i++)e_EA_A[i]=EA_H_A[i*n_a+i]/EA_U_A[i*n_a+i];
    for(int i=0;i<n_a;i++)e_IP_B[i]=IP_H_B[i*n_a+i]/IP_U_B[i*n_a+i];
    for(int i=0;i<n_a;i++)e_EA_B[i]=EA_H_B[i*n_a+i]/EA_U_B[i*n_a+i];
    
    fprintf(out_stream,"\n\nMPPT orbital enegries:\n");
    fprintf(out_stream,"IP(A):       ");fPrintMatr(out_stream,e_IP_A,1,n_a,0);
    fprintf(out_stream,"EA(A):       ");fPrintMatr(out_stream,e_EA_A,1,n_a,0);
    fprintf(out_stream,"IP(B):       ");fPrintMatr(out_stream,e_IP_B,1,n_a,0);
    fprintf(out_stream,"EA(B):       ");fPrintMatr(out_stream,e_EA_B,1,n_a,0);fprintf(out_stream,"_______________________________________________________________________\n\n\n");
    
    delete[] IP_H_A   ;
    delete[] EA_H_A   ;
    delete[] IP_H_B   ;
    delete[] EA_H_B   ;
    
    return 0;
    
}


int UPT_tensors::E2_calc_EE(){
    
    set_zero_matr(RF_P3_JK_A, n_a*n_a*n_a*n_a*n_a*n_a);
    set_zero_matr(RF_P3_JK_B, n_a*n_a*n_a*n_a*n_a*n_a);
    set_zero_matr(RF_P3_AB  , n_a*n_a*n_a*n_a*n_a*n_a);
    set_zero_matr(RF_P3_BA  , n_a*n_a*n_a*n_a*n_a*n_a);
    set_zero_matr(RF_PV_JK_A, n_a*n_a*n_a*n_a);
    set_zero_matr(RF_PV_JK_B, n_a*n_a*n_a*n_a);
    set_zero_matr(RF_PV_AB  , n_a*n_a*n_a*n_a);
    set_zero_matr(RF_PV_BA  , n_a*n_a*n_a*n_a);
    
    set_zero_matr(RF_PH_A, n_a*n_a);
    set_zero_matr(RF_PH_B, n_a*n_a);
    
//     set_zero_matr(RF_PS, N_fit);
    RF_PS=0;
    
    calc_EE_2_CCVV();printf_timer("CCVV res table");fprintf(out_stream,"\n");fflush(out_stream);
    calc_EE_2_CAVV();printf_timer("CAVV res table");fprintf(out_stream,"\n");fflush(out_stream);
    calc_EE_2_AAVV();printf_timer("AAVV res table");fprintf(out_stream,"\n");fflush(out_stream);
    calc_EE_2_CCAV();printf_timer("CCAV res table");fprintf(out_stream,"\n");fflush(out_stream);
    calc_EE_2_CCAA();printf_timer("CCAA res table");fprintf(out_stream,"\n");fflush(out_stream);
    calc_EE_2_CV  ();printf_timer("CV   res table");fprintf(out_stream,"\n");fflush(out_stream);
    calc_EE_2_AV  ();printf_timer("AV   res table");fprintf(out_stream,"\n");fflush(out_stream);//must be checked
    calc_EE_2_CA  ();printf_timer("CA   res table");fprintf(out_stream,"\n");fflush(out_stream);
    
    symmetrization(RF_PH_A, n_a);
    symmetrization(RF_PH_B, n_a);
    symmetrization(RF_PV_JK_A, n_a*n_a);
    symmetrization(RF_PV_JK_B, n_a*n_a);
    symmetrization(RF_PV_AB  , n_a*n_a);
    symmetrization(RF_PV_BA  , n_a*n_a);
    symmetrization(RF_P3_JK_A, n_a*n_a*n_a);
    symmetrization(RF_P3_JK_B, n_a*n_a*n_a);
    symmetrization(RF_P3_AB  , n_a*n_a*n_a);
    symmetrization(RF_P3_BA  , n_a*n_a*n_a);
    
    double * RF_P3_tmp = new double  [n_a*n_a*n_a*n_a*n_a*n_a];
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++)
    for(int v=0; v<n_a; v++)
    for(int w=0; w<n_a; w++)
    for(int x=0; x<n_a; x++)
    for(int y=0; y<n_a; y++){
        RF_P3_tmp[((((t*n_a+u)*n_a+v)*n_a+w)*n_a+x)*n_a+y]=    //AAAABB
        RF_P3_AB    [((((t*n_a+u)*n_a+x)*n_a+v)*n_a+w)*n_a+y]; //AABAAB
    }
    memcpy(RF_P3_AB,RF_P3_tmp,n_a*n_a*n_a*n_a*n_a*n_a*sizeof(double));
    
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++)
    for(int v=0; v<n_a; v++)
    for(int w=0; w<n_a; w++)
    for(int x=0; x<n_a; x++)
    for(int y=0; y<n_a; y++){
        RF_P3_tmp[((((t*n_a+u)*n_a+v)*n_a+w)*n_a+x)*n_a+y]=    //BBBBAA
        RF_P3_BA    [((((t*n_a+u)*n_a+x)*n_a+v)*n_a+w)*n_a+y]; //BBABBA
    }
    memcpy(RF_P3_BA,RF_P3_tmp,n_a*n_a*n_a*n_a*n_a*n_a*sizeof(double));
    
    delete[]RF_P3_tmp;
    
    return 0;
    
}

int UPT_tensors::calc_EE_2_CV(){
    
    double dE;
    double V;
    
    
    double * RF_MM_A ;    // (J-K) * (J-K)
    double * RF_MM_B ;    // (J-K) * (J-K)
    double * RF_JJ_A ;    //  J(aa)* J(aa)
    double * RF_JJ_B ;    //  J(aa)* J(aa)
    double * RF_JM_AB;    //  J(aa)*(J-K)
    double * RF_JM_BA;    //  J(bb)*(J-K)
    double * RF_MJ_AB;    // (J-K) * J(bb)
    double * RF_MJ_BA;    // (J-K) * J(aa)
    double * RF_ABBA ;    //  J(ab)* J(ba)
    double * RF_BAAB ;    //  J(ba)* J(ab)
    
    RF_MM_A  = new double  [n_a*n_a*n_a*n_a];
    RF_MM_B  = new double  [n_a*n_a*n_a*n_a];
    RF_JJ_A  = new double  [n_a*n_a*n_a*n_a];
    RF_JJ_B  = new double  [n_a*n_a*n_a*n_a];
    RF_JM_AB = new double  [n_a*n_a*n_a*n_a];
    RF_JM_BA = new double  [n_a*n_a*n_a*n_a];
    RF_MJ_AB = new double  [n_a*n_a*n_a*n_a];
    RF_MJ_BA = new double  [n_a*n_a*n_a*n_a];
    RF_ABBA  = new double  [n_a*n_a*n_a*n_a];
    RF_BAAB  = new double  [n_a*n_a*n_a*n_a];
    set_zero_matr(RF_MM_A , n_a*n_a*n_a*n_a);
    set_zero_matr(RF_MM_B , n_a*n_a*n_a*n_a);
    set_zero_matr(RF_JJ_A , n_a*n_a*n_a*n_a);
    set_zero_matr(RF_JJ_B , n_a*n_a*n_a*n_a);
    set_zero_matr(RF_JM_AB, n_a*n_a*n_a*n_a);
    set_zero_matr(RF_JM_BA, n_a*n_a*n_a*n_a);
    set_zero_matr(RF_MJ_AB, n_a*n_a*n_a*n_a);
    set_zero_matr(RF_MJ_BA, n_a*n_a*n_a*n_a);
    set_zero_matr(RF_ABBA , n_a*n_a*n_a*n_a);
    set_zero_matr(RF_BAAB , n_a*n_a*n_a*n_a);
    
    double * RF_HM_A;
    double * RF_MH_A;
    double * RF_HJ_A;
    double * RF_JH_A;
    RF_HM_A = new double  [n_a*n_a];
    RF_MH_A = new double  [n_a*n_a];
    RF_HJ_A = new double  [n_a*n_a];
    RF_JH_A = new double  [n_a*n_a];
    set_zero_matr(RF_HM_A, n_a*n_a);
    set_zero_matr(RF_MH_A, n_a*n_a);
    set_zero_matr(RF_HJ_A, n_a*n_a);
    set_zero_matr(RF_JH_A, n_a*n_a);
    
    double * RF_HM_B;
    double * RF_MH_B;
    double * RF_HJ_B;
    double * RF_JH_B;
    RF_HM_B = new double  [n_a*n_a];
    RF_MH_B = new double  [n_a*n_a];
    RF_HJ_B = new double  [n_a*n_a];
    RF_JH_B = new double  [n_a*n_a];
    set_zero_matr(RF_HM_B, n_a*n_a);
    set_zero_matr(RF_MH_B, n_a*n_a);
    set_zero_matr(RF_HJ_B, n_a*n_a);
    set_zero_matr(RF_JH_B, n_a*n_a);
    
    double RF_H=0;
//     RF_H = new double  [n_a];
//     set_zero_matr(RF_H, n_a);
    
    double * VCAA_AA;
    VCAA_AA = new double[n_c*n_v*n_a*n_a];
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        n_c*n_v,n_a*n_a,aux_n_ao,1.0,
                        RI_A->VC_RI_M,aux_n_ao,
                        RI_A->AA_RI_M,aux_n_ao,0.0,
                        VCAA_AA,n_a*n_a);
    
    double * VCAA_AB;
    VCAA_AB = new double[n_c*n_v*n_a*n_a];
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        n_c*n_v,n_a*n_a,aux_n_ao,1.0,
                        RI_A->VC_RI_M,aux_n_ao,
                        RI_B->AA_RI_M,aux_n_ao,0.0,
                        VCAA_AB,n_a*n_a);
    
    double * VCAA_BA;
    VCAA_BA = new double[n_c*n_v*n_a*n_a];
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        n_c*n_v,n_a*n_a,aux_n_ao,1.0,
                        RI_B->VC_RI_M,aux_n_ao,
                        RI_A->AA_RI_M,aux_n_ao,0.0,
                        VCAA_BA,n_a*n_a);
    
    double * VCAA_BB;
    VCAA_BB = new double[n_c*n_v*n_a*n_a];
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        n_c*n_v,n_a*n_a,aux_n_ao,1.0,
                        RI_B->VC_RI_M,aux_n_ao,
                        RI_B->AA_RI_M,aux_n_ao,0.0,
                        VCAA_BB,n_a*n_a);
    
    
    double * CAVA_AA;
    CAVA_AA = new double[n_a*n_c*n_a*n_v];
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        n_a*n_c,n_a*n_v,aux_n_ao,1.0,
                        RI_A->CA_RI_M,aux_n_ao,
                        RI_A->VA_RI_M,aux_n_ao,0.0,
                        CAVA_AA,n_a*n_v);
    
    double * CAVA_AB;
    CAVA_AB = new double[n_a*n_c*n_a*n_v];
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        n_a*n_c,n_a*n_v,aux_n_ao,1.0,
                        RI_A->CA_RI_M,aux_n_ao,
                        RI_B->VA_RI_M,aux_n_ao,0.0,
                        CAVA_AB,n_a*n_v);
    
    double * CAVA_BA;
    CAVA_BA = new double[n_a*n_c*n_a*n_v];
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        n_a*n_c,n_a*n_v,aux_n_ao,1.0,
                        RI_B->CA_RI_M,aux_n_ao,
                        RI_A->VA_RI_M,aux_n_ao,0.0,
                        CAVA_BA,n_a*n_v);
    
    double * CAVA_BB;
    CAVA_BB = new double[n_a*n_c*n_a*n_v];
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        n_a*n_c,n_a*n_v,aux_n_ao,1.0,
                        RI_B->CA_RI_M,aux_n_ao,
                        RI_B->VA_RI_M,aux_n_ao,0.0,
                        CAVA_BB,n_a*n_v);
//     double * VCAA_J;
//     VCAA_J = new double[n_c*n_v];
//     double * VCAA_K;
//     VCAA_K = new double[n_c*n_v];

    for(int i=0; i<n_c; i++)
    for(int u=0; u<n_a; u++)
    for(int a=0; a<n_v; a++)
    for(int t=0; t<n_a; t++){
        dE=e_c_A[i]+e_IP_B[t]-e_v_B[a]-e_EA_A[u];//is taken negative
        if(t==u)dE=e_c_A[i]-e_v_B[a];
        dE=dE/(dE*dE+edshift);
        
        for(int v=0; v<n_a; v++)
        for(int w=0; w<n_a; w++){
            //<i(a)t(b)|u(a)a(b)><v(a)a(b)|i(a)w(b)>=(iu|at)(iv|aw)
            V=CAVA_AB[((i*n_a+u)*n_v+a)*n_a+t]*
              CAVA_AB[((i*n_a+v)*n_v+a)*n_a+w];
              
            RF_BAAB[((t*n_a+u)*n_a+v)*n_a+w]+=V*dE;//BAAB
            
        }
    }
    for(int i=0; i<n_c; i++)
    for(int u=0; u<n_a; u++)
    for(int a=0; a<n_v; a++)
    for(int t=0; t<n_a; t++){
        dE=e_c_B[i]+e_IP_A[t]-e_v_A[a]-e_EA_B[u];//is taken negative
        if(t==u)dE=e_c_B[i]-e_v_A[a];
        dE=dE/(dE*dE+edshift);
        
        for(int v=0; v<n_a; v++)
        for(int w=0; w<n_a; w++){
            //<i(b)t(a)|u(b)a(a)><v(b)a(a)|i(b)w(a)>=(iu|at)(iv|aw)
            V=CAVA_BA[((i*n_a+u)*n_v+a)*n_a+t]*
              CAVA_BA[((i*n_a+v)*n_v+a)*n_a+w];
              
            RF_ABBA[((t*n_a+u)*n_a+v)*n_a+w]+=V*dE;//ABBA
            
        }
    }
    
    
    for(int i=0; i<n_c; i++)
    for(int a=0; a<n_v; a++)
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++){
        dE=e_c_A[i]+e_IP_A[t]-e_v_A[a]-e_EA_A[u];//is taken negative
        if(t==u)dE=e_c_A[i]-e_v_A[a];
        dE=dE/(dE*dE+edshift);
        
        for(int v=0; v<n_a; v++)
        for(int w=0; w<n_a; w++){
            //(<it|au>-<iu|at>)(<av|iw>-<av|wi>)
            V=(VCAA_AA[((a*n_c+i)*n_a+u)*n_a+t]-  //<it|au>=(ai|tu)
               CAVA_AA[((i*n_a+u)*n_v+a)*n_a+t])* //<it|ua>=(iu|at)
              (VCAA_AA[((a*n_c+i)*n_a+w)*n_a+v]-  //<av|iw>=(ai|wv)
               CAVA_AA[((i*n_a+v)*n_v+a)*n_a+w]); //<av|wi>=(iv|aw)
            
            RF_MM_A[((t*n_a+u)*n_a+v)*n_a+w]+=V*dE;//AAAA
        }
    }
    for(int i=0; i<n_c; i++)
    for(int a=0; a<n_v; a++)
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++){
        dE=e_c_B[i]+e_IP_B[t]-e_v_B[a]-e_EA_B[u];//is taken negative
        if(t==u)dE=e_c_B[i]-e_v_B[a];
        dE=dE/(dE*dE+edshift);
        
        for(int v=0; v<n_a; v++)
        for(int w=0; w<n_a; w++){
            //(<it|au>-<iu|at>)(<av|iw>-<av|wi>)
            V=(VCAA_BB[((a*n_c+i)*n_a+u)*n_a+t]-  //<it|au>=(ai|tu)
               CAVA_BB[((i*n_a+u)*n_v+a)*n_a+t])* //<it|ua>=(iu|at)
              (VCAA_BB[((a*n_c+i)*n_a+w)*n_a+v]-  //<av|iw>=(ai|wv)
               CAVA_BB[((i*n_a+v)*n_v+a)*n_a+w]); //<av|wi>=(iv|aw)
            
            RF_MM_B[((t*n_a+u)*n_a+v)*n_a+w]+=V*dE;//BBBB
        }
    }
    
    

    for(int i=0; i<n_c; i++)
    for(int a=0; a<n_v; a++)
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++){
        dE=e_c_A[i]+e_IP_A[t]-e_v_A[a]-e_EA_A[u];//is taken negative
        if(t==u)dE=e_c_A[i]-e_v_A[a];
        dE=dE/(dE*dE+edshift);
        
        for(int v=0; v<n_a; v++)
        for(int w=0; w<n_a; w++){
            //(<it|au>-<iu|at>)<av|iw>
            V=(VCAA_AA[((a*n_c+i)*n_a+u)*n_a+t]-  //<it|au>=(ai|tu)
               CAVA_AA[((i*n_a+u)*n_v+a)*n_a+t])* //<it|ua>=(iu|at)
               VCAA_AB[((a*n_c+i)*n_a+w)*n_a+v];  //<av|iw>=(ai|wv)
            
            RF_MJ_AB[((t*n_a+u)*n_a+v)*n_a+w]+=V*dE;//AABB
        }
    }
    for(int i=0; i<n_c; i++)
    for(int a=0; a<n_v; a++)
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++){
        dE=e_c_B[i]+e_IP_B[t]-e_v_B[a]-e_EA_B[u];//is taken negative
        if(t==u)dE=e_c_B[i]-e_v_B[a];
        dE=dE/(dE*dE+edshift);
        
        for(int v=0; v<n_a; v++)
        for(int w=0; w<n_a; w++){
            //(<it|au>-<iu|at>)<av|iw>
            V=(VCAA_BB[((a*n_c+i)*n_a+u)*n_a+t]-  //<it|au>=(ai|tu)
               CAVA_BB[((i*n_a+u)*n_v+a)*n_a+t])* //<it|ua>=(iu|at)
               VCAA_BA[((a*n_c+i)*n_a+w)*n_a+v];  //<av|iw>=(ai|wv)
            
            RF_MJ_BA[((t*n_a+u)*n_a+v)*n_a+w]+=V*dE;//BBAA
        }
    }
    
    
    for(int i=0; i<n_c; i++)
    for(int a=0; a<n_v; a++)
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++){
        dE=e_c_A[i]+e_IP_B[t]-e_v_A[a]-e_EA_B[u];
        if(t==u)dE=e_c_A[i]-e_v_A[a];
        dE=dE/(dE*dE+edshift);
        
        for(int v=0; v<n_a; v++)
        for(int w=0; w<n_a; w++){
            //<it|au>(<av|iw>-<av|wi>)
            V= VCAA_AB[((a*n_c+i)*n_a+u)*n_a+t]*  //<it|au>=(ai|tu)
              (VCAA_AA[((a*n_c+i)*n_a+w)*n_a+v]-  //<av|iw>=(ai|wv)
               CAVA_AA[((i*n_a+v)*n_v+a)*n_a+w]); //<av|wi>=(iv|aw)
            
            RF_JM_BA[((t*n_a+u)*n_a+v)*n_a+w]+=V*dE;//BBAA
        }
    }

    for(int i=0; i<n_c; i++)
    for(int a=0; a<n_v; a++)
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++){
        dE=e_c_B[i]+e_IP_A[t]-e_v_B[a]-e_EA_A[u];
        if(t==u)dE=e_c_B[i]-e_v_B[a];
        dE=dE/(dE*dE+edshift);
        
        for(int v=0; v<n_a; v++)
        for(int w=0; w<n_a; w++){
            //<it|au>(<av|iw>-<av|wi>)
            V= VCAA_BA[((a*n_c+i)*n_a+u)*n_a+t]*  //<it|au>=(ai|tu)
              (VCAA_BB[((a*n_c+i)*n_a+w)*n_a+v]-  //<av|iw>=(ai|wv)
               CAVA_BB[((i*n_a+v)*n_v+a)*n_a+w]); //<av|wi>=(iv|aw)
            
            RF_JM_AB[((t*n_a+u)*n_a+v)*n_a+w]+=V*dE;//AABB
        }
    }

    
    for(int i=0; i<n_c; i++)
    for(int a=0; a<n_v; a++)
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++){
        dE=e_c_B[i]+e_IP_A[t]-e_v_B[a]-e_EA_A[u];
        if(t==u)dE=e_c_B[i]-e_v_B[a];
        dE=dE/(dE*dE+edshift);
        
        for(int v=0; v<n_a; v++)
        for(int w=0; w<n_a; w++){
            //<it|au><av|iw>
            V=VCAA_BA[((a*n_c+i)*n_a+u)*n_a+t]*  //<it|au>=(ai|tu)
              VCAA_BA[((a*n_c+i)*n_a+w)*n_a+v];  //<av|iw>=(ai|wv)
            
            RF_JJ_A[((t*n_a+u)*n_a+v)*n_a+w]+=V*dE;//AAAA
        }
    }
    
    for(int i=0; i<n_c; i++)
    for(int a=0; a<n_v; a++)
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++){
        dE=e_c_A[i]+e_IP_B[t]-e_v_A[a]-e_EA_B[u];
        if(t==u)dE=e_c_A[i]-e_v_A[a];
        dE=dE/(dE*dE+edshift);
        
        for(int v=0; v<n_a; v++)
        for(int w=0; w<n_a; w++){
            //<it|au><av|iw>
            V=VCAA_AB[((a*n_c+i)*n_a+u)*n_a+t]*  //<it|au>=(ai|tu)
              VCAA_AB[((a*n_c+i)*n_a+w)*n_a+v];  //<av|iw>=(ai|wv)
            
            RF_JJ_B[((t*n_a+u)*n_a+v)*n_a+w]+=V*dE;//BBBB
        }
    }
//     printf_timer("opt table");
    
    
    for(int i=0; i<n_c; i++)
    for(int a=0; a<n_v; a++)
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++){
        //(<it|au>-<iu|at>)*H[i,a]
        V=(VCAA_AA[((a*n_c+i)*n_a+u)*n_a+t]-  //<it|au>=(ai|tu)
           CAVA_AA[((i*n_a+u)*n_v+a)*n_a+t])* //<it|ua>=(iu|at)
           H_CV_A[i*n_v+a]; 
        dE=e_c_A[i]+e_IP_A[t]-e_v_A[a]-e_EA_A[u];
        if(t==u)dE=e_c_A[i]-e_v_A[a];
        RF_MH_A[t*n_a+u]+=V*dE/(dE*dE+edshift);//AA
    }
    for(int i=0; i<n_c; i++)
    for(int a=0; a<n_v; a++)
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++){
        //(<it|au>-<iu|at>)*H[i,a]
        V=(VCAA_BB[((a*n_c+i)*n_a+u)*n_a+t]-  //<it|au>=(ai|tu)
           CAVA_BB[((i*n_a+u)*n_v+a)*n_a+t])* //<it|ua>=(iu|at)
           H_CV_B[i*n_v+a]; 
        dE=e_c_B[i]+e_IP_B[t]-e_v_B[a]-e_EA_B[u];
        if(t==u)dE=e_c_B[i]-e_v_B[a];
        RF_MH_B[t*n_a+u]+=V*dE/(dE*dE+edshift);//BB
    }
    
    
    for(int i=0; i<n_c; i++)
    for(int a=0; a<n_v; a++)
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++){
        //(<it|au>-<iu|at>)*H[i,a]
        V=(VCAA_AA[((a*n_c+i)*n_a+u)*n_a+t]-  //<it|au>=(ai|tu)
           CAVA_AA[((i*n_a+u)*n_v+a)*n_a+t])* //<it|ua>=(iu|at)
           H_CV_A [i*n_v+a]; 
        
        dE=e_c_A[i]-e_v_A[a];
        
        RF_HM_A[t*n_a+u]+=V*dE/(dE*dE+edshift);//AA
    }

    for(int i=0; i<n_c; i++)
    for(int a=0; a<n_v; a++)
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++){
        //(<it|au>-<iu|at>)*H[i,a]
        V=(VCAA_BB[((a*n_c+i)*n_a+u)*n_a+t]-  //<it|au>=(ai|tu)
           CAVA_BB[((i*n_a+u)*n_v+a)*n_a+t])* //<it|ua>=(iu|at)
           H_CV_B [i*n_v+a]; 
        
        dE=e_c_B[i]-e_v_B[a];
        
        RF_HM_B[t*n_a+u]+=V*dE/(dE*dE+edshift);//BB
    }

    for(int i=0; i<n_c; i++)
    for(int a=0; a<n_v; a++)
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++){
        //(<it|au>-<iu|at>)*H[i,a]
        V= VCAA_AB[((a*n_c+i)*n_a+u)*n_a+t]*  //<it|au>=(ai|tu)
           H_CV_A [i*n_v+a]; 
        dE=e_c_A[i]+e_IP_B[t]-e_v_A[a]-e_EA_B[u];
        if(t==u)dE=e_c_A[i]-e_v_A[a];
        RF_JH_B[t*n_a+u]+=V*dE/(dE*dE+edshift);//BB
        
    }
    
    for(int i=0; i<n_c; i++)
    for(int a=0; a<n_v; a++)
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++){
        //(<it|au>-<iu|at>)*H[i,a]
        V= VCAA_BA[((a*n_c+i)*n_a+u)*n_a+t]*  //<it|au>=(ai|tu)
           H_CV_B [i*n_v+a]; 
        dE=e_c_B[i]+e_IP_A[t]-e_v_B[a]-e_EA_A[u];
        if(t==u)dE=e_c_B[i]-e_v_B[a];
        RF_JH_A[t*n_a+u]+=V*dE/(dE*dE+edshift);//AA
        
    }
    
    for(int i=0; i<n_c; i++)
    for(int a=0; a<n_v; a++)
    for(int t=0  ;t<n_a;t++)
    for(int u=0  ;u<n_a;u++){
        //(<it|au>-<iu|at>)*H[i,a]
        V= VCAA_AB[((a*n_c+i)*n_a+u)*n_a+t]*  //<it|au>=(ai|tu)
           H_CV_A[i*n_v+a]; 
        dE=e_c_A[i]-e_v_A[a];
        RF_HJ_B[t*n_a+u]+=V*dE/(dE*dE+edshift);//BB
        
    }
    
    for(int i=0; i<n_c; i++)
    for(int a=0; a<n_v; a++)
    for(int t=0  ;t<n_a;t++)
    for(int u=0  ;u<n_a;u++){
        //(<it|au>-<iu|at>)*H[i,a]
        V= VCAA_BA[((a*n_c+i)*n_a+u)*n_a+t]*  //<it|au>=(ai|tu)
           H_CV_B [i*n_v+a]; 
        dE=e_c_B[i]-e_v_B[a];
        RF_HJ_A[t*n_a+u]+=V*dE/(dE*dE+edshift);//AA
        
    }
    
    for(int i=0; i<n_c; i++)
    for(int a=0; a<n_v; a++){
        //H[i,a]*H[i,a]
        V=H_CV_A[i*n_v+a]*H_CV_A[i*n_v+a]; 
        dE=e_c_A[i]-e_v_A[a];
        RF_H+=V*dE/(dE*dE+edshift);//A part
        
    }
    
    for(int i=0; i<n_c; i++)
    for(int a=0; a<n_v; a++){
        //H[i,a]*H[i,a]
        V=H_CV_B[i*n_v+a]*H_CV_B[i*n_v+a]; 
        dE=e_c_B[i]-e_v_B[a];
        RF_H+=V*dE/(dE*dE+edshift);//B part
        
    }
    
    
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++)
    for(int v=0; v<n_a; v++)
    for(int w=0; w<n_a; w++)
        RF_PV_JK_A[((u*n_a+t)*n_a+v)*n_a+w]+=+RF_MM_A[((t*n_a+v)*n_a+u)*n_a+w]
                                             +RF_JJ_A[((t*n_a+v)*n_a+u)*n_a+w]
                                             -RF_MM_A[((u*n_a+v)*n_a+t)*n_a+w]
                                             -RF_JJ_A[((u*n_a+v)*n_a+t)*n_a+w]
                                             -RF_MM_A[((t*n_a+w)*n_a+u)*n_a+v]
                                             -RF_JJ_A[((t*n_a+w)*n_a+u)*n_a+v]
                                             +RF_MM_A[((u*n_a+w)*n_a+t)*n_a+v]
                                             +RF_JJ_A[((u*n_a+w)*n_a+t)*n_a+v]
                                             ;
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++)
    for(int v=0; v<n_a; v++)
    for(int w=0; w<n_a; w++)
        RF_PV_JK_B[((u*n_a+t)*n_a+v)*n_a+w]+=+RF_MM_B[((t*n_a+v)*n_a+u)*n_a+w]
                                             +RF_JJ_B[((t*n_a+v)*n_a+u)*n_a+w]
                                             -RF_MM_B[((u*n_a+v)*n_a+t)*n_a+w]
                                             -RF_JJ_B[((u*n_a+v)*n_a+t)*n_a+w]
                                             -RF_MM_B[((t*n_a+w)*n_a+u)*n_a+v]
                                             -RF_JJ_B[((t*n_a+w)*n_a+u)*n_a+v]
                                             +RF_MM_B[((u*n_a+w)*n_a+t)*n_a+v]
                                             +RF_JJ_B[((u*n_a+w)*n_a+t)*n_a+v]
                                             ;
                                                       
    
    
    
//     set_zero_matr(RF_PV_AB, n_a*n_a*n_a*n_a*N_fit);
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++)
    for(int v=0; v<n_a; v++)
    for(int w=0; w<n_a; w++)
        //         A     B       A     B
        RF_PV_AB[((t*n_a+u)*n_a+v)*n_a+w]+=+RF_JM_AB[((t*n_a+v)*n_a+u)*n_a+w]//AABB
                                           +RF_MJ_AB[((t*n_a+v)*n_a+u)*n_a+w]//AABB
                                           +RF_JM_BA[((u*n_a+w)*n_a+t)*n_a+v]//BBAA
                                           +RF_MJ_BA[((u*n_a+w)*n_a+t)*n_a+v]//BBAA
                                           -RF_ABBA [((t*n_a+w)*n_a+u)*n_a+v]//ABBA
                                           -RF_BAAB [((u*n_a+v)*n_a+t)*n_a+w]//BAAB
                                           ;
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++)
    for(int v=0; v<n_a; v++)
    for(int w=0; w<n_a; w++)
        //         B     A       B     A
        RF_PV_BA[((t*n_a+u)*n_a+v)*n_a+w]+=+RF_JM_BA[((t*n_a+v)*n_a+u)*n_a+w]//BBAA
                                           +RF_MJ_BA[((t*n_a+v)*n_a+u)*n_a+w]//BBAA
                                           +RF_JM_AB[((u*n_a+w)*n_a+t)*n_a+v]//AABB
                                           +RF_MJ_AB[((u*n_a+w)*n_a+t)*n_a+v]//AABB
                                           -RF_BAAB [((t*n_a+w)*n_a+u)*n_a+v]//BAAB
                                           -RF_ABBA [((u*n_a+v)*n_a+t)*n_a+w]//ABBA
                                           ;
    
//     set_zero_matr(RF_PH,N_fit*n_a*n_a);
    for(int t=0; t<n_a; t++)
    for(int w=0; w<n_a; w++){
        RF_PH_A[t*n_a+w]+= RF_HM_A[w*n_a+t]
                          +RF_MH_A[t*n_a+w]
                          +RF_HJ_A[w*n_a+t]
                          +RF_JH_A[t*n_a+w];
        for(int u=0; u<n_a; u++)
             RF_PH_A[t*n_a+w]+= RF_MM_A[((t*n_a+u)*n_a+u)*n_a+w]
                               +RF_JJ_A[((t*n_a+u)*n_a+u)*n_a+w]
                               +RF_ABBA[((t*n_a+u)*n_a+u)*n_a+w];
    }
   
    for(int t=0; t<n_a; t++)
    for(int w=0; w<n_a; w++){
        RF_PH_B[t*n_a+w]+= RF_HM_B[w*n_a+t]
                          +RF_MH_B[t*n_a+w]
                          +RF_HJ_B[w*n_a+t]
                          +RF_JH_B[t*n_a+w];
        for(int u=0; u<n_a; u++)
             RF_PH_B[t*n_a+w]+= RF_MM_B[((t*n_a+u)*n_a+u)*n_a+w]
                               +RF_JJ_B[((t*n_a+u)*n_a+u)*n_a+w]
                               +RF_BAAB[((t*n_a+u)*n_a+u)*n_a+w];
    }
    RF_PS+=RF_H;//summ A+B
    
    delete[] RF_MM_A ;
    delete[] RF_MM_B ;
    delete[] RF_JJ_A ;
    delete[] RF_JJ_B ;
    delete[] RF_JM_AB;
    delete[] RF_JM_BA;
    delete[] RF_MJ_AB;
    delete[] RF_MJ_BA;
    delete[] RF_ABBA ;
    delete[] RF_BAAB ;
    delete[] RF_HM_A;
    delete[] RF_MH_A;
    delete[] RF_HJ_A;
    delete[] RF_JH_A;
    delete[] RF_HM_B;
    delete[] RF_MH_B;
    delete[] RF_HJ_B;
    delete[] RF_JH_B;
    delete[] VCAA_AA;
    delete[] VCAA_AB;
    delete[] VCAA_BA;
    delete[] VCAA_BB;
    delete[] CAVA_AA;
    delete[] CAVA_AB;
    delete[] CAVA_BA;
    delete[] CAVA_BB;
    
    return 0;
    
}

int UPT_tensors::calc_EE_2_AV(){
    
    
    
//     int aux_n_ao = aux_n_ao;
    double * VA_A=RI_A->VA_RI_M;
    double * CC_A=RI_A->CC_RI_M;
    double * VC_A=RI_A->VC_RI_M;
    double * CA_A=RI_A->CA_RI_M;
    double * VV_A=RI_A->AA_RI_M;
    
    double * VA_B=RI_B->VA_RI_M;
    double * CC_B=RI_B->CC_RI_M;
    double * VC_B=RI_B->VC_RI_M;
    double * CA_B=RI_B->CA_RI_M;
    double * VV_B=RI_B->AA_RI_M;
    
    
    double * VAAA_AA;
    double * VAAA_AB;
    double * VAAA_BA;
    double * VAAA_BB;
    VAAA_AA = new double[n_v*n_a*n_a*n_a];
    VAAA_AB = new double[n_v*n_a*n_a*n_a];
    VAAA_BA = new double[n_v*n_a*n_a*n_a];
    VAAA_BB = new double[n_v*n_a*n_a*n_a];
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        n_v*n_a,n_a*n_a,aux_n_ao,1.0,
                        VA_A,aux_n_ao,
                        VV_A,aux_n_ao,0.0,
                        VAAA_AA,n_a*n_a);
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        n_v*n_a,n_a*n_a,aux_n_ao,1.0,
                        VA_A,aux_n_ao,
                        VV_B,aux_n_ao,0.0,
                        VAAA_AB,n_a*n_a);
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        n_v*n_a,n_a*n_a,aux_n_ao,1.0,
                        VA_B,aux_n_ao,
                        VV_A,aux_n_ao,0.0,
                        VAAA_BA,n_a*n_a);
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        n_v*n_a,n_a*n_a,aux_n_ao,1.0,
                        VA_B,aux_n_ao,
                        VV_B,aux_n_ao,0.0,
                        VAAA_BB,n_a*n_a);
    
    double * RF_P3_JJ_ABAABA;
    double * RF_P3_JJ_BABBAB;
    double * RF_P3_JM_ABABBB;
    double * RF_P3_JM_BABAAA;
    double * RF_P3_MJ_AAABAB;
    double * RF_P3_MJ_BBBABA;
    double * RF_P3_HM_A;
    double * RF_P3_HM_B;
    double * RF_P3_MH_A;
    double * RF_P3_MH_B;
    double * RF_P3_HJ_ABAB;
    double * RF_P3_HJ_BABA;
    double * RF_P3_JH_ABAB;
    double * RF_P3_JH_BABA;
    
    RF_P3_JJ_ABAABA = new double[n_a*n_a*n_a*n_a*n_a*n_a];
    RF_P3_JJ_BABBAB = new double[n_a*n_a*n_a*n_a*n_a*n_a];
    RF_P3_JM_ABABBB = new double[n_a*n_a*n_a*n_a*n_a*n_a];
    RF_P3_JM_BABAAA = new double[n_a*n_a*n_a*n_a*n_a*n_a];
    RF_P3_MJ_AAABAB = new double[n_a*n_a*n_a*n_a*n_a*n_a];
    RF_P3_MJ_BBBABA = new double[n_a*n_a*n_a*n_a*n_a*n_a];
    
    RF_P3_HM_A = new double[n_a*n_a*n_a*n_a];
    RF_P3_HM_B = new double[n_a*n_a*n_a*n_a];
    RF_P3_MH_A = new double[n_a*n_a*n_a*n_a];
    RF_P3_MH_B = new double[n_a*n_a*n_a*n_a];
    RF_P3_HJ_ABAB = new double[n_a*n_a*n_a*n_a];
    RF_P3_HJ_BABA = new double[n_a*n_a*n_a*n_a];
    RF_P3_JH_ABAB = new double[n_a*n_a*n_a*n_a];
    RF_P3_JH_BABA = new double[n_a*n_a*n_a*n_a];
    
    set_zero_matr(RF_P3_JJ_ABAABA, n_a*n_a*n_a*n_a*n_a*n_a);
    set_zero_matr(RF_P3_JJ_BABBAB, n_a*n_a*n_a*n_a*n_a*n_a);
    set_zero_matr(RF_P3_JM_ABABBB, n_a*n_a*n_a*n_a*n_a*n_a);
    set_zero_matr(RF_P3_JM_BABAAA, n_a*n_a*n_a*n_a*n_a*n_a);
    set_zero_matr(RF_P3_MJ_AAABAB, n_a*n_a*n_a*n_a*n_a*n_a);
    set_zero_matr(RF_P3_MJ_BBBABA, n_a*n_a*n_a*n_a*n_a*n_a);
    
    set_zero_matr(RF_P3_HM_A, n_a*n_a*n_a*n_a);
    set_zero_matr(RF_P3_HM_B, n_a*n_a*n_a*n_a);
    set_zero_matr(RF_P3_MH_A, n_a*n_a*n_a*n_a);
    set_zero_matr(RF_P3_MH_B, n_a*n_a*n_a*n_a);
    set_zero_matr(RF_P3_HJ_ABAB, n_a*n_a*n_a*n_a);
    set_zero_matr(RF_P3_HJ_BABA, n_a*n_a*n_a*n_a);
    set_zero_matr(RF_P3_JH_ABAB, n_a*n_a*n_a*n_a);
    set_zero_matr(RF_P3_JH_BABA, n_a*n_a*n_a*n_a);
    
    
#pragma omp parallel
    {
    //MM == JK
        double V;
        double dE;
        double H1a;
        double H2a;
        
        int u,v;
        int th_id = omp_get_thread_num();
        
        for(int vu=th_id; vu<n_a*n_a; vu+=num_threads){
            v=vu/n_a;
            u=vu%n_a;
            for(int a=0; a<n_v; a++)
            for(int t=0; t<n_a; t++){
                dE=-e_v_A[a]+e_IP_A[t]+e_IP_A[u]-e_EA_A[v];
                if(v==t)dE=-e_v_A[a]+e_IP_A[u];
                if(v==u)dE=-e_v_A[a]+e_IP_A[t];
                dE=dE/(dE*dE+edshift);
                
                for(int y=0; y<n_a; y++)
                for(int x=0; x<n_a; x++)
                for(int w=0; w<n_a; w++){
                    //<aw|xy>-<aw|yx>=(yw|xa)-(xw|ya)
                    V =(VAAA_AA[((a*n_a+t)*n_a+v)*n_a+u]-
                        VAAA_AA[((a*n_a+u)*n_a+v)*n_a+t])*
                       (VAAA_AA[((a*n_a+x)*n_a+w)*n_a+y]-
                        VAAA_AA[((a*n_a+y)*n_a+w)*n_a+x]);
                    
                       RF_P3_JK_A[((((t*n_a+u)*n_a+v)*n_a+y)*n_a+x)*n_a+w]+=V*dE;
                    
                }
            }
            
            for(int a=0; a<n_v; a++)
            for(int t=0; t<n_a; t++){
                dE=-e_v_B[a]+e_IP_B[t]+e_IP_B[u]-e_EA_B[v];
                if(v==t)dE=-e_v_B[a]+e_IP_B[u];
                if(v==u)dE=-e_v_B[a]+e_IP_B[t];
                dE=dE/(dE*dE+edshift);
                
                for(int y=0; y<n_a; y++)
                for(int x=0; x<n_a; x++)
                for(int w=0; w<n_a; w++){
                    //<aw|xy>-<aw|yx>=(yw|xa)-(xw|ya)
                    V =(VAAA_BB[((a*n_a+t)*n_a+v)*n_a+u]-
                        VAAA_BB[((a*n_a+u)*n_a+v)*n_a+t])*
                       (VAAA_BB[((a*n_a+x)*n_a+w)*n_a+y]-
                        VAAA_BB[((a*n_a+y)*n_a+w)*n_a+x]);
                    
                       RF_P3_JK_B[((((t*n_a+u)*n_a+v)*n_a+y)*n_a+x)*n_a+w]+=V*dE;
                    
                }
            }
        
            //JJ
            for(int a=0; a<n_v; a++)
            for(int t=0; t<n_a; t++){
                dE=-e_v_A[a]+e_IP_B[t]+e_IP_A[u]-e_EA_B[v];
                if(v==t)dE=-e_v_A[a]+e_IP_A[u];
                if(v==u)dE=-e_v_A[a]+e_IP_B[t];
                dE=dE/(dE*dE+edshift);
                
                for(int y=0; y<n_a; y++)
                for(int x=0; x<n_a; x++)
                for(int w=0; w<n_a; w++){
                    //<a(a)w(b)|x(a)y(b)>=(yw|xa)
                    V=VAAA_AB[((a*n_a+u)*n_a+v)*n_a+t]*
                      VAAA_AB[((a*n_a+x)*n_a+w)*n_a+y];
                    
                    RF_P3_JJ_BABBAB[((((t*n_a+u)*n_a+v)*n_a+y)*n_a+x)*n_a+w]+=V*dE;//BABBAB
                    
                }
            }
            for(int a=0; a<n_v; a++)
            for(int t=0; t<n_a; t++){
                dE=-e_v_B[a]+e_IP_A[t]+e_IP_B[u]-e_EA_A[v];
                if(v==t)dE=-e_v_B[a]+e_IP_B[u];
                if(v==u)dE=-e_v_B[a]+e_IP_A[t];
                dE=dE/(dE*dE+edshift);
                
                for(int y=0; y<n_a; y++)
                for(int x=0; x<n_a; x++)
                for(int w=0; w<n_a; w++){
                    V=VAAA_BA[((a*n_a+u)*n_a+v)*n_a+t]*
                      VAAA_BA[((a*n_a+x)*n_a+w)*n_a+y];
                    
                    RF_P3_JJ_ABAABA[((((t*n_a+u)*n_a+v)*n_a+y)*n_a+x)*n_a+w]+=V*dE;//ABAABA
                    
                }
            }
            //MJ
            for(int a=0; a<n_v; a++)
            for(int t=0; t<n_a; t++){
                dE=-e_v_A[a]+e_IP_A[t]+e_IP_A[u]-e_EA_A[v];
                if(v==t)dE=-e_v_A[a]+e_IP_A[u];
                if(v==u)dE=-e_v_A[a]+e_IP_A[t];
                dE=dE/(dE*dE+edshift);
                
                for(int y=0; y<n_a; y++)
                for(int x=0; x<n_a; x++)
                for(int w=0; w<n_a; w++){
                    //<a(a)w(b)|x(a)y(b)>=(yw|xa)
                    V=(VAAA_AA[((a*n_a+t)*n_a+v)*n_a+u]-
                       VAAA_AA[((a*n_a+u)*n_a+v)*n_a+t])*
                       VAAA_AB[((a*n_a+x)*n_a+w)*n_a+y];
            
                    RF_P3_MJ_AAABAB[((((t*n_a+u)*n_a+v)*n_a+y)*n_a+x)*n_a+w]+=V*dE;//AAABAB
                    
                }
            }
            for(int a=0; a<n_v; a++)
            for(int t=0; t<n_a; t++){
                dE=-e_v_B[a]+e_IP_B[t]+e_IP_B[u]-e_EA_B[v];
                if(v==t)dE=-e_v_B[a]+e_IP_B[u];
                if(v==u)dE=-e_v_B[a]+e_IP_B[t];
                dE=dE/(dE*dE+edshift);
                
                for(int y=0; y<n_a; y++)
                for(int x=0; x<n_a; x++)
                for(int w=0; w<n_a; w++){
                    //<a(a)w(b)|x(a)y(b)>=(yw|xa)
                    V=(VAAA_BB[((a*n_a+t)*n_a+v)*n_a+u]-
                       VAAA_BB[((a*n_a+u)*n_a+v)*n_a+t])*
                       VAAA_BA[((a*n_a+x)*n_a+w)*n_a+y];
            
                    RF_P3_MJ_BBBABA[((((t*n_a+u)*n_a+v)*n_a+y)*n_a+x)*n_a+w]+=V*dE;//BBBABA
                    
                }
            }
            //JM
            for(int a=0; a<n_v; a++)
            for(int t=0; t<n_a; t++){
                dE=-e_v_A[a]+e_IP_B[t]+e_IP_A[u]-e_EA_B[v];
                if(v==t)dE=-e_v_A[a]+e_IP_A[u];
                if(v==u)dE=-e_v_A[a]+e_IP_B[t];
                dE=dE/(dE*dE+edshift);
                
                for(int y=0; y<n_a; y++)
                for(int x=0; x<n_a; x++)
                for(int w=0; w<n_a; w++){
                    //<aw|xy>-<aw|yx>=(yw|xa)-(xw|ya)
                    V= VAAA_AB[((a*n_a+u)*n_a+v)*n_a+t]*
                      (VAAA_BB[((a*n_a+x)*n_a+w)*n_a+y]-
                       VAAA_BB[((a*n_a+y)*n_a+w)*n_a+x]);
                    
                    RF_P3_JM_BABAAA[((((t*n_a+u)*n_a+v)*n_a+y)*n_a+x)*n_a+w]+=V*dE;//BABAAA
                    
                }
            }
            for(int a=0; a<n_v; a++)
            for(int t=0; t<n_a; t++){
                dE=-e_v_B[a]+e_IP_A[t]+e_IP_B[u]-e_EA_A[v];
                if(v==t)dE=-e_v_B[a]+e_IP_B[u];
                if(v==u)dE=-e_v_B[a]+e_IP_A[t];
                dE=dE/(dE*dE+edshift);
                
                for(int y=0; y<n_a; y++)
                for(int x=0; x<n_a; x++)
                for(int w=0; w<n_a; w++){
                    //<aw|xy>-<aw|yx>=(yw|xa)-(xw|ya)
                    V= VAAA_BA[((a*n_a+u)*n_a+v)*n_a+t]*
                      (VAAA_AA[((a*n_a+x)*n_a+w)*n_a+y]-
                       VAAA_AA[((a*n_a+y)*n_a+w)*n_a+x]);
                    
                    RF_P3_JM_ABABBB[((((t*n_a+u)*n_a+v)*n_a+y)*n_a+x)*n_a+w]+=V*dE;//ABABBB
                    
                }
            }
        }    
        //H
        for(int a=0; a<n_v; a++)
        for(int t=th_id; t<n_a; t+=num_threads){
            dE=-e_v_A[a]+e_IP_A[t];
            dE=dE/(dE*dE+edshift);
            
            for(int w=0; w<n_a; w++){
                H1a=H_AV_A[t*n_v+a];
                H2a=H_AV_A[w*n_v+a];
                            
                V=H1a*H2a;
                RF_PH_A[t*n_a+w]+=V*dE;
                
            }
        }
        for(int a=0; a<n_v; a++)
        for(int t=th_id; t<n_a; t+=num_threads){
            dE=-e_v_B[a]+e_IP_B[t];
            dE=dE/(dE*dE+edshift);
            
            for(int w=0; w<n_a; w++){
                H1a=H_AV_B[t*n_v+a];
                H2a=H_AV_B[w*n_v+a];
                            
                V=H1a*H2a;
                RF_PH_B[t*n_a+w]+=V*dE;
                
            }
        }
        //HM
        for(int a=0; a<n_v; a++)
        for(int t=th_id; t<n_a; t+=num_threads){
            dE=-e_v_A[a]+e_IP_A[t];
            dE=dE/(dE*dE+edshift);
            
            for(int y=0; y<n_a; y++)
            for(int x=0; x<n_a; x++)
            for(int w=0; w<n_a; w++){
                H1a=H_AV_A[t*n_v+a];
     
                //<aw|xy>-<aw|yx>=(yw|xa)-(xw|ya)
                V=H1a*
                  (VAAA_AA[((a*n_a+x)*n_a+w)*n_a+y]-
                   VAAA_AA[((a*n_a+y)*n_a+w)*n_a+x]);
                
                RF_P3_HM_A[((t*n_a+y)*n_a+x)*n_a+w]+=V*dE;
                
            }
        }
        for(int a=0; a<n_v; a++)
        for(int t=th_id; t<n_a; t+=num_threads){
            dE=-e_v_B[a]+e_IP_B[t];
            dE=dE/(dE*dE+edshift);
            
            for(int y=0; y<n_a; y++)
            for(int x=0; x<n_a; x++)
            for(int w=0; w<n_a; w++){
                H1a=H_AV_B[t*n_v+a];
     
                //<aw|xy>-<aw|yx>=(yw|xa)-(xw|ya)
                V=H1a*
                  (VAAA_BB[((a*n_a+x)*n_a+w)*n_a+y]-
                   VAAA_BB[((a*n_a+y)*n_a+w)*n_a+x]);
                
                RF_P3_HM_B[((t*n_a+y)*n_a+x)*n_a+w]+=V*dE;
                
            }
        }
        //MH
        for(int vu=th_id; vu<n_a*n_a; vu+=num_threads){
            v=vu/n_a;
            u=vu%n_a;
            for(int a=0; a<n_v; a++)
            for(int t=0; t<n_a; t++){
                dE=-e_v_A[a]+e_IP_A[t]+e_IP_A[u]-e_EA_A[v];
                if(v==t)dE=-e_v_A[a]+e_IP_A[u];
                if(v==u)dE=-e_v_A[a]+e_IP_A[t];
                dE=dE/(dE*dE+edshift);
                
                for(int w=0; w<n_a; w++){
                    H2a=H_AV_A[w*n_v+a];
                    //<aw|xy>-<aw|yx>=(yw|xa)-(xw|ya)
                    V=(VAAA_AA[((a*n_a+t)*n_a+v)*n_a+u]-
                       VAAA_AA[((a*n_a+u)*n_a+v)*n_a+t])*
                       H2a;
                    
                    RF_P3_MH_A[((t*n_a+u)*n_a+v)*n_a+w]+=V*dE;
                }
            }
        }
        for(int vu=th_id; vu<n_a*n_a; vu+=num_threads){
            v=vu/n_a;
            u=vu%n_a;
            for(int a=0; a<n_v; a++)
            for(int t=0; t<n_a; t++){
                dE=-e_v_B[a]+e_IP_B[t]+e_IP_B[u]-e_EA_B[v];
                if(v==t)dE=-e_v_B[a]+e_IP_B[u];
                if(v==u)dE=-e_v_B[a]+e_IP_B[t];
                dE=dE/(dE*dE+edshift);
                
                for(int w=0; w<n_a; w++){
                    H2a=H_AV_B[w*n_v+a];
                    //<aw|xy>-<aw|yx>=(yw|xa)-(xw|ya)
                    V=(VAAA_BB[((a*n_a+t)*n_a+v)*n_a+u]-
                       VAAA_BB[((a*n_a+u)*n_a+v)*n_a+t])*
                       H2a;
                    
                    RF_P3_MH_B[((t*n_a+u)*n_a+v)*n_a+w]+=V*dE;
                }
            }
        }
        //HJ
        for(int a=0; a<n_v; a++)
        for(int t=th_id; t<n_a; t+=num_threads){
            dE=-e_v_A[a]+e_IP_A[t];
            dE=dE/(dE*dE+edshift);
            
            for(int y=0; y<n_a; y++)
            for(int x=0; x<n_a; x++)
            for(int w=0; w<n_a; w++){
                H1a=H_AV_A[t*n_v+a];
                //<a(a)w(b)|x(a)y(b)>=(yw|xa)
                V=H1a*
                  VAAA_AB[((a*n_a+x)*n_a+w)*n_a+y];
                
                RF_P3_HJ_ABAB[((t*n_a+y)*n_a+x)*n_a+w]+=V*dE;//ABAB
            }
        }
        for(int a=0; a<n_v; a++)
        for(int t=th_id; t<n_a; t+=num_threads){
            dE=-e_v_B[a]+e_IP_B[t];
            dE=dE/(dE*dE+edshift);
            
            for(int y=0; y<n_a; y++)
            for(int x=0; x<n_a; x++)
            for(int w=0; w<n_a; w++){
                H1a=H_AV_B[t*n_v+a];
                //<a(a)w(b)|x(a)y(b)>=(yw|xa)
                V=H1a*
                  VAAA_BA[((a*n_a+x)*n_a+w)*n_a+y];
                
                RF_P3_HJ_BABA[((t*n_a+y)*n_a+x)*n_a+w]+=V*dE;//BABA
            }
        }
        //JH
        for(int vu=th_id; vu<n_a*n_a; vu+=num_threads){
            v=vu/n_a;
            u=vu%n_a;
            for(int a=0; a<n_v; a++)
            for(int t=0; t<n_a; t++){
                dE=-e_v_A[a]+e_IP_B[t]+e_IP_A[u]-e_EA_B[v];
                if(v==t)dE=-e_v_A[a]+e_IP_A[u];
                if(v==u)dE=-e_v_A[a]+e_IP_B[t];
                dE=dE/(dE*dE+edshift);
                
                for(int w=0; w<n_a; w++){
                    H2a=H_AV_A[w*n_v+a];     
                    //<a(a)w(b)|x(a)y(b)>=(yw|xa)
                    V=VAAA_AB[((a*n_a+u)*n_a+v)*n_a+t]*
                      H2a;

                    RF_P3_JH_BABA[((t*n_a+u)*n_a+v)*n_a+w]+=V*dE;//BABA
                }
            }
        }
        for(int vu=th_id; vu<n_a*n_a; vu+=num_threads){
            v=vu/n_a;
            u=vu%n_a;
            for(int a=0; a<n_v; a++)
            for(int t=0; t<n_a; t++){
                dE=-e_v_B[a]+e_IP_A[t]+e_IP_B[u]-e_EA_A[v];
                if(v==t)dE=-e_v_B[a]+e_IP_B[u];
                if(v==u)dE=-e_v_B[a]+e_IP_A[t];
                dE=dE/(dE*dE+edshift);
                
                for(int w=0; w<n_a; w++){
                    H2a=H_AV_B[w*n_v+a];     
                    //<a(a)w(b)|x(a)y(b)>=(yw|xa)
                    V=VAAA_BA[((a*n_a+u)*n_a+v)*n_a+t]*
                      H2a;

                    RF_P3_JH_ABAB[((t*n_a+u)*n_a+v)*n_a+w]+=V*dE;//ABAB
                }
            }
        }
    }
    //RF_PV_JK
    for(int w=0; w<n_a; w++)
    for(int t=0; t<n_a; t++)
    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++)
        RF_PV_JK_A[((w*n_a+t)*n_a+v)*n_a+u]+= RF_P3_HM_A[((t*n_a+u)*n_a+v)*n_a+w]
                                             -RF_P3_HM_A[((w*n_a+u)*n_a+v)*n_a+t]
                                             -RF_P3_MH_A[((t*n_a+w)*n_a+v)*n_a+u]
                                             +RF_P3_MH_A[((t*n_a+w)*n_a+u)*n_a+v];
    for(int w=0; w<n_a; w++)
    for(int t=0; t<n_a; t++)
    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++)
        RF_PV_JK_B[((w*n_a+t)*n_a+v)*n_a+u]+= RF_P3_HM_B[((t*n_a+u)*n_a+v)*n_a+w]
                                             -RF_P3_HM_B[((w*n_a+u)*n_a+v)*n_a+t]
                                             -RF_P3_MH_B[((t*n_a+w)*n_a+v)*n_a+u]
                                             +RF_P3_MH_B[((t*n_a+w)*n_a+u)*n_a+v];
    for(int w=0; w<n_a; w++)
    for(int t=0; t<n_a; t++)
    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++){
        //         A     B       A     B
        RF_PV_AB[((t*n_a+u)*n_a+v)*n_a+w]+= RF_P3_HJ_ABAB[((t*n_a+w)*n_a+v)*n_a+u]//ABAB
                                           +RF_P3_JH_BABA[((u*n_a+t)*n_a+w)*n_a+v]//BABA
                                           +RF_P3_HJ_BABA[((u*n_a+v)*n_a+w)*n_a+t]//BABA
                                           +RF_P3_JH_ABAB[((t*n_a+u)*n_a+v)*n_a+w]//ABAB
                                           ;//2nd!!!!!!!!!!!
    
        for(int x=0; x<n_a; x++)
            //         A     B       A     B
            RF_PV_AB[((t*n_a+u)*n_a+v)*n_a+w] += RF_P3_JJ_ABAABA[((((t*n_a+u)*n_a+x)*n_a+v)*n_a+w)*n_a+x]
                                                +RF_P3_JJ_BABBAB[((((u*n_a+t)*n_a+x)*n_a+w)*n_a+v)*n_a+x];
    }
    for(int w=0; w<n_a; w++)
    for(int t=0; t<n_a; t++)
    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++){
        //         B     A       B     A
        RF_PV_BA[((t*n_a+u)*n_a+v)*n_a+w]+= RF_P3_HJ_BABA[((t*n_a+w)*n_a+v)*n_a+u]//BABA
                                           +RF_P3_JH_ABAB[((u*n_a+t)*n_a+w)*n_a+v]//ABAB
                                           +RF_P3_HJ_ABAB[((u*n_a+v)*n_a+w)*n_a+t]//ABAB
                                           +RF_P3_JH_BABA[((t*n_a+u)*n_a+v)*n_a+w]//BABA
                                           ;//2nd!!!!!!!!!!!
    
        for(int x=0; x<n_a; x++)
            //         B     A       B     A
            RF_PV_BA[((t*n_a+u)*n_a+v)*n_a+w] += RF_P3_JJ_BABBAB[((((t*n_a+u)*n_a+x)*n_a+v)*n_a+w)*n_a+x]
                                                +RF_P3_JJ_ABAABA[((((u*n_a+t)*n_a+x)*n_a+w)*n_a+v)*n_a+x];
    }
    
    
    //RF_P3_AB - order differs with the original XMCQDPT code (res_fit.cpp)
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++)
    for(int x=0; x<n_a; x++)
    for(int v=0; v<n_a; v++)
    for(int w=0; w<n_a; w++)
    for(int y=0; y<n_a; y++)
        //           A     A      B      A      A      B
        RF_P3_AB[((((t*n_a+u)*n_a+x)*n_a+v)*n_a+w)*n_a+y]+=-RF_P3_MJ_AAABAB[((((t*n_a+u)*n_a+v)*n_a+y)*n_a+w)*n_a+x]//AAABAB
                                                           +RF_P3_MJ_AAABAB[((((t*n_a+u)*n_a+w)*n_a+y)*n_a+v)*n_a+x]//AAABAB
                                                           +RF_P3_JM_BABAAA[((((x*n_a+t)*n_a+y)*n_a+w)*n_a+v)*n_a+u]//BABAAA
                                                           -RF_P3_JM_BABAAA[((((x*n_a+u)*n_a+y)*n_a+w)*n_a+v)*n_a+t]//BABAAA
                                                           +RF_P3_JJ_ABAABA[((((t*n_a+x)*n_a+v)*n_a+w)*n_a+y)*n_a+u]//ABAABA
                                                           -RF_P3_JJ_ABAABA[((((t*n_a+x)*n_a+w)*n_a+v)*n_a+y)*n_a+u]//ABAABA
                                                           -RF_P3_JJ_ABAABA[((((u*n_a+x)*n_a+v)*n_a+w)*n_a+y)*n_a+t]//ABAABA
                                                           +RF_P3_JJ_ABAABA[((((u*n_a+x)*n_a+w)*n_a+v)*n_a+y)*n_a+t]//ABAABA
                                                           ;
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++)
    for(int x=0; x<n_a; x++)
    for(int v=0; v<n_a; v++)
    for(int w=0; w<n_a; w++)
    for(int y=0; y<n_a; y++)
        //           B     B      A      B      B      A
        RF_P3_BA[((((t*n_a+u)*n_a+x)*n_a+v)*n_a+w)*n_a+y]+=-RF_P3_MJ_BBBABA[((((t*n_a+u)*n_a+v)*n_a+y)*n_a+w)*n_a+x]//BBBABA
                                                           +RF_P3_MJ_BBBABA[((((t*n_a+u)*n_a+w)*n_a+y)*n_a+v)*n_a+x]//BBBABA
                                                           +RF_P3_JM_ABABBB[((((x*n_a+t)*n_a+y)*n_a+w)*n_a+v)*n_a+u]//ABABBB
                                                           -RF_P3_JM_ABABBB[((((x*n_a+u)*n_a+y)*n_a+w)*n_a+v)*n_a+t]//ABABBB
                                                           +RF_P3_JJ_BABBAB[((((t*n_a+x)*n_a+v)*n_a+w)*n_a+y)*n_a+u]//BABBAB
                                                           -RF_P3_JJ_BABBAB[((((t*n_a+x)*n_a+w)*n_a+v)*n_a+y)*n_a+u]//BABBAB
                                                           -RF_P3_JJ_BABBAB[((((u*n_a+x)*n_a+v)*n_a+w)*n_a+y)*n_a+t]//BABBAB
                                                           +RF_P3_JJ_BABBAB[((((u*n_a+x)*n_a+w)*n_a+v)*n_a+y)*n_a+t]//BABBAB
                                                           ;
                                    
    
    delete[] VAAA_AA;
    delete[] VAAA_AB;
    delete[] VAAA_BA;
    delete[] VAAA_BB;
    
    delete[] RF_P3_JJ_ABAABA;
    delete[] RF_P3_JJ_BABBAB;
    delete[] RF_P3_JM_ABABBB;
    delete[] RF_P3_JM_BABAAA;
    delete[] RF_P3_MJ_AAABAB;
    delete[] RF_P3_MJ_BBBABA;
    delete[] RF_P3_HM_A;
    delete[] RF_P3_HM_B;
    delete[] RF_P3_MH_A;
    delete[] RF_P3_MH_B;
    delete[] RF_P3_HJ_ABAB;
    delete[] RF_P3_HJ_BABA;
    delete[] RF_P3_JH_ABAB;
    delete[] RF_P3_JH_BABA;
    
    return 0;
    
}

int UPT_tensors::calc_EE_2_CA(){
    
    double dE;
    double H1a;
    double H2a;
    double V;
    
    
    double * CC_A=RI_A->CC_RI_M;
    double * CA_A=RI_A->CA_RI_M;
    double * AA_A=RI_A->AA_RI_M;
    double * CC_B=RI_B->CC_RI_M;
    double * CA_B=RI_B->CA_RI_M;
    double * AA_B=RI_B->AA_RI_M;

    
    double * CAAA_AA;
    double * CAAA_AB;
    double * CAAA_BA;
    double * CAAA_BB;
    CAAA_AA = new double[n_c*n_a*n_a*n_a];
    CAAA_AB = new double[n_c*n_a*n_a*n_a];
    CAAA_BA = new double[n_c*n_a*n_a*n_a];
    CAAA_BB = new double[n_c*n_a*n_a*n_a];
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        n_c*n_a,n_a*n_a,aux_n_ao,1.0,
                        CA_A,aux_n_ao,
                        AA_A,aux_n_ao,0.0,
                        CAAA_AA,n_a*n_a);
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        n_c*n_a,n_a*n_a,aux_n_ao,1.0,
                        CA_A,aux_n_ao,
                        AA_B,aux_n_ao,0.0,
                        CAAA_AB,n_a*n_a);
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        n_c*n_a,n_a*n_a,aux_n_ao,1.0,
                        CA_B,aux_n_ao,
                        AA_A,aux_n_ao,0.0,
                        CAAA_BA,n_a*n_a);
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        n_c*n_a,n_a*n_a,aux_n_ao,1.0,
                        CA_B,aux_n_ao,
                        AA_B,aux_n_ao,0.0,
                        CAAA_BB,n_a*n_a);
    
    double * RF_P3_VV_AAAAAA;
    double * RF_P3_VV_AAABBA;
    double * RF_P3_VV_AABAAB;
    double * RF_P3_VV_AABBBB;
    double * RF_P3_VV_BBAAAA;
    double * RF_P3_VV_BBABBA;
    double * RF_P3_VV_BBBAAB;
    double * RF_P3_VV_BBBBBB;
    double * RF_P3_VH_AAAA;
    double * RF_P3_VH_AABB;
    double * RF_P3_VH_BBAA;
    double * RF_P3_VH_BBBB;
    double * RF_P3_HV_AAAA;
    double * RF_P3_HV_ABBA;
    double * RF_P3_HV_BAAB;
    double * RF_P3_HV_BBBB;
    double * RF_P3_HH_A;
    double * RF_P3_HH_B;
    
    RF_P3_VV_AAAAAA = new double[n_a*n_a*n_a*n_a*n_a*n_a];
    RF_P3_VV_AAABBA = new double[n_a*n_a*n_a*n_a*n_a*n_a];
    RF_P3_VV_AABAAB = new double[n_a*n_a*n_a*n_a*n_a*n_a];
    RF_P3_VV_AABBBB = new double[n_a*n_a*n_a*n_a*n_a*n_a];
    RF_P3_VV_BBAAAA = new double[n_a*n_a*n_a*n_a*n_a*n_a];
    RF_P3_VV_BBABBA = new double[n_a*n_a*n_a*n_a*n_a*n_a];
    RF_P3_VV_BBBAAB = new double[n_a*n_a*n_a*n_a*n_a*n_a];
    RF_P3_VV_BBBBBB = new double[n_a*n_a*n_a*n_a*n_a*n_a];
    RF_P3_VH_AAAA   = new double[n_a*n_a*n_a*n_a        ];
    RF_P3_VH_AABB   = new double[n_a*n_a*n_a*n_a        ];
    RF_P3_VH_BBAA   = new double[n_a*n_a*n_a*n_a        ];
    RF_P3_VH_BBBB   = new double[n_a*n_a*n_a*n_a        ];
    RF_P3_HV_AAAA   = new double[n_a*n_a*n_a*n_a        ];
    RF_P3_HV_ABBA   = new double[n_a*n_a*n_a*n_a        ];
    RF_P3_HV_BAAB   = new double[n_a*n_a*n_a*n_a        ];
    RF_P3_HV_BBBB   = new double[n_a*n_a*n_a*n_a        ];
    RF_P3_HH_A      = new double[n_a*n_a                ];
    RF_P3_HH_B      = new double[n_a*n_a                ];
    
    set_zero_matr(RF_P3_VV_AAAAAA, n_a*n_a*n_a*n_a*n_a*n_a);
    set_zero_matr(RF_P3_VV_AAABBA, n_a*n_a*n_a*n_a*n_a*n_a);
    set_zero_matr(RF_P3_VV_AABAAB, n_a*n_a*n_a*n_a*n_a*n_a);
    set_zero_matr(RF_P3_VV_AABBBB, n_a*n_a*n_a*n_a*n_a*n_a);
    set_zero_matr(RF_P3_VV_BBAAAA, n_a*n_a*n_a*n_a*n_a*n_a);
    set_zero_matr(RF_P3_VV_BBABBA, n_a*n_a*n_a*n_a*n_a*n_a);
    set_zero_matr(RF_P3_VV_BBBAAB, n_a*n_a*n_a*n_a*n_a*n_a);
    set_zero_matr(RF_P3_VV_BBBBBB, n_a*n_a*n_a*n_a*n_a*n_a);
    set_zero_matr(RF_P3_VH_AAAA  , n_a*n_a*n_a*n_a        );
    set_zero_matr(RF_P3_VH_AABB  , n_a*n_a*n_a*n_a        );
    set_zero_matr(RF_P3_VH_BBAA  , n_a*n_a*n_a*n_a        );
    set_zero_matr(RF_P3_VH_BBBB  , n_a*n_a*n_a*n_a        );
    set_zero_matr(RF_P3_HV_AAAA  , n_a*n_a*n_a*n_a        );
    set_zero_matr(RF_P3_HV_ABBA  , n_a*n_a*n_a*n_a        );
    set_zero_matr(RF_P3_HV_BAAB  , n_a*n_a*n_a*n_a        );
    set_zero_matr(RF_P3_HV_BBBB  , n_a*n_a*n_a*n_a        );
    set_zero_matr(RF_P3_HH_A     , n_a*n_a                );
    set_zero_matr(RF_P3_HH_B     , n_a*n_a                );
    
    
    //VV
    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++)
    for(int i=0; i<n_c; i++)
    for(int t=0; t<n_a; t++){
        dE= e_c_A[i]+e_IP_A[t]-e_EA_A[u]-e_EA_A[v];
        if(t==u)dE= e_c_A[i]-e_EA_A[v];
        if(t==v)dE= e_c_A[i]-e_EA_A[u];
        dE=dE/(dE*dE+edshift);
        
        for(int y=0; y<n_a; y++)
        for(int x=0; x<n_a; x++)
        for(int w=0; w<n_a; w++){
            
            V=(CAAA_AA[((i*n_a+v)*n_a+u)*n_a+t]*
               CAAA_AA[((i*n_a+w)*n_a+x)*n_a+y]);
            
            RF_P3_VV_AAAAAA[((((t*n_a+u)*n_a+v)*n_a+y)*n_a+x)*n_a+w]+=V*dE;
            
        }
    }
    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++)
    for(int i=0; i<n_c; i++)
    for(int t=0; t<n_a; t++){
        dE= e_c_A[i]+e_IP_B[t]-e_EA_B[u]-e_EA_A[v];
        if(t==u)dE= e_c_A[i]-e_EA_A[v];
        if(t==v)dE= e_c_A[i]-e_EA_B[u];
        dE=dE/(dE*dE+edshift);
        
        for(int y=0; y<n_a; y++)
        for(int x=0; x<n_a; x++)
        for(int w=0; w<n_a; w++){
            
            V=(CAAA_AB[((i*n_a+v)*n_a+u)*n_a+t]*
               CAAA_AA[((i*n_a+w)*n_a+x)*n_a+y]);
            
            RF_P3_VV_BBAAAA[((((t*n_a+u)*n_a+v)*n_a+y)*n_a+x)*n_a+w]+=V*dE;
            
        }
    }
    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++)
    for(int i=0; i<n_c; i++)
    for(int t=0; t<n_a; t++){
        dE= e_c_A[i]+e_IP_A[t]-e_EA_A[u]-e_EA_A[v];
        if(t==u)dE= e_c_A[i]-e_EA_A[v];
        if(t==v)dE= e_c_A[i]-e_EA_A[u];
        dE=dE/(dE*dE+edshift);
        
        for(int y=0; y<n_a; y++)
        for(int x=0; x<n_a; x++)
        for(int w=0; w<n_a; w++){
            
            V=(CAAA_AA[((i*n_a+v)*n_a+u)*n_a+t]*
               CAAA_AB[((i*n_a+w)*n_a+x)*n_a+y]);
            
            RF_P3_VV_AAABBA[((((t*n_a+u)*n_a+v)*n_a+y)*n_a+x)*n_a+w]+=V*dE;
            
        }
    }
    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++)
    for(int i=0; i<n_c; i++)
    for(int t=0; t<n_a; t++){
        dE= e_c_A[i]+e_IP_B[t]-e_EA_B[u]-e_EA_A[v];
        if(t==u)dE= e_c_A[i]-e_EA_A[v];
        if(t==v)dE= e_c_A[i]-e_EA_B[u];
        dE=dE/(dE*dE+edshift);
        
        for(int y=0; y<n_a; y++)
        for(int x=0; x<n_a; x++)
        for(int w=0; w<n_a; w++){
            
            V=(CAAA_AB[((i*n_a+v)*n_a+u)*n_a+t]*
               CAAA_AB[((i*n_a+w)*n_a+x)*n_a+y]);
            
            RF_P3_VV_BBABBA[((((t*n_a+u)*n_a+v)*n_a+y)*n_a+x)*n_a+w]+=V*dE;
            
        }
    }
    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++)
    for(int i=0; i<n_c; i++)
    for(int t=0; t<n_a; t++){
        dE= e_c_B[i]+e_IP_A[t]-e_EA_A[u]-e_EA_B[v];
        if(t==u)dE= e_c_B[i]-e_EA_B[v];
        if(t==v)dE= e_c_B[i]-e_EA_A[u];
        dE=dE/(dE*dE+edshift);
        
        for(int y=0; y<n_a; y++)
        for(int x=0; x<n_a; x++)
        for(int w=0; w<n_a; w++){
            
            V=(CAAA_BA[((i*n_a+v)*n_a+u)*n_a+t]*
               CAAA_BA[((i*n_a+w)*n_a+x)*n_a+y]);
            
            RF_P3_VV_AABAAB[((((t*n_a+u)*n_a+v)*n_a+y)*n_a+x)*n_a+w]+=V*dE;
            
        }
    }
    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++)
    for(int i=0; i<n_c; i++)
    for(int t=0; t<n_a; t++){
        dE= e_c_B[i]+e_IP_B[t]-e_EA_B[u]-e_EA_B[v];
        if(t==u)dE= e_c_B[i]-e_EA_B[v];
        if(t==v)dE= e_c_B[i]-e_EA_B[u];
        dE=dE/(dE*dE+edshift);
        
        for(int y=0; y<n_a; y++)
        for(int x=0; x<n_a; x++)
        for(int w=0; w<n_a; w++){
            
            V=(CAAA_BB[((i*n_a+v)*n_a+u)*n_a+t]*
               CAAA_BA[((i*n_a+w)*n_a+x)*n_a+y]);
            
            RF_P3_VV_BBBAAB[((((t*n_a+u)*n_a+v)*n_a+y)*n_a+x)*n_a+w]+=V*dE;
            
        }
    }
    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++)
    for(int i=0; i<n_c; i++)
    for(int t=0; t<n_a; t++){
        dE= e_c_B[i]+e_IP_A[t]-e_EA_A[u]-e_EA_B[v];
        if(t==u)dE= e_c_B[i]-e_EA_B[v];
        if(t==v)dE= e_c_B[i]-e_EA_A[u];
        dE=dE/(dE*dE+edshift);
        
        for(int y=0; y<n_a; y++)
        for(int x=0; x<n_a; x++)
        for(int w=0; w<n_a; w++){
            
            V=(CAAA_BA[((i*n_a+v)*n_a+u)*n_a+t]*
               CAAA_BB[((i*n_a+w)*n_a+x)*n_a+y]);
            
            RF_P3_VV_AABBBB[((((t*n_a+u)*n_a+v)*n_a+y)*n_a+x)*n_a+w]+=V*dE;
            
        }
    }
    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++)
    for(int i=0; i<n_c; i++)
    for(int t=0; t<n_a; t++){
        dE= e_c_B[i]+e_IP_B[t]-e_EA_B[u]-e_EA_B[v];
        if(t==u)dE= e_c_B[i]-e_EA_B[v];
        if(t==v)dE= e_c_B[i]-e_EA_B[u];
        dE=dE/(dE*dE+edshift);
        
        for(int y=0; y<n_a; y++)
        for(int x=0; x<n_a; x++)
        for(int w=0; w<n_a; w++){
            
            V=(CAAA_BB[((i*n_a+v)*n_a+u)*n_a+t]*
               CAAA_BB[((i*n_a+w)*n_a+x)*n_a+y]);
            
            RF_P3_VV_BBBBBB[((((t*n_a+u)*n_a+v)*n_a+y)*n_a+x)*n_a+w]+=V*dE;
            
        }
    }
    

    
    
    
    //VH
    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++)
    for(int i=0; i<n_c; i++)
    for(int t=0; t<n_a; t++){
        dE= e_c_A[i]+e_IP_A[t]-e_EA_A[u]-e_EA_A[v];
        if(t==u)dE= e_c_A[i]-e_EA_A[v];
        if(t==v)dE= e_c_A[i]-e_EA_A[u];
        dE=dE/(dE*dE+edshift);
        for(int w=0; w<n_a; w++){
            H2a=H_CA_A[i*n_a+w];
            V=(CAAA_AA[((i*n_a+v)*n_a+u)*n_a+t]*
               H2a);
            RF_P3_VH_AAAA[((t*n_a+u)*n_a+v)*n_a+w]+=V*dE;
        }
    }

    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++)
    for(int i=0; i<n_c; i++)
    for(int t=0; t<n_a; t++){
        dE= e_c_A[i]+e_IP_B[t]-e_EA_B[u]-e_EA_A[v];
        if(t==u)dE= e_c_A[i]-e_EA_A[v];
        if(t==v)dE= e_c_A[i]-e_EA_B[u];
        dE=dE/(dE*dE+edshift);
        for(int w=0; w<n_a; w++){
            H2a=H_CA_A[i*n_a+w];
            V=(CAAA_AB[((i*n_a+v)*n_a+u)*n_a+t]*
               H2a);
            RF_P3_VH_BBAA[((t*n_a+u)*n_a+v)*n_a+w]+=V*dE;
        }
    }
    
    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++)
    for(int i=0; i<n_c; i++)
    for(int t=0; t<n_a; t++){
        dE= e_c_B[i]+e_IP_A[t]-e_EA_A[u]-e_EA_B[v];
        if(t==u)dE= e_c_B[i]-e_EA_B[v];
        if(t==v)dE= e_c_B[i]-e_EA_A[u];
        dE=dE/(dE*dE+edshift);
        for(int w=0; w<n_a; w++){
            H2a=H_CA_B[i*n_a+w];
            V=(CAAA_BA[((i*n_a+v)*n_a+u)*n_a+t]*
               H2a);
            RF_P3_VH_AABB[((t*n_a+u)*n_a+v)*n_a+w]+=V*dE;
        }
    }
    
    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++)
    for(int i=0; i<n_c; i++)
    for(int t=0; t<n_a; t++){
        dE= e_c_B[i]+e_IP_B[t]-e_EA_B[u]-e_EA_B[v];
        if(t==u)dE= e_c_B[i]-e_EA_B[v];
        if(t==v)dE= e_c_B[i]-e_EA_B[u];
        dE=dE/(dE*dE+edshift);
        for(int w=0; w<n_a; w++){
            H2a=H_CA_B[i*n_a+w];
            V=(CAAA_BB[((i*n_a+v)*n_a+u)*n_a+t]*
               H2a);
            RF_P3_VH_BBBB[((t*n_a+u)*n_a+v)*n_a+w]+=V*dE;
        }
    }

    //HV
    for(int i=0; i<n_c; i++)
    for(int t=0; t<n_a; t++){
        dE= e_c_A[i]-e_EA_A[t];
        dE=dE/(dE*dE+edshift);
        
        for(int y=0; y<n_a; y++)
        for(int x=0; x<n_a; x++)
        for(int w=0; w<n_a; w++){
            H1a=H_CA_A[i*n_a+t];
            V=(H1a*
               CAAA_AA[((i*n_a+w)*n_a+x)*n_a+y]);
            RF_P3_HV_AAAA[((t*n_a+y)*n_a+x)*n_a+w]+=V*dE;
        }
    }
    for(int i=0; i<n_c; i++)
    for(int t=0; t<n_a; t++){
        dE= e_c_A[i]-e_EA_A[t];
        dE=dE/(dE*dE+edshift);
        
        for(int y=0; y<n_a; y++)
        for(int x=0; x<n_a; x++)
        for(int w=0; w<n_a; w++){
            H1a=H_CA_A[i*n_a+t];
            V=(H1a*
               CAAA_AB[((i*n_a+w)*n_a+x)*n_a+y]);
            RF_P3_HV_ABBA[((t*n_a+y)*n_a+x)*n_a+w]+=V*dE;
        }
    }
    for(int i=0; i<n_c; i++)
    for(int t=0; t<n_a; t++){
        dE= e_c_B[i]-e_EA_B[t];
        dE=dE/(dE*dE+edshift);
        
        for(int y=0; y<n_a; y++)
        for(int x=0; x<n_a; x++)
        for(int w=0; w<n_a; w++){
            H1a=H_CA_B[i*n_a+t];
            V=(H1a*
               CAAA_BA[((i*n_a+w)*n_a+x)*n_a+y]);
            RF_P3_HV_BAAB[((t*n_a+y)*n_a+x)*n_a+w]+=V*dE;
        }
    }
    for(int i=0; i<n_c; i++)
    for(int t=0; t<n_a; t++){
        dE= e_c_B[i]-e_EA_B[t];
        dE=dE/(dE*dE+edshift);
        
        for(int y=0; y<n_a; y++)
        for(int x=0; x<n_a; x++)
        for(int w=0; w<n_a; w++){
            H1a=H_CA_B[i*n_a+t];
            V=(H1a*
               CAAA_BB[((i*n_a+w)*n_a+x)*n_a+y]);
            RF_P3_HV_BBBB[((t*n_a+y)*n_a+x)*n_a+w]+=V*dE;
        }
    }
    
    //HH
    for(int i=0; i<n_c; i++)
    for(int t=0; t<n_a; t++){
        dE= e_c_A[i]-e_EA_A[t];
        dE=dE/(dE*dE+edshift);
        
        for(int u=0; u<n_a; u++){
            H1a=H_CA_A[i*n_a+t];
            H2a=H_CA_A[i*n_a+u];
            V=H1a*H2a;
            RF_P3_HH_A[t*n_a+u]+=V*dE;
        }
    }
    for(int i=0; i<n_c; i++)
    for(int t=0; t<n_a; t++){
        dE= e_c_B[i]-e_EA_B[t];
        dE=dE/(dE*dE+edshift);
        
        for(int u=0; u<n_a; u++){
            H1a=H_CA_B[i*n_a+t];
            H2a=H_CA_B[i*n_a+u];
            V=H1a*H2a;
            RF_P3_HH_B[t*n_a+u]+=V*dE;
        }
    }
    
    
    for(int t=0; t<n_a; t++)
    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++)
    for(int w=0; w<n_a; w++)
    for(int x=0; x<n_a; x++)
    for(int y=0; y<n_a; y++)
        RF_P3_JK_A[((((t*n_a+u)*n_a+v)*n_a+y)*n_a+x)*n_a+w]+= RF_P3_VV_AAAAAA[((((t*n_a+v)*n_a+x)*n_a+y)*n_a+w)*n_a+u]
                                                             -RF_P3_VV_AAAAAA[((((u*n_a+v)*n_a+x)*n_a+y)*n_a+w)*n_a+t]
                                                             -RF_P3_VV_AAAAAA[((((t*n_a+v)*n_a+y)*n_a+x)*n_a+w)*n_a+u]
                                                             +RF_P3_VV_AAAAAA[((((u*n_a+v)*n_a+y)*n_a+x)*n_a+w)*n_a+t];
    for(int t=0; t<n_a; t++)
    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++)
    for(int w=0; w<n_a; w++)
    for(int x=0; x<n_a; x++)
    for(int y=0; y<n_a; y++)
        RF_P3_JK_B[((((t*n_a+u)*n_a+v)*n_a+y)*n_a+x)*n_a+w]+= RF_P3_VV_BBBBBB[((((t*n_a+v)*n_a+x)*n_a+y)*n_a+w)*n_a+u]
                                                             -RF_P3_VV_BBBBBB[((((u*n_a+v)*n_a+x)*n_a+y)*n_a+w)*n_a+t]
                                                             -RF_P3_VV_BBBBBB[((((t*n_a+v)*n_a+y)*n_a+x)*n_a+w)*n_a+u]
                                                             +RF_P3_VV_BBBBBB[((((u*n_a+v)*n_a+y)*n_a+x)*n_a+w)*n_a+t];
    //RF_PV_JK
    for(int t=0; t<n_a; t++)
    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++)
    for(int w=0; w<n_a; w++)
    for(int x=0; x<n_a; x++)
        RF_PV_JK_A[((u*n_a+t)*n_a+v)*n_a+w]+=-RF_P3_VV_AAAAAA[((((t*n_a+v)*n_a+x)*n_a+w)*n_a+x)*n_a+u]//A
                                             +RF_P3_VV_AAAAAA[((((u*n_a+v)*n_a+x)*n_a+w)*n_a+x)*n_a+t]
                                             +RF_P3_VV_AAAAAA[((((t*n_a+w)*n_a+x)*n_a+v)*n_a+x)*n_a+u]
                                             -RF_P3_VV_AAAAAA[((((u*n_a+w)*n_a+x)*n_a+v)*n_a+x)*n_a+t]
                                             -RF_P3_VV_AAAAAA[((((t*n_a+x)*n_a+v)*n_a+w)*n_a+u)*n_a+x]
                                             +RF_P3_VV_AAAAAA[((((t*n_a+v)*n_a+x)*n_a+w)*n_a+u)*n_a+x]
                                             +RF_P3_VV_AAAAAA[((((u*n_a+x)*n_a+v)*n_a+w)*n_a+t)*n_a+x]
                                             -RF_P3_VV_AAAAAA[((((u*n_a+v)*n_a+x)*n_a+w)*n_a+t)*n_a+x]
                                             +RF_P3_VV_AAAAAA[((((t*n_a+x)*n_a+w)*n_a+v)*n_a+u)*n_a+x]
                                             -RF_P3_VV_AAAAAA[((((t*n_a+w)*n_a+x)*n_a+v)*n_a+u)*n_a+x]
                                             -RF_P3_VV_AAAAAA[((((u*n_a+x)*n_a+w)*n_a+v)*n_a+t)*n_a+x]
                                             +RF_P3_VV_AAAAAA[((((u*n_a+w)*n_a+x)*n_a+v)*n_a+t)*n_a+x]
                                             +RF_P3_VV_AABAAB[((((t*n_a+v)*n_a+x)*n_a+w)*n_a+u)*n_a+x]//AB
                                             -RF_P3_VV_AABAAB[((((u*n_a+v)*n_a+x)*n_a+w)*n_a+t)*n_a+x]
                                             -RF_P3_VV_AABAAB[((((t*n_a+w)*n_a+x)*n_a+v)*n_a+u)*n_a+x]
                                             +RF_P3_VV_AABAAB[((((u*n_a+w)*n_a+x)*n_a+v)*n_a+t)*n_a+x];
    for(int t=0; t<n_a; t++)
    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++)
    for(int w=0; w<n_a; w++)
    for(int x=0; x<n_a; x++)
        RF_PV_JK_B[((u*n_a+t)*n_a+v)*n_a+w]+=-RF_P3_VV_BBBBBB[((((t*n_a+v)*n_a+x)*n_a+w)*n_a+x)*n_a+u]//B
                                             +RF_P3_VV_BBBBBB[((((u*n_a+v)*n_a+x)*n_a+w)*n_a+x)*n_a+t]
                                             +RF_P3_VV_BBBBBB[((((t*n_a+w)*n_a+x)*n_a+v)*n_a+x)*n_a+u]
                                             -RF_P3_VV_BBBBBB[((((u*n_a+w)*n_a+x)*n_a+v)*n_a+x)*n_a+t]
                                             -RF_P3_VV_BBBBBB[((((t*n_a+x)*n_a+v)*n_a+w)*n_a+u)*n_a+x]
                                             +RF_P3_VV_BBBBBB[((((t*n_a+v)*n_a+x)*n_a+w)*n_a+u)*n_a+x]
                                             +RF_P3_VV_BBBBBB[((((u*n_a+x)*n_a+v)*n_a+w)*n_a+t)*n_a+x]
                                             -RF_P3_VV_BBBBBB[((((u*n_a+v)*n_a+x)*n_a+w)*n_a+t)*n_a+x]
                                             +RF_P3_VV_BBBBBB[((((t*n_a+x)*n_a+w)*n_a+v)*n_a+u)*n_a+x]
                                             -RF_P3_VV_BBBBBB[((((t*n_a+w)*n_a+x)*n_a+v)*n_a+u)*n_a+x]
                                             -RF_P3_VV_BBBBBB[((((u*n_a+x)*n_a+w)*n_a+v)*n_a+t)*n_a+x]
                                             +RF_P3_VV_BBBBBB[((((u*n_a+w)*n_a+x)*n_a+v)*n_a+t)*n_a+x]
                                             +RF_P3_VV_BBABBA[((((t*n_a+v)*n_a+x)*n_a+w)*n_a+u)*n_a+x]//BA
                                             -RF_P3_VV_BBABBA[((((u*n_a+v)*n_a+x)*n_a+w)*n_a+t)*n_a+x]
                                             -RF_P3_VV_BBABBA[((((t*n_a+w)*n_a+x)*n_a+v)*n_a+u)*n_a+x]
                                             +RF_P3_VV_BBABBA[((((u*n_a+w)*n_a+x)*n_a+v)*n_a+t)*n_a+x];
    //RF_PV_JK
    for(int t=0; t<n_a; t++)
    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++)
    for(int w=0; w<n_a; w++)
        RF_PV_JK_A[((u*n_a+t)*n_a+v)*n_a+w]+=-RF_P3_VH_AAAA[((t*n_a+v)*n_a+w)*n_a+u]
                                             -RF_P3_HV_AAAA[((v*n_a+w)*n_a+u)*n_a+t]
                                             +RF_P3_VH_AAAA[((u*n_a+v)*n_a+w)*n_a+t]
                                             +RF_P3_HV_AAAA[((v*n_a+w)*n_a+t)*n_a+u]
                                             +RF_P3_VH_AAAA[((t*n_a+w)*n_a+v)*n_a+u]
                                             +RF_P3_HV_AAAA[((w*n_a+v)*n_a+u)*n_a+t]
                                             -RF_P3_VH_AAAA[((u*n_a+w)*n_a+v)*n_a+t]
                                             -RF_P3_HV_AAAA[((w*n_a+v)*n_a+t)*n_a+u];
    for(int t=0; t<n_a; t++)
    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++)
    for(int w=0; w<n_a; w++)
        RF_PV_JK_B[((u*n_a+t)*n_a+v)*n_a+w]+=-RF_P3_VH_AAAA[((t*n_a+v)*n_a+w)*n_a+u]
                                             -RF_P3_HV_AAAA[((v*n_a+w)*n_a+u)*n_a+t]
                                             +RF_P3_VH_AAAA[((u*n_a+v)*n_a+w)*n_a+t]
                                             +RF_P3_HV_AAAA[((v*n_a+w)*n_a+t)*n_a+u]
                                             +RF_P3_VH_AAAA[((t*n_a+w)*n_a+v)*n_a+u]
                                             +RF_P3_HV_AAAA[((w*n_a+v)*n_a+u)*n_a+t]
                                             -RF_P3_VH_AAAA[((u*n_a+w)*n_a+v)*n_a+t]
                                             -RF_P3_HV_AAAA[((w*n_a+v)*n_a+t)*n_a+u];

    
//  RF_P3_AB                 - order differs with the original XMCQDPT code (res_fit.cpp)
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++)
    for(int x=0; x<n_a; x++)
    for(int v=0; v<n_a; v++)
    for(int w=0; w<n_a; w++)
    for(int y=0; y<n_a; y++)
        //           A     A      B      A      A      B
        RF_P3_AB[((((t*n_a+u)*n_a+x)*n_a+v)*n_a+w)*n_a+y]+=-RF_P3_VV_AABAAB[((((t*n_a+v)*n_a+y)*n_a+w)*n_a+u)*n_a+x]//AABAAB
                                                           +RF_P3_VV_AABAAB[((((u*n_a+v)*n_a+y)*n_a+w)*n_a+t)*n_a+x]//AABAAB
                                                           +RF_P3_VV_AABAAB[((((t*n_a+w)*n_a+y)*n_a+v)*n_a+u)*n_a+x]//AABAAB
                                                           -RF_P3_VV_AABAAB[((((u*n_a+w)*n_a+y)*n_a+v)*n_a+t)*n_a+x]//AABAAB
                                                           -RF_P3_VV_AAABBA[((((t*n_a+v)*n_a+w)*n_a+y)*n_a+x)*n_a+u]//AAABBA
                                                           +RF_P3_VV_AAABBA[((((u*n_a+v)*n_a+w)*n_a+y)*n_a+x)*n_a+t]//AAABBA
                                                           +RF_P3_VV_AAABBA[((((t*n_a+w)*n_a+v)*n_a+y)*n_a+x)*n_a+u]//AAABBA
                                                           -RF_P3_VV_AAABBA[((((u*n_a+w)*n_a+v)*n_a+y)*n_a+x)*n_a+t]//AAABBA
                                                           -RF_P3_VV_BBAAAA[((((x*n_a+y)*n_a+v)*n_a+w)*n_a+u)*n_a+t]//BBAAAA
                                                           +RF_P3_VV_BBAAAA[((((x*n_a+y)*n_a+v)*n_a+w)*n_a+t)*n_a+u]//BBAAAA
                                                           +RF_P3_VV_BBAAAA[((((x*n_a+y)*n_a+w)*n_a+v)*n_a+u)*n_a+t]//BBAAAA
                                                           -RF_P3_VV_BBAAAA[((((x*n_a+y)*n_a+w)*n_a+v)*n_a+t)*n_a+u]//BBAAAA
                                                           ;
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++)
    for(int x=0; x<n_a; x++)
    for(int v=0; v<n_a; v++)
    for(int w=0; w<n_a; w++)
    for(int y=0; y<n_a; y++)
        //           B     B      A      B      B      A
        RF_P3_BA[((((t*n_a+u)*n_a+x)*n_a+v)*n_a+w)*n_a+y]+=-RF_P3_VV_BBABBA[((((t*n_a+v)*n_a+y)*n_a+w)*n_a+u)*n_a+x]//BBABBA
                                                           +RF_P3_VV_BBABBA[((((u*n_a+v)*n_a+y)*n_a+w)*n_a+t)*n_a+x]//BBABBA
                                                           +RF_P3_VV_BBABBA[((((t*n_a+w)*n_a+y)*n_a+v)*n_a+u)*n_a+x]//BBABBA
                                                           -RF_P3_VV_BBABBA[((((u*n_a+w)*n_a+y)*n_a+v)*n_a+t)*n_a+x]//BBABBA
                                                           -RF_P3_VV_BBBAAB[((((t*n_a+v)*n_a+w)*n_a+y)*n_a+x)*n_a+u]//BBBAAB
                                                           +RF_P3_VV_BBBAAB[((((u*n_a+v)*n_a+w)*n_a+y)*n_a+x)*n_a+t]//BBBAAB
                                                           +RF_P3_VV_BBBAAB[((((t*n_a+w)*n_a+v)*n_a+y)*n_a+x)*n_a+u]//BBBAAB
                                                           -RF_P3_VV_BBBAAB[((((u*n_a+w)*n_a+v)*n_a+y)*n_a+x)*n_a+t]//BBBAAB
                                                           -RF_P3_VV_AABBBB[((((x*n_a+y)*n_a+v)*n_a+w)*n_a+u)*n_a+t]//AABBBB
                                                           +RF_P3_VV_AABBBB[((((x*n_a+y)*n_a+v)*n_a+w)*n_a+t)*n_a+u]//AABBBB
                                                           +RF_P3_VV_AABBBB[((((x*n_a+y)*n_a+w)*n_a+v)*n_a+u)*n_a+t]//AABBBB
                                                           -RF_P3_VV_AABBBB[((((x*n_a+y)*n_a+w)*n_a+v)*n_a+t)*n_a+u]//AABBBB
                                                           ;

//     RF_PV_AB
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++)
    for(int v=0; v<n_a; v++)
    for(int w=0; w<n_a; w++)
    for(int x=0; x<n_a; x++)
        //         A     B       A     B
        RF_PV_AB[((t*n_a+u)*n_a+v)*n_a+w]+=-RF_P3_VV_AABAAB[((((t*n_a+x)*n_a+w)*n_a+v)*n_a+x)*n_a+u]// 1 - A?BA?B
                                           -RF_P3_VV_AAABBA[((((t*n_a+x)*n_a+v)*n_a+w)*n_a+u)*n_a+x]// 2 - A?ABB?
                                           +RF_P3_VV_AAABBA[((((t*n_a+v)*n_a+x)*n_a+w)*n_a+u)*n_a+x]// 3 - AA?BB? similar to line 10
                                           -RF_P3_VV_BBAAAA[((((u*n_a+w)*n_a+x)*n_a+v)*n_a+x)*n_a+t]// 4 - BB?A?A
                                           +RF_P3_VV_BBAAAA[((((u*n_a+w)*n_a+x)*n_a+v)*n_a+t)*n_a+x]// 5 - BB?AA? similar to line 8
                                           -RF_P3_VV_BBABBA[((((u*n_a+x)*n_a+v)*n_a+w)*n_a+x)*n_a+t]// 6 - B?AB&A
                                           -RF_P3_VV_BBBAAB[((((u*n_a+x)*n_a+w)*n_a+v)*n_a+t)*n_a+x]// 7 - B?BAA?
                                           +RF_P3_VV_BBBAAB[((((u*n_a+w)*n_a+x)*n_a+v)*n_a+t)*n_a+x]// 8 - BB?AA? similar to line 5
                                           -RF_P3_VV_AABBBB[((((t*n_a+v)*n_a+x)*n_a+w)*n_a+x)*n_a+u]// 9 - AA?B?B
                                           +RF_P3_VV_AABBBB[((((t*n_a+v)*n_a+x)*n_a+w)*n_a+u)*n_a+x]//10 - AA?BB? similar to line 3
                                           ;
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++)
    for(int v=0; v<n_a; v++)
    for(int w=0; w<n_a; w++)
    for(int x=0; x<n_a; x++)
        //         B     A       B     A
        RF_PV_BA[((t*n_a+u)*n_a+v)*n_a+w]+=-RF_P3_VV_BBABBA[((((t*n_a+x)*n_a+w)*n_a+v)*n_a+x)*n_a+u]// 1 - B?AB?A
                                           -RF_P3_VV_BBBAAB[((((t*n_a+x)*n_a+v)*n_a+w)*n_a+u)*n_a+x]// 2 - B?BAA?
                                           +RF_P3_VV_BBBAAB[((((t*n_a+v)*n_a+x)*n_a+w)*n_a+u)*n_a+x]// 3 - BB?AA? similar to line 10
                                           -RF_P3_VV_AABBBB[((((u*n_a+w)*n_a+x)*n_a+v)*n_a+x)*n_a+t]// 4 - AA?B?B
                                           +RF_P3_VV_AABBBB[((((u*n_a+w)*n_a+x)*n_a+v)*n_a+t)*n_a+x]// 5 - AA?BB? similar to line 8
                                           -RF_P3_VV_AABAAB[((((u*n_a+x)*n_a+v)*n_a+w)*n_a+x)*n_a+t]// 6 - A?BA&B
                                           -RF_P3_VV_AAABBA[((((u*n_a+x)*n_a+w)*n_a+v)*n_a+t)*n_a+x]// 7 - A?ABB?
                                           +RF_P3_VV_AAABBA[((((u*n_a+w)*n_a+x)*n_a+v)*n_a+t)*n_a+x]// 8 - AA?BB? similar to line 5
                                           -RF_P3_VV_BBAAAA[((((t*n_a+v)*n_a+x)*n_a+w)*n_a+x)*n_a+u]// 9 - BB?A?A
                                           +RF_P3_VV_BBAAAA[((((t*n_a+v)*n_a+x)*n_a+w)*n_a+u)*n_a+x]//10 - BB?AA? similar to line 3
                                           ;
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++)
    for(int v=0; v<n_a; v++)
    for(int w=0; w<n_a; w++)
        //         A     B       A     B
        RF_PV_AB[((t*n_a+u)*n_a+v)*n_a+w]+=-RF_P3_VH_BBAA[((u*n_a+w)*n_a+v)*n_a+t]//BBAA
                                           -RF_P3_HV_ABBA[((v*n_a+w)*n_a+u)*n_a+t]//ABBA
                                           -RF_P3_VH_AABB[((t*n_a+v)*n_a+w)*n_a+u]//AABB
                                           -RF_P3_HV_BAAB[((w*n_a+v)*n_a+t)*n_a+u]//BAAB
                                           ;
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++)
    for(int v=0; v<n_a; v++)
    for(int w=0; w<n_a; w++)
        //         B     A       B     A
        RF_PV_BA[((t*n_a+u)*n_a+v)*n_a+w]+=-RF_P3_VH_AABB[((u*n_a+w)*n_a+v)*n_a+t]//AABB
                                           -RF_P3_HV_BAAB[((v*n_a+w)*n_a+u)*n_a+t]//BAAB
                                           -RF_P3_VH_BBAA[((t*n_a+v)*n_a+w)*n_a+u]//BBAA
                                           -RF_P3_HV_ABBA[((w*n_a+v)*n_a+t)*n_a+u]//ABBA
                                           ;
//     set_zero_matr(RF_PH, n_a*n_a*N_fit);
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++)
    for(int v=0; v<n_a; v++)
    for(int y=0; y<n_a; y++)
        RF_PH_A[t*n_a+u]+=-RF_P3_VV_AAAAAA[((((t*n_a+y)*n_a+v)*n_a+u)*n_a+v)*n_a+y]//A??A??
                          +RF_P3_VV_AAAAAA[((((t*n_a+v)*n_a+y)*n_a+u)*n_a+v)*n_a+y]//A??A??
                          +RF_P3_VV_AABAAB[((((t*n_a+v)*n_a+y)*n_a+u)*n_a+v)*n_a+y]//A??A??
                          ;
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++)
    for(int v=0; v<n_a; v++)
    for(int y=0; y<n_a; y++)
        RF_PH_B[t*n_a+u]+=-RF_P3_VV_BBBBBB[((((t*n_a+y)*n_a+v)*n_a+u)*n_a+v)*n_a+y]//B??B??
                          +RF_P3_VV_BBBBBB[((((t*n_a+v)*n_a+y)*n_a+u)*n_a+v)*n_a+y]//B??B??
                          +RF_P3_VV_BBABBA[((((t*n_a+v)*n_a+y)*n_a+u)*n_a+v)*n_a+y]//B??B??
                          ;
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++)
    for(int y=0; y<n_a; y++)
        RF_PH_A[t*n_a+u]+=-RF_P3_VH_AAAA[((t*n_a+y)*n_a+u)*n_a+y]
                          +RF_P3_VH_AAAA[((t*n_a+u)*n_a+y)*n_a+y]
                          +RF_P3_VH_AABB[((t*n_a+u)*n_a+y)*n_a+y]
                          -RF_P3_HV_AAAA[((y*n_a+u)*n_a+y)*n_a+t]
                          +RF_P3_HV_AAAA[((y*n_a+u)*n_a+t)*n_a+y]
                          +RF_P3_HV_BAAB[((y*n_a+u)*n_a+t)*n_a+y];
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++)
    for(int y=0; y<n_a; y++)
        RF_PH_B[t*n_a+u]+=-RF_P3_VH_BBBB[((t*n_a+y)*n_a+u)*n_a+y]
                          +RF_P3_VH_BBBB[((t*n_a+u)*n_a+y)*n_a+y]
                          +RF_P3_VH_BBAA[((t*n_a+u)*n_a+y)*n_a+y]
                          -RF_P3_HV_BBBB[((y*n_a+u)*n_a+y)*n_a+t]
                          +RF_P3_HV_BBBB[((y*n_a+u)*n_a+t)*n_a+y]
                          +RF_P3_HV_ABBA[((y*n_a+u)*n_a+t)*n_a+y];
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++)
        RF_PH_A[t*n_a+u]+=-RF_P3_HH_A[u*n_a+t];
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++)
        RF_PH_B[t*n_a+u]+=-RF_P3_HH_A[u*n_a+t];
    
//     set_zero_matr(RF_PS, N_fit);
    for(int t=0; t<n_a; t++)
        RF_PS+=RF_P3_HH_A[t*n_a+t]+RF_P3_HH_B[t*n_a+t];
    

    delete[] CAAA_AA        ;
    delete[] CAAA_AB        ;
    delete[] CAAA_BA        ;
    delete[] CAAA_BB        ;
    delete[] RF_P3_VV_AAAAAA;
    delete[] RF_P3_VV_AAABBA;
    delete[] RF_P3_VV_AABAAB;
    delete[] RF_P3_VV_AABBBB;
    delete[] RF_P3_VV_BBAAAA;
    delete[] RF_P3_VV_BBABBA;
    delete[] RF_P3_VV_BBBAAB;
    delete[] RF_P3_VV_BBBBBB;
    delete[] RF_P3_VH_AAAA  ;
    delete[] RF_P3_VH_AABB  ;
    delete[] RF_P3_VH_BBAA  ;
    delete[] RF_P3_VH_BBBB  ;
    delete[] RF_P3_HV_AAAA  ;
    delete[] RF_P3_HV_ABBA  ;
    delete[] RF_P3_HV_BAAB  ;
    delete[] RF_P3_HV_BBBB  ;
    delete[] RF_P3_HH_A     ;
    delete[] RF_P3_HH_B     ;

    
    return 0;
}

int UPT_tensors::calc_EE_2_CCVV(){
    
    double * R;
    R = new double[num_threads];
    set_zero_matr(R,num_threads);
    double *V[num_threads];
    double *JK[num_threads];
    for(int i=0;i<num_threads;i++) V[i] = new double[n_c];
    for(int i=0;i<num_threads;i++)JK[i] = new double[n_c*n_c];
    //AA
    for(int a=0; a<n_v; a++){//fprintf(stderr,"CCVV a=%d\r",a);
#pragma omp parallel for
    for(int b=a; b<n_v; b++){/*fprintf(stderr,"CCVV a,b=%d,%d\r",a,b);*/
        
        int th_id = omp_get_thread_num();
        double C;
        
        int iajb=0;
        double dE;
        double Ep;
        double J;
        double K;
        
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        n_c,n_c,aux_n_ao,1.0,
                        RI_A->VC_RI_M+a*n_c*aux_n_ao,aux_n_ao,
                        RI_A->VC_RI_M+b*n_c*aux_n_ao,aux_n_ao,0.0,
                        JK[th_id],n_c);
        
        
        for(int i=0; i<n_c; i++){//fprintf(stderr,"i,a=%d,%d\r",i,a);
        for(int j=i; j<n_c; j++){
            J=JK[th_id][i*n_c+j];
            K=JK[th_id][j*n_c+i];
            C=J-K;
            V[th_id][j]=C*C;
        }
        
        Ep=e_c_A[i]-e_v_A[a]-e_v_A[b];
        for(int j=i; j<n_c; j++){
            dE=e_c_A[j]+Ep;
            R[th_id]+=V[th_id][j]*dE/(dE*dE+edshift);
            
        }
       
        }
    }
    }
    //BB
    for(int a=0; a<n_v; a++){//fprintf(stderr,"CCVV a=%d\r",a);
#pragma omp parallel for
    for(int b=a; b<n_v; b++){/*fprintf(stderr,"CCVV a,b=%d,%d\r",a,b);*/
        
        int th_id = omp_get_thread_num();
        double C;
        
        int iajb=0;
        double dE;
        double Ep;
        double J;
        double K;
        
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        n_c,n_c,aux_n_ao,1.0,
                        RI_B->VC_RI_M+a*n_c*aux_n_ao,aux_n_ao,
                        RI_B->VC_RI_M+b*n_c*aux_n_ao,aux_n_ao,0.0,
                        JK[th_id],n_c);
        
        
        for(int i=0; i<n_c; i++){//fprintf(stderr,"i,a=%d,%d\r",i,a);
        for(int j=i; j<n_c; j++){
            J=JK[th_id][i*n_c+j];
            K=JK[th_id][j*n_c+i];
            C=J-K;
            V[th_id][j]=C*C;
        }
        
        Ep=e_c_B[i]-e_v_B[a]-e_v_B[b];
        for(int j=i; j<n_c; j++){
            dE=e_c_B[j]+Ep;
            R[th_id]+=V[th_id][j]*dE/(dE*dE+edshift);
            
        }
       
        }
    }
    }
    //AB
    for(int a=0; a<n_v; a++){//fprintf(stderr,"CCVV a=%d\r",a);
#pragma omp parallel for
    for(int b=0; b<n_v; b++){/*fprintf(stderr,"CCVV a,b=%d,%d\r",a,b);*/
        
        int th_id = omp_get_thread_num();
        double C;
        
        int iajb=0;
        double dE;
        double Ep;
        double J;
        double K;
        
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        n_c,n_c,aux_n_ao,1.0,
                        RI_A->VC_RI_M+a*n_c*aux_n_ao,aux_n_ao,
                        RI_B->VC_RI_M+b*n_c*aux_n_ao,aux_n_ao,0.0,
                        JK[th_id],n_c);
        
        
        for(int i=0; i<n_c; i++){//fprintf(stderr,"i,a=%d,%d\r",i,a);
        for(int j=0; j<n_c; j++){
            J=JK[th_id][i*n_c+j];
            V[th_id][j]=J*J;
        }
        
        Ep=e_c_A[i]-e_v_A[a]-e_v_B[b];
        for(int j=0; j<n_c; j++){
            dE=e_c_B[j]+Ep;
            R[th_id]+=V[th_id][j]*dE/(dE*dE+edshift);
            
        }
       
        }
    }
    }
    for(int i=0;i<num_threads;i++)
        RF_PS+=R[i];
    
    delete[] R;
    for(int i=0;i<num_threads;i++) delete[] V [i];
    for(int i=0;i<num_threads;i++) delete[] JK[i];
    
    
    return 0;
}

int UPT_tensors::calc_EE_2_CAVV(){
    double **PH_th= new double * [num_threads];
    for(int i=0;i<num_threads;i++)PH_th[i]=new double[n_a*n_a];
    

    
    //AA+AB
    for(int i=0;i<num_threads;i++)set_zero_matr(PH_th[i],n_a*n_a);
    #pragma omp parallel
    {
        int nt = omp_get_thread_num();
        
        double *K;
        double *J;
        K= new double[n_a*n_c];
        J= new double[n_a*n_c];
        
        double dE;
        double V2;
        //AA
        for(int a=nt  ; a<n_v; a+=num_threads)
        for(int b=a  ; b<n_v; b++){//fprintf(stderr,"CAVVa,b=%d,%d\r",a,b);
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            n_c,n_a,aux_n_ao,1.0,
                            RI_A->VC_RI_M+a*n_c*aux_n_ao,aux_n_ao,
                            RI_A->VA_RI_M+b*n_a*aux_n_ao,aux_n_ao,0.0,
                            J,n_a);
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            n_c,n_a,aux_n_ao,1.0,
                            RI_A->VC_RI_M+b*n_c*aux_n_ao,aux_n_ao,
                            RI_A->VA_RI_M+a*n_a*aux_n_ao,aux_n_ao,0.0,
                            K,n_a);
            for(int i=0;i<n_c;i++)
            for(int t=0;t<n_a;t++)
            for(int v=0;v<n_a;v++){
                V2 =(J[i*n_a+t]-K[i*n_a+t])*(J[i*n_a+v]-K[i*n_a+v]);
//                 V2+=J[i*n_a+t]*J[i*n_a+v];
//                 if(a!=b)
//                     V2+=K[i*n_a+t]*K[i*n_a+v];
                dE=e_c_A[i]+e_IP_A[t]-e_v_A[a]-e_v_A[b];
                PH_th[nt][t*n_a+v]+=V2*dE/(dE*dE+edshift);
                
            }
        }
        //AB
        for(int a=nt  ; a<n_v; a+=num_threads)
        for(int b=0   ; b<n_v; b++){//fprintf(stderr,"CAVVa,b=%d,%d\r",a,b);
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            n_c,n_a,aux_n_ao,1.0,
                            RI_B->VC_RI_M+a*n_c*aux_n_ao,aux_n_ao,
                            RI_A->VA_RI_M+b*n_a*aux_n_ao,aux_n_ao,0.0,
                            J,n_a);
            for(int i=0;i<n_c;i++)
            for(int t=0;t<n_a;t++)
            for(int v=0;v<n_a;v++){
                V2 =J[i*n_a+t]*J[i*n_a+v];
                dE=e_c_B[i]+e_IP_A[t]-e_v_B[a]-e_v_A[b];
                PH_th[nt][t*n_a+v]+=V2*dE/(dE*dE+edshift);
                
            }
        }
        
        delete[] J;
        delete[] K;
    }
    
    for(long j=0; j<num_threads;j++)
    #pragma omp parallel for
        for(long i=0; i<n_a*n_a;i++)
            RF_PH_A[i]+=PH_th[j][i];
    //BB+BA
    for(int i=0;i<num_threads;i++)set_zero_matr(PH_th[i],n_a*n_a);
    #pragma omp parallel
    {
        int nt = omp_get_thread_num();
        
        double *K;
        double *J;
        K= new double[n_a*n_c];
        J= new double[n_a*n_c];
        
        double dE;
        double V2;
        //BB
        for(int a=nt  ; a<n_v; a+=num_threads)
        for(int b=a  ; b<n_v; b++){//fprintf(stderr,"CAVVa,b=%d,%d\r",a,b);
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            n_c,n_a,aux_n_ao,1.0,
                            RI_B->VC_RI_M+a*n_c*aux_n_ao,aux_n_ao,
                            RI_B->VA_RI_M+b*n_a*aux_n_ao,aux_n_ao,0.0,
                            J,n_a);
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            n_c,n_a,aux_n_ao,1.0,
                            RI_B->VC_RI_M+b*n_c*aux_n_ao,aux_n_ao,
                            RI_B->VA_RI_M+a*n_a*aux_n_ao,aux_n_ao,0.0,
                            K,n_a);
            for(int i=0;i<n_c;i++)
            for(int t=0;t<n_a;t++)
            for(int v=0;v<n_a;v++){
                V2 =(J[i*n_a+t]-K[i*n_a+t])*(J[i*n_a+v]-K[i*n_a+v]);
//                 V2+=J[i*n_a+t]*J[i*n_a+v];
//                 if(a!=b)
//                     V2+=K[i*n_a+t]*K[i*n_a+v];
                dE=e_c_B[i]+e_IP_B[t]-e_v_B[a]-e_v_B[b];
                PH_th[nt][t*n_a+v]+=V2*dE/(dE*dE+edshift);
                
            }
        }
        //BA
        for(int a=nt  ; a<n_v; a+=num_threads)
        for(int b=0   ; b<n_v; b++){//fprintf(stderr,"CAVVa,b=%d,%d\r",a,b);
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            n_c,n_a,aux_n_ao,1.0,
                            RI_A->VC_RI_M+a*n_c*aux_n_ao,aux_n_ao,
                            RI_B->VA_RI_M+b*n_a*aux_n_ao,aux_n_ao,0.0,
                            J,n_a);
            for(int i=0;i<n_c;i++)
            for(int t=0;t<n_a;t++)
            for(int v=0;v<n_a;v++){
                V2 =J[i*n_a+t]*J[i*n_a+v];
                dE=e_c_A[i]+e_IP_B[t]-e_v_A[a]-e_v_B[b];
                PH_th[nt][t*n_a+v]+=V2*dE/(dE*dE+edshift);
                
            }
        }
        
        delete[] J;
        delete[] K;
    }
    
    for(long j=0; j<num_threads;j++)
    #pragma omp parallel for
        for(long i=0; i<n_a*n_a;i++)
            RF_PH_B[i]+=PH_th[j][i];
        
    for(int i=0;i<num_threads;i++)delete[] PH_th[i];
    delete[] PH_th;
    
    
    return 0;
    
}

int UPT_tensors::calc_EE_2_AAVV(){
    
//     num_threads=1;
    double **JK_th= new double * [num_threads];
    for(int i=0;i<num_threads;i++)JK_th[i]=new double[n_a*n_a*n_a*n_a];

    double **AB_th= new double * [num_threads];
    for(int i=0;i<num_threads;i++)AB_th[i]=new double[n_a*n_a*n_a*n_a];

    //AA+AB
    for(int i=0;i<num_threads;i++)set_zero_matr(JK_th[i],n_a*n_a*n_a*n_a);
    for(int i=0;i<num_threads;i++)set_zero_matr(AB_th[i],n_a*n_a*n_a*n_a);
    #pragma omp parallel
    {
        int nt =omp_get_thread_num();
        
        double AB,JK,dE;
        
        double * Ja;
        
        Ja = new double[n_a*n_a];
        //AA
        for(int a=nt; a<n_v; a+=num_threads)
        for(int b=a; b<n_v; b++){//fprintf(stderr,"AAVV a,b=%3d,%3d  (thread %2d)    \r",a,b,nt);
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            n_a,n_a,aux_n_ao,1.0,
                            RI_A->VA_RI_M+a*n_a*aux_n_ao,aux_n_ao,
                            RI_A->VA_RI_M+b*n_a*aux_n_ao,aux_n_ao,0.0,
                            Ja,n_a);
            for(int t=0  ;t<n_a;t++)
            for(int u=0  ;u<n_a;u++)
            for(int v=0  ;v<n_a;v++)
            for(int w=0  ;w<n_a;w++){
                JK=(Ja[t*n_a+u]-Ja[u*n_a+t])*
                   (Ja[v*n_a+w]-Ja[w*n_a+v]);
                   
                dE=e_IP_A[t]+e_IP_A[u]-e_v_A[a]-e_v_A[b];
                JK_th[nt][((u*n_a+t)*n_a+v)*n_a+w]+=JK*dE/(dE*dE+edshift);///inverse t-u in res_fit_calc_2body_AA (?)
                
            }
                
            
        }
        //AB
        for(int a=nt; a<n_v; a+=num_threads)
        for(int b= 0; b<n_v; b++){//fprintf(stderr,"AAVV a,b=%3d,%3d  (thread %2d)    \r",a,b,nt);
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            n_a,n_a,aux_n_ao,1.0,
                            RI_A->VA_RI_M+a*n_a*aux_n_ao,aux_n_ao,
                            RI_B->VA_RI_M+b*n_a*aux_n_ao,aux_n_ao,0.0,
                            Ja,n_a);
            for(int t=0  ;t<n_a;t++)
            for(int u=0  ;u<n_a;u++)
            for(int v=0  ;v<n_a;v++)
            for(int w=0  ;w<n_a;w++){
                   
                AB         =Ja[t*n_a+u]*Ja[v*n_a+w];
                dE=e_IP_A[t]+e_IP_B[u]-e_v_A[a]-e_v_B[b];
                AB_th[nt][((t*n_a+u)*n_a+v)*n_a+w]+=AB*dE/(dE*dE+edshift);
                
            }
                
        }
        delete[] Ja;
    
    
    }
    
    for(long j=0; j<num_threads;j++)
    for(long i=0; i<n_a*n_a*n_a*n_a;i++)
        RF_PV_JK_A[i]+=JK_th[j][i];
    
    for(long j=0; j<num_threads;j++)
    for(long i=0; i<n_a*n_a*n_a*n_a;i++)
        RF_PV_AB[i]+=AB_th[j][i];
        
    //BB+BA
    for(int i=0;i<num_threads;i++)set_zero_matr(JK_th[i],n_a*n_a*n_a*n_a);
    for(int i=0;i<num_threads;i++)set_zero_matr(AB_th[i],n_a*n_a*n_a*n_a);
    #pragma omp parallel
    {
        int nt =omp_get_thread_num();
        
        double AB,JK,dE;
        
        double * Ja;
        
        Ja = new double[n_a*n_a];
        //BB
        for(int a=nt; a<n_v; a+=num_threads)
        for(int b=a; b<n_v; b++){//fprintf(stderr,"AAVV a,b=%3d,%3d  (thread %2d)    \r",a,b,nt);
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            n_a,n_a,aux_n_ao,1.0,
                            RI_B->VA_RI_M+a*n_a*aux_n_ao,aux_n_ao,
                            RI_B->VA_RI_M+b*n_a*aux_n_ao,aux_n_ao,0.0,
                            Ja,n_a);
            for(int t=0  ;t<n_a;t++)
            for(int u=0  ;u<n_a;u++)
            for(int v=0  ;v<n_a;v++)
            for(int w=0  ;w<n_a;w++){
                JK=(Ja[t*n_a+u]-Ja[u*n_a+t])*
                   (Ja[v*n_a+w]-Ja[w*n_a+v]);
                   
                dE=e_IP_B[t]+e_IP_B[u]-e_v_B[a]-e_v_B[b];
                JK_th[nt][((u*n_a+t)*n_a+v)*n_a+w]+=JK*dE/(dE*dE+edshift);///inverse t-u in res_fit_calc_2body_AA (?)
                
            }
                
            
        }
        //BA
        for(int a=nt; a<n_v; a+=num_threads)
        for(int b= 0; b<n_v; b++){//fprintf(stderr,"AAVV a,b=%3d,%3d  (thread %2d)    \r",a,b,nt);
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            n_a,n_a,aux_n_ao,1.0,
                            RI_B->VA_RI_M+a*n_a*aux_n_ao,aux_n_ao,
                            RI_A->VA_RI_M+b*n_a*aux_n_ao,aux_n_ao,0.0,
                            Ja,n_a);
            for(int t=0  ;t<n_a;t++)
            for(int u=0  ;u<n_a;u++)
            for(int v=0  ;v<n_a;v++)
            for(int w=0  ;w<n_a;w++){
                   
                AB         =Ja[t*n_a+u]*Ja[v*n_a+w];
                dE=e_IP_B[t]+e_IP_A[u]-e_v_B[a]-e_v_A[b];
                AB_th[nt][((t*n_a+u)*n_a+v)*n_a+w]+=AB*dE/(dE*dE+edshift);
                
            }
                
        }
        delete[] Ja;
    
    
    }
    
    for(long j=0; j<num_threads;j++)
    for(long i=0; i<n_a*n_a*n_a*n_a;i++)
        RF_PV_JK_B[i]+=JK_th[j][i];
    for(long j=0; j<num_threads;j++)
    for(long i=0; i<n_a*n_a*n_a*n_a;i++)
        RF_PV_BA[i]+=AB_th[j][i];
        
    for(int i=0;i<num_threads;i++)delete[] JK_th[i];
    delete[] JK_th;
    
    for(int i=0;i<num_threads;i++)delete[] AB_th[i];
    delete[] AB_th;
    return 0;
}

int UPT_tensors::calc_EE_2_CCAV(){
    
    
    
    
    double **RF_PH_th= new double * [num_threads];
    for(int i=0;i<num_threads;i++)RF_PH_th[i]   = new double[n_a*n_a];
    
    //AA+AB
    for(int i=0;i<num_threads;i++)set_zero_matr(RF_PH_th[i], n_a*n_a);
    
    #pragma omp parallel
    {
        int nt =omp_get_thread_num();
        
        double *K;
        K= new double[n_a*n_v];
        double *J;
        J= new double[n_a*n_v];
        
        double dE;
        double V2;
        //AA
        for(int i=nt; i<n_c; i+=num_threads)
        for(int j=i; j<n_c; j++){//fprintf(stderr," CCAV i,j=%d,%d\r",i,j);
            for(int a=0; a<n_v; a++)
            for(int t=0; t<n_a; t++){
                J[a*n_a+t]=cblas_ddot(aux_n_ao, RI_A->VC_RI_M+(a*n_c+i)*aux_n_ao, 1, RI_A->CA_RI_M+(j*n_a+t)*aux_n_ao, 1);
                K[a*n_a+t]=cblas_ddot(aux_n_ao, RI_A->VC_RI_M+(a*n_c+j)*aux_n_ao, 1, RI_A->CA_RI_M+(i*n_a+t)*aux_n_ao, 1);
            }
            
            for(int a=0;a<n_v;a++)
            for(int t=0;t<n_a;t++)
            for(int v=0;v<n_a;v++){
                V2 =(J[a*n_a+t]-K[a*n_a+t])*(J[a*n_a+v]-K[a*n_a+v]);
                dE=e_c_A[i]+e_c_A[j]-e_EA_A[t]-e_v_A[a];
                RF_PH_th[nt][v*n_a+t]-=V2*dE/(dE*dE+edshift);//sign changed!!!!
                
            }
        }
        //AB
        for(int i=nt; i<n_c; i+=num_threads)
        for(int j= 0; j<n_c; j++){//fprintf(stderr," CCAV i,j=%d,%d\r",i,j);
            for(int a=0; a<n_v; a++)
            for(int t=0; t<n_a; t++){
                J[a*n_a+t]=cblas_ddot(aux_n_ao, RI_B->VC_RI_M+(a*n_c+i)*aux_n_ao, 1, RI_A->CA_RI_M+(j*n_a+t)*aux_n_ao, 1);
            }
            
            for(int a=0;a<n_v;a++)
            for(int t=0;t<n_a;t++)
            for(int v=0;v<n_a;v++){
                V2 =J[a*n_a+t]*J[a*n_a+v];
                dE=e_c_B[i]+e_c_A[j]-e_EA_A[t]-e_v_B[a];
                RF_PH_th[nt][v*n_a+t]-=V2*dE/(dE*dE+edshift);//sign changed!!!!
                
            }
        }
        
        delete[] K;
        delete[] J;
    }
    
    for(long j=1; j<num_threads;j++)
    for(long i=0; i<n_a*n_a;i++)
        RF_PH_th[0][i]+=RF_PH_th[j][i];
    
    for(int t=0; t<n_a; t++)
        RF_PS-=RF_PH_th[0][t*n_a+t];//A part
    
    #pragma omp parallel for // parallelism is bad but i don't care
    for(int w=0; w<n_a; w++)
    for(int t=0; t<n_a; t++)
        RF_PH_A[w*n_a+t]+=RF_PH_th[0][w*n_a+t];
    
    
    
    //BB+BA
    for(int i=0;i<num_threads;i++)set_zero_matr(RF_PH_th[i], n_a*n_a);
    
    #pragma omp parallel
    {
        int nt =omp_get_thread_num();
        
        double *K;
        K= new double[n_a*n_v];
        double *J;
        J= new double[n_a*n_v];
        
        double dE;
        double V2;
        //BB
        for(int i=nt; i<n_c; i+=num_threads)
        for(int j=i; j<n_c; j++){//fprintf(stderr," CCAV i,j=%d,%d\r",i,j);
            for(int a=0; a<n_v; a++)
            for(int t=0; t<n_a; t++){
                J[a*n_a+t]=cblas_ddot(aux_n_ao, RI_B->VC_RI_M+(a*n_c+i)*aux_n_ao, 1, RI_B->CA_RI_M+(j*n_a+t)*aux_n_ao, 1);
                K[a*n_a+t]=cblas_ddot(aux_n_ao, RI_B->VC_RI_M+(a*n_c+j)*aux_n_ao, 1, RI_B->CA_RI_M+(i*n_a+t)*aux_n_ao, 1);
            }
            
            for(int a=0;a<n_v;a++)
            for(int t=0;t<n_a;t++)
            for(int v=0;v<n_a;v++){
                V2 =(J[a*n_a+t]-K[a*n_a+t])*(J[a*n_a+v]-K[a*n_a+v]);
                dE=e_c_B[i]+e_c_B[j]-e_EA_B[t]-e_v_B[a];
                RF_PH_th[nt][v*n_a+t]-=V2*dE/(dE*dE+edshift);//sign changed!!!!
                
            }
        }
        //BA
        for(int i=nt; i<n_c; i+=num_threads)
        for(int j= 0; j<n_c; j++){//fprintf(stderr," CCAV i,j=%d,%d\r",i,j);
            for(int a=0; a<n_v; a++)
            for(int t=0; t<n_a; t++){
                J[a*n_a+t]=cblas_ddot(aux_n_ao, RI_A->VC_RI_M+(a*n_c+i)*aux_n_ao, 1, RI_B->CA_RI_M+(j*n_a+t)*aux_n_ao, 1);
            }
            
            for(int a=0;a<n_v;a++)
            for(int t=0;t<n_a;t++)
            for(int v=0;v<n_a;v++){
                V2 =J[a*n_a+t]*J[a*n_a+v];
                dE=e_c_A[i]+e_c_B[j]-e_EA_B[t]-e_v_A[a];
                RF_PH_th[nt][v*n_a+t]-=V2*dE/(dE*dE+edshift);//sign changed!!!!
                
            }
        }
        
        delete[] K;
        delete[] J;
    }
    
    for(long j=1; j<num_threads;j++)
    for(long i=0; i<n_a*n_a;i++)
        RF_PH_th[0][i]+=RF_PH_th[j][i];
    
    for(int t=0; t<n_a; t++)
        RF_PS-=RF_PH_th[0][t*n_a+t];//B part
    
    #pragma omp parallel for // parallelism is bad but i don't care
    for(int w=0; w<n_a; w++)
    for(int t=0; t<n_a; t++)
        RF_PH_B[w*n_a+t]+=RF_PH_th[0][w*n_a+t];

    
    for(int i=0;i<num_threads;i++)delete[] RF_PH_th[i];
    delete[] RF_PH_th;
    
    
    return 0;
    
}

int UPT_tensors::calc_EE_2_CCAA(){
    
//     num_threads=1;
    
    double ** RF_PV_AB_th = new double * [num_threads];
    double ** RF_PV_JK_th = new double * [num_threads];
    
    for(int i=0;i<num_threads;i++)RF_PV_AB_th[i] = new double[n_a*n_a*n_a*n_a];
    for(int i=0;i<num_threads;i++)RF_PV_JK_th[i] = new double[n_a*n_a*n_a*n_a];
    
    
    //AA+AB
    for(int i=0;i<num_threads;i++)set_zero_matr(RF_PV_AB_th[i], n_a*n_a*n_a*n_a);
    for(int i=0;i<num_threads;i++)set_zero_matr(RF_PV_JK_th[i], n_a*n_a*n_a*n_a);
    
    
    #pragma omp parallel
    {
        int nt =omp_get_thread_num();
    
        double dE;
        double Ed;
        
        double V2;
        double AB;
        
        double * Ja;
        Ja = new double[n_a*n_a];
        //AA
        for(int i=nt; i<n_c; i+=num_threads)
        for(int j=i; j<n_c; j++){//fprintf(stderr,"CCAA i,j=%3d,%3d  (thread %2d)    \r",i,j,nt);
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            n_a,n_a,aux_n_ao,1.0,
                            RI_A->CA_RI_M+i*n_a*aux_n_ao,aux_n_ao,
                            RI_A->CA_RI_M+j*n_a*aux_n_ao,aux_n_ao,0.0,
                            Ja,n_a);
            for(int t=0  ;t<n_a;t++)
            for(int u=0  ;u<n_a;u++)
            for(int v=0  ;v<n_a;v++)
            for(int w=0  ;w<n_a;w++){
                
                V2=(Ja[t*n_a+u]-Ja[u*n_a+t])*
                   (Ja[v*n_a+w]-Ja[w*n_a+v]);
                
                dE=-e_EA_A[v]-e_EA_A[w]+e_c_A[i]+e_c_A[j];
                
                RF_PV_JK_th[nt][((u*n_a+t)*n_a+v)*n_a+w]+=V2*dE/(dE*dE+edshift);///inverse t-u in res_fit_calc_2body_AA (?)
                
            }
        }
        //AB
        for(int i=nt; i<n_c; i+=num_threads)
        for(int j= 0; j<n_c; j++){//fprintf(stderr,"CCAA i,j=%3d,%3d  (thread %2d)    \r",i,j,nt);
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            n_a,n_a,aux_n_ao,1.0,
                            RI_A->CA_RI_M+i*n_a*aux_n_ao,aux_n_ao,
                            RI_B->CA_RI_M+j*n_a*aux_n_ao,aux_n_ao,0.0,
                            Ja,n_a);
            for(int t=0  ;t<n_a;t++)
            for(int u=0  ;u<n_a;u++)
            for(int v=0  ;v<n_a;v++)
            for(int w=0  ;w<n_a;w++){
                
                AB         =Ja[t*n_a+u]*Ja[v*n_a+w];
                
                dE=-e_EA_A[v]-e_EA_B[w]+e_c_A[i]+e_c_B[j];
                
                RF_PV_AB_th[nt][((t*n_a+u)*n_a+v)*n_a+w]+=AB*dE/(dE*dE+edshift);
                
            }
        }
        
        delete[] Ja;
    }
    
    for(long j=1; j<num_threads;j++)
    for(long i=0; i<n_a*n_a*n_a*n_a;i++)
        RF_PV_JK_th[0][i]+=RF_PV_JK_th[j][i];

    for(long j=1; j<num_threads;j++)
    for(long i=0; i<n_a*n_a*n_a*n_a;i++)
        RF_PV_AB_th[0][i]+=RF_PV_AB_th[j][i];

    
    #pragma omp parallel for // parallelism is bad but i don't care
    for(long i=0;i<n_a*n_a*n_a*n_a;i++)RF_PV_JK_A[i]+=RF_PV_JK_th[0][i];
    
    #pragma omp parallel for // parallelism is bad but i don't care
    for(long i=0;i<n_a*n_a*n_a*n_a;i++)RF_PV_AB[i]+=RF_PV_AB_th[0][i];
    
    #pragma omp parallel for // parallelism is bad but i don't care
    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++)
    for(int t=0; t<n_a; t++)
        RF_PH_A[v*n_a+u]+=RF_PV_JK_th[0][((t*n_a+v)*n_a+t)*n_a+u]*0.25;
    
    #pragma omp parallel for // parallelism is bad but i don't care
    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++)
    for(int t=0; t<n_a; t++)
        RF_PH_A[v*n_a+t]-=RF_PV_JK_th[0][((u*n_a+v)*n_a+t)*n_a+u]*0.25;

    #pragma omp parallel for // parallelism is bad but i don't care
    for(int w=0; w<n_a; w++)
    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++)
        RF_PH_A[w*n_a+u]-=RF_PV_JK_th[0][((w*n_a+v)*n_a+v)*n_a+u]*0.25;

    #pragma omp parallel for // parallelism is bad but i don't care
    for(int w=0; w<n_a; w++)
    for(int u=0; u<n_a; u++)
    for(int t=0; t<n_a; t++)
        RF_PH_A[w*n_a+t]+=RF_PV_JK_th[0][((w*n_a+u)*n_a+t)*n_a+u]*0.25;

    #pragma omp parallel for // parallelism is bad but i don't care
    for(int w=0; w<n_a; w++)
    for(int u=0; u<n_a; u++)
    for(int t=0; t<n_a; t++)
        RF_PH_B[w*n_a+u]-=RF_PV_AB_th[0][((t*n_a+w)*n_a+t)*n_a+u]*0.5;

    #pragma omp parallel for // parallelism is bad but i don't care
    for(int w=0; w<n_a; w++)
    for(int u=0; u<n_a; u++)
    for(int t=0; t<n_a; t++)
        RF_PH_A[w*n_a+u]-=RF_PV_AB_th[0][((w*n_a+t)*n_a+u)*n_a+t]*0.5;



    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++)
        RF_PS+=RF_PV_JK_th[0][((u*n_a+v)*n_a+v)*n_a+u]*0.25;//A+B

    for(int v=0; v<n_a; v++)
    for(int w=0; w<n_a; w++)
        RF_PS-=RF_PV_JK_th[0][((w*n_a+v)*n_a+w)*n_a+v]*0.25;//A+B
    
    for(int u=0; u<n_a; u++)
    for(int t=0; t<n_a; t++)
        RF_PS+=RF_PV_AB_th[0][((t*n_a+u)*n_a+t)*n_a+u]*0.5;

    
    //AA+AB
    for(int i=0;i<num_threads;i++)set_zero_matr(RF_PV_AB_th[i], n_a*n_a*n_a*n_a);
    for(int i=0;i<num_threads;i++)set_zero_matr(RF_PV_JK_th[i], n_a*n_a*n_a*n_a);
    
    
    #pragma omp parallel
    {
        int nt =omp_get_thread_num();
    
        double dE;
        double Ed;
        
        double V2;
        double AB;
        
        double * Ja;
        Ja = new double[n_a*n_a];
        //BB
        for(int i=nt; i<n_c; i+=num_threads)
        for(int j=i; j<n_c; j++){//fprintf(stderr,"CCAA i,j=%3d,%3d  (thread %2d)    \r",i,j,nt);
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            n_a,n_a,aux_n_ao,1.0,
                            RI_B->CA_RI_M+i*n_a*aux_n_ao,aux_n_ao,
                            RI_B->CA_RI_M+j*n_a*aux_n_ao,aux_n_ao,0.0,
                            Ja,n_a);
            for(int t=0  ;t<n_a;t++)
            for(int u=0  ;u<n_a;u++)
            for(int v=0  ;v<n_a;v++)
            for(int w=0  ;w<n_a;w++){
                
                V2=(Ja[t*n_a+u]-Ja[u*n_a+t])*
                   (Ja[v*n_a+w]-Ja[w*n_a+v]);
                
                dE=-e_EA_B[v]-e_EA_B[w]+e_c_B[i]+e_c_B[j];
                
                RF_PV_JK_th[nt][((u*n_a+t)*n_a+v)*n_a+w]+=V2*dE/(dE*dE+edshift);///inverse t-u in res_fit_calc_2body_AA (?)
                
            }
        }
        //BA
        for(int i=nt; i<n_c; i+=num_threads)
        for(int j= 0; j<n_c; j++){//fprintf(stderr,"CCAA i,j=%3d,%3d  (thread %2d)    \r",i,j,nt);
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            n_a,n_a,aux_n_ao,1.0,
                            RI_B->CA_RI_M+i*n_a*aux_n_ao,aux_n_ao,
                            RI_A->CA_RI_M+j*n_a*aux_n_ao,aux_n_ao,0.0,
                            Ja,n_a);
            for(int t=0  ;t<n_a;t++)
            for(int u=0  ;u<n_a;u++)
            for(int v=0  ;v<n_a;v++)
            for(int w=0  ;w<n_a;w++){
                
                AB         =Ja[t*n_a+u]*Ja[v*n_a+w];
                
                dE=-e_EA_B[v]-e_EA_A[w]+e_c_B[i]+e_c_A[j];
                
                RF_PV_AB_th[nt][((t*n_a+u)*n_a+v)*n_a+w]+=AB*dE/(dE*dE+edshift);
                
            }
        }
        
        delete[] Ja;
    }
    
    for(long j=1; j<num_threads;j++)
    for(long i=0; i<n_a*n_a*n_a*n_a;i++)
        RF_PV_JK_th[0][i]+=RF_PV_JK_th[j][i];

    for(long j=1; j<num_threads;j++)
    for(long i=0; i<n_a*n_a*n_a*n_a;i++)
        RF_PV_AB_th[0][i]+=RF_PV_AB_th[j][i];

    
    #pragma omp parallel for // parallelism is bad but i don't care
    for(long i=0;i<n_a*n_a*n_a*n_a;i++)RF_PV_JK_B[i]+=RF_PV_JK_th[0][i];
    
    #pragma omp parallel for // parallelism is bad but i don't care
    for(long i=0;i<n_a*n_a*n_a*n_a;i++)RF_PV_BA[i]+=RF_PV_AB_th[0][i];
    
    #pragma omp parallel for // parallelism is bad but i don't care
    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++)
    for(int t=0; t<n_a; t++)
        RF_PH_B[v*n_a+u]+=RF_PV_JK_th[0][((t*n_a+v)*n_a+t)*n_a+u]*0.25;
    
    #pragma omp parallel for // parallelism is bad but i don't care
    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++)
    for(int t=0; t<n_a; t++)
        RF_PH_B[v*n_a+t]-=RF_PV_JK_th[0][((u*n_a+v)*n_a+t)*n_a+u]*0.25;

    #pragma omp parallel for // parallelism is bad but i don't care
    for(int w=0; w<n_a; w++)
    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++)
        RF_PH_B[w*n_a+u]-=RF_PV_JK_th[0][((w*n_a+v)*n_a+v)*n_a+u]*0.25;

    #pragma omp parallel for // parallelism is bad but i don't care
    for(int w=0; w<n_a; w++)
    for(int u=0; u<n_a; u++)
    for(int t=0; t<n_a; t++)
        RF_PH_B[w*n_a+t]+=RF_PV_JK_th[0][((w*n_a+u)*n_a+t)*n_a+u]*0.25;

    #pragma omp parallel for // parallelism is bad but i don't care
    for(int w=0; w<n_a; w++)
    for(int u=0; u<n_a; u++)
    for(int t=0; t<n_a; t++)
        RF_PH_A[w*n_a+u]-=RF_PV_AB_th[0][((t*n_a+w)*n_a+t)*n_a+u]*0.5;
    
    #pragma omp parallel for // parallelism is bad but i don't care
    for(int w=0; w<n_a; w++)
    for(int u=0; u<n_a; u++)
    for(int t=0; t<n_a; t++)
        RF_PH_B[w*n_a+u]-=RF_PV_AB_th[0][((w*n_a+t)*n_a+u)*n_a+t]*0.5;


    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++)
        RF_PS+=RF_PV_JK_th[0][((u*n_a+v)*n_a+v)*n_a+u]*0.25;//A+B

    for(int v=0; v<n_a; v++)
    for(int w=0; w<n_a; w++)
        RF_PS-=RF_PV_JK_th[0][((w*n_a+v)*n_a+w)*n_a+v]*0.25;//A+B
    
    for(int u=0; u<n_a; u++)
    for(int t=0; t<n_a; t++)
        RF_PS+=RF_PV_AB_th[0][((t*n_a+u)*n_a+t)*n_a+u]*0.5;

    
    
    for(int i=0;i<num_threads;i++)delete[] RF_PV_JK_th[i];
    for(int i=0;i<num_threads;i++)delete[] RF_PV_AB_th[i];
    
    
    delete[] RF_PV_JK_th;
    delete[] RF_PV_AB_th;
    
    return 0;
}


