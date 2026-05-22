# include <stdio.h>
# include "blas_link.h"
# include <omp.h>


# include "CI.h"
# include "matr.h"
# include "timer.h"

# include "from_hash.h"
# include "aldet.h"
# include "RI.h"
# include "PT_tensors_rel.h"

extern int num_threads;


inline double ED_with_shift(double E, double edshift){
    
    return E/(E*E+edshift);
}

PT_tensors_rel::PT_tensors_rel(){
            
    RF_PH_r = nullptr;
    RF_PH_i = nullptr;
    RF_PV   = nullptr;
    RF_PV_r = nullptr;
    RF_PV_i = nullptr;
    RF_P3_r = nullptr;
    RF_P3_i = nullptr;
    
    e_IP     = nullptr;
    e_EA     = nullptr;
    e_IP_2   = nullptr;
    e_EA_2   = nullptr;
    e_IP_AB  = nullptr;
    e_EA_AB  = nullptr;
    
    IP_U     = nullptr;
    EA_U     = nullptr;
    IP_Um    = nullptr;
    EA_Um    = nullptr;
    IP_Up    = nullptr;
    EA_Up    = nullptr;
     
    IP_2_U   = nullptr;
    EA_2_U   = nullptr;
    IP_2_Um  = nullptr;
    EA_2_Um  = nullptr;
    IP_2_Up  = nullptr;
    EA_2_Up  = nullptr;
    
    IP_AB_U  = nullptr;
    EA_AB_U  = nullptr;
    IP_AB_Um = nullptr;
    EA_AB_Um = nullptr;
    IP_AB_Up = nullptr;
    EA_AB_Up = nullptr;
    
    
    

}


int PT_tensors_rel::set_par(RI_complex_data * ext_RI,
                        double * ext_eps,
                        int ext_n_c, int ext_n_a, int ext_n_v, 
                        double * ext_H_AV, double * ext_H_CA, double * ext_H_CV,
                        double ext_edshift){
    
    RI      = ext_RI     ; 
    n_c     = ext_n_c    ; 
    n_a     = ext_n_a    ;
    n_v     = ext_n_v    ;
    H_AV_r    = ext_H_AV   ;
    H_CA_r    = ext_H_CA   ;
    H_CV_r    = ext_H_CV   ;
    
    H_AV_i    = ext_H_AV+2*n_v*2*n_a;
    H_CA_i    = ext_H_CA+2*n_c*2*n_a;
    H_CV_i    = ext_H_CV+2*n_c*2*n_v;
    
    edshift = ext_edshift;
        
    e_c     = ext_eps        ; 
    e_a     = ext_eps+2*n_c    ; 
    e_v     = ext_eps+2*n_c+2*n_a;
    
    
    
    
    e_IP    = new double[2*n_a];
    e_EA    = new double[2*n_a];
    IP_U    = new double[2*n_a*2*n_a];
    EA_U    = new double[2*n_a*2*n_a];
    IP_Um   = new double[2*n_a*2*n_a];
    EA_Um   = new double[2*n_a*2*n_a];
    IP_Up   = new double[2*n_a*2*n_a];
    EA_Up   = new double[2*n_a*2*n_a];
    e_IP_2  = new double[2*n_a*2*n_a];
    e_EA_2  = new double[2*n_a*2*n_a];
    IP_2_U  = new double[2*n_a*2*n_a*2*n_a*2*n_a];
    EA_2_U  = new double[2*n_a*2*n_a*2*n_a*2*n_a];
    IP_2_Um = new double[2*n_a*2*n_a*2*n_a*2*n_a];
    EA_2_Um = new double[2*n_a*2*n_a*2*n_a*2*n_a];
    IP_2_Up = new double[2*n_a*2*n_a*2*n_a*2*n_a];
    EA_2_Up = new double[2*n_a*2*n_a*2*n_a*2*n_a];

    e_IP_AB  = new double[2*n_a*2*n_a];
    e_EA_AB  = new double[2*n_a*2*n_a];
    IP_AB_U  = new double[2*n_a*2*n_a*2*n_a*2*n_a];
    EA_AB_U  = new double[2*n_a*2*n_a*2*n_a*2*n_a];
    IP_AB_Um = new double[2*n_a*2*n_a*2*n_a*2*n_a];
    EA_AB_Um = new double[2*n_a*2*n_a*2*n_a*2*n_a];
    IP_AB_Up = new double[2*n_a*2*n_a*2*n_a*2*n_a];
    EA_AB_Up = new double[2*n_a*2*n_a*2*n_a*2*n_a];
        

    
    RF_PH_r  = new double  [2*n_a*2*n_a                        ];
    RF_PH_i  = new double  [2*n_a*2*n_a                        ];
    RF_PV    = new double  [2*n_a*2*n_a*2*n_a*2*n_a            ];
    RF_PV_r  = new double  [2*n_a*2*n_a*2*n_a*2*n_a            ];
    RF_PV_i  = new double  [2*n_a*2*n_a*2*n_a*2*n_a            ];
    RF_P3_r  = new double  [2*n_a*2*n_a*2*n_a*2*n_a*2*n_a*2*n_a];
    RF_P3_i  = new double  [2*n_a*2*n_a*2*n_a*2*n_a*2*n_a*2*n_a];
    
    return 0;
}


int PT_tensors_rel::E2_calc_EE(){
    
    set_zero_matr(RF_P3_r, 2*n_a*2*n_a*2*n_a*2*n_a*2*n_a*2*n_a);
    set_zero_matr(RF_P3_i, 2*n_a*2*n_a*2*n_a*2*n_a*2*n_a*2*n_a);
    
    set_zero_matr(RF_PV_r, 2*n_a*2*n_a*2*n_a*2*n_a);
    set_zero_matr(RF_PV_i, 2*n_a*2*n_a*2*n_a*2*n_a);
    
    set_zero_matr(RF_PH_r, 2*n_a*2*n_a);
    set_zero_matr(RF_PH_i, 2*n_a*2*n_a);
    
//     set_zero_matr(RF_PS, N_fit);
    RF_PS=0;
    
    calc_EE_2_CCVV();printf_timer("CCVV res table");fprintf(out_stream,"\n");fflush(out_stream);//+
    calc_EE_2_CAVV();printf_timer("CAVV res table");fprintf(out_stream,"\n");fflush(out_stream);//??
    calc_EE_2_AAVV();printf_timer("AAVV res table");fprintf(out_stream,"\n");fflush(out_stream);//+
    calc_EE_2_CCAV();printf_timer("CCAV res table");fprintf(out_stream,"\n");fflush(out_stream);//+
    calc_EE_2_CCAA();printf_timer("CCAA res table");fprintf(out_stream,"\n");fflush(out_stream);//+
    calc_EE_2_CV  ();printf_timer("CV   res table");fprintf(out_stream,"\n");fflush(out_stream);//?
    calc_EE_2_AV  ();printf_timer("AV   res table");fprintf(out_stream,"\n");fflush(out_stream);//
    calc_EE_2_CA  ();printf_timer("CA   res table");fprintf(out_stream,"\n");fflush(out_stream);//?
//     
    printf("const=%e\n\n", RF_PS);
    printf("M1r:\n");
    PrintMatr(RF_PH_r,2*n_a,2*n_a,0);
    printf("M1i:\n");
    PrintMatr(RF_PH_i,2*n_a,2*n_a,0);
//     printf("M2_r:\n");
//     PrintMatr(RF_PV_r,2*n_a*2*n_a,2*n_a*2*n_a,0);
//     printf("M2_i:\n");
//     PrintMatr(RF_PV_i,2*n_a*2*n_a,2*n_a*2*n_a,0);
    
    
//     exit(0);
    
    
    symmetrization(RF_PH_r, 2*n_a);
    symmetrization(RF_PV_r, 2*n_a*2*n_a);
    symmetrization(RF_P3_r, 2*n_a*2*n_a*2*n_a);
    
    anti_symmetrization(RF_PH_i, 2*n_a);
    anti_symmetrization(RF_PV_i, 2*n_a*2*n_a);
    anti_symmetrization(RF_P3_i, 2*n_a*2*n_a*2*n_a);
    
//     double * RF_P3_tmp = new double  [n_a*n_a*n_a*n_a*n_a*n_a];
//     for(int t=0; t<n_a; t++)
//     for(int u=0; u<n_a; u++)
//     for(int v=0; v<n_a; v++)
//     for(int w=0; w<n_a; w++)
//     for(int x=0; x<n_a; x++)
//     for(int y=0; y<n_a; y++){
//         RF_P3_tmp[((((t*n_a+u)*n_a+v)*n_a+w)*n_a+x)*n_a+y]=
//         RF_P3_AB    [((((t*n_a+u)*n_a+x)*n_a+v)*n_a+w)*n_a+y];
//     }
//     memcpy(RF_P3_AB,RF_P3_tmp,n_a*n_a*n_a*n_a*n_a*n_a*sizeof(double));
//     
//     delete[]RF_P3_tmp;
    
    return 0;
    
}

int PT_tensors_rel::calc_EE_2_CV(){
    
    double dE;
    double Vr,Vi;
    double TUr,TUi,VWr,VWi;
    double Hr, Hi;
    
    double * RF_MMr;    // (J-K) * (J-K)
    double * RF_MMi;    // (J-K) * (J-K)
    
    RF_MMr = new double  [2*n_a*2*n_a*2*n_a*2*n_a];
    RF_MMi = new double  [2*n_a*2*n_a*2*n_a*2*n_a];
    set_zero_matr(RF_MMr, 2*n_a*2*n_a*2*n_a*2*n_a);
    set_zero_matr(RF_MMi, 2*n_a*2*n_a*2*n_a*2*n_a);
    
    double * RF_HMr;
    double * RF_MHr;
    double * RF_HMi;
    double * RF_MHi;
    RF_HMr = new double  [2*n_a*2*n_a];
    RF_MHr = new double  [2*n_a*2*n_a];
    RF_HMi = new double  [2*n_a*2*n_a];
    RF_MHi = new double  [2*n_a*2*n_a];
    set_zero_matr(RF_HMr, 2*n_a*2*n_a);
    set_zero_matr(RF_MHr, 2*n_a*2*n_a);
    set_zero_matr(RF_HMi, 2*n_a*2*n_a);
    set_zero_matr(RF_MHi, 2*n_a*2*n_a);
    
    double RF_Hr=0;
    double RF_Hi=0;
//     RF_H = new double  [2*n_a];
//     set_zero_matr(RF_H, 2*n_a);
    
    double * VCAAr;
    double * VCAAi;
    VCAAr = new double[2*n_c*2*n_v*2*n_a*2*n_a];
    VCAAi = new double[2*n_c*2*n_v*2*n_a*2*n_a];
    
    //<a*t*|iu>=(a*i|t*u)
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        2*n_c*2*n_v,2*n_a*2*n_a,RI->aux_n_ao,1.0,
                        RI->VC_RI_M_r,RI->aux_n_ao,
                        RI->AA_RI_M_r,RI->aux_n_ao,0.0,
                        VCAAr,2*n_a*2*n_a);
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        2*n_c*2*n_v,2*n_a*2*n_a,RI->aux_n_ao,-1.0,
                        RI->VC_RI_M_i,RI->aux_n_ao,
                        RI->AA_RI_M_i,RI->aux_n_ao,1.0,
                        VCAAr,2*n_a*2*n_a);
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        2*n_c*2*n_v,2*n_a*2*n_a,RI->aux_n_ao,1.0,
                        RI->VC_RI_M_r,RI->aux_n_ao,
                        RI->AA_RI_M_i,RI->aux_n_ao,0.0,
                        VCAAi,2*n_a*2*n_a);
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        2*n_c*2*n_v,2*n_a*2*n_a,RI->aux_n_ao,1.0,
                        RI->VC_RI_M_i,RI->aux_n_ao,
                        RI->AA_RI_M_r,RI->aux_n_ao,1.0,
                        VCAAi,2*n_a*2*n_a);
    double * CAVAr;
    double * CAVAi;
    CAVAr = new double[2*n_a*2*n_c*2*n_a*2*n_v];
    CAVAi = new double[2*n_a*2*n_c*2*n_a*2*n_v];
    
    //<a*t*|ui>=(a*u|it*)=\sum[(a*u|J)(i*t|J)*
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        2*n_a*2*n_c,2*n_a*2*n_v,RI->aux_n_ao,1.0,
                        RI->CA_RI_M_r,RI->aux_n_ao,
                        RI->VA_RI_M_r,RI->aux_n_ao,0.0,
                        CAVAr,2*n_a*2*n_v);
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        2*n_a*2*n_c,2*n_a*2*n_v,RI->aux_n_ao,1.0,
                        RI->CA_RI_M_i,RI->aux_n_ao,
                        RI->VA_RI_M_i,RI->aux_n_ao,1.0,
                        CAVAr,2*n_a*2*n_v);
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        2*n_a*2*n_c,2*n_a*2*n_v,RI->aux_n_ao,1.0,
                        RI->CA_RI_M_r,RI->aux_n_ao,
                        RI->VA_RI_M_i,RI->aux_n_ao,0.0,
                        CAVAi,2*n_a*2*n_v);
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        2*n_a*2*n_c,2*n_a*2*n_v,RI->aux_n_ao,-1.0,
                        RI->CA_RI_M_i,RI->aux_n_ao,
                        RI->VA_RI_M_r,RI->aux_n_ao,1.0,
                        CAVAi,2*n_a*2*n_v);
//     double * VCAA_J;
//     VCAA_J = new double[2*n_c*2*n_v];
//     double * VCAA_K;
//     VCAA_K = new double[2*n_c*2*n_v];

    
    for(int i=0; i<2*n_c; i++)
    for(int a=0; a<2*n_v; a++)
    for(int t=0; t<2*n_a; t++)
    for(int u=0; u<2*n_a; u++){
        dE=e_c[i]+e_a[t]-e_v[a]-e_a[u];//is taken negative
        dE=dE/(dE*dE+edshift);
        
        for(int v=0; v<2*n_a; v++)
        for(int w=0; w<2*n_a; w++){
            //(<i*t*|au>-<i*u*|at>)(<a*v*|iw>-<a*v*|wi>)
            TUr=(VCAAr[((a*2*n_c+i)*2*n_a+u)*2*n_a+t]-  //<it|au>=(ai|tu)
                 CAVAr[((i*2*n_a+u)*2*n_v+a)*2*n_a+t]); //<it|ua>=(iu|at)
            VWr=(VCAAr[((a*2*n_c+i)*2*n_a+w)*2*n_a+v]-  //<av|iw>=(ai|wv)
                 CAVAr[((i*2*n_a+v)*2*n_v+a)*2*n_a+w]); //<av|wi>=(iv|aw)
            
            TUi=(VCAAi[((a*2*n_c+i)*2*n_a+u)*2*n_a+t]-  //<it|au>=(ai|tu)
                 CAVAi[((i*2*n_a+u)*2*n_v+a)*2*n_a+t]); //<it|ua>=(iu|at)
            VWi=(VCAAi[((a*2*n_c+i)*2*n_a+w)*2*n_a+v]-  //<av|iw>=(ai|wv)
                 CAVAi[((i*2*n_a+v)*2*n_v+a)*2*n_a+w]); //<av|wi>=(iv|aw)
            
            
            RF_MMr[((t*2*n_a+u)*2*n_a+v)*2*n_a+w]+=(TUr*VWr-TUi*VWi)*dE;
            RF_MMi[((t*2*n_a+u)*2*n_a+v)*2*n_a+w]+=(TUi*VWr+TUr*VWi)*dE;
        }
    }
    

    
    for(int i=0; i<2*n_c; i++)
    for(int a=0; a<2*n_v; a++)
    for(int t=0; t<2*n_a; t++)
    for(int u=0; u<2*n_a; u++){
        //(<it|au>-<iu|at>)*H[i,a]
        TUr=(VCAAr[((a*2*n_c+i)*2*n_a+u)*2*n_a+t]-  //<it|au>=(ai|tu)
             CAVAr[((i*2*n_a+u)*2*n_v+a)*2*n_a+t]); //<it|ua>=(iu|at)
        Hr = H_CV_r[i*2*n_v+a]; 
        TUi=(VCAAi[((a*2*n_c+i)*2*n_a+u)*2*n_a+t]-  //<it|au>=(ai|tu)
             CAVAi[((i*2*n_a+u)*2*n_v+a)*2*n_a+t]); //<it|ua>=(iu|at)
        Hi = H_CV_i[i*2*n_v+a]; 
        dE=e_c[i]+e_a[t]-e_v[a]-e_a[u];
        
        RF_MHr[t*2*n_a+u]+=(TUr*Hr+TUi*Hi)*dE/(dE*dE+edshift);
        RF_MHi[t*2*n_a+u]+=(TUi*Hr-TUr*Hi)*dE/(dE*dE+edshift);
    }
    
    for(int i=0; i<2*n_c; i++)
    for(int a=0; a<2*n_v; a++)
    for(int t=0; t<2*n_a; t++)
    for(int u=0; u<2*n_a; u++){
        //(<it|au>-<iu|at>)*H[i,a]
        TUr=(VCAAr[((a*2*n_c+i)*2*n_a+u)*2*n_a+t]-  //<it|au>=(ai|tu)
             CAVAr[((i*2*n_a+u)*2*n_v+a)*2*n_a+t]); //<it|ua>=(iu|at)
        Hr = H_CV_r[i*2*n_v+a]; 
        TUi=(VCAAi[((a*2*n_c+i)*2*n_a+u)*2*n_a+t]-  //<it|au>=(ai|tu)
             CAVAi[((i*2*n_a+u)*2*n_v+a)*2*n_a+t]); //<it|ua>=(iu|at)
        Hi = H_CV_i[i*2*n_v+a]; 
        
        dE=e_c[i]-e_v[a];
        
        RF_HMr[t*2*n_a+u]+=(TUr*Hr+TUi*Hi)*dE/(dE*dE+edshift);
        RF_HMi[t*2*n_a+u]+=(TUr*Hi-TUi*Hr)*dE/(dE*dE+edshift);
    }

    
    for(int i=0; i<2*n_c; i++)
    for(int a=0; a<2*n_v; a++){
        //H[i,a]*H[i,a]
        Vr=H_CV_r[i*2*n_v+a]*H_CV_r[i*2*n_v+a]+
           H_CV_i[i*2*n_v+a]*H_CV_i[i*2*n_v+a]; 
        dE=e_c[i]-e_v[a];
        RF_Hr+=Vr*dE/(dE*dE+edshift);
        
    }
    
    for(int t=0; t<2*n_a; t++)
    for(int u=0; u<2*n_a; u++)
    for(int v=0; v<2*n_a; v++)
    for(int w=0; w<2*n_a; w++)
        RF_PV_r[((u*2*n_a+t)*2*n_a+v)*2*n_a+w]+=+RF_MMr[((t*2*n_a+v)*2*n_a+u)*2*n_a+w]
                                                -RF_MMr[((u*2*n_a+v)*2*n_a+t)*2*n_a+w]
                                                -RF_MMr[((t*2*n_a+w)*2*n_a+u)*2*n_a+v]
                                                +RF_MMr[((u*2*n_a+w)*2*n_a+t)*2*n_a+v]
                                                ;
    for(int t=0; t<2*n_a; t++)
    for(int u=0; u<2*n_a; u++)
    for(int v=0; v<2*n_a; v++)
    for(int w=0; w<2*n_a; w++)
        RF_PV_i[((u*2*n_a+t)*2*n_a+v)*2*n_a+w]+=+RF_MMi[((t*2*n_a+v)*2*n_a+u)*2*n_a+w]
                                                -RF_MMi[((u*2*n_a+v)*2*n_a+t)*2*n_a+w]
                                                -RF_MMi[((t*2*n_a+w)*2*n_a+u)*2*n_a+v]
                                                +RF_MMi[((u*2*n_a+w)*2*n_a+t)*2*n_a+v]
                                                ;
                                                       
    
    
//     set_zero_matr(RF_PH_r,N_fit*2*n_a*2*n_a);
    for(int t=0; t<2*n_a; t++)
    for(int w=0; w<2*n_a; w++){
        RF_PH_r[t*2*n_a+w]+= RF_HMr[w*2*n_a+t]
                            +RF_MHr[t*2*n_a+w];
        for(int u=0; u<2*n_a; u++)
             RF_PH_r[t*2*n_a+w]+= RF_MMr[((t*2*n_a+u)*2*n_a+u)*2*n_a+w];
    }
    
    for(int t=0; t<2*n_a; t++)
    for(int w=0; w<2*n_a; w++){
        RF_PH_i[t*2*n_a+w]+= RF_HMi[w*2*n_a+t]
                            +RF_MHi[t*2*n_a+w];
        for(int u=0; u<2*n_a; u++)
             RF_PH_i[t*2*n_a+w]+= RF_MMi[((t*2*n_a+u)*2*n_a+u)*2*n_a+w];
    }
    
    
    RF_PS+=RF_Hr;
    
    delete[] RF_MMr ;
    delete[] RF_MMi ;
    delete[] RF_HMr ;
    delete[] RF_MHr ;
    delete[] RF_HMi ;
    delete[] RF_MHi ;
    delete[] VCAAr  ;
    delete[] CAVAr  ;
    delete[] VCAAi  ;
    delete[] CAVAi  ;
//     delete[] VCAA_J;
//     delete[] VCAA_K;
    
    return 0;
    
}

int PT_tensors_rel::calc_EE_2_AV(){
    
    
    
    int aux_n_ao = RI->aux_n_ao;
    double * VA_r=RI->VA_RI_M_r;
    double * CC_r=RI->CC_RI_M_r;
    double * VC_r=RI->VC_RI_M_r;
    double * CA_r=RI->CA_RI_M_r;
    double * AA_r=RI->AA_RI_M_r;
    double * VA_i=RI->VA_RI_M_i;
    double * CC_i=RI->CC_RI_M_i;
    double * VC_i=RI->VC_RI_M_i;
    double * CA_i=RI->CA_RI_M_i;
    double * AA_i=RI->AA_RI_M_i;
    
    double * VAAA_r;
    double * VAAA_i;
    VAAA_r = new double[2*n_v*2*n_a*2*n_a*2*n_a];
    VAAA_i = new double[2*n_v*2*n_a*2*n_a*2*n_a];
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        2*n_v*2*n_a,2*n_a*2*n_a,aux_n_ao,1.0,
                        VA_r,aux_n_ao,
                        AA_r,aux_n_ao,0.0,
                        VAAA_r,2*n_a*2*n_a);
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        2*n_v*2*n_a,2*n_a*2*n_a,aux_n_ao,-1.0,
                        VA_i,aux_n_ao,
                        AA_i,aux_n_ao,1.0,
                        VAAA_r,2*n_a*2*n_a);
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        2*n_v*2*n_a,2*n_a*2*n_a,aux_n_ao,1.0,
                        VA_i,aux_n_ao,
                        AA_r,aux_n_ao,0.0,
                        VAAA_r,2*n_a*2*n_a);
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        2*n_v*2*n_a,2*n_a*2*n_a,aux_n_ao,1.0,
                        VA_r,aux_n_ao,
                        AA_i,aux_n_ao,1.0,
                        VAAA_r,2*n_a*2*n_a);
    
    
    double * RF_P3_HM_r;
    double * RF_P3_MH_r;
    double * RF_P3_HM_i;
    double * RF_P3_MH_i;
//     double * RF_P3_HJ;
//     double * RF_P3_JH;
    
//     RF_P3_JJ = new double[2*n_a*2*n_a*2*n_a*2*n_a*2*n_a*2*n_a];
//     RF_P3_JM = new double[2*n_a*2*n_a*2*n_a*2*n_a*2*n_a*2*n_a];
//     RF_P3_MJ = new double[2*n_a*2*n_a*2*n_a*2*n_a*2*n_a*2*n_a];
    
    RF_P3_HM_r = new double[2*n_a*2*n_a*2*n_a*2*n_a];
    RF_P3_MH_r = new double[2*n_a*2*n_a*2*n_a*2*n_a];
    RF_P3_HM_i = new double[2*n_a*2*n_a*2*n_a*2*n_a];
    RF_P3_MH_i = new double[2*n_a*2*n_a*2*n_a*2*n_a];
//     RF_P3_HJ = new double[2*n_a*2*n_a*2*n_a*2*n_a];
//     RF_P3_JH = new double[2*n_a*2*n_a*2*n_a*2*n_a];
    
//     set_zero_matr(RF_P3_JJ, 2*n_a*2*n_a*2*n_a*2*n_a*2*n_a*2*n_a);
//     set_zero_matr(RF_P3_JM, 2*n_a*2*n_a*2*n_a*2*n_a*2*n_a*2*n_a);
//     set_zero_matr(RF_P3_MJ, 2*n_a*2*n_a*2*n_a*2*n_a*2*n_a*2*n_a);
    set_zero_matr(RF_P3_HM_r, 2*n_a*2*n_a*2*n_a*2*n_a);
    set_zero_matr(RF_P3_MH_r, 2*n_a*2*n_a*2*n_a*2*n_a);
    set_zero_matr(RF_P3_HM_i, 2*n_a*2*n_a*2*n_a*2*n_a);
    set_zero_matr(RF_P3_MH_i, 2*n_a*2*n_a*2*n_a*2*n_a);
//     set_zero_matr(RF_P3_HJ, 2*n_a*2*n_a*2*n_a*2*n_a);
//     set_zero_matr(RF_P3_JH, 2*n_a*2*n_a*2*n_a*2*n_a);
    
    
#pragma omp parallel
    {
    //MM == JK
        double Vr, Vi;
        double Tr,Ti,Wr,Wi;
        double dE;
        double H1a_r,H1a_i;
        double H2a_r,H2a_i;
        
        int u,v;
        int th_id = omp_get_thread_num();
        
        for(int vu=th_id; vu<2*n_a*2*n_a; vu+=num_threads){
            v=vu/(2*n_a);
            u=vu%(2*n_a);
            for(int a=0; a<2*n_v; a++)
            for(int t=0; t<2*n_a; t++){
                dE=-e_v[a]+e_a[t]+e_a[u]-e_a[v];
                dE=dE/(dE*dE+edshift);
                
                for(int y=0; y<2*n_a; y++)
                for(int x=0; x<2*n_a; x++)
                for(int w=0; w<2*n_a; w++){
                    //<tu||av>=<tu|av>-<tu||va>=(ta|vu)-(tv|au)
                    //<xy||aw>=<xy|aw>-<xy||wa>=(xa|yw)-(xw|ay)
                    Tr= VAAA_r[((a*2*n_a+t)*2*n_a+v)*2*n_a+u]-
                        VAAA_r[((a*2*n_a+u)*2*n_a+v)*2*n_a+t];
                    Ti= VAAA_i[((a*2*n_a+t)*2*n_a+v)*2*n_a+u]-
                        VAAA_i[((a*2*n_a+u)*2*n_a+v)*2*n_a+t];
                    Wr= VAAA_r[((a*2*n_a+x)*2*n_a+w)*2*n_a+y]-
                        VAAA_r[((a*2*n_a+y)*2*n_a+w)*2*n_a+x];
                    Wi= VAAA_i[((a*2*n_a+x)*2*n_a+w)*2*n_a+y]-
                        VAAA_i[((a*2*n_a+y)*2*n_a+w)*2*n_a+x];
                    
                    RF_P3_r[((((t*2*n_a+u)*2*n_a+v)*2*n_a+y)*2*n_a+x)*2*n_a+w]+=(Tr*Wr+Ti*Wi)*dE;
                    RF_P3_i[((((t*2*n_a+u)*2*n_a+v)*2*n_a+y)*2*n_a+x)*2*n_a+w]+=(Ti*Wr-Tr*Wi)*dE;
                    
                }
            }
        
        }    
        //H
        for(int a=0; a<2*n_v; a++)
        for(int t=th_id; t<2*n_a; t+=num_threads){
            dE=-e_v[a]+e_a[t];
            dE=dE/(dE*dE+edshift);
            
            for(int w=0; w<2*n_a; w++){
                H1a_r=H_AV_r[t*2*n_v+a];
                H2a_r=H_AV_r[w*2*n_v+a];
                H1a_i=H_AV_i[t*2*n_v+a];
                H2a_i=H_AV_i[w*2*n_v+a];
                            
                RF_PH_r[t*2*n_a+w]+=(H1a_r*H2a_r+H1a_i*H2a_i)*dE;
                RF_PH_i[t*2*n_a+w]+=(H1a_i*H2a_r-H1a_r*H2a_i)*dE;
                
            }
        }
        //HM
        for(int a=0; a<2*n_v; a++)
        for(int t=th_id; t<2*n_a; t+=num_threads){
            dE=-e_v[a]+e_a[t];
            dE=dE/(dE*dE+edshift);
            
            for(int y=0; y<2*n_a; y++)
            for(int x=0; x<2*n_a; x++)
            for(int w=0; w<2*n_a; w++){
                H1a_r=H_AV_r[t*2*n_v+a];
                H1a_i=H_AV_i[t*2*n_v+a];
                //<aw|xy>-<aw|yx>=(yw|xa)-(xw|ya)
                Wr=VAAA_r[((a*2*n_a+x)*2*n_a+w)*2*n_a+y]-
                   VAAA_r[((a*2*n_a+y)*2*n_a+w)*2*n_a+x];
                Wi=VAAA_i[((a*2*n_a+x)*2*n_a+w)*2*n_a+y]-
                   VAAA_i[((a*2*n_a+y)*2*n_a+w)*2*n_a+x];
                
                RF_P3_HM_r[((t*2*n_a+y)*2*n_a+x)*2*n_a+w]+=(H1a_r*Wr+H1a_i*Wi)*dE;
                RF_P3_HM_i[((t*2*n_a+y)*2*n_a+x)*2*n_a+w]+=(H1a_i*Wr-H1a_r*Wi)*dE;
                
            }
        }
        //MH
        for(int vu=th_id; vu<2*n_a*2*n_a; vu+=num_threads){
            v=vu/(2*n_a);
            u=vu%(2*n_a);
            for(int a=0; a<2*n_v; a++)
            for(int t=0; t<2*n_a; t++){
                dE=-e_v[a]+e_a[t]+e_a[u]-e_a[v];
                dE=dE/(dE*dE+edshift);
                
                for(int w=0; w<2*n_a; w++){
                    H2a_r=H_AV_r[w*2*n_v+a];
                    H2a_i=H_AV_i[w*2*n_v+a];
                    //<aw|xy>-<aw|yx>=(yw|xa)-(xw|ya)
                    Tr= VAAA_r[((a*2*n_a+t)*2*n_a+v)*2*n_a+u]-
                        VAAA_r[((a*2*n_a+u)*2*n_a+v)*2*n_a+t];
                    Ti= VAAA_i[((a*2*n_a+t)*2*n_a+v)*2*n_a+u]-
                        VAAA_i[((a*2*n_a+u)*2*n_a+v)*2*n_a+t];
                    
                    RF_P3_MH_r[((t*2*n_a+u)*2*n_a+v)*2*n_a+w]+=(Tr*H2a_r+Ti*H2a_i)*dE;
                    RF_P3_MH_i[((t*2*n_a+u)*2*n_a+v)*2*n_a+w]+=(Ti*H2a_r-Tr*H2a_i)*dE;
                }
            }
        }
        
    }
    //RF_PV_JK
    for(int w=0; w<2*n_a; w++)
    for(int t=0; t<2*n_a; t++)
    for(int v=0; v<2*n_a; v++)
    for(int u=0; u<2*n_a; u++)
        RF_PV_r[((w*2*n_a+t)*2*n_a+v)*2*n_a+u]+= RF_P3_HM_r[((t*2*n_a+u)*2*n_a+v)*2*n_a+w]
                                                -RF_P3_HM_r[((w*2*n_a+u)*2*n_a+v)*2*n_a+t]
                                                -RF_P3_MH_r[((t*2*n_a+w)*2*n_a+v)*2*n_a+u]
                                                +RF_P3_MH_r[((t*2*n_a+w)*2*n_a+u)*2*n_a+v];
    for(int w=0; w<2*n_a; w++)
    for(int t=0; t<2*n_a; t++)
    for(int v=0; v<2*n_a; v++)
    for(int u=0; u<2*n_a; u++)
        RF_PV_i[((w*2*n_a+t)*2*n_a+v)*2*n_a+u]+= RF_P3_HM_i[((t*2*n_a+u)*2*n_a+v)*2*n_a+w]
                                                -RF_P3_HM_i[((w*2*n_a+u)*2*n_a+v)*2*n_a+t]
                                                -RF_P3_MH_i[((t*2*n_a+w)*2*n_a+v)*2*n_a+u]
                                                +RF_P3_MH_i[((t*2*n_a+w)*2*n_a+u)*2*n_a+v];
    
    delete[] VAAA_r;
    delete[] VAAA_i;
    delete[] RF_P3_HM_r;
    delete[] RF_P3_MH_r;
    delete[] RF_P3_HM_i;
    delete[] RF_P3_MH_i;
    
    return 0;
    
}

int PT_tensors_rel::calc_EE_2_CA(){
    
    double dE;
    double H1a_r,H1a_i;
    double H2a_r,H2a_i;
    double Vr,Vi;
    
    
    int aux_n_ao = RI->aux_n_ao;

    double * CC_r=RI->CC_RI_M_r;
    double * CC_i=RI->CC_RI_M_i;

    double * CA_r=RI->CA_RI_M_r;
    double * CA_i=RI->CA_RI_M_i;
    double * AA_r=RI->AA_RI_M_r;
    double * AA_i=RI->AA_RI_M_i;

    
    double * CAAA_r;
    double * CAAA_i;
    CAAA_r = new double[2*n_c*2*n_a*2*n_a*2*n_a];
    CAAA_i = new double[2*n_c*2*n_a*2*n_a*2*n_a];
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        2*n_c*2*n_a,2*n_a*2*n_a,aux_n_ao,1.0,
                        CA_r,aux_n_ao,
                        AA_r,aux_n_ao,0.0,
                        CAAA_r,2*n_a*2*n_a);
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        2*n_c*2*n_a,2*n_a*2*n_a,aux_n_ao,1.0,
                        CA_i,aux_n_ao,
                        AA_i,aux_n_ao,1.0,
                        CAAA_r,2*n_a*2*n_a);
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        2*n_c*2*n_a,2*n_a*2*n_a,aux_n_ao,-1.0,
                        CA_i,aux_n_ao,
                        AA_r,aux_n_ao,0.0,
                        CAAA_i,2*n_a*2*n_a);
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        2*n_c*2*n_a,2*n_a*2*n_a,aux_n_ao,1.0,
                        CA_r,aux_n_ao,
                        AA_i,aux_n_ao,1.0,
                        CAAA_i,2*n_a*2*n_a);
        
    
    double * RF_P3_VV_r;
    double * RF_P3_VH_r;
    double * RF_P3_HV_r;
    double * RF_P3_HH_r;
    double * RF_P3_VV_i;
    double * RF_P3_VH_i;
    double * RF_P3_HV_i;
    double * RF_P3_HH_i;
    
    
    RF_P3_VV_r = new double[2*n_a*2*n_a*2*n_a*2*n_a*2*n_a*2*n_a];
    RF_P3_VH_r = new double[2*n_a*2*n_a*2*n_a*2*n_a        ];
    RF_P3_HV_r = new double[2*n_a*2*n_a*2*n_a*2*n_a        ];
    RF_P3_HH_r = new double[2*n_a*2*n_a                ];
    RF_P3_VV_i = new double[2*n_a*2*n_a*2*n_a*2*n_a*2*n_a*2*n_a];
    RF_P3_VH_i = new double[2*n_a*2*n_a*2*n_a*2*n_a        ];
    RF_P3_HV_i = new double[2*n_a*2*n_a*2*n_a*2*n_a        ];
    RF_P3_HH_i = new double[2*n_a*2*n_a                ];
    
    set_zero_matr(RF_P3_VV_r, 2*n_a*2*n_a*2*n_a*2*n_a*2*n_a*2*n_a);
    set_zero_matr(RF_P3_VH_r, 2*n_a*2*n_a*2*n_a*2*n_a        );
    set_zero_matr(RF_P3_HV_r, 2*n_a*2*n_a*2*n_a*2*n_a        );
    set_zero_matr(RF_P3_HH_r, 2*n_a*2*n_a                );
    set_zero_matr(RF_P3_VV_i, 2*n_a*2*n_a*2*n_a*2*n_a*2*n_a*2*n_a);
    set_zero_matr(RF_P3_VH_i, 2*n_a*2*n_a*2*n_a*2*n_a        );
    set_zero_matr(RF_P3_HV_i, 2*n_a*2*n_a*2*n_a*2*n_a        );
    set_zero_matr(RF_P3_HH_i, 2*n_a*2*n_a                );
    
    
    
    //VV
    for(int v=0; v<2*n_a; v++)
    for(int u=0; u<2*n_a; u++)
    for(int i=0; i<2*n_c; i++)
    for(int t=0; t<2*n_a; t++){
        dE= e_c[i]+e_a[t]-e_a[u]-e_a[v];
        dE=dE/(dE*dE+edshift);
        
        for(int y=0; y<2*n_a; y++)
        for(int x=0; x<2*n_a; x++)
        for(int w=0; w<2*n_a; w++){
            
            Vr= CAAA_r[((i*2*n_a+v)*2*n_a+u)*2*n_a+t]*CAAA_r[((i*2*n_a+w)*2*n_a+x)*2*n_a+y]
               +CAAA_i[((i*2*n_a+v)*2*n_a+u)*2*n_a+t]*CAAA_i[((i*2*n_a+w)*2*n_a+x)*2*n_a+y];
            Vi= CAAA_i[((i*2*n_a+v)*2*n_a+u)*2*n_a+t]*CAAA_r[((i*2*n_a+w)*2*n_a+x)*2*n_a+y]
               -CAAA_r[((i*2*n_a+v)*2*n_a+u)*2*n_a+t]*CAAA_i[((i*2*n_a+w)*2*n_a+x)*2*n_a+y];
            
            RF_P3_VV_r[((((t*2*n_a+u)*2*n_a+v)*2*n_a+y)*2*n_a+x)*2*n_a+w]+=Vr*dE;
            RF_P3_VV_i[((((t*2*n_a+u)*2*n_a+v)*2*n_a+y)*2*n_a+x)*2*n_a+w]+=Vi*dE;
            
        }
    }
    //VH
    for(int v=0; v<2*n_a; v++)
    for(int u=0; u<2*n_a; u++)
    for(int i=0; i<2*n_c; i++)
    for(int t=0; t<2*n_a; t++){
        dE= e_c[i]+e_a[t]-e_a[u]-e_a[v];
        dE=dE/(dE*dE+edshift);
        for(int w=0; w<2*n_a; w++){
            H2a_r=H_CA_r[i*2*n_a+w];
            H2a_i=H_CA_i[i*2*n_a+w];
            Vr= CAAA_r[((i*2*n_a+v)*2*n_a+u)*2*n_a+t]*H2a_r
               -CAAA_i[((i*2*n_a+v)*2*n_a+u)*2*n_a+t]*H2a_i;
            Vi= CAAA_i[((i*2*n_a+v)*2*n_a+u)*2*n_a+t]*H2a_r
               +CAAA_r[((i*2*n_a+v)*2*n_a+u)*2*n_a+t]*H2a_i;
            
            RF_P3_VH_r[((t*2*n_a+u)*2*n_a+v)*2*n_a+w]+=Vr*dE;
            RF_P3_VH_i[((t*2*n_a+u)*2*n_a+v)*2*n_a+w]+=Vi*dE;
        }
    }

    //HV
    for(int i=0; i<2*n_c; i++)
    for(int t=0; t<2*n_a; t++){
        dE= e_c[i]-e_a[t];
        dE=dE/(dE*dE+edshift);
        
        for(int y=0; y<2*n_a; y++)
        for(int x=0; x<2*n_a; x++)
        for(int w=0; w<2*n_a; w++){
            H1a_r=H_CA_r[i*2*n_a+t];
            H1a_r=H_CA_r[i*2*n_a+t];
            
            Vr= H1a_r*CAAA_r[((i*2*n_a+w)*2*n_a+x)*2*n_a+y]
               -H1a_i*CAAA_i[((i*2*n_a+w)*2*n_a+x)*2*n_a+y];
            Vi=-H1a_i*CAAA_r[((i*2*n_a+w)*2*n_a+x)*2*n_a+y]
               -H1a_r*CAAA_i[((i*2*n_a+w)*2*n_a+x)*2*n_a+y];
            
               RF_P3_HV_r[((t*2*n_a+y)*2*n_a+x)*2*n_a+w]+=Vr*dE;
               RF_P3_HV_i[((t*2*n_a+y)*2*n_a+x)*2*n_a+w]+=Vi*dE;
        }
    }
    
    //HH
    for(int i=0; i<2*n_c; i++)
    for(int t=0; t<2*n_a; t++){
        dE= e_c[i]-e_a[t];
        dE=dE/(dE*dE+edshift);
        
        for(int u=0; u<2*n_a; u++){
            H1a_r=H_CA_r[i*2*n_a+t];
            H1a_i=H_CA_i[i*2*n_a+t];
            H2a_r=H_CA_r[i*2*n_a+u];
            H2a_i=H_CA_i[i*2*n_a+u];
//             printf("%e %e %e %e\n",H1a_r,H2a_r,H1a_i,H2a_i);
//             getchar();
            Vr= H1a_r*H2a_r+H1a_i*H2a_i;
            Vi=-H1a_i*H2a_r+H1a_r*H2a_i;
            RF_P3_HH_r[t*2*n_a+u]+=Vr*dE;
            RF_P3_HH_i[t*2*n_a+u]+=Vi*dE;
        }
    }
//     PrintMatr(RF_P3_HH_r,2*n_a,2*n_a,0);
//     exit(0);
    
    
    for(int t=0; t<2*n_a; t++)
    for(int v=0; v<2*n_a; v++)
    for(int u=0; u<2*n_a; u++)
    for(int w=0; w<2*n_a; w++)
    for(int x=0; x<2*n_a; x++)
    for(int y=0; y<2*n_a; y++)
        RF_P3_r[((((t*2*n_a+u)*2*n_a+v)*2*n_a+y)*2*n_a+x)*2*n_a+w]+= RF_P3_VV_r[((((t*2*n_a+v)*2*n_a+x)*2*n_a+y)*2*n_a+w)*2*n_a+u]
                                                                    -RF_P3_VV_r[((((u*2*n_a+v)*2*n_a+x)*2*n_a+y)*2*n_a+w)*2*n_a+t]
                                                                    -RF_P3_VV_r[((((t*2*n_a+v)*2*n_a+y)*2*n_a+x)*2*n_a+w)*2*n_a+u]
                                                                    +RF_P3_VV_r[((((u*2*n_a+v)*2*n_a+y)*2*n_a+x)*2*n_a+w)*2*n_a+t];
    for(int t=0; t<2*n_a; t++)
    for(int v=0; v<2*n_a; v++)
    for(int u=0; u<2*n_a; u++)
    for(int w=0; w<2*n_a; w++)
    for(int x=0; x<2*n_a; x++)
    for(int y=0; y<2*n_a; y++)
        RF_P3_i[((((t*2*n_a+u)*2*n_a+v)*2*n_a+y)*2*n_a+x)*2*n_a+w]+= RF_P3_VV_i[((((t*2*n_a+v)*2*n_a+x)*2*n_a+y)*2*n_a+w)*2*n_a+u]
                                                                    -RF_P3_VV_i[((((u*2*n_a+v)*2*n_a+x)*2*n_a+y)*2*n_a+w)*2*n_a+t]
                                                                    -RF_P3_VV_i[((((t*2*n_a+v)*2*n_a+y)*2*n_a+x)*2*n_a+w)*2*n_a+u]
                                                                    +RF_P3_VV_i[((((u*2*n_a+v)*2*n_a+y)*2*n_a+x)*2*n_a+w)*2*n_a+t];
    
    for(int t=0; t<2*n_a; t++)
    for(int v=0; v<2*n_a; v++)
    for(int u=0; u<2*n_a; u++)
    for(int w=0; w<2*n_a; w++)
    for(int x=0; x<2*n_a; x++)
        RF_PV_r[((u*2*n_a+t)*2*n_a+v)*2*n_a+w]+=-RF_P3_VV_r[((((t*2*n_a+v)*2*n_a+x)*2*n_a+w)*2*n_a+x)*2*n_a+u]
                                                +RF_P3_VV_r[((((u*2*n_a+v)*2*n_a+x)*2*n_a+w)*2*n_a+x)*2*n_a+t]
                                                +RF_P3_VV_r[((((t*2*n_a+w)*2*n_a+x)*2*n_a+v)*2*n_a+x)*2*n_a+u]
                                                -RF_P3_VV_r[((((u*2*n_a+w)*2*n_a+x)*2*n_a+v)*2*n_a+x)*2*n_a+t]
                                                -RF_P3_VV_r[((((t*2*n_a+x)*2*n_a+v)*2*n_a+w)*2*n_a+u)*2*n_a+x]
                                                +RF_P3_VV_r[((((t*2*n_a+v)*2*n_a+x)*2*n_a+w)*2*n_a+u)*2*n_a+x]
                                                +RF_P3_VV_r[((((u*2*n_a+x)*2*n_a+v)*2*n_a+w)*2*n_a+t)*2*n_a+x]
                                                -RF_P3_VV_r[((((u*2*n_a+v)*2*n_a+x)*2*n_a+w)*2*n_a+t)*2*n_a+x]
                                                +RF_P3_VV_r[((((t*2*n_a+x)*2*n_a+w)*2*n_a+v)*2*n_a+u)*2*n_a+x]
                                                -RF_P3_VV_r[((((t*2*n_a+w)*2*n_a+x)*2*n_a+v)*2*n_a+u)*2*n_a+x]
                                                -RF_P3_VV_r[((((u*2*n_a+x)*2*n_a+w)*2*n_a+v)*2*n_a+t)*2*n_a+x]
                                                +RF_P3_VV_r[((((u*2*n_a+w)*2*n_a+x)*2*n_a+v)*2*n_a+t)*2*n_a+x]
                                                +RF_P3_VV_r[((((t*2*n_a+v)*2*n_a+x)*2*n_a+w)*2*n_a+u)*2*n_a+x]
                                                -RF_P3_VV_r[((((u*2*n_a+v)*2*n_a+x)*2*n_a+w)*2*n_a+t)*2*n_a+x]
                                                -RF_P3_VV_r[((((t*2*n_a+w)*2*n_a+x)*2*n_a+v)*2*n_a+u)*2*n_a+x]
                                                +RF_P3_VV_r[((((u*2*n_a+w)*2*n_a+x)*2*n_a+v)*2*n_a+t)*2*n_a+x];
    for(int t=0; t<2*n_a; t++)
    for(int v=0; v<2*n_a; v++)
    for(int u=0; u<2*n_a; u++)
    for(int w=0; w<2*n_a; w++)
    for(int x=0; x<2*n_a; x++)
        RF_PV_i[((u*2*n_a+t)*2*n_a+v)*2*n_a+w]+=-RF_P3_VV_i[((((t*2*n_a+v)*2*n_a+x)*2*n_a+w)*2*n_a+x)*2*n_a+u]
                                                +RF_P3_VV_i[((((u*2*n_a+v)*2*n_a+x)*2*n_a+w)*2*n_a+x)*2*n_a+t]
                                                +RF_P3_VV_i[((((t*2*n_a+w)*2*n_a+x)*2*n_a+v)*2*n_a+x)*2*n_a+u]
                                                -RF_P3_VV_i[((((u*2*n_a+w)*2*n_a+x)*2*n_a+v)*2*n_a+x)*2*n_a+t]
                                                -RF_P3_VV_i[((((t*2*n_a+x)*2*n_a+v)*2*n_a+w)*2*n_a+u)*2*n_a+x]
                                                +RF_P3_VV_i[((((t*2*n_a+v)*2*n_a+x)*2*n_a+w)*2*n_a+u)*2*n_a+x]
                                                +RF_P3_VV_i[((((u*2*n_a+x)*2*n_a+v)*2*n_a+w)*2*n_a+t)*2*n_a+x]
                                                -RF_P3_VV_i[((((u*2*n_a+v)*2*n_a+x)*2*n_a+w)*2*n_a+t)*2*n_a+x]
                                                +RF_P3_VV_i[((((t*2*n_a+x)*2*n_a+w)*2*n_a+v)*2*n_a+u)*2*n_a+x]
                                                -RF_P3_VV_i[((((t*2*n_a+w)*2*n_a+x)*2*n_a+v)*2*n_a+u)*2*n_a+x]
                                                -RF_P3_VV_i[((((u*2*n_a+x)*2*n_a+w)*2*n_a+v)*2*n_a+t)*2*n_a+x]
                                                +RF_P3_VV_i[((((u*2*n_a+w)*2*n_a+x)*2*n_a+v)*2*n_a+t)*2*n_a+x]
                                                +RF_P3_VV_i[((((t*2*n_a+v)*2*n_a+x)*2*n_a+w)*2*n_a+u)*2*n_a+x]
                                                -RF_P3_VV_i[((((u*2*n_a+v)*2*n_a+x)*2*n_a+w)*2*n_a+t)*2*n_a+x]
                                                -RF_P3_VV_i[((((t*2*n_a+w)*2*n_a+x)*2*n_a+v)*2*n_a+u)*2*n_a+x]
                                                +RF_P3_VV_i[((((u*2*n_a+w)*2*n_a+x)*2*n_a+v)*2*n_a+t)*2*n_a+x];
    for(int t=0; t<2*n_a; t++)
    for(int v=0; v<2*n_a; v++)
    for(int u=0; u<2*n_a; u++)
    for(int w=0; w<2*n_a; w++)
        RF_PV_r[((u*2*n_a+t)*2*n_a+v)*2*n_a+w]+=-RF_P3_VH_r[((t*2*n_a+v)*2*n_a+w)*2*n_a+u]
                                                -RF_P3_HV_r[((v*2*n_a+w)*2*n_a+u)*2*n_a+t]
                                                +RF_P3_VH_r[((u*2*n_a+v)*2*n_a+w)*2*n_a+t]
                                                +RF_P3_HV_r[((v*2*n_a+w)*2*n_a+t)*2*n_a+u]
                                                +RF_P3_VH_r[((t*2*n_a+w)*2*n_a+v)*2*n_a+u]
                                                +RF_P3_HV_r[((w*2*n_a+v)*2*n_a+u)*2*n_a+t]
                                                -RF_P3_VH_r[((u*2*n_a+w)*2*n_a+v)*2*n_a+t]
                                                -RF_P3_HV_r[((w*2*n_a+v)*2*n_a+t)*2*n_a+u];
    for(int t=0; t<2*n_a; t++)
    for(int v=0; v<2*n_a; v++)
    for(int u=0; u<2*n_a; u++)
    for(int w=0; w<2*n_a; w++)
        RF_PV_i[((u*2*n_a+t)*2*n_a+v)*2*n_a+w]+=-RF_P3_VH_i[((t*2*n_a+v)*2*n_a+w)*2*n_a+u]
                                                -RF_P3_HV_i[((v*2*n_a+w)*2*n_a+u)*2*n_a+t]
                                                +RF_P3_VH_i[((u*2*n_a+v)*2*n_a+w)*2*n_a+t]
                                                +RF_P3_HV_i[((v*2*n_a+w)*2*n_a+t)*2*n_a+u]
                                                +RF_P3_VH_i[((t*2*n_a+w)*2*n_a+v)*2*n_a+u]
                                                +RF_P3_HV_i[((w*2*n_a+v)*2*n_a+u)*2*n_a+t]
                                                -RF_P3_VH_i[((u*2*n_a+w)*2*n_a+v)*2*n_a+t]
                                                -RF_P3_HV_i[((w*2*n_a+v)*2*n_a+t)*2*n_a+u];

    
//     set_zero_matr(RF_PH_r, 2*n_a*2*n_a*N_fit);
    for(int t=0; t<2*n_a; t++)
    for(int u=0; u<2*n_a; u++)
    for(int v=0; v<2*n_a; v++)
    for(int y=0; y<2*n_a; y++)
        RF_PH_r[t*2*n_a+u]+=-RF_P3_VV_r[((((t*2*n_a+y)*2*n_a+v)*2*n_a+u)*2*n_a+v)*2*n_a+y]
                            +RF_P3_VV_r[((((t*2*n_a+v)*2*n_a+y)*2*n_a+u)*2*n_a+v)*2*n_a+y];
    for(int t=0; t<2*n_a; t++)
    for(int u=0; u<2*n_a; u++)
    for(int v=0; v<2*n_a; v++)
    for(int y=0; y<2*n_a; y++)
        RF_PH_i[t*2*n_a+u]+=-RF_P3_VV_i[((((t*2*n_a+y)*2*n_a+v)*2*n_a+u)*2*n_a+v)*2*n_a+y]
                            +RF_P3_VV_i[((((t*2*n_a+v)*2*n_a+y)*2*n_a+u)*2*n_a+v)*2*n_a+y];
    for(int t=0; t<2*n_a; t++)
    for(int u=0; u<2*n_a; u++)
    for(int y=0; y<2*n_a; y++)
        RF_PH_r[t*2*n_a+u]+=-RF_P3_VH_r[((t*2*n_a+y)*2*n_a+u)*2*n_a+y]
                            +RF_P3_VH_r[((t*2*n_a+u)*2*n_a+y)*2*n_a+y]
                            -RF_P3_HV_r[((y*2*n_a+u)*2*n_a+y)*2*n_a+t]
                            +RF_P3_HV_r[((y*2*n_a+u)*2*n_a+t)*2*n_a+y];
    for(int t=0; t<2*n_a; t++)
    for(int u=0; u<2*n_a; u++)
    for(int y=0; y<2*n_a; y++)
        RF_PH_i[t*2*n_a+u]+=-RF_P3_VH_i[((t*2*n_a+y)*2*n_a+u)*2*n_a+y]
                            +RF_P3_VH_i[((t*2*n_a+u)*2*n_a+y)*2*n_a+y]
                            -RF_P3_HV_i[((y*2*n_a+u)*2*n_a+y)*2*n_a+t]
                            +RF_P3_HV_i[((y*2*n_a+u)*2*n_a+t)*2*n_a+y];
    for(int t=0; t<2*n_a; t++)
    for(int u=0; u<2*n_a; u++)
        RF_PH_r[t*2*n_a+u]+=-RF_P3_HH_r[u*2*n_a+t];
    for(int t=0; t<2*n_a; t++)
    for(int u=0; u<2*n_a; u++)
        RF_PH_i[t*2*n_a+u]+=-RF_P3_HH_i[u*2*n_a+t];
    
//     set_zero_matr(RF_PS, N_fit);
    for(int t=0; t<2*n_a; t++)
        RF_PS+=RF_P3_HH_r[t*2*n_a+t];
    

    delete[] CAAA_r  ;
    delete[] CAAA_i  ;


    delete[] RF_P3_VV_r;
    delete[] RF_P3_VH_r;
    delete[] RF_P3_HV_r;
    delete[] RF_P3_HH_r;
    delete[] RF_P3_VV_i;
    delete[] RF_P3_VH_i;
    delete[] RF_P3_HV_i;
    delete[] RF_P3_HH_i;

    
    return 0;
}

int PT_tensors_rel::calc_EE_2_CCVV(){
    
    double * R;
    R = new  double[num_threads];
    set_zero_matr(R,num_threads);
    double * V[  num_threads];
    double *JK[2*num_threads];
    for(int i=0;i<  num_threads;i++) V[i] = new double[2*n_c];
    for(int i=0;i<2*num_threads;i++)JK[i] = new double[4*n_c*n_c];
    
    for(int a=0; a<2*n_v; a++){//fprintf(stderr,"CCVV a=%d\r",a);
#pragma omp parallel for
    for(int b=a; b<2*n_v; b++){//fprintf(stderr,"CCVV a,b=%d,%d\r",a,b);
        
        int th_id = omp_get_thread_num();
        double Cr;
        double Ci;
        int iajb=0;
        double dE;
        double Ep;
        double Jr;
        double Kr;
        double Ji;
        double Ki;
        //<a*b*|ij>=(a*i|b*j)=\sum[(a*i|K)(b*j|K)]
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        2*n_c,2*n_c,RI->aux_n_ao,1.0,
                        RI->VC_RI_M_r+a*2*n_c*RI->aux_n_ao,RI->aux_n_ao,
                        RI->VC_RI_M_r+b*2*n_c*RI->aux_n_ao,RI->aux_n_ao,0.0,
                        JK[2*th_id],2*n_c);
        
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        2*n_c,2*n_c,RI->aux_n_ao,-1.0,
                        RI->VC_RI_M_i+a*2*n_c*RI->aux_n_ao,RI->aux_n_ao,
                        RI->VC_RI_M_i+b*2*n_c*RI->aux_n_ao,RI->aux_n_ao,1.0,
                        JK[2*th_id],2*n_c);
        
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        2*n_c,2*n_c,RI->aux_n_ao,1.0,
                        RI->VC_RI_M_i+a*2*n_c*RI->aux_n_ao,RI->aux_n_ao,
                        RI->VC_RI_M_r+b*2*n_c*RI->aux_n_ao,RI->aux_n_ao,0.0,
                        JK[2*th_id+1],2*n_c);
        
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        2*n_c,2*n_c,RI->aux_n_ao,1.0,
                        RI->VC_RI_M_r+a*2*n_c*RI->aux_n_ao,RI->aux_n_ao,
                        RI->VC_RI_M_i+b*2*n_c*RI->aux_n_ao,RI->aux_n_ao,1.0,
                        JK[2*th_id+1],2*n_c);
        
        
        for(int i=0; i<2*n_c; i++){//fprintf(stderr,"i,a=%d,%d\r",i,a);
        for(int j=i; j<2*n_c; j++){
            Jr=JK[2*th_id  ][i*2*n_c+j];
            Kr=JK[2*th_id  ][j*2*n_c+i];
            Ji=JK[2*th_id+1][i*2*n_c+j];
            Ki=JK[2*th_id+1][j*2*n_c+i];
            //<a*b*||ij>=<a*b*|ij>-<a*b*|ji>
            Cr=Jr-Kr;
            Ci=Ji-Ki;
            //<i*j*||ab><a*b*||ij>=|<a*b*||ij>|^2
            V[th_id][j]=(Cr*Cr+Ci*Ci);///complex conjugate????
            
        }
        
        Ep=e_c[i]-e_v[a]-e_v[b];
        for(int j=i; j<2*n_c; j++){
            dE=e_c[j]+Ep;
            R[th_id]+=V[th_id][j]*dE/(dE*dE+edshift);           
        }
       
        }
    }
    }
    
    for(int i=0;i<num_threads;i++){
        RF_PS  +=R[i];
    }
    
    delete[] R;
    for(int i=0;i<num_threads;i++) delete[] V [i];
    for(int i=0;i<num_threads;i++) delete[] JK[i];
    
    
    return 0;
}

int PT_tensors_rel::calc_EE_2_CAVV(){
    double **PH_th= new double * [2*num_threads];
    for(int i=0;i<2*num_threads;i++)PH_th[i]=new double[2*n_a*2*n_a];
    for(int i=0;i<2*num_threads;i++)set_zero_matr(PH_th[i],2*n_a*2*n_a);

    #pragma omp parallel
    {
        int nt = omp_get_thread_num();
        
        double *Kr;
        double *Jr;
        double *Ki;
        double *Ji;
        
        Kr= new double[2*n_a*2*n_c];
        Jr= new double[2*n_a*2*n_c];
        Ki= new double[2*n_a*2*n_c];
        Ji= new double[2*n_a*2*n_c];
        
        
        
        double dE;
        double V2r;
        double V2i;
        double Tr;
        double Ti;
        double Vr;
        double Vi;
        
//         fprintf(out_stream,"CAVV is not parallel\n");
        
        for(int a=nt  ; a<2*n_v; a+=num_threads)
        for(int b=a  ; b<2*n_v; b++){//fprintf(stderr,"CAVVa,b=%d,%d\r",a,b);
            //<a*b*|it>=(a*i|b*j)=\sum[(a*i|K)(b*t|K)]
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            2*n_c,2*n_a,RI->aux_n_ao,1.0,
                            RI->VC_RI_M_r+a*2*n_c*RI->aux_n_ao,RI->aux_n_ao,
                            RI->VA_RI_M_r+b*2*n_a*RI->aux_n_ao,RI->aux_n_ao,0.0,
                            Jr,2*n_a);
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            2*n_c,2*n_a,RI->aux_n_ao,-1.0,
                            RI->VC_RI_M_i+a*2*n_c*RI->aux_n_ao,RI->aux_n_ao,
                            RI->VA_RI_M_i+b*2*n_a*RI->aux_n_ao,RI->aux_n_ao,1.0,
                            Jr,2*n_a);
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            2*n_c,2*n_a,RI->aux_n_ao,1.0,
                            RI->VC_RI_M_r+a*2*n_c*RI->aux_n_ao,RI->aux_n_ao,
                            RI->VA_RI_M_i+b*2*n_a*RI->aux_n_ao,RI->aux_n_ao,0.0,
                            Ji,2*n_a);
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            2*n_c,2*n_a,RI->aux_n_ao,1.0,
                            RI->VC_RI_M_i+a*2*n_c*RI->aux_n_ao,RI->aux_n_ao,
                            RI->VA_RI_M_r+b*2*n_a*RI->aux_n_ao,RI->aux_n_ao,1.0,
                            Ji,2*n_a);
            //<b*a*|it>=(b*i|a*j)=\sum[(b*i|K)(a*t|K)]
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            2*n_c,2*n_a,RI->aux_n_ao,1.0,
                            RI->VC_RI_M_r+b*2*n_c*RI->aux_n_ao,RI->aux_n_ao,
                            RI->VA_RI_M_r+a*2*n_a*RI->aux_n_ao,RI->aux_n_ao,0.0,
                            Kr,2*n_a);
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            2*n_c,2*n_a,RI->aux_n_ao,-1.0,
                            RI->VC_RI_M_i+b*2*n_c*RI->aux_n_ao,RI->aux_n_ao,
                            RI->VA_RI_M_i+a*2*n_a*RI->aux_n_ao,RI->aux_n_ao,1.0,
                            Kr,2*n_a);
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            2*n_c,2*n_a,RI->aux_n_ao,1.0,
                            RI->VC_RI_M_i+b*2*n_c*RI->aux_n_ao,RI->aux_n_ao,
                            RI->VA_RI_M_r+a*2*n_a*RI->aux_n_ao,RI->aux_n_ao,0.0,
                            Ki,2*n_a);
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            2*n_c,2*n_a,RI->aux_n_ao,1.0,
                            RI->VC_RI_M_r+b*2*n_c*RI->aux_n_ao,RI->aux_n_ao,
                            RI->VA_RI_M_i+a*2*n_a*RI->aux_n_ao,RI->aux_n_ao,1.0,
                            Ki,2*n_a);
            
            for(int i=0;i<2*n_c;i++)
            for(int t=0;t<2*n_a;t++)
            for(int v=0;v<2*n_a;v++){
                Tr=Jr[i*2*n_a+t]-Kr[i*2*n_a+t];
                Ti=Ji[i*2*n_a+t]-Ki[i*2*n_a+t];
                Vr=Jr[i*2*n_a+v]-Kr[i*2*n_a+v];
                Vi=Ji[i*2*n_a+v]-Ki[i*2*n_a+v];
                //<i*v*||ab><a*b*||it>=<a*b*||iv>* <a*b*||ij>
                V2r =Tr*Vr+Ti*Vi;
                V2i =Ti*Vr-Tr*Vi;
//                 V2+=J[i*2*n_a+t]*J[i*2*n_a+v];
//                 if(a!=b)
//                     V2+=K[i*2*n_a+t]*K[i*2*n_a+v];
                dE=e_c[i]+e_a[t]-e_v[a]-e_v[b];
                PH_th[2*nt  ][t*2*n_a+v]+=V2r*dE/(dE*dE+edshift);
                PH_th[2*nt+1][t*2*n_a+v]+=V2i*dE/(dE*dE+edshift);
                
            }
        }
        
        delete[] Kr;
        delete[] Ki;
        delete[] Jr;
        delete[] Ji;
    }
    
    for(long j=0; j<num_threads;j++)
    #pragma omp parallel for
        for(long i=0; i<2*n_a*2*n_a;i++){
            RF_PH_r[i]+=PH_th[2*j  ][i];
            RF_PH_i[i]+=PH_th[2*j+1][i];
        }
        
    for(int i=0;i<2*num_threads;i++)delete[] PH_th[i];
    delete[] PH_th;
    
    
    return 0;
    
}

int PT_tensors_rel::calc_EE_2_AAVV(){
    
//     num_threads=1;
    double **Vr_th= new double * [num_threads];
    for(int i=0;i<num_threads;i++)Vr_th[i]=new double[2*n_a*2*n_a*2*n_a*2*n_a];
    for(int i=0;i<num_threads;i++)set_zero_matr(Vr_th[i],2*n_a*2*n_a*2*n_a*2*n_a);

    double **Vi_th= new double * [num_threads];
    for(int i=0;i<num_threads;i++)Vi_th[i]=new double[2*n_a*2*n_a*2*n_a*2*n_a];
    for(int i=0;i<num_threads;i++)set_zero_matr(Vi_th[i],2*n_a*2*n_a*2*n_a*2*n_a);

    
    #pragma omp parallel
    {
        int nt =omp_get_thread_num();
        
        double TUr, TUi,VWr, VWi, dE;
        
        double * Jr;
        double * Ji;
        
        Jr = new double[2*n_a*2*n_a];
        Ji = new double[2*n_a*2*n_a];
        
        for(int a=nt; a<2*n_v; a+=num_threads)
        for(int b=a ; b<2*n_v; b++){//fprintf(stderr,"AAVV a,b=%3d,%3d  (thread %2d)    \r",a,b,nt);
            //<a*b*|tu>
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            2*n_a,2*n_a,RI->aux_n_ao,1.0,
                            RI->VA_RI_M_r+a*2*n_a*RI->aux_n_ao,RI->aux_n_ao,
                            RI->VA_RI_M_r+b*2*n_a*RI->aux_n_ao,RI->aux_n_ao,0.0,
                            Jr,2*n_a);
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            2*n_a,2*n_a,RI->aux_n_ao,-1.0,
                            RI->VA_RI_M_i+a*2*n_a*RI->aux_n_ao,RI->aux_n_ao,
                            RI->VA_RI_M_i+b*2*n_a*RI->aux_n_ao,RI->aux_n_ao,1.0,
                            Jr,2*n_a);
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            2*n_a,2*n_a,RI->aux_n_ao,1.0,
                            RI->VA_RI_M_r+a*2*n_a*RI->aux_n_ao,RI->aux_n_ao,
                            RI->VA_RI_M_i+b*2*n_a*RI->aux_n_ao,RI->aux_n_ao,0.0,
                            Ji,2*n_a);
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            2*n_a,2*n_a,RI->aux_n_ao,1.0,
                            RI->VA_RI_M_i+a*2*n_a*RI->aux_n_ao,RI->aux_n_ao,
                            RI->VA_RI_M_r+b*2*n_a*RI->aux_n_ao,RI->aux_n_ao,1.0,
                            Ji,2*n_a);
            
            
            
            
            for(int t=0  ;t<2*n_a;t++)
            for(int u=0  ;u<2*n_a;u++)
            for(int v=0  ;v<2*n_a;v++)
            for(int w=0  ;w<2*n_a;w++){
                
                TUr=Jr[t*2*n_a+u]-Jr[u*2*n_a+t];
                TUi=Ji[t*2*n_a+u]-Ji[u*2*n_a+t];
                VWr=Jr[v*2*n_a+w]-Jr[w*2*n_a+v];
                VWi=Ji[v*2*n_a+w]-Ji[w*2*n_a+v];
                
                /*JK=(Ja[t*2*n_a+u]-Ja[u*2*n_a+t])*
                   (Ja[v*2*n_a+w]-Ja[w*2*n_a+v]);
                AB         =Ja[t*2*n_a+u]*Ja[v*2*n_a+w];
                if(a!=b)AB+=Ja[u*2*n_a+t]*Ja[w*2*n_a+v];
                */   
                dE=e_a[t]+e_a[u]-e_v[a]-e_v[b];
                Vr_th[nt][((u*2*n_a+t)*2*n_a+v)*2*n_a+w]+=(TUr*VWr+TUi*VWi)*dE/(dE*dE+edshift);///inverse t-u in res_fit_calc_2body_AA (?)
                Vi_th[nt][((t*2*n_a+u)*2*n_a+v)*2*n_a+w]+=(TUi*VWr-TUr*VWi)*dE/(dE*dE+edshift);//check sign!!!!!
                
                
            }
                
            
        }
    
       delete[] Jr;
       delete[] Ji;
    
    
    }
    
//     #pragma omp parallel for
    for(long j=0; j<num_threads;j++)
        for(long i=0; i<2*n_a*2*n_a*2*n_a*2*n_a;i++)
            RF_PV_r[i]+=Vr_th[j][i];
//     #pragma omp parallel for
    for(long j=0; j<num_threads;j++)
        for(long i=0; i<2*n_a*2*n_a*2*n_a*2*n_a;i++)
            RF_PV_i[i]+=Vi_th[j][i];
        
    for(int i=0;i<num_threads;i++)delete[] Vr_th[i];
    delete[] Vr_th;
    
    for(int i=0;i<num_threads;i++)delete[] Vi_th[i];
    delete[] Vi_th;
    return 0;
}

int PT_tensors_rel::calc_EE_2_CCAV(){
    

    
    
    double **RF_PH_th= new double * [2*num_threads];
    for(int i=0;i<2*num_threads;i++)RF_PH_th[i]   = new double[2*n_a*2*n_a];
    for(int i=0;i<2*num_threads;i++)set_zero_matr(RF_PH_th[i], 2*n_a*2*n_a);
    
    #pragma omp parallel
    {
        int nt =omp_get_thread_num();
        
        double *Kr;
        double *Ki;
        double *Jr;
        double *Ji;
        Kr= new double[2*n_a*2*n_v];
        Ki= new double[2*n_a*2*n_v];
        Jr= new double[2*n_a*2*n_v];
        Ji= new double[2*n_a*2*n_v];
        double Tr;
        double Ti;
        double Vr;
        double Vi;
        
        double dE;
        double V2r;
        double V2i;
        
        for(int i=nt; i<2*n_c; i+=num_threads)
        for(int j=i; j<2*n_c; j++){//fprintf(stderr," CCAV i,j=%d,%d\r",i,j);
            for(int a=0; a<2*n_v; a++)
            for(int t=0; t<2*n_a; t++){
                
                Jr[a*2*n_a+t]= cblas_ddot(RI->aux_n_ao, RI->VC_RI_M_r+(a*2*n_c+i)*RI->aux_n_ao, 1, RI->CA_RI_M_r+(j*2*n_a+t)*RI->aux_n_ao, 1)
                              +cblas_ddot(RI->aux_n_ao, RI->VC_RI_M_i+(a*2*n_c+i)*RI->aux_n_ao, 1, RI->CA_RI_M_i+(j*2*n_a+t)*RI->aux_n_ao, 1);
                Ji[a*2*n_a+t]=-cblas_ddot(RI->aux_n_ao, RI->VC_RI_M_r+(a*2*n_c+i)*RI->aux_n_ao, 1, RI->CA_RI_M_i+(j*2*n_a+t)*RI->aux_n_ao, 1)
                              +cblas_ddot(RI->aux_n_ao, RI->VC_RI_M_i+(a*2*n_c+i)*RI->aux_n_ao, 1, RI->CA_RI_M_r+(j*2*n_a+t)*RI->aux_n_ao, 1);
                Kr[a*2*n_a+t]= cblas_ddot(RI->aux_n_ao, RI->VC_RI_M_r+(a*2*n_c+j)*RI->aux_n_ao, 1, RI->CA_RI_M_r+(i*2*n_a+t)*RI->aux_n_ao, 1)
                              +cblas_ddot(RI->aux_n_ao, RI->VC_RI_M_i+(a*2*n_c+j)*RI->aux_n_ao, 1, RI->CA_RI_M_i+(i*2*n_a+t)*RI->aux_n_ao, 1);
                Ki[a*2*n_a+t]=-cblas_ddot(RI->aux_n_ao, RI->VC_RI_M_r+(a*2*n_c+j)*RI->aux_n_ao, 1, RI->CA_RI_M_i+(i*2*n_a+t)*RI->aux_n_ao, 1)
                              +cblas_ddot(RI->aux_n_ao, RI->VC_RI_M_i+(a*2*n_c+j)*RI->aux_n_ao, 1, RI->CA_RI_M_r+(i*2*n_a+t)*RI->aux_n_ao, 1);
            }
            
            for(int a=0;a<2*n_v;a++)
            for(int t=0;t<2*n_a;t++)
            for(int v=0;v<2*n_a;v++){
                Tr=Jr[a*2*n_a+t]-Kr[a*2*n_a+t];
                Ti=Ji[a*2*n_a+t]-Ki[a*2*n_a+t];
                Vr=Jr[a*2*n_a+v]-Kr[a*2*n_a+v];
                Vi=Ji[a*2*n_a+v]-Ki[a*2*n_a+v];
                //<i*v*||ab><a*b*||it>=<a*b*||iv>* <a*b*||ij>
                V2r =Tr*Vr+Ti*Vi;
                V2i =Ti*Vr-Tr*Vi;
                
//                 V2+=J[a*2*n_a+t]*J[a*2*n_a+v];    //i<j is filled by J
//                 if(i!=j)                      //i=j is filled only once J[i,i]=K[i,i]
//                     V2+=K[a*2*n_a+t]*K[a*2*n_a+v];//i>j is filled by K
                dE=e_c[i]+e_c[j]-e_a[t]-e_v[a];
                RF_PH_th[2*nt  ][v*2*n_a+t]-=V2r*dE/(dE*dE+edshift);//sign changed!!!!
                RF_PH_th[2*nt+1][v*2*n_a+t]-=V2i*dE/(dE*dE+edshift);//sign changed!!!!check im part!!!!
                
            }
        }
        
        delete[] Kr;
        delete[] Jr;
        delete[] Ki;
        delete[] Ji;
    }
    
    for(long j=1; j<num_threads;j++)
//     #pragma omp parallel for
        for(long i=0; i<2*n_a*2*n_a;i++){
            RF_PH_th[0][i]+=RF_PH_th[2*j  ][i];
            RF_PH_th[1][i]+=RF_PH_th[2*j+1][i];
        }
            
    
//     #pragma omp parallel for
    for(int t=0; t<2*n_a; t++)
        RF_PS-=RF_PH_th[0][t*2*n_a+t];
    
    
    #pragma omp parallel for // parallelism is bad but i don't care
    for(int w=0; w<2*n_a; w++)
    for(int t=0; t<2*n_a; t++){
        RF_PH_r[w*2*n_a+t]+=RF_PH_th[0][w*2*n_a+t];
        RF_PH_i[w*2*n_a+t]+=RF_PH_th[1][w*2*n_a+t];
    }
    
    for(int i=0;i<2*num_threads;i++)delete[] RF_PH_th[i];
    delete[] RF_PH_th;
    
    
    return 0;
    
}

int PT_tensors_rel::calc_EE_2_CCAA(){
    
//     num_threads=1;
    
    double ** Vr_th = new double * [num_threads];
    double ** Vi_th = new double * [num_threads];
    
    for(int i=0;i<num_threads;i++)Vr_th[i] = new double[2*n_a*2*n_a*2*n_a*2*n_a];
    for(int i=0;i<num_threads;i++)Vi_th[i] = new double[2*n_a*2*n_a*2*n_a*2*n_a];
    
    for(int i=0;i<num_threads;i++)set_zero_matr(Vr_th[i], 2*n_a*2*n_a*2*n_a*2*n_a);
    for(int i=0;i<num_threads;i++)set_zero_matr(Vi_th[i], 2*n_a*2*n_a*2*n_a*2*n_a);
    
    
    #pragma omp parallel
    {
        int nt =omp_get_thread_num();
    
        double dE;
//         double Ed;
        
        double TUr,TUi;
        double VWr,VWi;
        
        double * Jr;
        double * Ji;
        Jr = new double[2*n_a*2*n_a];
        Ji = new double[2*n_a*2*n_a];
        
        
        for(int i=nt; i<2*n_c; i+=num_threads)
        for(int j=i ; j<2*n_c; j++){//fprintf(stderr,"CCAA i,j=%3d,%3d  (thread %2d)    \r",i,j,nt);
            //<t*u*|ij>=<i*j*|tu>*
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            2*n_a,2*n_a,RI->aux_n_ao,1.0,
                            RI->CA_RI_M_r+i*2*n_a*RI->aux_n_ao,RI->aux_n_ao,
                            RI->CA_RI_M_r+j*2*n_a*RI->aux_n_ao,RI->aux_n_ao,0.0,
                            Jr,2*n_a);
            
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            2*n_a,2*n_a,RI->aux_n_ao,-1.0,
                            RI->CA_RI_M_i+i*2*n_a*RI->aux_n_ao,RI->aux_n_ao,
                            RI->CA_RI_M_i+j*2*n_a*RI->aux_n_ao,RI->aux_n_ao,1.0,
                            Jr,2*n_a);
            
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            2*n_a,2*n_a,RI->aux_n_ao,-1.0,
                            RI->CA_RI_M_i+i*2*n_a*RI->aux_n_ao,RI->aux_n_ao,
                            RI->CA_RI_M_r+j*2*n_a*RI->aux_n_ao,RI->aux_n_ao,0.0,
                            Ji,2*n_a);
            
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            2*n_a,2*n_a,RI->aux_n_ao,-1.0,
                            RI->CA_RI_M_r+i*2*n_a*RI->aux_n_ao,RI->aux_n_ao,
                            RI->CA_RI_M_i+j*2*n_a*RI->aux_n_ao,RI->aux_n_ao,1.0,
                            Ji,2*n_a);
            
            
            
            for(int t=0  ;t<2*n_a;t++)
            for(int u=0  ;u<2*n_a;u++)
            for(int v=0  ;v<2*n_a;v++)
            for(int w=0  ;w<2*n_a;w++){
                
                TUr=Jr[t*2*n_a+u]-Jr[u*2*n_a+t];
                TUi=Ji[t*2*n_a+u]-Ji[u*2*n_a+t];
                VWr=Jr[v*2*n_a+w]-Jr[w*2*n_a+v];
                VWi=Ji[v*2*n_a+w]-Ji[w*2*n_a+v];
                
                dE=-e_a[v]-e_a[w]+e_c[i]+e_c[j];
                
                Vr_th[nt][((u*2*n_a+t)*2*n_a+v)*2*n_a+w]+=(TUr*VWr+TUi*VWi)*dE/(dE*dE+edshift);///inverse t-u in res_fit_calc_2body_AA (?)
                Vi_th[nt][((t*2*n_a+u)*2*n_a+v)*2*n_a+w]+=(TUi*VWr-TUr*VWi)*dE/(dE*dE+edshift);//check sign
                
            }
        }
        
        delete[] Jr;
        delete[] Ji;
    }
    
    for(long j=1; j<num_threads;j++)
//     #pragma omp parallel for
        for(long i=0; i<2*n_a*2*n_a*2*n_a*2*n_a;i++)
            Vr_th[0][i]+=Vr_th[j][i];

    for(long j=1; j<num_threads;j++)
//     #pragma omp parallel for
        for(long i=0; i<2*n_a*2*n_a*2*n_a*2*n_a;i++)
            Vi_th[0][i]+=Vi_th[j][i];

    
    #pragma omp parallel for // parallelism is bad but i don't care
    for(long i=0;i<2*n_a*2*n_a*2*n_a*2*n_a;i++)RF_PV_r[i]+=Vr_th[0][i];
    
    #pragma omp parallel for // parallelism is bad but i don't care
    for(long i=0;i<2*n_a*2*n_a*2*n_a*2*n_a;i++)RF_PV_i[i]+=Vi_th[0][i];
    
    #pragma omp parallel for // parallelism is bad but i don't care
    for(int v=0; v<2*n_a; v++)
    for(int u=0; u<2*n_a; u++)
    for(int t=0; t<2*n_a; t++)
        RF_PH_r[v*2*n_a+u]+=Vr_th[0][((t*2*n_a+v)*2*n_a+t)*2*n_a+u]*0.25;
    
    #pragma omp parallel for // parallelism is bad but i don't care
    for(int v=0; v<2*n_a; v++)
    for(int u=0; u<2*n_a; u++)
    for(int t=0; t<2*n_a; t++)
        RF_PH_r[v*2*n_a+t]-=Vr_th[0][((u*2*n_a+v)*2*n_a+t)*2*n_a+u]*0.25;

    #pragma omp parallel for // parallelism is bad but i don't care
    for(int w=0; w<2*n_a; w++)
    for(int v=0; v<2*n_a; v++)
    for(int u=0; u<2*n_a; u++)
        RF_PH_r[w*2*n_a+u]-=Vr_th[0][((w*2*n_a+v)*2*n_a+v)*2*n_a+u]*0.25;

    #pragma omp parallel for // parallelism is bad but i don't care
    for(int w=0; w<2*n_a; w++)
    for(int u=0; u<2*n_a; u++)
    for(int t=0; t<2*n_a; t++)
        RF_PH_r[w*2*n_a+t]+=Vr_th[0][((w*2*n_a+u)*2*n_a+t)*2*n_a+u]*0.25;

    #pragma omp parallel for // parallelism is bad but i don't care
    for(int v=0; v<2*n_a; v++)
    for(int u=0; u<2*n_a; u++)
    for(int t=0; t<2*n_a; t++)
        RF_PH_i[v*2*n_a+u]+=Vi_th[0][((t*2*n_a+v)*2*n_a+t)*2*n_a+u]*0.25;
    
    #pragma omp parallel for // parallelism is bad but i don't care
    for(int v=0; v<2*n_a; v++)
    for(int u=0; u<2*n_a; u++)
    for(int t=0; t<2*n_a; t++)
        RF_PH_i[v*2*n_a+t]-=Vi_th[0][((u*2*n_a+v)*2*n_a+t)*2*n_a+u]*0.25;

    #pragma omp parallel for // parallelism is bad but i don't care
    for(int w=0; w<2*n_a; w++)
    for(int v=0; v<2*n_a; v++)
    for(int u=0; u<2*n_a; u++)
        RF_PH_i[w*2*n_a+u]-=Vi_th[0][((w*2*n_a+v)*2*n_a+v)*2*n_a+u]*0.25;

    #pragma omp parallel for // parallelism is bad but i don't care
    for(int w=0; w<2*n_a; w++)
    for(int u=0; u<2*n_a; u++)
    for(int t=0; t<2*n_a; t++)
        RF_PH_i[w*2*n_a+t]+=Vi_th[0][((w*2*n_a+u)*2*n_a+t)*2*n_a+u]*0.25;

//     #pragma omp parallel for // parallelism is bad but i don't care
//     for(int w=0; w<2*n_a; w++)
//     for(int u=0; u<2*n_a; u++)
//     for(int t=0; t<2*n_a; t++)
//         RF_PH_r[w*2*n_a+u]-=RF_PV_AB_th[0][((t*2*n_a+w)*2*n_a+t)*2*n_a+u];
// 


    for(int v=0; v<2*n_a; v++)
    for(int u=0; u<2*n_a; u++)
        RF_PS+=Vr_th[0][((u*2*n_a+v)*2*n_a+v)*2*n_a+u]*0.25;//A+B

    for(int v=0; v<2*n_a; v++)
    for(int w=0; w<2*n_a; w++)
        RF_PS-=Vr_th[0][((w*2*n_a+v)*2*n_a+w)*2*n_a+v]*0.25;//A+B
    
//     double im=0;
//     for(int v=0; v<2*n_a; v++)
//     for(int u=0; u<2*n_a; u++)
//         im+=Vi_th[0][((u*2*n_a+v)*2*n_a+v)*2*n_a+u]*0.5;//A+B
// 
//     for(int v=0; v<2*n_a; v++)
//     for(int w=0; w<2*n_a; w++)
//         im-=Vi_th[0][((w*2*n_a+v)*2*n_a+w)*2*n_a+v]*0.5;//A+B
//     
//     printf("im=%e\n",im);
//     
//     
//     for(int u=0; u<2*n_a; u++)
//     for(int t=0; t<2*n_a; t++)
//         RF_PS+=RF_PV_AB_th[0][((t*2*n_a+u)*2*n_a+t)*2*n_a+u];

    
    
    
    for(int i=0;i<num_threads;i++)delete[] Vr_th[i];
    for(int i=0;i<num_threads;i++)delete[] Vi_th[i];
    
    
    delete[] Vr_th;
    delete[] Vi_th;
    
    return 0;
}


PT_tensors_rel::~PT_tensors_rel(){
    
    if(RF_PH_r  != nullptr) delete[] RF_PH_r ;
    if(RF_PH_i  != nullptr) delete[] RF_PH_i ;
    if(RF_PV    != nullptr) delete[] RF_PV   ;
    if(RF_PV_r  != nullptr) delete[] RF_PV_r ;
    if(RF_PV_i  != nullptr) delete[] RF_PV_i ;
    if(RF_P3_r  != nullptr) delete[] RF_P3_r ;
    if(RF_P3_i  != nullptr) delete[] RF_P3_i ;
    if(e_IP     != nullptr) delete[] e_IP    ;
    if(e_EA     != nullptr) delete[] e_EA    ;
    if(e_IP_2   != nullptr) delete[] e_IP_2  ;
    if(e_EA_2   != nullptr) delete[] e_EA_2  ;
    if(e_IP_AB  != nullptr) delete[] e_IP_AB ;
    if(e_EA_AB  != nullptr) delete[] e_EA_AB ;
    if(IP_U     != nullptr) delete[] IP_U    ;
    if(EA_U     != nullptr) delete[] EA_U    ;
    if(IP_Um    != nullptr) delete[] IP_Um   ;
    if(EA_Um    != nullptr) delete[] EA_Um   ;
    if(IP_Up    != nullptr) delete[] IP_Up   ;
    if(EA_Up    != nullptr) delete[] EA_Up   ;
    if(IP_2_U   != nullptr) delete[] IP_2_U  ;
    if(EA_2_U   != nullptr) delete[] EA_2_U  ;
    if(IP_2_Um  != nullptr) delete[] IP_2_Um ;
    if(EA_2_Um  != nullptr) delete[] EA_2_Um ;
    if(IP_2_Up  != nullptr) delete[] IP_2_Up ;
    if(EA_2_Up  != nullptr) delete[] EA_2_Up ;
    if(IP_AB_U  != nullptr) delete[] IP_AB_U ;
    if(EA_AB_U  != nullptr) delete[] EA_AB_U ;
    if(IP_AB_Um != nullptr) delete[] IP_AB_Um;
    if(EA_AB_Um != nullptr) delete[] EA_AB_Um;
    if(IP_AB_Up != nullptr) delete[] IP_AB_Up;
    if(EA_AB_Up != nullptr) delete[] EA_AB_Up;
    
}


