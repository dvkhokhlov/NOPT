# include <stdio.h>
# include "blas_link.h"
# include <omp.h>


# include "CI.h"
# include "matr.h"
# include "timer.h"

# include "from_hash.h"
# include "aldet.h"
# include "RI.h"
# include "PT_tensors.h"
# include "libint_link.h"


# define max(a,b)  (((a)<(b))?(b):(a))
# define min(a,b)  (((a)>(b))?(b):(a))

extern int num_threads;


inline double ED_with_shift(double E, double edshift){
    
    return E/(E*E+edshift);
}


int PT_tensors::E2_calc_EE(){
    
    set_zero_matr(RF_P3_JK, n_a*n_a*n_a*n_a*n_a*n_a);
    set_zero_matr(RF_P3_AB, n_a*n_a*n_a*n_a*n_a*n_a);
    
    set_zero_matr(RF_PV_JK, n_a*n_a*n_a*n_a);
    set_zero_matr(RF_PV_AB, n_a*n_a*n_a*n_a);
    
    set_zero_matr(RF_PH, n_a*n_a);
    
//     set_zero_matr(RF_PS, N_fit);
    RF_PS=0;
    
    calc_EE_2_CCVV();printf_timer("CCVV res table");fprintf(out_stream,"\n");fflush(out_stream);//done
    calc_EE_2_CAVV();printf_timer("CAVV res table");fprintf(out_stream,"\n");fflush(out_stream);//done
    calc_EE_2_AAVV();printf_timer("AAVV res table");fprintf(out_stream,"\n");fflush(out_stream);//done1
    calc_EE_2_CCAV();printf_timer("CCAV res table");fprintf(out_stream,"\n");fflush(out_stream);//done
    calc_EE_2_CCAA();printf_timer("CCAA res table");fprintf(out_stream,"\n");fflush(out_stream);//done
    calc_EE_2_CV  ();printf_timer("CV   res table");fprintf(out_stream,"\n");fflush(out_stream);//done
    calc_EE_2_AV  ();printf_timer("AV   res table");fprintf(out_stream,"\n");fflush(out_stream);//done1
    calc_EE_2_CA  ();printf_timer("CA   res table");fprintf(out_stream,"\n");fflush(out_stream);//done
    
//     printf("const=%e\n\n", RF_PS);
//     printf("M1:\n");
//     PrintMatr(RF_PH,n_a,n_a,0);
//     printf("M2_JK:\n");
//     PrintMatr(RF_PV_JK,n_a*n_a,n_a*n_a,0);
//     printf("M2_AB:\n");
//     PrintMatr(RF_PV_AB,n_a*n_a,n_a*n_a,0);
    
    
//     exit(0);
    
    symmetrization(RF_PH, n_a);
    symmetrization(RF_PV_JK, n_a*n_a);
    symmetrization(RF_PV_AB, n_a*n_a);
    symmetrization(RF_P3_JK, n_a*n_a*n_a);
    symmetrization(RF_P3_AB, n_a*n_a*n_a);
    
    double * RF_P3_tmp = new double  [n_a*n_a*n_a*n_a*n_a*n_a];
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++)
    for(int v=0; v<n_a; v++)
    for(int w=0; w<n_a; w++)
    for(int x=0; x<n_a; x++)
    for(int y=0; y<n_a; y++){
        RF_P3_tmp[((((t*n_a+u)*n_a+v)*n_a+w)*n_a+x)*n_a+y]=
        RF_P3_AB    [((((t*n_a+u)*n_a+x)*n_a+v)*n_a+w)*n_a+y];
    }
    memcpy(RF_P3_AB,RF_P3_tmp,n_a*n_a*n_a*n_a*n_a*n_a*sizeof(double));
    
    delete[]RF_P3_tmp;
    
    return 0;
    
}

int PT_tensors::calc_EE_2_CV(){
    
    double dE;
    double V;
    
    
    double * RF_MM;    // (J-K) * (J-K)
    double * RF_JJ;    //  J(aa)* J(aa)
    double * RF_JM;    //  J(aa)*(J-K)
    double * RF_MJ;    // (J-K) * J(aa)
    double * RF_AB;    //  J(ab)* J(ab)
    
    RF_MM = new double  [n_a*n_a*n_a*n_a];
    RF_JJ = new double  [n_a*n_a*n_a*n_a];
    RF_JM = new double  [n_a*n_a*n_a*n_a];
    RF_MJ = new double  [n_a*n_a*n_a*n_a];
    RF_AB = new double  [n_a*n_a*n_a*n_a];
    set_zero_matr(RF_MM, n_a*n_a*n_a*n_a);
    set_zero_matr(RF_JJ, n_a*n_a*n_a*n_a);
    set_zero_matr(RF_JM, n_a*n_a*n_a*n_a);
    set_zero_matr(RF_MJ, n_a*n_a*n_a*n_a);
    set_zero_matr(RF_AB, n_a*n_a*n_a*n_a);
    
    double * RF_HM;
    double * RF_MH;
    double * RF_HJ;
    double * RF_JH;
    RF_HM = new double  [n_a*n_a];
    RF_MH = new double  [n_a*n_a];
    RF_HJ = new double  [n_a*n_a];
    RF_JH = new double  [n_a*n_a];
    set_zero_matr(RF_HM, n_a*n_a);
    set_zero_matr(RF_MH, n_a*n_a);
    set_zero_matr(RF_HJ, n_a*n_a);
    set_zero_matr(RF_JH, n_a*n_a);
    
    double RF_H=0;
//     RF_H = new double  [n_a];
//     set_zero_matr(RF_H, n_a);
    
    double * VCAA;
    VCAA = new double[n_c*n_v*n_a*n_a];
    double * CAVA;
    CAVA = new double[n_a*n_c*n_a*n_v];
    
    if(RI){
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        n_c*n_v,n_a*n_a,integrals->aux_n_ao,1.0,
                        integrals->VC_RI_M,integrals->aux_n_ao,
                        integrals->AA_RI_M,integrals->aux_n_ao,0.0,
                        VCAA,n_a*n_a);
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        n_a*n_c,n_a*n_v,integrals->aux_n_ao,1.0,
                        integrals->CA_RI_M,integrals->aux_n_ao,
                        integrals->VA_RI_M,integrals->aux_n_ao,0.0,
                        CAVA,n_a*n_v);
    }
    else{
        int n_ao = integrals->n_ao;
        int n_o  = n_c+n_a;//occupied
        double * CVEC = new double[(n_c+n_a)*n_ao];
        for(int i=0; i<n_ao ;i++){
            for(int j=0; j<n_c;j++) CVEC[i*(n_c+n_a)+j    ]=integrals->COR_CVEC[i*n_c+j];
            for(int j=0; j<n_a;j++) CVEC[i*(n_c+n_a)+j+n_c]=integrals->ACT_CVEC[i*n_a+j];
        }
        double * IIIE = new double[n_v*n_o*n_o*n_o];
        calc_2el_MO_INTS_XXXY(integrals->s, integrals->n_ao, IIIE, CVEC, integrals->VAC_CVEC, n_o, n_v);
        delete[] CVEC;
        
        for(int a=0;a<n_v;a++) //// ? ïî÷åìó â òàêîì ïîðÿäêå ñóììèðîâàíèå? âðîäå ÷èñëî àêòèâíûõ ñàìîå ìàëîå? è ýòî ñóììèðîâàíèå äîëæíî áûòü âíåøíåå? 
        for(int i=0;i<n_c;i++)
        for(int t=0;t<n_a;t++)
        for(int u=0;u<n_a;u++)
            VCAA[((a*n_c+i)*n_a+u)*n_a+t]=IIIE[(((t+n_c)*n_o+u+n_c)*n_o+i)*n_v+a];   //// ñðàâíèâàþ ñ cdas, ÷òî ýòî òàêîå ?? ÍÅÒ - <ia|tu>? èëè íàîáîðîò, double empty act act
                                                                                       //// òî÷íî, ýòî íàâåðíã (ia|tu)
        for(int i=0;i<n_c;i++)
        for(int u=0;u<n_a;u++)
        for(int a=0;a<n_v;a++)
        for(int t=0;t<n_a;t++)
            CAVA[((i*n_a+u)*n_v+a)*n_a+t]=IIIE[((i*n_o+u+n_c)*n_o+t+n_c)*n_v+a];  //// IIIE? inact inact inact act?
        
        
        
        delete[] IIIE;
    }
        
        
//     double * VCAA_J;
//     VCAA_J = new double[n_c*n_v];
//     double * VCAA_K;
//     VCAA_K = new double[n_c*n_v];

    for(int i=0; i<n_c; i++)
    for(int u=0; u<n_a; u++)
    for(int a=0; a<n_v; a++)
    for(int t=0; t<n_a; t++){
        dE=e_c[i]+e_IP[t]-e_v[a]-e_EA[u];//is taken negative
        if(t==u)dE=e_c[i]-e_v[a];
        dE=dE/(dE*dE+edshift);
        
        for(int v=0; v<n_a; v++)
        for(int w=0; w<n_a; w++){
            //<i(a)t(b)|u(a)a(b)><v(a)a(b)|i(a)w(b)>=(iu|at)(iv|aw)       //// òóò âðîäå îïå÷àòêà <it|ua>=(iu|ta)   ??? see res_fit.cpp
            V=CAVA[((i*n_a+u)*n_v+a)*n_a+t]*
              CAVA[((i*n_a+v)*n_v+a)*n_a+w];                 //// à òóò îòêóäà òàêèå êàóíòåðû? ÷òî ïðîïóñêàåòñÿ òàêèì îáðàçîì? ÿ äóìàë, ÷òî èíäåêñ áëîêà ìíîæèòñÿ íà ÷èñëî ÅÃÎ ýëåìåíòîâ
              
            RF_AB[((t*n_a+u)*n_a+v)*n_a+w]+=V*dE;
            
        }
    }
    
    for(int i=0; i<n_c; i++)
    for(int a=0; a<n_v; a++)
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++){
        dE=e_c[i]+e_IP[t]-e_v[a]-e_EA[u];//is taken negative
        if(t==u)dE=e_c[i]-e_v[a];
        dE=dE/(dE*dE+edshift);
        
        for(int v=0; v<n_a; v++)
        for(int w=0; w<n_a; w++){
            //(<it|au>-<iu|at>)(<av|iw>-<av|wi>)
            V=(VCAA[((a*n_c+i)*n_a+u)*n_a+t]-  //<it|au>=(ai|tu)
               CAVA[((i*n_a+u)*n_v+a)*n_a+t])* //<it|ua>=(iu|at)
              (VCAA[((a*n_c+i)*n_a+w)*n_a+v]-  //<av|iw>=(ai|wv)
               CAVA[((i*n_a+v)*n_v+a)*n_a+w]); //<av|wi>=(iv|aw)
            
            RF_MM[((t*n_a+u)*n_a+v)*n_a+w]+=V*dE;
        }
    }
    

    for(int i=0; i<n_c; i++)
    for(int a=0; a<n_v; a++)
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++){
        dE=e_c[i]+e_IP[t]-e_v[a]-e_EA[u];//is taken negative
        if(t==u)dE=e_c[i]-e_v[a];
        dE=dE/(dE*dE+edshift);
        
        for(int v=0; v<n_a; v++)
        for(int w=0; w<n_a; w++){
            //(<it|au>-<iu|at>)<av|iw>
            V=(VCAA[((a*n_c+i)*n_a+u)*n_a+t]-  //<it|au>=(ai|tu)         //// îïÿòü îùóùåíèå, ÷òî ñ ïîðÿäêîì ÷òî-òî íå òàê (ia|tu)
               CAVA[((i*n_a+u)*n_v+a)*n_a+t])* //<it|ua>=(iu|at)
               VCAA[((a*n_c+i)*n_a+w)*n_a+v];  //<av|iw>=(ai|wv)
            
            RF_MJ[((t*n_a+u)*n_a+v)*n_a+w]+=V*dE;
        }
    }
    
    
    for(int i=0; i<n_c; i++)
    for(int a=0; a<n_v; a++)
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++){
        dE=e_c[i]+e_IP[t]-e_v[a]-e_EA[u];
        if(t==u)dE=e_c[i]-e_v[a];
        dE=dE/(dE*dE+edshift);
        
        for(int v=0; v<n_a; v++)
        for(int w=0; w<n_a; w++){
            //<it|au>(<av|iw>-<av|wi>)
            V= VCAA[((a*n_c+i)*n_a+u)*n_a+t]*  //<it|au>=(ai|tu)
              (VCAA[((a*n_c+i)*n_a+w)*n_a+v]-  //<av|iw>=(ai|wv)
               CAVA[((i*n_a+v)*n_v+a)*n_a+w]); //<av|wi>=(iv|aw)
            
            RF_JM[((t*n_a+u)*n_a+v)*n_a+w]+=V*dE;
        }
    }

    
    for(int i=0; i<n_c; i++)
    for(int a=0; a<n_v; a++)
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++){
        dE=e_c[i]+e_IP[t]-e_v[a]-e_EA[u];
        if(t==u)dE=e_c[i]-e_v[a];
        dE=dE/(dE*dE+edshift);
        
        for(int v=0; v<n_a; v++)
        for(int w=0; w<n_a; w++){
            //<it|au><av|iw>
            V=VCAA[((a*n_c+i)*n_a+u)*n_a+t]*  //<it|au>=(ai|tu)
              VCAA[((a*n_c+i)*n_a+w)*n_a+v];  //<av|iw>=(ai|wv)
            
            RF_JJ[((t*n_a+u)*n_a+v)*n_a+w]+=V*dE;
        }
    }
//     printf_timer("opt table");
    
    
    for(int i=0; i<n_c; i++)
    for(int a=0; a<n_v; a++)
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++){
        //(<it|au>-<iu|at>)*H[i,a]
        V=(VCAA[((a*n_c+i)*n_a+u)*n_a+t]-  //<it|au>=(ai|tu)
           CAVA[((i*n_a+u)*n_v+a)*n_a+t])* //<it|ua>=(iu|at)
           H_CV[i*n_v+a]; 
        dE=e_c[i]+e_IP[t]-e_v[a]-e_EA[u];
        if(t==u)dE=e_c[i]-e_v[a];
        RF_MH[t*n_a+u]+=V*dE/(dE*dE+edshift);
    }
    
    for(int i=0; i<n_c; i++)
    for(int a=0; a<n_v; a++)
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++){
        //(<it|au>-<iu|at>)*H[i,a]
        V=(VCAA[((a*n_c+i)*n_a+u)*n_a+t]-  //<it|au>=(ai|tu)
           CAVA[((i*n_a+u)*n_v+a)*n_a+t])* //<it|ua>=(iu|at)
           H_CV[i*n_v+a]; 
        
//         dE=e_c[i]+e_IP[t]-e_v[a]-e_EA[u];
//         if(t==u)
            dE=e_c[i]-e_v[a];
        RF_HM[t*n_a+u]+=V*dE/(dE*dE+edshift);
    }

    for(int i=0; i<n_c; i++)
    for(int a=0; a<n_v; a++)
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++){
        //(<it|au>-<iu|at>)*H[i,a]
        V= VCAA[((a*n_c+i)*n_a+u)*n_a+t]*  //<it|au>=(ai|tu)
           H_CV[i*n_v+a]; 
        dE=e_c[i]+e_IP[t]-e_v[a]-e_EA[u];
        if(t==u)dE=e_c[i]-e_v[a];
        RF_JH[t*n_a+u]+=V*dE/(dE*dE+edshift);
        
    }
    
    for(int i=0; i<n_c; i++)
    for(int a=0; a<n_v; a++)
    for(int t=0  ;t<n_a;t++)
    for(int u=0  ;u<n_a;u++){
        //(<it|au>-<iu|at>)*H[i,a]
        V= VCAA[((a*n_c+i)*n_a+u)*n_a+t]*  //<it|au>=(ai|tu)
           H_CV[i*n_v+a]; 
//         dE=e_c[i]+e_IP[t]-e_v[a]-e_EA[u];
//         if(t==u)
            dE=e_c[i]-e_v[a];
        RF_HJ[t*n_a+u]+=V*dE/(dE*dE+edshift);
        
    }
    
    
    for(int i=0; i<n_c; i++)
    for(int a=0; a<n_v; a++){
        //H[i,a]*H[i,a]
        V=H_CV[i*n_v+a]*H_CV[i*n_v+a]; 
        dE=e_c[i]-e_v[a];
        RF_H+=V*dE/(dE*dE+edshift);
        
    }
    
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++)
    for(int v=0; v<n_a; v++)
    for(int w=0; w<n_a; w++)
        RF_PV_JK[((u*n_a+t)*n_a+v)*n_a+w]+=+RF_MM[((t*n_a+v)*n_a+u)*n_a+w]
                                           +RF_JJ[((t*n_a+v)*n_a+u)*n_a+w]
                                           -RF_MM[((u*n_a+v)*n_a+t)*n_a+w]
                                           -RF_JJ[((u*n_a+v)*n_a+t)*n_a+w]
                                           -RF_MM[((t*n_a+w)*n_a+u)*n_a+v]
                                           -RF_JJ[((t*n_a+w)*n_a+u)*n_a+v]
                                           +RF_MM[((u*n_a+w)*n_a+t)*n_a+v]
                                           +RF_JJ[((u*n_a+w)*n_a+t)*n_a+v]
                                           ;
                                                       
    
    
    
//     set_zero_matr(RF_PV_AB, n_a*n_a*n_a*n_a*N_fit);
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++)
    for(int v=0; v<n_a; v++)
    for(int w=0; w<n_a; w++)
        RF_PV_AB[((t*n_a+u)*n_a+v)*n_a+w]+=+RF_JM[((t*n_a+v)*n_a+u)*n_a+w]
                                           +RF_MJ[((t*n_a+v)*n_a+u)*n_a+w]
                                           +RF_JM[((u*n_a+w)*n_a+t)*n_a+v]
                                           +RF_MJ[((u*n_a+w)*n_a+t)*n_a+v]
                                           -RF_AB[((t*n_a+w)*n_a+u)*n_a+v]
                                           -RF_AB[((u*n_a+v)*n_a+t)*n_a+w];
    
//     set_zero_matr(RF_PH,N_fit*n_a*n_a);
    for(int t=0; t<n_a; t++)
    for(int w=0; w<n_a; w++){
        RF_PH[t*n_a+w]+= RF_HM[w*n_a+t]
                        +RF_MH[t*n_a+w]
                        +RF_HJ[w*n_a+t]
                        +RF_JH[t*n_a+w];
        for(int u=0; u<n_a; u++)
             RF_PH[t*n_a+w]+= RF_MM[((t*n_a+u)*n_a+u)*n_a+w]
                             +RF_JJ[((t*n_a+u)*n_a+u)*n_a+w]
                             +RF_AB[((t*n_a+u)*n_a+u)*n_a+w];
    }
    
    RF_PS+=2*RF_H;
    
    delete[] RF_MM ;
    delete[] RF_JJ ;
    delete[] RF_JM ;
    delete[] RF_MJ ;
    delete[] RF_AB ;
    delete[] RF_HM ;
    delete[] RF_MH ;
    delete[] RF_HJ ;
    delete[] RF_JH ;
    delete[] VCAA  ;
    delete[] CAVA  ;
//     delete[] VCAA_J;
//     delete[] VCAA_K;
    
    return 0;
    
}

int PT_tensors::calc_EE_2_AV(){   //// AV
    
    
    
    int aux_n_ao = integrals->aux_n_ao;
    double * VA=integrals->VA_RI_M;
    double * CC=integrals->CC_RI_M;
    double * VC=integrals->VC_RI_M;
    double * CA=integrals->CA_RI_M;
    double * AA=integrals->AA_RI_M;
    
    
    double * VAAA;
    VAAA = new double[n_v*n_a*n_a*n_a];
    if(RI)
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        n_v*n_a,n_a*n_a,aux_n_ao,1.0,
                        VA,aux_n_ao,
                        AA,aux_n_ao,0.0,
                        VAAA,n_a*n_a);                          //// ??? (at|K)(K|uv)  2el integrals??? 3rd term in dipole sum
    else{
        double * AAAV = new double[n_v*n_a*n_a*n_a];
        calc_2el_MO_INTS_XXXY(integrals->s, integrals->n_ao, AAAV, integrals->ACT_CVEC, integrals->VAC_CVEC, n_a, n_v);
        for(int a=0;a<n_v;a++)
        for(int t=0;t<n_a;t++)
        for(int u=0;u<n_a;u++)
        for(int v=0;v<n_a;v++)
            VAAA[((a*n_a+t)*n_a+u)*n_a+v]=AAAV[((u*n_a+v)*n_a+t)*n_v+a];
        delete[] AAAV;
    }

    
    double * RF_P3_JJ;
    double * RF_P3_JM;
    double * RF_P3_MJ;
    double * RF_P3_HM;
    double * RF_P3_MH;
    double * RF_P3_HJ;
    double * RF_P3_JH;
    
    RF_P3_JJ = new double[n_a*n_a*n_a*n_a*n_a*n_a];
    RF_P3_JM = new double[n_a*n_a*n_a*n_a*n_a*n_a];
    RF_P3_MJ = new double[n_a*n_a*n_a*n_a*n_a*n_a];
    
    RF_P3_HM = new double[n_a*n_a*n_a*n_a];
    RF_P3_MH = new double[n_a*n_a*n_a*n_a];
    RF_P3_HJ = new double[n_a*n_a*n_a*n_a];
    RF_P3_JH = new double[n_a*n_a*n_a*n_a];
    
    set_zero_matr(RF_P3_JJ, n_a*n_a*n_a*n_a*n_a*n_a);
    set_zero_matr(RF_P3_JM, n_a*n_a*n_a*n_a*n_a*n_a);
    set_zero_matr(RF_P3_MJ, n_a*n_a*n_a*n_a*n_a*n_a);
    set_zero_matr(RF_P3_HM, n_a*n_a*n_a*n_a);
    set_zero_matr(RF_P3_MH, n_a*n_a*n_a*n_a);
    set_zero_matr(RF_P3_HJ, n_a*n_a*n_a*n_a);
    set_zero_matr(RF_P3_JH, n_a*n_a*n_a*n_a);
    
    
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
            v=vu/n_a;     //// ???? çà÷åì íóæíû â òîì ÷èñëî îñòàòêè îò äåëåíèÿ
            u=vu%n_a;
            for(int a=0; a<n_v; a++)
            for(int t=0; t<n_a; t++){
                dE=-e_v[a]+e_IP[t]+e_IP[u]-e_EA[v];
                if(v==t)dE=-e_v[a]+e_IP[u];
                if(v==u)dE=-e_v[a]+e_IP[t];
                
                
                dE=dE/(dE*dE+edshift);
                
                for(int y=0; y<n_a; y++)
                for(int x=0; x<n_a; x++)
                for(int w=0; w<n_a; w++){
                    //<tu||av>=<tu|av>-<tu||va>=(ta|vu)-(tv|au)
                    //<xy||aw>=<xy|aw>-<xy||wa>=(xa|yw)-(xw|ay)
                    V =(VAAA[((a*n_a+t)*n_a+v)*n_a+u]-
                        VAAA[((a*n_a+u)*n_a+v)*n_a+t])*
                       (VAAA[((a*n_a+x)*n_a+w)*n_a+y]-
                        VAAA[((a*n_a+y)*n_a+w)*n_a+x]);
                       //a(w)a+(x)a+(y)a+(v)a(u)a(t)
                       RF_P3_JK[((((t*n_a+u)*n_a+v)*n_a+y)*n_a+x)*n_a+w]+=V*dE;
                    
                }
            }
        
            //JJ
            for(int a=0; a<n_v; a++)
            for(int t=0; t<n_a; t++){
                dE=-e_v[a]+e_IP[t]+e_IP[u]-e_EA[v];
                if(v==t)dE=-e_v[a]+e_IP[u];
                if(v==u)dE=-e_v[a]+e_IP[t];
                dE=dE/(dE*dE+edshift);
                
                for(int y=0; y<n_a; y++)
                for(int x=0; x<n_a; x++)
                for(int w=0; w<n_a; w++){
                    //<a(a)w(b)|x(a)y(b)>=(yw|xa)
                    V=VAAA[((a*n_a+u)*n_a+v)*n_a+t]*
                      VAAA[((a*n_a+x)*n_a+w)*n_a+y];
                    
                    RF_P3_JJ[((((t*n_a+u)*n_a+v)*n_a+y)*n_a+x)*n_a+w]+=V*dE;
                    
                }
            }
            //MJ
            for(int a=0; a<n_v; a++)
            for(int t=0; t<n_a; t++){
                dE=-e_v[a]+e_IP[t]+e_IP[u]-e_EA[v];
                if(v==t)dE=-e_v[a]+e_IP[u];
                if(v==u)dE=-e_v[a]+e_IP[t];
                dE=dE/(dE*dE+edshift);
                
                for(int y=0; y<n_a; y++)
                for(int x=0; x<n_a; x++)
                for(int w=0; w<n_a; w++){
                    //<a(a)w(b)|x(a)y(b)>=(yw|xa)
                    V=(VAAA[((a*n_a+t)*n_a+v)*n_a+u]-
                       VAAA[((a*n_a+u)*n_a+v)*n_a+t])*
                       VAAA[((a*n_a+x)*n_a+w)*n_a+y];
            
                    RF_P3_MJ[((((t*n_a+u)*n_a+v)*n_a+y)*n_a+x)*n_a+w]+=V*dE;
                    
                }
            }
            //JM
            for(int a=0; a<n_v; a++)
            for(int t=0; t<n_a; t++){
                dE=-e_v[a]+e_IP[t]+e_IP[u]-e_EA[v];
                if(v==t)dE=-e_v[a]+e_IP[u];
                if(v==u)dE=-e_v[a]+e_IP[t];
                dE=dE/(dE*dE+edshift);
                
                for(int y=0; y<n_a; y++)
                for(int x=0; x<n_a; x++)
                for(int w=0; w<n_a; w++){
                    //<aw|xy>-<aw|yx>=(yw|xa)-(xw|ya)
                    V= VAAA[((a*n_a+u)*n_a+v)*n_a+t]*
                      (VAAA[((a*n_a+x)*n_a+w)*n_a+y]-
                       VAAA[((a*n_a+y)*n_a+w)*n_a+x]);
                    
                    RF_P3_JM[((((t*n_a+u)*n_a+v)*n_a+y)*n_a+x)*n_a+w]+=V*dE;
                    
                }
            }
        }    
        //H
        for(int a=0; a<n_v; a++)
        for(int t=th_id; t<n_a; t+=num_threads){
            dE=-e_v[a]+e_IP[t];
            dE=dE/(dE*dE+edshift);
            
            for(int w=0; w<n_a; w++){
                H1a=H_AV[t*n_v+a];
                H2a=H_AV[w*n_v+a];
                            
                V=H1a*H2a;
                RF_PH[t*n_a+w]+=V*dE;
                
            }
        }
        //HM
        for(int a=0; a<n_v; a++)
        for(int t=th_id; t<n_a; t+=num_threads){
            dE=-e_v[a]+e_IP[t];
            dE=dE/(dE*dE+edshift);
            
            for(int y=0; y<n_a; y++)
            for(int x=0; x<n_a; x++)
            for(int w=0; w<n_a; w++){
                H1a=H_AV[t*n_v+a];
     
                //<aw|xy>-<aw|yx>=(yw|xa)-(xw|ya)
                V=H1a*
                  (VAAA[((a*n_a+x)*n_a+w)*n_a+y]-
                   VAAA[((a*n_a+y)*n_a+w)*n_a+x]);
                
                RF_P3_HM[((t*n_a+y)*n_a+x)*n_a+w]+=V*dE;
                
            }
        }
        //MH
        for(int vu=th_id; vu<n_a*n_a; vu+=num_threads){
            v=vu/n_a;
            u=vu%n_a;
            for(int a=0; a<n_v; a++)
            for(int t=0; t<n_a; t++){
                dE=-e_v[a]+e_IP[t]+e_IP[u]-e_EA[v];
                if(v==t)dE=-e_v[a]+e_IP[u];
                if(v==u)dE=-e_v[a]+e_IP[t];
                dE=dE/(dE*dE+edshift);
                
                for(int w=0; w<n_a; w++){
                    H2a=H_AV[w*n_v+a];
                    //<aw|xy>-<aw|yx>=(yw|xa)-(xw|ya)
                    V=(VAAA[((a*n_a+t)*n_a+v)*n_a+u]-
                       VAAA[((a*n_a+u)*n_a+v)*n_a+t])*
                       H2a;
                    
                    RF_P3_MH[((t*n_a+u)*n_a+v)*n_a+w]+=V*dE;
                }
            }
        }
        //HJ
        for(int a=0; a<n_v; a++)
        for(int t=th_id; t<n_a; t+=num_threads){
            dE=-e_v[a]+e_IP[t];
            dE=dE/(dE*dE+edshift);
            
            for(int y=0; y<n_a; y++)
            for(int x=0; x<n_a; x++)
            for(int w=0; w<n_a; w++){
                H1a=H_AV[t*n_v+a];
                //<a(a)w(b)|x(a)y(b)>=(yw|xa)
                V=H1a*
                  VAAA[((a*n_a+x)*n_a+w)*n_a+y];
                
                RF_P3_HJ[((t*n_a+y)*n_a+x)*n_a+w]+=V*dE;
            }
        }
        //JH
        for(int vu=th_id; vu<n_a*n_a; vu+=num_threads){
            v=vu/n_a;
            u=vu%n_a;
            for(int a=0; a<n_v; a++)
            for(int t=0; t<n_a; t++){
                dE=-e_v[a]+e_IP[t]+e_IP[u]-e_EA[v];
                if(v==t)dE=-e_v[a]+e_IP[u];
                if(v==u)dE=-e_v[a]+e_IP[t];
                dE=dE/(dE*dE+edshift);
                
                for(int w=0; w<n_a; w++){
                    H2a=H_AV[w*n_v+a];     
                    //<a(a)w(b)|x(a)y(b)>=(yw|xa)
                    V=VAAA[((a*n_a+u)*n_a+v)*n_a+t]*
                      H2a;

                    RF_P3_JH[((t*n_a+u)*n_a+v)*n_a+w]+=V*dE;
                }
            }
        }
    }
    //RF_PV_JK
    for(int w=0; w<n_a; w++)
    for(int t=0; t<n_a; t++)
    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++)
        RF_PV_JK[((w*n_a+t)*n_a+v)*n_a+u]+= RF_P3_HM[((t*n_a+u)*n_a+v)*n_a+w]
                                           -RF_P3_HM[((w*n_a+u)*n_a+v)*n_a+t]
                                           -RF_P3_MH[((t*n_a+w)*n_a+v)*n_a+u]
                                           +RF_P3_MH[((t*n_a+w)*n_a+u)*n_a+v];
    //RF_PV_AB
    for(int w=0; w<n_a; w++)
    for(int t=0; t<n_a; t++)
    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++){
        RF_PV_AB[((t*n_a+u)*n_a+v)*n_a+w]+= RF_P3_HJ[((t*n_a+w)*n_a+v)*n_a+u]
                                           +RF_P3_JH[((u*n_a+t)*n_a+w)*n_a+v]
                                           +RF_P3_HJ[((u*n_a+v)*n_a+w)*n_a+t]
                                           +RF_P3_JH[((t*n_a+u)*n_a+v)*n_a+w];
    
        for(int x=0; x<n_a; x++)
            RF_PV_AB[((t*n_a+u)*n_a+v)*n_a+w] += RF_P3_JJ[((((t*n_a+u)*n_a+x)*n_a+v)*n_a+w)*n_a+x]
                                                +RF_P3_JJ[((((u*n_a+t)*n_a+x)*n_a+w)*n_a+v)*n_a+x];
    }
    //RF_P3_AB - order differs with the original XMCQDPT code (res_fit.cpp)
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++)
    for(int x=0; x<n_a; x++)
    for(int v=0; v<n_a; v++)
    for(int w=0; w<n_a; w++)
    for(int y=0; y<n_a; y++)
        //           A     A      B      A      A      B
        RF_P3_AB[((((t*n_a+u)*n_a+x)*n_a+v)*n_a+w)*n_a+y]+=-RF_P3_MJ[((((t*n_a+u)*n_a+v)*n_a+y)*n_a+w)*n_a+x]
                                                           +RF_P3_MJ[((((t*n_a+u)*n_a+w)*n_a+y)*n_a+v)*n_a+x]
                                                           +RF_P3_JM[((((x*n_a+t)*n_a+y)*n_a+w)*n_a+v)*n_a+u]
                                                           -RF_P3_JM[((((x*n_a+u)*n_a+y)*n_a+w)*n_a+v)*n_a+t]
                                                           +RF_P3_JJ[((((t*n_a+x)*n_a+v)*n_a+w)*n_a+y)*n_a+u]
                                                           -RF_P3_JJ[((((t*n_a+x)*n_a+w)*n_a+v)*n_a+y)*n_a+u]
                                                           -RF_P3_JJ[((((u*n_a+x)*n_a+v)*n_a+w)*n_a+y)*n_a+t]
                                                           +RF_P3_JJ[((((u*n_a+x)*n_a+w)*n_a+v)*n_a+y)*n_a+t];
                                    
    
    delete[] VAAA;
    delete[] RF_P3_JJ;
    delete[] RF_P3_JM;
    delete[] RF_P3_MJ;
    delete[] RF_P3_HM;
    delete[] RF_P3_MH;
    delete[] RF_P3_HJ;
    delete[] RF_P3_JH;
    
    return 0;
    
}

int PT_tensors::calc_EE_2_CA(){
    
    double dE;
    double H1a;
    double H2a;
    double V;
    
    
    int aux_n_ao = integrals->aux_n_ao;

    double * CC=integrals->CC_RI_M;

    double * CA=integrals->CA_RI_M;
    double * AA=integrals->AA_RI_M;

    
    double * CAAA;
    CAAA = new double[n_c*n_a*n_a*n_a];
    
    if(RI)
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        n_c*n_a,n_a*n_a,aux_n_ao,1.0,
                        CA,aux_n_ao,
                        AA,aux_n_ao,0.0,
                        CAAA,n_a*n_a);    
    else{
        double * AAAC = new double[n_c*n_a*n_a*n_a];
        calc_2el_MO_INTS_XXXY(integrals->s, integrals->n_ao, AAAC, integrals->ACT_CVEC, integrals->COR_CVEC, n_a, n_c);
        for(int i=0;i<n_c;i++)
        for(int t=0;t<n_a;t++)
        for(int u=0;u<n_a;u++)
        for(int v=0;v<n_a;v++)
            CAAA[((i*n_a+t)*n_a+u)*n_a+v]=AAAC[((u*n_a+v)*n_a+t)*n_c+i];
            
        delete[] AAAC;
    }
        
    double * RF_P3_ss;
    double * RF_P3_os;
    double * RF_P3_JH;
    double * RF_P3_MH;
    double * RF_P3_HJ;
    double * RF_P3_HM;
    double * RF_P3_HH;
    
    RF_P3_ss = new double[n_a*n_a*n_a*n_a*n_a*n_a];
    RF_P3_os = new double[n_a*n_a*n_a*n_a*n_a*n_a];
    RF_P3_JH = new double[n_a*n_a*n_a*n_a        ];
    RF_P3_MH = new double[n_a*n_a*n_a*n_a        ];
    RF_P3_HJ = new double[n_a*n_a*n_a*n_a        ];
    RF_P3_HM = new double[n_a*n_a*n_a*n_a        ];
    RF_P3_HH = new double[n_a*n_a                ];
    
    set_zero_matr(RF_P3_ss, n_a*n_a*n_a*n_a*n_a*n_a);
    set_zero_matr(RF_P3_os, n_a*n_a*n_a*n_a*n_a*n_a);
    set_zero_matr(RF_P3_JH, n_a*n_a*n_a*n_a        );
    set_zero_matr(RF_P3_MH, n_a*n_a*n_a*n_a        );
    set_zero_matr(RF_P3_HJ, n_a*n_a*n_a*n_a        );
    set_zero_matr(RF_P3_HM, n_a*n_a*n_a*n_a        );
    set_zero_matr(RF_P3_HH, n_a*n_a                );
    
    
    
    //VV same spin
    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++)
    for(int i=0; i<n_c; i++)
    for(int t=0; t<n_a; t++){
        dE= e_c[i]+e_IP[t]-e_EA[u]-e_EA[v];
        if(t==u)dE= e_c[i]-e_EA[v];
        if(t==v)dE= e_c[i]-e_EA[u];
        dE=dE/(dE*dE+edshift);
        
        for(int y=0; y<n_a; y++)
        for(int x=0; x<n_a; x++)
        for(int w=0; w<n_a; w++){
            // i A
            // v A
            // u A
            // t A
            // w A
            // x A
            // y A
            V=(CAAA[((i*n_a+v)*n_a+u)*n_a+t]*
               CAAA[((i*n_a+w)*n_a+x)*n_a+y]);
            RF_P3_ss[((((t*n_a+u)*n_a+v)*n_a+y)*n_a+x)*n_a+w]+=V*dE;
            
        }
    }
    //VV - oposite spin
    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++)
    for(int i=0; i<n_c; i++)
    for(int t=0; t<n_a; t++){
        dE= e_c[i]+e_IP[t]-e_EA[u]-e_EA[v];
        if(t==u)dE= e_c[i]-e_EA[v];
        if(t==v)dE= e_c[i]-e_EA[u];
        dE=dE/(dE*dE+edshift);
        
        for(int y=0; y<n_a; y++)
        for(int x=0; x<n_a; x++)
        for(int w=0; w<n_a; w++){
            // i A
            // v A
            // u B
            // t B
            // w any
            // x any
            // y any
            V=(CAAA[((i*n_a+v)*n_a+u)*n_a+t]*
               CAAA[((i*n_a+w)*n_a+x)*n_a+y]);
            RF_P3_os[((((t*n_a+u)*n_a+v)*n_a+y)*n_a+x)*n_a+w]+=V*dE;
            
        }
    }
    
    
    //MH
    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++)
    for(int i=0; i<n_c; i++)
    for(int t=0; t<n_a; t++){
        dE= e_c[i]+e_IP[t]-e_EA[u]-e_EA[v];
        if(t==u)dE= e_c[i]-e_EA[v];
        if(t==v)dE= e_c[i]-e_EA[u];
        dE=dE/(dE*dE+edshift);
        for(int w=0; w<n_a; w++){
            H2a=H_CA[i*n_a+w];
            V=(CAAA[((i*n_a+v)*n_a+u)*n_a+t]
              -CAAA[((i*n_a+u)*n_a+v)*n_a+t])*
               H2a;
            RF_P3_MH[((t*n_a+u)*n_a+v)*n_a+w]+=V*dE;
        }
    }
    //JH
    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++)
    for(int i=0; i<n_c; i++)
    for(int t=0; t<n_a; t++){
        dE= e_c[i]+e_IP[t]-e_EA[u]-e_EA[v];
        if(t==u)dE= e_c[i]-e_EA[v];
        if(t==v)dE= e_c[i]-e_EA[u];
        dE=dE/(dE*dE+edshift);
        for(int w=0; w<n_a; w++){
            H2a=H_CA[i*n_a+w];
            V=(CAAA[((i*n_a+v)*n_a+u)*n_a+t]*
               H2a);
            RF_P3_JH[((t*n_a+u)*n_a+v)*n_a+w]+=V*dE;
        }
    }

    //HM
    for(int i=0; i<n_c; i++)
    for(int t=0; t<n_a; t++){
        dE= e_c[i]-e_EA[t];
        dE=dE/(dE*dE+edshift);
        
        for(int y=0; y<n_a; y++)
        for(int x=0; x<n_a; x++)
        for(int w=0; w<n_a; w++){
            H1a=H_CA[i*n_a+t];
            V= H1a*
              (CAAA[((i*n_a+w)*n_a+x)*n_a+y]
              -CAAA[((i*n_a+x)*n_a+w)*n_a+y]);
            RF_P3_HM[((t*n_a+y)*n_a+x)*n_a+w]+=V*dE;
        }
    }
    //HJ
    for(int i=0; i<n_c; i++)
    for(int t=0; t<n_a; t++){
        dE= e_c[i]-e_EA[t];
        dE=dE/(dE*dE+edshift);
        
        for(int y=0; y<n_a; y++)
        for(int x=0; x<n_a; x++)
        for(int w=0; w<n_a; w++){
            H1a=H_CA[i*n_a+t];
            V=(H1a*
               CAAA[((i*n_a+w)*n_a+x)*n_a+y]);
            RF_P3_HJ[((t*n_a+y)*n_a+x)*n_a+w]+=V*dE;
        }
    }
    
    //HH
    for(int i=0; i<n_c; i++)
    for(int t=0; t<n_a; t++){
        dE= e_c[i]-e_EA[t];
        dE=dE/(dE*dE+edshift);
        
        for(int u=0; u<n_a; u++){
            H1a=H_CA[i*n_a+t];
            H2a=H_CA[i*n_a+u];
            V=H1a*H2a;
            RF_P3_HH[t*n_a+u]+=V*dE;
        }
    }
    
    
    for(int t=0; t<n_a; t++)
    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++)
    for(int w=0; w<n_a; w++)
    for(int x=0; x<n_a; x++)
    for(int y=0; y<n_a; y++)
        //a(w)a+(x)a+(y)a+(v)a(u)a(t) - RF_P3_VV[((((t*n_a+v)*n_a+x)*n_a+y)*n_a+w)*n_a+u](t->u, x->y)
        RF_P3_JK[((((t*n_a+u)*n_a+v)*n_a+y)*n_a+x)*n_a+w]+= RF_P3_ss[((((t*n_a+v)*n_a+x)*n_a+y)*n_a+w)*n_a+u]
                                                           -RF_P3_ss[((((t*n_a+v)*n_a+y)*n_a+x)*n_a+w)*n_a+u]
                                                           -RF_P3_ss[((((u*n_a+v)*n_a+x)*n_a+y)*n_a+w)*n_a+t]
                                                           +RF_P3_ss[((((u*n_a+v)*n_a+y)*n_a+x)*n_a+w)*n_a+t];
    
    for(int t=0; t<n_a; t++)
    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++)
    for(int w=0; w<n_a; w++)
    for(int x=0; x<n_a; x++)
        RF_PV_JK[((u*n_a+t)*n_a+v)*n_a+w]+=-RF_P3_ss[((((t*n_a+v)*n_a+x)*n_a+w)*n_a+x)*n_a+u] //AAAAAA
                                           +RF_P3_ss[((((u*n_a+v)*n_a+x)*n_a+w)*n_a+x)*n_a+t] //AAAAAA
                                           +RF_P3_ss[((((t*n_a+w)*n_a+x)*n_a+v)*n_a+x)*n_a+u] //AAAAAA
                                           -RF_P3_ss[((((u*n_a+w)*n_a+x)*n_a+v)*n_a+x)*n_a+t] //AAAAAA
                                           -RF_P3_ss[((((t*n_a+x)*n_a+v)*n_a+w)*n_a+u)*n_a+x] //AAAAAA
                                           +RF_P3_ss[((((t*n_a+v)*n_a+x)*n_a+w)*n_a+u)*n_a+x] //AAAAAA 1
                                           +RF_P3_ss[((((u*n_a+x)*n_a+v)*n_a+w)*n_a+t)*n_a+x] //AAAAAA
                                           -RF_P3_ss[((((u*n_a+v)*n_a+x)*n_a+w)*n_a+t)*n_a+x] //AAAAAA 2
                                           +RF_P3_ss[((((t*n_a+x)*n_a+w)*n_a+v)*n_a+u)*n_a+x] //AAAAAA
                                           -RF_P3_ss[((((t*n_a+w)*n_a+x)*n_a+v)*n_a+u)*n_a+x] //AAAAAA 3
                                           -RF_P3_ss[((((u*n_a+x)*n_a+w)*n_a+v)*n_a+t)*n_a+x] //AAAAAA
                                           +RF_P3_ss[((((u*n_a+w)*n_a+x)*n_a+v)*n_a+t)*n_a+x] //AAAAAA 4
                                           +RF_P3_os[((((t*n_a+v)*n_a+x)*n_a+w)*n_a+u)*n_a+x] //AABAAB 1
                                           -RF_P3_os[((((u*n_a+v)*n_a+x)*n_a+w)*n_a+t)*n_a+x] //AABAAB 2
                                           -RF_P3_os[((((t*n_a+w)*n_a+x)*n_a+v)*n_a+u)*n_a+x] //AABAAB 3
                                           +RF_P3_os[((((u*n_a+w)*n_a+x)*n_a+v)*n_a+t)*n_a+x];//AABAAB 4
    for(int t=0; t<n_a; t++)
    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++)
    for(int w=0; w<n_a; w++)
        RF_PV_JK[((u*n_a+t)*n_a+v)*n_a+w]+=+RF_P3_MH[((t*n_a+w)*n_a+v)*n_a+u]
                                           +RF_P3_MH[((u*n_a+v)*n_a+w)*n_a+t]
                                           +RF_P3_HM[((v*n_a+w)*n_a+t)*n_a+u]
                                           +RF_P3_HM[((w*n_a+v)*n_a+u)*n_a+t];

    
//  RF_P3_AB                 - order differs with the original XMCQDPT code (res_fit.cpp)
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++)
    for(int x=0; x<n_a; x++)
    for(int v=0; v<n_a; v++)
    for(int w=0; w<n_a; w++)
    for(int y=0; y<n_a; y++)
        //           A     A      B      A      A      B
        RF_P3_AB[((((t*n_a+u)*n_a+x)*n_a+v)*n_a+w)*n_a+y]+=-RF_P3_os[((((t*n_a+v)*n_a+y)*n_a+w)*n_a+u)*n_a+x]
                                                           +RF_P3_os[((((u*n_a+v)*n_a+y)*n_a+w)*n_a+t)*n_a+x]
                                                           +RF_P3_os[((((t*n_a+w)*n_a+y)*n_a+v)*n_a+u)*n_a+x]
                                                           -RF_P3_os[((((u*n_a+w)*n_a+y)*n_a+v)*n_a+t)*n_a+x]
                                                           -RF_P3_ss[((((t*n_a+v)*n_a+w)*n_a+y)*n_a+x)*n_a+u]
                                                           +RF_P3_ss[((((u*n_a+v)*n_a+w)*n_a+y)*n_a+x)*n_a+t]
                                                           +RF_P3_ss[((((t*n_a+w)*n_a+v)*n_a+y)*n_a+x)*n_a+u]
                                                           -RF_P3_ss[((((u*n_a+w)*n_a+v)*n_a+y)*n_a+x)*n_a+t]
                                                           -RF_P3_os[((((x*n_a+y)*n_a+v)*n_a+w)*n_a+u)*n_a+t]
                                                           +RF_P3_os[((((x*n_a+y)*n_a+v)*n_a+w)*n_a+t)*n_a+u]
                                                           +RF_P3_os[((((x*n_a+y)*n_a+w)*n_a+v)*n_a+u)*n_a+t]
                                                           -RF_P3_os[((((x*n_a+y)*n_a+w)*n_a+v)*n_a+t)*n_a+u];

//     set_zero_matr(RF_PV_AB, n_a*n_a*n_a*n_a*N_fit);
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++)
    for(int v=0; v<n_a; v++)
    for(int w=0; w<n_a; w++)
    for(int x=0; x<n_a; x++)
        //         A     B      A      B
        RF_PV_AB[((t*n_a+u)*n_a+v)*n_a+w]+=-RF_P3_os[((((t*n_a+x)*n_a+w)*n_a+v)*n_a+x)*n_a+u]
                                           -RF_P3_ss[((((t*n_a+x)*n_a+v)*n_a+w)*n_a+u)*n_a+x]
                                           +RF_P3_ss[((((t*n_a+v)*n_a+x)*n_a+w)*n_a+u)*n_a+x] //1st - ss
                                           -RF_P3_ss[((((u*n_a+w)*n_a+x)*n_a+v)*n_a+x)*n_a+t] //2nd - ex ss
                                           +RF_P3_ss[((((u*n_a+w)*n_a+x)*n_a+v)*n_a+t)*n_a+x] //2nd - ss
                                           -RF_P3_os[((((u*n_a+x)*n_a+v)*n_a+w)*n_a+x)*n_a+t]
                                           -RF_P3_ss[((((u*n_a+x)*n_a+w)*n_a+v)*n_a+t)*n_a+x]
                                           +RF_P3_os[((((u*n_a+w)*n_a+x)*n_a+v)*n_a+t)*n_a+x] //2nd - os
                                           -RF_P3_ss[((((t*n_a+v)*n_a+x)*n_a+w)*n_a+x)*n_a+u] //1st - ex ss
                                           +RF_P3_os[((((t*n_a+v)*n_a+x)*n_a+w)*n_a+u)*n_a+x];//1st - os
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++)
    for(int v=0; v<n_a; v++)
    for(int w=0; w<n_a; w++)
        RF_PV_AB[((t*n_a+u)*n_a+v)*n_a+w]+=-RF_P3_JH[((u*n_a+w)*n_a+v)*n_a+t]
                                           -RF_P3_HJ[((v*n_a+w)*n_a+u)*n_a+t]
                                           -RF_P3_JH[((t*n_a+v)*n_a+w)*n_a+u]
                                           -RF_P3_HJ[((w*n_a+v)*n_a+t)*n_a+u];
    
//     set_zero_matr(RF_PH, n_a*n_a*N_fit);
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++)
    for(int v=0; v<n_a; v++)
    for(int y=0; y<n_a; y++)
        RF_PH[t*n_a+u]+=-RF_P3_ss[((((t*n_a+y)*n_a+v)*n_a+u)*n_a+v)*n_a+y]
                        +RF_P3_ss[((((t*n_a+v)*n_a+y)*n_a+u)*n_a+v)*n_a+y]
                        +RF_P3_os[((((t*n_a+v)*n_a+y)*n_a+u)*n_a+v)*n_a+y];
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++)
    for(int y=0; y<n_a; y++)
        RF_PH[t*n_a+u]+= RF_P3_MH[((t*n_a+u)*n_a+y)*n_a+y]
                        +RF_P3_JH[((t*n_a+u)*n_a+y)*n_a+y]
                        +RF_P3_HM[((y*n_a+u)*n_a+t)*n_a+y]
                        +RF_P3_HJ[((y*n_a+u)*n_a+t)*n_a+y];
    for(int t=0; t<n_a; t++)
    for(int u=0; u<n_a; u++)
        RF_PH[t*n_a+u]+=-RF_P3_HH[u*n_a+t];
    
//     set_zero_matr(RF_PS, N_fit);
    for(int t=0; t<n_a; t++)
        RF_PS+=RF_P3_HH[t*n_a+t]*2;
    

    delete[] CAAA  ;

    delete[] RF_P3_ss;
    delete[] RF_P3_os;
    delete[] RF_P3_JH;
    delete[] RF_P3_MH;
    delete[] RF_P3_HJ;
    delete[] RF_P3_HM;
    delete[] RF_P3_HH;

    
    return 0;
}

int PT_tensors::calc_EE_2_CCVV(){
    
    double * R;
    R = new double[num_threads];
    set_zero_matr(R,num_threads);
    double *V[num_threads];
    double *JK[num_threads];
    for(int i=0;i<num_threads;i++) V[i] = new double[n_c];
    for(int i=0;i<num_threads;i++)JK[i] = new double[n_c*n_c];
    
    double * IEIE;
    
    if(RI==0){
        int n_ao = integrals->n_ao;
        int n_o  = n_c+n_a;//occupied
        double * CVEC = new double[(n_c+n_a)*n_ao];
        for(int i=0; i<n_ao ;i++){
            for(int j=0; j<n_c;j++) CVEC[i*(n_c+n_a)+j    ]=integrals->COR_CVEC[i*n_c+j];
            for(int j=0; j<n_a;j++) CVEC[i*(n_c+n_a)+j+n_c]=integrals->ACT_CVEC[i*n_a+j];
        }
        IEIE = new double[n_v*n_o*n_v*n_o];
        calc_2el_MO_INTS_XYXY(integrals->s, n_ao, IEIE, CVEC, integrals->VAC_CVEC, n_o, n_v);
        
        delete[] CVEC;
        
    }
    
    
    for(int a=0; a<n_v; a++){//fprintf(stderr,"CCVV a=%d\r",a);
#pragma omp parallel for
    for(int b=a; b<n_v; b++){/*fprintf(stderr,"CCVV a,b=%d,%d\r",a,b);*/
        
        int th_id = omp_get_thread_num();
        double C;
        int n_o = n_c+n_a;//occupied
        int iajb=0;
        double dE;
        double Ep;
        double J;
        double K;
        
        if(RI)
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                        n_c,n_c,integrals->aux_n_ao,1.0,
                        integrals->VC_RI_M+a*n_c*integrals->aux_n_ao,integrals->aux_n_ao,
                        integrals->VC_RI_M+b*n_c*integrals->aux_n_ao,integrals->aux_n_ao,0.0,
                        JK[th_id],n_c);
        else{
            for(int i=0; i<n_c; i++)
            for(int j=i; j<n_c; j++){
                iajb = (i*n_o*n_v+j)*n_v;
                JK[th_id][i*n_c+j]=IEIE[iajb+a*n_o*n_v+b];
                JK[th_id][j*n_c+i]=IEIE[iajb+b*n_o*n_v+a];
            }
        }
            
            
        
        for(int i=0; i<n_c; i++){//fprintf(stderr,"i,a=%d,%d\r",i,a);
        for(int j=i; j<n_c; j++){
            J=JK[th_id][i*n_c+j];
            K=JK[th_id][j*n_c+i];
            V[th_id][j]=2*(J*J+K*K);
            
            if(a==b)V[th_id][j]=V[th_id][j]/2;
            if(i==j)V[th_id][j]=V[th_id][j]/2;
            C=J-K;
            V[th_id][j]+=2*C*C;
        }
        
        Ep=e_c[i]-e_v[a]-e_v[b];
        for(int j=i; j<n_c; j++){
            dE=e_c[j]+Ep;
            R[th_id]+=V[th_id][j]*dE/(dE*dE+edshift);
            
        }
       
        }
    }
    }
    if(RI==0)delete[] IEIE;
    
    for(int i=0;i<num_threads;i++)
        RF_PS+=R[i];
    
    delete[] R;
    for(int i=0;i<num_threads;i++) delete[] V [i];
    for(int i=0;i<num_threads;i++) delete[] JK[i];
    
    
    return 0;
}

int PT_tensors::calc_EE_2_CAVV(){
    double **PH_th= new double * [num_threads];
    for(int i=0;i<num_threads;i++)PH_th[i]=new double[n_a*n_a];
    for(int i=0;i<num_threads;i++)set_zero_matr(PH_th[i],n_a*n_a);
    
    int n_o=n_c+n_a;//occupied
    
    double * IEIE;
    
    if(RI==0){
        int n_ao = integrals->n_ao;
        double * CVEC = new double[(n_c+n_a)*n_ao];
        for(int i=0; i<n_ao ;i++){
            for(int j=0; j<n_c;j++) CVEC[i*(n_c+n_a)+j    ]=integrals->COR_CVEC[i*n_c+j];
            for(int j=0; j<n_a;j++) CVEC[i*(n_c+n_a)+j+n_c]=integrals->ACT_CVEC[i*n_a+j];
        }
        IEIE = new double[n_v*n_o*n_v*n_o];
        calc_2el_MO_INTS_XYXY(integrals->s, n_ao, IEIE, CVEC, integrals->VAC_CVEC, n_o, n_v);
        
        delete[] CVEC;
        
    }
    
    
    
    #pragma omp parallel
    {
        int nt = omp_get_thread_num();
        
        double *K;
        double *J;
        K= new double[n_a*n_c];
        J= new double[n_a*n_c];
        
        double dE;
        double V2;
        
//         fprintf(out_stream,"CAVV is not parallel\n");
        
        for(int a=nt  ; a<n_v; a+=num_threads)
        for(int b=a  ; b<n_v; b++){//fprintf(stderr,"CAVVa,b=%d,%d\r",a,b);
            
            if(RI){
                cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            n_c,n_a,integrals->aux_n_ao,1.0,
                            integrals->VC_RI_M+a*n_c*integrals->aux_n_ao,integrals->aux_n_ao,
                            integrals->VA_RI_M+b*n_a*integrals->aux_n_ao,integrals->aux_n_ao,0.0,
                            J,n_a);
                cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            n_c,n_a,integrals->aux_n_ao,1.0,
                            integrals->VC_RI_M+b*n_c*integrals->aux_n_ao,integrals->aux_n_ao,
                            integrals->VA_RI_M+a*n_a*integrals->aux_n_ao,integrals->aux_n_ao,0.0,
                            K,n_a);
            }
            else{
                for(int i=0  ; i<n_c; i++)
                for(int t=0  ; t<n_a; t++){
                    J[i*n_a+t]=IEIE[((i*n_v+a)*n_o+t+n_c)*n_v+b];
                    K[i*n_a+t]=IEIE[((i*n_v+b)*n_o+t+n_c)*n_v+a];
                }
            }
                
            
            for(int i=0;i<n_c;i++)
            for(int t=0;t<n_a;t++)
            for(int v=0;v<n_a;v++){
                V2 =(J[i*n_a+t]-K[i*n_a+t])*(J[i*n_a+v]-K[i*n_a+v]);
                V2+=J[i*n_a+t]*J[i*n_a+v];
                if(a!=b)
                    V2+=K[i*n_a+t]*K[i*n_a+v];
                dE=e_c[i]+e_IP[t]-e_v[a]-e_v[b];
                PH_th[nt][t*n_a+v]+=V2*dE/(dE*dE+edshift);
                
            }
        }
        
        delete[] K;
        delete[] J;
    }
    if(RI==0)delete[] IEIE;
    
    for(long j=0; j<num_threads;j++)
    #pragma omp parallel for
        for(long i=0; i<n_a*n_a;i++)
            RF_PH[i]+=PH_th[j][i];
        
    for(int i=0;i<num_threads;i++)delete[] PH_th[i];
    delete[] PH_th;
    
    
    return 0;
    
}

int PT_tensors::calc_EE_2_AAVV(){
    
//     num_threads=1;
    double **JK_th= new double * [num_threads];
    for(int i=0;i<num_threads;i++)JK_th[i]=new double[n_a*n_a*n_a*n_a];
    for(int i=0;i<num_threads;i++)set_zero_matr(JK_th[i],n_a*n_a*n_a*n_a);

    double **AB_th= new double * [num_threads];
    for(int i=0;i<num_threads;i++)AB_th[i]=new double[n_a*n_a*n_a*n_a];
    for(int i=0;i<num_threads;i++)set_zero_matr(AB_th[i],n_a*n_a*n_a*n_a);
    double * IEIE;
    int n_o  = n_c+n_a;//occupied
    
    if(RI==0){
        int n_ao = integrals->n_ao;
        double * CVEC = new double[(n_c+n_a)*n_ao];
        for(int i=0; i<n_ao ;i++){
            for(int j=0; j<n_c;j++) CVEC[i*(n_c+n_a)+j    ]=integrals->COR_CVEC[i*n_c+j];
            for(int j=0; j<n_a;j++) CVEC[i*(n_c+n_a)+j+n_c]=integrals->ACT_CVEC[i*n_a+j];
        }
        IEIE = new double[n_v*n_o*n_v*n_o];
        calc_2el_MO_INTS_XYXY(integrals->s, n_ao, IEIE, CVEC, integrals->VAC_CVEC, n_o, n_v);
        
        delete[] CVEC;
        
    }

    
    #pragma omp parallel
    {
        int nt =omp_get_thread_num();
        
        double AB,JK,dE;
        
        double * Ja;
        
        Ja = new double[n_a*n_a];
        
        for(int a=nt; a<n_v; a+=num_threads)
        for(int b=a; b<n_v; b++){//fprintf(stderr,"AAVV a,b=%3d,%3d  (thread %2d)    \r",a,b,nt);
            if(RI)
                cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            n_a,n_a,integrals->aux_n_ao,1.0,
                            integrals->VA_RI_M+a*n_a*integrals->aux_n_ao,integrals->aux_n_ao,
                            integrals->VA_RI_M+b*n_a*integrals->aux_n_ao,integrals->aux_n_ao,0.0,
                            Ja,n_a);
            else{
                for(int t=0  ;t<n_a;t++)
                for(int u=0  ;u<n_a;u++)
                    Ja[t*n_a+u]=IEIE[(((t+n_c)*n_v+a)*n_o+u+n_c)*n_v+b];
            }
            
            
            for(int t=0  ;t<n_a;t++)
            for(int u=0  ;u<n_a;u++)
            for(int v=0  ;v<n_a;v++)
            for(int w=0  ;w<n_a;w++){
                JK=(Ja[t*n_a+u]-Ja[u*n_a+t])*
                   (Ja[v*n_a+w]-Ja[w*n_a+v]);
                AB         =Ja[t*n_a+u]*Ja[v*n_a+w];
                if(a!=b)AB+=Ja[u*n_a+t]*Ja[w*n_a+v];
                   
                dE=e_IP[t]+e_IP[u]-e_v[a]-e_v[b];
                JK_th[nt][((u*n_a+t)*n_a+v)*n_a+w]+=JK*dE/(dE*dE+edshift);///inverse t-u in res_fit_calc_2body_AA (?)
                AB_th[nt][((t*n_a+u)*n_a+v)*n_a+w]+=AB*dE/(dE*dE+edshift);
                
                
            }
                
            
        }
    
       delete[] Ja;
    
    
    }
    if(RI==0)delete[] IEIE;
    
//     #pragma omp parallel for
    for(long j=0; j<num_threads;j++)
        for(long i=0; i<n_a*n_a*n_a*n_a;i++)
            RF_PV_JK[i]+=JK_th[j][i];
//     #pragma omp parallel for
    for(long j=0; j<num_threads;j++)
        for(long i=0; i<n_a*n_a*n_a*n_a;i++)
            RF_PV_AB[i]+=AB_th[j][i];
        
    for(int i=0;i<num_threads;i++)delete[] JK_th[i];
    delete[] JK_th;
    
    for(int i=0;i<num_threads;i++)delete[] AB_th[i];
    delete[] AB_th;
    return 0;
}

int PT_tensors::calc_EE_2_CCAV(){
    
    int n_o=n_c+n_a;//occupied
    
    
    double **RF_PH_th= new double * [num_threads];
    for(int i=0;i<num_threads;i++)RF_PH_th[i]   = new double[n_a*n_a];
    for(int i=0;i<num_threads;i++)set_zero_matr(RF_PH_th[i], n_a*n_a);
    
    
    double * IIIE;
    if(RI==0){
        IIIE = new double[n_v*n_o*n_o*n_o];
        
        int n_ao = integrals->n_ao;
        
        double * CVEC = new double[(n_c+n_a)*n_ao];
        for(int i=0; i<n_ao ;i++){
            for(int j=0; j<n_c;j++) CVEC[i*(n_c+n_a)+j    ]=integrals->COR_CVEC[i*n_c+j];
            for(int j=0; j<n_a;j++) CVEC[i*(n_c+n_a)+j+n_c]=integrals->ACT_CVEC[i*n_a+j];
        }
        calc_2el_MO_INTS_XXXY(integrals->s, integrals->n_ao, IIIE, CVEC, integrals->VAC_CVEC, n_o, n_v);
        delete[] CVEC;
        
    }
    
    
    
    #pragma omp parallel
    {
        int nt =omp_get_thread_num();
        
        double *K;
        K= new double[n_a*n_v];
        double *J;
        J= new double[n_a*n_v];
        
        double dE;
        double V2;
        
        for(int i=nt; i<n_c; i+=num_threads)
        for(int j=i; j<n_c; j++){//fprintf(stderr," CCAV i,j=%d,%d\r",i,j);
            if(RI){
                for(int a=0; a<n_v; a++)
                for(int t=0; t<n_a; t++){
                    J[a*n_a+t]=cblas_ddot(integrals->aux_n_ao, integrals->VC_RI_M+(a*n_c+i)*integrals->aux_n_ao, 1, integrals->CA_RI_M+(j*n_a+t)*integrals->aux_n_ao, 1);
                    K[a*n_a+t]=cblas_ddot(integrals->aux_n_ao, integrals->VC_RI_M+(a*n_c+j)*integrals->aux_n_ao, 1, integrals->CA_RI_M+(i*n_a+t)*integrals->aux_n_ao, 1);
                }
            }
            else{
                for(int a=0; a<n_v; a++)
                for(int t=0; t<n_a; t++){
                    J[a*n_a+t]=IIIE[((i*n_o+t+n_c)*n_o+j)*n_v+a];
                    K[a*n_a+t]=IIIE[((j*n_o+t+n_c)*n_o+i)*n_v+a];
                }
            }
                
            
            for(int a=0;a<n_v;a++)
            for(int t=0;t<n_a;t++)
            for(int v=0;v<n_a;v++){
                V2 =(J[a*n_a+t]-K[a*n_a+t])*(J[a*n_a+v]-K[a*n_a+v]);
                V2+=J[a*n_a+t]*J[a*n_a+v];    //i<j is filled by J
                if(i!=j)                      //i=j is filled only once J[i,i]=K[i,i]
                    V2+=K[a*n_a+t]*K[a*n_a+v];//i>j is filled by K
                dE=e_c[i]+e_c[j]-e_EA[t]-e_v[a];
                RF_PH_th[nt][v*n_a+t]-=V2*dE/(dE*dE+edshift);//sign changed!!!!
                
            }
        }
        
        delete[] K;
        delete[] J;
    }
    if(RI==0)delete[] IIIE;
    
    for(long j=1; j<num_threads;j++)
//     #pragma omp parallel for
        for(long i=0; i<n_a*n_a;i++)
            RF_PH_th[0][i]+=RF_PH_th[j][i];
    
//     #pragma omp parallel for
    for(int t=0; t<n_a; t++)
        RF_PS-=RF_PH_th[0][t*n_a+t]*2;
    
    #pragma omp parallel for // parallelism is bad but i don't care
    for(int w=0; w<n_a; w++)
    for(int t=0; t<n_a; t++)
        RF_PH[w*n_a+t]+=RF_PH_th[0][w*n_a+t];
    
    for(int i=0;i<num_threads;i++)delete[] RF_PH_th[i];
    delete[] RF_PH_th;
    
    
    return 0;
    
}

int PT_tensors::calc_EE_2_CCAA(){
    
//     num_threads=1;
    
    double ** RF_PV_AB_th = new double * [num_threads];
    double ** RF_PV_JK_th = new double * [num_threads];
    
    for(int i=0;i<num_threads;i++)RF_PV_AB_th[i] = new double[n_a*n_a*n_a*n_a];
    for(int i=0;i<num_threads;i++)RF_PV_JK_th[i] = new double[n_a*n_a*n_a*n_a];
    
    for(int i=0;i<num_threads;i++)set_zero_matr(RF_PV_AB_th[i], n_a*n_a*n_a*n_a);
    for(int i=0;i<num_threads;i++)set_zero_matr(RF_PV_JK_th[i], n_a*n_a*n_a*n_a);
    
    double * ACAC;
    
    if(RI==0){
        int n_ao = integrals->n_ao;
        ACAC = new double[n_a*n_c*n_a*n_c];
        calc_2el_MO_INTS_XYXY(integrals->s, n_ao, ACAC, integrals->ACT_CVEC, integrals->COR_CVEC, n_a, n_c);
    }
    
    
    
    
    #pragma omp parallel
    {
        int nt =omp_get_thread_num();
    
        double dE;
        double Ed;
        
        double V2;
        double AB;
        
        double * Ja;
        Ja = new double[n_a*n_a];
        
        
        for(int i=nt; i<n_c; i+=num_threads)
        for(int j=i; j<n_c; j++){//fprintf(stderr,"CCAA i,j=%3d,%3d  (thread %2d)    \r",i,j,nt);
            
            if(RI)
                cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                            n_a,n_a,integrals->aux_n_ao,1.0,
                            integrals->CA_RI_M+i*n_a*integrals->aux_n_ao,integrals->aux_n_ao,
                            integrals->CA_RI_M+j*n_a*integrals->aux_n_ao,integrals->aux_n_ao,0.0,
                            Ja,n_a);
            else{
                for(int t=0  ;t<n_a;t++)
                for(int u=0  ;u<n_a;u++)
                    Ja[t*n_a+u]=ACAC[((t*n_c+i)*n_a+u)*n_c+j];
            }
            
            
            for(int t=0  ;t<n_a;t++)
            for(int u=0  ;u<n_a;u++)
            for(int v=0  ;v<n_a;v++)
            for(int w=0  ;w<n_a;w++){
                
                V2=(Ja[t*n_a+u]-Ja[u*n_a+t])*
                   (Ja[v*n_a+w]-Ja[w*n_a+v]);
                AB         =Ja[t*n_a+u]*Ja[v*n_a+w];
                if(i!=j)AB+=Ja[u*n_a+t]*Ja[w*n_a+v];
                
                dE=-e_EA[v]-e_EA[w]+e_c[i]+e_c[j];
                
                RF_PV_JK_th[nt][((u*n_a+t)*n_a+v)*n_a+w]+=V2*dE/(dE*dE+edshift);///inverse t-u in res_fit_calc_2body_AA (?)
                RF_PV_AB_th[nt][((t*n_a+u)*n_a+v)*n_a+w]+=AB*dE/(dE*dE+edshift);
                
            }
        }
        
        delete[] Ja;
    }
    if(RI==0)delete[] ACAC;
    
    for(long j=1; j<num_threads;j++)
//     #pragma omp parallel for
        for(long i=0; i<n_a*n_a*n_a*n_a;i++)
            RF_PV_JK_th[0][i]+=RF_PV_JK_th[j][i];

    for(long j=1; j<num_threads;j++)
//     #pragma omp parallel for
        for(long i=0; i<n_a*n_a*n_a*n_a;i++)
            RF_PV_AB_th[0][i]+=RF_PV_AB_th[j][i];

    
    #pragma omp parallel for // parallelism is bad but i don't care
    for(long i=0;i<n_a*n_a*n_a*n_a;i++)RF_PV_JK[i]+=RF_PV_JK_th[0][i];
    
    #pragma omp parallel for // parallelism is bad but i don't care
    for(long i=0;i<n_a*n_a*n_a*n_a;i++)RF_PV_AB[i]+=RF_PV_AB_th[0][i];
    
    #pragma omp parallel for // parallelism is bad but i don't care
    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++)
    for(int t=0; t<n_a; t++)
        RF_PH[v*n_a+u]+=RF_PV_JK_th[0][((t*n_a+v)*n_a+t)*n_a+u]*0.25;
    
    #pragma omp parallel for // parallelism is bad but i don't care
    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++)
    for(int t=0; t<n_a; t++)
        RF_PH[v*n_a+t]-=RF_PV_JK_th[0][((u*n_a+v)*n_a+t)*n_a+u]*0.25;

    #pragma omp parallel for // parallelism is bad but i don't care
    for(int w=0; w<n_a; w++)
    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++)
        RF_PH[w*n_a+u]-=RF_PV_JK_th[0][((w*n_a+v)*n_a+v)*n_a+u]*0.25;

    #pragma omp parallel for // parallelism is bad but i don't care
    for(int w=0; w<n_a; w++)
    for(int u=0; u<n_a; u++)
    for(int t=0; t<n_a; t++)
        RF_PH[w*n_a+t]+=RF_PV_JK_th[0][((w*n_a+u)*n_a+t)*n_a+u]*0.25;

    #pragma omp parallel for // parallelism is bad but i don't care
    for(int w=0; w<n_a; w++)
    for(int u=0; u<n_a; u++)
    for(int t=0; t<n_a; t++)
        RF_PH[w*n_a+u]-=RF_PV_AB_th[0][((t*n_a+w)*n_a+t)*n_a+u];



    for(int v=0; v<n_a; v++)
    for(int u=0; u<n_a; u++)
        RF_PS+=RF_PV_JK_th[0][((u*n_a+v)*n_a+v)*n_a+u]*0.5;//A+B

    for(int v=0; v<n_a; v++)
    for(int w=0; w<n_a; w++)
        RF_PS-=RF_PV_JK_th[0][((w*n_a+v)*n_a+w)*n_a+v]*0.5;//A+B
    
    for(int u=0; u<n_a; u++)
    for(int t=0; t<n_a; t++)
        RF_PS+=RF_PV_AB_th[0][((t*n_a+u)*n_a+t)*n_a+u];

    
    
    
    for(int i=0;i<num_threads;i++)delete[] RF_PV_JK_th[i];
    for(int i=0;i<num_threads;i++)delete[] RF_PV_AB_th[i];
    
    
    delete[] RF_PV_JK_th;
    delete[] RF_PV_AB_th;
    
    return 0;
}

// int PT_tensors::calc_EE_2_VVVV(){
//     
//     
//     
//     int aux_n_ao = integrals->aux_n_ao;
//     double * VV=integrals->AA_RI_M;
//     
//     
//     double * VVVV;
//     VAAA = new double[n_v*n_a*n_a*n_a];
// //     if(RI)
//         cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
//                         n_v*n_a,n_a*n_a,aux_n_ao,1.0,
//                         VV,aux_n_ao,
//                         VV,aux_n_ao,0.0,
//                         VVVV,n_a*n_a);
// //     else{
// //         double * AAAV = new double[n_v*n_a*n_a*n_a];
// //         calc_2el_MO_INTS_XXXY(integrals->s, integrals->n_ao, AAAV, integrals->ACT_CVEC, integrals->VAC_CVEC, n_a, n_v);
// //         for(int a=0;a<n_v;a++)
// //         for(int t=0;t<n_a;t++)
// //         for(int u=0;u<n_a;u++)
// //         for(int v=0;v<n_a;v++)
// //             VAAA[((a*n_a+t)*n_a+u)*n_a+v]=AAAV[((u*n_a+v)*n_a+t)*n_v+a];
// //         delete[] AAAV;
// //     }
// // 
//     
//     
//     
//     double * RF_P4_ss = new double[n_a*n_a*n_a*n_a*n_a*n_a*n_a*n_a];
//     double * RF_P4_os = new double[n_a*n_a*n_a*n_a*n_a*n_a*n_a*n_a];
//     double * RF_P4_JH = new double[n_a*n_a*n_a*n_a*n_a*n_a        ];
//     double * RF_P4_MH = new double[n_a*n_a*n_a*n_a*n_a*n_a        ];
//     double * RF_P4_HJ = new double[n_a*n_a*n_a*n_a*n_a*n_a        ];
//     double * RF_P4_HM = new double[n_a*n_a*n_a*n_a*n_a*n_a        ];
//     double * RF_P4_HH = new double[n_a*n_a*n_a*n_a                ];
//     
//     set_zero_matr(RF_P3_ss, n_a*n_a*n_a*n_a*n_a*n_a*n_a*n_a);
//     set_zero_matr(RF_P3_os, n_a*n_a*n_a*n_a*n_a*n_a*n_a*n_a);
//     set_zero_matr(RF_P3_JH, n_a*n_a*n_a*n_a*n_a*n_a        );
//     set_zero_matr(RF_P3_MH, n_a*n_a*n_a*n_a*n_a*n_a        );
//     set_zero_matr(RF_P3_HJ, n_a*n_a*n_a*n_a*n_a*n_a        );
//     set_zero_matr(RF_P3_HM, n_a*n_a*n_a*n_a*n_a*n_a        );
//     set_zero_matr(RF_P3_HH, n_a*n_a*n_a*n_a                );
//     
//     
//     
//     //VV same spin
//     for(int v=0; v<n_a; v++)
//     for(int u=0; u<n_a; u++)
//     for(int w=0; w<n_a; w++)
//     for(int t=0; t<n_a; t++){
//         int m=0
//         dE= e_IP[t]+e_IP[u]-e_EA[v]-e_EA[w];
//         if(t==v){dE-= e_IP[t]-e_EA[v];m++}
//         if(u==v){dE-= e_IP[u]-e_EA[v];m++}
//         if(t==w){dE-= e_IP[t]-e_EA[w];m++}
//         if(u==w){dE-= e_IP[u]-e_EA[w];m++}
//         if(m>1) continue;
//         
//         
//         dE=dE/(dE*dE+edshift);
//         
//         for(int v1=0; v1<n_a; v1++)
//         for(int u1=0; u1<n_a; u1++)
//         for(int w1=0; w1<n_a; w1++)
//         for(int t1=0; t1<n_a; t1++){
//         
//             V=(VVVV[((t *n_a+v )*n_a+u )*n_a+w ]*
//                VVVV[((t1*n_a+v1)*n_a+u1)*n_a+w1]);
//             RF_P4_ss[((((t*n_a+u)*n_a+v)*n_a+w)*n_a+t1)*n_a+u1)*n_a+v1)*n_a+w]+=V*dE;
//             
//         }
//     }
//     //VV - oposite spin
//     for(int v=0; v<n_a; v++)
//     for(int u=0; u<n_a; u++)
//     for(int i=0; i<n_c; i++)
//     for(int t=0; t<n_a; t++){
//         dE= e_c[i]+e_IP[u]-e_EA[t]-e_EA[v];
//         if(t==u)dE= e_c[i]-e_EA[v];
//         if(v==u)dE= e_c[i]-e_EA[t];
//         dE=dE/(dE*dE+edshift);
//         
//         for(int y=0; y<n_a; y++)
//         for(int x=0; x<n_a; x++)
//         for(int w=0; w<n_a; w++){
//             // i A
//             // v A
//             // u B
//             // t B
//             // w any
//             // x any
//             // y any
//             V=(CAAA[((i*n_a+v)*n_a+u)*n_a+t]*
//                CAAA[((i*n_a+w)*n_a+x)*n_a+y]);
//             RF_P3_os[((((t*n_a+u)*n_a+v)*n_a+y)*n_a+x)*n_a+w]+=V*dE;
//             
//         }
//     }
//     
//     
//     //MH
//     for(int v=0; v<n_a; v++)
//     for(int u=0; u<n_a; u++)
//     for(int i=0; i<n_c; i++)
//     for(int t=0; t<n_a; t++){
//         dE= e_c[i]+e_IP[u]-e_EA[t]-e_EA[v];
//         if(t==u)dE= e_c[i]-e_EA[v];
//         if(v==u)dE= e_c[i]-e_EA[t];
//         dE=dE/(dE*dE+edshift);
//         for(int w=0; w<n_a; w++){
//             H2a=H_CA[i*n_a+w];
//             V=(CAAA[((i*n_a+v)*n_a+u)*n_a+t]
//               -CAAA[((i*n_a+u)*n_a+v)*n_a+t])*
//                H2a;
//             RF_P3_MH[((t*n_a+u)*n_a+v)*n_a+w]+=V*dE;
//         }
//     }
//     //JH
//     for(int v=0; v<n_a; v++)
//     for(int u=0; u<n_a; u++)
//     for(int i=0; i<n_c; i++)
//     for(int t=0; t<n_a; t++){
//         dE= e_c[i]+e_IP[u]-e_EA[t]-e_EA[v];
//         if(t==u)dE= e_c[i]-e_EA[v];
//         if(v==u)dE= e_c[i]-e_EA[t];
//         dE=dE/(dE*dE+edshift);
//         for(int w=0; w<n_a; w++){
//             H2a=H_CA[i*n_a+w];
//             V=(CAAA[((i*n_a+v)*n_a+u)*n_a+t]*
//                H2a);
//             RF_P3_JH[((t*n_a+u)*n_a+v)*n_a+w]+=V*dE;
//         }
//     }
// 
//     //HM
//     for(int i=0; i<n_c; i++)
//     for(int t=0; t<n_a; t++){
//         dE= e_c[i]-e_EA[t];
//         dE=dE/(dE*dE+edshift);
//         
//         for(int y=0; y<n_a; y++)
//         for(int x=0; x<n_a; x++)
//         for(int w=0; w<n_a; w++){
//             H1a=H_CA[i*n_a+t];
//             V= H1a*
//               (CAAA[((i*n_a+w)*n_a+x)*n_a+y]
//               -CAAA[((i*n_a+x)*n_a+w)*n_a+y]);
//             RF_P3_HM[((t*n_a+y)*n_a+x)*n_a+w]+=V*dE;
//         }
//     }
//     //HJ
//     for(int i=0; i<n_c; i++)
//     for(int t=0; t<n_a; t++){
//         dE= e_c[i]-e_EA[t];
//         dE=dE/(dE*dE+edshift);
//         
//         for(int y=0; y<n_a; y++)
//         for(int x=0; x<n_a; x++)
//         for(int w=0; w<n_a; w++){
//             H1a=H_CA[i*n_a+t];
//             V=(H1a*
//                CAAA[((i*n_a+w)*n_a+x)*n_a+y]);
//             RF_P3_HJ[((t*n_a+y)*n_a+x)*n_a+w]+=V*dE;
//         }
//     }
//     
//     //HH
//     for(int i=0; i<n_c; i++)
//     for(int t=0; t<n_a; t++){
//         dE= e_c[i]-e_EA[t];
//         dE=dE/(dE*dE+edshift);
//         
//         for(int u=0; u<n_a; u++){
//             H1a=H_CA[i*n_a+t];
//             H2a=H_CA[i*n_a+u];
//             V=H1a*H2a;
//             RF_P3_HH[t*n_a+u]+=V*dE;
//         }
//     }
//     
//     return 0;
//     
// }


int PT_tensors::calc_EE_1_AV(double* P_AV) {    //// AV block written by G.A.

    double dE;          
    double V;           

    double* VAAA;
    VAAA = new double[n_v * n_a * n_a * n_a];

    integrals->VAAA_calc(VAAA);                                  //// (at|uv) or (av|tu)

    set_zero_matr(RF_PV, n_a * n_a * n_a * n_a);

    for (int a = 0; a < n_v; a++)
        for (int t = 0; t < n_a; t++)
            for (int u = 0; u < n_a; u++)
                for (int v = 0; v < n_a; v++) {
                    dE = e_IP[t] + e_IP[u] - e_v[a] - e_EA[v];
                    dE = dE / (dE * dE + edshift);

                    for (int w = 0; w < n_a; w++) {
                        //<tu|va>*P[w,a];<tu|va>=(au|tv)
                        V = VAAA[((a * n_a + u) * n_a + t) * n_a + v] * P_AV[w * n_v + a];  //// 6th term
                        RF_PV[((t * n_a + u) * n_a + v) * n_a + w] += V * dE;

                    }
                }

    // #pragma omp parallel for
    for (int a = 0; a < n_v; a++)
        for (int t = 0; t < n_a; t++) {
            dE = e_IP[t] - e_v[a];
            dE = dE / (dE * dE + edshift);

            for (int w = 0; w < n_a; w++) {
                //H[t,a]*P[w,a]
                V = H_AV[t * n_v + a] * P_AV[w * n_v + a];    //// 4th   term
                RF_PH[t * n_a + w] += V * dE;

            }
        }
    //     }

    sym_JK_T(RF_PV_JK, RF_PV, n_a);     

    sym_AB_T(RF_PV_AB, RF_PV, n_a);     

    delete[] VAAA;

    return 0;
}


int PT_tensors::calc_EE_1_CA(double* P_CA) {                  //// CA block written by G.A.

    double dE;          
    double V;           

    double* CAAA;       ////  CAAA block <vw||ui>=<wv||iu>=<wv|iu>-<wv|ui>=(wi|vu)-(wu|vi)=(iw|uv)-(uw|iv)= (iw|uv)-(iv|uw)

    CAAA = new double[n_c * n_a * n_a * n_a];
    integrals->CAAA_calc(CAAA);                                 

    //// H1a = H_CA;    //// get externally F_core <i|F_core,j|t>
    double* RF_PH_tmp = new double[n_a * n_a];

    set_zero_matr(RF_PV, n_a * n_a * n_a * n_a);              
    set_zero_matr(RF_PH_tmp, n_a * n_a);


    for (int i = 0; i < n_c; i++)                                      //// i in C
        for (int t = 0; t < n_a; t++)                                  //// t in A 
            for (int u = 0; u < n_a; u++)                              //// u in A
                for (int v = 0; v < n_a; v++) {                        //// v in A
                    dE = e_c[i] + e_IP[t] - e_EA[u] - e_EA[v];         //// denominator OK
                    dE = dE / (dE * dE + edshift);

                    for (int w = 0; w < n_a; w++) {                    //// run over all active el -> last term -> rename ind as V_{tuvw}
                        //<it|uv>*P[w,a];<it|uv>=(iu|tv)               //// V_{tuvw} <m|d t+ u+ w v>, <i|d|t>*<vw||iu>=<i|d|t>[]
                        V = CAAA[((i * n_a + u) * n_a + t) * n_a + v] * P_CA[i * n_a + w];   
                        RF_PV[((t * n_a + u) * n_a + v) * n_a + w] += V * dE;  ////  term 9, need to change t and w, also u and v change positions

                    }
                }
    for (int i = 0; i < n_c; i++)
        for (int t = 0; t < n_a; t++) {
            dE = e_c[i] - e_EA[t];
            dE = dE / (dE * dE + edshift);

            for (int w = 0; w < n_a; w++) {                              //// sum over all Act, not only in vac
                //H[t,a]*P[w,a]
                V = H_CA[i * n_a + t] * P_CA[i * n_a + w];
                RF_PH_tmp[w * n_a + t] -= V * dE;                         //// term 7

            }
        }

    for (int t = 0; t < n_a; t++)
        RF_PS -= RF_PH_tmp[t * n_a + t] * 2;          //// second time minus sign, 1st term - to constant from diag part of F_core

    for (int t = 0; t < n_a; t++)
        for (int w = 0; w < n_a; w++)
            RF_PH[w * n_a + t] += RF_PH_tmp[w * n_a + t];   //// term 7 to D_{tu}

    for (int w = 0; w < n_a; w++)
        for (int t = 0; t < n_a; t++)
            for (int v = 0; v < n_a; v++)
                RF_PH[t * n_a + w] += 2 * RF_PV[((t * n_a + v) * n_a + w) * n_a + v]
                - RF_PV[((t * n_a + w) * n_a + v) * n_a + v];                             //// term 5

    sym_JK_CA_T(RF_PV_JK, RF_PV, n_a);         //// see PT_tensors_IPEA, for t=w or t=x and J-K, for 9th term

    sym_AB_CA_T(RF_PV_AB, RF_PV, n_a);         //// 9th term

    delete[] CAAA;
    delete[] RF_PH_tmp;

    return 0;
}

int PT_tensors::calc_EE_1_CV(double* P_CV) { //// CV block written by G.A.

    double V;
    double dE;         //// <au|it>=<au|it>-<au|ti>=(ai|ut)-(at|ui)=[real orbs]=??? (ia|tu)-(ui|at)=(ia|tu)-(ui|ta)

    double* VCAA;      //// (VC|AA)=(ai|tu)
    double* CAVA;      //// (CA|VA)=(iu|at)

    VCAA = new double[n_c * n_v * n_a * n_a]; 
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
        n_c * n_v, n_a * n_a, integrals->aux_n_ao, 1.0,
        integrals->VC_RI_M, integrals->aux_n_ao,
        integrals->AA_RI_M, integrals->aux_n_ao, 0.0,
        VCAA, n_a * n_a);                     //// multiplication of (ia|K) w/ (K|tu)
    CAVA = new double[n_a * n_c * n_a * n_v];
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
        n_a * n_c, n_a * n_v, integrals->aux_n_ao, 1.0,
        integrals->CA_RI_M, integrals->aux_n_ao,
        integrals->VA_RI_M, integrals->aux_n_ao, 0.0,
        CAVA, n_a * n_v);                     //// mult (ui|K) by (K|ta)

    for (int i = 0; i < n_c; i++)                                             //// i in C  
        for (int a = 0; a < n_v; a++)                                         //// a in V   //// first pair of indices defines block
            for (int t = 0; t < n_a; t++)                                     //// t in A
                for (int u = 0; u < n_a; u++) {                               //// u in A
                    //(<it|au>-<iu|at>)*P[i,a]                                //// new notation P[i,a]=<i|d|a>, P[i,a]*[(ia|tu)-(ui|ta)]
                    
                    V = (2 * VCAA[((a * n_c + i) * n_a + u) * n_a + t] -  //// !!! order of ind
                        CAVA[((i * n_a + u) * n_v + a) * n_a + t]) *      ////  2J[]-K[]  over space orbs when same-spin + opp-spin
                        P_CV[i * n_v + a];
                                        
                    dE = e_c[i] + e_IP[t] - e_v[a] - e_EA[u];             //// denominator is e_{t}-e_{u}+e_{i}-e_{a}
                    RF_PH[t * n_a + u] += V * dE / (dE * dE + edshift);   //// term (3)

                }

    for (int i = 0; i < n_c; i++)                                         //// i in C
        for (int a = 0; a < n_v; a++) {                                   //// a in V
            V = H_CV[i * n_v + a] * P_CV[i * n_v + a] * 2;                
            dE = e_c[i] - e_v[a];                                         
            RF_PS += V * dE / (dE * dE + edshift);

        }

    delete[] VCAA;
    delete[] CAVA;

    return 0;
}

int PT_tensors::P1_calc_EE(double* P1, aldet_data* CI, double* P_AV, double* P_CA, double* P_CV, int AV, int CA, int CV) { //// my take

         int n_s;   
         
         n_s = CI->n_states[0];
         
         aldet_data first_order_props;
         first_order_props.get_dim(CI->n_act, CI->na, CI->nb, 1, CI->mult, CI->print_number);
         first_order_props.copy_coef(0,CI,CI->n_states[0],0,1);

         RF_PS=0;
     
         set_zero_matr(RF_PH      , n_a*n_a        );
         set_zero_matr(RF_PV_JK   , n_a*n_a*n_a*n_a);
         set_zero_matr(RF_PV_AB   , n_a*n_a*n_a*n_a);
         
         //calculation of RF_PV_AB, RF_PV_JK and RF_PH
     
         if(CA)calc_EE_1_CA(P_CA); // CA
         if(CV)calc_EE_1_CV(P_CV); // CV
         if(AV)calc_EE_1_AV(P_AV); // AV
         
         
         printf_timer("calculation of 1st order correction of one-electron property (dipole moment d(1))");fflush(out_stream);
         
         first_order_props.simple_import_data(RF_PV_JK, RF_PV_AB, RF_PH, RF_PS);                //// ???
         first_order_props.H_calc(P1,CI->n_states[0]);

         char* timer_words = new char[256];
         if(CV==1)if(CA==1)if(AV==1)sprintf(timer_words, "calculation of first order properties\0");                      //// ??? for checking parts??
         if(CV==0)if(CA==1)if(AV==1)sprintf(timer_words, "calculation of CA and AV parts of first order properties\0");
         if(CV==1)if(CA==0)if(AV==1)sprintf(timer_words, "calculation of CV and AV parts of first order properties\0");
         if(CV==1)if(CA==1)if(AV==0)sprintf(timer_words, "calculation of CV and CA parts of first order properties\0");
         if(CV==1)if(CA==0)if(AV==0)sprintf(timer_words, "calculation of CV part of first order properties\0");
         if(CV==0)if(CA==1)if(AV==0)sprintf(timer_words, "calculation of CA part of first order properties\0");
         if(CV==0)if(CA==0)if(AV==1)sprintf(timer_words, "calculation of AV part of first order properties\0");
         if(CV==0)if(CA==0)if(AV==0)sprintf(timer_words, "calculation of first order properties\0");
         
         printf_timer(timer_words);fflush(out_stream);
         delete[] timer_words;
         
    return 0;

}

