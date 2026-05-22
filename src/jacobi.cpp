//C
#include "math.h"
//NOPT
#include "jacobi.h"
#include "matr.h"
#include "blas_link.h"

jacobi_mcscf_sd_engine::jacobi_mcscf_sd_engine(){
    
    rot_matr=nullptr;
    
}

int jacobi_mcscf_sd_engine::init(double * ext_G, double * ext_B, int ext_n_c, int ext_n_a, int ext_n_v, int ext_n_ao, int ext_n_s, double ext_x_max_jsd){
    
    n_c  = ext_n_c;
    n_a  = ext_n_a;
    n_v  = ext_n_v;
    n_s  = ext_n_s;
    n_ao = ext_n_ao;
    n_mo = n_c+n_a+n_v;
    dim  = n_c*n_a+n_a*n_v+n_c*n_v;
    
    x_max_jsd=ext_x_max_jsd;
    
    G=ext_G;
    
    B=ext_B;
    
    rot_matr=new double[n_mo*n_mo];
    
    return 0;
}

double jacobi_mcscf_sd_engine::find_max_el(){
    
    double max=0;
    
    for(int i=0;i<dim*n_s;i++)if(fabs(G[i])>max)max=fabs(G[i]);
    
    return max;
}

double rot_matr_scale(double * R, int n, double max){
    
    double x_max=0;
    for(int i=0  ; i<n; i++)
    for(int j=i+1; j<n; j++)
        if(fabs(R[i*n+j])>x_max)x_max=fabs(R[i*n+j]);
        
    if(x_max>max){
        double c= max/x_max;
        for(int i=0  ; i<n; i++)
        for(int j=0; j<n; j++)
            if(i!=j)R[i*n+j]=R[i*n+j]*c;
            
    }
    
    return x_max;
    
}

double jacobi_mcscf_sd_engine::step(double * VEC, double * BUF){
    
    //generation of rot matr
    set_zero_matr(rot_matr, n_mo*n_mo);
    for(int i=0;i<n_mo;i++) rot_matr[i*(n_mo+1)]=1.0;// 1 at diagonal
    
    //CA
    for(int i=0;i<n_c;i++)
    for(int t=0;t<n_a;t++){
//         printf("i=%d, t=%d\n",i,t);
        rot_matr[i*n_mo+ n_c+t ]= rotation(i*n_a+t,i,t+n_c); //upper corner
        rot_matr[i+n_mo*(n_c+t)]=-rot_matr[i*n_mo+ n_c+t ]; //lower corner
    }
    //CV
    for(int i=0;i<n_c;i++)
    for(int a=0;a<n_v;a++){
//         printf("i=%d, a=%d\n",i,a);
        rot_matr[i*n_mo+ n_c+n_a+a ]= rotation(i*n_v+a+n_c*n_a,i,a+n_c+n_a); //upper corner
        rot_matr[i+n_mo*(n_c+n_a+a)]=-rot_matr[i*n_mo+ n_c+n_a+a ]; //lower corner
    }
    //AV
    for(int t=0;t<n_a;t++)
    for(int a=0;a<n_v;a++){
//         printf("t=%d, a=%d\n",t,a);
        rot_matr[(n_c+t)*n_mo+ n_c+n_a+a ]= rotation(t*n_v+a+n_c*n_a+n_c*n_v,t+n_c,a+n_c+n_a); //upper corner
        rot_matr[(n_c+t)+n_mo*(n_c+n_a+a)]=-rot_matr[(n_c+t)*n_mo+ n_c+n_a+a ]; //lower corner
    }
    
//     fprintf(out_stream,"R\n");
//     PrintMatr(rot_matr,n_mo,n_mo,1);
    
    int scaling=0;
    
    double x_max = rot_matr_scale(rot_matr, n_mo, x_max_jsd) ;
    
    if(x_max > x_max_jsd){
        fprintf(out_stream,"SOSCF is scaling rotation angle matrix Xmax=%f     |\n",x_max);
        scaling=1;
    }
    
    ortogonalization(rot_matr, n_mo);//ortogonalization of rot matrix
    //rotation of vectors
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,n_mo,n_ao,n_mo,1.0, rot_matr  ,n_mo,VEC,n_ao,0.0,BUF,n_ao);//TMP = C * U^T
    for(int i=0;i<n_ao*n_ao;i++)VEC[i]=BUF[i];
    
    return x_max;
    
}

double jacobi_mcscf_sd_engine::rotation(int ab, int a, int b){
    
//     double c;
    double a2 =0;
    double b2 =0; 
    double aXb=0;
    
    for(int i_s=0; i_s<n_s; i_s++)
    {
        double g1=B[ab+i_s*dim];//B[b*n_ao+b+i_s*n_ao*n_ao] - B[a*n_ao+a+i_s*n_ao*n_ao];
        double g2=2*G[ab+i_s*dim];

        a2  += g1*g1;
        aXb += g1*g2;
        b2  += g2*g2;
//         printf("\n%e %e\n",B[ab+i_s*dim], G[ab+i_s*dim]);
//         printf("%e %e\n\n",G[ab+i_s*dim], 0.0);
//         getchar();
    }
    double t_on  = a2 - b2;
    double t_off = 2.0* aXb;
    
//     printf("of - %e on - %e\n",t_off,t_on + sqrt( t_on*t_on + t_off * t_off));
    
    double theta = 0.5* std::atan2(t_off,t_on + sqrt( t_on*t_on + t_off * t_off));
    double c = cos(theta);
    double s = sin(theta);
    
//     printf("c=%e s=%e\n",c,s);
//     if(fabs(s)>0.1)getchar();
    
//     D[i*dimM+j]+= G[0*2+0] * 2.0*s*(s*s-1.0)
//                             + G[0*2+1] * 2.0*c*(c*c-s*s)
//                             + G[1*2+1] * 4.0*c*c*s;
                
    return -s/c;
}

jacobi_mcscf_sd_engine::~jacobi_mcscf_sd_engine(){
    
//     if(B       !=nullptr)delete[] B       ;
    if(rot_matr!=nullptr)delete[] rot_matr;
    
}


