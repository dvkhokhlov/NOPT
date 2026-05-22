# include <stdio.h>
# include <math.h>
# include "l-bfgs_2_1.h"
# include "matr.h"
// # include "opt_link.h"

#ifndef LBFGS_RETURN_STEP
# define LBFGS_RETURN_X //can be changed on LBFGS_RETURN_STEP
# endif

l_bfgs_engine::l_bfgs_engine(){
    
    prev_grad_delta = NULL;
    prev_grad       = NULL;
    prev_step       = NULL;
    prev_skal_prod  = NULL;
    
}

int l_bfgs_engine::init(int n, int m){
    lbfgs_step_num =0;
    prev_grad_delta = new double[n*m];
    prev_grad = new double[n];
    for(int i=0;i<n;i++) prev_grad[i]=0;
    prev_step = new double[m*n];
    prev_skal_prod = new double[m];
    lbfgs_vec_num=m;
    return 0;
}

double l_bfgs_engine::step(double * x, double * grad, double * inv_h, int (*g_calc)(int, double *, double * ),int (*inv_h_calc)(int, double *, double * ), int n, double x_max){
    int vec_num=lbfgs_step_num;
    if (lbfgs_vec_num<lbfgs_step_num) vec_num=lbfgs_vec_num;
    int i,j;
    double * grad_transformed;
    grad_transformed = new double[n];
    double * s;
    s = new double [n];
    double * alpha;
    alpha =new double[vec_num];
    //taking grad
    g_calc(n,x,grad);
    //numgrad(calc_OF_by_par,x,grad,n);
    //calulation of gradient delta and skalar product prev_grad x pred_step
    if (lbfgs_vec_num>0)
    if (lbfgs_step_num>0){
//         printf("%d %d",(lbfgs_step_num-1)%lbfgs_vec_num,lbfgs_vec_num);
//         getchar();
        prev_skal_prod[(lbfgs_step_num-1)%lbfgs_vec_num]=0; ////CHECK IT!!!!!
        for(i=0; i<n;i++) {
            prev_grad_delta[i+(lbfgs_step_num-1)%lbfgs_vec_num*n]=grad[i]-prev_grad[i];
            prev_skal_prod[(lbfgs_step_num-1)%lbfgs_vec_num]+=prev_step[i+(lbfgs_step_num-1)%lbfgs_vec_num*n]*prev_grad_delta[i+(lbfgs_step_num-1)%lbfgs_vec_num*n];
        }
        prev_skal_prod[(lbfgs_step_num-1)%lbfgs_vec_num]=1/prev_skal_prod[(lbfgs_step_num-1)%lbfgs_vec_num];
    
    }
    //back up of gradient
    for(i=0; i<n;i++) grad_transformed[i]=grad[i];
    for(i=0; i<n;i++) prev_grad[i]=grad[i];
    
    //taking hess
#ifdef LBFGS_CALC_HESS
    inv_h_calc(n,x, inv_h);
#endif
    //first loop
    for(i=0;i<vec_num;i++){
        //alpha=rho*{prev_step x gr}
        alpha[i]=0;
        for(j=0;j<n;j++) alpha[i]+=prev_step[j+(lbfgs_step_num-i-1)%lbfgs_vec_num*n]*grad_transformed[j];
        alpha[i]=alpha[i]*prev_skal_prod[(lbfgs_step_num-i-1)%lbfgs_vec_num];
        //gr = gr - alpha * prev_grad
        for(j=0;j<n;j++)grad_transformed[j]-=alpha[i]*prev_grad_delta[j+(lbfgs_step_num-i-1)%lbfgs_vec_num*n];
    }
    //transformation throw diagonal approximate hess
    //Hess is not inverted!!!!!!
    for(i=0; i<n;i++) s[i]=grad_transformed[i]/inv_h[i];
    //second loop
    double betta;
    for(i=0;i<vec_num;i++){
        //alpha=rho*{prev_grad x s}
        betta=0;
        for(j=0;j<n;j++) betta+=prev_grad_delta[j+(lbfgs_step_num-vec_num+i)%lbfgs_vec_num*n]*s[j];
        betta=betta*prev_skal_prod[(lbfgs_step_num-vec_num+i)%lbfgs_vec_num];
        //s = s + (alpha - betta)* prev_step
        for(j=0;j<n;j++)s[j]+=(alpha[vec_num-i-1]-betta)*prev_step[j+(lbfgs_step_num-vec_num+i)%lbfgs_vec_num*n];
    }
    double zoom=-1;
    for (int i=0;i<n;i++){
      if(fabs(s[i])>x_max)if(x_max/fabs(s[i])<fabs(zoom))zoom = - x_max/fabs(s[i]);
//       printf("%d %lf\n",i,s[i]);
//       getchar();
    }
   
    for(i=0; i<n;i++){
#if defined LBFGS_RETURN_STEP
//         printf("s[%d] = %e*%e\n",i,s[i],zoom);
        x[i]=s[i]*zoom;
#elif defined LBFGS_RETURN_X
        x[i]+=s[i]*zoom;
#else
        printf("compillation option error:\nrecompile programm with\n");
        printf("define LBFGS_RETURN_STEP or LBFGS_RETURN_X\n");
        getchar();
#endif
    }
    //preparation to next step
    if (lbfgs_vec_num>0)
    for(i=0; i<n;i++){
        prev_step[i+(lbfgs_step_num)%lbfgs_vec_num*n]=s[i]*zoom;
    }
    lbfgs_step_num++;
    delete[] s;
    delete[] alpha;
    delete[] grad_transformed;
    return zoom;
}

l_bfgs_engine::~l_bfgs_engine(){
    
    if(prev_grad_delta!=NULL)delete[] prev_grad_delta;
    if(prev_grad      !=NULL)delete[] prev_grad      ;
    if(prev_step      !=NULL)delete[] prev_step      ;
    if(prev_skal_prod !=NULL)delete[] prev_skal_prod ;
    
}

bfgs_engine::bfgs_engine(){
    
    grad_delta = NULL;
    prev_grad  = NULL;
    prev_step  = NULL;
    inv_hess   = NULL;
    I_m_rsy    = NULL;
    B          = NULL;
    
//     prev_step  = NULL;
    
}

int bfgs_engine::init(int n){
    bfgs_step_num =0;
    dim=n;
    grad_delta = new double[n];
    prev_grad = new double[n];
    prev_step = new double[n];
    
    for(int i=0;i<n;i++) prev_grad[i]=0;
//     prev_step = new double[dim];
//     prev_skal_prod = new double[m];
//     lbfgs_vec_num=m;
    inv_hess= new double[n*n];
    I_m_rsy = new double[n*n];
    B       = new double[n*n];
    set_zero_matr(inv_hess,n*n);
    set_eq_diag(1.0);
    
    return 0;
}

int bfgs_engine::set_eq_diag(double c){
    
    for(int i=0;i<dim;i++)
        inv_hess[i*(dim+1)]=1.0/c;
    
    return 0;
    
}
    
int bfgs_engine::set_diag(double * c){
    
    for(int i=0;i<dim;i++)
        inv_hess[i*(dim+1)]=1.0/c[i];
    
    return 0;
    
}


double bfgs_engine::step(double * s, double * g, int (*g_calc)(int, double *, double * )){
    
    if(bfgs_step_num!=0){
        
        for(int i=0;i<dim;i++)grad_delta[i]=g[i]-prev_grad[i];
        
        double r=0;
        
        for(int i=0;i<dim;i++)r+=grad_delta[i]*prev_step[i];
        
        r=1.0/r;
        printf("%e\n",r);
//         getchar();
        
//         printf("Bk\n");
        
        for(int i=0;i<dim;i++)
        for(int j=0;j<dim;j++)
            I_m_rsy[i*dim+j]=-r*prev_step[i]*grad_delta[j];
//         PrintMatr(inv_hess,dim,dim,1);
        
        for(int i=0;i<dim;i++)
            I_m_rsy[i*(dim+1)]+=1.0;
        
//         printf("I-rsy\n");
//         PrintMatr(I_m_rsy,dim,dim,1);
        
        set_zero_matr(B,dim*dim);
        for(int i=0;i<dim;i++)
        for(int j=0;j<dim;j++)
        for(int k=0;k<dim;k++)
            B[i*dim+j]+=inv_hess[i*dim+k]*I_m_rsy[j*dim+k];
        
        set_zero_matr(inv_hess,dim*dim);
        for(int i=0;i<dim;i++)
        for(int j=0;j<dim;j++)
        for(int k=0;k<dim;k++)
            inv_hess[i*dim+j]+=I_m_rsy[i*dim+k]*B[k*dim+j];
        
//         printf("(1-r)Bk(1-r)\n");
//         PrintMatr(inv_hess,dim,dim,1);
        
        for(int i=0;i<dim;i++)
        for(int j=0;j<dim;j++)
            inv_hess[i*dim+j]+=r*prev_step[i]*prev_step[j];
        
//         printf("Bk+1\n");
//         PrintMatr(inv_hess,dim,dim,1);
    }
    
    bfgs_step_num++;
    for(int i=0;i<dim;i++){
        s[i]=0;
        for(int j=0;j<dim;j++)s[i]-=inv_hess[i*dim+j]*g[j];
    }
    
    double x_max=0.1;
    double zoom=1;
    for (int i=0;i<dim;i++){
      if(fabs(s[i])>x_max)if(x_max/fabs(s[i])<fabs(zoom))zoom = x_max/fabs(s[i]);
//       printf("%d %lf\n",i,s[i]);
//       getchar();
    }
    for (int i=0;i<dim;i++)s[i]=s[i]*zoom;
    for(int i=0;i<dim;i++)prev_grad[i]=g[i];
    for(int i=0;i<dim;i++)prev_step[i]=s[i];
        
    return 0.0;


    return 0.0;
}

bfgs_engine::~bfgs_engine(){
    
    if(grad_delta !=NULL)delete[] grad_delta ;
    if(prev_grad  !=NULL)delete[] prev_grad  ;
    if(prev_step  !=NULL)delete[] prev_step  ;
    if(inv_hess   !=NULL)delete[] inv_hess   ;
    if(I_m_rsy    !=NULL)delete[] I_m_rsy    ;
    if(B          !=NULL)delete[] B          ;
    
}


int do_nothing(int n, double * A, double * B){
    return 0;
}
