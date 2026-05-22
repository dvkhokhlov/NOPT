# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <etc.h>
# include <math.h>
#include "common_vars.h"


double * new_double_w_check(long long A, const char * stage){
    
    double * p =new (std::nothrow) double[1LL*A];
    
    if(p==nullptr){
        if     (fabs(A)>1024LL*1024*1024*1024)printf("ERROR: NOPT could not allocate %lld TB of memory\n",A/1024/1024/1024/1024);
        else if(fabs(A)>1024LL*1024*1024     )printf("ERROR: NOPT could not allocate %lld GB of memory\n",A/1024/1024/1024     );
        else if(fabs(A)>1024LL*1024          )printf("ERROR: NOPT could not allocate %lld MB of memory\n",A/1024/1024          );
        else if(fabs(A)>1024LL               )printf("ERROR: NOPT could not allocate %lld kB of memory\n",A/1024               );
        else                                  printf("ERROR: NOPT could not allocate %lld  B of memory\n",A                    );
        printf("       for %s\n",stage);
        exit(EXIT_FAILURE);
    }
    
    return p;

}

double dist(double x1, double y1, double z1,
            double x2, double y2, double z2,
            int control){
    double r;
    r=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2);
    if(control)if(r<1E-6){
        //fprintf(out_stream,"Error: zero distance\n\n");
        return -19;
    }
    return sqrt(r);
    
}

double angle(double x1, double y1, double z1,
             double x2, double y2, double z2,
             double x3, double y3, double z3,
             int control){
    double r1;
    r1=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2);
    if(control)if(r1<1E-6){
        //fprintf(out_stream,"Error: zero distance\n\n");
        return -19;
    }
    double r2;
    r2=(x3-x2)*(x3-x2)+(y3-y2)*(y3-y2)+(z3-z2)*(z3-z2);
    if(control)if(r2<1E-6){
        //fprintf(out_stream,"Error: zero distance\n\n");
        return -19;
    }
    double sc_pr;
    sc_pr=(x1-x2)*(x3-x2)+(y1-y2)*(y3-y2)+(z1-z2)*(z3-z2);
    
    double a;
    a= sc_pr/(sqrt(r1*r2));
    
    return acos(a);
    
}

double dihedral(double x1, double y1, double z1,
                double x2, double y2, double z2,
                double x3, double y3, double z3,
                double x4, double y4, double z4,
                int control){
    
    
    double r1212=(x1 - x2)*(x1 - x2)+(y1 - y2)*(y1 - y2)+(z1 - z2)*(z1 - z2);
    double r1223=(x1 - x2)*(x2 - x3)+(y1 - y2)*(y2 - y3)+(z1 - z2)*(z2 - z3);
    double r1234=(x1 - x2)*(x3 - x4)+(y1 - y2)*(y3 - y4)+(z1 - z2)*(z3 - z4);
    double r2323=(x2 - x3)*(x2 - x3)+(y2 - y3)*(y2 - y3)+(z2 - z3)*(z2 - z3);
    double r2334=(x2 - x3)*(x3 - x4)+(y2 - y3)*(y3 - y4)+(z2 - z3)*(z3 - z4);
    double r3434=(x3 - x4)*(x3 - x4)+(y3 - y4)*(y3 - y4)+(z3 - z4)*(z3 - z4);

    if(control)if((r1212<1E-6)||(r2323<1E-6)||(r3434<1E-6)){
        ////fprintf(out_stream,"Error: zero distance\n\n");
        return -19;
    }
    
    double q=-(r2323*r1234-r1223*r2334)/sqrt(r1212*r2323-r1223*r1223)/sqrt(r3434*r2323-r2334*r2334); // cos(theta)
    
    
    double t=acos(q); // absolute value of the torsion angle
    if (q> 1.001){fprintf(out_stream,"WARNING:cos(theta)>1. cos-1=%e\n",q-1);t=0.0;}
    if (q<-1.001){fprintf(out_stream,"WARNING:cos(theta)>1. cos-1=%e\n",q+1);t=M_PI;}
    
    // parameters of the plane ax+by+cz=d in which r1,r2,r3 are:
    double a=-y2*z1 + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3;
    double b= x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3;
    double c=-x2*y1 + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3;
    
    double d=-x3*y2*z1 + x2*y3*z1 + x3*y1*z2 - x1*y3*z2 - x2*y1*z3 + x1*y2*z3;
    
    double e=a*x4+b*y4+c*z4; // r4 above the plane or below the plane
    
    if(e>d){t=fabs(t);}else{t=-fabs(t);};
    
    return t;
    
}



// int write_file_name(char * source, int b, char * target, int e_l, int p){
//     
//     while(source[b]==' ')b++;
//     int end;
//     end = b;
//     if(strstr(source,".")==NULL){printf("'%s' - doesn't contain file name\n",source);return 1;}
//     while(source[end]!='.')end++;
//     end+=(e_l+1);
//     for(int i=b;i<end;i++)target[i-b]=source[i];
//     target[end-b]='\0';
// 
//     FILE * check;
//     check=fopen(target,"r");
//     if(p)if(check!=NULL){char tmpstr[255];fgets(tmpstr,255,check);printf(" first line: %s\n",tmpstr);}
//     if(check==NULL)return 1;
//     fclose(check);
//     return 0;
//     
// }

int compare_int_ar(int * A, int * B, int n){
    
    int d=0;
    for(int i =0; i<n;i++)
        if(A[i]!=B[i])d++;
    
    return d;
    
}

int D_change(char * line){
    
    int n=0;
    line = strstr(line,"D");
    while(line!=nullptr){
        line[0]='E';
        line = strstr(line,"D");
        n++;
    }
    
    return n;
}

int cut_zero(std::vector<double>* c){
    
    int l=0;
    
    for(int i=0;i<c->size();i++)
        if(fabs((*c)[i])>1e-8)l=i;
    
    c->resize(l+1);
    
    return 1;
}

int multi_iterator::set(int * max_number_inp, int dim_inp){
    
    dim = dim_inp;
    
    number = new int[dim];
    
    max_number = new int[dim];
    
    for(int i=0;i<dim;i++) max_number[i] = max_number_inp[i];
    
    not_finished=1;
    
    return 0;
}

int multi_iterator::zero(){
    
    for(int i=0;i<dim;i++) number[i] = 0;
    
    not_finished=1;
    
    return 0;
}

int multi_iterator::step_for_next(int n){
    
    if(n<dim){
        number[n]++;
        
        if(number[n]<max_number[n]) 
            return 1;
        else{
            number[n]=0;
            return step_for_next(n+1);
        }
        
    }
    else return 0;
    
}

int multi_iterator::next(){
    
    not_finished = step_for_next(0);
    
    return 0;
}

int multi_iterator::not_ended(){
    
    return not_finished;
}

multi_iterator::~multi_iterator(){
    
    delete[] number;
    delete[] max_number;
    
}
 
