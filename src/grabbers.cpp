#include <cstdio>
#include <vector>
#include <cstdlib>

int grab_E=0;
std::vector<double> E_grabbed;

int grab_D=0;
int grab_HS_data=0;
int grabber_n_states;
double *  H_grabbed=nullptr;
double *  S_grabbed=nullptr;
double * Dx_grabbed=nullptr;
double * Dy_grabbed=nullptr;
double * Dz_grabbed=nullptr;


FILE * e_scan=NULL;
FILE * d_scan=NULL;

int alloc_HS_grabber(int n, int m){
    
    if( H_grabbed!=nullptr)delete[] H_grabbed;
    if( S_grabbed!=nullptr)delete[] S_grabbed;
    if(Dx_grabbed!=nullptr)delete[]Dx_grabbed;
    if(Dy_grabbed!=nullptr)delete[]Dy_grabbed;
    if(Dz_grabbed!=nullptr)delete[]Dz_grabbed;
    
     H_grabbed=new double[n*m];
     S_grabbed=new double[n*m];
    Dx_grabbed=new double[n*m];
    Dy_grabbed=new double[n*m];
    Dz_grabbed=new double[n*m];
    
    grabber_n_states = m;
    
    return 0;
    
}

int alloc_D_grabber(int n, int m){
    
    if(Dx_grabbed!=nullptr)delete[]Dx_grabbed;
    if(Dy_grabbed!=nullptr)delete[]Dy_grabbed;
    if(Dz_grabbed!=nullptr)delete[]Dz_grabbed;
    
    Dx_grabbed=new double[n*m];
    Dy_grabbed=new double[n*m];
    Dz_grabbed=new double[n*m];
    
    grabber_n_states = m;
    
    return 0;
    
}

