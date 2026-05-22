#ifdef _XDR_FULL
#define _XDR_READ
#endif

#ifdef _XDR_READ
#include "rpc/xdr.h"
#endif


#include "matr.h"
#include <cstring>
#include <stdlib.h>

extern FILE* out_stream;

#define PASSWORD 13772497

int get_full_name(char ** f, char * p, const char * n){
    
    int l_n = strlen(n);
    
    if(p  ==NULL){
        *f = new char[l_n+1];
        sprintf(*f,"%s\0"     ,n);
        return -1;
    }
    
    int l_p = strlen(p);
    int l_s = 1;
    
    if(p[l_p-1]=='/') l_s=0;
    if(p[l_p-1]=='_') l_s=0;
    
    *f = new char[l_p+l_s+l_n+1];
    
    if(l_s==0   ){ sprintf(*f,"%s%s\0" ,p,n); return  0;}
    if(l_s==1   ){ sprintf(*f,"%s_%s\0",p,n); return  1;}
    
    return 0;
}

#ifdef _XDR_FULL


int xdr_array_write(double * M, const char * name, char * prefix, int n, int ask){
    
    char * full_name;
    
    get_full_name(&full_name, prefix, name);
    if(ask){
        fprintf(stderr,"data will be written to %s?\npress enter to continue\n", full_name);
        getchar();
    }
    
    FILE * data;
    data = fopen(full_name,"w");
    if(data==NULL){
		fprintf(out_stream,"ERROR: could not open output file %s\n\n", full_name);
	   	exit(EXIT_FAILURE);
    };
    
    XDR out_xdr;
    xdrstdio_create (&out_xdr, data, XDR_ENCODE);
    
    for(int ii=0; ii<n;ii++) xdr_double(&out_xdr,&M[ii]);
    
    xdr_destroy(&out_xdr);
    fclose(data);
    
    
    return 0;  
    
}

int xdr_VEC_write(const char * out_name, double * MO_VEC, double * orb_energy, int n_mo, int n_ao, double * MO_VEC_B, double * orb_energy_B){
    
    FILE * data;
    data = fopen(out_name,"w");
    XDR out_xdr;
    xdrstdio_create (&out_xdr, data, XDR_ENCODE);
    
    xdr_int(&out_xdr,&n_mo);
    xdr_int(&out_xdr,&n_ao);
    
    for(int ii=0; ii<n_mo*n_ao;ii++) xdr_double(&out_xdr,&MO_VEC[ii]);
    
    for(int ii=0; ii<n_mo;ii++) xdr_double(&out_xdr,&orb_energy[ii]);
    
    int flag;
    if(MO_VEC_B!=NULL){
        flag=1;
        xdr_int(&out_xdr,&flag);
        for(int ii=0; ii<n_mo*n_ao;ii++) xdr_double(&out_xdr,&MO_VEC_B[ii]);
        
        for(int ii=0; ii<n_mo;ii++) xdr_double(&out_xdr,&orb_energy_B[ii]);
    
    }
    else{
        flag=0;
        xdr_int(&out_xdr,&flag);
    }
    
    xdr_destroy(&out_xdr);
    fclose(data);
    
    return 0;
}

#endif

#ifdef _XDR_READ
int xdr_array_read(double * M, const char * name, char * prefix, int n){
    
    FILE * data;
    char * full_name;
    
    get_full_name(&full_name, prefix, name);
//     printf("s = %s\n l = %d\n", full_name, strlen(full_name));
//     getchar();
    fprintf(out_stream,"reading %s\n",full_name);
    data = fopen(full_name,"r");
    if(data==NULL){
        fprintf(out_stream,"WARNING: could not open backup file %s for input\n",full_name);
        fprintf(out_stream,"         data will be calculated\n\n",full_name);
        return 1;
    };
    
    
//     printf("s = %s\n l = %d\n", prefix, strlen(prefix));
//     getchar();
    

    
    set_zero_matr(M,n);
    XDR out_xdr;
    xdrstdio_create (&out_xdr, data, XDR_DECODE);
    
    for(int ii=0; ii<n;ii++) xdr_double(&out_xdr,&M[ii]);

    xdr_destroy(&out_xdr);
    fclose(data);
    
    return 0;  
    
}

int xdr_VEC_read(const char * inp_name, double * MO_VEC, double * orb_energy, int n_mo, int n_ao, double * MO_VEC_B, double * orb_energy_B){
    
    FILE * data;
    data = fopen(inp_name,"r");
    if(data==NULL){
        printf("ERROR: could not open output file %s\n\n",inp_name);
        exit(EXIT_FAILURE);
    };
    

    
    XDR out_xdr;
    xdrstdio_create (&out_xdr, data, XDR_DECODE);
    
    xdr_int(&out_xdr,&n_mo);
//     printf("n_mo = %d\n",n_mo);
//     getchar();
    xdr_int(&out_xdr,&n_ao);
//     printf("n_ao = %d\n",n_ao);
//     getchar();
    for(int ii=0; ii<n_mo*n_ao;ii++) xdr_double(&out_xdr,&MO_VEC[ii]);
    
    for(int ii=0; ii<n_mo;ii++) xdr_double(&out_xdr,&orb_energy[ii]);
    
    int flag;
    xdr_int(&out_xdr,&flag);
    if(flag){
        
        MO_VEC_B=new double[n_ao*n_ao];
        set_zero_matr(MO_VEC_B,n_ao*n_ao);
        for(int ii=0; ii<n_mo*n_ao;ii++) xdr_double(&out_xdr,&MO_VEC_B[ii]);
        
        orb_energy_B=new double[n_ao];
        set_zero_matr(orb_energy_B,n_ao);
        for(int ii=0; ii<n_mo;ii++) xdr_double(&out_xdr,&orb_energy_B[ii]);
        
    }
    else{
        MO_VEC_B=MO_VEC;
        orb_energy_B=orb_energy;
    }
    
    xdr_destroy(&out_xdr);
    fclose(data);
    
    return 0;
}

#endif

int binary_array_read(double * M, const char * name, char * prefix, int n){

#ifdef _XDR_FULL
    return xdr_array_read(M, name, prefix, n);
#else
    
    FILE * data;
    char * full_name;
    
    get_full_name(&full_name, prefix, name);
//     printf("s = %s\n l = %d\n", full_name, strlen(full_name));
//     getchar();
    fprintf(out_stream,"reading %s\n",full_name);
    data = fopen(full_name,"r");
    if(data==NULL){
        fprintf(out_stream,"WARNING: could not open backup file %s for input\n",full_name);
        fprintf(out_stream,"         data will be calculated\n\n",full_name);
        return 1;
    };
    int flag;//password
    
    fread(&flag,sizeof(int),1,data);
    if(flag!=PASSWORD){
        fprintf(out_stream,"WARNING: binary file %s has wrong format (possibly XDR)\n",full_name);
        fclose(data);
#ifdef _XDR_READ
        xdr_array_read(M, name, prefix, n);
#else
        exit(0);        
#endif
        return 0;
    }
    
    fread(&flag,sizeof(int),1,data);
    if(flag!=n){
        fprintf(out_stream,"ERROR: inconsistent number of values in file (%d) and calculation (%d)\n",flag,n);
        fprintf(out_stream,"       it means that file %s is generated for another basis or molecule\n",name);
        exit(0);
    }
    
    set_zero_matr(M,n);
    
    for(int ii=0; ii<n;ii++) fread(M+ii,sizeof(double),1,data);

    fclose(data);
    
    return 0;  
#endif
    
}

int binary_VEC_read(const char * inp_name, double * MO_VEC, double * orb_energy, int n_mo, int n_ao, double * MO_VEC_B, double * orb_energy_B){
    
#ifdef _XDR_FULL
    xdr_VEC_read(inp_name, MO_VEC, orb_energy, n_mo, n_ao, MO_VEC_B, orb_energy_B);
#else
    FILE * data;
    data = fopen(inp_name,"rb");
    if(data==NULL){
        printf("ERROR: could not open output file %s\n\n",inp_name);
        exit(EXIT_FAILURE);
    };
    
    int flag;//password
    
    fread(&flag,sizeof(int),1,data);
    if(flag!=PASSWORD){
        fprintf(out_stream,"WARNING: binary file %s has wrong format (possibly XDR)\n",inp_name);
        fclose(data);
#ifdef _XDR_READ
        xdr_VEC_read(inp_name, MO_VEC, orb_energy, n_mo, n_ao, MO_VEC_B, orb_energy_B);
#else
        exit(0);        
#endif
        return 0;
    }
    
    fread(&flag,sizeof(int),1,data);
    if(flag>n_mo){
        fprintf(out_stream,"ERROR: inconsistent number of MO in file (%d) and calculation (%d)\n",flag,n_mo);
        fprintf(out_stream,"       it means that file %s is generated for another basis or molecule\n",inp_name);
        exit(0);
    }
    n_mo=flag;
    fread(&flag,sizeof(int),1,data);
    if(flag!=n_ao){
        fprintf(out_stream,"ERROR: inconsistent number of AO in file (%d) and calculation (%d)\n",flag,n_ao);
        fprintf(out_stream,"       it means that file %s is generated for another basis or molecule\n",inp_name);
        exit(0);
    }
    for(int ii=0; ii<n_mo*n_ao;ii++) fread(MO_VEC+ii,sizeof(double),1,data);
    for(int ii=0; ii<n_mo     ;ii++) fread(orb_energy+ii,sizeof(double),1,data);
    
    
    fread(&flag,sizeof(int),1,data);
    if(flag){        
        MO_VEC_B=new double[n_ao*n_ao];
        set_zero_matr(MO_VEC_B,n_ao*n_ao);
        for(int ii=0; ii<n_mo*n_ao;ii++) fread(MO_VEC_B+ii,sizeof(double),1,data);
        
        orb_energy_B=new double[n_ao];
        set_zero_matr(orb_energy_B,n_ao);
        for(int ii=0; ii<n_mo     ;ii++) fread(orb_energy_B+ii,sizeof(double),1,data);
        
    }
    else{
        MO_VEC_B=MO_VEC;
        orb_energy_B=orb_energy;
    }
    
    fclose(data);
#endif   
    return n_mo;
}

int binary_array_write(double * M, const char * name, char * prefix, int n, int ask){

#ifdef _XDR_FULL
    xdr_array_write(M, name, prefix, n, ask);
#else
    
    char * full_name;
    
    get_full_name(&full_name, prefix, name);
    if(ask){
        fprintf(stderr,"data will be written to %s?\npress enter to continue\n", full_name);
        getchar();
    }
    
    FILE * data;
    data = fopen(full_name,"wb");
    if(data==NULL){
		fprintf(out_stream,"ERROR: could not open output file %s\n\n", full_name);
	   	exit(EXIT_FAILURE);
    };
    
    int flag =PASSWORD;//password
    
    fwrite(&flag,sizeof(int),1,data);
    fwrite(&n,sizeof(int),1,data);
    
    for(int ii=0; ii<n;ii++) fwrite(M+ii,sizeof(double),1,data);
    
    fclose(data);
    
#endif   
    return 0;  
}


int binary_VEC_write(const char * out_name, double * MO_VEC, double * orb_energy, int n_mo, int n_ao, double * MO_VEC_B, double * orb_energy_B){

#ifdef _XDR_FULL
    xdr_VEC_write(out_name, MO_VEC, orb_energy, n_mo, n_ao, MO_VEC_B, orb_energy_B);
#else
    
    FILE * data;
    data = fopen(out_name,"wb");
    
    int flag =PASSWORD;//password
    
    fwrite(&flag,sizeof(int),1,data);
    fwrite(&n_mo,sizeof(int),1,data);
    fwrite(&n_ao,sizeof(int),1,data);
    for(int ii=0; ii<n_mo*n_ao;ii++) fwrite(MO_VEC+ii,sizeof(double),1,data);
    for(int ii=0; ii<n_mo     ;ii++) fwrite(orb_energy+ii,sizeof(double),1,data);
    
    if(MO_VEC_B!=NULL){
        flag=1;//UHF flag
        fwrite(&flag,sizeof(int),1,data);
        for(int ii=0; ii<n_mo*n_ao;ii++) fwrite(MO_VEC_B+ii,sizeof(double),1,data);
        for(int ii=0; ii<n_mo     ;ii++) fwrite(orb_energy_B+ii,sizeof(double),1,data);
    }
    else{
        flag=0;//RHF flag
        fwrite(&flag,sizeof(int),1,data);
    }
    
    fclose(data);
    
#endif
    return 0;
}
