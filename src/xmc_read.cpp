////ALL XMC reading files once will be here
# include <stdio.h>
# include <string.h>
# include <vector>
# include <stdlib.h>

# include "matr.h"

#ifndef MAX_LINE_LEN
#define MAX_LINE_LEN 1000
#endif


std::vector<double> xmc_ifitd_read(char * xmc_name){
    
    FILE * XMCFile=fopen(xmc_name,"r");
    if(XMCFile==NULL){
        printf("ERROR: xmc_ifitd_read couldn't open file with xmc data\n");
        printf("       file: %s\n",xmc_name);
        exit(1);
    }
    char tmp[MAX_LINE_LEN]/*,tmp2[MAX_LINE_LEN]*/;
    
    std::vector<double> e;
    
    fgets(tmp,MAX_LINE_LEN,XMCFile);
    while (!feof(XMCFile)){
        fgets(tmp,MAX_LINE_LEN,XMCFile);
        if(strstr(tmp," ***  REDEFINED ENERGIES OF ACTIVE ORBITALS")!=0){
            fgets(tmp,MAX_LINE_LEN,XMCFile);
            fgets(tmp,MAX_LINE_LEN,XMCFile);
            fgets(tmp,MAX_LINE_LEN,XMCFile);
            while(strstr(tmp,"ORBITAL")){
                e.push_back(atof(tmp+43));
                fgets(tmp,MAX_LINE_LEN,XMCFile);
            }
            break;
        }
    }
    
    fclose(XMCFile);
    return e;
    
}

int H0PP_read(double * Hpp_XMC, int n_s, char * xmc_name){
    
    FILE * XMCFile=fopen(xmc_name,"r");
    
    char tmp[MAX_LINE_LEN];
    fgets(tmp,MAX_LINE_LEN,XMCFile);
    
    while (strstr(tmp,"*** MATRIX H(0)pp ***")==0) 
        fgets(tmp,MAX_LINE_LEN,XMCFile);

    int Num=0;
    
    for(int i=0; i<n_s; i+=4){
        fgets(tmp,MAX_LINE_LEN,XMCFile);
        fgets(tmp,MAX_LINE_LEN,XMCFile);
        fgets(tmp,MAX_LINE_LEN,XMCFile);
        for(int j=i;j<n_s;j++)
            {fscanf(XMCFile,"%d",&Num);
            Num-=i;
            if (Num>4) Num=4;
            for (int k=0; k<Num; k++)
                {
                fscanf(XMCFile,"%lf",&Hpp_XMC[j*n_s+i+k]);
                Hpp_XMC[(i+k)*n_s+j]=Hpp_XMC[j*n_s+i+k];
                }
            }
        }

//     printf("Hpp_XMC:\n");
//     PrintMatr(Hpp_XMC,n_s,n_s,0);

    return 0;
}

int H2_read(double * H2, int n_s, char * xmc_name){
    
    printf("\nreading H(2) from %s\n\n",xmc_name);
    FILE * xmc_inp;
    xmc_inp = fopen(xmc_name,"r");
    char * tmp_str= new char[MAX_LINE_LEN];
    fgets(tmp_str,MAX_LINE_LEN,xmc_inp);
    while(strstr(tmp_str,"###   XMC-QDPT2 RESULTS  ###")==NULL)fgets(tmp_str,MAX_LINE_LEN,xmc_inp);
    fgets(tmp_str,MAX_LINE_LEN,xmc_inp);
    fgets(tmp_str,MAX_LINE_LEN,xmc_inp);
    fgets(tmp_str,MAX_LINE_LEN,xmc_inp);
//     fgets(tmp_str,MAX_LINE_LEN,xmc_inp);
//     fgets(tmp_str,MAX_LINE_LEN,xmc_inp);
//     fgets(tmp_str,MAX_LINE_LEN,xmc_inp);
//     fgets(tmp_str,MAX_LINE_LEN,xmc_inp);
    int Num=0;
    for(int i=0; i<n_s; i+=4){
        fgets(tmp_str,MAX_LINE_LEN,xmc_inp);
        fgets(tmp_str,MAX_LINE_LEN,xmc_inp);
        fgets(tmp_str,MAX_LINE_LEN,xmc_inp);
        for(int j=i;j<n_s;j++)
            {fscanf(xmc_inp,"%d",&Num);
            Num-=i;
            if (Num>4) Num=4;
            for (int k=0; k<Num; k++)
                {
                fscanf(xmc_inp,"%lf",&H2[j*n_s+i+k]);
                H2[(i+k)*n_s+j]=H2[j*n_s+i+k];
                }
            }
        }
    
    
    printf("H2:\n");
    PrintMatr(H2,n_s,n_s,0);
    
//     printf("%s\n", tmp_str);
//     set_zero_matr(S0PP,n_s*n_s);
//     for(int i=0; i<n_s;i++)S[i*n_s+1]=1;

    
    
    return 0;
}


int read_H02_with_H0PP(double * H2, double * H0PP, int n_s, char * xmc_name){
    
    double * Hpp_XMC;
    
    Hpp_XMC = new double[n_s*n_s];
    
    H0PP_read(Hpp_XMC, n_s, xmc_name);
    
    printf("\nH0PP read from %s:\n",xmc_name);
    PrintMatr10(Hpp_XMC,n_s,n_s,0);
    
    printf("\n\n");
    double * sign = new double[n_s];
    sign[0]=1.0;
    for(int i=1; i<n_s; i++){
        printf("r[%d]=%.8e\n",i,H0PP[i]/Hpp_XMC[i]);
        if(H0PP[i]/Hpp_XMC[i]>0)sign[i]=1.0;
        else sign[i]=-1.0;
    }
    
    H2_read(H2, n_s, xmc_name);

//     printf("XMC_read:\n");
//     PrintMatr(H2,n_s,n_s,0);
    
    for(int i=0; i<n_s; i++)
    for(int j=0; j<n_s; j++)
        H2[i*n_s+j]=H2[i*n_s+j]*sign[i]*sign[j];
    
//     printf("XMC_recalc:\n");
//     PrintMatr(H2,n_s,n_s,0);
    
    delete[] Hpp_XMC;
    
//     exit(0);
    
    return 0;
}
