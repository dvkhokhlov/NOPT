#include "pseudo_potential.h"
#include <stdio.h>
#include <stdlib.h>

pseudo_potential::pseudo_potential(){
    n_pot=-1;
    n_oc =0;
    coeff.resize(0);
    alpha.resize(0);
    n.resize(0);
    out_core.resize(0);
    oc_2J.resize(0);
    oc_lib_shel_coeff.resize(0);
    oc_coeff.resize(0);
    oc_alpha.resize(0);
    oc_n.resize(0);
                    
}

int pseudo_potential::set_n_pot(){
    
    if(n_pot<1){
        printf("ERROR: can not create PP with %d terms\n",n_pot);
        exit(0);
    }
    
    coeff.resize(n_pot);
    alpha.resize(n_pot);
    n.resize(n_pot);
    l.resize(n_pot);
    
    return 0;
}

int pseudo_potential::set_pot_n_terms(int i, int n_terms){
    int m = coeff.size();
    if(i>=m){
        printf("ERROR: PP index is out of range %d/%d\n",i,m);
        exit(0);
    }
    
    
    coeff[i].resize(n_terms);
    alpha[i].resize(n_terms);
    n[i].resize(n_terms);
    
    return 0;
}

int pseudo_potential::print(FILE * stream, char * name, int type){
    
    if(type==0){//ECP
        fprintf(stream,"|Pseudopotential on atom %s(%3d)                                 |\n", name,atom);
        fprintf(stream,"|replacing  %3d electrons                                        |\n",n_el);
        fprintf(stream,"|________________________________________________________________|\n");
        fprintf(stream,"|              |    n    |       alpha       |         c         |\n");
        fprintf(stream,"|              |_________|___________________|___________________|\n");
    }
    else{//SO
        fprintf(stream,"|Pseudopotential on atom %s(%3d)                                                     |\n", name,atom);
        fprintf(stream,"|replacing  %3d electrons                                                            |\n",n_el);
        fprintf(stream,"|____________________________________________________________________________________|\n");
        fprintf(stream,"|              |    n    |       alpha       |         c         |     c(2L+1)/2     |\n");
        fprintf(stream,"|              |_________|___________________|___________________|___________________|\n");
    }
    
    
    for(int i=0;i<n_pot;i++)print_pot(i,stream,type);
    
    if(n_oc)if(n_pot){
        if(type==0)fprintf(stream,"|--------------|---------|-------------------|-------------------|\n");
        else       fprintf(stream,"|--------------|---------|-------------------|-------------------|-------------------|\n");
    }
    for(int i=0;i<n_oc;i++)print_oc(i,stream,type);
    
    
    return 0;
}
    

int pseudo_potential::print_pot(int i,FILE * stream, int type){
    
    if(i!=0){
        if(type==0)fprintf(stream,"|--------------|---------|-------------------|-------------------|\n");
        else       fprintf(stream,"|--------------|---------|-------------------|-------------------|-------------------|\n");
    }
    if(l[i]==-1)
        fprintf(stream,"| General term |         |                   |                   |\n");
    else{
        if(type==0)fprintf(stream,"| L =%3d       |         |                   |                   |\n",l[i]);
        else       fprintf(stream,"| L =%3d       |         |                   |                   |                   |\n",l[i]);
    }
    
    
    for(int i_g=0;i_g<n[i].size();i_g++){
        if(type==0)fprintf(stream,"|              |  %3d    |% 18.6f |% 18.6f |\n",n[i][i_g],alpha[i][i_g],coeff[i][i_g]);
        else       fprintf(stream,"|              |  %3d    |% 18.6f |% 18.6f |% 18.6f |\n",n[i][i_g],alpha[i][i_g],coeff[i][i_g],
                                  coeff[i][i_g]*(1.0*l[i]+0.5));
    }
    
    
    
    return 0;
}

int pseudo_potential::print_oc(int i, FILE * stream, int type){
    
    if(i!=0){
        if(type==0)fprintf(stream,"|--------------|---------|-------------------|-------------------|\n");
        else       fprintf(stream,"|--------------|---------|-------------------|-------------------|-------------------|\n");
    }
    
    int size=out_core[i].alpha.size();
    
    char x[30];
    if(type==0)sprintf(x,"");
    else       sprintf(x,"                   |");
    fprintf(stream,"| Out_core     |         |                   |                   |%s\n",x);
    fprintf(stream,"| orbital      |         |% 18.6f |% 18.6f |%s\n",out_core[i].alpha[0],oc_lib_shel_coeff[i][0],x);
    fprintf(stream,"| L = %3d      |         ",out_core[i].contr[0].l);
    if(size<2)
        fprintf(stream,"|                   |                   |%s\n",x);
    else
        fprintf(stream,"|% 18.6f |% 18.6f |%s\n",out_core[i].alpha[1],oc_lib_shel_coeff[i][1],x);
    fprintf(stream,"| J = %3d/2    |         ",oc_2J[i]);
    if(size<3)
        fprintf(stream,"|                   |                   |%s\n",x);
    else
        fprintf(stream,"|% 18.6f |% 18.6f |%s\n",out_core[i].alpha[2],oc_lib_shel_coeff[i][2],x);
    for(int i_s=3;i_s<out_core[i].alpha.size();i_s++)
        fprintf(stream,"|              |         |% 18.6f |% 18.6f |%s\n",out_core[i].alpha[i_s],oc_lib_shel_coeff[i][i_s],x);
    
    fprintf(stream,"| potential    |         |                   |                   |%s\n",x);
    for(int i_g=0;i_g<oc_n[i].size();i_g++)
        fprintf(stream,"|              |  %3d    |% 18.6f |% 18.6f |%s\n",oc_n[i][i_g],oc_alpha[i][i_g],oc_coeff[i][i_g],x);
    
    
    return 0;
}


#ifdef _USE_GRPP        
libgrpp_potential_t * pseudo_potential::make_grpp(int i){
    
    return libgrpp_new_potential(l[i],0,n[i].size(),n[i].data(),coeff[i].data(),alpha[i].data());
    
}

libgrpp_potential_t * pseudo_potential::make_oc_grpp(int i){
    
    return libgrpp_new_potential(out_core[i].contr[0].l,oc_2J[i],oc_n[i].size(),oc_n[i].data(),oc_coeff[i].data(),oc_alpha[i].data());
    
}


#endif
