# include <stdio.h>
# include <stdlib.h>
# include "matr.h"
# include "inp_out.h"
# include "nopa_pt.h"
# include "gv_solver.h"
# include "common_vars.h"
# include "grabbers.h"


int gv_solver::get_block_size_from_files(int ext_n_blocks, const char ** ext_names){
    
    CI_dim=0;
    n_blocks=ext_n_blocks;
    names = ext_names;
    dim = new int[n_blocks+1];
    block_dim = new int[n_blocks];
    
    int wrong_data=0;
    fprintf(out_stream,"      name | size |\n");
    fprintf(out_stream,"-----------|------|\n");
    
    for(int i=0; i<n_blocks; i++){
        block_dim[i] = read_block_size(names[i]);
        CI_dim+= block_dim[i];
        if(block_dim[i]==-1)wrong_data=1;
        fprintf(out_stream,"%10s | %4d |\n",names[i],block_dim[i]);
        
    }
    
    if(wrong_data){
        fprintf(out_stream,"There were some errors in diagonal blocks, program terminated.\n");
        fprintf(out_stream,"Non-diagonal blocks were not checked\n");
        exit(1);
    }
    
    return 0;
}

int gv_solver::get_block_size(int* ext_dim, const char** ext_names, int ext_n_blocks){
    
    CI_dim=0;
    n_blocks=ext_n_blocks;
    names = ext_names;
    dim = new int[n_blocks+1];
    block_dim = new int[n_blocks];
    
    fprintf(out_stream,"      name | size |\n");
    fprintf(out_stream,"-----------|------|\n");
    
    for(int i=0; i<n_blocks; i++){
        block_dim[i] = ext_dim[i];
        CI_dim+= block_dim[i];
        fprintf(out_stream,"%10s | %4d |\n",names[i],block_dim[i]);
    }
    
    return 0;
}



int gv_solver::alloc(){
    
    dim[0]=0;
    for(int i=0;i<n_blocks;i++)
        dim[i+1]=dim[i]+block_dim[i];
    
    
    ev = new double[CI_dim];
    H  = new double [CI_dim*CI_dim];set_zero_matr(H  , CI_dim*CI_dim);
    S  = new double [CI_dim*CI_dim];set_zero_matr(S  , CI_dim*CI_dim);
    T  = new double [CI_dim*CI_dim];set_zero_matr(T  , CI_dim*CI_dim);
    T1 = new double [CI_dim*CI_dim];set_zero_matr(T1 , CI_dim*CI_dim);
    d_x= new double [CI_dim*CI_dim];set_zero_matr(d_x, CI_dim*CI_dim);
    d_y= new double [CI_dim*CI_dim];set_zero_matr(d_y, CI_dim*CI_dim);
    d_z= new double [CI_dim*CI_dim];set_zero_matr(d_z, CI_dim*CI_dim);
    
    
    return 0;
}

int gv_solver::get_HS_data(){
    
    int wrong_data=0;
    for(int i=0; i<n_blocks; i++)
    for(int j=0; j<n_blocks; j++){
        if(read_HS_data(names[i],names[j],H  +dim[i]*CI_dim+dim[j],
                                          S  +dim[i]*CI_dim+dim[j],
                                          d_x+dim[i]*CI_dim+dim[j],
                                          d_y+dim[i]*CI_dim+dim[j],
                                          d_z+dim[i]*CI_dim+dim[j],
                        block_dim[i],block_dim[j],CI_dim))wrong_data=1;
        
    }
    
    if(wrong_data){
        fprintf(out_stream,"There were some errors in non-diagonal blocks, program terminated.\n");
        exit(1);
    }
    
    return 0;
}

int gv_solver::get_HS_from_global_memory(int i, int j){
    
    double * H_b = H  +dim[i]*CI_dim+dim[j];
    double * S_b = S  +dim[i]*CI_dim+dim[j];
    double * x_b = d_x+dim[i]*CI_dim+dim[j];
    double * y_b = d_y+dim[i]*CI_dim+dim[j];
    double * z_b = d_z+dim[i]*CI_dim+dim[j];
          
    for(int a=0; a<block_dim[i];a++)
    for(int b=0; b<block_dim[j];b++){
        H_b[a*CI_dim+b]= H_grabbed[a*grabber_n_states+b];
        S_b[a*CI_dim+b]= S_grabbed[a*grabber_n_states+b];
        x_b[a*CI_dim+b]=Dx_grabbed[a*grabber_n_states+b];
        y_b[a*CI_dim+b]=Dy_grabbed[a*grabber_n_states+b];
        z_b[a*CI_dim+b]=Dz_grabbed[a*grabber_n_states+b];
    }
    if(i==j) return 0;
    
    
    H_b = H  +dim[j]*CI_dim+dim[i];
    S_b = S  +dim[j]*CI_dim+dim[i];
    x_b = d_x+dim[j]*CI_dim+dim[i];
    y_b = d_y+dim[j]*CI_dim+dim[i];
    z_b = d_z+dim[j]*CI_dim+dim[i];
          
    for(int b=0; b<block_dim[j];b++)
    for(int a=0; a<block_dim[i];a++){
        H_b[b*CI_dim+a]= H_grabbed[a*grabber_n_states+b];
        S_b[b*CI_dim+a]= S_grabbed[a*grabber_n_states+b];
        x_b[b*CI_dim+a]=Dx_grabbed[a*grabber_n_states+b];
        y_b[b*CI_dim+a]=Dy_grabbed[a*grabber_n_states+b];
        z_b[b*CI_dim+a]=Dz_grabbed[a*grabber_n_states+b];
    }
    
    return 0;
    
}


int gv_solver::copy_HS(int i, int j, int k, int l){
    
    int warn=0;
    int d1 = block_dim[i];
    int d2 = block_dim[j];
    if(d1!=block_dim[k])warn=1;
    if(d2!=block_dim[l])warn=1;
    
    if(warn){
        printf("ERROR: copy_HS(%d,%d,%d,%d) can not copy block %dx%d to %dx%d\n", i,j,k,l,block_dim[k],block_dim[l],block_dim[i],block_dim[j]);
        exit(1);
    }
    
    double * H_bo = H  +dim[i]*CI_dim+dim[j];
    double * S_bo = S  +dim[i]*CI_dim+dim[j];
    double * x_bo = d_x+dim[i]*CI_dim+dim[j];
    double * y_bo = d_y+dim[i]*CI_dim+dim[j];
    double * z_bo = d_z+dim[i]*CI_dim+dim[j];
    double * H_bi = H  +dim[k]*CI_dim+dim[l];
    double * S_bi = S  +dim[k]*CI_dim+dim[l];
    double * x_bi = d_x+dim[k]*CI_dim+dim[l];
    double * y_bi = d_y+dim[k]*CI_dim+dim[l];
    double * z_bi = d_z+dim[k]*CI_dim+dim[l];
          
    for(int a=0; a<block_dim[i];a++)
    for(int b=0; b<block_dim[j];b++){
        H_bo[a*CI_dim+b]= H_bi[a*CI_dim+b];
        S_bo[a*CI_dim+b]= S_bi[a*CI_dim+b];
        x_bo[a*CI_dim+b]= x_bi[a*CI_dim+b];
        y_bo[a*CI_dim+b]= y_bi[a*CI_dim+b];
        z_bo[a*CI_dim+b]= z_bi[a*CI_dim+b];
    }
    if(i==j) return 0;
    
    
    H_bo = H  +dim[j]*CI_dim+dim[i];
    S_bo = S  +dim[j]*CI_dim+dim[i];
    x_bo = d_x+dim[j]*CI_dim+dim[i];
    y_bo = d_y+dim[j]*CI_dim+dim[i];
    z_bo = d_z+dim[j]*CI_dim+dim[i];
    H_bi = H  +dim[l]*CI_dim+dim[k];
    S_bi = S  +dim[l]*CI_dim+dim[k];
    x_bi = d_x+dim[l]*CI_dim+dim[k];
    y_bi = d_y+dim[l]*CI_dim+dim[k];
    z_bi = d_z+dim[l]*CI_dim+dim[k];
          
    for(int b=0; b<block_dim[j];b++)
    for(int a=0; a<block_dim[i];a++){
        H_bo[b*CI_dim+a]= H_bi[b*CI_dim+a];
        S_bo[b*CI_dim+a]= S_bi[b*CI_dim+a];
        x_bo[b*CI_dim+a]= x_bi[b*CI_dim+a];
        y_bo[b*CI_dim+a]= y_bi[b*CI_dim+a];
        z_bo[b*CI_dim+a]= z_bi[b*CI_dim+a];
    }
    
    return 0;
    
}

int gv_solver::set_zero_HS(int i, int j){
    
    
    double * H_b = H  +dim[i]*CI_dim+dim[j];
    double * S_b = S  +dim[i]*CI_dim+dim[j];
    double * x_b = d_x+dim[i]*CI_dim+dim[j];
    double * y_b = d_y+dim[i]*CI_dim+dim[j];
    double * z_b = d_z+dim[i]*CI_dim+dim[j];
          
    for(int a=0; a<block_dim[i];a++)
    for(int b=0; b<block_dim[j];b++){
        H_b[a*CI_dim+b]=0;
        S_b[a*CI_dim+b]=0;
        x_b[a*CI_dim+b]=0;
        y_b[a*CI_dim+b]=0;
        z_b[a*CI_dim+b]=0;
    }
    if(i==j) return 0;
    
    
    H_b = H  +dim[j]*CI_dim+dim[i];
    S_b = S  +dim[j]*CI_dim+dim[i];
    x_b = d_x+dim[j]*CI_dim+dim[i];
    y_b = d_y+dim[j]*CI_dim+dim[i];
    z_b = d_z+dim[j]*CI_dim+dim[i];
          
    for(int b=0; b<block_dim[j];b++)
    for(int a=0; a<block_dim[i];a++){
        H_b[b*CI_dim+a]=0;
        S_b[b*CI_dim+a]=0;
        x_b[b*CI_dim+a]=0;
        y_b[b*CI_dim+a]=0;
        z_b[b*CI_dim+a]=0;
    }
    
    return 0;
    
}


int gv_solver::print_init_data(int pH, int pS, int pD, int pE){
    
    if(pH){fprintf(out_stream,"H:\n");  fPrintMatr10(out_stream,H,CI_dim,CI_dim,0);}
    if(pS){fprintf(out_stream,"S:\n");  fPrintMatr10(out_stream,S,CI_dim,CI_dim,0);}
    if(pD){
        fprintf(out_stream,"d_x:\n");fPrintMatr10(out_stream,d_x,CI_dim,CI_dim,0);
        fprintf(out_stream,"d_y:\n");fPrintMatr10(out_stream,d_y,CI_dim,CI_dim,0);
        fprintf(out_stream,"d_z:\n");fPrintMatr10(out_stream,d_z,CI_dim,CI_dim,0);
    }
    
    if(pE){
        fprintf(out_stream,"\n\nInitial eigenvalues:\n");
        for(int i=0; i<n_blocks; i++){
            fprintf(out_stream,"\nblock %10s\n",names[i]);
            for(int j=0; j<block_dim[i]; j++)
                fprintf(out_stream,"E_old[%d] = %.10f (%.10e)\n",j,H[(dim[i]+j)*(CI_dim+1)]/S[(dim[i]+j)*(CI_dim+1)],S[(dim[i]+j)*(CI_dim+1)]);
        }
    }
    
    return 0;
}


int gv_solver::calc(){

    HC_SCE(H, S, T, T1, ev, CI_dim);
        
    fprintf(out_stream,"\n\nFinal eigenvalues:\n\n");
    PrintEnergy(ev, CI_dim,1);

    fprintf(out_stream,"\n\nEigenvectors:\n");
    for(int s=0;s<CI_dim;s++){
        fprintf(out_stream,"\nState %d:\n",s);
        for(int i=0; i<n_blocks; i++){
            fprintf(out_stream,"block %10s |",names[i]);
            for(int j=0; j<block_dim[i]; j++)
                fprintf(out_stream," % .10f",T[s*CI_dim+dim[i]+j]);
            fprintf(out_stream,"\n");
        }
    }
    
    transform(d_x,d_x,T,T1,CI_dim);
    transform(d_y,d_y,T,T1,CI_dim);
    transform(d_z,d_z,T,T1,CI_dim);
    
    fprintf(out_stream,"\n\nDipole:\n");
    PrintDipole(d_x, d_y, d_z, CI_dim);
    
    
    
    return 0;
}

