
#include <stdio.h>
#include <stdlib.h>
#include <cstring>

#include "inp_par_read.h"
#include "keywords.h"
#include "mol_link.h"
#include "etc.h"
#include "xmc_read.h"
#include "basis_lib_read.h"
#include "defaults.h"
#include "common_vars.h"
#include "inp_out.h"


const char*ny{"no yes"};

//RWC
backup_par::backup_par(){
    
    h1=1;
//     h2=1;
//     ri=1;
    prefix=NULL;
}

int backup_par::disable(){
    
    if(h1!=1)fprintf(out_stream,"\n\n\nWARNING: storing and reading 1-electron Hamiltonian is disabled\n\n\n");
    
    h1=1;
//     h2=1;
//     ri=1;
    return 0;
    
}


int backup_par::read_group(char * inp){
    
    recursive_file P;
    char line[BUF_LINE_LENGTH];
    
    P.r_open(inp);
        
    P.r_gets(line,BUF_LINE_LENGTH);;
    while((key_word_comp(line, backup_group_start)==0)&&(!P.r_eof()))P.r_gets(line,BUF_LINE_LENGTH);;
    if(key_word_comp(line, backup_group_start)==0){
        return 0;
    }
    
    while(!P.r_eof()){
        if(key_word_comp(line, h1_backup_kw))
            h1 = kw_to_i(line, h1_backup_kw, 1);
//     if(key_word_comp(line, h2_backup_kw))
//         h2 = kw_to_i(line, h2_backup_kw, 1);

//     if(key_word_comp(line, ri_backup_kw))
//         ri = kw_to_i(line, ri_backup_kw, 1);
    
        if(key_word_comp(    line, prefix_backup_kw))
            kw_to_s(&prefix, line, prefix_backup_kw);
        
        if(key_word_comp(line, backup_group_end))break;
        P.r_gets(line,BUF_LINE_LENGTH);;
    }
    
    
    
    
    return 0;
}

backup_par::~backup_par(){
    
    if(prefix!=NULL) delete[] prefix;
}

int backup_code_print(int a){
    
    if(a==0) fprintf(out_stream,"read (%d)\n",a);
    if(a==1) fprintf(out_stream,"calc, but not write (%d)\n",a);
    if(a==2) fprintf(out_stream,"calc and write (%d)\n",a);
    return 0;
}

//RHF
rhf_par::rhf_par(){
    
    y=0;
//     huckel_guess=0;
//     h_core_guess=0;
//     read_guess  =0;
    max_it = RHF_MAX_IT_DEFAULT;
    e_conv = RHF_EN_CON_DEFAULT;
    g_conv = RHF_GR_CON_DEFAULT;
    
};

int rhf_par::read_group(char * inp){
    
    recursive_file P;
    char line[BUF_LINE_LENGTH];
    
    P.r_open(inp);
        
    P.r_gets(line,BUF_LINE_LENGTH);;
    while((key_word_comp(line, rhf_group_start)==0)&&(!P.r_eof()))P.r_gets(line,BUF_LINE_LENGTH);;
    if(key_word_comp(line, rhf_group_start)==0){
        fprintf(out_stream,"ERROR: could not find $RHF group in file %s\n",inp);
        exit(1);
    }
    
    while(!P.r_eof()){
        read_line(line);
        if(key_word_comp(line, rhf_group_end))break;
        P.r_gets(line,BUF_LINE_LENGTH);;
    }
//     if((huckel_guess+h_core_guess+read_guess)==0){
//         huckel_guess=1;
//         h_core_guess=0;
//         read_guess  =0;
//     }
        
//     if((huckel_guess+h_core_guess+read_guess)!=1){
//         fprintf(out_stream,"ERROR: multiple choice of starting guess\n");
//         fprintf(out_stream,"       HUCKEL - %c%c%c\n",ny[3*huckel_guess],ny[3*huckel_guess+1],ny[3*huckel_guess+2]);
//         fprintf(out_stream,"       H_CORE - %c%c%c\n",ny[3*h_core_guess],ny[3*h_core_guess+1],ny[3*h_core_guess+2]);
//         fprintf(out_stream,"       READ   - %c%c%c\n",ny[3*  read_guess],ny[3*  read_guess+1],ny[3*  read_guess+2]);
//         fprintf(out_stream,"       you must choose only one.\n");
//         exit(1);
//     }
    
    
    
    
    return 0;
}

int rhf_par::read_line(char * inp){
    
    if(key_word_comp(inp, rhf_huckel)){
        fprintf(out_stream,"ERROR: while parsing $RHF group keyword HUCKEL is deprecated\n");
        exit(0);
    }

    if(key_word_comp(inp, rhf_h_core)){
        fprintf(out_stream,"ERROR: while parsing $RHF group keyword H_CORE is deprecated\n");
        exit(0);
    }
    
//     if(key_word_comp(inp, rhf_read))
//         read_guess=1;
    
    if(key_word_comp(inp, max_it_kw)){
        max_it = kw_to_i(inp, max_it_kw, RHF_MAX_IT_DEFAULT);
    }
    
    if(key_word_comp(inp, e_conv_kw)){
        e_conv = kw_to_f(inp, e_conv_kw, RHF_EN_CON_DEFAULT);
    }
    if(key_word_comp(inp, g_conv_kw)){
        g_conv = kw_to_f(inp, g_conv_kw, RHF_GR_CON_DEFAULT);
    }
    
    return 0;
}

int rhf_par::write_info(int n_cor){
    
    int s;
    fprintf(out_stream,"\n\n\n");
    fprintf(out_stream,"_______________________Starting_RHF_calculations_______________________\n\n");
    fprintf(out_stream,"Energy convergence:               %e\n",e_conv);
    fprintf(out_stream,"Orbital gradient convergence:     %e\n",g_conv);
    fprintf(out_stream,"Maximum number of iterations:     %d\n",max_it);
    fprintf(out_stream,"\n");
    fprintf(out_stream,"Number of core orbitals:          %d\n",n_cor);
    
//     fprintf(out_stream,"Starting guess:                   ");
//     if(huckel_guess)fprintf(out_stream,"HUCKEL\n");
//     if(h_core_guess)fprintf(out_stream,"H CORE\n");
//     if(  read_guess)fprintf(out_stream,"Read from file\n");
    
    fprintf(out_stream,"\n");
//     fprintf(out_stream,"Folder for output files: %s\n",out_folder_name);
    fprintf(out_stream,"_______________________________________________________________________\n\n\n");
    
    
    return 0;
    
}

rhf_par::~rhf_par(){
    
}

//CAS
cas_par::cas_par(){
    y=0;
    ci_solver = CISOLVER_ALDET;
    //convergence
    max_it = CAS_MAX_IT_DEFAULT;
    e_conv = CAS_EN_CON_DEFAULT;
    g_conv = CAS_GR_CON_DEFAULT;
    s_conv = CAS_STEP_CON_DEFAULT;
    x_max  = CAS_X_MAX_DEFAULT;
    method=0;
    
    //states
    n_s=1;
    w_state.resize(1);
    w_state[0]=1.0;
    
    track=0;
    rotate_orbs=1;
    w_state_type=-1;
    // Lambda=-1;

        
}

int cas_par::w_linear_read_line(char * inp){
    int found_L=0;
    int found_S=0;
    int found_W=0;
    int found_M=0;
    int L,S;
    double M;
    std::vector<double>W;
    
    
    if(key_word_comp(inp, w_linear_L_kw)){
        L = kw_to_i (inp, w_linear_L_kw, -18);
        if(L!=-18)found_L=1;
    }
    
    if(key_word_comp(inp, w_linear_sign_kw)){
        S = kw_to_i (inp, w_linear_sign_kw, 0);
        found_S=1;
    }
    
    if(key_word_comp(inp, mult_kw)){
        M = kw_to_f (inp, mult_kw, 0);
        found_M=1;
    }
    
    if(key_word_comp(inp, w_state_kw)){
        int n=kw_count(inp, w_state_kw,';');
        W.resize(n);
        kw_to_d_v(&W, inp, w_state_kw, n);
        cut_zero(&W);
        if(n) found_W=1;
    }
    
    if((found_L+found_M+found_S+found_W)==0)return 0;
    
    if(found_L*found_M*found_W*(found_S+L)==0){
        fprintf(out_stream,"ERROR: could not find w_state, multiplicity");
        if(L==0)fprintf(out_stream,", sign");
        fprintf(out_stream," and/or lambda while pharsing $W_STATE group. Line:\n");
        fprintf(out_stream,"       - %s\n",inp);
        exit(0);
    }
    if(L==0)if((S!=1)&&(S!=-1)){
        fprintf(out_stream,"ERROR: wrong value -- sign=%d\n",S);
        fprintf(out_stream,"       sign can be only +1 or -1 for Lambda=0\n");
        fprintf(out_stream,"       error is in line:\n%s",inp);
        exit(1);
    }
    
    
    w_state_by_rep.push_back(W);
    rep_lambda.push_back(L);
    rep_sign.push_back(S);
    rep_spin.push_back((M-1.0)/2.0);
    
    
    
    return 0;
}

int cas_par::w_state_by_rep_read(char * inp){
    
    recursive_file P;
    char line[BUF_LINE_LENGTH];
    if(LINEAR==0){
        nopt_printf("ERROR: IRREP_W_STATE is written only for group=LINEAR inp_par_read.cpp:225");exit(1);
    }
    
    P.r_open(inp);
        
    P.r_gets(line,BUF_LINE_LENGTH);;
    while((key_word_comp(line, w_state_group_start)==0)&&(!P.r_eof()))P.r_gets(line,BUF_LINE_LENGTH);;
    if(key_word_comp(line, w_state_group_start)==0){
        fprintf(out_stream,"ERROR: could not find $W_STATE group in file %s\n",inp);
        fprintf(out_stream,"       $W_STATE group is needed because w_state=ir.rep is set in $CAS group\n");
        exit(1);
    }
    P.r_gets(line,BUF_LINE_LENGTH);
    while(!P.r_eof()){
        w_linear_read_line(line);
        P.r_gets(line,BUF_LINE_LENGTH);;
        if(key_word_comp(line, w_state_group_end))break;
    }
    return 0;


}
int cas_par::read_group(char * inp){

    dav.read_group(inp);
    dmrg.read_group(inp);

    recursive_file P;
    char line[BUF_LINE_LENGTH];
    
    P.r_open(inp);
        
    P.r_gets(line,BUF_LINE_LENGTH);;
    while((key_word_comp(line, cas_group_start)==0)&&(!P.r_eof()))P.r_gets(line,BUF_LINE_LENGTH);;
    if(key_word_comp(line, cas_group_start)==0){
        fprintf(out_stream,"ERROR: could not find $CAS group in file %s\n",inp);
        exit(1);
    }
    
    while(!P.r_eof()){
        read_line(line);
        if(key_word_comp(line, cas_group_end))break;
        P.r_gets(line,BUF_LINE_LENGTH);;
    }
    
    if(dav.n_s==-1)
        dav.n_s = max(DAV_N_STATES_MULT_DEFAULT*n_s, DAV_N_STATES_DEFAULT);
    
    if(w_state_type==1){
        w_state.resize(n_s);
        for(int i=0;i<n_s;i++)w_state[i]=1.0;
    }
    if(w_state_type==2)w_state_by_rep_read(inp);
    
    if(n_s<w_state.size()){
        n_s=w_state.size();
        fprintf(out_stream,"\n\nWARNING: n_s was changed to %d due to large w_state length\n\n",n_s);
    }
    
    if(method>3 ){
        nopt_printf("ERROR: duplicate definition of CASSCF method\n");
        nopt_printf("       use SA or SM, but not both\n");
        exit(1);
    }
    method--;
    if(method==-1)method=CAS_METHOD_DEFAULT;

    if(ci_solver==CISOLVER_DMRG) dmrg.validate();

    return 0;
}

int cas_par::read_line(char * inp){
    
    if(key_word_comp(inp, max_it_kw)){
        max_it = kw_to_i(inp, max_it_kw, CAS_MAX_IT_DEFAULT);
    }
    
    
    
    if(key_word_comp(inp, cisolver_kw)){
        if      (kw_to_kw(inp, cisolver_kw, cisolver_aldet_kw)) ci_solver = CISOLVER_ALDET;
        else if (kw_to_kw(inp, cisolver_kw, cisolver_dmrg_kw )) ci_solver = CISOLVER_DMRG;
        else{
            fprintf(out_stream,"ERROR: unknown CISOLVER value; accepted: aldet, dmrg\n");
            exit(1);
        }
    }

    if(key_word_comp(inp, cas_track_kw)){
        track=1;
    }
    
    if(key_word_comp(inp, cas_rotate_orbs_kw)){
        rotate_orbs = kw_to_i(inp, cas_rotate_orbs_kw, 1);
    }
    
    if(key_word_comp(inp, e_conv_kw)){
        e_conv = kw_to_f(inp, e_conv_kw, CAS_EN_CON_DEFAULT);
    }
    
    if(key_word_comp(inp, g_conv_kw)){
        g_conv = kw_to_f(inp, g_conv_kw, CAS_GR_CON_DEFAULT);
    }
    
    if(key_word_comp(inp, s_conv_kw)){
        s_conv = kw_to_f(inp, s_conv_kw, CAS_STEP_CON_DEFAULT);
    }
    
    if(key_word_comp(inp, x_max_kw)){
        x_max = kw_to_f(inp, x_max_kw, CAS_X_MAX_DEFAULT);
    }
    
    if(key_word_comp(inp, cas_SA_kw)){
        method+=2;
    }
    if(key_word_comp(inp, cas_SM_kw)){
        method+=3;
    }
        
    
    if(key_word_comp(inp, num_state_kw)){
        n_s = kw_to_i(inp, num_state_kw, 1);
    }
    
    if(key_word_comp(inp, w_state_kw)){
        if(w_state_type!=-1){
            fprintf(out_stream,"ERROR: multiple definition of w_state\n");
            exit(1);
        }
        if(kw_to_kw(inp,w_state_kw,w_state_eq_kw)   ){w_state_type=1;goto ws_done;}
        if(kw_to_kw(inp,w_state_kw,w_state_irrep_kw)){w_state_type=2;goto ws_done;}
        w_state_type=0;
        int n=kw_count(inp, w_state_kw,';');
        w_state.resize(n);
        kw_to_d_v(&w_state, inp, w_state_kw, n);
        cut_zero(&w_state);
    }
    ws_done:
        
    return 0;
}

int cas_par::write_info(int n_a, int n_b, int n_o, int n_c, int mult){
    
    int s;
    fprintf(out_stream,"\n\n\n");
    fprintf(out_stream,"______________________Starting_CASSCF_calculations_____________________\n\n");
    
    fprintf(out_stream,"Number of states                  %d\n",n_s );
    fprintf(out_stream,"State multiplicity                %d " ,mult);
    if(mult==1) fprintf(out_stream,"(singlet)\n");
    if(mult==2) fprintf(out_stream,"(doublet)\n");
    if(mult==3) fprintf(out_stream,"(triplet)\n");
    if(mult==4) fprintf(out_stream,"(quadruplet)\n");
    if(mult >4) fprintf(out_stream,"(high-spin)\n");
    if(mult==0) fprintf(out_stream,"(all)\n");
    // if(LINEAR){
    // fprintf(out_stream,"State angle moment                ");
    // if(Lambda==-1)fprintf(out_stream,"(all)\n");
    // else          fprintf(out_stream,"%d\n",Lambda);
    // }
    fprintf(out_stream,"\n");
    fprintf(out_stream,"Active space:\n");
    fprintf(out_stream,"Number of alpha electrons         %d\n",n_a);
    fprintf(out_stream,"Number of beta  electrons         %d\n",n_b);
    fprintf(out_stream,"Number of active orbitals         %d\n",n_o);
    fprintf(out_stream,"\n");
    fprintf(out_stream,"Number of core orbitals:          %d\n",n_c);
    fprintf(out_stream,"\n");
    fprintf(out_stream,"Preparation of virtual orbitals:  %c%c%c\n",ny[3*rotate_orbs],ny[3*rotate_orbs+1],ny[3*rotate_orbs+2]);
    fprintf(out_stream,"State tracking:                   %c%c%c\n",ny[3*track      ],ny[3*track      +1],ny[3*track      +2]);
    if(method==1)fprintf(out_stream,"State averaging:\n");
    if(method==2)fprintf(out_stream,"Separate minimization:\n");
    if(w_state_type==0){
        fprintf(out_stream,"Weights:\t");
        for(auto &w: w_state) fprintf(out_stream,"% .2e ",w);
        fprintf(out_stream,"\n");
    }
    if(w_state_type==1)fprintf(out_stream,"Weights:\tall equal to 1\n");
    if(w_state_type==2){
        for(int i=0;i<w_state_by_rep.size();i++){
            fprintf(out_stream,"%4d) ",i);
            fprintf(out_stream,"L =%3d ",rep_lambda[i]);
            if(rep_lambda[i]==0)if(rep_sign[i]== 1)fprintf(out_stream,"(+), ");
            if(rep_lambda[i]==0)if(rep_sign[i]==-1)fprintf(out_stream,"(-), ");
            if(rep_lambda[i]!=0)                   fprintf(out_stream,"   , ");
            fprintf(out_stream,"S=%4.1f, ",rep_spin[i]);
            fprintf(out_stream,"Weights:\t");
            for(auto &w: w_state_by_rep[i]) fprintf(out_stream,"% .2e ",w);
            fprintf(out_stream,";\n");
        }
    }
    
    fprintf(out_stream,"\n");
    fprintf(out_stream,"SCF settings:\n");
    fprintf(out_stream,"Energy convergence:               %e\n",e_conv);
    fprintf(out_stream,"Orbital gradient convergence:     %e\n",g_conv);
    fprintf(out_stream,"Rotation matrix convergence :     %e\n",s_conv);
    fprintf(out_stream,"Maximum number of iterations:     %d\n",max_it);
    fprintf(out_stream,"Maximum SOSCF step          :     %e\n",x_max);
    fprintf(out_stream,"\n");
    if      (ci_solver==CISOLVER_ALDET){
        fprintf(out_stream,"CI solver:                        determinant CI (aldet)\n");
        dav.write_info();
    }
    else if (ci_solver==CISOLVER_DMRG){
        fprintf(out_stream,"CI solver:                        DMRG\n");
        dmrg.write_info();
    }
    else{
        fprintf(out_stream,"CI solver:                        unknown (%d)\n",ci_solver);
    }
    
//     fprintf(out_stream,"Folder for output files: %s\n",out_folder_name);
    fprintf(out_stream,"_______________________________________________________________________\n\n\n");
    
    
    return 0;
    
}

cas_par::~cas_par(){
    
}

//-----------------------------------------------------------------------------------------------------

//CIS
cis_par::cis_par(){
    y=0;

    n_s = 1;
    n_f = 0;
    method = 1;
    nat_orb = 0;		// новое

}

int cis_par::read_group(char * inp){
    
    dav.read_group(inp);
    
    recursive_file P;
    char line[BUF_LINE_LENGTH];
    
    P.r_open(inp);
        
    P.r_gets(line,BUF_LINE_LENGTH);;
    while((key_word_comp(line, cis_group_start)==0)&&(!P.r_eof()))P.r_gets(line,BUF_LINE_LENGTH);;
    if(key_word_comp(line, cis_group_start)==0){
        fprintf(out_stream,"ERROR: could not find $CIS group in file %s\n",inp);
        exit(1);
    }
    
    while(!P.r_eof()){
        read_line(line);
        if(key_word_comp(line, cis_group_end))break;
        P.r_gets(line,BUF_LINE_LENGTH);;
    }

    if ( (method!=1)&&(method!=2) )
    {
    fprintf(out_stream, "WARNING: Method = %d is not available. \n", method);
    fprintf(out_stream, "         Method = 1 is used \n");
    method = 1;
    }

    if(dav.n_s==-1)
        dav.n_s = max(DAV_N_STATES_MULT_DEFAULT*n_s, DAV_N_STATES_DEFAULT);
    
    
    
    if(n_s<w_state.size())
    {
        n_s=w_state.size();
        fprintf(out_stream,"\n\n WARNING: n_s was changed to %d due to large w_state length\n\n", n_s);
    }
        
    return 0;
}

int cis_par::read_line(char * inp){
    
    if(key_word_comp(inp, num_frozen_orb_kw)){
        n_f = kw_to_i(inp, num_frozen_orb_kw, 0);
    }
    
    if(key_word_comp(inp, num_state_kw)){
        n_s = kw_to_i(inp, num_state_kw, 1);
    }
    
    if(key_word_comp(inp, method_kw)){
        method = kw_to_i(inp, method_kw, 1);
    }

    if(key_word_comp(inp, nat_orb_kw)){
        nat_orb = kw_to_i(inp, nat_orb_kw, 0);
    }

    if(key_word_comp(inp, w_state_kw)){
        int n=kw_count(inp, w_state_kw,';');
        w_state.resize(n);
        kw_to_d_v(&w_state, inp, w_state_kw, n);
        cut_zero(&w_state);
    }
            
    return 0;
}

int cis_par::write_info(){
    
    fprintf(out_stream,"\n\n\n");
    fprintf(out_stream,"______________________Starting_CIS_calculations________________________\n\n");
    
    fprintf(out_stream,"Number of states                  %d\n",n_s );
    fprintf(out_stream,"Number of frozen orbitals         %d\n",n_f );
    fprintf(out_stream,"Method used                       CIS_%d\n",method);
    if (nat_orb){fprintf(out_stream,"Make Natural Orbitals             yes");}
    else{fprintf(out_stream,"Make Natural Orbitals             no");}
    
    fprintf(out_stream,"\n");
    dav.write_info();
    
//     fprintf(out_stream,"Folder for output files: %s\n",out_folder_name);
    fprintf(out_stream,"_______________________________________________________________________\n\n\n");
    
    
    return 0;
    
}

cis_par::~cis_par(){
    
}


//-----------------------------------------------------------------------------------------------------

dav_par::dav_par(){
    
    n_s    = -1;
    n_bf   = DAV_N_BF_DEFAULT;
    max_it = DAV_MAX_IT_DEFAULT;
    e_conv = DAV_EN_CON_DEFAULT;
    r_conv = DAV_R_CON_DEFAULT;
    se_min = DAV_SE_MIN_DEFAULT;
    edshift= DAV_ED_SHIFT_DEFAULT;
    sparsed_Hc = 0;
    
}

int dav_par::read_group(char * inp){
    
    recursive_file P;
    char line[BUF_LINE_LENGTH];
    
    P.r_open(inp);
        
    P.r_gets(line,BUF_LINE_LENGTH);;
    while((key_word_comp(line, dav_group_start)==0)&&(!P.r_eof()))P.r_gets(line,BUF_LINE_LENGTH);;
    if(key_word_comp(line, dav_group_start)==0){
        return 0;
    }
    
    while(!P.r_eof()){
        read_line(line);
        if(key_word_comp(line, dav_group_end))break;
        P.r_gets(line,BUF_LINE_LENGTH);;
    }
    
    
    
    return 0;
}

int dav_par::read_line(char * inp){
    
    if(key_word_comp(inp, num_state_kw)){
        n_s = kw_to_f(inp, num_state_kw, -1);
    }
    
    if(key_word_comp(inp, dav_n_bf_kw)){
        n_bf = kw_to_i(inp, dav_n_bf_kw, DAV_N_BF_DEFAULT);
    }
    
    if(key_word_comp(inp, max_it_kw)){
        max_it = kw_to_i(inp, max_it_kw, DAV_MAX_IT_DEFAULT);
    }
    
    if(key_word_comp(inp, e_conv_kw)){
        e_conv = kw_to_f(inp, e_conv_kw, DAV_EN_CON_DEFAULT);
    }
    
    if(key_word_comp(inp, r_conv_kw)){
        r_conv = kw_to_f(inp, r_conv_kw, DAV_R_CON_DEFAULT);
    }
    
    if(key_word_comp(inp, dav_se_kw)){
        se_min = kw_to_f(inp, dav_se_kw, DAV_SE_MIN_DEFAULT);
    }
    
    if(key_word_comp(inp, edshift_kw)){
        edshift = kw_to_f(inp, edshift_kw, DAV_ED_SHIFT_DEFAULT);
    }

    if(key_word_comp(inp, dav_sparsed_HC_kw)){
        sparsed_Hc = kw_to_i(inp, dav_sparsed_HC_kw, 0);
    }

    
    
    return 0;
}

int dav_par::write_info(){
    
    fprintf(out_stream,"Davidson diagonalization settings:\n");
    fprintf(out_stream,"Number of correlation vectors:    %d\n",n_s    );
    fprintf(out_stream,"Number of basis functions:        %d\n",n_bf   );
    fprintf(out_stream,"Maximum number of iterations:     %d\n",max_it );
    fprintf(out_stream,"Energy convergence:               %e\n",e_conv );
    fprintf(out_stream,"Residue convergence:              %e\n",r_conv );
    fprintf(out_stream,"Orthogonality cut-off:            %e\n",se_min );
    fprintf(out_stream,"Energy denominator shift:         %e\n",edshift);
    fprintf(out_stream,"Sparsed multiplication:           %s\n",sparsed_Hc?"yes":"no");
    fprintf(out_stream,"\n");
    
    return 0;
}

dav_par::~dav_par(){

}

//-----------------------------------------------------------------------------------------------------

//DMRG
dmrg_par::dmrg_par(){

    m         = DMRG_M_DEFAULT;          // 0 = unset; required when CISOLVER=dmrg
    sweeps    = DMRG_SWEEPS_DEFAULT;
    sweep_tol = DMRG_SWEEP_TOL_DEFAULT;
    hf_occ    = DMRG_HF_OCC_INTEGRAL;
    schedule  = DMRG_SCHED_DEFAULT;
    save_dir  = DMRG_SAVE_DIR_DEFAULT;
    memory    = DMRG_MEMORY_DEFAULT;
    localize  = DMRG_LOC_OFF;
    dump_loc_orbs = 0;
    loc_order = DMRG_LOCORDER_FIEDLER;
    warm_start       = DMRG_WARM_START_DEFAULT;
    warm_sweeps      = DMRG_WARM_SWEEPS_DEFAULT;
    rot_m            = DMRG_ROT_M_DEFAULT;
    rot_steps        = DMRG_ROT_STEPS_DEFAULT;
    warm_start_after = DMRG_WARM_START_AFTER_DEFAULT;
    warm_rotate      = DMRG_WARM_ROTATE_DEFAULT;
    print_dets       = DMRG_PRINT_DETS_DEFAULT;
    det_rot_m        = DMRG_DET_ROT_M_DEFAULT;
    det_rot_steps    = DMRG_DET_ROT_STEPS_DEFAULT;
    extract_m        = DMRG_EXTRACT_M_DEFAULT;
    extract_cutoff   = DMRG_EXTRACT_CUTOFF_DEFAULT;

}

int dmrg_par::read_group(char * inp){

    recursive_file P;
    char line[BUF_LINE_LENGTH];

    P.r_open(inp);

    P.r_gets(line,BUF_LINE_LENGTH);;
    while((key_word_comp(line, dmrg_group_start)==0)&&(!P.r_eof()))P.r_gets(line,BUF_LINE_LENGTH);;
    if(key_word_comp(line, dmrg_group_start)==0){
        return 0;
    }

    while(!P.r_eof()){
        read_line(line);
        if(key_word_comp(line, dmrg_group_end))break;
        P.r_gets(line,BUF_LINE_LENGTH);;
    }

    return 0;
}

int dmrg_par::read_line(char * inp){

    if(key_word_comp(inp, dmrg_m_kw)){
        m = kw_to_i(inp, dmrg_m_kw, DMRG_M_DEFAULT);
    }

    if(key_word_comp(inp, dmrg_sweeps_kw)){
        sweeps = kw_to_i(inp, dmrg_sweeps_kw, DMRG_SWEEPS_DEFAULT);
    }

    if(key_word_comp(inp, dmrg_sweep_tol_kw)){
        sweep_tol = kw_to_f(inp, dmrg_sweep_tol_kw, DMRG_SWEEP_TOL_DEFAULT);
    }

    if(key_word_comp(inp, dmrg_hf_occ_kw)){
        if(kw_to_kw(inp, dmrg_hf_occ_kw, dmrg_hf_occ_integral_kw)) hf_occ = DMRG_HF_OCC_INTEGRAL;
        else                                                       hf_occ = DMRG_HF_OCC_UNKNOWN;
    }

    if(key_word_comp(inp, dmrg_schedule_kw)){
        if(kw_to_kw(inp, dmrg_schedule_kw, dmrg_schedule_default_kw)) schedule = DMRG_SCHED_DEFAULT;
        else                                                          schedule = DMRG_SCHED_UNKNOWN;
    }

    if(key_word_comp(inp, dmrg_save_dir_kw)){
        char* tmp=nullptr;
        kw_to_s(&tmp, inp, dmrg_save_dir_kw);
        if(tmp){ save_dir=tmp; delete[] tmp; }
    }

    if(key_word_comp(inp, dmrg_memory_kw)){
        memory = kw_to_f(inp, dmrg_memory_kw, DMRG_MEMORY_DEFAULT);
    }

    if(key_word_comp(inp, dmrg_localize_kw)){
        if     (kw_to_kw(inp, dmrg_localize_kw, dmrg_localize_off_kw))  localize = DMRG_LOC_OFF;
        else if(kw_to_kw(inp, dmrg_localize_kw, dmrg_localize_pm_kw))   localize = DMRG_LOC_PM;
        else if(kw_to_kw(inp, dmrg_localize_kw, dmrg_localize_boys_kw)) localize = DMRG_LOC_BOYS;
        else                                                           localize = DMRG_LOC_UNKNOWN;
    }

    if(key_word_comp(inp, dmrg_dump_loc_kw)){
        dump_loc_orbs = 1;
    }

    if(key_word_comp(inp, dmrg_loc_order_kw)){
        if     (kw_to_kw(inp, dmrg_loc_order_kw, dmrg_loc_order_fiedler_kw)) loc_order = DMRG_LOCORDER_FIEDLER;
        else if(kw_to_kw(inp, dmrg_loc_order_kw, dmrg_loc_order_gaopt_kw))   loc_order = DMRG_LOCORDER_GAOPT;
        else if(kw_to_kw(inp, dmrg_loc_order_kw, dmrg_loc_order_none_kw))    loc_order = DMRG_LOCORDER_NONE;
        else                                                                loc_order = DMRG_LOCORDER_UNKNOWN;
    }

    if(key_word_comp(inp, dmrg_warm_start_kw)){
        if     (kw_to_kw(inp, dmrg_warm_start_kw, dmrg_warm_off_kw)) warm_start = DMRG_WARM_OFF;
        else if(kw_to_kw(inp, dmrg_warm_start_kw, dmrg_warm_on_kw))  warm_start = DMRG_WARM_ON;
        else                                                         warm_start = DMRG_WARM_UNKNOWN;
    }

    if(key_word_comp(inp, dmrg_warm_sweeps_kw)){
        warm_sweeps = kw_to_i(inp, dmrg_warm_sweeps_kw, DMRG_WARM_SWEEPS_DEFAULT);
    }

    if(key_word_comp(inp, dmrg_rot_m_kw)){
        rot_m = kw_to_i(inp, dmrg_rot_m_kw, DMRG_ROT_M_DEFAULT);
    }

    if(key_word_comp(inp, dmrg_rot_steps_kw)){
        rot_steps = kw_to_i(inp, dmrg_rot_steps_kw, DMRG_ROT_STEPS_DEFAULT);
    }

    if(key_word_comp(inp, dmrg_warm_start_after_kw)){
        warm_start_after = kw_to_i(inp, dmrg_warm_start_after_kw, DMRG_WARM_START_AFTER_DEFAULT);
    }

    if(key_word_comp(inp, dmrg_warm_rotate_kw)){
        if     (kw_to_kw(inp, dmrg_warm_rotate_kw, dmrg_warm_off_kw)) warm_rotate = DMRG_WARM_OFF;
        else if(kw_to_kw(inp, dmrg_warm_rotate_kw, dmrg_warm_on_kw))  warm_rotate = DMRG_WARM_ON;
        else                                                          warm_rotate = DMRG_WARM_UNKNOWN;
    }

    if(key_word_comp(inp, dmrg_print_dets_kw)){
        if     (kw_to_kw(inp, dmrg_print_dets_kw, dmrg_warm_off_kw)) print_dets = DMRG_WARM_OFF;
        else if(kw_to_kw(inp, dmrg_print_dets_kw, dmrg_warm_on_kw))  print_dets = DMRG_WARM_ON;
        else                                                         print_dets = DMRG_WARM_UNKNOWN;
    }

    if(key_word_comp(inp, dmrg_det_rot_m_kw)){
        det_rot_m = kw_to_i(inp, dmrg_det_rot_m_kw, DMRG_DET_ROT_M_DEFAULT);
    }

    if(key_word_comp(inp, dmrg_det_rot_steps_kw)){
        det_rot_steps = kw_to_i(inp, dmrg_det_rot_steps_kw, DMRG_DET_ROT_STEPS_DEFAULT);
    }

    if(key_word_comp(inp, dmrg_extract_m_kw)){
        extract_m = kw_to_i(inp, dmrg_extract_m_kw, DMRG_EXTRACT_M_DEFAULT);
    }

    if(key_word_comp(inp, dmrg_extract_cutoff_kw)){
        extract_cutoff = kw_to_f(inp, dmrg_extract_cutoff_kw, DMRG_EXTRACT_CUTOFF_DEFAULT);
    }

    return 0;
}

int dmrg_par::validate(){

    int ok=1;

    if(m<=0){
        fprintf(out_stream,"ERROR: $DMRG bond dimension m=%d must be > 0 (set 'm' in $DMRG)\n",m);
        ok=0;
    }
    if(sweeps<=0){
        fprintf(out_stream,"ERROR: $DMRG sweeps=%d must be > 0\n",sweeps);
        ok=0;
    }
    if(warm_sweeps<=0){ warm_sweeps = sweeps/2; if(warm_sweeps<1) warm_sweeps = 1; } // auto = sweeps/2
    if(sweep_tol<=0){
        fprintf(out_stream,"ERROR: $DMRG sweep_tol=%e must be > 0\n",sweep_tol);
        ok=0;
    }
    if(hf_occ==DMRG_HF_OCC_UNKNOWN){
        fprintf(out_stream,"ERROR: $DMRG unknown hf_occ value; accepted: integral\n");
        ok=0;
    }
    if(schedule==DMRG_SCHED_UNKNOWN){
        fprintf(out_stream,"ERROR: $DMRG unknown schedule value; accepted: default\n");
        ok=0;
    }
    if(localize==DMRG_LOC_UNKNOWN){
        fprintf(out_stream,"ERROR: $DMRG unknown localize value; accepted: off, pm, boys\n");
        ok=0;
    }
    if(localize==DMRG_LOC_BOYS){
        fprintf(out_stream,"ERROR: $DMRG localize=boys not implemented yet; accepted: off, pm\n");
        ok=0;
    }
    if(loc_order==DMRG_LOCORDER_UNKNOWN){
        fprintf(out_stream,"ERROR: $DMRG unknown loc_order value; accepted: fiedler, none\n");
        ok=0;
    }
    if(loc_order==DMRG_LOCORDER_GAOPT){
        fprintf(out_stream,"ERROR: $DMRG loc_order=gaopt not implemented yet; accepted: fiedler, none\n");
        ok=0;
    }
    if(save_dir.empty()){
        fprintf(out_stream,"ERROR: $DMRG save_dir must not be empty\n");
        ok=0;
    }
    if(memory<=0){
        fprintf(out_stream,"ERROR: $DMRG memory=%g must be > 0 (block2 double-stack size in GB)\n",memory);
        ok=0;
    }
    if(warm_start==DMRG_WARM_UNKNOWN){
        fprintf(out_stream,"ERROR: $DMRG unknown warm_start value; accepted: off, on\n");
        ok=0;
    }
    if(warm_rotate==DMRG_WARM_UNKNOWN){
        fprintf(out_stream,"ERROR: $DMRG unknown warm_rotate value; accepted: off, on\n");
        ok=0;
    }
    if(print_dets==DMRG_WARM_UNKNOWN){
        fprintf(out_stream,"ERROR: $DMRG unknown print_dets value; accepted: off, on\n");
        ok=0;
    }
    if(det_rot_m<=0){ det_rot_m = 2*m < 1500 ? 2*m : 1500; } // auto = min(2m, 1500)
    if(det_rot_steps<=0){
        fprintf(out_stream,"ERROR: $DMRG det_rot_steps=%d must be > 0\n",det_rot_steps);
        ok=0;
    }
    if(extract_m<0){
        fprintf(out_stream,"ERROR: $DMRG extract_m=%d must be >= 0 (0 = no compression)\n",extract_m);
        ok=0;
    }
    if(extract_cutoff<=0){
        fprintf(out_stream,"ERROR: $DMRG extract_cutoff=%e must be > 0\n",extract_cutoff);
        ok=0;
    }
    if(warm_start==DMRG_WARM_ON){
        if(rot_m<0){
            fprintf(out_stream,"ERROR: $DMRG rot_m=%d must be >= 0 (0 = use m)\n",rot_m);
            ok=0;
        }
        if(rot_steps<=0){
            fprintf(out_stream,"ERROR: $DMRG rot_steps=%d must be > 0\n",rot_steps);
            ok=0;
        }
        if(warm_start_after<0){
            fprintf(out_stream,"ERROR: $DMRG warm_start_after=%d must be >= 0\n",warm_start_after);
            ok=0;
        }
    }

    if(!ok)exit(1);

    return 0;
}

int dmrg_par::write_info(){

    fprintf(out_stream,"DMRG settings:\n");
    fprintf(out_stream,"Bond dimension (m):               %d\n",m);
    fprintf(out_stream,"Maximum number of sweeps:         %d\n",sweeps);
    fprintf(out_stream,"Sweep energy convergence:         %e\n",sweep_tol);
    if(hf_occ==DMRG_HF_OCC_INTEGRAL)
        fprintf(out_stream,"Initial occupancy:                integral\n");
    if(schedule==DMRG_SCHED_DEFAULT)
        fprintf(out_stream,"Sweep schedule:                   default\n");
    if(localize==DMRG_LOC_OFF)
        fprintf(out_stream,"Active-space localization:        off\n");
    if(localize==DMRG_LOC_PM)
        fprintf(out_stream,"Active-space localization:        Pipek-Mezey\n");
    if(dump_loc_orbs)
        fprintf(out_stream,"Dump localized orbitals:          yes\n");
    if(loc_order==DMRG_LOCORDER_FIEDLER)
        fprintf(out_stream,"DMRG orbital ordering:            Fiedler\n");
    if(loc_order==DMRG_LOCORDER_NONE)
        fprintf(out_stream,"DMRG orbital ordering:            none (input order)\n");
    fprintf(out_stream,"Scratch directory (save_dir):     %s\n",save_dir.c_str());
    fprintf(out_stream,"Memory (block2 double stack):     %g GB\n",memory);
    if(warm_start==DMRG_WARM_ON){
        fprintf(out_stream,"MPS warm-start:                   on (after %d cold iter)\n",warm_start_after);
        fprintf(out_stream,"Warm re-solve sweeps:             %d\n",warm_sweeps);
        fprintf(out_stream,"Rotate reused MPS:                %s\n",warm_rotate==DMRG_WARM_ON?"yes":"no (reuse-only)");
        fprintf(out_stream,"MPS-rotation bond dim (rot_m):    %d\n",rot_m==0?m:rot_m);
        if(warm_rotate==DMRG_WARM_ON)
            fprintf(out_stream,"MPS-rotation TE steps (rot_steps):%d\n",rot_steps);
    }
    else
        fprintf(out_stream,"MPS warm-start:                   off\n");
    if(print_dets==DMRG_WARM_ON){
        fprintf(out_stream,"Leading determinants:             on\n");
        fprintf(out_stream,"  read-out rotation bond dim:     %d (%d TE steps)\n",det_rot_m,det_rot_steps);
        if(extract_m>0)
            fprintf(out_stream,"  extraction compression:         m=%d, cutoff=%e\n",extract_m,extract_cutoff);
        else
            fprintf(out_stream,"  extraction:                     no compression, cutoff=%e\n",extract_cutoff);
    }
    else
        fprintf(out_stream,"Leading determinants:             off\n");
    fprintf(out_stream,"\n");

    return 0;
}

dmrg_par::~dmrg_par(){

}

//XMC
xmc_par::xmc_par(){
    y=0;
    edshift=0.02;
    have_ifitd=0;
    d_only=0;
    pt1_d=1;
//     mult=0;
    gamess_file_name = NULL;
    
    n_fit    =XMC_N_FIT    ;
    n_fit_pol=XMC_N_FIT_POL;
    
}

int xmc_par::read_group(char * inp, cas_par * ext_cas){
    
    cas = ext_cas;
    
    recursive_file P;
    char line[BUF_LINE_LENGTH];
    
    P.r_open(inp);
        
    P.r_gets(line,BUF_LINE_LENGTH);;
    while((key_word_comp(line, xmc_group_start)==0)&&(!P.r_eof()))P.r_gets(line,BUF_LINE_LENGTH);;
    if(key_word_comp(line, xmc_group_start)==0){
        fprintf(out_stream,"ERROR: could not find $XMC group in file %s\n",inp);
        exit(1);
    }
    
    while(!P.r_eof()){
        read_line(line);
        if(key_word_comp(line, xmc_group_end))break;
        P.r_gets(line,BUF_LINE_LENGTH);;
    }
    if(d_only)if(gamess_file_name==NULL){
        fprintf(out_stream,"ERROR: you must use INP=... with D_ONLY=1\n");
        exit(0);
    }
    if(d_only)if(gamess_file_name!=NULL){
        FILE * XMCFile=fopen(gamess_file_name,"r");
        if(XMCFile==NULL){
            fprintf(out_stream,"ERROR: xmc_par::read_group couldn't open file with xmc data for reading energy\n");
            fprintf(out_stream,"       file: %s\n",gamess_file_name);
            exit(1);
        }
        fclose(XMCFile);
    
    }
    

    if(have_ifitd){
//         fprintf(out_stream,"%s\n",gamess_file_name);
//         getchar();
        if(gamess_file_name==NULL){
            fprintf(out_stream,"ERROR: you must use INP=... with IFITD=1\n");
            exit(0);
        }
        ifitd_energy = xmc_ifitd_read(gamess_file_name);
    }
    
    
        
    return 0;
}

int xmc_par::read_line(char * inp){
    
    if(key_word_comp(inp, edshift_kw))
        edshift = kw_to_f(inp, edshift_kw, 0.02);
    
    if(key_word_comp(inp, xmc_d_only_kw)){
        d_only = kw_to_i(inp, xmc_d_only_kw, 0);
    }
    
    if(key_word_comp(inp, xmc_n_fit_pol_kw)){
        n_fit_pol = kw_to_i(inp, xmc_n_fit_pol_kw, XMC_N_FIT_POL);
    }
    
    if(key_word_comp(inp, xmc_n_fit_kw)){
        n_fit = kw_to_i(inp, xmc_n_fit_kw, XMC_N_FIT);
    }
        
    if(key_word_comp(inp, inp_kw)){
        kw_to_s(&gamess_file_name, inp, inp_kw);
//         fprintf(out_stream,"%s\n",gamess_file_name);
    }
    
    
//     fprintf(out_stream,"ed = %e\n",edshift);
//     getchar();
    
    if(key_word_comp(inp, ifitd_kw)){
        have_ifitd=kw_to_i(inp, ifitd_kw, 0);
    }
    
    if(key_word_comp(inp, pt1_dipole_kw))
        pt1_d = kw_to_i(inp, pt1_dipole_kw,1);

//     getchar();



    return 0;
}

int xmc_par::write_info(int n_a, int n_b, int n_o, int mult){
    
    int s;
    fprintf(out_stream,"\n\n\n");
    fprintf(out_stream,"_____________________Starting_XMCQDPT2_calculations____________________\n\n");
    
    fprintf(out_stream,"Number of states                  %d\n",cas->n_s );
    fprintf(out_stream,"State multiplicity                %d " ,mult);
    if(mult==1) fprintf(out_stream,"(singlet)\n");
    if(mult==2) fprintf(out_stream,"(doublet)\n");
    if(mult==3) fprintf(out_stream,"(triplet)\n");
    if(mult==4) fprintf(out_stream,"(quadruplet)\n");
    if(mult >4) fprintf(out_stream,"(high-spin)\n");
    if(mult==0) fprintf(out_stream,"(all)\n");
    fprintf(out_stream,"\n");
    fprintf(out_stream,"Active space:\n");
    fprintf(out_stream,"Number of alpha electrons         %d\n",n_a);
    fprintf(out_stream,"Number of beta  electrons         %d\n",n_b);
    fprintf(out_stream,"Number of active orbitals         %d\n",n_o);
    fprintf(out_stream,"\n");
    fprintf(out_stream,"State averaging:\n");
    fprintf(out_stream,"Weights:\t");
    for(auto &w: cas->w_state) fprintf(out_stream,"% .2e ",w);
    fprintf(out_stream,"\n");
    fprintf(out_stream,"\n");
    fprintf(out_stream,"Energy Denominator SHIFT          %e (a.u.)\n", edshift);
    fprintf(out_stream,"\n");
    fprintf(out_stream,"Resolvent fitting:\n");
    fprintf(out_stream,"grid size                         %d\n", n_fit);
    fprintf(out_stream,"polynomial degree                 %d\n", n_fit_pol*2+1);
    fprintf(out_stream,"\n");
    fprintf(out_stream,"Additional options:\n");
    fprintf(out_stream,"skip energy calculation an read from file        -     %c%c%c\n",ny[3*d_only]
                                                                            ,ny[3*d_only+1]
                                                                            ,ny[3*d_only+2]);
    fprintf(out_stream,"read fitting energies from file                  -     %c%c%c\n",ny[3*have_ifitd]
                                                                            ,ny[3*have_ifitd+1]
                                                                            ,ny[3*have_ifitd+2]);
    if((have_ifitd!=0)||(d_only!=0))
        fprintf(out_stream,"File with external data: %s\n",gamess_file_name);
//     dav.write_info();
    fprintf(out_stream,"\nPT first order term: %s\n",pt1_d?"yes":"no");
    fprintf(out_stream,"_______________________________________________________________________\n\n\n");
    
    
    return 0;
    
}


xmc_par::~xmc_par(){
    
    if(gamess_file_name != NULL)
        delete[] gamess_file_name;
        
}

cdas_par::cdas_par(){
    
    y=0;
    edshift=0.00;
    have_eps=0;
    gamess_file_name=NULL;
    IPEA=0;
    MPPT=0;
    HOMO=0;
    actual=0;
    sing_e=0;
    mult_e=0;
    orb_e=0;
    fit_e=0;
    rotate_orbs=1;//default - perform rotation
    pt1_d=1;

}

int cdas_par::read_group(char * inp, cas_par * ext_cas){
    
    cas=ext_cas;
    
    recursive_file P;
    char line[BUF_LINE_LENGTH];
    
    P.r_open(inp);
        
    P.r_gets(line,BUF_LINE_LENGTH);;
    while((key_word_comp(line, cdas_group_start)==0)&&(!P.r_eof()))P.r_gets(line,BUF_LINE_LENGTH);;
    if(key_word_comp(line, cdas_group_start)==0){
        fprintf(out_stream,"ERROR: could not find $CDAS group in file %s\n",inp);
        exit(1);
    }
    
    while(!P.r_eof()){
        read_line(line);
        if(key_word_comp(line, cdas_group_end))break;
        P.r_gets(line,BUF_LINE_LENGTH);;
    }
    
    if((IPEA+MPPT+HOMO+actual+sing_e+mult_e+orb_e+fit_e)!=1){
        fprintf(out_stream,"ERROR: incorrect choice of energy scheme\n");
        fprintf(out_stream,"       IPEA                    - %c%c%c\n",ny[3*IPEA  ],ny[3*IPEA  +1],ny[3*IPEA  +2]);
        fprintf(out_stream,"       MPPT                    - %c%c%c\n",ny[3*MPPT  ],ny[3*MPPT  +1],ny[3*MPPT  +2]);
        fprintf(out_stream,"       HOMO                    - %c%c%c\n",ny[3*HOMO  ],ny[3*HOMO  +1],ny[3*HOMO  +2]);
        fprintf(out_stream,"       ACTUAL                  - %c%c%c\n",ny[3*actual],ny[3*actual+1],ny[3*actual+2]);
        fprintf(out_stream,"       ENERGY                  - %c%c%c\n",ny[3*sing_e],ny[3*sing_e+1],ny[3*sing_e+2]);
        fprintf(out_stream,"       ENERGIES                - %c%c%c\n",ny[3*mult_e],ny[3*mult_e+1],ny[3*mult_e+2]);
        fprintf(out_stream,"       USE_ORB_FOR_ENERGY      - %c%c%c\n",ny[3*orb_e ],ny[3*orb_e +1],ny[3*orb_e +2]);
        fprintf(out_stream,"       USE_FIREFLY_FIT_ENERGY  - %c%c%c\n",ny[3*fit_e ],ny[3*fit_e +1],ny[3*fit_e +2]);
        fprintf(out_stream,"       you must choose only one.\n");
        exit(1);
    }
        
    
    
    return 0;
}

int cdas_par::read_line(char * inp){
        
    if(key_word_comp(inp, edshift_kw))
        edshift = kw_to_f(inp, edshift_kw, 0.02);
    
    if(key_word_comp(inp, cdas_ipea_kw)){
        IPEA=1;
    }
    
    if(key_word_comp(inp, cdas_mppt_kw)){
        MPPT=1;
    }
    
    if(key_word_comp(inp, cdas_homo_kw)){
        HOMO=1;
    }
    
    if(key_word_comp(inp, cdas_actual_kw)){
        actual=1;
    }
    
    if(key_word_comp(inp, cdas_sing_en_kw)){
        sing_e=1;
        have_eps=0;
        eps.resize(1);
        eps[0]=kw_to_f(inp, cdas_sing_en_kw, 0.0);
    }
    
    if(key_word_comp(inp, cdas_mult_en_kw)){
        mult_e=1;
        have_eps=1;
        int n=kw_count(inp, cdas_mult_en_kw,';');
        kw_to_d_v(&eps, inp, cdas_mult_en_kw, n);
    }
    if(key_word_comp(inp, cdas_orb_en_kw)){
        orb_e=1;
        n_orb=kw_to_i(inp, cdas_orb_en_kw, -1);
    }
    
    if(key_word_comp(inp, cdas_fit_en_kw)){
        fit_e=1;
        have_eps=1;
        kw_to_s(&gamess_file_name, inp, cdas_fit_en_kw);
        eps=xmc_ifitd_read(gamess_file_name);
    }
    
    if(key_word_comp(inp, cdas_rot_orbs))
        rotate_orbs = kw_to_i(inp, cdas_rot_orbs,1);
    
    if(key_word_comp(inp, pt1_dipole_kw))
        pt1_d = kw_to_i(inp, pt1_dipole_kw,1);
    
    
    
    return 0;
}

int cdas_par::write_info(int n_a, int n_b, int n_o, int mult){
    
    int s;
    fprintf(out_stream,"\n\n\n");
    fprintf(out_stream,"_____________________Starting_CDAS_PT2_calculations____________________\n\n");
    
    fprintf(out_stream,"Number of states                  %d\n",cas->n_s );
    fprintf(out_stream,"State multiplicity                %d " ,mult);
    if(mult==1) fprintf(out_stream,"(singlet)\n");
    if(mult==2) fprintf(out_stream,"(doublet)\n");
    if(mult==3) fprintf(out_stream,"(triplet)\n");
    if(mult==4) fprintf(out_stream,"(quadruplet)\n");
    if(mult >4) fprintf(out_stream,"(high-spin)\n");
    if(mult==0) fprintf(out_stream,"(all)\n");
    fprintf(out_stream,"\n");
    fprintf(out_stream,"Active space:\n");
    fprintf(out_stream,"Number of alpha electrons         %d\n",n_a);
    fprintf(out_stream,"Number of beta  electrons         %d\n",n_b);
    fprintf(out_stream,"Number of active orbitals         %d\n",n_o);
    fprintf(out_stream,"\n");
    fprintf(out_stream,"State averaging:\n");
    fprintf(out_stream,"Weights:\t");
    for(auto &w: cas->w_state) fprintf(out_stream,"% .2e ",w);
    fprintf(out_stream,"\n");
    fprintf(out_stream,"\n");
    fprintf(out_stream,"Energy Denominator SHIFT          %e (a.u.)\n", edshift);
    fprintf(out_stream,"\n");
    fprintf(out_stream,"Energy scheme                                    -     ");
    if(IPEA  )fprintf(out_stream,"IPEA\n");
    if(MPPT  )fprintf(out_stream,"MPPT\n");
    if(HOMO  )fprintf(out_stream,"HOMO\n");
    if(actual)fprintf(out_stream,"ACTUAL\n");
    if(sing_e)fprintf(out_stream,"ENERGY\n");
    if(mult_e)fprintf(out_stream,"ENERGIES\n");
    if(orb_e )fprintf(out_stream,"USE_ORB_FOR_ENERGY\n");
    if(fit_e )fprintf(out_stream,"USE_FIREFLY_FIT_ENERGY\n");
    if(fit_e!=0)
        fprintf(out_stream,"File with external data: %s\n",gamess_file_name);
    // fprintf(out_stream,"\nOrbital rotation: ");
         // if(rotate_orbs==1)fprintf(out_stream,"1 - use semi-canonical orbitals\n");
    // else if(rotate_orbs==0)fprintf(out_stream,"0 - use initial orbitals\n");
    // else if(rotate_orbs==2)fprintf(out_stream,"2 - use semi-canonic inactive and initial active orbitals\n");
    // else{
        // nopt_printf("\n\nERROR: invalid value of cdas.rotate_orbs=%d\n", rotate_orbs);
        // exit(0);
    // }
    fprintf(out_stream,"\nPT first order term: %s\n",pt1_d?"yes":"no");
    

    
//     dav.write_info();
    
    fprintf(out_stream,"_______________________________________________________________________\n\n\n");
    
    
    return 0;
    
}


cdas_par::~cdas_par(){
    
    if(gamess_file_name!=NULL)delete[] gamess_file_name;
    
}

int print_dipole     = PRINT_DIPOLE;
int print_dispersion = PRINT_DISP;
int print_quadrupole = PRINT_QUADRUPOLE;
int print_mulliken   = PRINT_MULLIKEN;

int read_prop_line(char * inp){
    
    if(key_word_comp(inp, print_dipole_kw))
        print_dipole = kw_to_i(inp, print_dipole_kw, 1);
    
    if(key_word_comp(inp, print_dispersion_kw))
        print_dispersion = kw_to_i(inp, print_dispersion_kw, 0);
    
    if(key_word_comp(inp, print_quadrupole_kw))
        print_quadrupole = kw_to_i(inp, print_quadrupole_kw, 0);
    
    if(key_word_comp(inp, print_mulliken_kw))
        print_mulliken = kw_to_i(inp, print_mulliken_kw, 0);
    
    return 0;
}

int read_prop_group(char * inp){
    
    recursive_file P;
    char line[BUF_LINE_LENGTH];
    
    P.r_open(inp);
    if(P.file[0]==NULL){
        fprintf(out_stream,"ERROR: could not open %s\n",inp);
        exit(1);
    }
//     fprintf(out_stream,"could open %s\n",inp);
        
    P.r_gets(line,BUF_LINE_LENGTH);;
    while((key_word_comp(line, prop_group)==0)&&(!P.r_eof()))P.r_gets(line,BUF_LINE_LENGTH);;
    if(key_word_comp(line, prop_group)==0){
        return 0;
    }
    
    
    while(true){
        read_prop_line(line);
        if(key_word_comp(line, prop_group_end))break;
        if(P.r_eof()) break;
        P.r_gets(line,BUF_LINE_LENGTH);;
    }
    
//     
//     if(P.r_eof()){
//         return 0;
//     }
    
    
    
    return 0;
}

inp_par::inp_par(){
    inp_name = new char[BUF_LINE_LENGTH];
    ecp=0;
    num_threads=1;
    doci_dec=0;
    out_folder_defined=0;
    out_folder_name=NULL;
    point_group=NULL;
    
    job_name = NULL;
    
    
}

int inp_par::write_info(){
    
    int s;
    fprintf(out_stream,"\n\n\n");
    fprintf(out_stream,"_____________________________Input_summary_____________________________\n\n");
    
    fprintf(out_stream,"\n");
    fprintf(out_stream,"Number of threads: % d\n",num_threads);
    fprintf(out_stream,"\n");
    fprintf(out_stream,"Job name: %s\n",job_name);
    fprintf(out_stream,"\n");
    fprintf(out_stream,"Methods:\n");
    fprintf(out_stream,"Restricted Hartree-Fock                          -     %c%c%c\n",ny[3*rhf.y],ny[3*rhf.y+1],ny[3*rhf.y+2]);
    fprintf(out_stream,"Complete Active Space SCF                        -     %c%c%c\n",ny[3*cas.y],ny[3*cas.y+1],ny[3*cas.y+2]);
    fprintf(out_stream,"eXtended MultiConfiguration QuasiDegenerate PT   -     %c%c%c\n",ny[3*xmc.y],ny[3*xmc.y+1],ny[3*xmc.y+2]);
    fprintf(out_stream,"Complete Degenerate Active Space PT              -     %c%c%c\n",ny[3*cdas.y],ny[3*cdas.y+1],ny[3*cdas.y+2]);
    fprintf(out_stream,"Settings:\n");
    fprintf(out_stream,"D5                                               -     %c%c%c\n",ny[3*D5],ny[3*D5+1],ny[3*D5+2]);
    fprintf(out_stream,"Resolution of Identity                           -     %c%c%c\n",ny[3*RI],ny[3*RI+1],ny[3*RI+2]);
    fprintf(out_stream,"Spin-orbit                                       -     %c%c%c\n",ny[3*SO],ny[3*SO+1],ny[3*SO+2]);
    fprintf(out_stream,"Symmetry point group                             -     ");if(IS_SYM)fprintf(out_stream,"%s\n",point_group);
                                                                      else      fprintf(out_stream,"no(C1)\n");
    fprintf(out_stream,"Printing:\n");
    fprintf(out_stream,"dipole moment                                    -     %c%c%c\n",ny[3*print_dipole    ],ny[3*print_dipole    +1],ny[3*print_dipole    +2]);
    fprintf(out_stream,"electron dispersion                              -     %c%c%c\n",ny[3*print_dispersion],ny[3*print_dispersion+1],ny[3*print_dispersion+2]);
    fprintf(out_stream,"quadrupole moment                                -     %c%c%c\n",ny[3*print_quadrupole],ny[3*print_quadrupole+1],ny[3*print_quadrupole+2]);
    fprintf(out_stream,"Mulliken charges                                 -     %c%c%c\n",ny[3*print_mulliken  ],ny[3*print_mulliken  +1],ny[3*print_mulliken  +2]);
    
    if(backup.h1!=1){
    fprintf(out_stream,"\n");
    fprintf(out_stream,"Backup settings:\n");
    fprintf(out_stream,"- 1el data: "); backup_code_print(backup.h1);
//     fprintf(out_stream,"- 2el data: "); backup_code_print(backup.h2);
//     fprintf(out_stream,"- RI  data: "); backup_code_print(backup.ri);
    fprintf(out_stream,"- prefix  : %s\n", backup.prefix);
    }

    
    fprintf(out_stream,"_______________________________________________________________________\n\n\n");
    
//     exit(0);
    return 0;
    
}

int inp_par::read_p_group(char * inp){
    
    recursive_file P;
    char line[BUF_LINE_LENGTH];
    
    P.r_open(inp);
    if(P.file[0]==NULL){
        fprintf(out_stream,"ERROR: could not open %s\n",inp);
        exit(1);
    }
//     fprintf(out_stream,"could open %s\n",inp);
        
    P.r_gets(line,BUF_LINE_LENGTH);
    while((key_word_comp(line, par_group)==0)&&(!P.r_eof()))P.r_gets(line,BUF_LINE_LENGTH);;
    if(P.r_eof()){
        fprintf(out_stream,"ERROR: could not find $PAR group in file %s\n",inp);
        exit(1);
    }
    
    while(!P.r_eof()){
        read_p_line(line);
        if(key_word_comp(line, par_group_end))break;
        P.r_gets(line,BUF_LINE_LENGTH);;
    }
    
    if(job_name==NULL){
        job_name = new char[BUF_LINE_LENGTH];
        sprintf(job_name,"%s\0",DEFAULT_JOB_NAME);
    }
    
    
    return 0;
}

std::vector<int> inp_par::count_groups(char * inp, std::vector<const char *>start_kw, std::vector<const char *>stop_kw){
    
    
    recursive_file P;
    char line[BUF_LINE_LENGTH];
    
    std::vector<int> out;
    out.resize(0);
    
    P.r_open(inp);
    if(P.file[0]==NULL){
        fprintf(out_stream,"ERROR: could not open %s\n",inp);
        exit(1);
    }
    int line_num=0;
    P.r_gets(line,BUF_LINE_LENGTH);
    while(!P.r_eof()){
        if(key_word_comp(line, start_kw)){
            out.push_back(line_num);
            while(key_word_comp(line, stop_kw)==0){
                P.r_gets(line,BUF_LINE_LENGTH);line_num++;//printf("A:%d - %s\n", line_num,line);
                if(P.r_eof()){
                    fprintf(out_stream,"ERROR: could not find end of the group started at line %d\n",out.back());
                    exit(1);
                }
            }
            P.r_gets(line,BUF_LINE_LENGTH);line_num++;//printf("B:%d - %s\n", line_num,line);
        }
        else{
            P.r_gets(line,BUF_LINE_LENGTH);line_num++;//printf("C:%d - %s\n", line_num,line);
        }
    }
    
//     for(int i=0;i<out.size();i++)
//         nopt_printf("%d\n",out[i]);
    
    
    return out;
}

int inp_par::read_symm_group(char * inp){
    
    recursive_file P;
    char line[BUF_LINE_LENGTH];
    
    P.r_open(inp);
    if(P.file[0]==NULL){
        fprintf(out_stream,"ERROR: could not open %s\n",inp);
        exit(1);
    }
//     fprintf(out_stream,"could open %s\n",inp);
        
    P.r_gets(line,BUF_LINE_LENGTH);;
    while((key_word_comp(line, symm_group_start)==0)&&(!P.r_eof()))P.r_gets(line,BUF_LINE_LENGTH);;
    if(P.r_eof()){
        IS_SYM=0;
        return 0;
    }
    
    while(!P.r_eof()){
        read_symm_line(line);
        if(key_word_comp(line, symm_group_end))break;
        P.r_gets(line,BUF_LINE_LENGTH);;
    }
    
    
    
    
    return 0;
}



int inp_par::read(char * ext_inp){
    
    sprintf(inp_name,"%s",ext_inp);
    
    read_p_group(inp_name);
    
    mol_line       = count_groups(inp_name,      mol_group      ,      mol_group_end);
    
    int mc=0;
    if(cas.y)  mc=1;
    if(xmc.y)  mc=1;
    if(cdas.y) mc=1;
    if((cas.y+xmc.y+cdas.y+rhf.y)==0){
        mc=1;
        nopt_printf("WARNING: no QM method chosen\n");
        nopt_printf("         probably nothing to do\n");
    }
    
    
    if(mc){
        act_space_line = count_groups(inp_name,act_space_group_start,act_space_group_end);
    
        if(mol_line.size()!=act_space_line.size()){
            nopt_printf("ERROR: inconsistent number of\n");
            nopt_printf("       $MOL       groups (%3d) and\n",      mol_line.size());
            nopt_printf("       $ACT_SPACE groups (%3d) \n"   ,act_space_line.size());
            exit(0);
        }
    }
    
    read_prop_group(inp_name);
    
    read_symm_group(inp_name);
    
    backup.read_group(inp_name);
    
    if(out_folder_defined==0){
        out_folder_name=new char[3];
        out_folder_name[0]='.';
        out_folder_name[1]='/';
        out_folder_name[2]='\0';
    }
    
    if(xmc.y) cas.y=1;
    if(cdas.y)cas.y=1;
        
    if(rhf.y)  rhf.read_group(inp_name);
    if(cas.y)  cas.read_group(inp_name);
    if(xmc.y)  xmc.read_group(inp_name,&cas);
    if(cdas.y)cdas.read_group(inp_name,&cas);
    if(cis.y)  cis.read_group(inp_name);

    
    
    
    return 0;
    
}

int inp_par::read_p_line(char * inp){
    
    if(key_word_comp(inp, num_threads_kw))
        num_threads = kw_to_i(inp, num_threads_kw, 1);
    
    
    if(key_word_comp(inp, doci_dec_kw))
        doci_dec=kw_to_i(inp, doci_dec_kw,0);
    
    if(key_word_comp(inp, out_folder_kw)){
        out_folder_defined=1;
        kw_to_s(&out_folder_name, inp, out_folder_kw);
    }
    
    if(key_word_comp(inp, d5_kw))
        D5 = kw_to_i(inp, d5_kw, 0);

    if(key_word_comp(inp, ri_kw)){
        RI = kw_to_i(inp, ri_kw, 0);
        if(RI_flag!=-1){
            if(RI_flag!=RI){
                fprintf(stdout,"ERROR: redeclaration of RI flag\n");
                fprintf(stdout,"       one or more of your .inp files containes both\n");
                fprintf(stdout,"       RI=0 and RI=1\n");
                exit(0);
            }
        }
        else RI_flag=RI;
    }
    
    if(key_word_comp(inp, so_kw))
        SO = kw_to_i(inp, so_kw, 0);
    
    
    if(key_word_comp(inp, rhf_kw))
        rhf.y=kw_to_i(inp, rhf_kw,0);
    
    if(key_word_comp(inp, cas_kw))
        cas.y=kw_to_i(inp, cas_kw,0);
    
    if(key_word_comp(inp, xmc_kw))
        xmc.y=kw_to_i(inp, xmc_kw,0);

    if(key_word_comp(inp, cdas_kw))
        cdas.y=kw_to_i(inp, cdas_kw,0);

    if(key_word_comp(inp, cis_kw))
        cis.y=kw_to_i(inp, cis_kw,0);

    if(key_word_comp(inp, name_kw)){
//         out_folder_defined=1;
        kw_to_s(&(job_name), inp, name_kw);
    }
    
//     fprintf(out_stream,"nt=%d\n",num_threads);
//     char* nt;
//     nt = key_word_find(inp,num_threads_kw);
//     if(nt==NULL)fprintf(out_stream,"not found\n");
//     if(nt!=NULL)fprintf(out_stream,"%d\n",atoi(strstr(nt,"=")+1));
    
    
    return 0;
}

int inp_par::read_symm_line(char * inp){
    
    if(key_word_comp(inp, symm_point_group_kw)){
        kw_to_s(&point_group, inp, symm_point_group_kw);
        IS_SYM=1;
        if(strstr(point_group,"C1"))IS_SYM=0;
        if(key_word_comp(point_group,linear_kw))LINEAR=1;
    }
    
    
    return 0;
}

inp_par::~inp_par(){
    
    if(out_folder_name!=NULL)delete[] out_folder_name;
    if(point_group    !=NULL)delete[] point_group    ;
    if(inp_name       !=NULL)delete[] inp_name       ;
    if(job_name       !=NULL)delete[] job_name       ;

}





std::vector<int> molecule_read_by_inp_par(molecule * M, inp_par * P){
    
    char * report_name = new char[BUF_LINE_LENGTH];
    // sprintf(report_name,"%s_out_stream.txt\0",P->job_name);
    // mol_report=stdout;//fopen(report_name,"w");
    
    std::vector<int> return_info;
    
    int s = P->mol_line.size();
    
    int mc=0;
    if( P->cas.y||P->xmc.y||P->cdas.y)mc=1;
    if((P->cas.y +P->xmc.y +P->cdas.y+P->rhf.y)==0)mc=1;
        
        
    molecule *Mf;
    Mf = new molecule[s];
        
    fprintf(out_stream,"Reading molecule:\n");
    
    
    for(int i=0;i<s;i++){

        fprintf(out_stream,"-----------------------------------------------------------------------\n");
        
        sprintf(Mf[i].source_name,"%s",P->inp_name);
        Mf[i].inp_mol_line       = P->mol_line      [i];
        if(mc) Mf[i].inp_act_space_line = P->act_space_line[i];
        
        
        Mf[i].geom_read();
//         Mf[i].geom_print();
        
        fprintf(out_stream,"                      charge   = %d\n",Mf[i].mol_charge);
        
        Mf[i].basis_name_read();
        fprintf(out_stream,"                      basis    = %s",Mf[i].basis_name);
        if(Mf[i].basis_in_lib)fprintf(out_stream," (library)\n");
        else                  fprintf(out_stream," (current folder)\n");
        
        if((Mf[i].ecp_name!=Mf[i].basis_name)||SO){
            fprintf(out_stream,"                      ecp      = %s",Mf[i].ecp_name);
            if(Mf[i].ecp_in_lib)fprintf(out_stream," (library)\n");
            else                fprintf(out_stream," (current folder)\n");
        }
        
        if(SO){
            fprintf(out_stream,"                      so       = %s",Mf[i].so_name);
            if(Mf[i].so_in_lib)fprintf(out_stream," (library)\n");
            else               fprintf(out_stream," (current folder)\n");
        }
        
        
        if(Mf[i].basis_format==0)
            Mf[i].s=basis_lib_read_gbs(Mf+i,Mf[i].basis_name,
                                       Mf[i].basis_in_lib,Mf[i].basis_start_line,
                                       &(Mf[i].lib_coef),&(Mf[i].shell_center),false, nullptr, nullptr);
        if(Mf[i].basis_format==1)
            Mf[i].s=basis_lib_read_exp(Mf+i,Mf[i].basis_name,
                                       Mf[i].basis_in_lib,Mf[i].basis_start_line,
                                       &(Mf[i].lib_coef),&(Mf[i].shell_center),false);
            
            
        if(Mf[i].ecp_format==0)ecp_lib_read_gbs (Mf+i,Mf[i].ecp_in_lib,Mf[i].ecp_start_line,Mf[i].ecp_name,"ECP");
        if(Mf[i].ecp_format==1)ecp_lib_read_exp (Mf+i,Mf[i].ecp_in_lib,Mf[i].ecp_start_line,Mf[i].ecp_name,"ECP");
        if(Mf[i].ecp_format==2)ecp_lib_read_grpp(Mf+i,Mf[i].ecp_in_lib,Mf[i].ecp_start_line,Mf[i].ecp_name,"ECP");
        if(SO){
            if(Mf[i].so_format==0)ecp_lib_read_gbs (Mf+i,Mf[i].so_in_lib,Mf[i].so_start_line,Mf[i].so_name,"SO");
            if(Mf[i].so_format==1)ecp_lib_read_exp (Mf+i,Mf[i].so_in_lib,Mf[i].so_start_line,Mf[i].so_name,"SO");
            if(Mf[i].so_format==2)ecp_lib_read_grpp(Mf+i,Mf[i].so_in_lib,Mf[i].so_start_line,Mf[i].so_name,"SO");
        }
        
        
        Mf[i].n_ao=0;
        for(auto &s: Mf[i].s)
            Mf[i].n_ao+=s.size();
        Mf[i].n_mo=Mf[i].n_ao;
        if(D5){
            Mf[i].n_mo=0;
            for(auto &s: Mf[i].s)
                Mf[i].n_mo+=s.contr[0].l*2+1;
        }

        return_info.push_back(Mf[i].s.size());
        
        if(RI){
            if(Mf[i].ri_basis_format==0)
                Mf[i].ri_s=basis_lib_read_gbs(Mf+i,Mf[i].ri_basis_name,
                                              Mf[i].ri_basis_in_lib,Mf[i].ri_basis_start_line,
                                              &(Mf[i].lib_coef_ri),&(Mf[i].ri_shell_center), true, nullptr, nullptr);
                                              
            if(Mf[i].ri_basis_format==1)
                Mf[i].ri_s=basis_lib_read_exp(Mf+i,Mf[i].ri_basis_name,
                                              Mf[i].ri_basis_in_lib,Mf[i].ri_basis_start_line,
                                              &(Mf[i].lib_coef_ri),&(Mf[i].ri_shell_center), true);
            fprintf(out_stream,"                      RI basis = %s",Mf[i].ri_basis_name);
            if(Mf[i].ri_basis_in_lib)fprintf(out_stream," (library)\n");
            else                     fprintf(out_stream," (current folder)\n");
        
        }
        
        Mf[i].MO_gen();
        
        if(mc)Mf[i].active_space_read(P->rhf.y,1-P->cas.y);
        else
            Mf[i].STATES_set_zero();
        
    }
    fprintf(out_stream,"_______________________________________________________________________\n\n");

    mol_array_link(M, Mf, s, 1);
    
    if(IS_SYM)M->S.init(P->point_group);
    
    delete[] Mf;
    
    if(LINEAR)M->lin_prepare();
    //fclose(mol_report);
    // mol_report=stdout;
    
    
    M->geom_print();
    
    M->basis_print(0);
    if(RI)M->basis_print(1);
    M->PP_print();
    
    
    
    
    fflush(out_stream);
    
    return return_info;
    
}


