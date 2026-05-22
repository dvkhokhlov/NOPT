#include <cstring>
#include <cstdio>
#include <vector>
#include <stdlib.h>
#include "common_vars.h"

using namespace std;

vector<const char *>rhf_kw{{"rhf"},
                                 {"RHF"}};

vector<const char *>rhf_group_start{{"$rhf"},
                                    {"$RHF"}};

vector<const char *>rhf_group_end{{"$rhfend"},
                                  {"$RHFEND"},
                                  {"$end"},
                                  {"$END"}};

vector<const char *>rhf_huckel{{"huckel"},
                               {"Huckel"},
                               {"HUCKEL"}};
                         
vector<const char *>rhf_h_core{{"h_core"},
                               {"H_CORE"}};
                  
vector<const char *>rhf_read{{"read"},
                             {"READ"},
                             {"moread"},
                             {"MOREAD"},
                             {"readmo"},
                             {"READMO"},
                             {"read_orb"},
                             {"READ_ORB"},
                             {"read_orbs"},
                             {"READ_ORBS"}};                         
                         
vector<const char *>max_it_kw{{"maxit"},
                              {"MAXIT"},
                              {"max_iter"},
                              {"MAX_ITER"},
                              {"max_iters"},
                              {"MAX_ITERS"}};

vector<const char *>e_conv_kw{{"dE"},
                              {"DE"},
                              {"deltaE"},
                              {"DELTAE"},
                              {"delta_E"},
                              {"DELTA_E"}};
                        
vector<const char *>g_conv_kw{{"grad"},
                              {"GRAD"},
                              {"orb_grad"},
                              {"ORB_GRAD"}};                        
                              
vector<const char *>s_conv_kw{{"step"},
                              {"STEP"}};                        
                              
vector<const char *>x_max_kw{{"x_max"},
                             {"X_MAX"}};                        

                              
vector<const char *>r_conv_kw{{"dR"},
                              {"DR"},
                              {"deltaR"},
                              {"DELTAR"},
                              {"delta_R"},
                              {"DELTA_R"}};
                        
vector<const char *>cas_kw{{"cas"},
                           {"CAS"},
                           {"casscf"},
                           {"CASSCF"}};

vector<const char *>cas_SA_kw{{"SA"},
                              {"averaged"},
                              {"AVERAGED"}};
                           
vector<const char *>cas_DA_kw{{"DA"},
                              {"density_averaged"},
                              {"DENSITY_AVERAGED"}};
                           
vector<const char *>cas_SM_kw{{"SM"},
                              {"separete"},
                              {"SEPARATE"}};
                           
                           
vector<const char *>cas_group_start{{"$cas"},
                                    {"$CAS"},
                                    {"$casscf"},
                                    {"$CASSCF"}};                        

vector<const char *>cas_group_end{{"cas_end"},
                                  {"CAS_END"},
                                  {"casend"},
                                  {"CASEND"},
                                  {"casscf_end"},
                                  {"CASSCF_END"},
                                  {"casscfend"},
                                  {"CASSCFEND"},
                                  {"$end"},
                                  {"$END"}};
                        
vector<const char *>num_state_kw{{"n_s"},
                                 {"N_S"},
                                 {"numstate"},
                                 {"NUMSTATE"},
                                 {"num_state"},
                                 {"NUM_STATE"}};

vector<const char *>num_frozen_orb_kw{{"n_f"},
                                      {"N_F"},
                                      {"num_frozen_orb"},
                                      {"NUM_FROZEN_ORB"},
                                      {"num_frozen_orbitals"},
                                      {"NUM_FROZEN_ORBITALS"}};

vector<const char *>method_kw{{"method"},
                              {"METHOD"},
                              {"Method"}};

vector<const char *>nat_orb_kw{{"no"},
			       {"NO"},
			       {"nat_orb"},
                               {"NAT_ORB"},
                               {"natural_orbitals"},
                               {"NATURAL_ORBITALS"}};

vector<const char *>cis_kw{{"cis"},
                           {"CIS"}};

vector<const char *>cis_group_start{{"$cis"},
                                    {"$CIS"}};

vector<const char *>cis_group_end{{"cis_end"},
                                  {"CIS_END"},
                                  {"cisend"},
                                  {"CISEND"},
                                  {"$end"},
                                  {"$END"}};
                           
vector<const char *>w_state_kw{{"wstate"},
                               {"WSTATE"},
                               {"w_state"},
                               {"W_STATE"}};
                               
vector<const char *>w_state_eq_kw{{"all_1"},
                                  {"ALL_1"},
                                  {"eq"},
                                  {"EQ"},
                                  {"equal"},
                                  {"EQUAL"}};
                                  
vector<const char *>w_state_irrep_kw{{"rep"},
                                     {"REP"},
                                     {"ir.rep"},
                                     {"IR.REP"},
                                     {"irrep"},
                                     {"IRREP"}};
                                  
vector<const char *>cas_track_kw{{"track"},
                                 {"TRACK"},
                                 {"tracking"},
                                 {"TRACKING"}};

vector<const char *>cas_rotate_orbs_kw{{"rotate"},
                                       {"ROTATE"}};

vector<const char *>w_state_group_start{{"$state"},
                                        {"$WSTATE"},
                                        {"$w_state"},
                                        {"$W_STATE"}};

vector<const char *>w_state_group_end{{"$stateend"},
                                      {"$WSTATEEND"},
                                      {"$w_state_end"},
                                      {"$W_STATE_END"},
                                      {"$end"},
                                      {"$END"}};

vector<const char *>w_linear_L_kw{{"lambda"},
                                  {"Lambda"},
                                  {"LAMBDA"}};

vector<const char *>w_linear_sign_kw{{"sign"},
                                     {"Sign"},
                                     {"SIGN"}};

                                 
vector<const char *>mult_kw{{"mult"},
                            {"MULT"}};

vector<const char *>print_number_kw{{"p_n"},
                                    {"P_N"},
                                    {"p_num"},
                                    {"P_NUM"},
                                    {"print_number"},
                                    {"PRINT_NUMBER"}};

vector<const char *>reorder_kw{{"reorder"},
                               {"REORDER"},
                               {"norder"},
                               {"NORDER"}};
                         
vector<const char *>orbitals_kw{{"orbitals"},
                                {"ORBITALS"},
                                {"act_orb"},
                                {"ACT_ORB"},
                                {"active"},
                                {"ACTIVE"}};
                         
                              
vector<const char *>dav_group_start{{"$dav"},
                                    {"$DAV"}};

                              
// vector<const char *>max_dav_it_kw{{"dav_it"},
//                                   {"DAV_IT"},
//                                   {"davit"},
//                                   {"DAVIT"}};

vector<const char *>dav_se_kw{{"s_ev"},
                              {"s_eval"},
                              {"S_EV"},
                              {"S_EVAL"}};

vector<const char *>dav_n_bf_kw{{"n_bf"},
                                {"N_BF"}};

vector<const char *>dav_sparsed_HC_kw{{"sparsed_Hc"}};
                          
vector<const char *>dav_group_end{{"dav_end"},
                                  {"DAV_END"},
                                  {"davend"},
                                  {"DAVEND"},
                                  {"$end"},
                                  {"$END"}};
                                  
vector<const char *>MO_group_start{{"$VEC1"},
                                   {"$vec1"},
                                   {"$VEC"},
                                   {"$vec"},
                                   {"$VEC_A_1"},
                                   {"$vec_a_1"},
                                   {"$VEC_A"},
                                   {"$vec_a"},
                                   {"_MO"},
                                   {"_mo"},
                                   {"[MO]"}};
                                   
vector<const char *>MO_text_group_start{{"$VEC1"},
                                        {"$vec1"},
                                        {"$VEC"},
                                        {"$vec"},
                                        {"$VEC_A_1"},
                                        {"$vec_a_1"},
                                        {"$VEC_A"},
                                        {"$vec_a"},
                                        {"[MO]"},};


vector<const char *>MO_group_b_start{{"$VEC_B_1"},
                                     {"$vec_b_1"},
                                     {"$MO_B"},
                                     {"$mo_b"},
                                     {"$VEC_B"},
                                     {"$vec_b"}};

vector<const char *>MO_group_end{{"$vecend1"},
                                 {"$VECEND1"},
                                 {"$vecend"},
                                 {"$VECEND"},
                                 {"$END"},
                                 {"$end"}};

vector<const char *>basis_kw{{"_basis"},
                             {"_BASIS"}};
                       
vector<const char *>basis_end_kw{{"basisend"},
                                 {"BASISEND"},
                                 {"basis_end"},
                                 {"BASIS_END"}};

vector<const char *>ri_basis_kw{{"_ribasis"},
                                {"_RIBASIS"},
                                {"_ri_basis"},
                                {"_RI_BASIS"},
                                {"_auxbasis"},
                                {"_AUXBASIS"},
                                {"_aux_basis"},
                                {"_AUX_BASIS"}};
                                
vector<const char *>auto_aux_kw{{"auto_aux"},
                                {"AUTO_AUX"},
                                {"auto"},
                                {"AUTO"},
                                {"ri_auto"},
                                {"RI_AUTO"}};

vector<const char *>ecp_kw{{"_ecp"},
                           {"_ECP"}};
                                
vector<const char *>so_kw{{"SO"},
                          {"so"},
                          {"spin_orbit"},
                          {"SPIN_ORBIT"}};

vector<const char *>so_lib_kw{{"_SO"},
                              {"_so"},
                              {"_spin_orbit"},
                              {"_SPIN_ORBIT"}};

                          
vector<const char *>ecp_so_end_kw{{"ecp_end"},
                                  {"ECP_END"},
                                  {"so_end"},
                                  {"SO_END"},
                                  {"ecpend"},
                                  {"ECPEND"},
                                  {"soend"},
                                  {"SOEND"}};

vector<const char *>read_basis_kw{{"_fbasis"},
                                  {"_FBASIS"},
                                  {"_f_basis"},
                                  {"_F_BASIS"}};

                          
                          
vector<const char *>geom_group_start{{"$geom1"},
                                     {"$GEOM1"},
                                     {"$geom"},
                                     {"$GEOM"},
                                     {"$geometry"},
                                     {"$GEOMETRY"},
                                     {"_geom1"},
                                     {"_GEOM1"},
                                     {"_geom"},
                                     {"_GEOM"},
                                     {"_geometry"},
                                     {"_GEOMETRY"},
                                     {"_xyz"},
                                     {"_XYZ"},
                                     {"_zmat"},
                                     {"_z-mat"},
                                     {"_Zmat"},
                                     {"_Z-mat"},
                                     {"_ZMAT"},
                                     {"_Z-MAT"},
                                     {"[Atoms]"}};    
                               
vector<const char *>z_mat_kw{{"zmat"},
                             {"z-mat"},
                             {"Zmat"},
                             {"Z-mat"},
                             {"ZMAT"},
                             {"Z-MAT"}};    
                             
vector<const char *>act_space_group_start{{"$act_space"},
                                          {"$ACT_SPACE"},
                                          {"$active_space"},
                                          {"$ACTIVE_SPACE"},
                                          {"$states1"},
                                          {"$STATES1"},
                                          {"$states"},
                                          {"$STATES"}};    
                             
vector<const char *>act_space_group_end{{"$act_space_end"},
                                        {"$ACT_SPACE_END"},
                                        {"$active_space_end"},
                                        {"$ACTIVE_SPACE_END"},
                                        {"$statesend1"},
                                        {"$STATESEND1"},
                                        {"$statesend"},
                                        {"$STATESEND"},
                                        {"$END"},
                                        {"$end"}};

vector<const char *>geom_dim_unit{{"_unit"}, 
                                  {"_units"}, 
                                  {"_UNIT"},
                                  {"_UNITS"}};

                               
vector<const char *>geom_dim_au{{"au"}, 
                                {"AU"}, 
                                {"bohr"},
                                {"BOHR"}};


vector<const char *>geom_dim_nm{{"nm"}, 
                                {"NM"}, 
                                {"nanometer"},
                                {"NANOMETER"}};

vector<const char *>geom_dim_ang{{"ang"}, 
                                 {"ANG"}, 
                                 {"angstrom"},
                                 {"ANGSTROM"}};
                           
vector<const char *>wf_type{{"$wf_type"}, 
                            {"$WF_TYPE"}};

vector<const char *>wf_type_end{{"$wf_type_end"}, 
                                {"$WF_TYPE_END"}};

vector<const char *>uhf{{"uhf"}, 
                        {"UHF"}};

vector<const char *>wide{{"wide"}, 
                         {"WIDE"}};

vector<const char *>energy_group_start{{"$ENERGY1"},
                                       {"$energy1"},
                                       {"$ENERGY"},
                                       {"$energy"},
                                       {"$ENERGY_A_1"},
                                       {"$energy_a_1"},
                                       {"$ENERGY_A"},
                                       {"$energy_a"}};

vector<const char *>energy_group_b_start{{"$ENERGY_B_1"},
                                         {"$energy_b_1"},
                                         {"$ENERGY_B"},
                                         {"$energy_b"}};

vector<const char *>energy_group_end{{"$energy_end"},
                                     {"$ENERGY_END"},
                                     {"$energyend"},
                                     {"$ENERGYEND"},
                                     {"$END"},
                                     {"$end"}};
                               
                               
vector<const char *>par_group{{"$par"},
                              {"$PAR"}};
                    
vector<const char *>mol_group{{"$mol"},
                              {"$MOL"}};
                        
vector<const char *>par_group_end{{"$par_end"},
                                  {"$PAR_END"},
                                  {"$parend"},
                                  {"$PAREND"},
                                  {"$end"},
                                  {"$END"}};
                            
vector<const char *>mol_group_end{{"$molend"},
                                  {"$molend1"},
                                  {"$MOLEND"},
                                  {"$MOLEND1"},
                                  {"$end"},
                                  {"$end1"},
                                  {"$END"},
                                  {"$END1"}};

vector<const char *>num_threads_kw{{"num_threads"},
                                   {"NUM_THREAD"},
                                   {"numthreads"},
                                   {"NUMTHREAD"},
                                   {"n_threads"},
                                   {"N_THREAD"},
                                   {"nthreads"},
                                   {"NTHREAD"}};
                       
vector<const char *>n_mol_kw{{"n_mol"},
                             {"N_MOL"},
                             {"num_mol"},
                             {"NUM_MOL"},
                             {"nmol"},
                             {"NMOL"},
                             {"nummol"},
                             {"NUMMOL"},
                             {"n_mols"},
                             {"N_MOLS"},
                             {"num_mols"},
                             {"NUM_MOLS"},
                             {"nmols"},
                             {"NMOLS"},
                             {"nummols"},
                             {"NUMMOLS"}};
                       
vector<const char *>doci_dec_kw{{"dec"},
                                {"DEC"}};
                       
vector<const char *>out_folder_kw{{"outfolder"},
                                  {"OUTFOLDER"},
                                  {"out_folder"},
                                  {"OUT_FOLDER"}};
                            
vector<const char *>n_frags_kw{{"n_frag"},
                               {"N_FRAG"},
                               {"num_frag"},
                               {"NUM_FRAG"},
                               {"nfrag"},
                               {"NFRAG"},
                               {"numfrag"},
                               {"NUMFRAG"},
                               {"n_frags"},
                               {"N_FRAGS"},
                               {"num_frags"},
                               {"NUM_FRAGS"},
                               {"nfrags"},
                               {"NFRAGS"},
                               {"numfrags"},
                               {"NUMFRAGS"}};
                         
vector<const char *>inp_kw{{"INPS"},
                           {"inps"},
                           {"INP"}};

vector<const char *>charge_kw{{"_charge"},
                              {"_CHARGE"},
                              {"_charges"},
                              {"_CHARGES"}};
                        
vector<const char *>name_kw{{"name"},
                            {"NAME"}};

vector<const char *>occ_group_start{{"$occ1"},
                                    {"$OCC1"},
                                    {"$occ"},
                                    {"$OCC"}};
                              
vector<const char *>occ_group_end{{"$occend1"},
                                  {"$OCCEND1"},
                                  {"$occend"},
                                  {"$OCCEND"},
                                  {"$END"},
                                  {"$end"}};
                      
vector<const char *>xmc_kw{{"xmc"},
                           {"XMC"},
                           {"xmcqdpt"},
                           {"XMCQDPT"}};
                     
vector<const char *>xmc_group_start{{"$xmc"},
                                    {"$XMC"},
                                    {"$xmcqdpt"},
                                    {"$XMCQDPT"}};

vector<const char *>xmc_group_end{{"$xmcend"},
                                  {"$XMCEND"},
                                  {"$xmcqdptend"},
                                  {"$XMCQDPTEND"},
                                  {"$END"},
                                  {"$end"}};
                            
vector<const char *>edshift_kw{{"edshift"},
                               {"EDSHIFT"}};

vector<const char *>avecoe_kw{{"avecoe"},
                              {"AVECOE"}};

vector<const char *>ifitd_kw{{"ifitd"},
                             {"IFITD"}};

vector<const char *>xmc_d_only_kw{{"d_only"},
                                  {"D_ONLY"}};
                            
vector<const char *>xmc_n_fit_kw{{"n_fit"},
                                 {"N_FIT"}};
                            
vector<const char *>xmc_n_fit_pol_kw{{"n_fit_pol"},
                                     {"N_FIT_POL"}};
                            
vector<const char *>cdas_kw{{"cdas"},
                            {"CDAS"},
                            {"cqd"},
                            {"CQD"},
                            {"cqdas"},
                            {"CQDAS"}};
                     
vector<const char *>cdas_group_start{{"$cdas"},
                                     {"$CDAS"},
                                     {"$cqd"},
                                     {"$CQD"},
                                     {"$cqdas"},
                                     {"$CQDAS"}};

vector<const char *>cdas_group_end{{"$cdas_end"},
                                   {"$CDAS_END"},
                                   {"$cqd_end"},
                                   {"$CQD_END"},
                                   {"$cqdas_end"},
                                   {"$CQDAS_END"},
                                   {"$END"},
                                   {"$end"}};
                            
vector<const char *>cdas_ipea_kw{{"IPEA"},
                                 {"ipea"}};
                          
vector<const char *>cdas_mppt_kw{{"MPPT"},
                                 {"mppt"}};
                          
vector<const char *>cdas_homo_kw{{"HOMO"},
                                 {"homo"}};
                          
vector<const char *>cdas_actual_kw{{"ACTUAL"},
                                   {"actual"}};
                            
vector<const char *>cdas_sing_en_kw{{"ENERGY"},
                                    {"energy"}};
                            
vector<const char *>cdas_mult_en_kw{{"ENERGIES"},                            
                                    {"energies"}};
                             
vector<const char *>cdas_orb_en_kw{{"USE_ORB_FOR_ENERGY"},
                                   {"use_orb_for_energy"}};
                            
vector<const char *>cdas_fit_en_kw{{"USE_FIREFLY_FIT_ENERGY"},
                                   {"use_firefly_fit_energy"}};
                                   
vector<const char *>cdas_rot_orbs{{"rotate"},
                                  {"ROTATE"}};
                            
vector<const char *>pt1_dipole_kw{{"pt1_dipole"},
                                  {"PT1_dipole"},
                                  {"PT1_DIPOLE"}};

vector<const char *>d5_kw{{"d5"},
                          {"D5"}};
                            
vector<const char *>ri_kw{{"ri"},
                          {"RI"}};

vector<const char *>n_st_kw{{"n_st"},
                            {"N_ST"},
                            {"num_st"},
                            {"NUM_ST"},
                            {"nst"},
                            {"NST"},
                            {"numst"},
                            {"NUMST"},
                            {"n_sts"},
                            {"N_STS"},
                            {"num_sts"},
                            {"NUM_STS"},
                            {"nsts"},
                            {"NSTS"},
                            {"numsts"},
                            {"NUMSTS"}};
                        
vector<const char *>n_act_kw{{"n_act"},
                             {"N_ACT"},
                             {"nact"},
                             {"NACT"},
                             {"n_val"},
                             {"N_VAL"},
                             {"nval"},
                             {"NVAL"}};
                        
vector<const char *>n_alp_kw{{"n_alp"},
                             {"N_ALP"},
                             {"nalp"},
                             {"NALP"}};                        


vector<const char *>n_bet_kw{{"n_bet"},
                             {"N_BET"},
                             {"nbet"},
                             {"NBET"}};                        

vector<const char *>backup_group_start{{"$backup"},
                                       {"$BACKUP"}};

vector<const char *>backup_group_end{{"$backupend"},
                                     {"$BACKUPEND"},
                                     {"$END"},
                                     {"$end"}};

vector<const char *>h1_backup_kw{{"h1"},
                                 {"H1"}};
                       
vector<const char *>h2_backup_kw{{"h2"},
                                 {"H2"}};
                        
vector<const char *>ri_backup_kw{{"ri"},
                                 {"RI"}};


                        
vector<const char *>prefix_backup_kw{{"prefix"},
                                     {"PREFIX"}};

vector<const char *>n_bfgs_vec_kw{{"n_bfgs_vec"},
                                  {"N_BFGS_VEC"},
                                  {"n_bfgs"},
                                  {"N_BFGS"},
                                  {"num_bfgs_vec"},
                                  {"NUM_BFGS_VEC"},
                                  {"num_bfgs"},
                                  {"NUM_BFGS"}};

vector<const char *>prop_group{{"$prop"},
                               {"$PROP"}};
                         
vector<const char *>prop_group_end{{"$propend"},
                                   {"$PROPEND"},
                                   {"$prop_end"},
                                   {"$PROP_END"},
                                   {"$end"},
                                   {"$END"}};
                             
vector<const char *>print_dipole_kw{{"dipole"},
                                    {"DIPOLE"}};
                         
// vector<const char *>first_order_property_kw{{"dipole1"},
//                                             {"DIPOLE1"}};
//                                             {"CDAS_FO"}};
//                                             {"DIPOLE"}};

vector<const char *>print_dispersion_kw{{"disp"},
                                        {"DISP"}};

vector<const char *>print_quadrupole_kw{{"quad"},
                                        {"QUAD"}};

vector<const char *>print_mulliken_kw{{"mulliken"},
                                      {"MULLIKEN"}};
                                
vector<const char *>symm_group_start{{"$symm"},
                                     {"$SYMM"},
                                     {"$sym"},
                                     {"$SYM"}
};
                                
vector<const char *>symm_group_end{{"$symmend"},
                                   {"$SYMMEND"},
                                   {"$symm_end"},
                                   {"$SYMM_END"},
                                   {"$end"},
                                   {"$END"}};
                             
vector<const char *>symm_point_group_kw{{"group"},
                                        {"GROUP"}};

vector<const char *>linear_kw{{"lin"},
                              {"LIN"},
                              {"linear"},
                              {"LINEAR"}};


vector<const char *>include_kw{{"#include"}};



vector<const char *>lib_kw{{"lib"},
                           {"LIB"},
                           {"library"},
                           {"LIBRARY"}};
                              
                            
int failed_find_symbol(char * inp, vector<const char *> keywords, char s){
    
    fprintf(out_stream,"ERROR: failed reading one of the keywords:\n");
    for(const auto&k:keywords)
    fprintf(out_stream,"       - %s\n",k);
    fprintf(out_stream,"       Could not find \"%c\" in line:\n       %s\n",s,inp);
    exit(0);
    
    return 1;
}

int is_splitter(char a){
    
    if(a==' ')return 1;
    if(a==',')return 1;
    if(a==';')return 1;
    if(a=='\n')return 1;
    if(a=='=')return 1;

//     if(a==" ")return 1;
    return 0;
}


int valid_pointer(char * a, char *b){
    if(a==nullptr) return 0;
    if(a==b) return 1;
    return is_splitter(a[-1]);
}

char *key_word_find(char * str, vector<const char *> keywords){
    
    char *f;
    std::vector<char *> vf;
    vf.resize(0);
    
    for(const auto&k:keywords){
        f=strstr(str,k);
        if(valid_pointer(f,str)){
            if(vf.size()==0) vf.push_back(f);
            else if(f!=vf[0])vf.push_back(f);
        }
    }
    // fprintf(out_stream,"%d\n",vf.size());
    if(vf.size()==1)return vf[0];
    if(vf.size()>1){
        fprintf(out_stream,"ERROR: in line:\n");
        fprintf(out_stream,"       %s\n",str);
        fprintf(out_stream,"       found more than 1 keyword form list:\n");
        for(const auto&k:keywords)
            fprintf(out_stream,"       - %s\n",k);
        exit(1);
    }
    
    
    return nullptr;
    
}


int key_word_comp(char * str, vector<const char *> keywords){
    
    if(key_word_find(str, keywords)==nullptr) return 0;
    else                                      return 1;
    
//     for(const auto&k:keywords)
//         if(strstr(str,k)!=NULL){ 
//             return 1;
//         }
//     
//     return 0;
    
} 


int key_word_eq(char * str, vector<const char *> keywords){
    
    
    while(is_splitter(str[0]))str++;
    char *f;
    
    for(const auto&k:keywords){
        f=strstr(str,k);
        if(f==str) return 1;
    }
    return 0;
    
}


int kw_to_i(char * inp, vector<const char *> keywords, int def_val){
    
    char* f;
    f = key_word_find(inp,keywords);
    if(f!=NULL){
        if(strstr(f,"=")==NULL)failed_find_symbol(inp,keywords,'=');
        return atoi(strstr(f,"=")+1);
    }
    return def_val;
}

float kw_to_f(char * inp, vector<const char *> keywords, float def_val){
    
    char* f;
    f = key_word_find(inp,keywords);
    if(f!=NULL){
        if(strstr(f,"=")==NULL)failed_find_symbol(inp,keywords,'=');            
        return atof(strstr(f,"=")+1);
    }
    return def_val;
}

int lenth(char* s){
//     fprintf(out_stream,"s = %s\n",s);
        
    int i=0;
    int r=0;
    while(is_splitter(s[i]))i++;
    while(is_splitter(s[i])==0){
//         fprintf(out_stream,"%c\n",s[i]);
        i++;
        r++;
    }
    return r;
}

char* kw_cpy(char * out, char * inp, int l){
    
    int i=0;
    while(is_splitter(inp[i]))i++;
    for(int j=0;j<l;i++,j++)
        out[j]=inp[i];
    out[l]='\0';
    
    return (inp+i);
}

int kw_to_kw(char * inp, vector<const char *> keywords, vector<const char *> keywords_2){
    
    char* f;
    char* s;
    f = key_word_find(inp,keywords);
//     if(nt==NULL)fprintf(out_stream,"not found\n");
    if(f!=NULL){
        if(strstr(f,"=")==NULL)failed_find_symbol(inp,keywords,'=');
        s=strstr(f,"=");while(is_splitter(s[0]))s++;
        if(key_word_eq(s,keywords_2))return 1;
    }
    return 0;
}


int kw_to_s(char ** name, char * inp, vector<const char *> keywords){
    
    char* f;
    char* s;
    f = key_word_find(inp,keywords);
    if(f!=NULL){
        if(strstr(f,"=")==NULL)failed_find_symbol(inp,keywords,'=');
        s=strstr(f,"=");while(is_splitter(s[0]))s++;
        int l = lenth(s);
        *name=new char[l+1];
        kw_cpy(*name,s,l);
//         getchar();
    }
    return 1;
}

int kw_to_s_w_def(char ** name, char * inp, vector<const char *> keywords, char * def_val){
    
    char* f;
    char* s;
    int r_val=1;
    f = key_word_find(inp,keywords);
    if(f!=NULL){
        if(strstr(f,"=")==NULL) {s=def_val;r_val=0;}
        else s=strstr(f,"=");while(is_splitter(s[0]))s++;
        
        int l = lenth(s);
        *name=new char[l+1];
        kw_cpy(*name,s,l);
//         getchar();
    }
    return r_val;
}

int kw_to_s_v(std::vector<char*>* name, char * inp, vector<const char *> keywords, int n){
    
    char* f;
    char* s;
    char* s_tmp;
    f = key_word_find(inp,keywords);
    if(f!=NULL){
        s=strstr(f,"=");while(is_splitter(s[0]))s++;
        for(int i=0;i<n;i++){
            int l = lenth(s);
            (*name)[i]=new char[l+1];
            s_tmp = kw_cpy((*name)[i],s,l);
            s=s_tmp;
        }

//         getchar();
    }
    return 1;
}

int kw_to_i_v(std::vector<int>* c, char * inp, vector<const char *> keywords, int n){
    
    std::vector<char*> c_char;
    c_char.resize(n);
    kw_to_s_v(&c_char,inp,keywords,n);
    for(int i=0;i<n;i++){
        (*c)[i]=atoi(c_char[i]);
        delete[]c_char[i];
    }
    
    return 1;
}

int kw_to_d_v(std::vector<double>* c, char * inp, vector<const char *> keywords, int n){
    
    c->resize(n);
    std::vector<char*> c_char;
    c_char.resize(n);
    kw_to_s_v(&c_char,inp,keywords,n);
    for(int i=0;i<n;i++){
        (*c)[i]=atof(c_char[i]);
        delete[]c_char[i];
    }
    
    return 1;
}


int kw_count(char * inp, vector<const char *> keywords, char end){
    char* f;
    char* s;
    char* s_tmp;
    int c=0;
    f = key_word_find(inp,keywords);
    if(f!=NULL){
        c=1;
        if(strchr(f,end)==NULL)failed_find_symbol(inp,keywords,end);
        if(strstr(f,"=")==NULL)failed_find_symbol(inp,keywords,'=');
        s=strstr(f,"=");while(is_splitter(s[0]))s++;
        int i=0;
        while(s[i]!=end){
            if(is_splitter(s[i]))c++;
            i++;
        }
    }
    
    return c;
}
