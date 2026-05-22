//C++ headers
// # include <vector>

//user
# include "chem_data.h"
# include "molecule.h"
# include "matr.h"
# include "common_vars.h"

int * eq_atoms_A;
int * eq_atoms_B;


int STATES_link(molecule * AB, molecule * A, molecule * B){
    
    if(nQM){
        
        AB->n_frag         = A->n_frag        + B->n_frag;
        AB->n_cor_orb = A->n_cor_orb + B->n_cor_orb;
        
        AB->n_states     = new int[AB->n_frag];
        AB->n_act_el_alp = new int[AB->n_frag];
        AB->n_act_el_bet = new int[AB->n_frag];
        AB->n_act_orb    = new int[AB->n_frag];
        AB->n_cor_orb_f  = new int[AB->n_frag];
        AB->CI           = new aldet_data[AB->n_frag];
        
//         AB->csf_occup[0]    = new std::vector<std::vector<std::vector<int>>>  [AB->n_frag];
//         AB->csf_coef [0]    = new std::vector<std::vector<double>> [AB->n_frag];
        
        for(int i_f=0; i_f<A->n_frag;i_f++){
            AB->n_states[i_f]     = A->n_states[i_f]    ;
            AB->n_act_el_alp[i_f] = A->n_act_el_alp[i_f];
            AB->n_act_el_bet[i_f] = A->n_act_el_bet[i_f];
            AB->n_act_orb[i_f]    = A->n_act_orb[i_f]   ;
            AB->n_cor_orb_f[i_f]  = A->n_cor_orb_f[i_f] ;
            
            aldet_copy(AB->CI+i_f,A->CI+i_f);
            
        }
        
        for(int i_f=0; i_f<B->n_frag;i_f++){
            int i_fp = i_f + A->n_frag;
            AB->n_act_el_alp[i_fp] = B->n_act_el_alp[i_f];
            AB->n_act_el_bet[i_fp] = B->n_act_el_bet[i_f];
            AB->n_act_orb[i_fp]    = B->n_act_orb[i_f]   ;
            AB->n_cor_orb_f[i_fp]  = B->n_cor_orb_f[i_f]   ;
            AB->n_states[i_fp]     = B->n_states[i_f]    ;
            
            aldet_copy(AB->CI+i_fp,B->CI+i_f);
        }
    }
    else{
        
        AB->n_frag         = 1;
        AB->n_cor_orb = A->n_cor_orb + B->n_cor_orb;
        
        AB->n_states     = new int[1];
        AB->n_act_el_alp = new int[1];
        AB->n_act_el_bet = new int[1];
        AB->n_act_orb    = new int[1];
        AB->n_cor_orb_f  = new int[1];
        AB->CI           = new aldet_data[1];
        AB->n_states    [0]=0;
        AB->n_act_el_alp[0]=0;
        AB->n_act_el_bet[0]=0;
        AB->n_act_orb   [0]=0;
        AB->n_cor_orb_f [0]=0;
        
        
        for(int i_f=0; i_f<A->n_frag;i_f++){
            AB->n_act_el_alp[0] += A->n_act_el_alp[i_f];
            AB->n_act_el_bet[0] += A->n_act_el_bet[i_f];
            AB->n_act_orb   [0] += A->n_act_orb[i_f]   ;
            AB->n_cor_orb_f [0] += A->n_cor_orb_f[i_f] ;
            AB->n_states    [0] += A->n_states[i_f]    ;
        }
        
        for(int i_f=0; i_f<B->n_frag;i_f++){
//             int i_fp = i_f + A->n_frag;
            AB->n_act_el_alp[0] += B->n_act_el_alp[i_f];
            AB->n_act_el_bet[0] += B->n_act_el_bet[i_f];
            AB->n_act_orb   [0] += B->n_act_orb[i_f]   ;
            AB->n_cor_orb_f [0] += B->n_cor_orb_f[i_f]   ;
            AB->n_states    [0] += B->n_states[i_f]    ;
        }
        
        AB->CI[0].get_dim(AB->n_act_orb[0],
                          AB->n_act_el_alp[0], 
                          AB->n_act_el_bet[0], 2, 
                          0,//mult, 
                          10//print_number
                          );//two sets are needed for CASSCF
        AB->CI[0].n_states[0]=0;
        AB->CI[0].n_states[1]=0;
    
        
    }
  
    return 0;
    
}

int MO_link(double * AB, double * A, double * B, int dimA, int dimB, 
                                                 int corA, int corB, 
                                                 int actA, int actB, 
                                                 int virtA, int virtB){
    
//     printf("dim: %d + %d\n",dimA, dimB);
//     printf("cor: %d + %d\n",corA, corB);
//     printf("act: %d + %d\n",actA, actB);
//     printf("virt: %d + %d\n",virtA, virtB);
    
//     getchar();
    
    
    
    set_zero_matr(AB,(dimA+dimB)*(dimA+dimB));
    
    int voidA = dimA-corA-actA-virtA;
    int voidB = dimB-corB-actB-virtB;
    
    //corA
    for(int i=0; i<corA; i++)
        for(int j=0; j<dimA; j++)
            AB[i*(dimA+dimB)+j]=A[i*dimA+j];
    //corB
    for(int i=0; i<corB; i++)
        for(int j=0; j<dimB; j++)
            AB[(i+corA)*(dimA+dimB)+j+dimA]=B[i*dimB+j];
    //actA
    for(int i=0; i<actA; i++)
        for(int j=0; j<dimA; j++)
            AB[(i+corA+corB)*(dimA+dimB)+j]=A[(i+corA)*dimA+j];
    //actB
    for(int i=0; i<actB; i++)
        for(int j=0; j<dimB; j++)
            AB[(i+corA+corB+actA)*(dimA+dimB)+j+dimA]=B[(i+corB)*dimB+j];
    //virtA
    for(int i=0; i<virtA; i++)
        for(int j=0; j<dimA; j++)
            AB[(i+corA+corB+actA+actB)*(dimA+dimB)+j]=A[(i+corA+actA)*dimA+j];
    //virtB
    for(int i=0; i<virtB; i++)
        for(int j=0; j<dimB; j++)
            AB[(i+corA+corB+actA+actB+virtA)*(dimA+dimB)+j+dimA]=B[(i+corB+actB)*dimB+j];
    
    //voidA
    for(int i=0; i<voidA; i++)
        for(int j=0; j<dimA; j++)
            AB[(i+corA+corB+actA+actB+virtA+virtB)*(dimA+dimB)+j]=A[(i+corA+actA+virtA)*dimA+j];
    //voidB
    for(int i=0; i<voidB; i++)
        for(int j=0; j<dimB; j++)
            AB[(i+corA+corB+actA+actB+virtA+virtB+voidA)*(dimA+dimB)+j+dimA]=B[(i+corB+actB+virtB)*dimB+j];
        
    
//     printf("AB_VEC\n");
//     PrintMatr(AB,(dimA+dimB),(dimA+dimB),1);
        
    return 0;
    
}

int orb_E_link(double * AB, double * A, double * B, int dimA, int dimB, 
                                                    int corA, int corB, 
                                                    int actA, int actB, 
                                                    int virtA, int virtB){
    
    set_zero_matr(AB,dimA+dimB);
    
    int voidA = dimA-corA-actA-virtA;
    int voidB = dimB-corB-actB-virtB;
    
    //corA
    for(int i=0; i<corA; i++)
        AB[i]=A[i];
    //corB
    for(int i=0; i<corB; i++)
        AB[i+corA]=B[i];
    //actA
    for(int i=0; i<actA; i++)
        AB[i+corA+corB]=A[i+corA];
    //actB
    for(int i=0; i<actB; i++)
        AB[i+corA+corB+actA]=B[i+corB];
    //virtA
    for(int i=0; i<virtA; i++)
        AB[i+corA+corB+actA+actB]=A[i+corA+actA];
    //virtB
    for(int i=0; i<virtB; i++)
        AB[i+corA+corB+actA+actB+virtA]=B[i+corB+actB];
    
    //voidA
    for(int i=0; i<voidA; i++)
        AB[i+corA+corB+actA+actB+virtA+virtB]=A[i+corA+actA+virtA];
    //voidB
    for(int i=0; i<voidB; i++)
        AB[i+corA+corB+actA+actB+virtA+virtB+voidA]=B[i+corB+actB+virtB];
        
    return 0;
    
}



int MO_link_2(double * AB, double * A, double * B, int dimA, int dimB, int corA, int corB, int actA, int actB){
    
    set_zero_matr(AB,(dimA+dimB)*(dimA+dimB));
    
    int virtA = dimA-corA-actA;
    int virtB = dimB-corB-actB;
    
    
    for(int i=0; i<dimA; i++)
        for(int j=0; j<dimA; j++)
            AB[i*(dimA+dimB)+j]=A[i*dimA+j];
    
    for(int i=0; i<dimB; i++)
        for(int j=0; j<dimB; j++)
            AB[(i+dimA)*(dimA+dimB)+j+dimA]=B[i*dimB+j];
           
    
//     printf("AB_VEC\n");
//     PrintMatr(AB,(dimA+dimB),(dimA+dimB),1);
        
    return 0;
    
}

int geom_link(molecule * AB, molecule * A, molecule * B){

  AB->au_coef=1.0; // dim = Bohr!!!!
  
  if(eq_atoms_A==NULL){
      eq_atoms_A = new int[A->n_atoms];
      eq_atoms_B = new int[B->n_atoms];
      AB->n_atoms = A->n_atoms + B->n_atoms;
      printf("START geom_link\n");
      for(int i=0; i<A->n_atoms; i++)eq_atoms_A[i]=-1;
      for(int i=0; i<B->n_atoms; i++)eq_atoms_B[i]=-1;
      for(int i=0; i<A->n_atoms; i++)
      for(int j=0; j<B->n_atoms; j++){
          if((A->atom_coord[3*i  ]-B->atom_coord[3*j  ])*(A->atom_coord[3*i  ]-B->atom_coord[3*j  ])+
             (A->atom_coord[3*i+1]-B->atom_coord[3*j+1])*(A->atom_coord[3*i+1]-B->atom_coord[3*j+1])+
             (A->atom_coord[3*i+2]-B->atom_coord[3*j+2])*(A->atom_coord[3*i+2]-B->atom_coord[3*j+2])<1E-4){
              printf("WARNING %d %d\n",i,j);
              eq_atoms_A[i]=j;
              eq_atoms_B[j]=i;
              AB->n_atoms--;
              
          }
      }
  }
  
  if(AB->n_atoms != (A->n_atoms + B->n_atoms)){
      printf("ERROR: NOPT doesn't support linking of the molecules with common atoms\n");
      exit(0);
  }
  
  
  AB->atom_coord        = new double[AB->n_atoms*3];
  AB->nucl_charges_full = new double[AB->n_atoms];
  AB->nucl_charges_calc = new double[AB->n_atoms];
  AB->n_ecp_electrons   = new int   [AB->n_atoms];
  AB->atom_names        = new char* [AB->n_atoms];
  AB->atom_mass         = new double[AB->n_atoms];
  AB->atom_is_ghost     = new int   [AB->n_atoms];
  
//   for(int i=0;i<AB->n_ion_states;i++)AB->basis[i].n_atoms=AB->n_atoms;
//   for(int i=0;i<AB->n_ion_states;i++)AB->basis[i].coords=AB->atom_coord;
  
  AB->n_el_full = 0;
  AB->n_el_calc = 0;
  
  //copy A
  for(int i=0;i<A->n_atoms;i++){
    if(eq_atoms_A[i]!=-1){
        if(A->nucl_charges_full[i]!=0){
            AB->nucl_charges_full[i]=A->nucl_charges_full[i];
            AB->nucl_charges_calc[i]=A->nucl_charges_calc[i];
            AB->n_ecp_electrons  [i]=A->n_ecp_electrons  [i];
            AB->atom_names   [i]    = new char[3];
            AB->atom_names   [i][0] = A->atom_names   [i][0];
            AB->atom_names   [i][1] = A->atom_names   [i][1];
            AB->atom_names   [i][2] = A->atom_names   [i][2];
            AB->atom_mass    [i]    = A->atom_mass    [i];
            AB->atom_is_ghost[i]    = A->atom_is_ghost[i];
        }
        else if(B->nucl_charges_full[eq_atoms_A[i]]!=0){
            AB->nucl_charges_full[i]=B->nucl_charges_full[eq_atoms_A[i]];
            AB->nucl_charges_calc[i]=B->nucl_charges_calc[eq_atoms_A[i]];
            AB->n_ecp_electrons  [i]=B->n_ecp_electrons  [eq_atoms_A[i]];
            AB->atom_names   [i]    = new char[3];
            AB->atom_names   [i][0] = B->atom_names   [eq_atoms_A[i]][0];
            AB->atom_names   [i][1] = B->atom_names   [eq_atoms_A[i]][1];
            AB->atom_names   [i][2] = B->atom_names   [eq_atoms_A[i]][2];
            AB->atom_mass    [i]    = B->atom_mass    [eq_atoms_A[i]];
            AB->atom_is_ghost[i]    = B->atom_is_ghost[eq_atoms_A[i]];
        }
        else{
            printf("ERROR: both A(%d) and B(%d) are linkers\n",i,eq_atoms_A[i]);
            exit(1);
        }
    }
    else{
        AB->nucl_charges_full[i]  = A->nucl_charges_full[i];
        AB->nucl_charges_calc[i]  = A->nucl_charges_calc[i];
        AB->n_ecp_electrons  [i]  = A->n_ecp_electrons  [i];
        AB->atom_names   [i]    = new char[3];
        AB->atom_names   [i][0] = A->atom_names   [i][0];
        AB->atom_names   [i][1] = A->atom_names   [i][1];
        AB->atom_names   [i][2] = A->atom_names   [i][2];
        AB->atom_mass    [i]    = A->atom_mass    [i];
        AB->atom_is_ghost[i]    = A->atom_is_ghost[i];
        
        
    }
//     printf("WARNING: number of electrons must be tested (mol_link.cpp:311)\n");
    
    AB->n_el_full=A->n_el_full+B->n_el_full;
    AB->n_el_calc=A->n_el_calc+B->n_el_calc;
//     printf("%d %d %d\n",AB->n_el_calc,A->n_el_calc,B->n_el_calc);
//     getchar();

//     printf("atom is %c%c\n",atom_names[i][0],atom_names[i][1]);
//     getchar();
    AB->atom_coord[i*3+0]=A->atom_coord[i*3+0];
    AB->atom_coord[i*3+1]=A->atom_coord[i*3+1];
    AB->atom_coord[i*3+2]=A->atom_coord[i*3+2];
  }
  //copy B
  for(int i=0, i_AB=A->n_atoms;i< B->n_atoms;i++){
    if(eq_atoms_B[i]==-1){
//         AB->n_el_full+=B->nucl_charges_full[i];
//         AB->n_el_calc+=B->nucl_charges_calc[i];
        AB->nucl_charges_full[i_AB]    = B->nucl_charges_full[i];
        AB->nucl_charges_calc[i_AB]    = B->nucl_charges_calc[i];
        AB->n_ecp_electrons  [i_AB]    = B->n_ecp_electrons  [i];
        AB->atom_names   [i_AB]    = new char[3];
        AB->atom_names   [i_AB][0] = B->atom_names   [i][0];
        AB->atom_names   [i_AB][1] = B->atom_names   [i][1];
        AB->atom_names   [i_AB][2] = B->atom_names   [i][2];
        AB->atom_mass    [i_AB]    = B->atom_mass    [i];
        AB->atom_is_ghost[i_AB]    = B->atom_is_ghost[i];
//         printf("atom is %c%c%d\n",AB->atom_names[i][0],AB->atom_names[i][1],i);
//         getchar();
        AB->atom_coord[(i_AB)*3+0]=B->atom_coord[i*3+0];
        AB->atom_coord[(i_AB)*3+1]=B->atom_coord[i*3+1];
        AB->atom_coord[(i_AB)*3+2]=B->atom_coord[i*3+2];
        eq_atoms_B[i]=i_AB;
        i_AB++;
    }
  }
//   for(int i=0; i<A->n_atoms; i++)printf("eq_atoms_A[%d] = %d\n",i,eq_atoms_A[i]);
//   for(int i=0; i<B->n_atoms; i++)printf("eq_atoms_B[%d] = %d\n",i,eq_atoms_B[i]);
  
//   AB->n_el_full-=AB->mol_charge;
//   AB->n_el_calc-=AB->mol_charge;
  AB->z.init(AB->n_atoms);
  AB->z.calc_Z(AB->atom_coord);
  
  return 0;
    
}

int shell_link(std::vector<Shell> * AB, std::vector<Shell> * A, std::vector<Shell> * B){
    
//     printf("WARNING: shell_link was not tested, lib_coef and shell_center were not coppied (see mol_link.cpp:236)\n\n");
    
    int size_A=A->size();
    int size_B=B->size();
    
    AB->resize(size_A+size_B);
    
    for(int i=0; i<size_A;i++)
        AB[0][i]=A[0][i];
    
    for(int i=0; i<size_B;i++)
        AB[0][i+size_A]=B[0][i];
    
    return 0;
}

int pp_link(std::vector<pseudo_potential> * AB, std::vector<pseudo_potential> * A, std::vector<pseudo_potential> * B){
    
    printf("WARNING: pp_link was not tested(see mol_link.cpp:363)\n\n");
    
//     int size_A=A->size();
//     int size_B=B->size();
//     
//     AB->resize(size_A+size_B);
//     
//     for(int i=0; i<size_A;i++)
//         AB[0][i]=A[0][i];
//     
//     for(int i=0; i<size_B;i++)
//         AB[0][i+size_A]=B[0][i];
//     
    return 0;
}


int shell_center_link(std::vector<int> * AB, std::vector<int> * A, std::vector<int> * B, int n_A){
    
//     printf("WARNING: shell_link was not tested, lib_coef and shell_center were not coppied (see mol_link.cpp:236)\n\n");
    
    int size_A=A->size();
    int size_B=B->size();
    
    AB->resize(size_A+size_B);
    
    for(int i=0; i<size_A;i++)
        AB[0][i]=A[0][i];
    
    for(int i=0; i<size_B;i++)
        AB[0][i+size_A]=B[0][i]+n_A;
    
    return 0;
}

int shell_lib_coef_link(std::vector<std::vector<double>> * AB, std::vector<std::vector<double>> * A, std::vector<std::vector<double>> * B){
    
//     printf("WARNING: shell_link was not tested, lib_coef and shell_center were not coppied (see mol_link.cpp:236)\n\n");
    
    int size_A=A->size();
    int size_B=B->size();
    
    AB->resize(size_A+size_B);
    
    for(int i=0; i<size_A;i++)
        AB[0][i]=A[0][i];
    
    for(int i=0; i<size_B;i++)
        AB[0][i+size_A]=B[0][i];
    
    return 0;
}



int mol_link(molecule * AB,molecule * A, molecule * B, int gen_states){
    
//         geom_link
    AB->mol_charge=A->mol_charge+B->mol_charge;
    
    geom_link(AB,A,B);
    
    int n_act_A=0;
    int n_act_B=0;
    
    for(int i=0;i<A->n_frag;i++)n_act_A+=A->n_act_orb[i];
    for(int i=0;i<B->n_frag;i++)n_act_B+=B->n_act_orb[i];
    
    AB->n_wf=A->n_wf*B->n_wf;
            
    //shell_link
    shell_link(&(AB->s   ),&(A->s   ),&(B->s   ));
    shell_link(&(AB->ri_s),&(A->ri_s),&(B->ri_s));
    shell_center_link  (&(AB->shell_center),&(A->shell_center),&(B->shell_center),A->n_atoms);
    shell_center_link  (&(AB->ri_shell_center),&(A->ri_shell_center),&(B->ri_shell_center),A->n_atoms);
    shell_lib_coef_link(&(AB->lib_coef    ),&(A->lib_coef    ),&(B->lib_coef    ));
    shell_lib_coef_link(&(AB->lib_coef_ri ),&(A->lib_coef_ri ),&(B->lib_coef_ri ));
    pp_link(&(AB->PP),&(A->PP),&(B->PP));

//     MO_link
    AB->n_ao=A->n_ao+B->n_ao;
    AB->n_ro=A->n_ao+B->n_ao;
    AB->n_mo=A->n_mo+B->n_mo;
//     printf("AO: %d=%d+%d\n",AB->n_ao,A->n_ao,B->n_ao);
//     printf("MO: %d=%d+%d\n",AB->n_mo,A->n_mo,B->n_mo);
//     getchar();
    
    int n_virt_A = A->n_mo - A->n_cor_orb - n_act_A;
    int n_virt_B = B->n_mo - B->n_cor_orb - n_act_B;
    
    AB->MO_VEC  = new double[(A->n_ao+B->n_ao)*(A->n_ao+B->n_ao)];
    AB->MO_VEC_R= new double[(A->n_ao+B->n_ao)*(A->n_ao+B->n_ao)];
    AB->NO_VEC  = new double[(A->n_ao+B->n_ao)*(A->n_ao+B->n_ao)];
    AB->NO_VEC_B= new double[(A->n_ao+B->n_ao)*(A->n_ao+B->n_ao)];
    MO_link(AB->MO_VEC,A->MO_VEC,B->MO_VEC, A->n_ao, B->n_ao,
                                            A->n_cor_orb, B->n_cor_orb,
                                            n_act_A     , n_act_B,      
                                            n_virt_A     , n_virt_B);
    
    AB->rep_num=new int[A->n_ao+B->n_ao];
    for(int i=0;i<A->n_ao+B->n_ao;i++)AB->rep_num[i]=-1;
    AB->MO_VEC_B=AB->MO_VEC;
    
    AB->orb_energy  = new double[A->n_ao+B->n_ao];
    orb_E_link(AB->orb_energy,
                A->orb_energy,
                B->orb_energy, 
                A->n_ao     , B->n_ao,
                A->n_cor_orb, B->n_cor_orb,
                n_act_A     , n_act_B,      
                n_virt_A     , n_virt_B);
    
    
    //TODO orb_energy, nat_orb_occ
//     printf("AB.MO_VEC:\n");
//     PrintMatr(AB->MO_VEC,AB->n_ao,AB->n_ao,1);
    
//     states_link
    if(gen_states)
    STATES_link(AB,A,B);
    
    
    delete[] eq_atoms_A;
    delete[] eq_atoms_B;
    
    eq_atoms_A=NULL;
    eq_atoms_B=NULL;
    
    return 0;
    
}

int geom_cpy(molecule * O, molecule * IN){

    O->au_coef=IN->au_coef; 
    
    if(O->n_atoms==0){
        O->n_atoms=IN->n_atoms;
        O->atom_coord       =  new double[O->n_atoms*3];
        O->nucl_charges_full = new double[O->n_atoms];
        O->nucl_charges_calc = new double[O->n_atoms];
        O->n_ecp_electrons   = new int   [O->n_atoms];
        O->atom_names        = new char* [O->n_atoms];
        O->atom_mass         = new double[O->n_atoms];
        O->atom_is_ghost     = new int   [O->n_atoms];
        
//         for(int i=0;i<O->n_ion_states;i++)O->basis[i].n_atoms=O->n_atoms;
//         for(int i=0;i<O->n_ion_states;i++)O->basis[i].coords=O->atom_coord;
        
        O->n_el_full = 0;
        O->n_el_calc = 0;
    
        //copy A
        for(int i=0;i<IN->n_atoms;i++){
            O->nucl_charges_full[i]=IN->nucl_charges_full[i];
            O->nucl_charges_calc[i]=IN->nucl_charges_calc[i];
            O->n_ecp_electrons  [i]=IN->n_ecp_electrons  [i];
//             O->n_el_full+=O->nucl_charges_full[i];
//             O->n_el_calc+=O->nucl_charges_calc[i];
            O->atom_names   [i]    = new char[3];
            O->atom_names   [i][0] = IN->atom_names   [i][0];
            O->atom_names   [i][1] = IN->atom_names   [i][1];
            O->atom_names   [i][2] = IN->atom_names   [i][2];
            O->atom_mass    [i]    = IN->atom_mass    [i]   ;
            O->atom_is_ghost[i]    = IN->atom_is_ghost[i]   ;
//     printf("atom is %c%c\n",atom_names[i][0],atom_names[i][1]);
//     getchar();
            O->atom_coord[i*3+0]=IN->atom_coord[i*3+0];
            O->atom_coord[i*3+1]=IN->atom_coord[i*3+1];
            O->atom_coord[i*3+2]=IN->atom_coord[i*3+2];
        }
        O->n_el_full=IN->n_el_full;
        O->n_el_calc=IN->n_el_calc;
        
        O->z = IN->z;
//         for(int i=0;i<IN->n_ion_states;i++)printf("n_el[%d] = %d\n",i+start,O->n_el[i+start]);
//         getchar();
    }
    else{
        int comp_not_passed=0;
        if(O->n_atoms!=IN->n_atoms){
            printf("ERROR: mol_cpy can not copy molecule - number of atoms is different\n");
            exit(1);
        }
        
        //compare geometry
        O->n_el_full = 0;
        O->n_el_calc = 0;
        for(int i=0;i<IN->n_atoms;i++){
            O->n_el_full+=O->nucl_charges_full[i];
            O->n_el_calc+=O->nucl_charges_calc[i];
            if(fabs(O->nucl_charges_full[i]-IN->nucl_charges_full[i])<1e-6 )comp_not_passed=1;
            if(fabs(O->nucl_charges_calc[i]-IN->nucl_charges_calc[i])<1e-6 )comp_not_passed=1;//must be checked
            if(O->n_ecp_electrons  [i]  !=IN->n_ecp_electrons  [i] )comp_not_passed=1;//must be checked
            if(O->atom_coord[i*3+0]!=IN->atom_coord[i*3+0])comp_not_passed=1;
            if(O->atom_coord[i*3+1]!=IN->atom_coord[i*3+1])comp_not_passed=1;
            if(O->atom_coord[i*3+2]!=IN->atom_coord[i*3+2])comp_not_passed=1;
        }
        O->n_el_full-=O->mol_charge;
        O->n_el_calc-=O->mol_charge;
//         for(int i=0;i<IN->n_ion_states;i++)printf("n_el[%d] = %d\n",i+start,O->n_el[i+start]);
//         getchar();
        if(comp_not_passed){
            printf("ERROR: mol_cpy can not copy molecule - geometry is different\n");
            exit(1);
        }
    }
     
    return 0;
    
}

int MO_cpy(double * O, double * IN, int dim_O, int dim_I/*, int * index*/){
    
    set_zero_matr(O,dim_O*dim_O);
    
//     int virtA = dimA-corA-actA;
//     int virtB = dimB-corB-actB;
    
    //corA
    for(int i=0; i<dim_O; i++)
    for(int j=0; j<dim_I; j++)
        O[/*index[*/i/*]*/*dim_O+/*index[*/j/*]*/]+=IN[i*dim_I+j];

    return 0;
    
}

int STATES_cpy(molecule * O, molecule * IN){

    O->n_frag    = IN->n_frag   ;
    O->n_cor_orb = IN->n_cor_orb;
    
    O->n_states     = new int[O->n_frag];
    O->n_act_el_alp = new int[O->n_frag];
    O->n_act_el_bet = new int[O->n_frag];
    O->n_act_orb    = new int[O->n_frag];
    O->n_cor_orb_f  = new int[O->n_frag];
    O->CI           = new aldet_data[O->n_frag];
//     O->csf_occup[0]    = new std::vector<std::vector<std::vector<int>>>  [O->n_frag];
//     O->csf_coef[0]     = new std::vector<std::vector<double>> [O->n_frag];
//     printf("I\n");
//     IN->CI[0].print_states(0,IN->n_states[0]);
    for(int i_f=0; i_f<IN->n_frag;i_f++){
        O->n_act_el_alp[i_f] = IN->n_act_el_alp[i_f];
        O->n_act_el_bet[i_f] = IN->n_act_el_bet[i_f];
        O->n_act_orb[i_f]    = IN->n_act_orb[i_f]   ;
        O->n_cor_orb_f[i_f]  = IN->n_cor_orb_f[i_f] ;
        O->n_states[i_f]     = IN->n_states[i_f]    ;
        if(IN->CI!=NULL){
            aldet_copy(O->CI+i_f,IN->CI+i_f);
        }
        
    }
//     printf("O\n");
//     O->CI[0].print_states(0,O->n_states[0]);
  
    return 0;
    
}


int mol_cpy(molecule *  O, molecule * IN , int gen_states){
        
    O->mol_charge=IN->mol_charge;
    
    geom_cpy(O,IN);
//     int * index;
    O->n_wf=IN->n_wf;
//         basis_link
   O->s=IN->s;
   O->lib_coef    =IN->lib_coef    ;
   O->lib_coef_ri =IN->lib_coef_ri ;
   O->shell_center=IN->shell_center;
   
   O->ri_s=IN->ri_s;
   O->ri_shell_center=IN->ri_shell_center;
   O->PP=IN->PP;
   O->SO_PP=IN->SO_PP;
   
   
//    MO_link
   O->n_mo=IN->n_mo;
   O->n_ao=IN->n_ao;
   O->n_ro=IN->n_ao;
   O->reorder=IN->reorder;
   O-> a_num=IN-> a_num;
   O->cv_num=IN->cv_num;
   O->MO_VEC  = new double[O->n_ao*O->n_ao];
   O->MO_VEC_R= new double[O->n_ao*O->n_ao];
   O->NO_VEC  = new double[O->n_ao*O->n_ao];
   O->NO_VEC_B= new double[O->n_ao*O->n_ao];
   MO_cpy(O->MO_VEC,IN->MO_VEC,O->n_ao,IN->n_ao/*, index*/);
   O->rep_num=new int[O->n_ao];
   for(int i=0;i<O->n_ao;i++)O->rep_num[i]=IN->rep_num[i];
    
   delete[] O->MO_VEC_B;
   O->MO_VEC_B=O->MO_VEC;
   O->orb_energy   = new double[O->n_ao];
   O->nat_orb_occ  = new double[O->n_ao];
   O->nat_orb_occ_B= new double[O->n_ao];
//    for(int i=0;i<O->n_ao;i++)
//    set_zero_matr(O->orb_energy,O->n_ao);
   memcpy(O->orb_energy   , IN->orb_energy   , O->n_ao*sizeof(double));
   if(IN->nat_orb_occ  !=nullptr)memcpy(O->nat_orb_occ  , IN->nat_orb_occ  , O->n_ao*sizeof(double));
   if(IN->nat_orb_occ_B!=nullptr)memcpy(O->nat_orb_occ_B, IN->nat_orb_occ_B, O->n_ao*sizeof(double));

//         states_link
   if(gen_states) STATES_cpy(O, IN);
    
    return 0;
    
}

int mol_array_link(molecule * O, molecule * IN, int n, int gen_states){
    
    molecule * T;
    int m=0;
    if(n>1){
        T = new molecule[n-1];
//         getchar();
        mol_link(T,IN,IN+1,gen_states);
        for(int i=1; i<n-1; i++)mol_link(T+i,T+i-1,IN+i+1,gen_states);
        m=n-2;
    }
    else T = IN;
    
//     T[m].GAMESS_type_out_print("AB.out",0,-1);
    mol_cpy(O,T+m, gen_states);
    if(n>1){//check it
        delete[] T;
    }
    
    return 0;
}
