# ifdef LIBINT
# include <vector>
# include <libint2/atom.h>
using libint2::Atom;
# include <libint2/shell.h>
using libint2::Shell;
#endif
#include "pseudo_potential.h"
#include "aldet.h"
#ifndef __molecule
#define __molecule
#include "symmetry.h"
#include "z_matrix.h"

class molecule;

class molecule 
{
    public:
        int n_atoms;
        int n_frag;
        
        double *atom_coord;
//         double *mcs_rot_matr;
        double *nucl_charges_full;//real charge for input/output
        double *nucl_charges_calc;//effective charge for calculation 
        int * n_ecp_electrons  ;//number of electrons frozen by ECP
        
        double *atom_mass;
        char ** atom_names;
        int * atom_is_ghost;
        Z_matrix z;
        
        double au_coef;
        
        //for each ionic state
//         AO_basis * basis; //to be deleted!!! cannot be deleted because is needed by:
                          //molecule_1el_data::RxP_calc
                          //molecule_1el_data::P_calc
                          //smth else
        
        std::vector<Shell> s;
        std::vector<std::vector<double>> lib_coef;
        std::vector<int> shell_center;
        
        std::vector<Shell> ri_s;
        std::vector<std::vector<double>> lib_coef_ri;
        std::vector<int> ri_shell_center;
        
        std::vector<Shell> read_s;
        
        std::vector<pseudo_potential> PP;
        
        std::vector<pseudo_potential> SO_PP;
        
        
        
        int n_ao; // dimension of the AO basis
        int n_ro; // dimension of the AO basis for reading
        int n_mo; // number of MO can be less then n_ao
        double * MO_VEC;
        double * MO_VEC_B;
        double * MO_BUF;
        double * MO_VEC_R;//vectors for reading with additional basis
        double * orb_energy;
        double * orb_energy_B;
        int mol_charge;
        int n_wf;
        int n_el_full;//full number
        int n_el_calc;//effective number
        int * n_act_el_alp;
        int * n_act_el_bet;
        int * n_act_orb;
        int * n_states;
        int n_cor_orb;
        int * n_cor_orb_f;
        
        int reorder;
        std::vector<int> a_num; //active orbitals numbers for reordering
        std::vector<int> cv_num;//core and virtual orbitals numbers for reordering
        
        //natural orbitals
        double * nat_orbs_m;
        double * NO_VEC;
        double * NO_VEC_B;
        double * nat_orb_occ;
        double * nat_orb_occ_B;
        
        //symmetry
        symmetry S;
        int * rep_num;
        int * rep_num_backup;
        int * rep_AO_num;
        //linear
        double * Lambda_AO;
        
        //CI
        int * read_CI;
        char ** CI_file_name;
        int n_CI;
        aldet_data * CI;
        
        char * source_name;
        int inp_mol_line;
        int inp_act_space_line;
        
        char * basis_name;
        int basis_in_lib;
        int basis_format;
        int basis_start_line;
        
        char * ri_basis_name;
        int ri_basis_in_lib;
        int ri_basis_format;
        int ri_basis_start_line;
        
        char * read_basis_name;
        int read_basis_in_lib;
        int read_basis_format;
        int read_basis_start_line;
        
        char * ecp_name;
        int ecp_in_lib;
        int ecp_format;
        int ecp_start_line;
        
        char * so_name;
        int so_in_lib;
        int so_format;
        int so_start_line;
        
        
        char * MO_name;
        
        double *states_energies;
        
        int ecp;
        
        int n_ao_d5;    //number of molecular orbitals excluding d6-artifact orbitals
        int have_ecp;   //flag about presence of effective core potential
        double * S_AO ; //orbital overlap matrix can be in atomic basis
        double * S_MO ; //orbital overlap matrix can be in molecular basis (check the code to put correctly)
        double * S_M05;
        double * H_AO ; //1-el Hamiltonian in atomic basis
        double * Dx_AO; //1-el dipole_x in atomic basis
        double * Dy_AO; //1-el dipole_y in atomic basis
        double * Dz_AO; //1-el dipole_z in atomic basis
        double * Qxx_AO;//1-el quadrupole_xx in atomic basis
        double * Qxy_AO;//1-el quadrupole_xy in atomic basis
        double * Qxz_AO;//1-el quadrupole_xz in atomic basis
        double * Qyy_AO;//1-el quadrupole_yy in atomic basis
        double * Qyz_AO;//1-el quadrupole_yz in atomic basis
        double * Qzz_AO;//1-el quadrupole_zz in atomic basis
        double * P_AO;  //1-el impulse (3 projections) in atomic basis
        double * RxP_AO;//1-el coordinate-impulse vector product (3 projections) in atomic basis
        double * BUF;
        
        double * V_D6;  //                              d6-artifact orbitals
        double * T_D5;  //rotation matrix for excluding d6-artifact orbitals
        
//         std::vector<Atom> atoms;//atomic coordinates in the libint format
        double V_nuc;               //V_nn - nuclear interaction energy
        double Dx_nuc,Dy_nuc,Dz_nuc;//dipole components caused by nuclear charge
        char * prefix;  //prefix for names of input/output files
        int h1_rwc;
        
        int mc;//flag for RI calculation of F operator -- calc_F_AO()
        //mc=0 - single reference
        //mc=1 - multireference SCF calculations
        //mc=2 - multireference non-SCF calculations
        
        double F_x, F_y, F_z; //// electric field intensity
        
        
        molecule();
        
        int geom_read();
        int geom_print();
        int basis_read();
        int basis_name_read();
        int basis_print(int t);
        int PP_print();
        
        int active_space_read(int warning, int read_states, int ci_alloc = ALDET_ALLOC_FULL);
        int reorder_orbitals();
        int STATES_set_zero();
        int active_space_print(int f);
        int MO_gen();
        int huckel_guess();
        int project_to_full_basis();
        int UHF_MO_read();
        int molden_read();
//         int MO_read_wn(int i);
        int SCF_type_def();
        int MO_print(const char * out_name);
//         int MO_cpy(int i,int j, const char* type);
        int write_orb_index();
        int NO_print(const char * name, char S);
        // dump the active orbitals rotated by U (n_act x n_act, [a*n_act+p]) as a GAMESS
        // visualization, e.g. localized orbitals. MO_VEC is restored on return (transient).
        int LOC_print(const char * name, const double * U);
        int GAMESS_type_out_print(const char* out_name, int n_print_mo);
        int GAMESS_geom_print(const char* out_name);
        double center_of_mass(int xyz);
        int inertion_tensor(double * T);
//         int move_to_mcs(int calc_rm, int inv);//molecular coordinate system;
        # ifdef LIBINT
        int add_libint_atoms(std::vector<Atom> *);
        int MO_libint_format();//from MO_VEC_R to MO_VEC_R
        int MO_gamess_format();//from MO_VEC   to MO_VEC_R
        std::vector<std::pair<double,std::array<double, 3>>> libint_point_charges();
        #endif
        
        int efrag_read();
        
        int nat_orb_calc(double * gamma, int renorm, char S);
        
        int make_symmetry();
        
        int lin_prepare();
        
        int check_orb_symmetry();
        
        int sort_orbs(char t);
        
        int sort_orbs_by_rep(int first, int last);
        
        int MO_backup();
        
        int MO_restore();
        
        int sync_shell_to_geom();
        
        int resize_CI(int n);
        
        int alloc_basic(); //allocation of basic memory
        int alloc_Q_AO();     //allocation of special memory fo quadrupole operators
        int clear_Q_AO();     //cleaning of special memory fo quadrupole operators
        
        int gen_1el_data();//set parameters --- old function from NOPT-2.x.x
        
        int calc_S_AO();
        int calc_SM05();
        int MO_orth();
        int calc_H_AO();
        
        int write_H_AO();
        int read_H_AO();
        int calc_F_AO  (double * F               , double * DM                 , double c);
        int calc_U_F_AO(double * F_A,double * F_B, double * DM_A, double * DM_B, double c);
        int calc_D_AO();
        int calc_Q_AO();
        int RxP_calc();//written with old AO_basis sturcture
        int P_calc();  //written with old AO_basis sturcture
        
        int diag_X_AO_in_MO(double * X);
        int diag_X_MO_block(double * X, int n0, int dim, double * U);
        
        int calc_d5_n_ao();
        int make_D5_matr();
        
        
        ~molecule();
        //functions not transfered from doCI
//         int geom_shift(double dx, double dy, double dz);
//         int geom_rotate(double fi,
//                         double  a, double  b, double  c,
//                         double xc, double yc, double zc);

    
};


int VEC_energy_cpy(double * O_V,double * I_V, double * o_e, double * i_e, int n);

int GAMESS_VEC_data_print(FILE * out, double * V, int n_mo, int n_ao);

int mol2_fprint(const char* out_name, molecule * a, molecule * b);

double E_nuc_calc(molecule *a);

double V_nuc_calc(molecule *a, molecule *b);

double D_nuc_calc(molecule * a, int xyz);


#endif

int transform_from_col_MO(double * M_O, double * M, int n, double * V1, int n1, double * V2, int n2);
int transform_from_col_MO_compl(double * M_O, double * M, int n, double * V1, int n1, double * V2, int n2);

int change_orbs(double * V, int i, int j, int dim);
