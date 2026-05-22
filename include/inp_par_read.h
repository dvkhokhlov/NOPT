#ifndef IPR_H
#define IPR_H

# include <vector>

class backup_par // read/write/calc
{
    public:
        int h1;// 0 = read, 1 = calc, 2 = calc and write
//         int h2;
//         int ri;
        
        char * prefix;
        
        backup_par();
        int read_group(char * inp); 
        int disable(); 
//         int read_line(char * inp); 
        ~backup_par();
};

class rhf_par
{
    public:
        int y;
        //guess
//         int huckel_guess;
//         int h_core_guess;
//         int read_guess;
        //convergence
        int max_it;
        double e_conv;
        double g_conv;
        
        rhf_par();
        int read_group(char * inp); 
        int read_line(char * inp); 
        int write_info(int n_cor);
        ~rhf_par();
    
};

class dav_par
{
    public:
        int n_s;
        int n_bf;
        int max_it;
        double e_conv;
        double r_conv;
        double se_min;
        double edshift;
        int sparsed_Hc;

        
        
        dav_par();
        int read_group(char * inp); 
        int read_line(char * inp); 
        int write_info();
        ~dav_par();
    
};


class cas_par
{
    public:
        int y;
        //convergence
        int max_it;
        double e_conv;
        double g_conv;
        double s_conv;
        
        double x_max;
        
        int method;
        int SA;
        int SM;
        int DA;
        
        //states
        int n_s;
        
        std::vector<double> w_state;
        int w_state_type;// 0 - read, 1 - all eq to 1, 2 - separated by symmetry ir.rep.
        std::vector<std::vector<double>> w_state_by_rep;
        std::vector<int> rep_num;
        std::vector<int> rep_lambda;
        std::vector<int> rep_sign;
        std::vector<double> rep_spin;
        
        //tracking
        int track;
        
        //rotate orbitals before cas_par
        int rotate_orbs;
        
        //davidson
        dav_par dav;
        
        cas_par();
        int read_group(char * inp); 
        int read_line(char * inp); 
        int w_linear_read_line(char * line);
        int w_state_by_rep_read(char * inp);
        int write_info(int n_a, int n_b, int n_o, int n_c, int mult);
        ~cas_par();
    
};

class cis_par
{
    public:
        int y;
        //convergence
        
        //states
        int n_s;
	int n_f;
	int method;
	int nat_orb;
	
        std::vector<double> w_state;
        
        //tracking
        int track;
        
        //davidson
        dav_par dav;
        
        cis_par();
        int read_group(char * inp); 
        int read_line(char * inp); 
        int write_info();
        ~cis_par();
    
};

class xmc_par
{
    public:
        int y;
//         int n_s;
        double edshift;
//         std::vector<double> avecoe;
        cas_par * cas;
        int have_ifitd;
        std::vector<double> ifitd_energy;
        int d_only;
        char* gamess_file_name;
        int n_fit;
        int n_fit_pol;
        int pt1_d;
        
        xmc_par();
        int read_group(char * inp, cas_par * ext_cas); 
        int read_line(char * inp); 
        int write_info(int n_a, int n_b, int n_o, int mult);
//         int write_ci_info(int n_a, int n_b, int n_o, int mult);
        ~xmc_par();
        
    
};

class cdas_par
{
    public:
        int y;
//         int n_s;
        double edshift;
//         std::vector<double> avecoe;
        cas_par * cas;
        int have_eps;
        std::vector<double> eps;
        char* gamess_file_name;
//         int mult;
        
        int IPEA;
        int MPPT;
        int HOMO;
        int actual;
        int sing_e;
        int mult_e;
        int orb_e;
        int fit_e;
        int n_orb;
        int rotate_orbs;
        int pt1_d;
        
        cdas_par();
        int read_group(char * inp, cas_par * ext_cas); 
        int read_line(char * inp); 
        int write_info(int n_a, int n_b, int n_o, int mult);
//         int write_ci_info(int n_a, int n_b, int n_o, int mult);
        ~cdas_par();
        
    
};

# include "molecule.h"

class inp_par
{
    public:
        char* inp_name;
        char* job_name;
        std::vector<int> mol_line;
        std::vector<int> act_space_line;
        
        char* out_folder_name;
        int doci_dec;
        
        char* point_group;
        
        rhf_par rhf;
        cas_par cas;
        cis_par cis;
        xmc_par xmc;
        cdas_par cdas;
        
//         int charge;
//         int parsing_done;
        
//         int non_diag=0;
        int ecp;
        int num_threads;
        int out_folder_defined;
        
        inp_par();
        int write_info();
        int read(char * inp_name);
        int read_p_group(char * inp);
        int read_p_line(char * inp); 
        std::vector<int> count_groups(char * inp, std::vector<const char *>start_kw, std::vector<const char *>stop_kw);
        int read_m_line(char * inp, int n); 
        int read_symm_group(char * inp);
        int read_symm_line(char * inp); 
        ~inp_par();
        
};

std::vector<int> molecule_read_by_inp_par(molecule * M, inp_par * P);



#endif
