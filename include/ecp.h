#ifndef __ecp_H
#define __ecp_H
#ifdef _USE_GRPP
# include <libint2/shell.h>
extern "C"{
    # include <libgrpp.h>
}
# include <molecule.h>

using libint2::Shell;

class nopt_grpp;

class nopt_grpp
{
    public:
        //grpp-formated shells
        int n_sh;
        int n_ao;
        libgrpp_shell_t ** s;
        
        //potentials
        int n_pot;
        libgrpp_potential_t ** p;
        double * origin;
        
        //SO potentials
        int n_so;
        libgrpp_potential_t ** so_p;
        double * so_origin;
        
        //out_core potentials ECP
        int n_oc;
        libgrpp_shell_t ** oc_s;
        libgrpp_potential_t ** oc_p;
        double * oc_origin;
        
        //out_core potentials SO
        int n_oc_so;
        libgrpp_shell_t ** oc_s_so;
        libgrpp_potential_t ** oc_p_so;
        double * oc_origin_so;
        
        
        
        int init_done;
        
        nopt_grpp();
        int set_mol(molecule * M);
        int refresh(molecule * M);
        int calc_H_AO(double * m);
        int calc_SO_AO(double * Sx, double * Sy, double * Sz);
        ~nopt_grpp();
   
};

extern nopt_grpp grpp_engine;
#endif
#endif
