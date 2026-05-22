#ifndef __pseudo_potential
#define __pseudo_potential

#include <vector>
# include <libint2/shell.h>
using libint2::Shell;


#ifdef _USE_GRPP
extern "C"{
    # include <libgrpp.h>
}
#endif

class pseudo_potential
{
    public:
        int n_pot;
        int n_el;// not yet taken into account
        std::vector<std::vector<double>> coeff;
        std::vector<std::vector<double>> alpha;
        std::vector<std::vector<int   >> n;
        
        int n_oc;
        std::vector<Shell> out_core;
        std::vector<int> oc_2J;//2*J of the outcore orbital
        std::vector<std::vector<double>> oc_lib_shel_coeff;
        std::vector<std::vector<double>> oc_coeff;
        std::vector<std::vector<double>> oc_alpha;
        std::vector<std::vector<int   >> oc_n;
        
        std::vector<int> l;
        //l[0]=-1 is general term
        //l[1]=0 is S correction, l[2]=1 - P
        int atom;
        
        pseudo_potential();
        int set_n_pot();
        int set_pot_n_terms(int i, int n);
        int print(FILE * stream, char * name, int type);
        int print_pot(int i, FILE * stream, int type);
        int print_oc (int i, FILE * stream, int type);
#ifdef _USE_GRPP        
        libgrpp_potential_t *make_grpp(int i);
        libgrpp_potential_t *make_oc_grpp(int i);
#endif

}; 

#endif
