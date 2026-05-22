#ifndef __symmetry
#define __symmetry


class symmetry;

class symmetry 
{
    public:
        char * group;
        char T;
        int n;
        char t;
        //symmetry operations
        int n_op;
        double ** Op_Space;
        double ** Op_Orb;
        //irreducible representations
        int n_rep;
        double ** ch;
        char ** rep_name;
        double * dim;
        double ** P;
        //for linear
        double ** T_D5;
        int * n_ao_d5;
        
        //molecule
        int * at_refl;
        
        
        symmetry();
        int alloc();
        int init(char* ext_name);
        int print_rep(FILE *, int n);
        ~symmetry();
        
        
};

int OxR(double * OR, double * O, double *R);

int at_search(double * R, double * L, int n, double eps);

#endif
