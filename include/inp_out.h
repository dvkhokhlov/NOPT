#ifndef INP_OUT_H
#define INP_OUT_H


#include <vector>

class recursive_file;

class recursive_file{
    
    public:
        std::vector<FILE *> file;
        int level;
        
        recursive_file();
        int r_open(char * name);
        int go_to_begining();
        int r_gets(char * out, int l);
        int r_smart_fgets(char * l, int * n, char* fn, const char * job);
        int r_eof();
        int r_close();
        
        ~recursive_file();
};

int smart_fgets(char * l, int * n, FILE* f, char* fn, const char * job);

int print_file_content(char * inp);

int read_block_size(const char * B);

int read_HS_data(const char * B, const char * K, double * H, double * S, double * x, double * y, double * z, int ny, int nx, int ld);

int read_H_data(const char * B, const double * H, int nx);

int read_gv_eigenvalue(double *C, char * ev_name, int * dim, int w, int s);

int ReadMatr10(double * M, int n, int m, char * f_name, const char * keyword);

int PrintDipoleTable(FILE * stream, int dim);

int PrintEnergyTable(FILE * stream, int dim);

int ReadDipole(double * x, double * y, double * z, int n, char * f_name, const char * keyword);

int PrintDipole(double * x, double * y, double * z, int n);

int PrintDispersion(double * xx, double * yy, double * zz, int n);

int PrintQuadrupole(double * xx, double * yy, double * zz,
                    double * xy, double * xz, double * yz, int n);

int PrintEnergy(double * ev, int dim, int scan);

int PrintEnergy_derivative_corrected(double * E, double * D, int dim, int change);

int PrintDer(double * D, int n_s, int n_d);

int PrintRowVec(double * V, int n1, int n2);

std::vector<double> read_E_XMC(const char* B, const char * E);

int print_CI_results(FILE * target, int i, int j, double S, double V_nuc, double H_1el, double H_2el, double d_x, double d_y, double d_z);
#endif
