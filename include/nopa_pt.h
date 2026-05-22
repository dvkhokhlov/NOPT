#ifndef NOPA_PT_H
#define NOPA_PT_H

# include <vector>


class nopa_v_data;

class nopa_v_data
{
public:
    
    int n_block;
    
    double * V;
    
    int CI_dim;
    int * dim;
    int * block_dim;
    double * H;
    double * S;
    
    
    std::vector<double>*E_PT;
    double ** V_A;
    double ** S_A;
    double ** V_B;
    double ** S_B;
    double ** V_A_0;
    double ** V_B_0;
    double * T2;//bufer
    char * v_name;
    
    int * n_a;
    int * n_b;
    
    double * S_ort;
    
    
    int H_update_v_ort(double * H, double * S, int argc, char ** argv, int CI_dim, int * dim, int * block_dim);
    int H_update_simple(double * H, double * S, int argc, char ** argv, int CI_dim, int * dim, int * block_dim);
    int H_update_simple_with_non_diag_V(double * H, double * S, int argc, char ** argv, int CI_dim, int * dim, int * block_dim);
    int H_update_leow(double * H, double * S, int argc, char ** argv, int CI_dim, int * dim, int * block_dim);
    int H_update_leow_with_non_diag_V(double * H, double * S, int argc, char ** argv, int CI_dim, int * dim, int * block_dim);
    int H_update_leow_G(double * ext_H, double * ext_S, int ext_CI_dim, int * ext_dim, int * ext_block_dim);
    int alloc(int n);
    int alloc2(int n, int ext_CI_dim);
    int get_V_from_global_memory(int i, int a);
    int V_read(int i, char * name);
    int copy_V(int i, int j);
    int V_read_all(char ** name);
    int V_block_recalc(int i);
    int V_block_recalc_all();
    int V_orthogonal_calc();
    
    int read_1el(int a, int b, char ** name);
    int read_1el_all(char ** name);
    
    int H_update();
    int mulliken_update(double * M);
    int leowdin_update (double * M, int inv);
        
};

#endif
