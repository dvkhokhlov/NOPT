#ifndef GV_SOL_H
#define GV_SOL_H

class gv_solver;

class gv_solver 
{
    public:
        int n_blocks;
        int CI_dim;
        const char **names;
        int * dim;
        int * block_dim;
        double * ev;
        double * H ;
        double * S ;
        double * T ;
        double * T1;
        double * T2;
        double * d_x;
        double * d_y;
        double * d_z;
    
        
        
        int get_block_size_from_files(int ext_n_blocks, const char ** names);
        int get_block_size           (int* ext_dim, const char** ext_names, int ext_n_blocks);
        int alloc();
        int get_HS_data();
        int get_HS_from_global_memory(int i, int j);
        int copy_HS(int i, int j, int k, int l);
        int set_zero_HS(int i, int j);
        int print_init_data(int pH, int pS, int pD, int pE);
        int calc();
        
        
};


#endif
