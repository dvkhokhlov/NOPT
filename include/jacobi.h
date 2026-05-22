class jacobi_mcscf_sd_engine;

class jacobi_mcscf_sd_engine{
    public:
        
        int n_c;
        int n_a; 
        int n_v; 
        int n_s; 
        int n_ao;
        int n_mo;
        int dim;
        
        double x_max_jsd;
        
        double * G;
        double * B;
        double * rot_matr;
        
        jacobi_mcscf_sd_engine();
        int init(double * ext_G, double * ext_B, int ext_n_c, int ext_n_a, int ext_n_v, int ext_n_ao, int ext_n_s, double ext_x_max_jsd);
        double find_max_el();
        double step(double * VEC, double * BUF);
        double rotation(int ab, int a, int b);
        ~jacobi_mcscf_sd_engine();
        
};
