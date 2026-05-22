#ifdef _USE_GRPP
#include "ecp.h"
#include "matr.h"
# include "common_vars.h"
# include <omp.h>

nopt_grpp grpp_engine;

#define MAX_BUF 10000

void add_block_to_matrix(int dim_1, int dim_2, double *matrix,
                         int block_dim_1, int block_dim_2, double *block,
                         int col_offset, int row_offset, double factor)
{
    for (int i = 0; i < block_dim_1; i++) {
        for (int j = 0; j < block_dim_2; j++) {
            double element = block[i * block_dim_2 + j];
            matrix[(col_offset + i) * dim_2 + (row_offset + j)] += factor * element;
        }
    }
}


int norm(int num_shells, libgrpp_shell_t **shell_list){
    
    double buf[MAX_BUF];

    for (int ishell = 0; ishell < num_shells; ishell++) {

        libgrpp_shell_t *shell = shell_list[ishell];
        libgrpp_overlap_integrals(shell, shell, buf);
        double norm = 1.0/sqrt(buf[0]);
        
        for(int i_p=0;i_p<shell->num_primitives;i_p++)
            shell->coeffs[i_p]=shell->coeffs[i_p]*norm;
    }
    
    return 0;
}

int print_norm(int num_shells, libgrpp_shell_t **shell_list){
    
    double buf[MAX_BUF];

    for (int ishell = 0; ishell < num_shells; ishell++) {

        libgrpp_shell_t *shell = shell_list[ishell];
        libgrpp_overlap_integrals(shell, shell, buf);
        printf("n(%d) = %e\n",ishell, buf[0]);
        getchar();
    }
    
    return 0;
}




nopt_grpp::nopt_grpp(){
    
    s           =nullptr;
    p           =nullptr;
    origin      =nullptr;
    so_p        =nullptr;
    so_origin   =nullptr;
    oc_s        =nullptr;
    oc_p        =nullptr;
    oc_origin   =nullptr;
    oc_s_so     =nullptr;
    oc_p_so     =nullptr;
    oc_origin_so=nullptr;

    init_done=0;
    
}



int nopt_grpp::set_mol(molecule * M){
    
    
    if(init_done){
        refresh(M);
        return 0;
    }
    init_done=1;
    
    //allocate shells
    n_sh = M->s.size();
    
    s = (libgrpp_shell_t **) malloc(sizeof(libgrpp_shell_t *) * n_sh);
    
    //allocate ECP(1)
    n_pot=0;
    
    for(int i=0;i<M->PP.size();i++) n_pot+=M->PP[i].n_pot;
    
    p = (libgrpp_potential_t **) malloc(sizeof(libgrpp_potential_t *) * n_pot);
    origin = new double[3*n_pot];

    //allocate SO
    n_so=0;
    
    for(int i=0;i<M->SO_PP.size();i++) n_so+=M->SO_PP[i].n_pot;
    
    so_p = (libgrpp_potential_t **) malloc(sizeof(libgrpp_potential_t *) * n_so);
    so_origin = new double[3*n_so];
    
    //allocate out_core
    n_oc=0;
    for(int i=0;i<M->PP.size();i++) n_oc+=M->PP[i].n_oc;
    oc_s = (libgrpp_shell_t **) malloc(sizeof(libgrpp_shell_t *) * n_oc);
    oc_p = (libgrpp_potential_t **) malloc(sizeof(libgrpp_potential_t *) * n_oc);
    oc_origin = new double[3*n_oc];

    //allocate out_core SO
    n_oc_so=0;
    for(int i=0;i<M->SO_PP.size();i++) n_oc_so+=M->SO_PP[i].n_oc;
    oc_s_so = (libgrpp_shell_t **) malloc(sizeof(libgrpp_shell_t *) * n_oc_so);
    oc_p_so = (libgrpp_potential_t **) malloc(sizeof(libgrpp_potential_t *) * n_oc_so);
    oc_origin_so = new double[3*n_oc_so];

    //generate shells
    n_ao=M->n_ao;
    for(int i=0;i<n_sh;i++){
        s[i]=libgrpp_new_shell(M->s[i].O.data(), 
                               M->s[i].contr[0].l, 
                               M->s[i].alpha.size(), 
                               M->lib_coef[i].data(), 
                               M->s[i].alpha.data());        
    }
    norm(n_sh, s);
    
    
    //generate ECP(1)
    for(int i=0,i_p=0;i<M->PP.size()  ;i++)
    for(int j=0           ;j<M->PP[i].n_pot;j++){
        p[i_p]=M->PP[i].make_grpp(j);
        origin[3*i_p  ]=M->atom_coord[3*M->PP[i].atom  ];
        origin[3*i_p+1]=M->atom_coord[3*M->PP[i].atom+1];
        origin[3*i_p+2]=M->atom_coord[3*M->PP[i].atom+2];
        i_p++;
    }
    //generate SO
    for(int i=0,i_p=0;i<M->SO_PP.size()  ;i++)
    for(int j=0           ;j<M->SO_PP[i].n_pot;j++){
        so_p[i_p]=M->SO_PP[i].make_grpp(j);
        so_origin[3*i_p  ]=M->atom_coord[3*M->SO_PP[i].atom  ];
        so_origin[3*i_p+1]=M->atom_coord[3*M->SO_PP[i].atom+1];
        so_origin[3*i_p+2]=M->atom_coord[3*M->SO_PP[i].atom+2];
        i_p++;
    }
    //generate out_core
    for(int i=0,i_sh=0;i<M->PP.size() ;i++)
    for(int j=0       ;j<M->PP[i].n_oc;j++){
        oc_s[i_sh]=libgrpp_new_shell(M->PP[i].out_core[j].O.data(), 
                                     M->PP[i].out_core[j].contr[0].l, 
                                     M->PP[i].out_core[j].alpha.size(), 
                                     M->PP[i].oc_lib_shel_coeff[j].data(), 
                                     M->PP[i].out_core[j].alpha.data());
        
        oc_p[i_sh]=M->PP[i].make_oc_grpp(j);
        oc_origin[3*i_sh  ]=M->atom_coord[3*M->PP[i].atom  ];
        oc_origin[3*i_sh+1]=M->atom_coord[3*M->PP[i].atom+1];
        oc_origin[3*i_sh+2]=M->atom_coord[3*M->PP[i].atom+2];
        i_sh++;
    }
    //generate out_core SO
    for(int i=0,i_sh=0;i<M->SO_PP.size() ;i++)
    for(int j=0       ;j<M->SO_PP[i].n_oc;j++){
        oc_s_so[i_sh]=libgrpp_new_shell(M->SO_PP[i].out_core[j].O.data(), 
                                        M->SO_PP[i].out_core[j].contr[0].l, 
                                        M->SO_PP[i].out_core[j].alpha.size(), 
                                        M->SO_PP[i].oc_lib_shel_coeff[j].data(), 
                                        M->SO_PP[i].out_core[j].alpha.data());
        
        oc_p_so[i_sh]=M->SO_PP[i].make_oc_grpp(j);
        oc_origin_so[3*i_sh  ]=M->atom_coord[3*M->SO_PP[i].atom  ];
        oc_origin_so[3*i_sh+1]=M->atom_coord[3*M->SO_PP[i].atom+1];
        oc_origin_so[3*i_sh+2]=M->atom_coord[3*M->SO_PP[i].atom+2];
        i_sh++;
    }
    
//     print_norm(n_oc, oc_s);
    
//     norm(n_oc, oc_s); must not be normalized???????
    
    return 0;
}
    
int nopt_grpp::refresh(molecule * M){
    
    // printf("nopt_grpp::refresh() is not written\n");
    // exit(1);
    
    
    for(int i=0;i<n_sh;i++){
        s[i]->origin[0]=M->s[i].O[0];
        s[i]->origin[1]=M->s[i].O[1];
        s[i]->origin[2]=M->s[i].O[2];
    }
    norm(n_sh, s);
    
    
    //generate ECP(1)
    for(int i=0,i_p=0;i<M->PP.size()  ;i++)
    for(int j=0           ;j<M->PP[i].n_pot;j++){
        // p[i_p]=M->PP[i].make_grpp(j);
        origin[3*i_p  ]=M->atom_coord[3*M->PP[i].atom  ];
        origin[3*i_p+1]=M->atom_coord[3*M->PP[i].atom+1];
        origin[3*i_p+2]=M->atom_coord[3*M->PP[i].atom+2];
        i_p++;
    }
    //generate SO
    for(int i=0,i_p=0;i<M->SO_PP.size()  ;i++)
    for(int j=0           ;j<M->SO_PP[i].n_pot;j++){
        // so_p[i_p]=M->SO_PP[i].make_grpp(j);
        so_origin[3*i_p  ]=M->atom_coord[3*M->SO_PP[i].atom  ];
        so_origin[3*i_p+1]=M->atom_coord[3*M->SO_PP[i].atom+1];
        so_origin[3*i_p+2]=M->atom_coord[3*M->SO_PP[i].atom+2];
        i_p++;
    }
    //generate out_core
    for(int i=0,i_sh=0;i<M->PP.size() ;i++)
    for(int j=0       ;j<M->PP[i].n_oc;j++){
        oc_s[i_sh]->origin[0]=M->PP[i].out_core[j].O[0];
        oc_s[i_sh]->origin[1]=M->PP[i].out_core[j].O[1];
        oc_s[i_sh]->origin[2]=M->PP[i].out_core[j].O[2];
        
        // oc_p[i_sh]=M->PP[i].make_oc_grpp(j);
        oc_origin[3*i_sh  ]=M->atom_coord[3*M->PP[i].atom  ];
        oc_origin[3*i_sh+1]=M->atom_coord[3*M->PP[i].atom+1];
        oc_origin[3*i_sh+2]=M->atom_coord[3*M->PP[i].atom+2];
        i_sh++;
    }
    //generate out_core SO
    for(int i=0,i_sh=0;i<M->SO_PP.size() ;i++)
    for(int j=0       ;j<M->SO_PP[i].n_oc;j++){
        oc_s_so[i_sh]->origin[0]=M->SO_PP[i].out_core[j].O[0];
        oc_s_so[i_sh]->origin[1]=M->SO_PP[i].out_core[j].O[1];
        oc_s_so[i_sh]->origin[2]=M->SO_PP[i].out_core[j].O[2];
        
        // oc_p_so[i_sh]=M->SO_PP[i].make_oc_grpp(j);
        oc_origin_so[3*i_sh  ]=M->atom_coord[3*M->SO_PP[i].atom  ];
        oc_origin_so[3*i_sh+1]=M->atom_coord[3*M->SO_PP[i].atom+1];
        oc_origin_so[3*i_sh+2]=M->atom_coord[3*M->SO_PP[i].atom+2];
        i_sh++;
    }
    
    
    return 0;
}

int nopt_grpp::calc_H_AO(double * m){
    
    double * MP;
    MP = new double[ num_threads*(n_ao*n_ao+64)];
    set_zero_matr(MP,num_threads*(n_ao*n_ao+64));
    
    
    
#pragma omp parallel
//     for(int thread_id=0;thread_id<num_threads;thread_id++)
    {
        double buf[MAX_BUF];
        double buf_x[MAX_BUF];
        double buf_y[MAX_BUF];
        double buf_z[MAX_BUF];
        
        
        double * __restrict__ MT;//trying restrict??????
        int thread_id = omp_get_thread_num();
        MT = MP + thread_id*(n_ao*n_ao+64);
        

        int ioffset = 0;
        for (int i_s=0;i_s<n_sh;i_s++) {
            libgrpp_shell_t *bra = s[i_s];
            int bra_dim = libgrpp_get_shell_size(bra);
        
            int joffset = 0;
            for (int j_s=0;j_s<n_sh;j_s++){
                
                libgrpp_shell_t *ket = s[j_s];
                int ket_dim = libgrpp_get_shell_size(ket);
                if (((i_s*n_sh+j_s) % num_threads) == thread_id){
                    
                    
                    for(int i_p=0;i_p<n_pot;i_p++){

                        if(p[i_p][0].L==-1)
                            libgrpp_type1_integrals (bra, ket, origin+3*i_p, p[i_p], buf);
                        else
                            libgrpp_type2_integrals (bra, ket, origin+3*i_p, p[i_p], buf);
                        add_block_to_matrix(n_ao, n_ao, MT, bra_dim, ket_dim, buf, ioffset, joffset, 1.0);
                    }
                    
                    libgrpp_outercore_potential_integrals(bra, ket, oc_origin, n_oc, oc_p, oc_s, buf, buf_x, buf_y, buf_z);
                    add_block_to_matrix(n_ao, n_ao, MT, bra_dim, ket_dim, buf, ioffset, joffset, 1.0);
                        
                }
                
                joffset += ket_dim;
            }
        
            ioffset += bra_dim;
        }
//         printf("H:\n");
//         PrintMatr(m,n_ao,n_ao,1);
    }
    
    for(int i=0; i<n_ao*n_ao; i++){
        for(int j=0; j<num_threads; j++)
            m[i]+=MP[j*(n_ao*n_ao+64)+i];
        
    }
    
    delete [] MP;
    
    return 0;
}


int nopt_grpp::calc_SO_AO(double * Sx, double * Sy, double * Sz){
    
    double * SxP;
    SxP = new double[ num_threads*(n_ao*n_ao+64)];
    set_zero_matr(SxP,num_threads*(n_ao*n_ao+64));
    
    double * SyP;
    SyP = new double[ num_threads*(n_ao*n_ao+64)];
    set_zero_matr(SyP,num_threads*(n_ao*n_ao+64));
    
    double * SzP;
    SzP = new double[ num_threads*(n_ao*n_ao+64)];
    set_zero_matr(SzP,num_threads*(n_ao*n_ao+64));
    
    
    
#pragma omp parallel
//     for(int thread_id=0;thread_id<num_threads;thread_id++)
    {
        double buf_x[MAX_BUF];
        double buf_y[MAX_BUF];
        double buf_z[MAX_BUF];
        double buf_n[MAX_BUF];
        
        double * __restrict__ SxT;//trying restrict??????
        double * __restrict__ SyT;//trying restrict??????
        double * __restrict__ SzT;//trying restrict??????
        int thread_id = omp_get_thread_num();
        SxT = SxP + thread_id*(n_ao*n_ao+64);
        SyT = SyP + thread_id*(n_ao*n_ao+64);
        SzT = SzP + thread_id*(n_ao*n_ao+64);
        

        int ioffset = 0;
        for (int i_s=0;i_s<n_sh;i_s++) {
            
            libgrpp_shell_t *bra = s[i_s];
            int bra_dim = libgrpp_get_shell_size(bra);
        
            int joffset = 0;
            for (int j_s=0;j_s<n_sh;j_s++){
                
                libgrpp_shell_t *ket = s[j_s];
                int ket_dim = libgrpp_get_shell_size(ket);
                if (((i_s*n_sh+j_s) % num_threads) == thread_id){
                    
                    
                    for(int i_p=0;i_p<n_so;i_p++){
                        libgrpp_spin_orbit_integrals(bra, ket, so_origin+3*i_p, so_p[i_p], buf_x, buf_y, buf_z);
                        
                        add_block_to_matrix(n_ao, n_ao, SxT, bra_dim, ket_dim, buf_x, ioffset, joffset, 1.0);
                        add_block_to_matrix(n_ao, n_ao, SyT, bra_dim, ket_dim, buf_y, ioffset, joffset, 1.0);
                        add_block_to_matrix(n_ao, n_ao, SzT, bra_dim, ket_dim, buf_z, ioffset, joffset, 1.0);
                    }
                    libgrpp_outercore_potential_integrals(bra, ket, oc_origin_so, n_oc_so, oc_p_so, oc_s_so, buf_n, buf_x, buf_y, buf_z);
                    
                    add_block_to_matrix(n_ao, n_ao, SxT, bra_dim, ket_dim, buf_x, ioffset, joffset, 1.0);
                    add_block_to_matrix(n_ao, n_ao, SyT, bra_dim, ket_dim, buf_y, ioffset, joffset, 1.0);
                    add_block_to_matrix(n_ao, n_ao, SzT, bra_dim, ket_dim, buf_z, ioffset, joffset, 1.0);
//                     }
                    
                }
                
                joffset += ket_dim;
            }
        
            ioffset += bra_dim;
        }
//         printf("H:\n");
//         PrintMatr(m,n_ao,n_ao,1);
    }
    
    for(int i=0; i<n_ao*n_ao; i++){
        for(int j=0; j<num_threads; j++)
            Sx[i]+=SxP[j*(n_ao*n_ao+64)+i];
        
    }
    for(int i=0; i<n_ao*n_ao; i++){
        for(int j=0; j<num_threads; j++)
            Sy[i]+=SyP[j*(n_ao*n_ao+64)+i];
        
    }
    for(int i=0; i<n_ao*n_ao; i++){
        for(int j=0; j<num_threads; j++)
            Sz[i]+=SzP[j*(n_ao*n_ao+64)+i];
        
    }
    
    delete [] SxP;
    delete [] SyP;
    delete [] SzP;
    
    return 0;
}



nopt_grpp::~nopt_grpp(){
    
    if(s!=nullptr){
        for(int i=0;i<n_sh;i++)libgrpp_delete_shell(s[i]);
        free(s);
    }
    if(oc_s!=nullptr){
        for(int i=0;i<n_oc;i++)libgrpp_delete_shell(oc_s[i]);
        free(oc_s);
    }
    if(oc_s_so!=nullptr){
        for(int i=0;i<n_oc_so;i++)libgrpp_delete_shell(oc_s_so[i]);
        free(oc_s_so);
    }
    
    if(p!=nullptr){
        for(int i=0;i<n_pot;i++)libgrpp_delete_potential(p[i]);
        free(p);
    }
    
    if(so_p!=nullptr){
        for(int i=0;i<n_so;i++)libgrpp_delete_potential(so_p[i]);
        free(so_p);
    }
    
    if(oc_p!=nullptr){
        for(int i=0;i<n_oc;i++)libgrpp_delete_potential(oc_p[i]);
        free(oc_p);
    }
    if(oc_p_so!=nullptr){
        for(int i=0;i<n_oc_so;i++)libgrpp_delete_potential(oc_p_so[i]);
        free(oc_p_so);
    }
    
    if(      origin!=nullptr) delete[]       origin;
    if(   so_origin!=nullptr) delete[]    so_origin;
    if(   oc_origin!=nullptr) delete[]    oc_origin;
    if(oc_origin_so!=nullptr) delete[] oc_origin_so;
    
}


    
#endif



