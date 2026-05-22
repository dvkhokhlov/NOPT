//standart
# include <stdio.h>
# include <omp.h>

//user
# include "blas_link.h"
# include "aldet_rel.h"
# include "matr.h"
# include "timer.h"
# include  "from_hash.h"
# include "common_vars.h"
# include "version.h"
# include "defaults.h"

# define max(a,b)  (((a)<(b))?(b):(a))
# define min(a,b)  (((a)>(b))?(b):(a))

#define NUM_AVX 1

#define TEST_ALDET

aldet_rel_data::aldet_rel_data(){
    
    import_done=0;
    
    n_act  = 0;
    n_el   = 0;
    n_sets = 0;
    print_number = CAS_PRINT_NUMBER_DEFAULT;
    
    do_PT = 0;
    
    Nd = 0;
    f = NULL;
    vec = NULL;
    
    bit = NULL;
    buf = NULL;//more space for +a-b and -a+b ????
    
    
    n_states = NULL;
    coef_r   = NULL;
    coef_i   = NULL;
    
    E_states  = NULL;
    S2        = NULL;
    L2        = NULL;
    P         = NULL;
    E_act     = NULL;
    
    H_diag    = NULL;
    H_diag_appr=NULL;
    Hv_r      = NULL;
    Hv_i      = NULL;
    Hv_buf_r  = NULL;
    Hv_buf_i  = NULL;
    

    J_act  = NULL;
    K_act  = NULL;
    
//     sym_ab = NULL;
    
    e1_ind    = NULL;
    e1_orbs   = NULL;
    e1_sign   = NULL;
    
    n_e1_so   = 0;;
    e1_ind_so = NULL;  
    e1_h_so   = NULL;  
        
    
    e2_orbs   = NULL;
    e2_sign   = NULL;
    e2_ind    = NULL;  
    e2_V      = NULL;  
    
    e3_ind    = NULL;  
    e3_V      = NULL;  
    
    spin_sign = NULL;
    
    F_act_r     = NULL;
    F_act_i     = NULL;
//     F_act_B     = NULL;
    act_INTS_AA = NULL;
    act_INTS_AB = NULL;
    act_INTS_BB = NULL;
    T3_AAA     = NULL;
    T3_AAB     = NULL;
    T3_BBA     = NULL;
    T3_BBB     = NULL;
    T2_AA      = NULL;
    T2_AB      = NULL;
    T2_BB      = NULL;
    T1_A       = NULL;
    T1_B       = NULL;
    
    
}


int aldet_rel_data::get_dim(int ext_n_act, int ext_na, int ext_nb, int ext_n_sets, int ext_print_number){
    
    int status =0;
    // 0 - normal
    // 1 - negative n_el
    // 2 - excessive n_el
    
    n_act        = ext_n_act       ;
    n_el         = ext_na + ext_nb ;
    n_sets       = ext_n_sets      ;
    print_number = ext_print_number;
    
    if(n_el<0    )  status+=1;
    if(n_el>2*n_act)status+=2;
    
    
    if(status)return status;
        
    Nd = (int) std::lround(tgammal(2*n_act+1) / tgammal(n_el+1) / tgammal(2*n_act-n_el+1));
        
    
    if(print_number>Nd)
        print_number=Nd;
    
    f  = new int [n_el * 2*n_act];
    vec   = new int [Nd*(n_el+1)];
    get_factorials(n_el, 2*n_act, f);
    
    get_vec(n_el, 2*n_act, Nd, vec);
    
    bit = new int [2*n_act];
    
    buf = new int [n_el + 2];//more space for +a-b and -a+b ????
    buf[0]=0;
    
    n_states = new int     [n_sets];
    coef_r   = new double *[n_sets];for(int i=0;i<n_sets;i++)coef_r  [i]=NULL;
    coef_i   = new double *[n_sets];for(int i=0;i<n_sets;i++)coef_i  [i]=NULL;
    E_states = new double *[n_sets];for(int i=0;i<n_sets;i++)E_states[i]=NULL;
    S2       = new double *[n_sets];for(int i=0;i<n_sets;i++)S2      [i]=NULL;
    L2       = new double *[n_sets];for(int i=0;i<n_sets;i++)L2      [i]=NULL;
    P        = new double *[n_sets];for(int i=0;i<n_sets;i++)P       [i]=NULL;
    E_act    = new double *[n_sets];for(int i=0;i<n_sets;i++)E_act   [i]=NULL;
    
//     if(na==0)if(nb==0)return 0;
    
    J_act = new double[n_act*n_act*Nd];//????????
    K_act = new double[n_act*n_act*Nd];//????????
    
#ifdef TEST_ALDET
#ifdef DEBUG_WARN
    fprintf(out_stream,"\n\nWARNING current version (%s) has non-optimal ALDet CI engine - diagonal elements are not skiped!\n\n\n", VERSION);
#endif
    n_e1 = n_el*(2*n_act-n_el+1);
#endif
#ifndef TEST_ALDET    
    n_e1 = na_el*(2*n_act-n_el);
#endif
    
    if(n_e1%NUM_AVX)n_e1 = n_e1+NUM_AVX-n_e1%NUM_AVX;
    
    act_INTS_AA = new double[n_act*n_act*n_act*n_act];
    act_INTS_AB = new double[n_act*n_act*n_act*n_act];
    act_INTS_BB = new double[n_act*n_act*n_act*n_act];
    F_act_r     = new double[n_act*n_act*4          ];
    F_act_i     = new double[n_act*n_act*4          ];
//     F_act_B     = F_act_r+n_act*n_act;
    
    e1_ind  = new int[Nd*n_e1];
    e1_orbs = new int[Nd*n_e1];
    e1_sign = new double[Nd*n_e1];
    
//     e1_asm_ints_b = new int[2*Nb*n_e1_b];
#ifdef TEST_ALDET
    n_e2 = n_el*(n_el-1)*(2*n_act-n_el+2)*(2*n_act-n_el+1)/4;
#endif
#ifndef TEST_ALDET
    n_e2 = n_el*(n_el-1)*(2*n_act-n_el)*(2*n_act-n_el-1)/4;
#endif
    
    
    e2_ind  = new int   [Nd*n_e2];
    e2_V    = new double[Nd*n_e2];
    e2_sign = new double[Nd*n_e2];
    e2_orbs = new int   [Nd*n_e2];
    
    
//     a_spin_sign = new double[Na*n_act];
//     b_spin_sign = new double[Nb*n_act];
//     
//     for(int i_CI = 0; i_CI<Na; i_CI++){
//         if(na>0)a_spin_sign[i_CI*n_act+vec_a[i_CI*(na+1)+na]-1]=1;
//         for (int k=na-1; k>0; k--){
//             a_spin_sign[i_CI*n_act+vec_a[i_CI*(na+1)+k]-1]=-a_spin_sign[i_CI*n_act+vec_a[i_CI*(na+1)+k+1]-1];
//         }
//         
//     }
//         
//     for(int j_CI = 0; j_CI<Nb; j_CI++){
//         memset(bit_b,0,n_act*sizeof(int));
//         for (int k=1; k<nb+1; k++) bit_b[vec_b[j_CI*(nb+1)+k]-1] = 1;
//         double sign=1;
//         for(int t=0  ;t<n_act;t++){
//             if(bit_b[t]==0){
//                 b_spin_sign[j_CI*n_act+t]=sign;
//             }
//             else
//                 sign=-sign;
//         }
//         
//     }
    
    return 0;
}

int aldet_rel_data::simple_import_data(double * ext_act_INTS,
                                   double * ext_act_INTS_AB,
                                   double * ext_F_act,
                                   double   ext_E_core){
    
    if(n_act==0){
        printf("ERROR: can not import arrays to aldet_rel_data wth n_act=0\n");
        printf("       it means that you use corrupted version of NOPT\n");
        exit(0);
    }
    
    memcpy(act_INTS_AA, ext_act_INTS   , sizeof(double)*n_act*n_act*n_act*n_act);
    memcpy(act_INTS_AB, ext_act_INTS_AB, sizeof(double)*n_act*n_act*n_act*n_act);
    memcpy(act_INTS_BB, ext_act_INTS   , sizeof(double)*n_act*n_act*n_act*n_act);
    set_zero_matr(F_act_r,n_act*n_act*4);
    for(int i=0;i<n_act;i++)
    for(int j=0;j<n_act;j++)
        F_act_r[i*2*n_act+j]=ext_F_act[i*n_act+j];
    for(int i=0;i<n_act;i++)
    for(int j=0;j<n_act;j++)
        F_act_r[(i+n_act)*2*n_act+j+n_act]=ext_F_act[i*n_act+j];
    
    set_zero_matr(F_act_i,n_act*n_act*4);

    
    E_core      = ext_E_core    ;
    
    import_done=1;
    
    return 0;
    
}

int aldet_rel_data::U_simple_import_data(double * ext_act_INTS_AA,
                                     double * ext_act_INTS_AB,
                                     double * ext_act_INTS_BB,
                                     double * ext_F_act_A,
                                     double * ext_F_act_B,
                                     double   ext_E_core){
    
//     if(n_act==0){
        printf("ERROR: U import doesn't work for rel.");//can not import arrays to aldet_rel_data wth n_act=0\n");
//         printf("       it means that you use corrupted version of NOPT\n");
        exit(0);
//     }
    
    memcpy(act_INTS_AA, ext_act_INTS_AA, sizeof(double)*n_act*n_act*n_act*n_act);
    memcpy(act_INTS_AB, ext_act_INTS_AB, sizeof(double)*n_act*n_act*n_act*n_act);
    memcpy(act_INTS_BB, ext_act_INTS_BB, sizeof(double)*n_act*n_act*n_act*n_act);
    set_zero_matr(F_act_r,n_act*n_act*4);
    for(int i=0;i<n_act;i++)
    for(int j=0;j<n_act;j++)
        F_act_r[i*2*n_act+j]=ext_F_act_A[i*n_act+j];
    for(int i=0;i<n_act;i++)
    for(int j=0;j<n_act;j++)
        F_act_r[(i+n_act)*2*n_act+j+n_act]=ext_F_act_B[i*n_act+j];
    
    
    set_zero_matr(F_act_i,n_act*n_act*4);

    E_core      = ext_E_core    ;
    
    import_done=1;
    
    
    return 0;
    
}

int aldet_rel_data::PT2_alloc(){
    
    if(T3_AAA==nullptr) T3_AAA = new double[n_act*n_act*n_act*n_act*n_act*n_act];
    if(T3_AAB==nullptr) T3_AAB = new double[n_act*n_act*n_act*n_act*n_act*n_act];
    if(T3_BBA==nullptr) T3_BBA = new double[n_act*n_act*n_act*n_act*n_act*n_act];
    if(T3_BBB==nullptr) T3_BBB = new double[n_act*n_act*n_act*n_act*n_act*n_act];
    if(T2_AA ==nullptr) T2_AA  = new double[n_act*n_act*n_act*n_act            ];
    if(T2_AB ==nullptr) T2_AB  = new double[n_act*n_act*n_act*n_act            ];
    if(T2_BB ==nullptr) T2_BB  = new double[n_act*n_act*n_act*n_act            ];
    if(T1_A  ==nullptr) T1_A   = new double[n_act*n_act                        ];
    if(T1_B  ==nullptr) T1_B   = new double[n_act*n_act                        ];
    
    return 0;
    
}



int aldet_rel_data::PT2_import_data(double * ext_T3,
                                double * ext_T3_AB,
                                double * ext_T2,
                                double * ext_T2_AB,
                                double * ext_T1,
                                double   ext_T0){
    
    if(n_act==0){
        printf("ERROR: can not import arrays to aldet_rel_data wth n_act=0\n");
        printf("       it means that you use corrupted version of NOPT\n");
        exit(0);
    }
    
    PT2_alloc();
    
    memcpy(T3_AAA, ext_T3   , sizeof(double)*n_act*n_act*n_act*n_act*n_act*n_act);
    memcpy(T3_AAB, ext_T3_AB, sizeof(double)*n_act*n_act*n_act*n_act*n_act*n_act);
    memcpy(T3_BBA, ext_T3_AB, sizeof(double)*n_act*n_act*n_act*n_act*n_act*n_act);
    memcpy(T3_BBB, ext_T3   , sizeof(double)*n_act*n_act*n_act*n_act*n_act*n_act);
    memcpy(T2_AA , ext_T2   , sizeof(double)*n_act*n_act*n_act*n_act            );
    memcpy(T2_AB , ext_T2_AB, sizeof(double)*n_act*n_act*n_act*n_act            );
    memcpy(T2_BB , ext_T2   , sizeof(double)*n_act*n_act*n_act*n_act            );
    memcpy(T1_A  , ext_T1   , sizeof(double)*n_act*n_act                        );
    memcpy(T1_B  , ext_T1   , sizeof(double)*n_act*n_act                        );
    T0    = ext_T0   ;
    
    do_PT = 1;
    
    return 0;
    
}

int aldet_rel_data::UPT2_import_data(double * ext_T3_AAA,
                                 double * ext_T3_AAB,
                                 double * ext_T3_BBA,
                                 double * ext_T3_BBB,
                                 double * ext_T2_AA,
                                 double * ext_T2_AB,
                                 double * ext_T2_BB,
                                 double * ext_T1_A,
                                 double * ext_T1_B,
                                 double   ext_T0){
    
//     if(n_act==0){
        printf("ERROR: U import doesn't work for rel.");//can not import arrays to aldet_rel_data wth n_act=0\n");
//         printf("       it means that you use corrupted version of NOPT\n");
        exit(0);
//     }
    
    PT2_alloc();
    
    memcpy(T3_AAA, ext_T3_AAA , sizeof(double)*n_act*n_act*n_act*n_act*n_act*n_act);
    memcpy(T3_AAB, ext_T3_AAB , sizeof(double)*n_act*n_act*n_act*n_act*n_act*n_act);
    memcpy(T3_BBA, ext_T3_BBA , sizeof(double)*n_act*n_act*n_act*n_act*n_act*n_act);
    memcpy(T3_BBB, ext_T3_BBB , sizeof(double)*n_act*n_act*n_act*n_act*n_act*n_act);
    memcpy(T2_AA , ext_T2_AA  , sizeof(double)*n_act*n_act*n_act*n_act            );
    memcpy(T2_AB , ext_T2_AB  , sizeof(double)*n_act*n_act*n_act*n_act            );
    memcpy(T2_BB , ext_T2_BB  , sizeof(double)*n_act*n_act*n_act*n_act            );
    memcpy(T1_A  , ext_T1_A   , sizeof(double)*n_act*n_act                        );
    memcpy(T1_B  , ext_T1_B   , sizeof(double)*n_act*n_act                        );
    T0    = ext_T0      ;
    
    do_PT = 1;
    
    return 0;
    
}

int aldet_rel_data::PT_update(){
//     if(na==0)if(nb==0)return 0;
#ifndef TEST_ALDET
    exit(0);
#endif

    for(int i=0; i<n_act;i++)
    for(int j=0; j<n_act;j++)
    for(int k=0; k<n_act;k++)
    for(int l=0; l<n_act;l++)
    for(int m=0; m<n_act;m++)
        T2_AA [((  i*n_act+j         )*n_act+k)*n_act+l]+=
        T3_AAA[((((j*n_act+i)*n_act+m)*n_act+l)*n_act+k)*n_act+m];

    if(T2_BB!=T2_AA)
    for(int i=0; i<n_act;i++)
    for(int j=0; j<n_act;j++)
    for(int k=0; k<n_act;k++)
    for(int l=0; l<n_act;l++)
    for(int m=0; m<n_act;m++)
        T2_BB [((  i*n_act+j         )*n_act+k)*n_act+l]+=
        T3_BBB[((((j*n_act+i)*n_act+m)*n_act+l)*n_act+k)*n_act+m];

    
    if(act_INTS_AB==act_INTS_AA){
        printf("ERROR: use of the external references in aldet_data is deprecated check act_INT_AA and act_INTS_AB\n");
        exit(0);
//         act_INTS_AB_alloced_inside = new double[n_act*n_act*n_act*n_act];
//         act_INTS_AB = act_INTS_AB_alloced_inside;
//         memcpy(act_INTS_AB,act_INTS_AA, n_act*n_act*n_act*n_act*sizeof(double));
    }
    for(int i=0; i<n_act;i++)
    for(int j=0; j<n_act;j++)
    for(int k=0; k<n_act;k++)
    for(int l=0; l<n_act;l++)
        act_INTS_AB[((i*n_act+j)*n_act+k)*n_act+l]=//AABB
        act_INTS_AB[((i*n_act+j)*n_act+k)*n_act+l]+//AABB
        T2_AB     [((i*n_act+k)*n_act+j)*n_act+l];//ABAB
    //AAAA etc. are updated in gen_ext_ind_PT()
    
    gen_ext_ind_PT();
    
    for(int i=0;i<n_act;i++)
    for(int j=0;j<n_act;j++)
        F_act_r[i*2*n_act+j]+=T1_A[i*n_act+j];
    for(int i=0;i<n_act;i++)
    for(int j=0;j<n_act;j++)
        F_act_r[(i+n_act)*2*n_act+j+n_act]+=T1_B[i*n_act+j];
    
        
//     H_diag_calc_PT();
    
//     for(int i=0;i<Nd;i++)
//         H_diag[i]+=T0;
    
    E_core+=T0;
    
    return 0;
}

int aldet_rel_data::SO_update(double * Sx, double * Sy, double * Sz){
    
    for(int i=0;i<n_act;i++)
    for(int j=0;j<n_act;j++)
        F_act_r[i*2*n_act+j+n_act]+= Sy[i*n_act+j]*0.5;
    
    for(int i=0;i<n_act;i++)
    for(int j=0;j<n_act;j++)
        F_act_r[(i+n_act)*2*n_act+j]-=Sy[i*n_act+j]*0.5;
    
    
    for(int i=0;i<n_act;i++)
    for(int j=0;j<n_act;j++)
        F_act_i[i*2*n_act+j] += Sz[i*n_act+j]*0.5;
    
    for(int i=0;i<n_act;i++)
    for(int j=0;j<n_act;j++)
        F_act_i[(i+n_act)*2*n_act+j+n_act]-=Sz[i*n_act+j]*0.5;
    
    
    for(int i=0;i<n_act;i++)
    for(int j=0;j<n_act;j++)
        F_act_i[i*2*n_act+j+n_act]+= Sx[i*n_act+j]*0.5;
    
    for(int i=0;i<n_act;i++)
    for(int j=0;j<n_act;j++)
        F_act_i[(i+n_act)*2*n_act+j]+=Sx[i*n_act+j]*0.5;
    
    
    return 0;
}


int aldet_rel_data::gen_ext_ind(){
    
    if(n_el==0)return 0;
    
    buf[0]=0;
    
    double sign1;
    double sign2;
    double sign3;
    double sign4;
    int i_ext;
    
    if(import_done==0){
        fprintf(out_stream,"ERROR: in function gen_ext_ind()\n");
        fprintf(out_stream,"       external references must be set by aldet_rel_data impport function\n");
        exit(0);
    }
    
    set_zero_matr(J_act , n_act*n_act*Nd);
    set_zero_matr(K_act , n_act*n_act*Nd);
    
    //AA
    for(int i_CI =0; i_CI<Nd; i_CI++){
        i_ext=0;
        memset(bit,0,2*n_act*sizeof(int));
        for (int k=1; k<n_el+1; k++) bit[vec[i_CI*(n_el+1)+k]-1] = 1;
//         printf_occ(i_CI);
//         printf("\n");
        sign1=1.0;
        for(int t=0  ;t<2*n_act;t++)if(bit[t]!=0){
            sign2=1.0;
            bit[t]=0;
            for(int v=0  ;v<2*n_act;v++){
                if(bit[v]==0){
                    bit[v]=1;
                    e1_sign[i_CI*n_e1+i_ext] = sign1*sign2;
                    e1_ind [i_CI*n_e1+i_ext] = get_ind_from_ON(bit, 2*n_act, n_el, f, buf);
                    e1_orbs[i_CI*n_e1+i_ext] = 2*t*n_act+v;
                    bit[v]=0;
//                     printf_occ(e1_ind [i_CI*n_e1+i_ext]);
//                     printf("\n%d %d (%d %d) - %e\n", i_CI, e1_ind [i_CI*n_e1+i_ext],t,v, F_act_r[e1_orbs[i_CI*n_e1+i_ext]]);
                    i_ext++;
                }
                else
                    sign2=-sign2;
            }
            bit[t]=1;
            sign1=-sign1;
        }
        while(i_ext!=n_e1){
            e1_sign[i_CI*n_e1+i_ext] = 0;
            e1_ind [i_CI*n_e1+i_ext] = -1;
            e1_orbs[i_CI*n_e1+i_ext] = e1_orbs    [i_CI*n_e1+i_ext-1];
            i_ext++;
        }
    }
    
//     getchar();
    //AAAA
    for(int i_CI =0; i_CI<Nd; i_CI++){
        
        i_ext=0;
        memset(bit,0,2*n_act*sizeof(int));
        for (int k=1; k<n_el+1; k++) bit[vec[i_CI*(n_el+1)+k]-1] = 1;
        
        for(int t=0  ;t<n_act;t++)if(bit[t]!=0){
            sign1=1.0;
            bit[t]=0;
            for(int u=t+1;u<n_act;u++){
                if(bit[u]!=0){
                    bit[u]=0;
                    for(int v=0  ;v<n_act;v++){
                        if(bit[v]==0){
                            sign2=1.0;
                            bit[v]=1;
                            for(int w=v+1;w<n_act;w++){
                                if(bit[w]==0){
                                    bit[w]=1;
#ifndef TEST_ALDET
                                    if(t!=v)
                                    if(t!=w)
                                    if(u!=v)
                                    if(u!=w)
#endif
                                    {
                                        e2_orbs[i_CI*n_e2+i_ext] =((t*n_act+u)*n_act+v)*n_act+w;
                                        e2_ind [i_CI*n_e2+i_ext] = get_ind_from_ON(bit, 2*n_act, n_el, f, buf);
                                        e2_V   [i_CI*n_e2+i_ext] =(act_INTS_AA[((t*n_act+v)*n_act+u)*n_act+w]-
                                                                   act_INTS_AA[((t*n_act+w)*n_act+u)*n_act+v])*
                                                                   sign1*sign2;
                                        e2_sign[i_ext]= sign1*sign2;
                                        i_ext++;
                                    }
                                    bit[w]=0;
                                }
                                else
                                    sign2=-sign2;
                            }
                            bit[v]=0;
                        }
                    }
                    bit[u]=1;
                    sign1=-sign1;
                }
            }
            bit[t]=1;
        }
        for(int t=0  ;t<n_act;t++)if(bit[t+n_act]!=0){
            sign1=1.0;
            bit[t+n_act]=0;
            for(int u=t+1;u<n_act;u++){
                if(bit[u+n_act]!=0){
                    bit[u+n_act]=0;
                    for(int v=0  ;v<n_act;v++){
                        if(bit[v+n_act]==0){
                            sign2=1.0;
                            bit[v+n_act]=1;
                            for(int w=v+1;w<n_act;w++){
                                if(bit[w+n_act]==0){
                                    bit[w+n_act]=1;
#ifndef TEST_ALDET
                                    if(t!=v)
                                    if(t!=w)
                                    if(u!=v)
                                    if(u!=w)
#endif
                                    {
                                        e2_orbs[i_CI*n_e2+i_ext] =((t*n_act+u)*n_act+v)*n_act+w;
                                        e2_ind [i_CI*n_e2+i_ext] = get_ind_from_ON(bit, 2*n_act, n_el, f, buf);
                                        e2_V   [i_CI*n_e2+i_ext] =(act_INTS_BB[((t*n_act+v)*n_act+u)*n_act+w]-
                                                                   act_INTS_BB[((t*n_act+w)*n_act+u)*n_act+v])*
                                                                   sign1*sign2;
                                        e2_sign[i_ext]= sign1*sign2;
                                        i_ext++;
                                    }
                                    bit[w+n_act]=0;
                                }
                                else
                                    sign2=-sign2;
                            }
                            bit[v+n_act]=0;
                        }
                    }
                    bit[u+n_act]=1;
                    sign1=-sign1;
                }
            }
            bit[t+n_act]=1;
        }
        sign1=1.0;
        for(int t=0  ;t<n_act;t++)if(bit[t]!=0){
            bit[t]=0;
            sign2=1.0;
            for(int u=0;u<n_act;u++){
                if(bit[u+n_act]!=0){
                    bit[u+n_act]=0;
                    sign3=1.0;
                    for(int v=0  ;v<n_act;v++){
                        if(bit[v]==0){
                            sign4=1.0;
                            bit[v]=1;
                            for(int w=0;w<n_act;w++){
                                if(bit[w+n_act]==0){
                                    bit[w+n_act]=1;
#ifndef TEST_ALDET
                                    if(t!=v)
                                    if(t!=w)
                                    if(u!=v)
                                    if(u!=w)
#endif
                                    {
                                        e2_orbs[i_CI*n_e2+i_ext] =((t*n_act+u)*n_act+v)*n_act+w;
                                        e2_ind [i_CI*n_e2+i_ext] = get_ind_from_ON(bit, 2*n_act, n_el, f, buf);
                                        e2_V   [i_CI*n_e2+i_ext] =(act_INTS_AB[((t*n_act+v)*n_act+u)*n_act+w])*
                                                                   sign1*sign2*sign3*sign4;
                                        e2_sign[i_ext]= sign1*sign2*sign3*sign4;
                                        i_ext++;
                                    }
                                    bit[w+n_act]=0;
                                }
                                else
                                    sign4=-sign4;
                            }
                            bit[v]=0;
                        }
                        else
                            sign3=-sign3;
                    }
                    bit[u+n_act]=1;
                    sign2=-sign2;
                }
            }
            bit[t]=1;
            sign1=-sign1;
        }
        while(i_ext<n_e2){
            e2_sign[i_CI*n_e2+i_ext] = 0;
            e2_ind [i_CI*n_e2+i_ext] = -1;
            e2_orbs[i_CI*n_e2+i_ext] = -1;
            e2_V   [i_CI*n_e2+i_ext] = 0;
            i_ext++;
        }
    }
    
    return 0;
    
}

int aldet_rel_data::gen_ext_ind_PT(){
    
    buf[0]=0;
    
    double sign1;
    double sign2;
    double sign3;
    double sign4;
    double t_sign;
    
    int i_ext;
    
//     
//     set_zero_matr(J_act_a , n_act*n_act*Na);
//     set_zero_matr(J_act_b , n_act*n_act*Nb);
//     set_zero_matr(K_act_a , n_act*n_act*Na);
//     set_zero_matr(K_act_b , n_act*n_act*Nb);
//     
    //this will be optimized
    n_e1 = n_el*(2*n_act-n_el+1);
    if(n_e1%NUM_AVX)n_e1 = n_e1+NUM_AVX-n_e1%NUM_AVX;
    
    if(e1_ind  !=NULL) delete[] e1_ind ;
    if(e1_orbs !=NULL) delete[] e1_orbs;
    if(e1_sign !=NULL) delete[] e1_sign;
    
    e1_ind  = new int[Nd*n_e1];
    e1_orbs = new int[Nd*n_e1];
    e1_sign = new double[Nd*n_e1];
    
    n_e2 = n_el*(n_el-1)*(2*n_act-n_el+2)*(2*n_act-n_el+1)/4;
    if(n_e2%NUM_AVX)n_e2 = n_e2+NUM_AVX-n_e2%NUM_AVX;
    
    if(e2_sign  !=NULL) delete[] e2_sign ;
    if(e2_orbs  !=NULL) delete[] e2_orbs ;
    if(e2_ind != NULL) delete[] e2_ind;
    if(e2_V   != NULL) delete[] e2_V  ;
    
    
    e2_sign = new double[Nd*n_e2];
    e2_orbs = new int[Nd*n_e2];
    e2_ind  = new int[Nd*n_e2];
    e2_V = new double[Nd*n_e2];
    
    n_e3 = n_el*(n_el-1)*(n_el-2)*(2*n_act-n_el+3)*(2*n_act-n_el+2)*(2*n_act-n_el+1)/36;
    
    if(e3_ind != NULL) delete[] e3_ind;
    if(e3_V   != NULL) delete[] e3_V  ;
    
    
    e3_ind  = new int[Nd*n_e3];
    e3_V = new double[Nd*n_e3];
    
    //AA
    for(int i_CI =0; i_CI<Nd; i_CI++){
        i_ext=0;
        memset(bit,0,2*n_act*sizeof(int));
        for (int k=1; k<n_el+1; k++) bit[vec[i_CI*(n_el+1)+k]-1] = 1;
//         printf_occ(i_CI);
//         printf("\n");
        sign1=1.0;
        for(int t=0  ;t<2*n_act;t++)if(bit[t]!=0){
            sign2=1.0;
            bit[t]=0;
            for(int v=0  ;v<2*n_act;v++){
                if(bit[v]==0){
                    bit[v]=1;
                    e1_sign[i_CI*n_e1+i_ext] = sign1*sign2;
                    e1_ind [i_CI*n_e1+i_ext] = get_ind_from_ON(bit, 2*n_act, n_el, f, buf);
                    e1_orbs[i_CI*n_e1+i_ext] = t*2*n_act+v;
                    bit[v]=0;
//                     printf_occ(e1_ind [i_CI*n_e1+i_ext]);
//                     printf("\n%d %d (%d %d) - %e\n", i_CI, e1_ind [i_CI*n_e1+i_ext],t,v, F_act_r[e1_orbs[i_CI*n_e1+i_ext]]);
                    i_ext++;
                }
                else
                    sign2=-sign2;
            }
            bit[t]=1;
            sign1=-sign1;
        }
        while(i_ext!=n_e1){
            e1_sign[i_CI*n_e1+i_ext] = 0;
            e1_ind [i_CI*n_e1+i_ext] = -1;
            e1_orbs[i_CI*n_e1+i_ext] = e1_orbs    [i_CI*n_e1+i_ext-1];
            i_ext++;
        }
    }
    //AAAA
    for(int i_CI =0; i_CI<Nd; i_CI++){
        
        i_ext=0;
        memset(bit,0,2*n_act*sizeof(int));
        for (int k=1; k<n_el+1; k++) bit[vec[i_CI*(n_el+1)+k]-1] = 1;
        
        for(int t=0  ;t<n_act;t++)if(bit[t]!=0){
            sign1=1.0;
            bit[t]=0;
            for(int u=t+1;u<n_act;u++){
                if(bit[u]!=0){
                    bit[u]=0;
                    for(int v=0  ;v<n_act;v++){
                        if(bit[v]==0){
                            sign2=1.0;
                            bit[v]=1;
                            for(int w=v+1;w<n_act;w++){
                                if(bit[w]==0){
                                    bit[w]=1;
#ifndef TEST_ALDET
                                    if(t!=v)
                                    if(t!=w)
                                    if(u!=v)
                                    if(u!=w)
#endif
                                    {
                                        e2_orbs[i_CI*n_e2+i_ext] =((t*n_act+u)*n_act+v)*n_act+w;
                                        e2_ind [i_CI*n_e2+i_ext] = get_ind_from_ON(bit, 2*n_act, n_el, f, buf);
                                        e2_V   [i_CI*n_e2+i_ext] =(act_INTS_AA[((t*n_act+v)*n_act+u)*n_act+w]-
                                                                   act_INTS_AA[((t*n_act+w)*n_act+u)*n_act+v]+
                                                                   T2_AA     [((u*n_act+t)*n_act+v)*n_act+w])*
                                                                   sign1*sign2;
                                        e2_sign[i_ext]= sign1*sign2;
                                        i_ext++;
                                    }
                                    bit[w]=0;
                                }
                                else
                                    sign2=-sign2;
                            }
                            bit[v]=0;
                        }
                    }
                    bit[u]=1;
                    sign1=-sign1;
                }
            }
            bit[t]=1;
        }
        for(int t=0  ;t<n_act;t++)if(bit[t+n_act]!=0){
            sign1=1.0;
            bit[t+n_act]=0;
            for(int u=t+1;u<n_act;u++){
                if(bit[u+n_act]!=0){
                    bit[u+n_act]=0;
                    for(int v=0  ;v<n_act;v++){
                        if(bit[v+n_act]==0){
                            sign2=1.0;
                            bit[v+n_act]=1;
                            for(int w=v+1;w<n_act;w++){
                                if(bit[w+n_act]==0){
                                    bit[w+n_act]=1;
#ifndef TEST_ALDET
                                    if(t!=v)
                                    if(t!=w)
                                    if(u!=v)
                                    if(u!=w)
#endif
                                    {
                                        e2_orbs[i_CI*n_e2+i_ext] =((t*n_act+u)*n_act+v)*n_act+w;
                                        e2_ind [i_CI*n_e2+i_ext] = get_ind_from_ON(bit, 2*n_act, n_el, f, buf);
                                        e2_V   [i_CI*n_e2+i_ext] =(act_INTS_BB[((t*n_act+v)*n_act+u)*n_act+w]-
                                                                   act_INTS_BB[((t*n_act+w)*n_act+u)*n_act+v]+
                                                                   T2_BB     [((u*n_act+t)*n_act+v)*n_act+w])*
                                                                   sign1*sign2;
                                        e2_sign[i_ext]= sign1*sign2;
                                        i_ext++;
                                    }
                                    bit[w+n_act]=0;
                                }
                                else
                                    sign2=-sign2;
                            }
                            bit[v+n_act]=0;
                        }
                    }
                    bit[u+n_act]=1;
                    sign1=-sign1;
                }
            }
            bit[t+n_act]=1;
        }
        sign1=1.0;
        for(int t=0  ;t<n_act;t++)if(bit[t]!=0){
            bit[t]=0;
            sign2=1.0;
            for(int u=0;u<n_act;u++){
                if(bit[u+n_act]!=0){
                    bit[u+n_act]=0;
                    sign3=1.0;
                    for(int v=0  ;v<n_act;v++){
                        if(bit[v]==0){
                            sign4=1.0;
                            bit[v]=1;
                            for(int w=0;w<n_act;w++){
                                if(bit[w+n_act]==0){
                                    bit[w+n_act]=1;
#ifndef TEST_ALDET
                                    if(t!=v)
                                    if(t!=w)
                                    if(u!=v)
                                    if(u!=w)
#endif
                                    {
                                        e2_orbs[i_CI*n_e2+i_ext] =((t*n_act+u)*n_act+v)*n_act+w;
                                        e2_ind [i_CI*n_e2+i_ext] = get_ind_from_ON(bit, 2*n_act, n_el, f, buf);
                                        e2_V   [i_CI*n_e2+i_ext] =(act_INTS_AB[((t*n_act+v)*n_act+u)*n_act+w])*
                                                                   sign1*sign2*sign3*sign4;
                                        e2_sign[i_ext]= sign1*sign2*sign3*sign4;
                                        i_ext++;
                                    }
                                    bit[w+n_act]=0;
                                }
                                else
                                    sign4=-sign4;
                            }
                            bit[v]=0;
                        }
                        else
                            sign3=-sign3;
                    }
                    bit[u+n_act]=1;
                    sign2=-sign2;
                }
            }
            bit[t]=1;
            sign1=-sign1;
        }
        while(i_ext<n_e2){
            e2_sign[i_CI*n_e2+i_ext] = 0;
            e2_ind [i_CI*n_e2+i_ext] = -1;
            e2_orbs[i_CI*n_e2+i_ext] = -1;
            e2_V   [i_CI*n_e2+i_ext] = 0;
            i_ext++;
        }
    }
//     //AAAAAA
//     i_ext=0;
    for(int i_CI =0; i_CI<Nd; i_CI++){//fprintf(stderr,"i_CI=%d  (thread %d)        \r",i_CI, i_th);
        
        i_ext=0;
        memset(bit,0,2*n_act*sizeof(int));
        for (int k=1; k<n_el+1; k++) bit[vec[i_CI*(n_el+1)+k]-1] = 1;
        
        double t_sign=1;
        for(int t=0  ;t<n_act;t++){
            if(bit[t]!=0){
                bit[t]=0;
                sign1=t_sign;
                for(int u=t+1;u<n_act;u++)if(bit[u]!=0){
                    bit[u]=0;
                            sign3=sign1;//1;
                            for(int w=u+1/*0*/;w<n_act;w++)if(bit[w]!=0){
                                bit[w]=0;
                                
                                
                    sign2=1;
                    for(int v=0;v<n_act;v++){
                        if(bit[v]==0){
                            bit[v]=1;
                                for(int x=v+1/*0*/;x<n_act;x++)if(bit[x]==0){
                                    sign4=1;
                                    for(int y=x+1;y<n_act;y++){
                                        if(bit[y]==0){
                                            bit[x]=1;
                                            bit[y]=1;
                                            // change order from a--+++- to ---+++
                                            // regular order a+(v)a+(x)a+(y)a(w)a(u)a(t);v>x>y,t>u>v
                                            e3_ind [i_CI*n_e3+i_ext] = get_ind_from_ON(bit, 2*n_act, n_el, f, buf);
                                            e3_V   [i_CI*n_e3+i_ext] =(-T3_AAA[((((t*n_act+u)*n_act+v)*n_act+y)*n_act+x)*n_act+w]//;
                                                                       +T3_AAA[((((t*n_act+w)*n_act+v)*n_act+y)*n_act+x)*n_act+u]
                                                                       -T3_AAA[((((u*n_act+w)*n_act+v)*n_act+y)*n_act+x)*n_act+t]
                                                                       +T3_AAA[((((t*n_act+u)*n_act+x)*n_act+y)*n_act+v)*n_act+w]
                                                                       -T3_AAA[((((t*n_act+w)*n_act+x)*n_act+y)*n_act+v)*n_act+u]
                                                                       +T3_AAA[((((u*n_act+w)*n_act+x)*n_act+y)*n_act+v)*n_act+t]
                                                                       -T3_AAA[((((t*n_act+u)*n_act+y)*n_act+x)*n_act+v)*n_act+w]
                                                                       +T3_AAA[((((t*n_act+w)*n_act+y)*n_act+x)*n_act+v)*n_act+u]
                                                                       -T3_AAA[((((u*n_act+w)*n_act+y)*n_act+x)*n_act+v)*n_act+t])*
                                                                       t_sign*sign1*sign2*sign3*sign4;
                                            i_ext++;
                                            bit[x]=0;
                                            bit[y]=0;    
                                        }
                                        else
                                            sign4=-sign4;
                                    }
                                }
                            bit[v]=0;                            }
                        else
                            sign2=-sign2;
                    }
                                bit[w]=1;
                                sign3=-sign3;                                }
                    bit[u]=1;
                    sign1=-sign1;
                }
                bit[t]=1;
                t_sign=-t_sign;
            }
        }
        for(int t=0  ;t<n_act;t++){
            if(bit[t+n_act]!=0){
                bit[t+n_act]=0;
                sign1=t_sign;
                for(int u=t+1;u<n_act;u++)if(bit[u+n_act]!=0){
                    bit[u+n_act]=0;
                            sign3=sign1;//1;////HERE!!!!!!!!!!!!!!!!!!!!
                            for(int w=u+1/*0*/;w<n_act;w++)if(bit[w+n_act]!=0){
                                bit[w+n_act]=0;
                                
                                
                    sign2=1;
                    for(int v=0;v<n_act;v++){
                        if(bit[v+n_act]==0){
                            bit[v+n_act]=1;
                                for(int x=v+1/*0*/;x<n_act;x++)if(bit[x+n_act]==0){
                                    sign4=1;
                                    for(int y=x+1;y<n_act;y++){
                                        if(bit[y+n_act]==0){
                                            bit[x+n_act]=1;
                                            bit[y+n_act]=1;
                                            // change order from a--+++- to ---+++
                                            // regular order a+(v)a+(x)a+(y)a(w)a(u)a(t);v>x>y,t>u>v
                                            e3_ind [i_CI*n_e3+i_ext] = get_ind_from_ON(bit, 2*n_act, n_el, f, buf);
                                            e3_V   [i_CI*n_e3+i_ext] =(-T3_AAA[((((t*n_act+u)*n_act+v)*n_act+y)*n_act+x)*n_act+w]//;
                                                                       +T3_AAA[((((t*n_act+w)*n_act+v)*n_act+y)*n_act+x)*n_act+u]
                                                                       -T3_AAA[((((u*n_act+w)*n_act+v)*n_act+y)*n_act+x)*n_act+t]
                                                                       +T3_AAA[((((t*n_act+u)*n_act+x)*n_act+y)*n_act+v)*n_act+w]
                                                                       -T3_AAA[((((t*n_act+w)*n_act+x)*n_act+y)*n_act+v)*n_act+u]
                                                                       +T3_AAA[((((u*n_act+w)*n_act+x)*n_act+y)*n_act+v)*n_act+t]
                                                                       -T3_AAA[((((t*n_act+u)*n_act+y)*n_act+x)*n_act+v)*n_act+w]
                                                                       +T3_AAA[((((t*n_act+w)*n_act+y)*n_act+x)*n_act+v)*n_act+u]
                                                                       -T3_AAA[((((u*n_act+w)*n_act+y)*n_act+x)*n_act+v)*n_act+t])*
                                                                       t_sign*sign1*sign2*sign3*sign4;
                                            i_ext++;
                                            bit[x+n_act]=0;
                                            bit[y+n_act]=0;    
                                        }
                                        else
                                            sign4=-sign4;
                                    }
                                }
                            bit[v+n_act]=0;                            }
                        else
                            sign2=-sign2;
                    }
                                bit[w+n_act]=1;
                                sign3=-sign3;                                }
                    bit[u+n_act]=1;
                    sign1=-sign1;
                }
                bit[t+n_act]=1;
                t_sign=-t_sign;
            }
        }
        for(int t=0  ;t<n_act;t++)if(bit[t]!=0){
            sign1=1.0;
            bit[t]=0;
            for(int u=t+1;u<n_act;u++){
                if(bit[u]!=0){
                    bit[u]=0;
                    for(int v=0  ;v<n_act;v++){
                        if(bit[v]==0){
                            sign2=1.0;
                            bit[v]=1;
                            for(int w=v+1;w<n_act;w++){
                                if(bit[w]==0){
                                    bit[w]=1;
                                    sign3=1.0;
                                    for(int x=0  ;x<n_act;x++){
                                        if(bit[x+n_act]!=0){
                                            sign4=1.0;
                                            bit[x+n_act]=0;
                                            for(int y=0;y<n_act;y++){
                                                if(bit[y+n_act]==0){
                                                    bit[y+n_act]=1;
                                                    e3_ind [i_CI*n_e3+i_ext] = get_ind_from_ON(bit, 2*n_act, n_el, f, buf);
                                                    e3_V   [i_CI*n_e3+i_ext] =T3_AAB[((((t*n_act+u)*n_act+v)*n_act+w)*n_act+x)*n_act+y]*
                                                                              sign1*sign2*sign3*sign4;
                                                    i_ext++;
                                                    bit[y+n_act]=0;
                                                }
                                                else
                                                    sign4=-sign4;
                                            }
                                            bit[x+n_act]=1;
                                            sign3=-sign3;
                                        }
                                    }
                                    bit[w]=0;
                                }
                                else
                                    sign2=-sign2;
                            }
                            bit[v]=0;
                        }
                    }
                    bit[u]=1;
                    sign1=-sign1;
                }
            }
            bit[t]=1;
        }
        for(int t=0  ;t<n_act;t++)if(bit[t+n_act]!=0){
            sign1=1.0;
            bit[t+n_act]=0;
            for(int u=t+1;u<n_act;u++){
                if(bit[u+n_act]!=0){
                    bit[u+n_act]=0;
                    for(int v=0  ;v<n_act;v++){
                        if(bit[v+n_act]==0){
                            sign2=1.0;
                            bit[v+n_act]=1;
                            for(int w=v+1;w<n_act;w++){
                                if(bit[w+n_act]==0){
                                    bit[w+n_act]=1;
                                    sign3=1.0;
                                    for(int x=0  ;x<n_act;x++){
                                        if(bit[x]!=0){
                                            sign4=1.0;
                                            bit[x]=0;
                                            for(int y=0;y<n_act;y++){
                                                if(bit[y]==0){
                                                    bit[y]=1;
                                                    e3_ind [i_CI*n_e3+i_ext] = get_ind_from_ON(bit, 2*n_act, n_el, f, buf);
                                                    e3_V   [i_CI*n_e3+i_ext] =T3_AAB[((((t*n_act+u)*n_act+v)*n_act+w)*n_act+x)*n_act+y]*
                                                                              sign1*sign2*sign3*sign4;
                                                    i_ext++;
                                                    bit[y]=0;
                                                }
                                                else
                                                    sign4=-sign4;
                                            }
                                            bit[x]=1;
                                            sign3=-sign3;
                                        }
                                    }
                                    bit[w+n_act]=0;
                                }
                                else
                                    sign2=-sign2;
                            }
                            bit[v+n_act]=0;
                        }
                    }
                    bit[u+n_act]=1;
                    sign1=-sign1;
                }
            }
            bit[t+n_act]=1;
        }
        
        while(i_ext<n_e3){
            e3_V  [i_CI*n_e3+i_ext] = 0;
            e3_ind[i_CI*n_e3+i_ext] = -1;
            i_ext++;
        }
    }
    
//     fprintf(out_stream,"bbb%d %d\n", i_ext, n_e3_b*Nb);
//     getchar();
    
    return 0;
    
}

int aldet_rel_data::gen_ext_ind_SO(){
    
//     buf[0]=0;
//     
//     double sign1;
//     double sign2;
//     int i_ext;
//     
// //     
// //     set_zero_matr(J_act_a , n_act*n_act*Na);
// //     set_zero_matr(J_act_b , n_act*n_act*Nb);
// //     set_zero_matr(K_act_a , n_act*n_act*Na);
// //     set_zero_matr(K_act_b , n_act*n_act*Nb);
// //     
//     //this will be optimized
//     n_e1_so = n_el*n_act;
//     if(n_e1_so%NUM_AVX)n_e1_so = n_e1_so+NUM_AVX-n_e1_so%NUM_AVX;
//     
//     if(e1_ind_so != NULL) delete[]e1_ind_so     ;
//     if(e1_h_so   != NULL) delete[]e1_h_so       ;
//     
//     e1_ind_so = new int[Nd*n_e1_so];
//     e1_h_so   = new double[Nd*n_e1_so];
//     
//     //AA
//         for(int i_CI =0; i_CI<Nd; i_CI++){
//         i_ext=0;
//         memset(bit,0,2*n_act*sizeof(int));
//         for (int k=1; k<n_el+1; k++) bit[vec[i_CI*(n_el+1)+k]-1] = 1;
// //         printf_occ(i_CI);
// //         printf("\n");
//         sign1=1.0;
//         for(int t=0  ;t<2*n_act;t++)if(bit[t]!=0){
//             sign2=1.0;
//             bit[t]=0;
//             for(int v=0  ;v<2*n_act;v++){
//                 if(bit[v]==0){
//                     bit[v]=1;
//                     e1_sign[i_CI*n_e1+i_ext] = sign1*sign2;
//                     e1_ind [i_CI*n_e1+i_ext] = get_ind_from_ON(bit, 2*n_act, n_el, f, buf);
//                     e1_orbs[i_CI*n_e1+i_ext] = 2*t*n_act+v;
//                     bit[v]=0;
// //                     printf_occ(e1_ind [i_CI*n_e1+i_ext]);
// //                     printf("\n%d %d (%d %d) - %e\n", i_CI, e1_ind [i_CI*n_e1+i_ext],t,v, F_act_r[e1_orbs[i_CI*n_e1+i_ext]]);
//                     i_ext++;
//                 }
//                 else
//                     sign2=-sign2;
//             }
//             bit[t]=1;
//             sign1=-sign1;
//         }
//         while(i_ext!=n_e1){
//             e1_sign[i_CI*n_e1+i_ext] = 0;
//             e1_ind [i_CI*n_e1+i_ext] = -1;
//             e1_orbs[i_CI*n_e1+i_ext] = e1_orbs    [i_CI*n_e1+i_ext-1];
//             i_ext++;
//         }
//     }
//     
    return 0;
    
}

int aldet_rel_data::init_zero_vec(int n_s, int i_set){
    
    if(n_s==-1) n_s=Nd;
    n_states[i_set] = n_s;
    
    if(coef_r  [i_set]==NULL) coef_r  [i_set] = new double[Nd*n_s];
    if(coef_i  [i_set]==NULL) coef_i  [i_set] = new double[Nd*n_s];
    if(E_states[i_set]==NULL) E_states[i_set] = new double[   n_s];
    if(S2      [i_set]==NULL) S2      [i_set] = new double[   n_s];
    if(L2      [i_set]==NULL) L2      [i_set] = new double[   n_s];
    if(P       [i_set]==NULL) P       [i_set] = new double[   n_s];
//     if(sym_ab  [i_set]==NULL) sym_ab  [i_set] = new double[   n_s];
    
    
    set_zero_matr(coef_r[i_set], Nd*n_s);
    set_zero_matr(coef_i[i_set], Nd*n_s);
    
    
    return 0;
}




int aldet_rel_data::H_full_calc_and_diag(int i_set){
    
    if(n_states[i_set] != Nd){
        if(coef_r  [i_set]!=NULL) delete[] coef_r  [i_set];
        if(coef_i  [i_set]!=NULL) delete[] coef_i  [i_set];
        if(E_states[i_set]!=NULL) delete[] E_states[i_set];
        if(S2      [i_set]!=NULL) delete[] S2      [i_set];
        if(L2      [i_set]!=NULL) delete[] L2      [i_set];
        if(P       [i_set]!=NULL) delete[] P       [i_set];
        
        coef_r  [i_set] = new double[Nd*Nd];
        coef_i  [i_set] = new double[Nd*Nd];
        E_states[i_set] = new double[   Nd];
        S2      [i_set] = new double[   Nd];
        L2      [i_set] = new double[   Nd];
        P       [i_set] = new double[   Nd];
        
        n_states[i_set] = Nd;
    }
    
    set_zero_matr(coef_r[i_set], Nd*Nd);
    set_zero_matr(coef_i[i_set], Nd*Nd);
    
    H_full_calc  (coef_r[i_set]);
    H_full_calc_i(coef_i[i_set]);

    lapack_herm_diag(coef_r[i_set],coef_i[i_set], E_states[i_set], Nd);
    
    transpose(coef_r[i_set],Nd,Nd);
    transpose(coef_i[i_set],Nd,Nd);
    
    
    return 0;
}

int aldet_rel_data::H_full_calc(double *H){
    
    set_zero_matr(H,Nd*Nd);
    
    int j_CI;
    double s;
    double h;
    
    for(int i_CI=0;i_CI<Nd;i_CI++)H[i_CI*Nd+i_CI]+=E_core;
    
    for(int i_CI=0;i_CI<Nd;i_CI++)
    for(int k = 0;k<n_e1;k++){
        j_CI = e1_ind[i_CI*n_e1+k];
        s = e1_sign[i_CI*n_e1+k];
        h = F_act_r[e1_orbs[i_CI*n_e1+k]];
        
        if(j_CI!=-1)H[i_CI*Nd+j_CI]+=s*h;    
        
    }
    
    for(int i_CI=0;i_CI<Nd;i_CI++)
    for(int k = 0;k<n_e2;k++){
        j_CI = e2_ind[i_CI*n_e2+k];
        h = e2_V[i_CI*n_e2+k];
        
        if(j_CI!=-1)H[i_CI*Nd+j_CI]+=h;    
        
    }
    if(do_PT)
    for(int i_CI=0;i_CI<Nd;i_CI++)
    for(int k = 0;k<n_e3;k++){
        j_CI = e3_ind[i_CI*n_e3+k];
        h = e3_V[i_CI*n_e3+k];
        
        if(j_CI!=-1)H[i_CI*Nd+j_CI]+=h;    
        
    }
    
    
    
    
    return 0;
    
    
}

int aldet_rel_data::H_full_calc_i(double *Hi){
    
    set_zero_matr(Hi,Nd*Nd);
    
    int j_CI;
    double s;
    double h;
    
    
    for(int i_CI=0;i_CI<Nd;i_CI++)
    for(int k = 0;k<n_e1;k++){
        j_CI = e1_ind[i_CI*n_e1+k];
        s = e1_sign[i_CI*n_e1+k];
        h = F_act_i[e1_orbs[i_CI*n_e1+k]];
        
        if(j_CI!=-1)Hi[i_CI*Nd+j_CI]+=s*h;    
        
    }
    
    return 0;
    
}


int aldet_rel_data::printf_occ(int i_CI){
    
    int * bit2 = new int[2*n_act];
    memset(bit2,0,2*n_act*sizeof(int));
    
    for (int k=1; k<n_el+1; k++) bit2[vec[i_CI*(n_el+1)+k]-1] = 1;
    for (int t=0; t<n_act;t++)
        if(bit2[t])fprintf(out_stream,"1");
        else        fprintf(out_stream,"0");
    fprintf(out_stream," | ");
    for (int t=0; t<n_act;t++)
        if(bit2[t+n_act])fprintf(out_stream,"1");
        else            fprintf(out_stream,"0");
    
    delete[] bit2;
        
    return 0;
}

int aldet_rel_data::calc_spin_matr(double * S2, double * Sz, double * C_r, double * C_i,int n_s, int ld){
    
    double * S_M_c_r = new double[Nd*n_s];
    double * S_M_c_i = new double[Nd*n_s];
    set_zero_matr(S_M_c_r, Nd*n_s);
    set_zero_matr(S_M_c_i, Nd*n_s);
    
    
    
    int * bit_a=bit;
    int * bit_b=bit+n_act;
    double sign;
    for(int i_CI = 0; i_CI<Nd; i_CI++){
        memset(bit,0,2*n_act*sizeof(int));
        for (int k=1; k<n_el+1; k++) bit[vec[i_CI*(n_el+1)+k]-1] = 1;
        
        int Sz_i=0;
        for(int t=0  ;t<n_act;t++)if(bit_a[t])Sz_i++;
        for(int t=0  ;t<n_act;t++)if(bit_b[t])Sz_i--;
        double Sz_2=Sz_i;
        Sz_2=Sz_2/2.0;
        Sz_2=Sz_2*(Sz_2-1.0);
        for(int i_s=0;i_s<n_s;i_s++)
            S_M_c_r[i_CI*n_s+i_s]+= C_r[i_CI*Nd+i_s]*Sz_2;
        for(int i_s=0;i_s<n_s;i_s++)
            S_M_c_i[i_CI*n_s+i_s]+= C_i[i_CI*Nd+i_s]*Sz_2;
                    
    }
    
    
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                         n_s,n_s,Nd,1.0,
                         C_r,Nd,
                         S_M_c_r,n_s,0.0,
                         S2    ,n_s);
    
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                         n_s,n_s,Nd,1.0,
                         C_i,Nd,
                         S_M_c_i,n_s,1.0,
                         S2    ,n_s);
    
    
    set_zero_matr(S_M_c_r, Nd*n_s);
    set_zero_matr(S_M_c_i, Nd*n_s);
    
    for(int i_CI = 0; i_CI<Nd; i_CI++){
        
        memset(bit,0,2*n_act*sizeof(int));
        for (int k=1; k<n_el+1; k++) bit[vec[i_CI*(n_el+1)+k]-1] = 1;
                
        for(int t=0  ;t<n_act;t++)if(bit_a[t]!=0)if(bit_b[t]==0){
            bit_a[t]=0;
            bit_b[t]=1;
            sign=1.0;
            for(int u=t+1;u<t+n_act;u++)if(bit[u])sign=-sign;
            auto i_CIext = get_ind_from_ON(bit, 2*n_act, n_el, f, buf);
            for(int i_s=0;i_s<n_s;i_s++)
                S_M_c_r[i_CIext*n_s+i_s]+= C_r[i_CI*ld+i_s]*sign;
            for(int i_s=0;i_s<n_s;i_s++)
                S_M_c_i[i_CIext*n_s+i_s]+= C_i[i_CI*ld+i_s]*sign;
            bit_a[t]=1;
            bit_b[t]=0;
            
        }
            
    }
    
    
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                         n_s,n_s,Nd,1.0,
                         S_M_c_r,n_s,
                         S_M_c_r,n_s,1.0,
                         S2    ,n_s);
    
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                         n_s,n_s,Nd,1.0,
                         S_M_c_i,n_s,
                         S_M_c_i,n_s,1.0,
                         S2    ,n_s);
    
//     PrintMatr(O, n_s,n_s,0);
//     for(int i_s=0;i_s<n_s;i_s++)
//         fprintf(out_stream,"%d - %e\n",i_s,O[i_s*(n_s+1)]);
//     getchar();
//     delete[] faM   ;
//     delete[] fbM   ;
    delete[] S_M_c_r ;
    delete[] S_M_c_i ;
    
    return 0;
}


int aldet_rel_data::print_states(int i_s, int n_s, int print){
    
    print_number=3;
    double * S2_loc = new double[n_s*n_s];
    double * Sz_loc = new double[n_s*n_s];
//     double * L2_loc;
//     double * P_loc;
    double * c_2 = new double[Nd];
    int * max_coef = new int[12];
    calc_spin_matr(S2_loc,Sz_loc,coef_r[i_s],coef_i[i_s],n_s,n_states[i_s]);
//     if(LINEAR){
//         L2_loc = new double[n_s*n_s];
//         P_loc  = new double[n_s*n_s];
//         calc_L2_matr(L2_loc,i_s,n_s,n_states[i_s]);
//         calc_P_matr ( P_loc,i_s,n_s,n_states[i_s]);
//     }
        
    
    for(int i=0; i<n_s; i++){
//         S2[i_s][i]=S2_loc[i*n_s+i];
        if(print)fprintf(out_stream,"State %d  E  = % 18.10f S^2 = %.2f",i,E_states[i_s][i],S2_loc[i*n_s+i]);
//         if(LINEAR){
//             L2[i_s][i]=L2_loc[i*n_s+i];
//             P [i_s][i]= P_loc[i*n_s+i];
//             if(print)fprintf(out_stream," L^2 = %.2f",L2_loc[i*n_s+i]);
//             if(print)if(L2_loc[i*n_s+i]<1e-3){
//                 if     (fabs(P_loc[i*n_s+i]-1.0)<1e-8)printf(" (+)");
//                 else if(fabs(P_loc[i*n_s+i]+1.0)<1e-8)printf(" (-)");
//                 else                                  printf(" (?)");
//             }
//         }
        
        for(int j=0; j<Nd;j++)c_2[j]=coef_r[i_s][j*n_states[i_s]+i]*
                                     coef_r[i_s][j*n_states[i_s]+i]+
                                     coef_i[i_s][j*n_states[i_s]+i]*
                                     coef_i[i_s][j*n_states[i_s]+i];
        if(print){
            fprintf(out_stream,":\n");
            find_max_coef(max_coef, print_number, c_2, Nd, 0,1);
            for(int j=0;j<print_number;j++){
                fprintf(out_stream,"% .10e| % .10e I | ",coef_r[i_s][max_coef[j]*n_states[i_s]+i],coef_i[i_s][max_coef[j]*n_states[i_s]+i]);
                printf_occ(max_coef[j]);
                fprintf(out_stream,"\n");
            }
        }
    }
    
    delete[] max_coef;
//     delete[] S2_loc;
//     if(LINEAR)delete[] L2_loc;
//     if(LINEAR)delete[]  P_loc;
    
    return 0;
}

int sum_index(int i, int j, aldet_data * CI, aldet_rel_data * CI_r){
    
    
    memset(CI_r->bit,0,2*CI->n_act*sizeof(int));
    for (int k=1; k<CI->na+1; k++) CI_r->bit[CI->vec_a[i*(CI->na+1)+k]-1          ] = 1;
    for (int k=1; k<CI->nb+1; k++) CI_r->bit[CI->vec_b[j*(CI->nb+1)+k]-1+CI->n_act] = 1;
    
    
    int ij = get_ind_from_ON(CI_r->bit, 2*CI_r->n_act, CI_r->n_el, CI_r->f, CI_r->buf);
//     CI->printf_occ_a(i);printf("|");CI->printf_occ_b(j);printf("\n");
//     CI_r->printf_occ(ij);printf("\n");
//     getchar();
    return ij;
}

int compress_H_full_to_H_ab(double * H, double * Hf, aldet_data * CI, aldet_rel_data * CI_r){
    
    set_zero_matr(H,CI->Nd*CI->Nd);
    
    for(int i=0; i<CI->Na;i++)
    for(int j=0; j<CI->Nb;j++){
        int ij = sum_index(i,j,CI,CI_r);
        for(int k=0; k<CI->Na;k++)
        for(int l=0; l<CI->Nb;l++){
            int kl = sum_index(k,l,CI,CI_r);
            H[(i*CI->Nb+j)*CI->Nd+k*CI->Nb+l]=Hf[ij*CI_r->Nd+kl];
        }
    }
    
//     PrintMatr(H,CI->Nd,CI->Nd,0);
//     exit(0);
    
    
    return 0;
    
    
}

aldet_rel_data::~aldet_rel_data(){
    
    if(f   != NULL)delete[] f;
    if(vec != NULL)delete[] vec;
    
    if(bit != NULL)delete[] bit;
    if(buf != NULL)delete[] buf;//more space for +a-b and -a+b ????
    
    if(n_states != NULL)delete[] n_states;

    if(coef_r  !=NULL)for(int i=0;i<n_sets;i++)if(coef_r  [i]!=NULL)delete[] coef_r  [i];delete[] coef_r  ;
    if(coef_i  !=NULL)for(int i=0;i<n_sets;i++)if(coef_i  [i]!=NULL)delete[] coef_i  [i];delete[] coef_i  ;
    if(E_states!=NULL)for(int i=0;i<n_sets;i++)if(E_states[i]!=NULL)delete[] E_states[i];delete[] E_states;
    if(S2      !=NULL)for(int i=0;i<n_sets;i++)if(S2      [i]!=NULL)delete[] S2      [i];delete[] S2      ;
    if(L2      !=NULL)for(int i=0;i<n_sets;i++)if(L2      [i]!=NULL)delete[] L2      [i];delete[] L2      ;
    if(P       !=NULL)for(int i=0;i<n_sets;i++)if(P       [i]!=NULL)delete[] P       [i];delete[] P       ;
    if(E_act   !=NULL)for(int i=0;i<n_sets;i++)if(E_act   [i]!=NULL)delete[] E_act   [i];delete[] E_act   ;
    

    
    if(H_diag      !=NULL) delete[] H_diag     ;
    if(H_diag_appr !=NULL) delete[] H_diag_appr;
    
    if(Hv_r      != NULL) delete[] Hv_r      ;
    if(Hv_i      != NULL) delete[] Hv_i      ;
    if(Hv_buf_r  != NULL) delete[] Hv_buf_r  ;
    if(Hv_buf_i  != NULL) delete[] Hv_buf_i  ;
    
//     if(Ia         != NULL) delete[] Ia        ;
//     if(Ia_friends != NULL) delete[] Ia_friends;
//     if(Ib         != NULL) delete[] Ib        ;
//     if(Ib_friends != NULL) delete[] Ib_friends;
//     
    if(J_act  != NULL) delete[] J_act;
    if(K_act  != NULL) delete[] K_act;
    
    if(e1_ind    != NULL) delete[]e1_ind      ;
    if(e1_orbs   != NULL) delete[]e1_orbs     ;
    if(e1_sign   != NULL) delete[]e1_sign     ;

    if(e1_ind_so != NULL) delete[]e1_ind_so     ;
    if(e1_h_so   != NULL) delete[]e1_h_so       ;

    
    if(e2_sign   !=NULL) delete[] e2_sign ;
    if(e2_orbs   !=NULL) delete[] e2_orbs ;
    if(e2_ind    != NULL) delete[] e2_ind;  
    if(e2_V      != NULL) delete[] e2_V  ;  
    
    if(e3_ind    != NULL) delete[] e3_ind;
    if(e3_V      != NULL) delete[] e3_V  ;
    
    if(spin_sign != NULL) delete[] spin_sign;
    
    if(F_act_r    !=NULL) delete[] F_act_r    ;
    if(F_act_i    !=NULL) delete[] F_act_i    ;
//     if(F_act_B    !=NULL) delete[] F_act_B    ;
    if(act_INTS_AA!=NULL) delete[] act_INTS_AA;
    if(act_INTS_AB!=NULL) delete[] act_INTS_AB;
    if(act_INTS_BB!=NULL) delete[] act_INTS_BB;
    if(T3_AAA    !=NULL) delete[] T3_AAA    ;
    if(T3_AAB    !=NULL) delete[] T3_AAB    ;
    if(T3_BBA    !=NULL) delete[] T3_BBA    ;
    if(T3_BBB    !=NULL) delete[] T3_BBB    ;
    if(T2_AA     !=NULL) delete[] T2_AA     ;
    if(T2_AB     !=NULL) delete[] T2_AB     ;
    if(T2_BB     !=NULL) delete[] T2_BB     ;
    if(T1_A      !=NULL) delete[] T1_A      ;
    if(T1_B      !=NULL) delete[] T1_B      ;

    
}
