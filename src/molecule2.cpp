# include <stdio.h>
# include "blas_link.h"
# include "molecule.h"
# include "doCI_data.h"
# include "matr.h"
# include "timer.h"
# include "common_vars.h"
# include "RI.h"
// 
// 
# ifdef LIBINT
# include <vector>
# include <libint2/atom.h>
# include <libint2/engine.h>
using libint2::Atom;
using libint2::Shell;
using libint2::Engine;
using libint2::Operator;
# include "libint_link.h"
#endif


int calc_2el_CI(molecule * M, doCI_data * D){
    
    if(RI==0)
        calc_2el_MO_INTS(M->s, M->n_ao, D->DM_C, D->J, D->K, D->act_INTS, D->R_ACT, D->L_ACT, D->n_act[0]);        
        
    if(RI==1)if(nQM==1)
        calc_2el_MO_INTS_RI(M->n_ao, D->DM_C, D->J, D->K, D->act_INTS, D->VEC[0], D->R_ACT, D->L_ACT, M->n_cor_orb, D->n_act[0],1);
    
    if(RI==1)if(nQM==0){
        double * MO = D->VEC[0] + M->n_cor_orb*M->n_ao;
        calc_2el_MO_INTS_RI(M->n_ao, D->DM_C, D->J, D->K, D->act_INTS, D->VEC[0], MO, MO, M->n_cor_orb, D->n_act[0],0);
    }
    
    
    return 0;
}

int calc_2el_CI_1st_order(molecule * M, doCI_data * D, int nshA, int nshB){
    
    std::vector<Engine> e_engines(num_threads);
//     e_engines.resize(num_threads);
    e_engines[0] = Engine(Operator::coulomb,GTO_MAX,L_MAX);
    for (size_t i = 0; i != num_threads; ++i) {
        e_engines[i] = e_engines[0];
    }
    
    
    double *J[2];
//     J[0]=new double[M->n_ao*M->n_ao];set_zero_matr(J[0],M->n_ao*M->n_ao);
//     J[1]=new double[M->n_ao*M->n_ao];set_zero_matr(J[1],M->n_ao*M->n_ao);
//     
//     double *K[2];
//     K[0]=new double[M->n_ao*M->n_ao];set_zero_matr(K[0],M->n_ao*M->n_ao);
//     K[1]=new double[M->n_ao*M->n_ao];set_zero_matr(K[1],M->n_ao*M->n_ao);
//     
    D->K_R_1st=new double[M->n_ao*M->n_ao];set_zero_matr(D->K_R_1st,M->n_ao*M->n_ao);
    D->K_L_1st=new double[M->n_ao*M->n_ao];set_zero_matr(D->K_L_1st,M->n_ao*M->n_ao);
    
//     printf("R_O:\n");
//     PrintMatr(D->R_ACT,M->n_ao,D->n_act[0],1);
//     printf("L_O_1st:\n");
//     PrintMatr(D->L_ACT_1st_order,M->n_ao,D->n_act[0],1);
    
    std::vector<Shell> sh_A;
    std::vector<Shell> sh_B;
//     sh_A.resize(0);
//     sh_B.resize(0);
    Shell tmp_s;
    printf("number of shells - A=%d, B=%d, %d\n",nshA, nshB,M->s.size());
//     getchar();
    
    for(int i=0     ;i<nshA;i++){
//         std::cerr << M->s[i]<<std::endl;
//         tmp_s = Shell(M->s[i].alpha,{{M->s[i].contr[0].l, false, M->s[i].contr[0].coeff},},{{ M->s[i].O[0], M->s[i].O[1], M->s[i].O[2]}});
        sh_A.push_back(M->s[i]);
        
    }
    for(int i=nshA;i<nshA+nshB;i++)sh_B.push_back(M->s[i]);

    int n_ao_A = 0;
    for(int i=0;i<nshA;i++)n_ao_A+=sh_A[i].size();
    
    int n_ao_B = 0;
    for(int i=0;i<nshB;i++)n_ao_B+=sh_B[i].size();
    
    
    int n_act_A = D->n_act_f[0][0];
    int n_act_B = D->n_act_f[0][1];
    double * L_ACT_A = new double[n_act_A*n_ao_A];
    double * R_ACT_A = new double[n_act_A*n_ao_A];
    double * L_ACT_B = new double[n_act_B*n_ao_B];
    double * R_ACT_B = new double[n_act_B*n_ao_B];
    double * L_ACT_A_1st = new double[n_act_A*n_ao_B];
    double * R_ACT_A_1st = new double[n_act_A*n_ao_B];
    double * L_ACT_B_1st = new double[n_act_B*n_ao_A];
    double * R_ACT_B_1st = new double[n_act_B*n_ao_A];
    
    for(int j=0;j<n_ao_A ;j++)
    for(int i=0;i<n_act_A;i++)
        L_ACT_A[j*n_act_A+i]=D->L_ACT[j*D->n_act[0]+i];
    
//     printf("L_A:\n");
//     PrintMatr(L_ACT_A,n_ao_A,n_act_A,1);
    
    for(int j=0;j<n_ao_B ;j++)
    for(int i=0;i<n_act_B;i++)
        L_ACT_B[j*n_act_B+i]=D->L_ACT[(j+n_ao_A)*D->n_act[0]+i+n_act_A];
    
//     printf("L_B:\n");
//     PrintMatr(L_ACT_B,n_ao_A,n_act_A,1);
    
    for(int j=0;j<n_ao_A ;j++)
    for(int i=0;i<n_act_A;i++)
        R_ACT_A[j*n_act_A+i]=D->R_ACT[j*D->n_act[0]+i];
    
//     printf("R_A:\n");
//     PrintMatr(R_ACT_A,n_ao_A,n_act_A,1);
    
    for(int j=0;j<n_ao_B ;j++)
    for(int i=0;i<n_act_B;i++)
        R_ACT_B[j*n_act_B+i]=D->R_ACT[(j+n_ao_A)*D->n_act[0]+i+n_act_A];
    
//     printf("R_B:\n");
//     PrintMatr(R_ACT_B,n_ao_A,n_act_A,1);
    
    for(int j=0;j<n_ao_B ;j++)
    for(int i=0;i<n_act_A;i++)
        L_ACT_A_1st[j*n_act_A+i]=D->L_ACT_1st_order[(j+n_ao_A)*D->n_act[0]+i];
    
//     printf("L_A:\n");
//     PrintMatr(L_ACT_A_1st,n_ao_A,n_act_B,1);
    
    for(int j=0;j<n_ao_A ;j++)
    for(int i=0;i<n_act_B;i++)
        L_ACT_B_1st[j*n_act_B+i]=D->L_ACT_1st_order[j*D->n_act[0]+i+n_act_A];
    
//     printf("L_B:\n");
//     PrintMatr(L_ACT_B_1st,n_ao_B,n_act_A,1);
    
    for(int j=0;j<n_ao_B ;j++)
    for(int i=0;i<n_act_A;i++)
        R_ACT_A_1st[j*n_act_A+i]=D->R_ACT_1st_order[(j+n_ao_A)*D->n_act[0]+i];
    
//     printf("R_A:\n");
//     PrintMatr(R_ACT_A_1st,n_ao_A,n_act_B,1);
    
    for(int j=0;j<n_ao_A ;j++)
    for(int i=0;i<n_act_B;i++)
        R_ACT_B_1st[j*n_act_B+i]=D->R_ACT_1st_order[j*D->n_act[0]+i+n_act_A];
    
//     printf("R_B:\n");
//     PrintMatr(R_ACT_B_1st,n_ao_B,n_act_A,1);
    
    if(n_act_A!=n_act_B){
        printf("realized only for n_act_A==n_act_B");
        exit(0);
    }
    
    double *RO_INTS = new double[D->n_act[0]*D->n_act[0]*D->n_act[0]*D->n_act[0]];
    double *LO_INTS = new double[D->n_act[0]*D->n_act[0]*D->n_act[0]*D->n_act[0]];
    
    double *RO_INTS_J = new double[n_act_A*n_act_A*n_act_A*n_act_A];
    double *LO_INTS_J = new double[n_act_A*n_act_A*n_act_A*n_act_A];
    double *RO_INTS_K = new double[n_act_A*n_act_A*n_act_A*n_act_A];
    double *LO_INTS_K = new double[n_act_A*n_act_A*n_act_A*n_act_A];
    
    
    
    double ** DM_K = new double*[3];
    DM_K[0]=D->DM_C;
    DM_K[1]=D->DM_L_1st;
    DM_K[2]=D->DM_R_1st;
    
    double ** K_calc = new double*[3];
    K_calc[0]=D->K;
    K_calc[1]=D->K_L_1st;
    K_calc[2]=D->K_R_1st;
    
    set_zero_matr(D->J,M->n_ao*M->n_ao);
    set_zero_matr(D->K,M->n_ao*M->n_ao);
    
//     calc_2el_MO_INTS(M->s, M->n_ao, D->DM_C, D->J, D->K, D->act_INTS, D->R_ACT, D->L_ACT, D->n_act[0],&(e_engines));

    calc_2el_MO_INTS_4sets(sh_A, n_ao_A, NULL,NULL,0, &(D->DM_R_1st),&(D->K_R_1st),0, RO_INTS_J, R_ACT_B_1st, R_ACT_A, L_ACT_A, L_ACT_A, n_act_A,&(e_engines));

    calc_2el_MO_INTS_4sets(sh_A, n_ao_A, NULL,NULL,0, &(D->DM_R_1st),&(D->K_R_1st),0, LO_INTS_J, R_ACT_A, R_ACT_A, L_ACT_B_1st, L_ACT_A, n_act_A,&(e_engines));

    calc_2el_MO_INTS_4sets_2shells(sh_B, sh_A, n_ao_B, n_ao_A, NULL,NULL,0, &(D->DM_R_1st),&(D->K_R_1st),0, LO_INTS_K, R_ACT_B, R_ACT_A, L_ACT_A_1st, L_ACT_A, n_act_A,&(e_engines));

    calc_2el_MO_INTS_4sets_2shells(sh_B, sh_A, n_ao_B, n_ao_A, NULL,NULL,0, &(D->DM_R_1st),&(D->K_R_1st),0, RO_INTS_K, R_ACT_A_1st, R_ACT_A, L_ACT_B, L_ACT_A, n_act_A,&(e_engines));

    calc_2el_MO_INTS_multi_JK(M->s, M->n_ao, &(D->DM_C), &(D->J), 1, DM_K, K_calc, 3, D->act_INTS, D->R_ACT, D->L_ACT, D->n_act[0],&(e_engines));
    
    
    
//     DM_to_F_transform(&(D->J), K_calc,
//                       &(D->DM_C), 1,
//                       DM_K, 3, 
//                       M->s, M->n_ao, 
//                       &(e_engines),num_threads);
    
//     DM_to_F_transform(NULL, &(D->K_R_1st),
//                       NULL, 0,
//                       &(D->DM_R_1st), 1, 
//                       M->s, M->n_ao, 
//                       &(e_engines),num_threads);
    
//     DM_to_F_transform(NULL, &(D->K_L_1st),
//                       NULL, 0,
//                       &(D->DM_L_1st), 1, 
//                       M->s, M->n_ao, 
//                       &(e_engines),num_threads);
    
    
//     calc_2el_MO_INTS_4sets(M->s, M->n_ao, D->DM_C_F,J,0, &(D->DM_R_1st),&(D->K_R_1st),1, RO_INTS, D->R_ACT_1st_order, D->R_ACT, D->L_ACT, D->L_ACT, 0,&(e_engines)); //RO - (i j^(1) |kl)
    
//     calc_2el_MO_INTS_4sets(M->s, M->n_ao, D->DM_C_F+1,J+1,0, &(D->DM_L_1st),&(D->K_L_1st),1, LO_INTS, D->R_ACT, D->R_ACT, D->L_ACT_1st_order, D->L_ACT, 0,&(e_engines)); //RO - (i^(1) j|kl)
    
    
    D->act_INTS_1st_order = new double[D->n_act[0]*D->n_act[0]*D->n_act[0]*D->n_act[0]];
    set_zero_matr(D->act_INTS_1st_order,D->n_act[0]*D->n_act[0]*D->n_act[0]*D->n_act[0]);
    
//     printf("LO_INTS:\n");
//     PrintMatr(LO_INTS,D->n_act[0]*D->n_act[0],D->n_act[0]*D->n_act[0],1);
//     printf("LO_INTS_K:\n");
//     PrintMatr(LO_INTS_K,n_act_B*n_act_B,n_act_B*n_act_B,1);
//     exit(0);
    
    int num;
    int num2;
    int num3;
    
    int ** n_act_f = D->n_act_f;
    int * n_act = D->n_act;
    
    for(int i=0;i<D->n_act[0]*D->n_act[0]*D->n_act[0]*D->n_act[0];i++)
        D->act_INTS_1st_order[i]=D->act_INTS[i];
    
//     printf("act_INTS:\n");
//     PrintMatr(D->act_INTS,n_act[0]*n_act[0],n_act[0]*n_act[0],1); 
//     printf("act_INTS_1st_order:\n");
//     PrintMatr(D->act_INTS_1st_order,n_act[0]*n_act[0],n_act[0]*n_act[0],1);    
//     exit(0);
    
    for(int i=0            ;i<D->n_act_f[0][0]              ;i++)
    for(int j=0            ;j<D->n_act_f[0][0]              ;j++)
    for(int k=0            ;k<D->n_act_f[0][0]              ;k++)
    for(int l=n_act_f[0][0];l<D->n_act_f[0][0]+n_act_f[0][1];l++){
        num =((i*n_act[0]+j)*n_act[0]+k)*n_act[0]+l;
        num2=((k*n_act[0]+l)*n_act[0]+i)*n_act[0]+j;
        num3=((k*n_act_A+l-n_act_A)*n_act_A+i)*n_act_A+j;
        //AAAB working!!!! in ab|a+b-
        //(ij|kl)^(0+1)=(ij|kl)^(0)+(i j|k^(1) l)+(i j|k l^(1))
        D->act_INTS_1st_order[num]+=RO_INTS_J[num3]+LO_INTS_K[num3];
    }
    
    for(int i=0            ;i<D->n_act_f[0][0]              ;i++)
    for(int j=0            ;j<D->n_act_f[0][0]              ;j++)
    for(int k=n_act_f[0][0];k<D->n_act_f[0][0]+n_act_f[0][1];k++)
    for(int l=0            ;l<D->n_act_f[0][0]              ;l++){
        num =((i*n_act[0]+j)*n_act[0]+k)*n_act[0]+l;
        num2=((k*n_act[0]+l)*n_act[0]+i)*n_act[0]+j;
        num3=(((k-n_act_A)*n_act_A+l)*n_act_A+i)*n_act_A+j;
        //AABA working in a+b-|ab
        //(ij|kl)^(0+1)=(ij|kl)^(0)+(i j|k^(1) l)+(i j|k l^(1))
        D->act_INTS_1st_order[num]+=RO_INTS_K[num3]+LO_INTS_J[num3];
    }
    
    for(int i=0            ;i<D->n_act_f[0][0]              ;i++)
    for(int j=n_act_f[0][0];j<D->n_act_f[0][0]+n_act_f[0][1];j++)
    for(int k=0            ;k<D->n_act_f[0][0]              ;k++)
    for(int l=0            ;l<D->n_act_f[0][0]              ;l++){
        num=((i*n_act[0]+j)*n_act[0]+k)*n_act[0]+l;
        num3=((i*n_act_A+j-n_act_A)*n_act_A+k)*n_act_A+l;
        //ABAA
        //(ij|kl)^(0+1)=(ij|kl)^(0)+(i^(1) j|k l)+(i j^(1)|k l)
        D->act_INTS_1st_order[num]+=RO_INTS_J[num3]+LO_INTS_K[num3];
    }
    
    for(int i=n_act_f[0][0];i<D->n_act_f[0][0]+n_act_f[0][1];i++)
    for(int j=0            ;j<D->n_act_f[0][0]              ;j++)
    for(int k=0            ;k<D->n_act_f[0][0]              ;k++)
    for(int l=0            ;l<D->n_act_f[0][0]              ;l++){
        int num=((i*n_act[0]+j)*n_act[0]+k)*n_act[0]+l;
        num3=(((i-n_act_A)*n_act_A+j)*n_act_A+k)*n_act_A+l;
        //BAAA
        //(ij|kl)^(0+1)=(ij|kl)^(0)+(i^(1) j|k l)+(i j^(1)|k l)
        D->act_INTS_1st_order[num]+=RO_INTS_K[num3]+LO_INTS_J[num3];
    }

    calc_2el_MO_INTS_4sets(sh_B, n_ao_B, NULL,NULL,0, &(D->DM_R_1st),&(D->K_R_1st),0, RO_INTS_J, R_ACT_A_1st, R_ACT_B, L_ACT_B, L_ACT_B, n_act_B,&(e_engines));

    calc_2el_MO_INTS_4sets(sh_B, n_ao_B, NULL,NULL,0, &(D->DM_R_1st),&(D->K_R_1st),0, LO_INTS_J, R_ACT_B, R_ACT_B, L_ACT_A_1st, L_ACT_B, n_act_B,&(e_engines));
    
    calc_2el_MO_INTS_4sets_2shells(sh_A, sh_B, n_ao_A, n_ao_B, NULL,NULL,0, &(D->DM_R_1st),&(D->K_R_1st),0, LO_INTS_K, R_ACT_A, R_ACT_B, L_ACT_B_1st, L_ACT_B, n_act_A,&(e_engines));

    calc_2el_MO_INTS_4sets_2shells(sh_A, sh_B, n_ao_A, n_ao_B, NULL,NULL,0, &(D->DM_R_1st),&(D->K_R_1st),0, RO_INTS_K, R_ACT_B_1st, R_ACT_B, L_ACT_A, L_ACT_B, n_act_A,&(e_engines));
    
    
    
    
    //3B-1A
    for(int i=n_act_f[0][0];i<D->n_act_f[0][0]+n_act_f[0][1];i++)
    for(int j=n_act_f[0][0];j<D->n_act_f[0][0]+n_act_f[0][1];j++)
    for(int k=n_act_f[0][0];k<D->n_act_f[0][0]+n_act_f[0][1];k++)
    for(int l=0            ;l<D->n_act_f[0][0]              ;l++){
        num =((i*n_act[0]+j)*n_act[0]+k)*n_act[0]+l;
//         num2=((k*n_act[0]+l)*n_act[0]+i)*n_act[0]+j;
        num3=(((k-n_act_A)*n_act_B+l)*n_act_B+i-n_act_A)*n_act_B+j-n_act_A;
        //BBBA
        //(ij|kl)^(0+1)=(ij|kl)^(0)+(i j|k^(1) l)+(i j|k l^(1))
        D->act_INTS_1st_order[num]+=RO_INTS_J[num3]+LO_INTS_K[num3];
    }
    
    for(int i=n_act_f[0][0];i<D->n_act_f[0][0]+n_act_f[0][1];i++)
    for(int j=n_act_f[0][0];j<D->n_act_f[0][0]+n_act_f[0][1];j++)
    for(int k=0            ;k<D->n_act_f[0][0]              ;k++)
    for(int l=n_act_f[0][0];l<D->n_act_f[0][0]+n_act_f[0][1];l++){
        num =((i*n_act[0]+j)*n_act[0]+k)*n_act[0]+l;
//         num2=((k*n_act[0]+l)*n_act[0]+i)*n_act[0]+j;
        num3=((k*n_act_B+l-n_act_A)*n_act_B+i-n_act_A)*n_act_B+j-n_act_A;
        //BBAB
        //(ij|kl)^(0+1)=(ij|kl)^(0)+(i j|k^(1) l)+(i j|k l^(1))
        D->act_INTS_1st_order[num]+=RO_INTS_K[num3]+LO_INTS_J[num3];
    }
    
    for(int i=n_act_f[0][0];i<D->n_act_f[0][0]+n_act_f[0][1];i++)
    for(int j=0            ;j<D->n_act_f[0][0]              ;j++)
    for(int k=n_act_f[0][0];k<D->n_act_f[0][0]+n_act_f[0][1];k++)
    for(int l=n_act_f[0][0];l<D->n_act_f[0][0]+n_act_f[0][1];l++){
        num=((i*n_act[0]+j)*n_act[0]+k)*n_act[0]+l;
        num3=(((i-n_act_A)*n_act_B+j)*n_act_B+k-n_act_A)*n_act_B+l-n_act_A;
        //BABB
        //(ij|kl)^(0+1)=(ij|kl)^(0)+(i^(1) j|k l)+(i j^(1)|k l)
        D->act_INTS_1st_order[num]+=RO_INTS_J[num3]+LO_INTS_K[num3];
    }
    
    for(int i=0            ;i<D->n_act_f[0][0]              ;i++)
    for(int j=n_act_f[0][0];j<D->n_act_f[0][0]+n_act_f[0][1];j++)
    for(int k=n_act_f[0][0];k<D->n_act_f[0][0]+n_act_f[0][1];k++)
    for(int l=n_act_f[0][0];l<D->n_act_f[0][0]+n_act_f[0][1];l++){
        num=((i*n_act[0]+j)*n_act[0]+k)*n_act[0]+l;
        num3=((i*n_act_B+j-n_act_A)*n_act_B+k-n_act_A)*n_act_B+l-n_act_A;
        //ABBB
        //(ij|kl)^(0+1)=(ij|kl)^(0)+(i^(1) j|k l)+(i j^(1)|k l)
        D->act_INTS_1st_order[num]+=RO_INTS_K[num3]+LO_INTS_J[num3];
    }
    
    
//     for(int i=0; i<M->n_ao*M->n_ao;i++)
//         D->J[i]=J[0][i]+J[1][i];
//     for(int i=0; i<M->n_ao*M->n_ao;i++)
//         D->K[i]=K[0][i]+K[1][i];
//     
    
    
    
    
    return 0;
}
