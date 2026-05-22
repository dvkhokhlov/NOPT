#include <stdio.h>
#include <cstring>
// #include <stdlib.h>
#include "matr.h"
#include "trcamm.h"
// # include "libint_link.h"

extern int num_threads;

int gen_atomic_DM(double * A_DM, double * DM, int atom, std::vector<Shell> * s, std::vector<int> * shell_center, int n_ao){
    
//     for(int i=0; i<n_ao; i++){
//         printf("C[%d] = %d\n", i, ao_center[i]);
//         getchar();
//     }
    
    set_zero_matr(A_DM,n_ao*n_ao);
    
    int ao_num=0;
    for(int sh_i=0;sh_i<s->size();sh_i++){
        int sh_s=s[0][sh_i].contr[0].size();
        if(shell_center[0][sh_i]==atom){
            for(int i=0;i<sh_s;i++)
                memcpy(A_DM+(ao_num+i)*n_ao,DM+(ao_num+i)*n_ao,n_ao*sizeof(double));
        
        }
        ao_num+=s[0][sh_i].contr[0].size();
    }
    
    symmetrization_with_scaling(A_DM, n_ao,1.0);
//     PrintMatr(A_DM,n_ao,n_ao,1);
//     exit(1);
    
    
    return 0;
}

int gen_atomic_DM_anti(double * A_DM, double * DM, int atom, int * ao_center, int n_ao){
    
    
    set_zero_matr(A_DM,n_ao*n_ao);
    
    for(int i=0; i<n_ao; i++)if(ao_center[i]==atom)
        memcpy(A_DM+i*n_ao,DM+i*n_ao,n_ao*sizeof(double));
            
    anti_symmetrization_with_scaling(A_DM, n_ao,1.0);
   
    
    return 0;
}

// int transform_to_libint_mat(double * M, AO_basis * b, int n_ao){
//     
//     
//     int k;
//     double c;
//     for(int j=0; j<n_ao;j++){
//         k=0;
//         while(k<b[0].n_ao){
//             if(b[0].l_ao[k]==0){k++;continue;}
//             if(b[0].l_ao[k]==1){k+=3;continue;}
//             if(b[0].l_ao[k]==2){
//                 // 0=0
//                 // 1<->3
//                 // 2<-4<-5<-2
//                 c            =M[j*n_ao+k+1];
//                 M[j*n_ao+k+1]=M[j*n_ao+k+3]/sqrt(3.0);
//                 M[j*n_ao+k+3]=c;
//                 c            =M[j*n_ao+k+2];
//                 M[j*n_ao+k+2]=M[j*n_ao+k+4]/sqrt(3.0);
//                 M[j*n_ao+k+4]=M[j*n_ao+k+5]/sqrt(3.0);
//                 M[j*n_ao+k+5]=c;
// 
//                 k+=6;
//                 continue;
//             }
//             if(b[0].l_ao[k]==3){
//                 // 0=0
//                 // 1<-3<-5<-7<-6<-1
//                 // 2<-4<-9<-2
//                 // 8=8
//                 c            =M[j*n_ao+k+1]          ;
//                 M[j*n_ao+k+1]=M[j*n_ao+k+3]/sqrt( 5.0);
//                 M[j*n_ao+k+3]=M[j*n_ao+k+5]/sqrt( 5.0);
//                 M[j*n_ao+k+5]=M[j*n_ao+k+7]/sqrt( 5.0);
//                 M[j*n_ao+k+7]=M[j*n_ao+k+6]/sqrt( 5.0);
//                 M[j*n_ao+k+6]=c;
//                 c            =M[j*n_ao+k+2];
//                 M[j*n_ao+k+2]=M[j*n_ao+k+4]/sqrt( 5.0);
//                 M[j*n_ao+k+4]=M[j*n_ao+k+9]/sqrt(15.0);
//                 M[j*n_ao+k+9]=c;
//                 M[j*n_ao+k+8]=M[j*n_ao+k+8]/sqrt( 5.0);
//                 k+=10;
//                 continue;
//             }
//         }
//     }
//     
//     for(int j=0; j<n_ao;j++){
//         k=0;
//         while(k<b[0].n_ao){
//             if(b[0].l_ao[k]==0){k++;continue;}
//             if(b[0].l_ao[k]==1){k+=3;continue;}
//             if(b[0].l_ao[k]==2){
//                 // 0=0
//                 // 1<->3
//                 // 2<-4<-5<-2
//                 c              =M[j+n_ao*(k+1)];
//                 M[j+n_ao*(k+1)]=M[j+n_ao*(k+3)]/sqrt(3.0);
//                 M[j+n_ao*(k+3)]=c;
//                 c              =M[j+n_ao*(k+2)];
//                 M[j+n_ao*(k+2)]=M[j+n_ao*(k+4)]/sqrt(3.0);
//                 M[j+n_ao*(k+4)]=M[j+n_ao*(k+5)]/sqrt(3.0);
//                 M[j+n_ao*(k+5)]=c;
// 
//                 k+=6;
//                 continue;
//             }
//             if(b[0].l_ao[k]==3){
//                 // 0=0
//                 // 1<-3<-5<-7<-6<-1
//                 // 2<-4<-9<-2
//                 // 8=8
//                 c              =M[j+n_ao*(k+1)]          ;
//                 M[j+n_ao*(k+1)]=M[j+n_ao*(k+3)]/sqrt( 5.0);
//                 M[j+n_ao*(k+3)]=M[j+n_ao*(k+5)]/sqrt( 5.0);
//                 M[j+n_ao*(k+5)]=M[j+n_ao*(k+7)]/sqrt( 5.0);
//                 M[j+n_ao*(k+7)]=M[j+n_ao*(k+6)]/sqrt( 5.0);
//                 M[j+n_ao*(k+6)]=c;
//                 c              =M[j+n_ao*(k+2)];
//                 M[j+n_ao*(k+2)]=M[j+n_ao*(k+4)]/sqrt( 5.0);
//                 M[j+n_ao*(k+4)]=M[j+n_ao*(k+9)]/sqrt(15.0);
//                 M[j+n_ao*(k+9)]=c;
//                 M[j+n_ao*(k+8)]=M[j+n_ao*(k+8)]/sqrt( 5.0);
//                 k+=10;
//                 continue;
//             }
//         }
//     }
//     
//     
//     
//     return 1;
// }



// int calc_TrCAMM_by_DM(molecule * A, double ** DM, double * S_WF, int n_s){
// 
//     libint2::initialize();
// 
//     
//     int n_ao = A->basis[0].n_ao;
//     printf("n_ao = %d\n", n_ao);
//     
//     molecule_1el_data cD;
//     cD.gen_start_data(A, 0, NULL, num_threads, NULL);
//     cD.alloc_Q_AO();
// 
//     cD.calc_S_AO();
//     cD.calc_D_AO();
//     cD.calc_Q_AO();
//     
//     cD.P_calc();
//     cD.RxP_calc();
//     
// //     PrintMatr(cD.P_AO, n_ao, n_ao,0);
//     
// //     exit(1);
// 
//     
//     printf("# TrCAMM\n");
//     printf("%d %d\n", n_s, A->n_atoms);
// 
//     double S=0;
//     double * DM_ij;
//     double * A_DM;
//     A_DM = new double[n_ao*n_ao];
//     double x,y,z;
//     double c, d_x, d_y, d_z;
//     double q_xx, q_xy, q_xz, q_yy, q_yz, q_zz;
//     
//     double * gr_c    =new double[A->n_atoms];
//     double * gr_d_x  =new double[A->n_atoms];
//     double * gr_d_y  =new double[A->n_atoms];
//     double * gr_d_z  =new double[A->n_atoms];
//     double * gr_q_xx =new double[A->n_atoms];
//     double * gr_q_xy =new double[A->n_atoms];
//     double * gr_q_xz =new double[A->n_atoms];
//     double * gr_q_yy =new double[A->n_atoms];
//     double * gr_q_yz =new double[A->n_atoms];
//     double * gr_q_zz =new double[A->n_atoms];
//     
//     printf("\n\nZero diagonal element is written as it is, other diagonal elements are given as a difference with zero-element\n\n");
//     printf("\n\nbug with transposed DM is fixed!!!\n\n");
//     
//     
//     for(int i_s=0;i_s<n_s;i_s++)
//     for(int j_s=0;j_s<n_s;j_s++){
//         
//         printf("[%d %d]\n",i_s,j_s);
//         DM_ij = DM[i_s*n_s+j_s];
// //         if(i_s==j_s)S=1;
// //         else        S=0;
//         S=S_WF[i_s*n_s+j_s];
//         
//         
//         for(int i_a=0;i_a<A->n_atoms;i_a++){
//             x=A->atom_coord[3*i_a+0];
//             y=A->atom_coord[3*i_a+1];
//             z=A->atom_coord[3*i_a+2];
//             
//             //by the 273 standard TrCAMM is taken in the electron population way, but not in the charge way 
//             c   = calc_atom_multipole(A_DM, DM_ij, cD.S_MO  , i_a, A->basis[0].ao_center, n_ao,0);//S_MO! there is no error
//             d_x = calc_atom_multipole(A_DM, DM_ij, cD.Dx_AO , i_a, A->basis[0].ao_center, n_ao,0);
//             d_y = calc_atom_multipole(A_DM, DM_ij, cD.Dy_AO , i_a, A->basis[0].ao_center, n_ao,0);
//             d_z = calc_atom_multipole(A_DM, DM_ij, cD.Dz_AO , i_a, A->basis[0].ao_center, n_ao,0);
//             q_xx= calc_atom_multipole(A_DM, DM_ij, cD.Qxx_AO, i_a, A->basis[0].ao_center, n_ao,0);
//             q_xy= calc_atom_multipole(A_DM, DM_ij, cD.Qxy_AO, i_a, A->basis[0].ao_center, n_ao,0);
//             q_xz= calc_atom_multipole(A_DM, DM_ij, cD.Qxz_AO, i_a, A->basis[0].ao_center, n_ao,0);
//             q_yy= calc_atom_multipole(A_DM, DM_ij, cD.Qyy_AO, i_a, A->basis[0].ao_center, n_ao,0);
//             q_yz= calc_atom_multipole(A_DM, DM_ij, cD.Qyz_AO, i_a, A->basis[0].ao_center, n_ao,0);
//             q_zz= calc_atom_multipole(A_DM, DM_ij, cD.Qzz_AO, i_a, A->basis[0].ao_center, n_ao,0);
//             
//             q_xx=q_xx-d_x*x-x*d_x+c*x*x;
//             q_xy=q_xy-d_x*y-x*d_y+c*x*y;
//             q_xz=q_xz-d_x*z-x*d_z+c*x*z;
//             q_yy=q_yy-d_y*y-y*d_y+c*y*y;
//             q_yz=q_yz-d_y*z-y*d_z+c*y*z;
//             q_zz=q_zz-d_z*z-z*d_z+c*z*z;
//             
//             d_x=d_x-c*x;
//             d_y=d_y-c*y;
//             d_z=d_z-c*z;
//             
//             c=c-S*A->nucl_charges[i_a];
//             
//             if((i_s*n_s+j_s)==0){
//                 gr_c   [i_a] = c   ;
//                 gr_d_x [i_a] = d_x ;
//                 gr_d_y [i_a] = d_y ;
//                 gr_d_z [i_a] = d_z ;
//                 gr_q_xx[i_a] = q_xx;
//                 gr_q_xy[i_a] = q_xy;
//                 gr_q_xz[i_a] = q_xz;
//                 gr_q_yy[i_a] = q_yy;
//                 gr_q_yz[i_a] = q_yz;
//                 gr_q_zz[i_a] = q_zz;
//             }
//             if(i_s==j_s)if(i_s!=0){
//                 c   -= gr_c   [i_a];
//                 d_x -= gr_d_x [i_a];
//                 d_y -= gr_d_y [i_a];
//                 d_z -= gr_d_z [i_a];
//                 q_xx-= gr_q_xx[i_a];
//                 q_xy-= gr_q_xy[i_a];
//                 q_xz-= gr_q_xz[i_a];
//                 q_yy-= gr_q_yy[i_a];
//                 q_yz-= gr_q_yz[i_a];
//                 q_zz-= gr_q_zz[i_a];
//             }
//                 
//             
//             
//             printf("%.10e %.10e %.10e  %.10e  %.10e %.10e %.10e  %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e \n",
//                    A->atom_coord[i_a*3  ],
//                    A->atom_coord[i_a*3+1],
//                    A->atom_coord[i_a*3+2],
//                    c,
//                    d_x, d_y, d_z,
//                    q_xx, q_xy, q_xz, q_xy, q_yy, q_yz, q_xz, q_yz, q_zz);
//             
//         }
//     }
//     
//     double * Px_AO = new double[n_ao*n_ao];
//     double * Py_AO = new double[n_ao*n_ao];
//     double * Pz_AO = new double[n_ao*n_ao];
//     for(int i2=0;i2<n_ao*n_ao;i2++){
//         Px_AO[i2]=cD.P_AO[3*i2+0];
//         Py_AO[i2]=cD.P_AO[3*i2+1];
//         Pz_AO[i2]=cD.P_AO[3*i2+2];
//     }
//     
//     double * MUx_AO = new double[n_ao*n_ao];
//     double * MUy_AO = new double[n_ao*n_ao];
//     double * MUz_AO = new double[n_ao*n_ao];
//     for(int i2=0;i2<n_ao*n_ao;i2++){
//         MUx_AO[i2]=cD.RxP_AO[9*i2+3*1+2]-cD.RxP_AO[9*i2+3*2+1];
//         MUy_AO[i2]=cD.RxP_AO[9*i2+3*2+0]-cD.RxP_AO[9*i2+3*0+2];
//         MUz_AO[i2]=cD.RxP_AO[9*i2+3*0+1]-cD.RxP_AO[9*i2+3*1+0];
//     }
//     
//     
//     
//     transform_to_libint_mat(Px_AO,A->basis, n_ao);
//     transform_to_libint_mat(Py_AO,A->basis, n_ao);
//     transform_to_libint_mat(Pz_AO,A->basis, n_ao);
//     transform_to_libint_mat(MUx_AO,A->basis, n_ao);
//     transform_to_libint_mat(MUy_AO,A->basis, n_ao);
//     transform_to_libint_mat(MUz_AO,A->basis, n_ao);
//     
//     
//     double p_x, p_y, p_z;
// //     double mu_x, mu_y, mu_z;
//     
// //     double * gr_p_x  =new double[A->n_atoms];
// //     double * gr_p_y  =new double[A->n_atoms];
// //     double * gr_p_z  =new double[A->n_atoms];
//     double * mu_x =new double[A->n_atoms];
//     double * mu_y =new double[A->n_atoms];
//     double * mu_z =new double[A->n_atoms];
//     
//     
//     printf("# Magnetic TrCAMM\n");
//     for(int i_s=0;i_s<n_s;i_s++)
//     for(int j_s=0;j_s<n_s;j_s++){
//         
//         printf("[%d %d]\n",i_s,j_s);
//         DM_ij = DM[i_s*n_s+j_s];
//         
//         printf("# P_i\n");
//         for(int i_a=0;i_a<A->n_atoms;i_a++){
//             x=A->atom_coord[3*i_a+0];
//             y=A->atom_coord[3*i_a+1];
//             z=A->atom_coord[3*i_a+2];
//             // multiply on -1 !!!!!
//             p_x =  - calc_atom_multipole(A_DM, DM_ij, Px_AO , i_a, A->basis[0].ao_center, n_ao,1);
//             p_y =  - calc_atom_multipole(A_DM, DM_ij, Py_AO , i_a, A->basis[0].ao_center, n_ao,1);
//             p_z =  - calc_atom_multipole(A_DM, DM_ij, Pz_AO , i_a, A->basis[0].ao_center, n_ao,1);
//             // multiply on -1 !!!!
//             mu_x[i_a]= - calc_atom_multipole(A_DM, DM_ij, MUx_AO , i_a, A->basis[0].ao_center, n_ao,1);
//             mu_y[i_a]= - calc_atom_multipole(A_DM, DM_ij, MUy_AO , i_a, A->basis[0].ao_center, n_ao,1);
//             mu_z[i_a]= - calc_atom_multipole(A_DM, DM_ij, MUz_AO , i_a, A->basis[0].ao_center, n_ao,1);
//             
//             mu_x[i_a]+=-y*p_z+z*p_y;
//             mu_y[i_a]+=-z*p_x+x*p_z;
//             mu_z[i_a]+=-x*p_y+y*p_x;
//             
// //             if((i_s*n_s+j_s)==0){
// //                 gr_p_x [i_a] = p_x ;
// //                 gr_p_y [i_a] = p_y ;
// //                 gr_p_z [i_a] = p_z ;
// //             }
// //             if(i_s==j_s){
// //                 p_x -= gr_p_x [i_a];
// //                 p_y -= gr_p_y [i_a];
// //                 p_z -= gr_p_z [i_a];
// //             }
//                 
//             
//             
//             printf("%.10e %.10e %.10e \n",p_x, p_y, p_z);
//             
//         }
//         printf("# Mu_i\n");
//         for(int i_a=0;i_a<A->n_atoms;i_a++){
//             printf("%.10e %.10e %.10e \n",mu_x[i_a], mu_y[i_a], mu_z[i_a]);
//         }
//         
//     }
//     
//     
//     
//     
//     cD.clear_Q_AO();
//     
//     
//     
//     delete[] A_DM;
//     delete[]Px_AO;
//     delete[]Py_AO;
//     delete[]Pz_AO;
//     delete[]MUx_AO;
//     delete[]MUy_AO;
//     delete[]MUz_AO;
//     
//     
//     delete[] gr_c   ; 
//     delete[] gr_d_x ; 
//     delete[] gr_d_y ; 
//     delete[] gr_d_z ; 
//     delete[] gr_q_xx; 
//     delete[] gr_q_xy; 
//     delete[] gr_q_xz; 
//     delete[] gr_q_yy; 
//     delete[] gr_q_yz; 
//     delete[] gr_q_zz; 
//     delete[] mu_x;
//     delete[] mu_y;
//     delete[] mu_z;
//     
//     
//     return 0;
// }
// 
// int calc_TrCAMM_by_DM_2(molecule * A, double ** DM, double * S_WF, int n_s1, int n_s2){
// 
//     libint2::initialize();
// 
//     
//     int n_ao = A->basis[0].n_ao;
//     printf("n_ao = %d\n", n_ao);
//     
//     molecule_1el_data cD;
//     cD.gen_start_data(A, 0, NULL, num_threads, NULL);
//     cD.alloc_Q_AO();
// 
//     cD.calc_S_AO();
//     cD.calc_D_AO();
//     cD.calc_Q_AO();
//     
//     cD.P_calc();
//     cD.RxP_calc();
//     
// //     PrintMatr(cD.P_AO, n_ao, n_ao,0);
//     
// //     exit(1);
// 
//     
//     printf("# TrCAMM for nondiagonal part %dx%d\n",n_s1,n_s2);
//     printf("%d %d\n", n_s1, A->n_atoms);
// 
//     double S=0;
//     double * DM_ij;
//     double * A_DM;
//     A_DM = new double[n_ao*n_ao];
//     double x,y,z;
//     double c, d_x, d_y, d_z;
//     double q_xx, q_xy, q_xz, q_yy, q_yz, q_zz;
//     
//     double * gr_c    =new double[A->n_atoms];
//     double * gr_d_x  =new double[A->n_atoms];
//     double * gr_d_y  =new double[A->n_atoms];
//     double * gr_d_z  =new double[A->n_atoms];
//     double * gr_q_xx =new double[A->n_atoms];
//     double * gr_q_xy =new double[A->n_atoms];
//     double * gr_q_xz =new double[A->n_atoms];
//     double * gr_q_yy =new double[A->n_atoms];
//     double * gr_q_yz =new double[A->n_atoms];
//     double * gr_q_zz =new double[A->n_atoms];
//     
//     printf("\n\nZero diagonal element is written as it is, other diagonal elements are given as a difference with zero-element\n\n");
//     printf("\n\nbug with transposed DM is fixed!!!\n\n");
//     
//     
//     for(int i_s=0;i_s<n_s1;i_s++)
//     for(int j_s=0;j_s<n_s2;j_s++){
//         
//         printf("[%d %d]\n",i_s,j_s);
//         DM_ij = DM[i_s*n_s2+j_s];
// //         if(i_s==j_s)S=1;
// //         else        S=0;
//         S=S_WF[i_s*n_s2+j_s];
//         
//         
//         for(int i_a=0;i_a<A->n_atoms;i_a++){
//             x=A->atom_coord[3*i_a+0];
//             y=A->atom_coord[3*i_a+1];
//             z=A->atom_coord[3*i_a+2];
//             
//             //by the 273 standard TrCAMM is taken in the electron population way, but not in the charge way 
//             c   = calc_atom_multipole(A_DM, DM_ij, cD.S_MO  , i_a, A->basis[0].ao_center, n_ao,0);//S_MO! there is no error
//             d_x = calc_atom_multipole(A_DM, DM_ij, cD.Dx_AO , i_a, A->basis[0].ao_center, n_ao,0);
//             d_y = calc_atom_multipole(A_DM, DM_ij, cD.Dy_AO , i_a, A->basis[0].ao_center, n_ao,0);
//             d_z = calc_atom_multipole(A_DM, DM_ij, cD.Dz_AO , i_a, A->basis[0].ao_center, n_ao,0);
//             q_xx= calc_atom_multipole(A_DM, DM_ij, cD.Qxx_AO, i_a, A->basis[0].ao_center, n_ao,0);
//             q_xy= calc_atom_multipole(A_DM, DM_ij, cD.Qxy_AO, i_a, A->basis[0].ao_center, n_ao,0);
//             q_xz= calc_atom_multipole(A_DM, DM_ij, cD.Qxz_AO, i_a, A->basis[0].ao_center, n_ao,0);
//             q_yy= calc_atom_multipole(A_DM, DM_ij, cD.Qyy_AO, i_a, A->basis[0].ao_center, n_ao,0);
//             q_yz= calc_atom_multipole(A_DM, DM_ij, cD.Qyz_AO, i_a, A->basis[0].ao_center, n_ao,0);
//             q_zz= calc_atom_multipole(A_DM, DM_ij, cD.Qzz_AO, i_a, A->basis[0].ao_center, n_ao,0);
//             
//             q_xx=q_xx-d_x*x-x*d_x+c*x*x;
//             q_xy=q_xy-d_x*y-x*d_y+c*x*y;
//             q_xz=q_xz-d_x*z-x*d_z+c*x*z;
//             q_yy=q_yy-d_y*y-y*d_y+c*y*y;
//             q_yz=q_yz-d_y*z-y*d_z+c*y*z;
//             q_zz=q_zz-d_z*z-z*d_z+c*z*z;
//             
//             d_x=d_x-c*x;
//             d_y=d_y-c*y;
//             d_z=d_z-c*z;
//             
//             c=c-S*A->nucl_charges[i_a];
//             
//             printf("%.10e %.10e %.10e  %.10e  %.10e %.10e %.10e  %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e \n",
//                    A->atom_coord[i_a*3  ],
//                    A->atom_coord[i_a*3+1],
//                    A->atom_coord[i_a*3+2],
//                    c,
//                    d_x, d_y, d_z,
//                    q_xx, q_xy, q_xz, q_xy, q_yy, q_yz, q_xz, q_yz, q_zz);
//             
//         }
//     }
//     
//     double * Px_AO = new double[n_ao*n_ao];
//     double * Py_AO = new double[n_ao*n_ao];
//     double * Pz_AO = new double[n_ao*n_ao];
//     for(int i2=0;i2<n_ao*n_ao;i2++){
//         Px_AO[i2]=cD.P_AO[3*i2+0];
//         Py_AO[i2]=cD.P_AO[3*i2+1];
//         Pz_AO[i2]=cD.P_AO[3*i2+2];
//     }
//     
//     double * MUx_AO = new double[n_ao*n_ao];
//     double * MUy_AO = new double[n_ao*n_ao];
//     double * MUz_AO = new double[n_ao*n_ao];
//     for(int i2=0;i2<n_ao*n_ao;i2++){
//         MUx_AO[i2]=cD.RxP_AO[9*i2+3*1+2]-cD.RxP_AO[9*i2+3*2+1];
//         MUy_AO[i2]=cD.RxP_AO[9*i2+3*2+0]-cD.RxP_AO[9*i2+3*0+2];
//         MUz_AO[i2]=cD.RxP_AO[9*i2+3*0+1]-cD.RxP_AO[9*i2+3*1+0];
//     }
//     
//     
//     
//     transform_to_libint_mat(Px_AO,A->basis, n_ao);
//     transform_to_libint_mat(Py_AO,A->basis, n_ao);
//     transform_to_libint_mat(Pz_AO,A->basis, n_ao);
//     transform_to_libint_mat(MUx_AO,A->basis, n_ao);
//     transform_to_libint_mat(MUy_AO,A->basis, n_ao);
//     transform_to_libint_mat(MUz_AO,A->basis, n_ao);
//     
//     
//     double p_x, p_y, p_z;
// //     double mu_x, mu_y, mu_z;
//     
// //     double * gr_p_x  =new double[A->n_atoms];
// //     double * gr_p_y  =new double[A->n_atoms];
// //     double * gr_p_z  =new double[A->n_atoms];
//     double * mu_x =new double[A->n_atoms];
//     double * mu_y =new double[A->n_atoms];
//     double * mu_z =new double[A->n_atoms];
//     
//     
//     printf("# Magnetic TrCAMM\n");
//     for(int i_s=0;i_s<n_s1;i_s++)
//     for(int j_s=0;j_s<n_s2;j_s++){
//         
//         printf("[%d %d]\n",i_s,j_s);
//         DM_ij = DM[i_s*n_s2+j_s];
//         
//         printf("# P_i\n");
//         for(int i_a=0;i_a<A->n_atoms;i_a++){
//             x=A->atom_coord[3*i_a+0];
//             y=A->atom_coord[3*i_a+1];
//             z=A->atom_coord[3*i_a+2];
//             // multiply on -1 !!!!!
//             p_x =  - calc_atom_multipole(A_DM, DM_ij, Px_AO , i_a, A->basis[0].ao_center, n_ao,1);
//             p_y =  - calc_atom_multipole(A_DM, DM_ij, Py_AO , i_a, A->basis[0].ao_center, n_ao,1);
//             p_z =  - calc_atom_multipole(A_DM, DM_ij, Pz_AO , i_a, A->basis[0].ao_center, n_ao,1);
//             // multiply on -1 !!!!
//             mu_x[i_a]= - calc_atom_multipole(A_DM, DM_ij, MUx_AO , i_a, A->basis[0].ao_center, n_ao,1);
//             mu_y[i_a]= - calc_atom_multipole(A_DM, DM_ij, MUy_AO , i_a, A->basis[0].ao_center, n_ao,1);
//             mu_z[i_a]= - calc_atom_multipole(A_DM, DM_ij, MUz_AO , i_a, A->basis[0].ao_center, n_ao,1);
//             
//             mu_x[i_a]+=-y*p_z+z*p_y;
//             mu_y[i_a]+=-z*p_x+x*p_z;
//             mu_z[i_a]+=-x*p_y+y*p_x;
//             
// //             if((i_s*n_s+j_s)==0){
// //                 gr_p_x [i_a] = p_x ;
// //                 gr_p_y [i_a] = p_y ;
// //                 gr_p_z [i_a] = p_z ;
// //             }
// //             if(i_s==j_s){
// //                 p_x -= gr_p_x [i_a];
// //                 p_y -= gr_p_y [i_a];
// //                 p_z -= gr_p_z [i_a];
// //             }
//                 
//             
//             
//             printf("%.10e %.10e %.10e \n",p_x, p_y, p_z);
//             
//         }
//         printf("# Mu_i\n");
//         for(int i_a=0;i_a<A->n_atoms;i_a++){
//             printf("%.10e %.10e %.10e \n",mu_x[i_a], mu_y[i_a], mu_z[i_a]);
//         }
//         
//     }
//     
//     
//     
//     
//     cD.clear_Q_AO();
//     
//     
//     
//     delete[] A_DM;
//     delete[]Px_AO;
//     delete[]Py_AO;
//     delete[]Pz_AO;
//     delete[]MUx_AO;
//     delete[]MUy_AO;
//     delete[]MUz_AO;
//     
//     
//     delete[] gr_c   ; 
//     delete[] gr_d_x ; 
//     delete[] gr_d_y ; 
//     delete[] gr_d_z ; 
//     delete[] gr_q_xx; 
//     delete[] gr_q_xy; 
//     delete[] gr_q_xz; 
//     delete[] gr_q_yy; 
//     delete[] gr_q_yz; 
//     delete[] gr_q_zz; 
//     delete[] mu_x;
//     delete[] mu_y;
//     delete[] mu_z;
//     
//     
//     return 0;
// }
