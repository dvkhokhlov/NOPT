# include <vector>
# include <libint2/shell.h>
# include <libint2/engine.h>
# include "matr.h"
# include <libint2.hpp>
# include <omp.h>
# include <chrono>
# include "blas_link.h"
# include "timer.h"
#include <unordered_map>
#include <mutex>
# include "libint_functions.h"
# include "ecp.h"
# include "common_vars.h"

using libint2::Shell;
using libint2::Engine;
using libint2::Operator;

std::vector<Engine> s_g_engines;//overlap
std::vector<Engine> t_g_engines;//kinetic
std::vector<Engine> v_g_engines;//NAI
// std::vector<Engine> m_g_engines;//MM-NAI
std::vector<Engine> e_g_engines;//ERI
std::vector<Engine> d_g_engines;//dipole
std::vector<Engine> q_g_engines;//quadrupole
std::vector<Engine> r_g_engines;//resolution of identity
std::vector<Engine> r2g_engines;//resolution of identity


int AO_1el_from_2shells(double * M, std::vector<Shell> s1, std::vector<Shell> s2, int dim1, int dim2, char E, int r){
    
    std::vector<Engine> * engines;
    
    if(E=='s')engines=&s_g_engines;
    if(E=='t')engines=&t_g_engines;
    if(E=='v')engines=&v_g_engines;
    if(E=='m')//engines=&m_g_engines;
    {
        printf("m-engine is deprecated\n");
        exit(1);
    }
    if(E=='d')engines=&d_g_engines;
    if(E=='q')engines=&q_g_engines;
    
    
    double * MP;
    MP = new double[num_threads*(dim1*dim2+64)];
    set_zero_matr(MP,num_threads*(dim1*dim2+64));
    
#pragma omp parallel 
    {
        int fi_num,fj_num;
        double * __restrict__ MT;//trying restrict??????
        int thread_id = omp_get_thread_num();
        MT = MP + thread_id*(dim1*dim2+64);
        const auto& buf_vec = engines->at(thread_id).results();
        fi_num=0;
        for(int i=0;i<s1.size();i++){
            fj_num=0;
            //fprintf(stderr,"AO[%7d/%7d]  threadid=%7d\r",i,s1.size(),thread_id);
//             getchar();
            for(int j=0;j<s2.size();j++){
                if (((i*s2.size()+j) % num_threads) == thread_id){
//                     auto deg = (i == j) ? 0.5 : 1.0;
                    engines->at(thread_id).compute(s1[i],s2[j]);
                    const auto* tmp_S = buf_vec[r];
//                     if (tmp_S == nullptr)
//                         continue;
                    if (tmp_S != nullptr)
                    for(int k=0;k<s1[i].contr[0].size();k++)
                    for(int l=0;l<s2[j].contr[0].size();l++){
                        MT[(fi_num+k)*dim2+fj_num+l]=tmp_S[k*(s2[j].contr[0].size())+l];
//                         MT[(fi_num+k)*dim1+fj_num+l]+=MT[(fi_num+k)*dim1+fj_num+l]
                    }
                    
                }
                fj_num+=s2[j].contr[0].size();
            }
            fi_num+=s1[i].contr[0].size();
        }
    }
    
    for(int i=0; i<dim1*dim2; i++){
        for(int j=0; j<num_threads; j++)
            M[i]+=MP[j*(dim1*dim2+64)+i];
        
    }
    delete[] MP;
    return 0;
}

int AO_1el_from_shell(double * M, std::vector<Shell> s, int dim, char E, int r){
    
    AO_1el_from_2shells(M, s, s, dim, dim, E, r);
    return 0;
}

int add_matr_from_2shells(double * M, std::vector<Shell> s1, std::vector<Shell> s2, int dim1, int dim2, std::vector<libint2::Engine> * engines, int r){
    
    double * MP;
    MP = new double[num_threads*dim1*dim2];
    set_zero_matr(MP,num_threads*dim1*dim2);
    
    auto lambda = [&](const int thread_id) {
        int fi_num,fj_num;
        double * MT;
        MT = MP + thread_id*dim1*dim2;
        const auto& buf_vec = engines->at(thread_id).results();
        fi_num=0;
        for(int i=0;i<s1.size();i++){
            fj_num=0;
            //fprintf(stderr,"AO[%7d/%7d]  threadid=%7d\r",i,s1.size(),thread_id);
//             getchar();
            for(int j=0;j<s2.size();j++){
                if (((i*s2.size()+j) % num_threads) == thread_id){
//                     auto deg = (i == j) ? 0.5 : 1.0;
                    engines->at(thread_id).compute(s1[i],s2[j]);
                    const auto* tmp_S = buf_vec[r];
//                     if (tmp_S == nullptr)
//                         continue;
                    if (tmp_S != nullptr)
                    for(int k=0;k<s1[i].contr[0].size();k++)
                    for(int l=0;l<s2[j].contr[0].size();l++){
                        MT[(fi_num+k)*dim2+fj_num+l]=tmp_S[k*(s2[j].contr[0].size())+l];
//                         MT[(fi_num+k)*dim1+fj_num+l]+=MT[(fi_num+k)*dim1+fj_num+l]
                    }
                    
                }
                fj_num+=s2[j].contr[0].size();
            }
            fi_num+=s1[i].contr[0].size();
        }
    };
    libint2_parallel_do(lambda);
    
    for(int i=0; i<dim1*dim2; i++){
        for(int j=0; j<num_threads; j++)
            M[i]+=MP[j*dim1*dim2+i];
        
    }
    delete[] MP;
    return 0;
}

int add_matr_from_shell(double * M, std::vector<Shell> s, int dim, std::vector<libint2::Engine> * engines, int r){
    
    add_matr_from_2shells(M, s, s, dim, dim, engines, r);
    return 0;
}

int double_array_clone(double *** S_th,double **S,int nt,int n_a,int dim){
    
    S_th[0]=S;   
    
    for(int i=1;i<nt;i++){
        S_th[i]=new double *[n_a];
        for(int j=0;j<n_a;j++)S_th[i][j]=new double[dim];
        for(int j=0;j<n_a;j++)set_zero_matr(S_th[i][j],dim);
    }
    
    return 0;
}

int clone_sum_and_clear(double *** S,int nt,int n_a,int dim){
    
    for(int i=1;i<nt ;i++)
    for(int j=0;j<n_a;j++)
    for(int k=0;k<dim;k++)
            S[0][j][k]+=S[i][j][k];
    
    for(int i=1;i<nt ;i++)
    for(int j=0;j<n_a;j++)
        delete[] S[i][j];
    
    for(int i=1;i<nt ;i++)
        delete[] S[i];
    
    
    return 1;
}

int gen_max_DM(const double * DM_I, int n_ao,double * DM_O, std::vector<size_t> n_bf,std::vector<size_t> map){
        
    int n = n_bf.size();
//     double d;
    double * t;
    int nbf1;
    int nbf2;
    t =new double[50*50];
    set_zero_matr(DM_O,n*n);
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++){
            nbf1 = n_bf[i];
            nbf2 = n_bf[j];
            set_zero_matr(t,50*50);
            for(int ii=0;ii<nbf1;ii++)
                for(int jj=0;jj<nbf2;jj++)
                    t[ii*50+jj]=DM_I[(map[i]+ii)*n_ao+map[j]+jj];
//                     for(int c = 0; c < s;c++){
//                         d = CI[c].max_el(map[i]+ii,map[j]+jj);
//                         if(t[ii*50+jj]<d)
//                             t[ii*50+jj]=d;
//                     }
            DM_O[i*n+j] = find_max_abs_value(t, nbf1, nbf2, 50);
        }
    
//     printf("max DM:\n");
//     PrintMatr(M,n,n,1);
    delete[] t;
    return 0;
    
}

int DM_to_F_transform(double ** J, double ** K,
                      double ** DM_J, int n_j,
                      double ** DM_K, int n_k, 
                      std::vector<Shell> shells, int n_ao){
  
    int nsh = shells.size();
    double ** Schwarz;
    
    std::vector<size_t> shell2bf[num_threads];
    shell2bf[0] = map_shell_to_basis_function(shells);
    for(int i=1;i<num_threads;i++)shell2bf[i]=shell2bf[0];
    
    std::vector<size_t> n_bf;
    n_bf = n_basis_function_in_shell(shells);
//     double ** DM_m;
//     DM_m = new double*[num_threads];
//     for(int i=0;i<num_threads;i++) DM_m[i] = new double[nsh*nsh];
    
//     gen_max_DM(DM,n_ao,DM_m[0],n_bf,shell2bf[0]);
//     for(int i=1;i<num_threads;i++)for(int j=0;j<nsh*nsh;j++) DM_m[i][j]=DM_m[0][j];
    
    double eps = 0;
    std::tie(obs_shellpair_list, obs_shellpair_data) = compute_shellpairs(shells,eps);
    
    
    double ** J_th[num_threads];
    double_array_clone(J_th,J,num_threads,n_j,n_ao*n_ao);
    double ** K_th[num_threads];
    double_array_clone(K_th,K,num_threads,n_k,n_ao*n_ao);
    
    
#pragma omp parallel 
    {        
        int thread_id=omp_get_thread_num();;
//         double * int_tr_tens_n=int_tr_tens[thread_id];
        int num_ints_computed = 0;
//         auto engine = engines->at(thread_id);
        const auto& buf = e_g_engines.at(thread_id).results();
        long s1234;
        int s1;
        for(s1=0, s1234=0; s1!=nsh; ++s1) {
          //fprintf(stderr,"AO[%7d/%7d]  threadid=%7d\r",s1,nsh,thread_id);
          auto bf1_first = shell2bf[thread_id][s1]; // first basis function in this shell
          auto n1 = shells[s1].size();   // number of basis functions in this shell
//           auto sp12_iter = obs_shellpair_data.at(s1).begin();
//           for (const auto s2 : obs_shellpair_list[s1]){
          for (auto s2=0;s2!=nsh;++s2){
            if (s2 > s1)
                  break;
            
            auto bf2_first = shell2bf[thread_id][s2];
            auto n2 = shells[s2].size();
            for(auto s3=0; s3<=s1; ++s3) {
//             for(auto s3=0; s3!=nsh; ++s3) {
            
              auto bf3_first = shell2bf[thread_id][s3];
              auto n3 = shells[s3].size();
//               auto sp34_iter = obs_shellpair_data.at(s3).begin();
              const auto s4_max = (s1 == s3) ? s2 : s3;
              for (const auto s4 : obs_shellpair_list[s3]) {
//               for (auto s4=0;s4!=nsh;++s4){
                if (s4 > s4_max)
                  break;
                  
                  
                if ((s1234++) % num_threads != thread_id) continue;
                auto bf4_first = shell2bf[thread_id][s4];
                auto n4 = shells[s4].size();
                // compute the permutational degeneracy (i.e. # of equivalents) of the given shell set
                auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
                auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
                auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
                auto s1234_deg = s12_deg * s34_deg * s12_34_deg/8;
              
        
//                 double DM_max = std::max(DM_m[thread_id][s1*nsh+s2],
//                                          DM_m[thread_id][s1*nsh+s3]);
//                 DM_max = std::max(DM_max,DM_m[thread_id][s1*nsh+s4]);
//                 DM_max = std::max(DM_max,DM_m[thread_id][s2*nsh+s1]);
//                 DM_max = std::max(DM_max,DM_m[thread_id][s2*nsh+s3]);
//                 DM_max = std::max(DM_max,DM_m[thread_id][s2*nsh+s4]);
//                 DM_max = std::max(DM_max,DM_m[thread_id][s3*nsh+s1]);
//                 DM_max = std::max(DM_max,DM_m[thread_id][s3*nsh+s2]);
//                 DM_max = std::max(DM_max,DM_m[thread_id][s3*nsh+s4]);
//                 DM_max = std::max(DM_max,DM_m[thread_id][s4*nsh+s1]);
//                 DM_max = std::max(DM_max,DM_m[thread_id][s4*nsh+s2]);
//                 DM_max = std::max(DM_max,DM_m[thread_id][s4*nsh+s3]);
              
        
//                 if (DM_max*Schwarz[thread_id][s1*nsh + s2] * Schwarz[thread_id][s3*nsh + s4] < 1E-11){
//                     delete[]  TEMP_Tensor[thread_id];
//                     continue;
//                 }
                e_g_engines.at(thread_id).compute(shells[s1], shells[s2], shells[s3], shells[s4]);
                const auto* buf_1234 = buf[0];
                
                
                if (buf_1234 == nullptr){
                    continue;
                } // if all integrals screened out, skip to next quartet
//                 num_ints_computed += n1 * n2 * n3 * n4;
        
        
            for(auto f1=0, f1234=0; f1!=n1; ++f1) {
              const auto bf1 = f1 + bf1_first;
              for(auto f2=0; f2!=n2; ++f2) {
                  const auto bf2 = f2 + bf2_first;
                  for(auto f3=0; f3!=n3; ++f3) {
                  const auto bf3 = f3 + bf3_first;
                      for(auto f4=0; f4!=n4; ++f4, ++f1234) {
                          const auto bf4 = f4 + bf4_first;
                          const auto value = buf_1234[f1234];
                          const auto value_scal_by_deg = value * s1234_deg;//<13|24>
                          
                          for(int i_j=0;i_j<n_j;i_j++){
                              J_th[thread_id][i_j][bf1*n_ao+bf2]+=value_scal_by_deg*DM_J[i_j][bf3*n_ao+bf4];
                              J_th[thread_id][i_j][bf1*n_ao+bf2]+=value_scal_by_deg*DM_J[i_j][bf4*n_ao+bf3];
                              J_th[thread_id][i_j][bf2*n_ao+bf1]+=value_scal_by_deg*DM_J[i_j][bf3*n_ao+bf4];
                              J_th[thread_id][i_j][bf2*n_ao+bf1]+=value_scal_by_deg*DM_J[i_j][bf4*n_ao+bf3];
                              J_th[thread_id][i_j][bf3*n_ao+bf4]+=value_scal_by_deg*DM_J[i_j][bf1*n_ao+bf2];
                              J_th[thread_id][i_j][bf4*n_ao+bf3]+=value_scal_by_deg*DM_J[i_j][bf1*n_ao+bf2];
                              J_th[thread_id][i_j][bf3*n_ao+bf4]+=value_scal_by_deg*DM_J[i_j][bf2*n_ao+bf1];
                              J_th[thread_id][i_j][bf4*n_ao+bf3]+=value_scal_by_deg*DM_J[i_j][bf2*n_ao+bf1];
                          }
                          for(int i_k=0;i_k<n_k;i_k++){////CHECK!!!! EXCHANGE may be wrong! previous version 1234->1432, current 1234->1324
                              K_th[thread_id][i_k][bf1*n_ao+bf3]+=value_scal_by_deg*DM_K[i_k][bf2*n_ao+bf4];
                              K_th[thread_id][i_k][bf1*n_ao+bf4]+=value_scal_by_deg*DM_K[i_k][bf2*n_ao+bf3];
                              K_th[thread_id][i_k][bf2*n_ao+bf3]+=value_scal_by_deg*DM_K[i_k][bf1*n_ao+bf4];
                              K_th[thread_id][i_k][bf2*n_ao+bf4]+=value_scal_by_deg*DM_K[i_k][bf1*n_ao+bf3];
                              K_th[thread_id][i_k][bf3*n_ao+bf1]+=value_scal_by_deg*DM_K[i_k][bf4*n_ao+bf2];
                              K_th[thread_id][i_k][bf4*n_ao+bf1]+=value_scal_by_deg*DM_K[i_k][bf3*n_ao+bf2];
                              K_th[thread_id][i_k][bf3*n_ao+bf2]+=value_scal_by_deg*DM_K[i_k][bf4*n_ao+bf1];
                              K_th[thread_id][i_k][bf4*n_ao+bf2]+=value_scal_by_deg*DM_K[i_k][bf3*n_ao+bf1];
                          }
                          
                      }
                    }
                  }
                }
        
              }
            }
          }
        }
    }
//     fprintf(stderr,"\n");
    //fprintf(stderr,"\r                           \r");
    clone_sum_and_clear(J_th, num_threads, n_j,n_ao*n_ao);
    clone_sum_and_clear(K_th, num_threads, n_k,n_ao*n_ao);
    
//     for(int i=0;i<num_threads;i++)delete[] Schwarz[i];
//     delete[] Schwarz;
    
    
    return 0;
}

double calc_2el_MO_INTS(std::vector<Shell> shells, int n_ao, const double * __restrict__ DM_D, double * J, double * K, double * MO_INTS, const double * __restrict__ R_VEC, const double * __restrict__ L_VEC, int n_v){

    
  int nsh = shells.size();
//   double ** Schwarz;
//   Schwarz = new double *[num_threads];
//   for(int i=0;i<num_threads;i++) Schwarz[i] = new double[nsh*nsh];
//   set_zero_matr(Schwarz[0],nsh*nsh);
//   compute_schwarz_matr(shells, num_threads, Schwarz[0]);
//   for(int i=1;i<num_threads;i++)for(int j=0;j<nsh*nsh;j++) Schwarz[i][j]=Schwarz[0][j];
  
  std::vector<size_t> shell2bf[num_threads];
  shell2bf[0] = map_shell_to_basis_function(shells);
  for(int i=1;i<num_threads;i++)shell2bf[i]=shell2bf[0];
  
  std::vector<size_t> n_bf;
  n_bf = n_basis_function_in_shell(shells);
  double ** DM_m;
  DM_m = new double*[num_threads];
  for(int i=0;i<num_threads;i++) DM_m[i] = new double[nsh*nsh];
  
  gen_max_DM(DM_D,n_ao,DM_m[0],n_bf,shell2bf[0]);
  
  for(int i=1;i<num_threads;i++)for(int j=0;j<nsh*nsh;j++) DM_m[i][j]=DM_m[0][j];
  
  double eps = 0;
  std::tie(obs_shellpair_list, obs_shellpair_data) = compute_shellpairs(shells,eps);
// printf("screened\n");
//   double * res;
//   res = new double[num_threads];
  // loop over permutationally-unique set of shells
  
  double * mo_ints[num_threads];
  mo_ints[0] = MO_INTS;
  for(int i=1;i<num_threads;i++)mo_ints[i]  = new double[n_v*n_v*n_v*n_v];
  for(int i=0;i<num_threads;i++)set_zero_matr(mo_ints[i],n_v*n_v*n_v*n_v);
  double * J_th[num_threads];
  J_th[0] = J;
  for(int i=1;i<num_threads;i++)J_th[i]  = new double[n_ao*n_ao];
  for(int i=0;i<num_threads;i++)set_zero_matr(J_th[i],n_ao*n_ao);
  double * K_th[num_threads];
  K_th[0] = K;
  for(int i=1;i<num_threads;i++)K_th[i]  = new double[n_ao*n_ao];
  for(int i=0;i<num_threads;i++)set_zero_matr(K_th[i],n_ao*n_ao);
  
  
  auto lambda = [&](const int thread_id) {

//       res[thread_id]=0;
      double * u4;
      double * v4;
      double * uv34;
      double * m_u;
      double * m_v;

      int mg = L_MAX*(L_MAX+1);
      int mg2 = mg*mg;
      int mg3 = mg2*mg;
      
      u4   = new double[n_v*        mg3];
      v4   = new double[n_v*        mg3];
      uv34 = new double[n_v*n_v*    mg2];
      m_u  = new double[n_v*n_v*n_v*mg ];
      m_v  = new double[n_v*n_v*n_v*mg ];
//       double * int_tr_tens_n=int_tr_tens[thread_id];
      int num_ints_computed = 0;
//       auto engine = engines->at(thread_id);
      const auto& buf = e_g_engines.at(thread_id).results();
      long s1234;
      int s1;
      for(s1=0, s1234=0; s1!=nsh; ++s1) {
        //fprintf(stderr,"AO[%7d/%7d]  threadid=%7d\r",s1,nsh,thread_id);
//         fflush(stderr);
//         if (s1 % num_threads != thread_id) continue;
        auto bf1_first = shell2bf[thread_id][s1]; // first basis function in this shell
        auto n1 = shells[s1].size();   // number of basis functions in this shell
        auto sp12_iter = obs_shellpair_data.at(s1).begin();
        set_zero_matr(m_u, n_v*n_v*n_v*mg);
        set_zero_matr(m_v, n_v*n_v*n_v*mg);
        for (const auto s2 : obs_shellpair_list[s1]){
          if (s2 > s1)
                  break;  
//           set_zero_matr(u4v3, n_v*n_v*mg2);
          set_zero_matr(uv34, n_v*n_v*mg2);
          auto bf2_first = shell2bf[thread_id][s2];
          auto n2 = shells[s2].size();
          for(auto s3=0; s3<=s1; ++s3) {
//             printf("s3 = %d\n",s3);
            auto bf3_first = shell2bf[thread_id][s3];
            auto n3 = shells[s3].size();
            auto sp34_iter = obs_shellpair_data.at(s3).begin();
//             printf("n2 = %d\n",n2);
            const auto s4_max = (s1 == s3) ? s2 : s3;
            set_zero_matr(u4,n_v*mg3);
            set_zero_matr(v4,n_v*mg3);
            for (const auto s4 : obs_shellpair_list[s3]) {
              if (s4 > s4_max)
                  break;

    
              if ((s1234++) % num_threads != thread_id) continue;
              auto bf4_first = shell2bf[thread_id][s4];
              auto n4 = shells[s4].size();
    
              // compute the permutational degeneracy (i.e. # of equivalents) of the given shell set
              auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
              auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
              auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
              auto s1234_deg = s12_deg * s34_deg * s12_34_deg/8;
              e_g_engines.at(thread_id).compute(shells[s1], shells[s2], shells[s3], shells[s4]);
              const auto* buf_1234 = buf[0]; //(12|34)
              
              
              if (buf_1234 == nullptr){
                  continue;
              } // if all integrals screened out, skip to next quartet
//               num_ints_computed += n1 * n2 * n3 * n4;


          for(auto f1=0, f1234=0, f123=0; f1!=n1; ++f1) {
            const auto bf1 = f1 + bf1_first;
            for(auto f2=0; f2!=n2; ++f2) {
                const auto bf2 = f2 + bf2_first;
                for(auto f3=0; f3!=n3; ++f3,++f123) {
                    const auto bf3 = f3 + bf3_first;
                    for(auto f4=0; f4!=n4; ++f4, ++f1234) {
                        const auto bf4 = f4 + bf4_first;
                        const auto value = buf_1234[f1234];
                        const auto value_scal_by_deg = value * s1234_deg;//<13|24>
//                          printf("i = %e\n",value_scal_by_deg);
//                          if(fabs(value_scal_by_deg)>1e-10)
                        double f_int=0;
                        
//                         f_int=2*(DM_D[bf1*n_ao+bf2]+DM_D[bf2*n_ao+bf1])*
//                                 (S*(DM_D[bf3*n_ao+bf4]+DM_D[bf4*n_ao+bf3])+DM_A[bf3*n_ao+bf4]+DM_A[bf4*n_ao+bf3]+DM_B[bf4*n_ao+bf3]+DM_B[bf3*n_ao+bf4])+
//                               2*(DM_D[bf3*n_ao+bf4]+DM_D[bf4*n_ao+bf3])*
//                                 (S*(DM_D[bf1*n_ao+bf2]+DM_D[bf2*n_ao+bf1])+DM_A[bf1*n_ao+bf2]+DM_A[bf2*n_ao+bf1]+DM_B[bf1*n_ao+bf2]+DM_B[bf2*n_ao+bf1])
//                                  -DM_D[bf1*n_ao+bf4]*(S*DM_D[bf3*n_ao+bf2]+DM_A[bf3*n_ao+bf2]+DM_B[bf3*n_ao+bf2])
//                                  -DM_D[bf1*n_ao+bf3]*(S*DM_D[bf4*n_ao+bf2]+DM_A[bf4*n_ao+bf2]+DM_B[bf4*n_ao+bf2])
//                                  -DM_D[bf2*n_ao+bf3]*(S*DM_D[bf4*n_ao+bf1]+DM_A[bf4*n_ao+bf1]+DM_B[bf4*n_ao+bf1])
//                                  -DM_D[bf2*n_ao+bf4]*(S*DM_D[bf3*n_ao+bf1]+DM_B[bf3*n_ao+bf1]+DM_A[bf3*n_ao+bf1])
//                                  -DM_D[bf3*n_ao+bf1]*(S*DM_D[bf2*n_ao+bf4]+DM_A[bf2*n_ao+bf4]+DM_B[bf2*n_ao+bf4])
//                                  -DM_D[bf3*n_ao+bf2]*(S*DM_D[bf1*n_ao+bf4]+DM_A[bf1*n_ao+bf4]+DM_B[bf1*n_ao+bf4])
//                                  -DM_D[bf4*n_ao+bf1]*(S*DM_D[bf2*n_ao+bf3]+DM_A[bf2*n_ao+bf3]+DM_B[bf2*n_ao+bf3])
//                                  -DM_D[bf4*n_ao+bf2]*(S*DM_D[bf1*n_ao+bf3]+DM_A[bf1*n_ao+bf3]+DM_B[bf1*n_ao+bf3]);
//                         for(int i_j=0;i_j<n_j;i_j++){
                            J_th[thread_id]/*[i_j]*/[bf1*n_ao+bf2]+=value_scal_by_deg*DM_D/*[i_j]*/[bf3*n_ao+bf4];
                            J_th[thread_id]/*[i_j]*/[bf1*n_ao+bf2]+=value_scal_by_deg*DM_D/*[i_j]*/[bf4*n_ao+bf3];
                            J_th[thread_id]/*[i_j]*/[bf2*n_ao+bf1]+=value_scal_by_deg*DM_D/*[i_j]*/[bf3*n_ao+bf4];
                            J_th[thread_id]/*[i_j]*/[bf2*n_ao+bf1]+=value_scal_by_deg*DM_D/*[i_j]*/[bf4*n_ao+bf3];
                            J_th[thread_id]/*[i_j]*/[bf3*n_ao+bf4]+=value_scal_by_deg*DM_D/*[i_j]*/[bf1*n_ao+bf2];
                            J_th[thread_id]/*[i_j]*/[bf4*n_ao+bf3]+=value_scal_by_deg*DM_D/*[i_j]*/[bf1*n_ao+bf2];
                            J_th[thread_id]/*[i_j]*/[bf3*n_ao+bf4]+=value_scal_by_deg*DM_D/*[i_j]*/[bf2*n_ao+bf1];
                            J_th[thread_id]/*[i_j]*/[bf4*n_ao+bf3]+=value_scal_by_deg*DM_D/*[i_j]*/[bf2*n_ao+bf1];
//                         }
//                         for(int i_k=0;i_k<n_k;i_k++){
                            K_th[thread_id]/*[i_k]*/[bf1*n_ao+bf3]+=value_scal_by_deg*DM_D/*[i_k]*/[bf2*n_ao+bf4];
                            K_th[thread_id]/*[i_k]*/[bf1*n_ao+bf4]+=value_scal_by_deg*DM_D/*[i_k]*/[bf2*n_ao+bf3];
                            K_th[thread_id]/*[i_k]*/[bf2*n_ao+bf3]+=value_scal_by_deg*DM_D/*[i_k]*/[bf1*n_ao+bf4];
                            K_th[thread_id]/*[i_k]*/[bf2*n_ao+bf4]+=value_scal_by_deg*DM_D/*[i_k]*/[bf1*n_ao+bf3];
                            K_th[thread_id]/*[i_k]*/[bf3*n_ao+bf1]+=value_scal_by_deg*DM_D/*[i_k]*/[bf4*n_ao+bf2];
                            K_th[thread_id]/*[i_k]*/[bf4*n_ao+bf1]+=value_scal_by_deg*DM_D/*[i_k]*/[bf3*n_ao+bf2];
                            K_th[thread_id]/*[i_k]*/[bf3*n_ao+bf2]+=value_scal_by_deg*DM_D/*[i_k]*/[bf4*n_ao+bf1];
                            K_th[thread_id]/*[i_k]*/[bf4*n_ao+bf2]+=value_scal_by_deg*DM_D/*[i_k]*/[bf3*n_ao+bf1];
//                         }         


                        for(int l_vo=0;l_vo<n_v;l_vo++){
                            v4[f123*n_v+l_vo]+=R_VEC[bf4*n_v+l_vo]*value_scal_by_deg;
                            u4[f123*n_v+l_vo]+=L_VEC[bf4*n_v+l_vo]*value_scal_by_deg;
                        }
//                         for(int l_vo=0;l_vo<n_v;l_vo++){
//                             u4[f123*n_v+l_vo]+=R_VEC[bf3*n_v+l_vo]*value_scal_by_deg;
//                         }
//                         res[thread_id] += value_scal_by_deg*f_int;
                      
                    }
                  }
                }
              }

            }
            for(auto f1=0, f12=0, f123=0; f1!=n1; ++f1) {
                const auto bf1 = f1 + bf1_first;
                for(auto f2=0; f2!=n2; ++f2,++f12) {
                    const auto bf2 = f2 + bf2_first;
                    for(auto f3=0; f3!=n3; ++f3,++f123) {
                        const auto bf3 = f3 + bf3_first;
                        for(int k_vo=0;k_vo<n_v;k_vo++)
                        for(int l_vo=0;l_vo<n_v;l_vo++){
                            uv34[(f12*n_v+k_vo)*n_v+l_vo]+=L_VEC[bf3*n_v+k_vo]*v4[f123*n_v+l_vo]
                                                         + R_VEC[bf3*n_v+l_vo]*u4[f123*n_v+k_vo];//!!!!!!!!!!!!!
                        }

                    }
                }
            }
            
          }
          for(auto f1=0, f12=0; f1!=n1; ++f1) {
              const auto bf1 = f1 + bf1_first;
              for(auto f2=0; f2!=n2; ++f2,++f12) {
                  const auto bf2 = f2 + bf2_first;
                  for(int j_vo=0;j_vo<n_v;j_vo++)
                  for(int k_vo=0;k_vo<n_v;k_vo++)
                  for(int l_vo=0;l_vo<n_v;l_vo++){
//                       fprintf(stderr,"%d\n",((f1*n_v+j_vo)*n_v+k_vo)*n_v+l_vo);
                      m_u[((f1*n_v+j_vo)*n_v+k_vo)*n_v+l_vo]+=R_VEC[bf2*n_v+j_vo]*uv34[(f12*n_v+k_vo)*n_v+l_vo];
                      m_v[((f1*n_v+j_vo)*n_v+k_vo)*n_v+l_vo]+=L_VEC[bf2*n_v+j_vo]*uv34[(f12*n_v+k_vo)*n_v+l_vo];
                      }

                  }
              }
          }
        for(auto f1=0; f1!=n1; ++f1) {
              const auto bf1 = f1 + bf1_first;
              for(int i_vo=0;i_vo<n_v;i_vo++)
              for(int j_vo=0;j_vo<n_v;j_vo++)
              for(int k_vo=0;k_vo<n_v;k_vo++)
              for(int l_vo=0;l_vo<n_v;l_vo++){
                  mo_ints[thread_id][((i_vo*n_v+j_vo)*n_v+k_vo)*n_v+l_vo]+=L_VEC[bf1*n_v+i_vo]*m_u[((f1*n_v+j_vo)*n_v+k_vo)*n_v+l_vo]
                                                                         + L_VEC[bf1*n_v+k_vo]*m_u[((f1*n_v+l_vo)*n_v+i_vo)*n_v+j_vo]
                                                                         + R_VEC[bf1*n_v+j_vo]*m_v[((f1*n_v+i_vo)*n_v+k_vo)*n_v+l_vo]
                                                                         + R_VEC[bf1*n_v+l_vo]*m_v[((f1*n_v+k_vo)*n_v+i_vo)*n_v+j_vo];
                  
                  }
              }                
      }
      
      delete[] u4  ;
      delete[] v4  ;
      delete[] uv34;
      delete[] m_u ;
      delete[] m_v ;
  };

//   libint2_parallel_do(lambda1);
//    printf_timer("after lambda");
  libint2_parallel_do(lambda);
  //fprintf(stderr,"\r                                        \r");
  printf_timer("after AO_ints");
//   printf("lambda done\n");
  double result = 0.0;
//   for(int i=0;i<num_threads;i++)result+=res[i];
//   printf("get thread sum moints\n");
  for(int i=1;i<num_threads;i++)
  for(int j=0;j<n_v*n_v*n_v*n_v;j++)
      mo_ints[0][j]+=mo_ints[i][j];
  for(int i=1;i<num_threads;i++)
  for(int j=0;j<n_ao*n_ao;j++)
      J_th[0][j]+=J_th[i][j];
  for(int i=1;i<num_threads;i++)
  for(int j=0;j<n_ao*n_ao;j++)
      K_th[0][j]+=K_th[i][j];
//   printf("get thread sum moints - done\n");
//   for(int a=0;a<n_v    ;a++)
//   for(int b=0;b<n_v    ;b++)
//   for(int k=0;k<n_v*n_v;k++)
//       TH_AB[(a*n_ao+b)*n_v*n_v+k]+=mo_ints[0][(a*n_v+b)*n_v*n_v+k]+
//                                    mo_ints[0][(b*n_v+a)*n_v*n_v+k];
  
  for(int i=1;i<num_threads;i++)delete[] mo_ints[i];
  for(int i=1;i<num_threads;i++)delete[] J_th[i];
  for(int i=1;i<num_threads;i++)delete[] K_th[i];
//   printf("free moints done\n");
//   for(int i=0;i<num_threads;i++)
//   for(int k_vo=0;k_vo<n_v ;k_vo++)
//   for(int l_vo=0;l_vo<n_v ;l_vo++)
//   for(int c_ao=0;c_ao<n_ao;c_ao++)
//   for(int d_ao=0;d_ao<n_ao;d_ao++)
//             result+=int_tr_tens[i][((c_ao*n_ao+d_ao)*n_v+k_vo)*n_v+l_vo]*R_VEC[c_ao*n_v+k_vo]*R_VEC[d_ao*n_v+l_vo];
  

//   for(int i=0;i<num_threads;i++)delete[] Schwarz[i];
//   delete[] Schwarz;
  for(int i=0;i<num_threads;i++) delete[] DM_m[i];
  delete[] DM_m;
//   printf("exit 2el block\n");
  
  return result;
}

double calc_2el_MO_INTS_multi_JK(std::vector<Shell> shells, int n_ao, double ** __restrict__ DM_J, double ** J, int n_j, double ** __restrict__ DM_K, double ** K, int n_k, double * MO_INTS, const double * __restrict__ R_VEC, const double * __restrict__ L_VEC, int n_v, std::vector<libint2::Engine> * engines){

  
  int nsh = shells.size();
  double ** Schwarz;
  Schwarz = new double *[num_threads];
  for(int i=0;i<num_threads;i++) Schwarz[i] = new double[nsh*nsh];
  set_zero_matr(Schwarz[0],nsh*nsh);
  compute_schwarz_matr(shells, num_threads, Schwarz[0]);
  for(int i=1;i<num_threads;i++)for(int j=0;j<nsh*nsh;j++) Schwarz[i][j]=Schwarz[0][j];
  
  std::vector<size_t> shell2bf[num_threads];
  shell2bf[0] = map_shell_to_basis_function(shells);
  for(int i=1;i<num_threads;i++)shell2bf[i]=shell2bf[0];
  
  std::vector<size_t> n_bf;
  n_bf = n_basis_function_in_shell(shells);
//   double ** DM_m;
//   DM_m = new double*[num_threads];
//   for(int i=0;i<num_threads;i++) DM_m[i] = new double[nsh*nsh];
//   
//   gen_max_DM(DM_D,n_ao,DM_m[0],n_bf,shell2bf[0]);
//   
//   for(int i=1;i<num_threads;i++)for(int j=0;j<nsh*nsh;j++) DM_m[i][j]=DM_m[0][j];
  
  double eps = 0;
  std::tie(obs_shellpair_list, obs_shellpair_data) = compute_shellpairs(shells,eps);
// printf("screened\n");
//   double * res;
//   res = new double[num_threads];
  // loop over permutationally-unique set of shells
  
  double * mo_ints[num_threads];
  mo_ints[0] = MO_INTS;
  for(int i=1;i<num_threads;i++)mo_ints[i]  = new double[n_v*n_v*n_v*n_v];
  for(int i=0;i<num_threads;i++)set_zero_matr(mo_ints[i],n_v*n_v*n_v*n_v);
  double ** J_th[num_threads];
  double_array_clone(J_th,J,num_threads,n_j,n_ao*n_ao);
  double ** K_th[num_threads];
  double_array_clone(K_th,K,num_threads,n_k,n_ao*n_ao);
  
  
  auto lambda = [&](const int thread_id) {

//       res[thread_id]=0;
      double * u4;
      double * v4;
      double * uv34;
      double * m_u;
      double * m_v;

      int mg = L_MAX*(L_MAX+1);
      int mg2 = mg*mg;
      int mg3 = mg2*mg;
      
      u4   = new double[n_v*        mg3];
      v4   = new double[n_v*        mg3];
      uv34 = new double[n_v*n_v*    mg2];
      m_u  = new double[n_v*n_v*n_v*mg ];
      m_v  = new double[n_v*n_v*n_v*mg ];
//       double * int_tr_tens_n=int_tr_tens[thread_id];
      int num_ints_computed = 0;
//       auto engine = engines->at(thread_id);
      const auto& buf = engines->at(thread_id).results();
      long s1234;
      int s1;
      for(s1=0, s1234=0; s1!=nsh; ++s1) {
        //fprintf(stderr,"AO[%7d/%7d]  threadid=%7d\r",s1,nsh,thread_id);
//         fflush(stderr);
//         if (s1 % num_threads != thread_id) continue;
        auto bf1_first = shell2bf[thread_id][s1]; // first basis function in this shell
        auto n1 = shells[s1].size();   // number of basis functions in this shell
        auto sp12_iter = obs_shellpair_data.at(s1).begin();
        set_zero_matr(m_u, n_v*n_v*n_v*mg);
        set_zero_matr(m_v, n_v*n_v*n_v*mg);
        for (const auto s2 : obs_shellpair_list[s1]){
          if (s2 > s1)
                  break;  
//           set_zero_matr(u4v3, n_v*n_v*mg2);
          set_zero_matr(uv34, n_v*n_v*mg2);
          auto bf2_first = shell2bf[thread_id][s2];
          auto n2 = shells[s2].size();
          for(auto s3=0; s3<=s1; ++s3) {
//             printf("s3 = %d\n",s3);
            auto bf3_first = shell2bf[thread_id][s3];
            auto n3 = shells[s3].size();
            auto sp34_iter = obs_shellpair_data.at(s3).begin();
//             printf("n2 = %d\n",n2);
            const auto s4_max = (s1 == s3) ? s2 : s3;
            set_zero_matr(u4,n_v*mg3);
            set_zero_matr(v4,n_v*mg3);
            for (const auto s4 : obs_shellpair_list[s3]) {
              if (s4 > s4_max)
                  break;

    
              if ((s1234++) % num_threads != thread_id) continue;
              auto bf4_first = shell2bf[thread_id][s4];
              auto n4 = shells[s4].size();
    
              // compute the permutational degeneracy (i.e. # of equivalents) of the given shell set
              auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
              auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
              auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
              auto s1234_deg = s12_deg * s34_deg * s12_34_deg/8;
              engines->at(thread_id).compute(shells[s1], shells[s2], shells[s3], shells[s4]);
              const auto* buf_1234 = buf[0]; //(12|34)
              
              
              if (buf_1234 == nullptr){
                  continue;
              } // if all integrals screened out, skip to next quartet
//               num_ints_computed += n1 * n2 * n3 * n4;


          for(auto f1=0, f1234=0, f123=0; f1!=n1; ++f1) {
            const auto bf1 = f1 + bf1_first;
            for(auto f2=0; f2!=n2; ++f2) {
                const auto bf2 = f2 + bf2_first;
                for(auto f3=0; f3!=n3; ++f3,++f123) {
                    const auto bf3 = f3 + bf3_first;
                    for(auto f4=0; f4!=n4; ++f4, ++f1234) {
                        const auto bf4 = f4 + bf4_first;
                        const auto value = buf_1234[f1234];
                        const auto value_scal_by_deg = value * s1234_deg;//<13|24>
//                          printf("i = %e\n",value_scal_by_deg);
//                          if(fabs(value_scal_by_deg)>1e-10)
                        double f_int=0;
                        
//                         f_int=2*(DM_D[bf1*n_ao+bf2]+DM_D[bf2*n_ao+bf1])*
//                                 (S*(DM_D[bf3*n_ao+bf4]+DM_D[bf4*n_ao+bf3])+DM_A[bf3*n_ao+bf4]+DM_A[bf4*n_ao+bf3]+DM_B[bf4*n_ao+bf3]+DM_B[bf3*n_ao+bf4])+
//                               2*(DM_D[bf3*n_ao+bf4]+DM_D[bf4*n_ao+bf3])*
//                                 (S*(DM_D[bf1*n_ao+bf2]+DM_D[bf2*n_ao+bf1])+DM_A[bf1*n_ao+bf2]+DM_A[bf2*n_ao+bf1]+DM_B[bf1*n_ao+bf2]+DM_B[bf2*n_ao+bf1])
//                                  -DM_D[bf1*n_ao+bf4]*(S*DM_D[bf3*n_ao+bf2]+DM_A[bf3*n_ao+bf2]+DM_B[bf3*n_ao+bf2])
//                                  -DM_D[bf1*n_ao+bf3]*(S*DM_D[bf4*n_ao+bf2]+DM_A[bf4*n_ao+bf2]+DM_B[bf4*n_ao+bf2])
//                                  -DM_D[bf2*n_ao+bf3]*(S*DM_D[bf4*n_ao+bf1]+DM_A[bf4*n_ao+bf1]+DM_B[bf4*n_ao+bf1])
//                                  -DM_D[bf2*n_ao+bf4]*(S*DM_D[bf3*n_ao+bf1]+DM_B[bf3*n_ao+bf1]+DM_A[bf3*n_ao+bf1])
//                                  -DM_D[bf3*n_ao+bf1]*(S*DM_D[bf2*n_ao+bf4]+DM_A[bf2*n_ao+bf4]+DM_B[bf2*n_ao+bf4])
//                                  -DM_D[bf3*n_ao+bf2]*(S*DM_D[bf1*n_ao+bf4]+DM_A[bf1*n_ao+bf4]+DM_B[bf1*n_ao+bf4])
//                                  -DM_D[bf4*n_ao+bf1]*(S*DM_D[bf2*n_ao+bf3]+DM_A[bf2*n_ao+bf3]+DM_B[bf2*n_ao+bf3])
//                                  -DM_D[bf4*n_ao+bf2]*(S*DM_D[bf1*n_ao+bf3]+DM_A[bf1*n_ao+bf3]+DM_B[bf1*n_ao+bf3]);
                        for(int i_j=0;i_j<n_j;i_j++){
                            J_th[thread_id][i_j][bf1*n_ao+bf2]+=value_scal_by_deg*DM_J[i_j][bf3*n_ao+bf4];
                            J_th[thread_id][i_j][bf1*n_ao+bf2]+=value_scal_by_deg*DM_J[i_j][bf4*n_ao+bf3];
                            J_th[thread_id][i_j][bf2*n_ao+bf1]+=value_scal_by_deg*DM_J[i_j][bf3*n_ao+bf4];
                            J_th[thread_id][i_j][bf2*n_ao+bf1]+=value_scal_by_deg*DM_J[i_j][bf4*n_ao+bf3];
                            J_th[thread_id][i_j][bf3*n_ao+bf4]+=value_scal_by_deg*DM_J[i_j][bf1*n_ao+bf2];
                            J_th[thread_id][i_j][bf4*n_ao+bf3]+=value_scal_by_deg*DM_J[i_j][bf1*n_ao+bf2];
                            J_th[thread_id][i_j][bf3*n_ao+bf4]+=value_scal_by_deg*DM_J[i_j][bf2*n_ao+bf1];
                            J_th[thread_id][i_j][bf4*n_ao+bf3]+=value_scal_by_deg*DM_J[i_j][bf2*n_ao+bf1];
                        }
                        for(int i_k=0;i_k<n_k;i_k++){
                            K_th[thread_id][i_k][bf1*n_ao+bf3]+=value_scal_by_deg*DM_K[i_k][bf2*n_ao+bf4];
                            K_th[thread_id][i_k][bf1*n_ao+bf4]+=value_scal_by_deg*DM_K[i_k][bf2*n_ao+bf3];
                            K_th[thread_id][i_k][bf2*n_ao+bf3]+=value_scal_by_deg*DM_K[i_k][bf1*n_ao+bf4];
                            K_th[thread_id][i_k][bf2*n_ao+bf4]+=value_scal_by_deg*DM_K[i_k][bf1*n_ao+bf3];
                            K_th[thread_id][i_k][bf3*n_ao+bf1]+=value_scal_by_deg*DM_K[i_k][bf4*n_ao+bf2];
                            K_th[thread_id][i_k][bf4*n_ao+bf1]+=value_scal_by_deg*DM_K[i_k][bf3*n_ao+bf2];
                            K_th[thread_id][i_k][bf3*n_ao+bf2]+=value_scal_by_deg*DM_K[i_k][bf4*n_ao+bf1];
                            K_th[thread_id][i_k][bf4*n_ao+bf2]+=value_scal_by_deg*DM_K[i_k][bf3*n_ao+bf1];
                        }         


                        for(int l_vo=0;l_vo<n_v;l_vo++){
                            v4[f123*n_v+l_vo]+=R_VEC[bf4*n_v+l_vo]*value_scal_by_deg;
                            u4[f123*n_v+l_vo]+=L_VEC[bf4*n_v+l_vo]*value_scal_by_deg;
                        }
//                         for(int l_vo=0;l_vo<n_v;l_vo++){
//                             u4[f123*n_v+l_vo]+=R_VEC[bf3*n_v+l_vo]*value_scal_by_deg;
//                         }
//                         res[thread_id] += value_scal_by_deg*f_int;
                      
                    }
                  }
                }
              }

            }
            for(auto f1=0, f12=0, f123=0; f1!=n1; ++f1) {
                const auto bf1 = f1 + bf1_first;
                for(auto f2=0; f2!=n2; ++f2,++f12) {
                    const auto bf2 = f2 + bf2_first;
                    for(auto f3=0; f3!=n3; ++f3,++f123) {
                        const auto bf3 = f3 + bf3_first;
                        for(int k_vo=0;k_vo<n_v;k_vo++)
                        for(int l_vo=0;l_vo<n_v;l_vo++){
                            uv34[(f12*n_v+k_vo)*n_v+l_vo]+=L_VEC[bf3*n_v+k_vo]*v4[f123*n_v+l_vo]
                                                         + R_VEC[bf3*n_v+l_vo]*u4[f123*n_v+k_vo];//!!!!!!!!!!!!!
                        }

                    }
                }
            }
            
          }
          for(auto f1=0, f12=0; f1!=n1; ++f1) {
              const auto bf1 = f1 + bf1_first;
              for(auto f2=0; f2!=n2; ++f2,++f12) {
                  const auto bf2 = f2 + bf2_first;
                  for(int j_vo=0;j_vo<n_v;j_vo++)
                  for(int k_vo=0;k_vo<n_v;k_vo++)
                  for(int l_vo=0;l_vo<n_v;l_vo++){
                      m_u[((f1*n_v+j_vo)*n_v+k_vo)*n_v+l_vo]+=R_VEC[bf2*n_v+j_vo]*uv34[(f12*n_v+k_vo)*n_v+l_vo];
                      m_v[((f1*n_v+j_vo)*n_v+k_vo)*n_v+l_vo]+=L_VEC[bf2*n_v+j_vo]*uv34[(f12*n_v+k_vo)*n_v+l_vo];
                      }

                  }
              }
          }
        for(auto f1=0; f1!=n1; ++f1) {
              const auto bf1 = f1 + bf1_first;
              for(int i_vo=0;i_vo<n_v;i_vo++)
              for(int j_vo=0;j_vo<n_v;j_vo++)
              for(int k_vo=0;k_vo<n_v;k_vo++)
              for(int l_vo=0;l_vo<n_v;l_vo++){
                  mo_ints[thread_id][((i_vo*n_v+j_vo)*n_v+k_vo)*n_v+l_vo]+=L_VEC[bf1*n_v+i_vo]*m_u[((f1*n_v+j_vo)*n_v+k_vo)*n_v+l_vo]
                                                                         + L_VEC[bf1*n_v+k_vo]*m_u[((f1*n_v+l_vo)*n_v+i_vo)*n_v+j_vo]
                                                                         + R_VEC[bf1*n_v+j_vo]*m_v[((f1*n_v+i_vo)*n_v+k_vo)*n_v+l_vo]
                                                                         + R_VEC[bf1*n_v+l_vo]*m_v[((f1*n_v+k_vo)*n_v+i_vo)*n_v+j_vo];
                  
                  }
              }                
      }
      
      delete[] u4  ;
      delete[] v4  ;
      delete[] uv34;
      delete[] m_u ;
      delete[] m_v ;
  };

//   libint2_parallel_do(lambda1);
//    printf_timer("after lambda");
  libint2_parallel_do(lambda);
  //fprintf(stderr,"\r                                        \r");
  printf_timer("after AO_ints");
//   printf("lambda done\n");
  double result = 0.0;
//   for(int i=0;i<num_threads;i++)result+=res[i];
//   printf("get thread sum moints\n");
  for(int i=1;i<num_threads;i++)
  for(int j=0;j<n_v*n_v*n_v*n_v;j++)
      mo_ints[0][j]+=mo_ints[i][j];
  clone_sum_and_clear(J_th, num_threads, n_j,n_ao*n_ao);
  clone_sum_and_clear(K_th, num_threads, n_k,n_ao*n_ao);
  
  for(int i=1;i<num_threads;i++)delete[] mo_ints[i];
  

  for(int i=0;i<num_threads;i++)delete[] Schwarz[i];
  delete[] Schwarz;
//   for(int i=0;i<num_threads;i++) delete[] DM_m[i];
//   delete[] DM_m;
//   printf("exit 2el block\n");
  
  return result;
}


double calc_2el_MO_INTS_4sets(std::vector<Shell> shells, int n_ao, double ** __restrict__ DM_J, double ** J, int n_j, double ** __restrict__ DM_K, double ** K, int n_k, double * MO_INTS, const double * __restrict__ R_VEC1, const double * __restrict__ R_VEC2, const double * __restrict__ L_VEC1, const double * __restrict__ L_VEC2, int n_v, std::vector<libint2::Engine> * engines){

  
  int nsh = shells.size();
  double ** Schwarz;
  Schwarz = new double *[num_threads];
  for(int i=0;i<num_threads;i++) Schwarz[i] = new double[nsh*nsh];
  set_zero_matr(Schwarz[0],nsh*nsh);
  compute_schwarz_matr(shells, num_threads, Schwarz[0]);
  for(int i=1;i<num_threads;i++)for(int j=0;j<nsh*nsh;j++) Schwarz[i][j]=Schwarz[0][j];
  
  std::vector<size_t> shell2bf[num_threads];
  shell2bf[0] = map_shell_to_basis_function(shells);
  for(int i=1;i<num_threads;i++)shell2bf[i]=shell2bf[0];
  
  std::vector<size_t> n_bf;
  n_bf = n_basis_function_in_shell(shells);
//   double ** DM_m;
//   DM_m = new double*[num_threads];
//   for(int i=0;i<num_threads;i++) DM_m[i] = new double[nsh*nsh];
  
//   gen_max_DM(DM_D,n_ao,DM_m[0],n_bf,shell2bf[0]);
  
//   for(int i=1;i<num_threads;i++)for(int j=0;j<nsh*nsh;j++) DM_m[i][j]=DM_m[0][j];
  
  double eps = 0;
  std::tie(obs_shellpair_list, obs_shellpair_data) = compute_shellpairs(shells,eps);
// printf("screened\n");
//   double * res;
//   res = new double[num_threads];
  // loop over permutationally-unique set of shells
  
  double * mo_ints[num_threads];
  mo_ints[0] = MO_INTS;
  for(int i=1;i<num_threads;i++)mo_ints[i]  = new double[n_v*n_v*n_v*n_v];
  for(int i=0;i<num_threads;i++)set_zero_matr(mo_ints[i],n_v*n_v*n_v*n_v);

  double ** J_th[num_threads];
  double_array_clone(J_th,J,num_threads,n_j,n_ao*n_ao);
  double ** K_th[num_threads];
  double_array_clone(K_th,K,num_threads,n_k,n_ao*n_ao);
    
  
  auto lambda = [&](const int thread_id) {

//       res[thread_id]=0;
      double * t4;
      double * u4;
      double * v4;
      double * w4;
      double * tu34;
      double * vw34;
      double * m_t;
      double * m_u;
      double * m_v;
      double * m_w;

      int mg = L_MAX*(L_MAX+1);
      int mg2 = mg*mg;
      int mg3 = mg2*mg;
      
      t4   = new double[n_v*        mg3];
      u4   = new double[n_v*        mg3];
      v4   = new double[n_v*        mg3];
      w4   = new double[n_v*        mg3];
      
      tu34 = new double[n_v*n_v*    mg2];
      vw34 = new double[n_v*n_v*    mg2];
      
      m_t  = new double[n_v*n_v*n_v*mg/*2*/];
      m_u  = new double[n_v*n_v*n_v*mg/*2*/];
      m_v  = new double[n_v*n_v*n_v*mg/*2*/];
      m_w  = new double[n_v*n_v*n_v*mg/*2*/];
      
//       double * int_tr_tens_n=int_tr_tens[thread_id];
      int num_ints_computed = 0;
//       auto engine = engines->at(thread_id);
      const auto& buf = engines->at(thread_id).results();
      long s1234;
      int s1;
      for(s1=0, s1234=0; s1!=nsh; ++s1) {
        //fprintf(stderr,"AO[%7d/%7d]  threadid=%7d\r",s1,nsh,thread_id);
//         fflush(stderr);
//         if (s1 % num_threads != thread_id) continue;
        auto bf1_first = shell2bf[thread_id][s1]; // first basis function in this shell
        auto n1 = shells[s1].size();   // number of basis functions in this shell
        auto sp12_iter = obs_shellpair_data.at(s1).begin();
        set_zero_matr(m_t, n_v*n_v*n_v*mg/*2*/);
        set_zero_matr(m_u, n_v*n_v*n_v*mg/*2*/);
        set_zero_matr(m_v, n_v*n_v*n_v*mg/*2*/);
        set_zero_matr(m_w, n_v*n_v*n_v*mg/*2*/);
        for (const auto s2 : obs_shellpair_list[s1]){
          if (s2 > s1)
                  break;  
//           set_zero_matr(u4v3, n_v*n_v*mg2);
          set_zero_matr(tu34, n_v*n_v*mg2);
          set_zero_matr(vw34, n_v*n_v*mg2);
          auto bf2_first = shell2bf[thread_id][s2];
          auto n2 = shells[s2].size();
          for(auto s3=0; s3<=s1; ++s3) {
//             printf("s3 = %d\n",s3);
            auto bf3_first = shell2bf[thread_id][s3];
            auto n3 = shells[s3].size();
            auto sp34_iter = obs_shellpair_data.at(s3).begin();
//             printf("n2 = %d\n",n2);
            const auto s4_max = (s1 == s3) ? s2 : s3;
            set_zero_matr(t4,n_v*mg3);
            set_zero_matr(u4,n_v*mg3);
            set_zero_matr(v4,n_v*mg3);
            set_zero_matr(w4,n_v*mg3);
            for (const auto s4 : obs_shellpair_list[s3]) {
              if (s4 > s4_max)
                  break;

    
              if ((s1234++) % num_threads != thread_id) continue;
              auto bf4_first = shell2bf[thread_id][s4];
              auto n4 = shells[s4].size();
    
              // compute the permutational degeneracy (i.e. # of equivalents) of the given shell set
              auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
              auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
              auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
              auto s1234_deg = s12_deg * s34_deg * s12_34_deg/8;
              engines->at(thread_id).compute(shells[s1], shells[s2], shells[s3], shells[s4]);
              const auto* buf_1234 = buf[0]; //(12|34)
              
              
              if (buf_1234 == nullptr){
                  continue;
              } // if all integrals screened out, skip to next quartet
//               num_ints_computed += n1 * n2 * n3 * n4;


          for(auto f1=0, f1234=0, f123=0; f1!=n1; ++f1) {
            const auto bf1 = f1 + bf1_first;
            for(auto f2=0; f2!=n2; ++f2) {
                const auto bf2 = f2 + bf2_first;
                for(auto f3=0; f3!=n3; ++f3,++f123) {
                    const auto bf3 = f3 + bf3_first;
                    for(auto f4=0; f4!=n4; ++f4, ++f1234) {
                        const auto bf4 = f4 + bf4_first;
                        const auto value = buf_1234[f1234];
                        const auto value_scal_by_deg = value * s1234_deg;//<13|24>
//                          printf("i = %e\n",value_scal_by_deg);
//                          if(fabs(value_scal_by_deg)>1e-10)
                        double f_int=0;
                        
//                         f_int=2*(DM_D[bf1*n_ao+bf2]+DM_D[bf2*n_ao+bf1])*
//                                 (S*(DM_D[bf3*n_ao+bf4]+DM_D[bf4*n_ao+bf3])+DM_A[bf3*n_ao+bf4]+DM_A[bf4*n_ao+bf3]+DM_B[bf4*n_ao+bf3]+DM_B[bf3*n_ao+bf4])+
//                               2*(DM_D[bf3*n_ao+bf4]+DM_D[bf4*n_ao+bf3])*
//                                 (S*(DM_D[bf1*n_ao+bf2]+DM_D[bf2*n_ao+bf1])+DM_A[bf1*n_ao+bf2]+DM_A[bf2*n_ao+bf1]+DM_B[bf1*n_ao+bf2]+DM_B[bf2*n_ao+bf1])
//                                  -DM_D[bf1*n_ao+bf4]*(S*DM_D[bf3*n_ao+bf2]+DM_A[bf3*n_ao+bf2]+DM_B[bf3*n_ao+bf2])
//                                  -DM_D[bf1*n_ao+bf3]*(S*DM_D[bf4*n_ao+bf2]+DM_A[bf4*n_ao+bf2]+DM_B[bf4*n_ao+bf2])
//                                  -DM_D[bf2*n_ao+bf3]*(S*DM_D[bf4*n_ao+bf1]+DM_A[bf4*n_ao+bf1]+DM_B[bf4*n_ao+bf1])
//                                  -DM_D[bf2*n_ao+bf4]*(S*DM_D[bf3*n_ao+bf1]+DM_B[bf3*n_ao+bf1]+DM_A[bf3*n_ao+bf1])
//                                  -DM_D[bf3*n_ao+bf1]*(S*DM_D[bf2*n_ao+bf4]+DM_A[bf2*n_ao+bf4]+DM_B[bf2*n_ao+bf4])
//                                  -DM_D[bf3*n_ao+bf2]*(S*DM_D[bf1*n_ao+bf4]+DM_A[bf1*n_ao+bf4]+DM_B[bf1*n_ao+bf4])
//                                  -DM_D[bf4*n_ao+bf1]*(S*DM_D[bf2*n_ao+bf3]+DM_A[bf2*n_ao+bf3]+DM_B[bf2*n_ao+bf3])
//                                  -DM_D[bf4*n_ao+bf2]*(S*DM_D[bf1*n_ao+bf3]+DM_A[bf1*n_ao+bf3]+DM_B[bf1*n_ao+bf3]);
                        for(int i_j=0;i_j<n_j;i_j++){
                            J_th[thread_id][i_j][bf1*n_ao+bf2]+=value_scal_by_deg*DM_J[i_j][bf3*n_ao+bf4];
                            J_th[thread_id][i_j][bf1*n_ao+bf2]+=value_scal_by_deg*DM_J[i_j][bf4*n_ao+bf3];
                            J_th[thread_id][i_j][bf2*n_ao+bf1]+=value_scal_by_deg*DM_J[i_j][bf3*n_ao+bf4];
                            J_th[thread_id][i_j][bf2*n_ao+bf1]+=value_scal_by_deg*DM_J[i_j][bf4*n_ao+bf3];
                            J_th[thread_id][i_j][bf3*n_ao+bf4]+=value_scal_by_deg*DM_J[i_j][bf1*n_ao+bf2];
                            J_th[thread_id][i_j][bf4*n_ao+bf3]+=value_scal_by_deg*DM_J[i_j][bf1*n_ao+bf2];
                            J_th[thread_id][i_j][bf3*n_ao+bf4]+=value_scal_by_deg*DM_J[i_j][bf2*n_ao+bf1];
                            J_th[thread_id][i_j][bf4*n_ao+bf3]+=value_scal_by_deg*DM_J[i_j][bf2*n_ao+bf1];
                        }
                        for(int i_k=0;i_k<n_k;i_k++){
                            K_th[thread_id][i_k][bf1*n_ao+bf3]+=value_scal_by_deg*DM_K[i_k][bf2*n_ao+bf4];
                            K_th[thread_id][i_k][bf1*n_ao+bf4]+=value_scal_by_deg*DM_K[i_k][bf2*n_ao+bf3];
                            K_th[thread_id][i_k][bf2*n_ao+bf3]+=value_scal_by_deg*DM_K[i_k][bf1*n_ao+bf4];
                            K_th[thread_id][i_k][bf2*n_ao+bf4]+=value_scal_by_deg*DM_K[i_k][bf1*n_ao+bf3];
                            K_th[thread_id][i_k][bf3*n_ao+bf1]+=value_scal_by_deg*DM_K[i_k][bf4*n_ao+bf2];
                            K_th[thread_id][i_k][bf4*n_ao+bf1]+=value_scal_by_deg*DM_K[i_k][bf3*n_ao+bf2];
                            K_th[thread_id][i_k][bf3*n_ao+bf2]+=value_scal_by_deg*DM_K[i_k][bf4*n_ao+bf1];
                            K_th[thread_id][i_k][bf4*n_ao+bf2]+=value_scal_by_deg*DM_K[i_k][bf3*n_ao+bf1];
                        }         


                        for(int l_vo=0;l_vo<n_v;l_vo++){
                            t4[f123*n_v+l_vo]+=L_VEC1[bf4*n_v+l_vo]*value_scal_by_deg;//l<->i
                            u4[f123*n_v+l_vo]+=R_VEC1[bf4*n_v+l_vo]*value_scal_by_deg;//l<->j
                            v4[f123*n_v+l_vo]+=L_VEC2[bf4*n_v+l_vo]*value_scal_by_deg;//l<->k
                            w4[f123*n_v+l_vo]+=R_VEC2[bf4*n_v+l_vo]*value_scal_by_deg;//l<->l
                        }
//                         for(int l_vo=0;l_vo<n_v;l_vo++){
//                             u4[f123*n_v+l_vo]+=R_VEC1[bf3*n_v+l_vo]*value_scal_by_deg;
//                         }
//                         res[thread_id] += value_scal_by_deg*f_int;
                      
                    }
                  }
                }
              }

            }
            for(auto f1=0, f12=0, f123=0; f1!=n1; ++f1) {
                const auto bf1 = f1 + bf1_first;
                for(auto f2=0; f2!=n2; ++f2,++f12) {
                    const auto bf2 = f2 + bf2_first;
                    for(auto f3=0; f3!=n3; ++f3,++f123) {
                        const auto bf3 = f3 + bf3_first;
                        for(int k_vo=0;k_vo<n_v;k_vo++)
                        for(int l_vo=0;l_vo<n_v;l_vo++){
                            tu34[(f12*n_v+k_vo)*n_v+l_vo]+=L_VEC1[bf3*n_v+k_vo]*u4[f123*n_v+l_vo]
                                                         + R_VEC1[bf3*n_v+l_vo]*t4[f123*n_v+k_vo];//k<->i,l<->j
                            vw34[(f12*n_v+k_vo)*n_v+l_vo]+=L_VEC2[bf3*n_v+k_vo]*w4[f123*n_v+l_vo]
                                                         + R_VEC2[bf3*n_v+l_vo]*v4[f123*n_v+k_vo];//k<->k,l<->l
                        }

                    }
                }
            }
            
          }
          for(auto f1=0, f12=0; f1!=n1; ++f1) {
              const auto bf1 = f1 + bf1_first;
              for(auto f2=0; f2!=n2; ++f2,++f12) {
                  const auto bf2 = f2 + bf2_first;
                  for(int j_vo=0;j_vo<n_v;j_vo++)
                  for(int k_vo=0;k_vo<n_v;k_vo++)
                  for(int l_vo=0;l_vo<n_v;l_vo++){
                      m_t[((f1*n_v+j_vo)*n_v+k_vo)*n_v+l_vo]+=R_VEC1[bf2*n_v+j_vo]*vw34[(f12*n_v+k_vo)*n_v+l_vo];//k<->k,l<->l,j<->j
                      m_u[((f1*n_v+j_vo)*n_v+k_vo)*n_v+l_vo]+=L_VEC1[bf2*n_v+j_vo]*vw34[(f12*n_v+k_vo)*n_v+l_vo];//k<->k,l<->l,j<->i
                      m_v[((f1*n_v+j_vo)*n_v+k_vo)*n_v+l_vo]+=R_VEC2[bf2*n_v+j_vo]*tu34[(f12*n_v+k_vo)*n_v+l_vo];//k<->i,l<->j,j<->l
                      m_w[((f1*n_v+j_vo)*n_v+k_vo)*n_v+l_vo]+=L_VEC2[bf2*n_v+j_vo]*tu34[(f12*n_v+k_vo)*n_v+l_vo];//k<->i,l<->j,j<->k
                      }

                  }
              }
          }
        for(auto f1=0; f1!=n1; ++f1) {
              const auto bf1 = f1 + bf1_first;
              for(int i_vo=0;i_vo<n_v;i_vo++)
              for(int j_vo=0;j_vo<n_v;j_vo++)
              for(int k_vo=0;k_vo<n_v;k_vo++)
              for(int l_vo=0;l_vo<n_v;l_vo++){//(L1 R1|L2 R2)
                  mo_ints[thread_id][((i_vo*n_v+j_vo)*n_v+k_vo)*n_v+l_vo]+=L_VEC1[bf1*n_v+i_vo]*m_t[((f1*n_v+j_vo)*n_v+k_vo)*n_v+l_vo]
                                                                         + R_VEC1[bf1*n_v+j_vo]*m_u[((f1*n_v+i_vo)*n_v+k_vo)*n_v+l_vo]
                                                                         + L_VEC2[bf1*n_v+k_vo]*m_v[((f1*n_v+l_vo)*n_v+i_vo)*n_v+j_vo]
                                                                         + R_VEC2[bf1*n_v+l_vo]*m_w[((f1*n_v+k_vo)*n_v+i_vo)*n_v+j_vo];
                  
                  }
              }                
      }
      
      delete[] t4  ;
      delete[] u4  ;
      delete[] v4  ;
      delete[] w4  ;
      delete[] tu34;
      delete[] vw34;
      delete[] m_t ;
      delete[] m_u ;
      delete[] m_v ;
      delete[] m_w ;
  };

//   libint2_parallel_do(lambda1);
//    printf_timer("after lambda");
  libint2_parallel_do(lambda);
  //fprintf(stderr,"\r                                        \r");
  printf_timer("after AO_ints");
//   printf("lambda done\n");
  double result = 0.0;
//   for(int i=0;i<num_threads;i++)result+=res[i];
//   printf("get thread sum moints\n");
  for(int i=1;i<num_threads;i++)
  for(int j=0;j<n_v*n_v*n_v*n_v;j++)
      mo_ints[0][j]+=mo_ints[i][j];
  
  clone_sum_and_clear(J_th, num_threads, n_j,n_ao*n_ao);
  clone_sum_and_clear(K_th, num_threads, n_k,n_ao*n_ao);
  
  for(int i=1;i<num_threads;i++)delete[] mo_ints[i];

  for(int i=0;i<num_threads;i++)delete[] Schwarz[i];
  delete[] Schwarz;
//   for(int i=0;i<num_threads;i++) delete[] DM_m[i];
//   delete[] DM_m;
//   printf("exit 2el block\n");
  
  return result;
}

double calc_2el_MO_INTS_4sets_2shells(std::vector<Shell> shells1,std::vector<Shell> shells2, int n_ao1, int n_ao2, double ** __restrict__ DM_J, double ** J, int n_j, double ** __restrict__ DM_K, double ** K, int n_k, double * MO_INTS, const double * __restrict__ R_VEC1, const double * __restrict__ R_VEC2, const double * __restrict__ L_VEC1, const double * __restrict__ L_VEC2, int n_v, std::vector<libint2::Engine> * engines){

  
  int nsh1 = shells1.size();
  int nsh2 = shells2.size();
  
  double ** Schwarz;
  Schwarz = new double *[num_threads];
  for(int i=0;i<num_threads;i++) Schwarz[i] = new double[nsh1*nsh1];
  set_zero_matr(Schwarz[0],nsh1*nsh1);
  compute_schwarz_matr(shells1, num_threads, Schwarz[0]);
  for(int i=1;i<num_threads;i++)for(int j=0;j<nsh1*nsh1;j++) Schwarz[i][j]=Schwarz[0][j];
  
  std::vector<size_t> shell2bf1[num_threads];
  shell2bf1[0] = map_shell_to_basis_function(shells1);
  for(int i=1;i<num_threads;i++)shell2bf1[i]=shell2bf1[0];
  
  std::vector<size_t> shell2bf2[num_threads];
  shell2bf2[0] = map_shell_to_basis_function(shells2);
  for(int i=1;i<num_threads;i++)shell2bf2[i]=shell2bf2[0];
  
  
//   std::vector<size_t> n_bf;
//   n_bf = n_basis_function_in_shell(shells);
//   double ** DM_m;
//   DM_m = new double*[num_threads];
//   for(int i=0;i<num_threads;i++) DM_m[i] = new double[nsh*nsh];
  
//   gen_max_DM(DM_D,n_ao,DM_m[0],n_bf,shell2bf[0]);
  
//   for(int i=1;i<num_threads;i++)for(int j=0;j<nsh*nsh;j++) DM_m[i][j]=DM_m[0][j];
  
  double eps = 0;
  std::tie(obs_shellpair_list1, obs_shellpair_data1) = compute_shellpairs(shells1,eps);
  std::tie(obs_shellpair_list2, obs_shellpair_data2) = compute_shellpairs(shells2,eps);
// printf("screened\n");
//   double * res;
//   res = new double[num_threads];
  // loop over permutationally-unique set of shells
  
  double * mo_ints[num_threads];
  mo_ints[0] = MO_INTS;
  for(int i=1;i<num_threads;i++)mo_ints[i]  = new double[n_v*n_v*n_v*n_v];
  for(int i=0;i<num_threads;i++)set_zero_matr(mo_ints[i],n_v*n_v*n_v*n_v);

  double ** J_th[num_threads];
  double_array_clone(J_th,J,num_threads,n_j,n_ao1*n_ao2);
  double ** K_th[num_threads];
  double_array_clone(K_th,K,num_threads,n_k,n_ao1*n_ao2);
    
  
  auto lambda = [&](const int thread_id) {

//       res[thread_id]=0;
      double * t4;
      double * u4;
      double * v4;
      double * w4;
      double * tu34;
      double * vw34;
      double * m_t;
      double * m_u;
      double * m_v;
      double * m_w;

      int mg = L_MAX*(L_MAX+1);
      int mg2 = mg*mg;
      int mg3 = mg2*mg;
      
//       t4   = new double[n_v*        mg3];
//       u4   = new double[n_v*        mg3];
      v4   = new double[n_v*        mg3];
      w4   = new double[n_v*        mg3];
      
//       tu34 = new double[n_v*n_v*    mg2];
      vw34 = new double[n_v*n_v*    mg2];
      
      m_t  = new double[n_v*n_v*n_v*mg/*2*/];
      m_u  = new double[n_v*n_v*n_v*mg/*2*/];
//       m_v  = new double[n_v*n_v*n_v*mg/*2*/];
//       m_w  = new double[n_v*n_v*n_v*mg/*2*/];
      
//       double * int_tr_tens_n=int_tr_tens[thread_id];
      int num_ints_computed = 0;
//       auto engine = engines->at(thread_id);
      const auto& buf = engines->at(thread_id).results();
      long s1234;
      int s1;
      for(s1=0, s1234=0; s1!=nsh1; ++s1) {
        //fprintf(stderr,"AO[%7d/%7d]  threadid=%7d\r",s1,nsh1,thread_id);
//         fflush(stderr);
//         if (s1 % num_threads != thread_id) continue;
        auto bf1_first = shell2bf1[thread_id][s1]; // first basis function in this shell
        auto n1 = shells1[s1].size();   // number of basis functions in this shell
        auto sp12_iter = obs_shellpair_data1.at(s1).begin();
        set_zero_matr(m_t, n_v*n_v*n_v*mg/*2*/);
        set_zero_matr(m_u, n_v*n_v*n_v*mg/*2*/);
//         set_zero_matr(m_v, n_v*n_v*n_v*mg/*2*/);
//         set_zero_matr(m_w, n_v*n_v*n_v*mg/*2*/);
        for (const auto s2 : obs_shellpair_list1[s1]){
          if (s2 > s1)
                  break;  

//           set_zero_matr(tu34, n_v*n_v*mg2);
          set_zero_matr(vw34, n_v*n_v*mg2);
          auto bf2_first = shell2bf1[thread_id][s2];
          auto n2 = shells1[s2].size();
          for(auto s3=0; s3!=nsh2; ++s3) {
//             printf("s3 = %d\n",s3);
            auto bf3_first = shell2bf2[thread_id][s3];
            auto n3 = shells2[s3].size();
            auto sp34_iter = obs_shellpair_data2.at(s3).begin();
//             printf("n2 = %d\n",n2);
//             const auto s4_max = (s1 == s3) ? s2 : s3;
//             set_zero_matr(t4,n_v*mg3);
//             set_zero_matr(u4,n_v*mg3);
            set_zero_matr(v4,n_v*mg3);
            set_zero_matr(w4,n_v*mg3);
            for (const auto s4 : obs_shellpair_list2[s3]) {
              if (s4 > s3)
                  break;

    
              if ((s1234++) % num_threads != thread_id) continue;
              auto bf4_first = shell2bf2[thread_id][s4];
              auto n4 = shells2[s4].size();
    
              // compute the permutational degeneracy (i.e. # of equivalents) of the given shell set
              auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
              auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
//               auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
              auto s1234_deg = s12_deg * s34_deg /4;
              engines->at(thread_id).compute(shells1[s1], shells1[s2], shells2[s3], shells2[s4]);
              const auto* buf_1234 = buf[0]; //(12|34)
              
              
              if (buf_1234 == nullptr){
                  continue;
              } // if all integrals screened out, skip to next quartet
//               num_ints_computed += n1 * n2 * n3 * n4;


          for(auto f1=0, f1234=0, f123=0; f1!=n1; ++f1) {
            const auto bf1 = f1 + bf1_first;
            for(auto f2=0; f2!=n2; ++f2) {
                const auto bf2 = f2 + bf2_first;
                for(auto f3=0; f3!=n3; ++f3,++f123) {
                    const auto bf3 = f3 + bf3_first;
                    for(auto f4=0; f4!=n4; ++f4, ++f1234) {
                        const auto bf4 = f4 + bf4_first;
                        const auto value = buf_1234[f1234];
                        const auto value_scal_by_deg = value * s1234_deg;//<13|24>
//                          printf("i = %e\n",value_scal_by_deg);
//                          if(fabs(value_scal_by_deg)>1e-10)
                        double f_int=0;
                        
//                         f_int=2*(DM_D[bf1*n_ao+bf2]+DM_D[bf2*n_ao+bf1])*
//                                 (S*(DM_D[bf3*n_ao+bf4]+DM_D[bf4*n_ao+bf3])+DM_A[bf3*n_ao+bf4]+DM_A[bf4*n_ao+bf3]+DM_B[bf4*n_ao+bf3]+DM_B[bf3*n_ao+bf4])+
//                               2*(DM_D[bf3*n_ao+bf4]+DM_D[bf4*n_ao+bf3])*
//                                 (S*(DM_D[bf1*n_ao+bf2]+DM_D[bf2*n_ao+bf1])+DM_A[bf1*n_ao+bf2]+DM_A[bf2*n_ao+bf1]+DM_B[bf1*n_ao+bf2]+DM_B[bf2*n_ao+bf1])
//                                  -DM_D[bf1*n_ao+bf4]*(S*DM_D[bf3*n_ao+bf2]+DM_A[bf3*n_ao+bf2]+DM_B[bf3*n_ao+bf2])
//                                  -DM_D[bf1*n_ao+bf3]*(S*DM_D[bf4*n_ao+bf2]+DM_A[bf4*n_ao+bf2]+DM_B[bf4*n_ao+bf2])
//                                  -DM_D[bf2*n_ao+bf3]*(S*DM_D[bf4*n_ao+bf1]+DM_A[bf4*n_ao+bf1]+DM_B[bf4*n_ao+bf1])
//                                  -DM_D[bf2*n_ao+bf4]*(S*DM_D[bf3*n_ao+bf1]+DM_B[bf3*n_ao+bf1]+DM_A[bf3*n_ao+bf1])
//                                  -DM_D[bf3*n_ao+bf1]*(S*DM_D[bf2*n_ao+bf4]+DM_A[bf2*n_ao+bf4]+DM_B[bf2*n_ao+bf4])
//                                  -DM_D[bf3*n_ao+bf2]*(S*DM_D[bf1*n_ao+bf4]+DM_A[bf1*n_ao+bf4]+DM_B[bf1*n_ao+bf4])
//                                  -DM_D[bf4*n_ao+bf1]*(S*DM_D[bf2*n_ao+bf3]+DM_A[bf2*n_ao+bf3]+DM_B[bf2*n_ao+bf3])
//                                  -DM_D[bf4*n_ao+bf2]*(S*DM_D[bf1*n_ao+bf3]+DM_A[bf1*n_ao+bf3]+DM_B[bf1*n_ao+bf3]);
//                         for(int i_j=0;i_j<n_j;i_j++){
//                             J_th[thread_id][i_j][bf1*n_ao2+bf2]+=value_scal_by_deg*DM_J[i_j][bf3*n_ao2+bf4];
//                             J_th[thread_id][i_j][bf1*n_ao2+bf2]+=value_scal_by_deg*DM_J[i_j][bf4*n_ao2+bf3];
//                             J_th[thread_id][i_j][bf2*n_ao2+bf1]+=value_scal_by_deg*DM_J[i_j][bf3*n_ao2+bf4];
//                             J_th[thread_id][i_j][bf2*n_ao2+bf1]+=value_scal_by_deg*DM_J[i_j][bf4*n_ao2+bf3];
//                             J_th[thread_id][i_j][bf3*n_ao2+bf4]+=value_scal_by_deg*DM_J[i_j][bf1*n_ao2+bf2];
//                             J_th[thread_id][i_j][bf4*n_ao2+bf3]+=value_scal_by_deg*DM_J[i_j][bf1*n_ao2+bf2];
//                             J_th[thread_id][i_j][bf3*n_ao2+bf4]+=value_scal_by_deg*DM_J[i_j][bf2*n_ao2+bf1];
//                             J_th[thread_id][i_j][bf4*n_ao2+bf3]+=value_scal_by_deg*DM_J[i_j][bf2*n_ao2+bf1];
//                         }
//                         for(int i_k=0;i_k<n_k;i_k++){
//                             K_th[thread_id][i_k][bf1*n_ao2+bf3]+=value_scal_by_deg*DM_K[i_k][bf2*n_ao2+bf4];
//                             K_th[thread_id][i_k][bf1*n_ao2+bf4]+=value_scal_by_deg*DM_K[i_k][bf2*n_ao2+bf3];
//                             K_th[thread_id][i_k][bf2*n_ao2+bf3]+=value_scal_by_deg*DM_K[i_k][bf1*n_ao2+bf4];
//                             K_th[thread_id][i_k][bf2*n_ao2+bf4]+=value_scal_by_deg*DM_K[i_k][bf1*n_ao2+bf3];
//                             K_th[thread_id][i_k][bf3*n_ao2+bf1]+=value_scal_by_deg*DM_K[i_k][bf4*n_ao2+bf2];
//                             K_th[thread_id][i_k][bf4*n_ao2+bf1]+=value_scal_by_deg*DM_K[i_k][bf3*n_ao2+bf2];
//                             K_th[thread_id][i_k][bf3*n_ao2+bf2]+=value_scal_by_deg*DM_K[i_k][bf4*n_ao2+bf1];
//                             K_th[thread_id][i_k][bf4*n_ao2+bf2]+=value_scal_by_deg*DM_K[i_k][bf3*n_ao2+bf1];
//                         }         
// 

                        for(int l_vo=0;l_vo<n_v;l_vo++){
//                             t4[f123*n_v+l_vo]+=L_VEC1[bf4*n_v+l_vo]*value_scal_by_deg;//l<->i
//                             u4[f123*n_v+l_vo]+=R_VEC1[bf4*n_v+l_vo]*value_scal_by_deg;//l<->j
                            v4[f123*n_v+l_vo]+=L_VEC2[bf4*n_v+l_vo]*value_scal_by_deg;//l<->k
                            w4[f123*n_v+l_vo]+=R_VEC2[bf4*n_v+l_vo]*value_scal_by_deg;//l<->l
                        }
//                         for(int l_vo=0;l_vo<n_v;l_vo++){
//                             u4[f123*n_v+l_vo]+=R_VEC1[bf3*n_v+l_vo]*value_scal_by_deg;
//                         }
//                         res[thread_id] += value_scal_by_deg*f_int;
                      
                    }
                  }
                }
              }

            }
            for(auto f1=0, f12=0, f123=0; f1!=n1; ++f1) {
                const auto bf1 = f1 + bf1_first;
                for(auto f2=0; f2!=n2; ++f2,++f12) {
                    const auto bf2 = f2 + bf2_first;
                    for(auto f3=0; f3!=n3; ++f3,++f123) {
                        const auto bf3 = f3 + bf3_first;
                        for(int k_vo=0;k_vo<n_v;k_vo++)
                        for(int l_vo=0;l_vo<n_v;l_vo++){
//                             tu34[(f12*n_v+k_vo)*n_v+l_vo]+=L_VEC1[bf3*n_v+k_vo]*u4[f123*n_v+l_vo]
//                                                          + R_VEC1[bf3*n_v+l_vo]*t4[f123*n_v+k_vo];//k<->i,l<->j
                            vw34[(f12*n_v+k_vo)*n_v+l_vo]+=L_VEC2[bf3*n_v+k_vo]*w4[f123*n_v+l_vo]
                                                         + R_VEC2[bf3*n_v+l_vo]*v4[f123*n_v+k_vo];//k<->k,l<->l
                        }

                    }
                }
            }
            
          }
          for(auto f1=0, f12=0; f1!=n1; ++f1) {
              const auto bf1 = f1 + bf1_first;
              for(auto f2=0; f2!=n2; ++f2,++f12) {
                  const auto bf2 = f2 + bf2_first;
                  for(int j_vo=0;j_vo<n_v;j_vo++)
                  for(int k_vo=0;k_vo<n_v;k_vo++)
                  for(int l_vo=0;l_vo<n_v;l_vo++){
                      m_t[((f1*n_v+j_vo)*n_v+k_vo)*n_v+l_vo]+=R_VEC1[bf2*n_v+j_vo]*vw34[(f12*n_v+k_vo)*n_v+l_vo];//k<->k,l<->l,j<->j
                      m_u[((f1*n_v+j_vo)*n_v+k_vo)*n_v+l_vo]+=L_VEC1[bf2*n_v+j_vo]*vw34[(f12*n_v+k_vo)*n_v+l_vo];//k<->k,l<->l,j<->i
//                       m_v[((f1*n_v+j_vo)*n_v+k_vo)*n_v+l_vo]+=R_VEC2[bf2*n_v+j_vo]*tu34[(f12*n_v+k_vo)*n_v+l_vo];//k<->i,l<->j,j<->l
//                       m_w[((f1*n_v+j_vo)*n_v+k_vo)*n_v+l_vo]+=L_VEC2[bf2*n_v+j_vo]*tu34[(f12*n_v+k_vo)*n_v+l_vo];//k<->i,l<->j,j<->k
                      }

                  }
              }
          }
        for(auto f1=0; f1!=n1; ++f1) {
              const auto bf1 = f1 + bf1_first;
              for(int i_vo=0;i_vo<n_v;i_vo++)
              for(int j_vo=0;j_vo<n_v;j_vo++)
              for(int k_vo=0;k_vo<n_v;k_vo++)
              for(int l_vo=0;l_vo<n_v;l_vo++){//(L1 R1|L2 R2)
                  mo_ints[thread_id][((i_vo*n_v+j_vo)*n_v+k_vo)*n_v+l_vo]+=L_VEC1[bf1*n_v+i_vo]*m_t[((f1*n_v+j_vo)*n_v+k_vo)*n_v+l_vo]
                                                                         + R_VEC1[bf1*n_v+j_vo]*m_u[((f1*n_v+i_vo)*n_v+k_vo)*n_v+l_vo]
                                                                         /*+ L_VEC2[bf1*n_v+k_vo]*m_v[((f1*n_v+l_vo)*n_v+i_vo)*n_v+j_vo]*/
                                                                         /*+ R_VEC2[bf1*n_v+l_vo]*m_w[((f1*n_v+k_vo)*n_v+i_vo)*n_v+j_vo]*/;
                  
                  }
              }                
      }
      
//       delete[] t4  ;
//       delete[] u4  ;
      delete[] v4  ;
      delete[] w4  ;
//       delete[] tu34;
      delete[] vw34;
      delete[] m_t ;
      delete[] m_u ;
//       delete[] m_v ;
//       delete[] m_w ;
  };

//   libint2_parallel_do(lambda1);
//    printf_timer("after lambda");
  libint2_parallel_do(lambda);
  //fprintf(stderr,"\r                                        \r");
  printf_timer("after AO_ints");
//   printf("lambda done\n");
  double result = 0.0;
//   for(int i=0;i<num_threads;i++)result+=res[i];
//   printf("get thread sum moints\n");
  for(int i=1;i<num_threads;i++)
  for(int j=0;j<n_v*n_v*n_v*n_v;j++)
      mo_ints[0][j]+=mo_ints[i][j];
  
  clone_sum_and_clear(J_th, num_threads, n_j,n_ao1*n_ao2);
  clone_sum_and_clear(K_th, num_threads, n_k,n_ao1*n_ao2);
  
  for(int i=1;i<num_threads;i++)delete[] mo_ints[i];

  for(int i=0;i<num_threads;i++)delete[] Schwarz[i];
  delete[] Schwarz;
//   for(int i=0;i<num_threads;i++) delete[] DM_m[i];
//   delete[] DM_m;
//   printf("exit 2el block\n");
  
  return result;
}

double calc_2el_MO_INTS_XYXY(std::vector<Shell>  shells, int n_ao, double * XYXY_INTS, const double * __restrict__ X_VEC, const double * __restrict__ Y_VEC, int n_x, int n_y){
  
  //nthreads = nt;
  
  int nsh = shells.size();
//   double ** Schwarz;
//   Schwarz = new double *[num_threads];
//   for(int i=0;i<num_threads;i++) Schwarz[i] = new double[nsh*nsh];
//   set_zero_matr(Schwarz[0],nsh*nsh);
//   compute_schwarz_matr(shells, num_threads, Schwarz[0]);
//   for(int i=1;i<num_threads;i++)for(int j=0;j<nsh*nsh;j++) Schwarz[i][j]=Schwarz[0][j];
  
  std::vector<size_t> shell2bf[num_threads];
  shell2bf[0] = map_shell_to_basis_function(shells);
  for(int i=1;i<num_threads;i++)shell2bf[i]=shell2bf[0];
  
  std::vector<size_t> n_bf;
  n_bf = n_basis_function_in_shell(shells);
//   double ** DM_m;
//   DM_m = new double*[num_threads];
//   for(int i=0;i<num_threads;i++) DM_m[i] = new double[nsh*nsh];
  
//   gen_max_DM(DM_D,n_ao,DM_m[0],n_bf,shell2bf[0]);
  
//   for(int i=1;i<num_threads;i++)for(int j=0;j<nsh*nsh;j++) DM_m[i][j]=DM_m[0][j];
  
  double eps = 0;
  std::tie(obs_shellpair_list, obs_shellpair_data) = compute_shellpairs(shells,eps);
  
  double * mo_ints[num_threads];
  mo_ints[0] = XYXY_INTS;
  for(int i=1;i<num_threads;i++)mo_ints[i]  = new double[n_x*n_x*n_y*n_y];
  for(int i=0;i<num_threads;i++)set_zero_matr(mo_ints[i],n_x*n_x*n_y*n_y);
  
  auto lambda = [&](const int thread_id) {

//       res[thread_id]=0;
      
      std::vector<Shell> shells_local = shells;
      
      double * __restrict__ x4;
      double * __restrict__ xy34;
      double * __restrict__ m_y;

      int mg = L_MAX*(L_MAX+1);
      int mg2 = mg*mg;
      int mg3 = mg2*mg;
      
      x4   = new double[n_x*        mg3];
      xy34 = new double[n_x*n_x*    mg2];
      m_y  = new double[n_x*n_y*n_x*mg];
      int num_ints_computed = 0;

      const auto& buf = e_g_engines.at(thread_id).results();
      long s1234;
      int s1;
      for(s1=0, s1234=0; s1!=nsh; ++s1) {
        //fprintf(stderr,"AO[%7d/%7d]  threadid=%7d\r",s1,nsh,thread_id);
        if (s1 % num_threads != thread_id) continue;
        auto bf1_first = shell2bf[thread_id][s1]; // first basis function in this shell
        auto n1 = shells_local[s1].size();   // number of basis functions in this shell
        auto sp12_iter = obs_shellpair_data.at(s1).begin();
        set_zero_matr(m_y, n_x*n_y*n_x*mg);
        for(auto s2=0; s2<=s1; ++s2) {
          set_zero_matr(xy34, n_x*n_x*mg2);
          auto bf2_first = shell2bf[thread_id][s2];
          auto n2 = shells_local[s2].size();
          for (const auto s3 : obs_shellpair_list[s1]){
            auto bf3_first = shell2bf[thread_id][s3];
            auto n3 = shells_local[s3].size();
            auto sp34_iter = obs_shellpair_data.at(s3).begin();
            const auto s4_max = (s1 == s3) ? s2 : s3;
            set_zero_matr(x4,n_x*mg3);
            for (const auto s4 : obs_shellpair_list[s2]) {
              
              auto bf4_first = shell2bf[thread_id][s4];
              auto n4 = shells_local[s4].size();
    
              auto s12_deg = (s1 == s2) ? 0.5 : 1.0;

              e_g_engines.at(thread_id).compute(shells_local[s1], shells_local[s3], shells_local[s2], shells_local[s4]);
              const auto* buf_1234 = buf[0]; //(12|34)
              
              
              if (buf_1234 == nullptr){
                  continue;
              } // if all integrals screened out, skip to next quartet

          for(auto f1=0, f1234=0, f123=0; f1!=n1; ++f1) {
            const auto bf1 = f1 + bf1_first;
            for(auto f3=0; f3!=n3; ++f3) {
                const auto bf3 = f3 + bf3_first;
                for(auto f2=0; f2!=n2; ++f2,++f123) {
                    const auto bf2 = f2 + bf2_first;
                    for(auto f4=0; f4!=n4; ++f4, ++f1234) {
                        const auto bf4 = f4 + bf4_first;
                        const auto value = buf_1234[f1234];
                        const auto value_scal_by_deg = value * s12_deg;//<13|24>


                        for(int l_x=0;l_x<n_x;l_x++){
                            x4[f123*n_x+l_x]+=X_VEC[bf4*n_x+l_x]*value_scal_by_deg;
                        }
                      
                    }
                  }
                }
              }

            }
            for(auto f1=0, f12=0, f123=0; f1!=n1; ++f1) {
                const auto bf1 = f1 + bf1_first;
                for(auto f3=0; f3!=n3; ++f3) {
                    const auto bf3 = f3 + bf3_first;
                    for(auto f2=0; f2!=n2; ++f2,++f123) {
                        const auto bf2 = f2 + bf2_first;
                        for(int k_x=0;k_x<n_x;k_x++)
                        for(int l_x=0;l_x<n_x;l_x++){
                            f12=f1*n2+f2;
                            xy34[(f12*n_x+k_x)*n_x+l_x]+=X_VEC[bf3*n_x+k_x]*x4[f123*n_x+l_x];
                        }

                    }
                }
            }
            
          }
          for(auto f1=0, f12=0; f1!=n1; ++f1) {
              const auto bf1 = f1 + bf1_first;
              for(auto f2=0; f2!=n2; ++f2,++f12) {
                  const auto bf2 = f2 + bf2_first;
                  for(int j_y=0;j_y<n_y;j_y++)
                  for(int k_x=0;k_x<n_x;k_x++)
                  for(int l_x=0;l_x<n_x;l_x++){
                      m_y[((f1*n_y+j_y)*n_x+k_x)*n_x+l_x]+=Y_VEC[bf2*n_y+j_y]*xy34[(f12*n_x+k_x)*n_x+l_x];
                      }

                  }
              }
          }
        for(auto f1=0; f1!=n1; ++f1) {
              const auto bf1 = f1 + bf1_first;
              for(int i_x=0;i_x<n_x;i_x++)
              for(int j_y=0;j_y<n_y;j_y++)
              for(int k_x=0;k_x<n_x;k_x++)
              for(int l_y=0;l_y<n_y;l_y++){
                  mo_ints[thread_id][((i_x*n_y+j_y)*n_x+k_x)*n_y+l_y]+=Y_VEC[bf1*n_y+j_y]*m_y[((f1*n_y+l_y)*n_x+i_x)*n_x+k_x]
                                                                         + Y_VEC[bf1*n_y+l_y]*m_y[((f1*n_y+j_y)*n_x+k_x)*n_x+i_x];
                  }
              }                
      }
      
//       delete[] e4  ;
      delete[] x4  ;
      delete[] xy34;
//       delete[] m_d ;
      delete[] m_y ;
  };

  libint2_parallel_do(lambda);
  //fprintf(stderr,"\r                                        \r");
  printf_timer("after XYXY transformation of AO_ints");

  double result = 0.0;

  for(int i=1;i<num_threads;i++)
  for(int j=0;j<n_x*n_y*n_x*n_y;j++)
      mo_ints[0][j]+=mo_ints[i][j];
  
  for(int i=1;i<num_threads;i++)delete[] mo_ints[i];
  
//   for(int i=0;i<num_threads;i++)delete[] Schwarz[i];
//   delete[] Schwarz;
//   for(int i=0;i<num_threads;i++) delete[] DM_m[i];
//   delete[] DM_m;
//   printf("exit 2el block\n");
  
  return result;
}

double calc_2el_MO_INTS_XXXY(std::vector<Shell> shells, int n_ao, double * XXXY_INTS, const double * __restrict__ X_VEC, const double * __restrict__ Y_VEC, int n_x, int n_y){

  //nthreads = nt;
  
  
  int nsh = shells.size();
//   double ** Schwarz;
//   Schwarz = new double *[num_threads];
//   for(int i=0;i<num_threads;i++) Schwarz[i] = new double[nsh*nsh];
//   set_zero_matr(Schwarz[0],nsh*nsh);
//   compute_schwarz_matr(shells, num_threads, Schwarz[0]);
//   for(int i=1;i<num_threads;i++)for(int j=0;j<nsh*nsh;j++) Schwarz[i][j]=Schwarz[0][j];
  
  std::vector<size_t> shell2bf[num_threads];
  shell2bf[0] = map_shell_to_basis_function(shells);
  for(int i=1;i<num_threads;i++)shell2bf[i]=shell2bf[0];
  
  std::vector<size_t> n_bf;
  n_bf = n_basis_function_in_shell(shells);
//   double ** DM_m;
//   DM_m = new double*[num_threads];
//   for(int i=0;i<num_threads;i++) DM_m[i] = new double[nsh*nsh];
  
//   gen_max_DM(DM_D,n_ao,DM_m[0],n_bf,shell2bf[0]);
  
//   for(int i=1;i<num_threads;i++)for(int j=0;j<nsh*nsh;j++) DM_m[i][j]=DM_m[0][j];
  
  double eps = 0;
  std::tie(obs_shellpair_list, obs_shellpair_data) = compute_shellpairs(shells,eps);
  
  double * mo_ints[num_threads];
  mo_ints[0] = XXXY_INTS;
  for(int i=1;i<num_threads;i++)mo_ints[i]  = new double[n_x*n_x*n_x*n_y];
  for(int i=0;i<num_threads;i++)set_zero_matr(mo_ints[i],n_x*n_x*n_x*n_y);
  double * J_th[num_threads];
  auto lambda = [&](const int thread_id) {

      double * __restrict__ x4;
      double * __restrict__ xx34;
      double * __restrict__ m_y;

      int mg = L_MAX*(L_MAX+1);
      int mg2 = mg*mg;
      int mg3 = mg2*mg;
      
      x4   = new double[n_x*        mg3];
      xx34 = new double[n_x*n_x*    mg2];
      m_y  = new double[n_x*n_x*n_x*mg ];

      int num_ints_computed = 0;

      const auto& buf = e_g_engines.at(thread_id).results();
      long s1234;
      int s1;
      for(s1=0, s1234=0; s1!=nsh; ++s1) {
        //fprintf(stderr,"AO[%7d/%7d]  threadid=%7d\r",s1,nsh,thread_id);

        if (s1 % num_threads != thread_id) continue;
        auto bf1_first = shell2bf[thread_id][s1]; // first basis function in this shell
        auto n1 = shells[s1].size();   // number of basis functions in this shell
        auto sp12_iter = obs_shellpair_data.at(s1).begin();
        set_zero_matr(m_y , n_x*n_x*n_x*mg);
        for (const auto s2 : obs_shellpair_list[s1]){
          
          set_zero_matr(xx34, n_x*n_x*mg2);
          
          auto bf2_first = shell2bf[thread_id][s2];
          auto n2 = shells[s2].size();
          for(auto s3=0; s3!=nsh; ++s3) {
            auto bf3_first = shell2bf[thread_id][s3];
            auto n3 = shells[s3].size();
            auto sp34_iter = obs_shellpair_data.at(s3).begin();

            const auto s4_max = (s1 == s3) ? s2 : s3;
            set_zero_matr(x4,n_x*mg3);
//             set_zero_matr(e4,n_y*mg3);
            for (const auto s4 : obs_shellpair_list[s3]) {
              if (s4 > s3)
                  break; 
//               if ((s1234++) % num_threads != thread_id) continue;
              auto bf4_first = shell2bf[thread_id][s4];
              auto n4 = shells[s4].size();
    
              // compute the permutational degeneracy (i.e. # of equivalents) of the given shell set
//               auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
              auto s34_deg = (s3 == s4) ? 0.5 : 1.0;
//               auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
//               auto s1234_deg = 1.0;//s12_deg * s34_deg * s12_34_deg/8;
              e_g_engines.at(thread_id).compute(shells[s1], shells[s2], shells[s3], shells[s4]);
              const auto* buf_1234 = buf[0]; //(12|34)
              
              
              if (buf_1234 == nullptr){
                  continue;
              } // if all integrals screened out, skip to next quartet
//               num_ints_computed += n1 * n2 * n3 * n4;


          for(auto f1=0, f1234=0, f123=0; f1!=n1; ++f1) {
            const auto bf1 = f1 + bf1_first;
            for(auto f2=0; f2!=n2; ++f2) {
                const auto bf2 = f2 + bf2_first;
                for(auto f3=0; f3!=n3; ++f3,++f123) {
                    const auto bf3 = f3 + bf3_first;
                    for(auto f4=0; f4!=n4; ++f4, ++f1234) {
                        const auto bf4 = f4 + bf4_first;
                        const auto value = buf_1234[f1234];
                        const auto value_scal_by_deg = value * s34_deg;//<13|24>

//                         for(int l_y=0;l_y<n_y;l_y++){
//                             e4[f123*n_y+l_y]+=Y_VEC[bf4*n_y+l_y]*value_scal_by_deg;
//                         }
                        for(int l_x=0;l_x<n_x;l_x++){
                            x4[f123*n_x+l_x]+=X_VEC[bf4*n_x+l_x]*value_scal_by_deg;
                        }
                      
                    }
                  }
                }
              }

            }
            for(auto f1=0, f12=0, f123=0; f1!=n1; ++f1) {
                const auto bf1 = f1 + bf1_first;
                for(auto f2=0; f2!=n2; ++f2,++f12) {
                    const auto bf2 = f2 + bf2_first;
                    for(auto f3=0; f3!=n3; ++f3,++f123) {
                        const auto bf3 = f3 + bf3_first;
                        for(int k_x=0;k_x<n_x;k_x++)
                        for(int l_x=0;l_x<n_x;l_x++)
                            xx34[(f12*n_x+k_x)*n_x+l_x]+=X_VEC[bf3*n_x+k_x]*x4[f123*n_x+l_x]
                                                         + X_VEC[bf3*n_x+l_x]*x4[f123*n_x+k_x];//!!!!!!!!!!!!!
                        
                        

                    }
                }
            }
            
          }
          for(auto f1=0, f12=0; f1!=n1; ++f1) {
              const auto bf1 = f1 + bf1_first;
              for(auto f2=0; f2!=n2; ++f2,++f12) {
                  const auto bf2 = f2 + bf2_first;
                  for(int j_x=0;j_x<n_x;j_x++)
                  for(int k_x=0;k_x<n_x;k_x++)
                  for(int l_x=0;l_x<n_x;l_x++){
                      m_y[((f1*n_x+j_x)*n_x+k_x)*n_x+l_x]+=X_VEC[bf2*n_x+j_x]*xx34[(f12*n_x+k_x)*n_x+l_x];
                      }


                  }
              }
          }
        for(auto f1=0; f1!=n1; ++f1) {
              const auto bf1 = f1 + bf1_first;
              for(int i_x=0;i_x<n_x;i_x++)
              for(int j_x=0;j_x<n_x;j_x++)
              for(int k_x=0;k_x<n_x;k_x++)
              for(int l_y=0;l_y<n_y;l_y++){
                  mo_ints[thread_id][((i_x*n_x+j_x)*n_x+k_x)*n_y+l_y]+=//X_VEC[bf1*n_x+i_x]*m_d1[((f1*n_x+j_x)*n_x+k_x)*n_y+l_y];
//                                                                          + X_VEC[bf1*n_x+j_x]*m_d1[((f1*n_x+i_x)*n_x+k_x)*n_y+l_y]
//                                                                          + X_VEC[bf1*n_x+k_x]*m_d2[((f1*n_y+l_y)*n_x+i_x)*n_x+j_x]
                                                                          Y_VEC[bf1*n_y+l_y]*m_y [((f1*n_x+k_x)*n_x+i_x)*n_x+j_x];
                  
                  }
              }                
      }
      
//       delete[] e4  ;
      delete[] x4  ;
//       delete[] xy34;
      delete[] xx34;
//       delete[] m_d1;
//       delete[] m_d2;
      delete[] m_y ;
  };

  libint2_parallel_do(lambda);
  //fprintf(stderr,"\r                                        \r");
  printf_timer("after AO_ints");

  for(int i=1;i<num_threads;i++)
  for(int j=0;j<n_x*n_x*n_x*n_y;j++)
      mo_ints[0][j]+=mo_ints[i][j];

  for(int i=1;i<num_threads;i++)delete[] mo_ints[i];


//   for(int i=0;i<num_threads;i++)delete[] Schwarz[i];
//   delete[] Schwarz;
//   for(int i=0;i<num_threads;i++) delete[] DM_m[i];
//   delete[] DM_m;
//   printf("exit 2el block\n");
  
  return 0;
}

int calc_2el_MP2(std::vector<Shell> shells, int n_ao, int n_mo, double * T, double * V, std::vector<libint2::Engine> * engines){

  
  int nsh = shells.size();
  double ** Schwarz;
  Schwarz = new double *[num_threads];
  for(int i=0;i<num_threads;i++) Schwarz[i] = new double[nsh*nsh];
  set_zero_matr(Schwarz[0],nsh*nsh);
  compute_schwarz_matr(shells, num_threads, Schwarz[0]);
  for(int i=1;i<num_threads;i++)for(int j=0;j<nsh*nsh;j++) Schwarz[i][j]=Schwarz[0][j];
  
  std::vector<size_t> shell2bf[num_threads];
  shell2bf[0] = map_shell_to_basis_function(shells);
  for(int i=1;i<num_threads;i++)shell2bf[i]=shell2bf[0];
  
  std::vector<size_t> n_bf;
  n_bf = n_basis_function_in_shell(shells);
//   double ** DM_m;
//   DM_m = new double*[num_threads];
//   for(int i=0;i<num_threads;i++) DM_m[i] = new double[nsh*nsh];
//   
//   gen_max_DM(DM_D,n_ao,DM_m[0],n_bf,shell2bf[0]);
//   
//   for(int i=1;i<num_threads;i++)for(int j=0;j<nsh*nsh;j++) DM_m[i][j]=DM_m[0][j];
  
  double eps = 0;
  std::tie(obs_shellpair_list, obs_shellpair_data) = compute_shellpairs(shells,eps);
// printf("screened\n");
  double * res;
  res = new double[num_threads];
  // loop over permutationally-unique set of shells
  
//   double * int_tr_tens[num_threads];
//   for(int i=0;i<num_threads;i++){
//       int_tr_tens[i] = new double[n_ao*n_ao*n_v*n_v];
//       set_zero_matr(int_tr_tens[i],n_ao*n_ao*n_v*n_v);
//   }
  
  
  auto lambda = [&](const int thread_id) {

      res[thread_id]=0;
//       double * int_tr_tens_n=int_tr_tens[thread_id];
      int num_ints_computed = 0;
//       auto engine = engines->at(thread_id);
      const auto& buf = engines->at(thread_id).results();
      long s1234;
      int s1;
      for(s1=0, s1234=0; s1!=nsh; ++s1) {
        //fprintf(stderr,"AO[%7d/%7d]  threadid=%7d\r",s1,nsh,thread_id);
//         fflush(stderr);
//         if (s1 % num_threads != thread_id) continue;
        auto bf1_first = shell2bf[thread_id][s1]; // first basis function in this shell
        auto n1 = shells[s1].size();   // number of basis functions in this shell
        auto sp12_iter = obs_shellpair_data.at(s1).begin();
        for (const auto s2 : obs_shellpair_list[s1]){
//           if (s2 > s1)
//                   break;  
//           const auto* sp12 = sp12_iter->get();
//           ++sp12_iter;
    
          auto bf2_first = shell2bf[thread_id][s2];
          auto n2 = shells[s2].size();
//           printf("n2 = %d\n",n2);
          for(auto s3=0; s3!=nsh/*s3<=s1*/; ++s3) {
//             printf("s3 = %d\n",s3);
            auto bf3_first = shell2bf[thread_id][s3];
            auto n3 = shells[s3].size();
            auto sp34_iter = obs_shellpair_data.at(s3).begin();
//             printf("n2 = %d\n",n2);
            const auto s4_max = (s1 == s3) ? s2 : s3;
            for (const auto s4 : obs_shellpair_list[s3]) {
//               if (s4 > s4_max)
//                   break;
//               printf("ijkl = %d, %d, %d,%d\n",s1,s2,s3,s4);
//               const auto* sp34 = sp34_iter->get();
//               ++sp34_iter;
//               printf("n2 = %d\n",n2);
    
              if ((s1234++) % num_threads != thread_id) continue;
              auto bf4_first = shell2bf[thread_id][s4];
              auto n4 = shells[s4].size();
    
              // compute the permutational degeneracy (i.e. # of equivalents) of the given shell set
//               auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
//               auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
//               auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
//               auto s1234_deg = s12_deg * s34_deg * s12_34_deg/8;
//               double DM_max = std::max(DM_m[thread_id][s1*nsh+s2],
//                                        DM_m[thread_id][s1*nsh+s3]);
//               DM_max = std::max(DM_max,DM_m[thread_id][s1*nsh+s4]);
//               DM_max = std::max(DM_max,DM_m[thread_id][s2*nsh+s1]);
//               DM_max = std::max(DM_max,DM_m[thread_id][s2*nsh+s3]);
//               DM_max = std::max(DM_max,DM_m[thread_id][s2*nsh+s4]);
//               DM_max = std::max(DM_max,DM_m[thread_id][s3*nsh+s1]);
//               DM_max = std::max(DM_max,DM_m[thread_id][s3*nsh+s2]);
//               DM_max = std::max(DM_max,DM_m[thread_id][s3*nsh+s4]);
//               DM_max = std::max(DM_max,DM_m[thread_id][s4*nsh+s1]);
//               DM_max = std::max(DM_max,DM_m[thread_id][s4*nsh+s2]);
//               DM_max = std::max(DM_max,DM_m[thread_id][s4*nsh+s3]);
            
    //           const auto tstart = std::chrono::high_resolution_clock::now();
    
//               if (DM_max*Schwarz[thread_id][s1*nsh + s2] * Schwarz[thread_id][s3*nsh + s4] < 1E-11){
//                   delete[]  TEMP_Tensor[thread_id];
//                   continue;
//               }
              engines->at(thread_id).compute(shells[s1], shells[s2], shells[s3], shells[s4]);
              const auto* buf_1234 = buf[0];
              
              
              if (buf_1234 == nullptr){
                  continue;
              } // if all integrals screened out, skip to next quartet
//               num_ints_computed += n1 * n2 * n3 * n4;


          for(auto f1=0, f1234=0; f1!=n1; ++f1) {
            const auto bf1 = f1 + bf1_first;
            for(auto f2=0; f2!=n2; ++f2) {
                const auto bf2 = f2 + bf2_first;
                for(auto f3=0; f3!=n3; ++f3) {
                const auto bf3 = f3 + bf3_first;
                    for(auto f4=0; f4!=n4; ++f4, ++f1234) {
                        const auto bf4 = f4 + bf4_first;
                        const auto value = buf_1234[f1234];
//                         const auto value_scal_by_deg = value * s1234_deg;//<13|24>
//                          printf("i = %e\n",value_scal_by_deg);
//                          if(fabs(value_scal_by_deg)>1e-10)
                        
                        int nklI=((bf3*n_ao+bf2)*n_ao+bf4)*n_mo;//n is bf3, k is bf2,l is bf4
//                         printf("v[%d,%d,%d,%d] = %e\n",bf1,bf2,bf3,bf4,value);
                        for(int i=0   ;i<n_mo;i++/*)
                        for(int j=0   ;j<n_mo;j++*/,nklI++){
                            T[nklI]+=value*V[i*n_ao+bf1];
    
                        }
                      
                    }
                  }
                }
              }

            }
          }
        }
      }
  };
  //fprintf(stderr,"\n");
//   libint2_parallel_do(lambda1);
//    printf_timer("after lambda");
  libint2_parallel_do(lambda);
// //   printf_timer("after lambda");
// //   printf("lambda done\n");
//   double result = 0.0;
//   for(int i=0;i<num_threads;i++)result+=res[i];
  
//   for(int i=0;i<num_threads;i++)
//   for(int k_vo=0;k_vo<n_v ;k_vo++)
//   for(int l_vo=0;l_vo<n_v ;l_vo++)
//   for(int c_ao=0;c_ao<n_ao;c_ao++)
//   for(int d_ao=0;d_ao<n_ao;d_ao++)
//             result+=int_tr_tens[i][((c_ao*n_ao+d_ao)*n_v+k_vo)*n_v+l_vo]*R_VEC[c_ao*n_v+k_vo]*R_VEC[d_ao*n_v+l_vo];
  

  for(int i=0;i<num_threads;i++)delete[] Schwarz[i];
  delete[] Schwarz;
//   for(int i=0;i<num_threads;i++) delete[] DM_m[i];
//   delete[] DM_m;
  
  
  return 0.0;
}

int RI_V_calc(std::vector<Shell> shells, int n_ao, double * V){
  
    int nsh = shells.size();
    
    std::vector<size_t> shell2bf[num_threads];
    shell2bf[0] = map_shell_to_basis_function(shells);
    for(int i=1;i<num_threads;i++)shell2bf[i]=shell2bf[0];
    
    std::vector<size_t> n_bf;
    n_bf = n_basis_function_in_shell(shells);
    
    double eps = 0;
    
    
#pragma omp parallel 
    {
        int thread_id = omp_get_thread_num();
        
        const auto& buf = r2g_engines[thread_id].results();
        long s12;
        int s1;
        for(s1=0, s12=0; s1!=nsh; ++s1) {
            //fprintf(stderr,"AO[%7d/%7d]  threadid=%7d\r",s1,nsh,thread_id);
            auto bf1_first = shell2bf[thread_id][s1]; // first basis function in this shell
            auto n1 = shells[s1].size();   // number of basis functions in this shell
            for (auto s2=0;s2!=nsh;++s2){
                if ((s12++) % num_threads != thread_id) continue;
                
                auto bf2_first = shell2bf[thread_id][s2];
                auto n2 = shells[s2].size();
                
                r2g_engines[thread_id].compute(shells[s1], shells[s2]);
                const auto* buf_1234 = buf[0];
                
                if (buf_1234 == nullptr){
                    continue;
                } // if all integrals screened out, skip to next
                
                for(auto f1=0, f12=0; f1!=n1; ++f1) {
                    const auto bf1 = f1 + bf1_first;
                    for(auto f2=0; f2!=n2; ++f2,++f12){
                        const auto bf2 = f2 + bf2_first;
                        const auto value = buf_1234[f12];
                        V[bf1*n_ao+bf2]=value;
                    }
                }
            }
        }
    }
    //fprintf(stderr,"\r                           \r");
//     for(int i=1;i<num_threads       ;i++)
//     for(int j=0;j<n_ao*n_ao;j++){
//         V[j]+=V_th[i][j];
//     }
//     for(int i=1;i<num_threads       ;i++)delete[]V_th[i];
//     delete[] V_th;
    
    return 0;
}

double calc_2el_MO_INTS_XY_RI(std::vector<Shell> shells, int n_ao, std::vector<Shell> aux_shells, int aux_n_ao,double * XY_RI_INTS, const double * __restrict__ X_VEC, const double * __restrict__ Y_VEC, int n_x, int n_y, std::vector<libint2::Engine> * engines){


  int nsh = shells.size();
  int aux_nsh = aux_shells.size();
//   double ** Schwarz;
//   Schwarz = new double *[num_threads];
//   for(int i=0;i<num_threads;i++) Schwarz[i] = new double[nsh*nsh];
//   set_zero_matr(Schwarz[0],nsh*nsh);
//   compute_schwarz_matr(shells, num_threads, Schwarz[0]);
//   for(int i=1;i<num_threads;i++)for(int j=0;j<nsh*nsh;j++) Schwarz[i][j]=Schwarz[0][j];
  
  std::vector<size_t> shell2bf[num_threads];
  shell2bf[0] = map_shell_to_basis_function(shells);
  for(int i=1;i<num_threads;i++)shell2bf[i]=shell2bf[0];
  
  std::vector<size_t> aux_shell2bf[num_threads];
  aux_shell2bf[0] = map_shell_to_basis_function(aux_shells);
  for(int i=1;i<num_threads;i++)aux_shell2bf[i]=aux_shell2bf[0];
  
  
//   std::vector<size_t> n_bf;
//   n_bf = n_basis_function_in_shell(shells);
//   double ** DM_m;
//   DM_m = new double*[num_threads];
//   for(int i=0;i<num_threads;i++) DM_m[i] = new double[nsh*nsh];
  
//   gen_max_DM(DM_D,n_ao,DM_m[0],n_bf,shell2bf[0]);
  
//   for(int i=1;i<num_threads;i++)for(int j=0;j<nsh*nsh;j++) DM_m[i][j]=DM_m[0][j];
  
//   double eps = 0;
//   std::tie(obs_shellpair_list, obs_shellpair_data) = compute_shellpairs(shells,eps);
  
  double * mo_ints[num_threads];
  mo_ints[0] = XY_RI_INTS;
  for(int i=1;i<num_threads;i++)mo_ints[i]  = new double[aux_n_ao*n_x*n_y];
  for(int i=0;i<num_threads;i++)set_zero_matr(mo_ints[i],aux_n_ao*n_x*n_y);
  
  auto lambda = [&](const int thread_id) {

      double * __restrict__ x3;
      double * __restrict__ xy23;
      
      int mg = L_MAX*(L_MAX+1);
      int mg2 = mg*mg;
//       int mg3 = mg2*mg;
      
      x3   = new double[n_x*        mg2];
      xy23 = new double[n_x*n_y*    mg];

      int num_ints_computed = 0;

      const auto& buf = engines->at(thread_id).results();
      long s1234;
      int s1;
      for(s1=0, s1234=0; s1!=aux_nsh; ++s1) {
        //fprintf(stderr,"RI[%7d/%7d]  threadid=%7d\r",s1,aux_nsh,thread_id);

        if (s1 % num_threads != thread_id) continue;
        auto bf1_first = aux_shell2bf[thread_id][s1]; // first basis function in this shell
        auto n1 = aux_shells[s1].size();   // number of basis functions in this shell
//         auto sp12_iter = obs_shellpair_data.at(s1).begin();
//         set_zero_matr(m_d1, n_y*n_x*n_x*mg);
//         set_zero_matr(m_d2, n_y*n_x*n_x*mg);
        set_zero_matr(xy23, n_x*n_y*mg);

        for (int s2=0; s2<nsh; s2++){
          
          set_zero_matr(x3,n_x*mg2);
          
          auto bf2_first = shell2bf[thread_id][s2];
          auto n2 = shells[s2].size();
          for(int s3=0; s3<nsh; s3++) {
              auto bf3_first = shell2bf[thread_id][s3];
              auto n3 = shells[s3].size();
//               auto sp34_iter = obs_shellpair_data.at(s3).begin();
              
//               const auto s4_max = (s1 == s3) ? s2 : s3;
//               set_zero_matr(e4,n_y*mg3);
              engines->at(thread_id).compute(aux_shells[s1], shells[s2], shells[s3]);
              const auto* buf_123 = buf[0]; //(1|23)
                
              if (buf_123 == nullptr){
                  continue;
              } // if all integrals screened out, skip to next quartet
              
              for(auto f1=0, f123=0, f12=0; f1!=n1; ++f1) {
                  const auto bf1 = f1 + bf1_first;
                  for(auto f2=0; f2!=n2; ++f2,++f12) {
                      const auto bf2 = f2 + bf2_first;
                      for(auto f3=0; f3!=n3; ++f3,++f123) {
                          const auto bf3 = f3 + bf3_first;
                          const auto value = buf_123[f123];
                          const auto value_scal_by_deg = value;//<13|24>
                          for(int l_x=0;l_x<n_x;l_x++){
                              x3[f12*n_x+l_x]+=X_VEC[bf3*n_x+l_x]*value_scal_by_deg;
                          }
                      }
                  }
              }
          }
          
          for(auto f1=0, f12=0; f1!=n1; ++f1) {
              const auto bf1 = f1 + bf1_first;
              for(auto f2=0; f2!=n2; ++f2,++f12) {
                  const auto bf2 = f2 + bf2_first;
                  for(int k_y=0;k_y<n_y;k_y++)
                  for(int l_x=0;l_x<n_x;l_x++){
                      xy23[(f1*n_y+k_y)*n_x+l_x]+=Y_VEC[bf2*n_y+k_y]*x3[f12*n_x+l_x];
//                       printf("%e\n", xy23[(f1*n_y+k_y)*n_x+l_x]);
//                   getchar();
                  }
              }
          }
        }
        
        for(auto f1=0; f1!=n1; ++f1) {
            const auto bf1 = f1 + bf1_first;
//             for(int i_xo=0;i_xo<aux_n_ao;i_xo++)
            for(int k_y=0;k_y<n_y;k_y++)
            for(int l_x=0;l_x<n_x;l_x++){
                mo_ints[thread_id][(bf1*n_y+k_y)*n_x+l_x]=xy23[(f1*n_y+k_y)*n_x+l_x];
//                 printf("%e\n", mo_ints[thread_id][(bf1*n_y+k_y)*n_x+l_x]);
//                 getchar();
            }
        }
      }
      delete[] x3  ;
      delete[] xy23;
  };

  libint2_parallel_do(lambda);
  //fprintf(stderr,"\r                                        \r");
  printf_timer("after AO_ints");

  for(int i=1;i<num_threads;i++)
  for(int j=0;j<aux_n_ao*n_x*n_y;j++)
      mo_ints[0][j]+=mo_ints[i][j];

  for(int i=1;i<num_threads;i++)delete[] mo_ints[i];


//   for(int i=0;i<num_threads;i++)delete[] Schwarz[i];
//   delete[] Schwarz;
//   for(int i=0;i<num_threads;i++) delete[] DM_m[i];
//   delete[] DM_m;
//   printf("exit 2el block\n");
  
  return 0;
}

int calc_2el_AO_INTS_RI(std::vector<Shell> shells, int n_ao, std::vector<Shell> aux_shells, int aux_n_ao, double * AO_RI_INTS){

  int nsh = shells.size();
  int aux_nsh = aux_shells.size();
  
  std::vector<size_t> shell2bf[num_threads];
  shell2bf[0] = map_shell_to_basis_function(shells);
  for(int i=1;i<num_threads;i++)shell2bf[i]=shell2bf[0];
  
  std::vector<size_t> aux_shell2bf[num_threads];
  aux_shell2bf[0] = map_shell_to_basis_function(aux_shells);
  for(int i=1;i<num_threads;i++)aux_shell2bf[i]=aux_shell2bf[0];
//   printf_timer("prepare");
//   double * ao_ints[num_threads];
//   ao_ints[0] = AO_RI_INTS;
//   for(int i=1;i<num_threads;i++)ao_ints[i]  = new double[aux_n_ao*n_ao*n_ao];
//   for(int i=0;i<num_threads;i++)set_zero_matr(ao_ints[i],aux_n_ao*n_ao*n_ao);
//   printf_timer("alloc");
  auto lambda = [&](const int thread_id) {

      int num_ints_computed = 0;

      const auto& buf = r_g_engines.at(thread_id).results();
      long s1234;
      int s1;
      for(s1=0, s1234=0; s1!=nsh; ++s1) {
        //fprintf(stderr,"RI[%7d/%7d]  threadid=%7d\r",s1,nsh,thread_id);

        if (s1 % num_threads != thread_id) continue;
        auto bf1_first = shell2bf[thread_id][s1]; // first basis function in this shell
        auto n1 = shells[s1].size();   // number of basis functions in this shell
//         auto sp12_iter = obs_shellpair_data.at(s1).begin();
        for (int s2=0; s2<nsh; s2++){
          
          auto bf2_first = shell2bf[thread_id][s2];
          auto n2 = shells[s2].size();
          for(int s3=0; s3<aux_nsh; s3++) {
              auto bf3_first = aux_shell2bf[thread_id][s3];
              auto n3 = aux_shells[s3].size();
//               auto sp34_iter = obs_shellpair_data.at(s3).begin();
              
//               const auto s4_max = (s1 == s3) ? s2 : s3;
//               set_zero_matr(e4,n_v*mg3);
              r_g_engines.at(thread_id).compute(aux_shells[s3], shells[s2], shells[s1]);
              const auto* buf_123 = buf[0]; //(1|23)
                
              if(buf_123 == nullptr)
                  for(auto f3=0, f123=0, f12=0; f3!=n3; ++f3) {
                      const auto bf3 = f3 + bf3_first;
                      for(auto f2=0; f2!=n2; ++f2,++f12) {
                          const auto bf2 = f2 + bf2_first;
                          for(auto f1=0; f1!=n1; ++f1,++f123) {
                              const auto bf1 = f1 + bf1_first;
                              AO_RI_INTS[(bf1*n_ao+bf2)*aux_n_ao+bf3] =0;
                              
                          }
                      }
                  }// if all integrals screened out, skip to next quartet
              else
                  for(auto f3=0, f123=0, f12=0; f3!=n3; ++f3) {
                      const auto bf3 = f3 + bf3_first;
                      for(auto f2=0; f2!=n2; ++f2,++f12) {
                          const auto bf2 = f2 + bf2_first;
                          for(auto f1=0; f1!=n1; ++f1,++f123) {
                              const auto bf1 = f1 + bf1_first;
                              const auto value = buf_123[f123];
//                               const auto value_scal_by_deg = value;//<13|24>
                              AO_RI_INTS[(bf1*n_ao+bf2)*aux_n_ao+bf3] =value/*_scal_by_deg*/;
                              
                          }
                      }
                  }
          }
        }
      }
  };
//   printf_timer("come to AO_ints");
  libint2_parallel_do(lambda);
  //fprintf(stderr,"\n                                        \r");
//   printf_timer("after AO_ints");

//   for(int i=1;i<num_threads;i++)
//   for(int j=0;j<aux_n_ao*n_ao*n_ao;j++)
//       ao_ints[0][j]+=ao_ints[i][j];

//   for(int i=1;i<num_threads;i++)delete[] ao_ints[i];
  
  
  return 0;
}

int calc_2el_AO_INTS_RI_1atom(std::vector<Shell> shells, int n_ao, std::vector<Shell> aux_shells, int aux_n_ao, double * AO_RI_INTS){

  
  int nsh = shells.size();
  int aux_nsh = aux_shells.size();
  
  std::vector<size_t> shell2bf[num_threads];
  shell2bf[0] = map_shell_to_basis_function(shells);
  for(int i=1;i<num_threads;i++)shell2bf[i]=shell2bf[0];
  
  std::vector<size_t> aux_shell2bf[num_threads];
  aux_shell2bf[0] = map_shell_to_basis_function(aux_shells);
  for(int i=1;i<num_threads;i++)aux_shell2bf[i]=aux_shell2bf[0];
//   printf_timer("prepare");
//   double * ao_ints[num_threads];
//   ao_ints[0] = AO_RI_INTS;
//   for(int i=1;i<num_threads;i++)ao_ints[i]  = new double[aux_n_ao*n_ao*n_ao];
  set_zero_matr(AO_RI_INTS,aux_n_ao*n_ao*n_ao);
//   printf_timer("alloc");
  auto lambda = [&](const int thread_id) {

      int num_ints_computed = 0;

      const auto& buf = r_g_engines.at(thread_id).results();
      long s1234;
      int s1;
      for(s1=0, s1234=0; s1!=nsh; ++s1) {
        //fprintf(stderr,"RI[%7d/%7d]  threadid=%7d\r",s1,nsh,thread_id);

        if (s1 % num_threads != thread_id) continue;
        auto bf1_first = shell2bf[thread_id][s1]; // first basis function in this shell
        auto n1 = shells[s1].size();   // number of basis functions in this shell
//         auto sp12_iter = obs_shellpair_data.at(s1).begin();
        for (int s2=0; s2<nsh; s2++){
          
          if((shells[s1].O[0]-shells[s2].O[0])*(shells[s1].O[0]-shells[s2].O[0])+
             (shells[s1].O[1]-shells[s2].O[1])*(shells[s1].O[1]-shells[s2].O[1])+
             (shells[s1].O[2]-shells[s2].O[2])*(shells[s1].O[2]-shells[s2].O[2])>1e-6)continue;
          auto bf2_first = shell2bf[thread_id][s2];
          auto n2 = shells[s2].size();
          for(int s3=0; s3<aux_nsh; s3++) {
              auto bf3_first = aux_shell2bf[thread_id][s3];
              auto n3 = aux_shells[s3].size();
//               auto sp34_iter = obs_shellpair_data.at(s3).begin();
              
//               const auto s4_max = (s1 == s3) ? s2 : s3;
//               set_zero_matr(e4,n_v*mg3);
              r_g_engines.at(thread_id).compute(aux_shells[s3], shells[s2], shells[s1]);
              const auto* buf_123 = buf[0]; //(1|23)
                
              if (buf_123 == nullptr){
                  continue;
              } // if all integrals screened out, skip to next quartet
              
              for(auto f3=0, f123=0, f12=0; f3!=n3; ++f3) {
                  const auto bf3 = f3 + bf3_first;
                  for(auto f2=0; f2!=n2; ++f2,++f12) {
                      const auto bf2 = f2 + bf2_first;
                      for(auto f1=0; f1!=n1; ++f1,++f123) {
                          const auto bf1 = f1 + bf1_first;
                          const auto value = buf_123[f123];
//                           const auto value_scal_by_deg = value;//<13|24>
                          AO_RI_INTS[(bf1*n_ao+bf2)*aux_n_ao+bf3] =value/*_scal_by_deg*/;
                          
                      }
                  }
              }
          }
        }
      }
  };
//   printf_timer("come to AO_ints");
  libint2_parallel_do(lambda);
  //fprintf(stderr,"\n                                        \r");
//   printf_timer("after AO_ints");

//   for(int i=1;i<num_threads;i++)
//   for(int j=0;j<aux_n_ao*n_ao*n_ao;j++)
//       ao_ints[0][j]+=ao_ints[i][j];

//   for(int i=1;i<num_threads;i++)delete[] ao_ints[i];
  
  
  return 0;
}

int calc_2el_AO_INTS_RI_aux_norm(std::vector<Shell> shells, int n_ao, std::vector<Shell> * aux_shells, int aux_n_ao){

  
  int nsh = shells.size();
  int aux_nsh = aux_shells[0].size();
  
  std::vector<size_t> shell2bf[num_threads];
  shell2bf[0] = map_shell_to_basis_function(shells);
  for(int i=1;i<num_threads;i++)shell2bf[i]=shell2bf[0];
  
  std::vector<size_t> aux_shell2bf[num_threads];
  aux_shell2bf[0] = map_shell_to_basis_function(aux_shells[0]);
  for(int i=1;i<num_threads;i++)aux_shell2bf[i]=aux_shell2bf[0];
//   printf_timer("prepare");
//   double * ao_ints[num_threads];
//   ao_ints[0] = AO_RI_INTS;
//   for(int i=1;i<num_threads;i++)ao_ints[i]  = new double[aux_n_ao*n_ao*n_ao];
//   printf_timer("alloc");
      double * V_max = new double[aux_nsh];
  {
      int thread_id=0;

      int num_ints_computed = 0;
      
      set_zero_matr(V_max,aux_nsh);

      const auto& buf = r_g_engines.at(thread_id).results();
      long s1234;
      int s1;
      for(s1=0, s1234=0; s1!=nsh; ++s1) {
        //fprintf(stderr,"RI[%7d/%7d]  threadid=%7d\r",s1,nsh,thread_id);

//         if (s1 % num_threads != thread_id) continue;
        auto bf1_first = shell2bf[thread_id][s1]; // first basis function in this shell
        auto n1 = shells[s1].size();   // number of basis functions in this shell
//         auto sp12_iter = obs_shellpair_data.at(s1).begin();
        for (int s2=0; s2<nsh; s2++){
          
//           if((shells[s1].O[0]-shells[s2].O[0])*(shells[s1].O[0]-shells[s2].O[0])+
//              (shells[s1].O[1]-shells[s2].O[1])*(shells[s1].O[1]-shells[s2].O[1])+
//              (shells[s1].O[2]-shells[s2].O[2])*(shells[s1].O[2]-shells[s2].O[2])>1e-6)continue;
          auto bf2_first = shell2bf[thread_id][s2];
          auto n2 = shells[s2].size();
          for(int s3=0; s3<aux_nsh; s3++) {
              auto bf3_first = aux_shell2bf[thread_id][s3];
              auto n3 = aux_shells[0][s3].size();
//               auto sp34_iter = obs_shellpair_data.at(s3).begin();
              
//               const auto s4_max = (s1 == s3) ? s2 : s3;
//               set_zero_matr(e4,n_v*mg3);
              r_g_engines.at(thread_id).compute(aux_shells[0][s3], shells[s2], shells[s1]);
              const auto* buf_123 = buf[0]; //(1|23)
                
              if (buf_123 == nullptr){
                  continue;
              } // if all integrals screened out, skip to next quartet
              
              for(auto f3=0, f123=0, f12=0; f3!=n3; ++f3) {
                  const auto bf3 = f3 + bf3_first;
                  for(auto f2=0; f2!=n2; ++f2,++f12) {
                      const auto bf2 = f2 + bf2_first;
                      for(auto f1=0; f1!=n1; ++f1,++f123) {
                          const auto bf1 = f1 + bf1_first;
                          const auto value = fabs(buf_123[f123]);
//                           const auto value_scal_by_deg = value;//<13|24>
                          if(V_max[s3]<value)V_max[s3]=value;
                          
                      }
                  }
              }
          }
        }
      }
  }
//   for(int s3=0; s3<aux_nsh; s3++)
//       printf("maxI(%d) =%e\n",s3,V_max[s3]);
//   getchar();
  for(int s3=0; s3<aux_nsh; s3++)
  for(int i=0;i<aux_shells[0][s3].contr.size();i++)
  for(int j=0;j<aux_shells[0][s3].contr[i].coeff.size();j++)
      aux_shells[0][s3].contr[i].coeff[j]=aux_shells[0][s3].contr[i].coeff[j]/V_max[s3];
  
  //fprintf(stderr,"\n                                        \r");
  
  
  return 0;
}


int calc_2el_AO_INTS_RI_aux_transform(std::vector<Shell> shells, int n_ao, std::vector<Shell> aux_shells, int aux_n_ao, double * __restrict__ T, double * __restrict__ AO_RI_INTS){

  
  int nsh = shells.size();
  int aux_nsh = aux_shells.size();
  
  std::vector<size_t> shell2bf[num_threads];
  shell2bf[0] = map_shell_to_basis_function(shells);
  for(int i=1;i<num_threads;i++)shell2bf[i]=shell2bf[0];
  
  std::vector<size_t> aux_shell2bf[num_threads];
  aux_shell2bf[0] = map_shell_to_basis_function(aux_shells);
  for(int i=1;i<num_threads;i++)aux_shell2bf[i]=aux_shell2bf[0];
  
  for(int i=0;i<num_threads+1;i++)
      nopt_parallel_index[i]=(i*nsh)/num_threads;

  
// #ifdef _OPENBLAS
//     int ntb = openblas_get_num_threads();
//     openblas_set_num_threads(1);
// #endif
  
  set_zero_matr(AO_RI_INTS,n_ao*n_ao*aux_n_ao);
#pragma omp parallel
  {
      int thread_id = omp_get_thread_num();
      
//       double * __restrict__ AA_ints = new double[aux_n_ao*50*50];
      
      const auto& buf = r_g_engines.at(thread_id).results();
      
      int start  =0; nopt_parallel_index[thread_id  ];
      int finish =nsh; nopt_parallel_index[thread_id+1];

      for(int s1=start; s1 < finish; s1++) {
        //fprintf(stderr,"RI[%7d/%7d]  threadid=%7d\r",s1,nsh,thread_id);

        if (s1 % num_threads != thread_id) continue;
        
        auto bf1_first = shell2bf[thread_id][s1]; // first basis function in this shell
        auto n1 = shells[s1].size();   // number of basis functions in this shell
//         auto sp12_iter = obs_shellpair_data.at(s1).begin();
        for (int s2=0; s2<nsh; s2++){
          
          auto bf2_first = shell2bf[thread_id][s2];
          auto n2 = shells[s2].size();
          for(int s3=0; s3<aux_nsh; s3++) {
              
              auto bf3_first = aux_shell2bf[thread_id][s3];
              auto n3 = aux_shells[s3].size();
              
              r_g_engines.at(thread_id).compute(aux_shells[s3], shells[s2], shells[s1]);
              
              const auto* buf_123 = buf[0]; //(1|23)
                
              if (buf_123 == nullptr){
//                   for(auto f3=0, f123=0, f12=0; f3!=n3; ++f3) {
//                       const auto bf3 = f3 + bf3_first;
//                       for(auto f2=0; f2!=n2; ++f2,++f12) {
//                           const auto bf2 = f2 + bf2_first;
//                           for(auto f1=0; f1!=n1; ++f1,++f123) {
//                               const auto bf1 = f1 + bf1_first;
//                               AA_ints[(f1*n2+f2)*aux_n_ao+bf3]=0;/*_scal_by_deg*/;
//                               
//                           }
//                       }
//                   }
                  
                  continue;
              } // if all integrals screened out, skip to next quartet
              
              for(auto f3=0, f123=0, f12=0; f3!=n3; ++f3) {
                  const auto bf3 = f3 + bf3_first;
                  for(auto f2=0; f2!=n2; ++f2,++f12) {
                      const auto bf2 = f2 + bf2_first;
                      for(auto f1=0; f1!=n1; ++f1,++f123) {
                          const auto bf1 = f1 + bf1_first;
                          const auto value = buf_123[f123];
//                           const auto value_scal_by_deg = value;//<13|24>
//                           AA_ints[(f1*n2+f2)*aux_n_ao+bf3]=value;/*_scal_by_deg*/;
                            for(int bf4=0;bf4<aux_n_ao;bf4++)
                                AO_RI_INTS[(bf1*n_ao+bf2)*aux_n_ao+bf4]+=value*T[bf3*aux_n_ao+bf4];/*_scal_by_deg*/;
                          
                      }
                  }
              }
          }
//           for(auto f2=0; f2!=n2; ++f2) {
//               const auto bf2 = f2 + bf2_first;
//               for(auto f1=0; f1!=n1; ++f1) {
//                   const auto bf1 = f1 + bf1_first;
//                   for(int bf3=0;bf3<aux_n_ao;bf3++)
//                   for(int bf4=0;bf4<aux_n_ao;bf4++)
//                       AO_RI_INTS[(bf1*n_ao+bf2)*aux_n_ao+bf4]+=AA_ints[(f1*n2+f2)*aux_n_ao+bf3]*T[bf3*aux_n_ao+bf4];/*_scal_by_deg*/;
//               }
//           }
        }
      }
  };
  
//  #ifdef _OPENBLAS
//     openblas_set_num_threads(ntb);
// #endif

  
  printf_timer("come to AO_ints");
  //fprintf(stderr,"\n                                        \r");
  printf_timer("after AO_ints");

 
  
  return 0;
}





