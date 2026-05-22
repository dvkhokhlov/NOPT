# include <vector>
# include <libint2/shell.h>
# include <libint2/engine.h>
# include <libint2.hpp>
# include <omp.h>
# include <chrono>
#include <unordered_map>
#include <mutex>

using libint2::Shell;
using libint2::Engine;
using libint2::Operator;

extern const double max_engine_precision;

size_t nbasis(const std::vector<libint2::Shell>& shells);


std::vector<size_t> map_shell_to_basis_function(const std::vector<libint2::Shell>& shells) ;

std::vector<size_t> n_basis_function_in_shell(const std::vector<libint2::Shell>& shells) ;

using shellpair_list_t = std::unordered_map<size_t, std::vector<size_t>>;
extern shellpair_list_t obs_shellpair_list;  // shellpair list for OBS
using shellpair_data_t = std::vector<std::vector<std::shared_ptr<libint2::ShellPair>>>;  // in same order as shellpair_list_t
extern shellpair_data_t obs_shellpair_data;

extern shellpair_list_t obs_shellpair_list1;  // shellpair list for OBS
extern shellpair_data_t obs_shellpair_data1;

extern shellpair_list_t obs_shellpair_list2;  // shellpair list for OBS
extern shellpair_data_t obs_shellpair_data2;


template <typename Lambda>
void libint2_parallel_do(Lambda& lambda) {
#ifdef _OPENMP
#pragma omp parallel 
  {
    auto thread_id = omp_get_thread_num();
    lambda(thread_id);
   
  }
#else  // use C++11 threads
  std::vector<std::thread> threads;
  for (int thread_id = 0; thread_id != libint2::nthreads; ++thread_id) {
    if (thread_id != nthreads - 1)
      threads.push_back(std::thread(lambda, thread_id));
    else
      lambda(thread_id);
  }  // threads_id
  for (int thread_id = 0; thread_id < nthreads - 1; ++thread_id)
    threads[thread_id].join();
#endif
}


int compute_schwarz_matr(std::vector<Shell> shells, int nt, double * K) ;


std::tuple<shellpair_list_t,shellpair_data_t>
compute_shellpairs(std::vector<Shell> shells,
                   const double threshold) ;

