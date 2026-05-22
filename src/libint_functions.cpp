# include <vector>
# include <libint2/shell.h>
# include <libint2/engine.h>
# include <libint2.hpp>
# include <omp.h>
# include <chrono>
#include <unordered_map>
#include <mutex>

# include "common_vars.h"
# include "matr.h"
# include "libint_functions.h"

using libint2::Shell;
using libint2::Engine;
using libint2::Operator;

const double max_engine_precision = std::numeric_limits<double>::epsilon() / 1e10;

size_t nbasis(const std::vector<libint2::Shell>& shells) {
  size_t n = 0;
  for (const auto& shell: shells)
    n += shell.size();
  return n;
}


std::vector<size_t> map_shell_to_basis_function(const std::vector<libint2::Shell>& shells) {
  std::vector<size_t> result;
//   result.reserve(shells.size());

  size_t n = 0;
  for (auto shell: shells) {
    result.push_back(n);
    n += shell.size();
  }

  return result;
}

std::vector<size_t> n_basis_function_in_shell(const std::vector<libint2::Shell>& shells) {
  std::vector<size_t> result;
//   result.reserve(shells.size());

  size_t n = 0;
  for (auto shell: shells) {
    n = shell.size();
    result.push_back(n);
  }

  return result;
}

int nthreads;

/// fires off \c nthreads instances of lambda in parallel
// template <typename Lambda>


int compute_schwarz_matr(std::vector<Shell> shells, int nt, double * K) {
  int nsh1 = shells.size();
  
  // construct the 2-electron repulsion integrals engine
  using libint2::Engine;
  
  nthreads = nt;
  std::vector<Engine> engines(nthreads);

  // !!! very important: cannot screen primitives in Schwarz computation !!!
  //auto epsilon = 0.;
  engines[0] = Engine(Operator::coulomb,GTO_MAX,L_MAX);
    
  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }

  
  auto compute = [&](int thread_id) {

    const auto& buf = engines[thread_id].results();

    // loop over permutationally-unique set of shells
    for (auto s1 = 0l, s12 = 0l; s1 != nsh1; ++s1) {
      
      auto n1 = shells[s1].size();  // number of basis functions in this shell

//       auto s2_max = shells1_equiv_shells2 ? s1 : nsh2 - 1;
      for (auto s2 = 0; s2 <= s1; ++s2, ++s12) {
        if (s12 % nthreads != thread_id) continue;

        auto n2 = shells[s2].size();
        auto n12 = n1 * n2;

        engines[thread_id].compute(shells[s1], shells[s2],shells[s1], shells[s2]);
//         assert(buf[0] != nullptr &&
//                "to compute Schwarz ints turn off primitive screening");
        if(buf[0]==nullptr){
            K[s1*nsh1 + s2] = 0.0;
            K[s2*nsh1 + s1] = 0.0;
        }
        // to apply Schwarz inequality to individual integrals must use the diagonal elements
        // to apply it to sets of functions (e.g. shells) use the whole shell-set of ints here
//         Eigen::Map<const Matrix> buf_mat(buf[0], n12, n12);
        else{
            double norm2 = find_max_abs_value(buf[0], n12, 1, n12);
            K[s1*nsh1 + s2] = std::sqrt(norm2);
            K[s2*nsh1 + s1] = K[s1*nsh1 + s2];//if (shells_equiv_shells) K(s2, s1) = K(s1, s2);
        }
      }
    }
  };  // thread lambda

  libint2_parallel_do(compute);

//   timer.stop(0);
//   std::cout << "done (" << timer.read(0) << " s)" << std::endl;

  return 0;
}

// using shellpair_list_t = std::unordered_map<size_t, std::vector<size_t>>;
shellpair_list_t obs_shellpair_list;  // shellpair list for OBS
// using shellpair_data_t = std::vector<std::vector<std::shared_ptr<libint2::ShellPair>>>;  // in same order as shellpair_list_t
shellpair_data_t obs_shellpair_data;  // shellpair data for OBS

shellpair_list_t obs_shellpair_list1;  // shellpair list for OBS
shellpair_data_t obs_shellpair_data1;  // shellpair data for OBS
shellpair_list_t obs_shellpair_list2;  // shellpair list for OBS
shellpair_data_t obs_shellpair_data2;  // shellpair data for OBS


std::vector<Engine> s_h_engines;

std::tuple<shellpair_list_t,shellpair_data_t>
compute_shellpairs(std::vector<Shell> shells,
                   const double threshold) {
  const auto nsh1 = shells.size();
//   printf("n_ao_sh = %d\n",nsh1);
  

  // construct the 2-electron repulsion(?????????) integrals engine
  
  shellpair_list_t splist;

  std::mutex mx;

  auto compute = [&](int thread_id) {

    auto& engine = s_h_engines[thread_id];
    const auto& buf = engine.results();

    // loop over permutationally-unique set of shells
    for (auto s1 = 0l, s12 = 0l; s1 != nsh1; ++s1) {
//       printf("start s1=%d\n",s1);
      mx.lock();
      if (splist.find(s1) == splist.end())
        splist.insert(std::make_pair(s1, std::vector<size_t>()));
      mx.unlock();
      

      auto n1 = shells[s1].size();  // number of basis functions in this shell

//       auto s2_max = bs1_equiv_bs2 ? s1 : nsh2 - 1;
      for (auto s2 = 0; s2 != nsh1; ++s2, ++s12) {
//         printf("start s2=%d\n",s2);
        if (s12 % num_threads != thread_id) continue;
//         printf("s1 =%d, s2 = %d\n",s1,s2);
    
        auto on_same_center = (shells[s1].O == shells[s2].O);
//         std::cout << (shells[s1].O == shells[s2].O) << std::endl; getchar();
        bool significant = on_same_center;
        if (not on_same_center) {
          auto n2 = shells[s2].size();
          s_h_engines[thread_id].compute(shells[s1], shells[s2]);
//           Eigen::Map<const Matrix> buf_mat(buf[0], n1, n2);
          auto norm = find_max_abs_value(buf[0], n1*n2, 1, 1);
//           printf("m = %e\n",norm);
//           norm = sqrt(cblas_ddot(n1*n2,buf[0],1,buf[0],1));
//           printf("l = %e\n",norm);
//           getchar();
          significant = (norm >= threshold);
        }

        if (significant) {
          mx.lock();
          splist[s1].emplace_back(s2);
          mx.unlock();
        }
      }
    }
  };  // end of compute

  libint2_parallel_do(compute);

  // resort shell list in increasing order, i.e. splist[s][s1] < splist[s][s2] if s1 < s2
  // N.B. only parallelized over 1 shell index
  auto sort = [&](int thread_id) {
    for (auto s1 = 0l; s1 != nsh1; ++s1) {
      if (s1 % num_threads == thread_id) {
        auto& list = splist[s1];
        std::sort(list.begin(), list.end());
      }
    }
  };  // end of sort

  libint2_parallel_do(sort);

  // compute shellpair data assuming that we are computing to default_epsilon
  // N.B. only parallelized over 1 shell index
  const auto ln_max_engine_precision = std::log(max_engine_precision);
  shellpair_data_t spdata(splist.size());
  auto make_spdata = [&](int thread_id) {
    for (auto s1 = 0l; s1 != nsh1; ++s1) {
      if (s1 % num_threads == thread_id) {
        for(const auto& s2 : splist[s1]) {
          spdata[s1].emplace_back(std::make_shared<libint2::ShellPair>(shells[s1],shells[s2],ln_max_engine_precision));
        }
      }
    }
  };  // end of make_spdata

  libint2_parallel_do(make_spdata);

//   timer.stop(0);
//   std::cout << "done (" << timer.read(0) << " s)" << std::endl;

  return std::make_tuple(splist,spdata);
}


