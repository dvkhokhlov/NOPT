#pragma once
//
// Internal state + shared helpers for the block2 DMRG backend. Heavy block2 headers live
// here, so include this ONLY from the block2 backend TUs (block2_casci_wrap / block2_dmrg /
// block2_mps_to_det). The block2-free public interface stays in block2_casci_wrap.h.

#include <cstdint>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>
#include <omp.h>

#include "block2_casci_wrap.h"   // block2_casci_wrap class + dmrg_par (via inp_par_read.h)
#include "common_vars.h"         // out_stream
#include "blas_link.h"           // openblas_set_num_threads

#include "block2_core.hpp"
#include "block2_dmrg.hpp"

using namespace block2;

// -------------------------- engine: all block2 state ----------------------------------
struct dmrgci_engine {
    dmrg_par cfg;
    // twos = 2S (total spin) drives the spin-adapted SU2 solve; twosz = 2*M_S = na-nb is the
    // projection the determinant read-out expands to, so it matches the native (aldet) M_S.
    int n_act, n_elec, twos, twosz, mult, n_s, print_number;
    SU2 target;                       // active-space symmetry sector (C1: pg = 0)
    std::vector<uint8_t> orbsym;      // per-orbital irrep; C1 -> all 0 (set in set_act_rep_num)
    std::vector<double> E_states;     // per-state energies (returned by E_states_ptr)
    bool storage_ready = false;       // has init_state_storage run?

    // Localization. Empty/off => solve in the given basis.
    std::vector<double> U_loc;        // n_act x n_act, [a*n_act+p]; valid when localize_on
    bool localize_on = false;
    std::vector<double> F_loc, g_loc; // integrals rotated into the localized basis (when localize_on)

    // Active-block canonicalization (eigenvectors of the active Fock, [a*n_act+p]) supplied by the
    // host before the final print; rotates the reported leading determinants into the canonical basis.
    std::vector<double> U_canon;
    bool have_canon = false;

    // Warm-start (localization-rotation MPS reuse). R_active takes the previous macro-iteration's
    // active basis to the current one; valid only when have_rotation. Consumed by the next solve().
    std::vector<double> R_active;     // n_act x n_act, [a*n_act+p]
    bool have_rotation = false;

    // block2 objects for the current macro-iteration (rebuilt each import_integrals)
    std::shared_ptr<FCIDUMP<double>> fcidump;
    std::shared_ptr<HamiltonianQC<SU2, double>> hamil;
    std::shared_ptr<MPO<SU2, double>> mpo;
    std::shared_ptr<MultiMPSInfo<SU2>> mps_info;  // persists solve -> RDM read-out
    std::shared_ptr<MultiMPS<SU2, double>> mps;   // the converged (state-averaged) wavefunction
    std::vector<double> d2_cache;                 // per-state block2 2-RDM: n_s blocks of n_act^4
    bool d2_valid = false;                        // is d2_cache current for this solve?
    std::vector<double> dmfull_cache;             // full n_s x n_s spin-summed 1-RDM (properties), delocalized
    bool dmfull_valid = false;                     // is dmfull_cache current for this solve?
    int solve_count = 0;                        // macro-iteration index -> unique MPS tag
    int last_n_sweeps = 0;                      // sweeps actually run in the last solve
    double last_sweep_dE = 0.0;                 // |dE| between the final two sweeps (achieved convergence)
    bool last_hit_max = false;                  // last solve used its full sweep budget with dE > sweep_tol
    std::vector<uint16_t> reorder_perm;         // DMRG lattice order (Fiedler); empty => input order

    dmrgci_engine(int n_act_, int n_elec_, int twos_, int twosz_, int mult_, int n_s_,
                  int print_number_, const dmrg_par &c)
        : cfg(c), n_act(n_act_), n_elec(n_elec_), twos(twos_), twosz(twosz_), mult(mult_),
          n_s(n_s_), print_number(print_number_), target(n_elec_, twos_, 0), orbsym(n_act_, 0),
          E_states(n_s_, 0.0) {}
};

// -------------------------- shared block2-backend helpers ----------------------------------
namespace nopt_block2 {

// Pin block2 to the host OpenMP/BLAS thread count on entering a block2 region; restore on exit.
struct host_threads_guard {
    int omp_saved;
    int blas_saved;
    host_threads_guard() : omp_saved(omp_get_max_threads()), blas_saved(openblas_get_num_threads()) {}
    ~host_threads_guard() {
        omp_set_num_threads(omp_saved);
        openblas_set_num_threads(blas_saved); // restore the host's OpenBLAS count independently of OMP
    }
    host_threads_guard(const host_threads_guard &) = delete;
    host_threads_guard &operator=(const host_threads_guard &) = delete;
};

void ensure_block2_runtime(const std::string &save_dir_root, int n_threads);
void remove_tag_files(const std::string &tag);
void assert_stack_clean(const char *where);

// Fit a lower-bond-dim copy of an MPS (identity-MPO Linear) — a cheaper TRIE for the read-out.
// Defined engine-side; called from the read-out TU.
std::shared_ptr<MPS<SU2, double>>
compress_single_mps(dmrgci_engine &e, const std::shared_ptr<MPS<SU2, double>> &ket,
                    int target_m, const std::string &ctag);

} // namespace nopt_block2
