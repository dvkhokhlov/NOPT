// block2_casci_wrap — casci_solver backend over the external block2 DMRG library.

#include "block2_casci_wrap.h"

// Heavy block2 headers — confined to this cpp
#include "block2_core.hpp"
#include "block2_dmrg.hpp"

#include <cstdio>
#include <cstdlib>
#include <mutex>
#include <string>
#include <unistd.h>     // getpid
#include <omp.h>

#include "common_vars.h"   // out_stream

using namespace block2;

// ---- process-global block2 runtime ---------------------------------------------------
namespace {

std::once_flag g_b2_once;

void ensure_block2_runtime(const std::string &save_dir_root, int n_threads) {
    std::call_once(g_b2_once, [&] {
        // Per-process scratch subdir under the configured root (default /dev/shm), so
        // concurrent/!repeated runs never share block2's renormalized-operator files.
        const std::string scratch =
            save_dir_root + "/nopt_dmrg_" + std::to_string((long)getpid());

        Random::rand_seed(0);
        // isize/dsize are BYTE sizes of the integer/double stacks 
        frame_<double>() = std::make_shared<DataFrame<double>>(
            (size_t)1 << 24, (size_t)1 << 30, scratch);
        frame_<double>()->use_main_stack = false;
        frame_<double>()->minimal_disk_usage = true;
        frame_<double>()->minimal_memory_usage = false;

        threading_() = std::make_shared<Threading>(
            ThreadingTypes::OperatorBatchedGEMM | ThreadingTypes::Global,
            n_threads, n_threads, 1);
        threading_()->seq_type = SeqTypes::Tasked;

        fprintf(out_stream,
                "block2 runtime initialized (scratch=%s, threads=%d)\n",
                scratch.c_str(), n_threads);
    });
}

} // namespace

// -------------------------- engine: all block2 state ----------------------------------
struct dmrgci_engine {
    dmrg_par cfg;
    explicit dmrgci_engine(const dmrg_par &c) : cfg(c) {}
};

// ---- P2.0 not-implemented sentinel ---------------------------------------------------
[[noreturn]] static void nyi(const char *what) {
    fprintf(out_stream,
            "ERROR: block2_casci_wrap::%s not implemented yet (P2.0 skeleton)\n", what);
    exit(0);
}

// ---- ctor / dtor ---------------------------------------------------------------------
block2_casci_wrap::block2_casci_wrap(const dmrg_par &cfg)
    : impl_(std::make_unique<dmrgci_engine>(cfg)) {
    int nthr = omp_get_max_threads();
    if (nthr < 1) nthr = 1;
    ensure_block2_runtime(cfg.save_dir, nthr);
}

block2_casci_wrap::~block2_casci_wrap() = default;

// ----------------------- casci_solver contract ----------------------------------------
void block2_casci_wrap::init_state_storage(int, int) { nyi("init_state_storage"); }
bool block2_casci_wrap::has_coef(int) const          { nyi("has_coef"); }
void block2_casci_wrap::set_act_rep_num(int *)        { nyi("set_act_rep_num"); }
void block2_casci_wrap::import_integrals(double *, double *, double) { nyi("import_integrals"); }
int  block2_casci_wrap::solve(int, int, bool)         { nyi("solve"); }
void block2_casci_wrap::calc_DM_diag(double *, int)   { nyi("calc_DM_diag"); }
void block2_casci_wrap::G_calc(double *)              { nyi("G_calc"); }
void block2_casci_wrap::calc_DMA(double *, int, int)  { nyi("calc_DMA"); }
void block2_casci_wrap::calc_DMB(double *, int, int)  { nyi("calc_DMB"); }
int    block2_casci_wrap::n_states()      const { nyi("n_states"); }
int    block2_casci_wrap::mult()          const { nyi("mult"); }
double block2_casci_wrap::E_state(int)    const { nyi("E_state"); }
double block2_casci_wrap::S2_state(int)   const { nyi("S2_state"); }
double block2_casci_wrap::L2_state(int)   const { nyi("L2_state"); }
double block2_casci_wrap::P_state(int)    const { nyi("P_state"); }
double *block2_casci_wrap::E_states_ptr() const { nyi("E_states_ptr"); }
void block2_casci_wrap::gen_ext_ind()                 { nyi("gen_ext_ind"); }
void block2_casci_wrap::print_states(int, int, int)   { nyi("print_states"); }
void block2_casci_wrap::write_civec(int, char *)      { nyi("write_civec"); }
