// block2_casci_wrap — casci_solver backend over the external block2 DMRG library.

#include "block2_casci_wrap.h"

// Heavy block2 headers — confined to this cpp
#include "block2_core.hpp"
#include "block2_dmrg.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <map>
#include <numeric>
#include <string>
#include <system_error>
#include <vector>
#include <unistd.h>     // getpid
#include <omp.h>

#include "common_vars.h"      // out_stream
#include "localized_dmrg.h"   // rotate1/rotate2 (active-space basis transforms)
#include "blas_link.h"        // cblas_dgemm (warm-start rotation regularization)
#include "mps_rotation.h"     // evolve_sa_multimps (multi-root SA MPS rotation)

using namespace block2;

// ---- process-global block2 runtime ---------------------------------------------------
namespace {

void report_dmrg_orbital_map(std::FILE *out, int n_act,
                             const std::vector<int> &canon_of_lattice,
                             const std::vector<double> &lattice_sign,
                             const char *basis_label) {
    std::fprintf(out, "DMRG determinant read-out orbital map (lattice k -> %s p, sign):\n", basis_label);
    std::fprintf(out, "  %-9s:", basis_label);
    for (int k = 0; k < n_act; k++) std::fprintf(out, " %d", canon_of_lattice[k]);
    std::fprintf(out, "\n  sign     :");
    for (int k = 0; k < n_act; k++) std::fprintf(out, " %c", lattice_sign[k] < 0.0 ? '-' : '+');
    std::fprintf(out, "\n");
}

// Process-global block2 runtime: owns the frame_/threading_ globals and the per-process scratch
// dir, and tears all of it down (RAII) at program exit.
struct Block2Runtime {
    std::string scratch;

    Block2Runtime(const std::string &save_dir_root, int n_threads) {
        // Per-process scratch subdir under the configured root (default /dev/shm), so
        // concurrent/repeated runs never share block2's renormalized-operator files.
        scratch = save_dir_root + "/nopt_dmrg_" + std::to_string((long)getpid());

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
    }

    ~Block2Runtime() {
        // Release block2's hold on the DataFrame, then remove the scratch
        if (frame_<double>() != nullptr) {
            frame_<double>()->activate(0);
            frame_<double>() = nullptr;
        }
        threading_() = nullptr;
        std::error_code ec;
        std::filesystem::remove_all(scratch, ec); // best-effort; never throw from a dtor
    }

    Block2Runtime(const Block2Runtime &) = delete;
    Block2Runtime &operator=(const Block2Runtime &) = delete;
};

// Build the runtime once, on first use; its destructor runs at program exit (see above).
void ensure_block2_runtime(const std::string &save_dir_root, int n_threads) {
    static Block2Runtime runtime(save_dir_root, n_threads);
    (void)runtime;
}

// Best-effort removal of all scratch files belonging to one MPS tag. Correctness does not depend
// on this — each solve writes a fresh unique tag, so a stale file can never be read back — but it
// keeps the per-process scratch dir from growing across macro-iterations. block2 embeds the tag,
// always trailing-'.'-delimited, in every per-tag filename:
//   <prefix>.MPS.<tag>.<i>        single-MPS site tensors       (mps.hpp)
//   <prefix>.MMPS.<tag>.<i>       MultiMPS site tensors + data  (state_averaged.hpp; SA solve/extract)
//   <prefix>.MMPS-WFN.<tag>.<i>   MultiMPS per-root wavefunctions
//   <prefix>.MPS.INFO.<tag>.*     (Multi)MPSInfo bond dims      (MultiMPSInfo inherits MPSInfo)
//   <tag>-mps_info.bin            serialized (Multi)MPSInfo
// The trailing '.' keeps tag "work_1" from matching "work_10"; ".MMPS." never prefix-matches ".MPS.".
void remove_tag_files(const std::string &tag) {
    auto fr = frame_<double>();
    if (fr == nullptr)
        return;
    std::error_code ec;
    const std::string t_mps  = fr->prefix + ".MPS." + tag + ".";
    const std::string t_mmps = fr->prefix + ".MMPS." + tag + ".";
    const std::string t_wfn  = fr->prefix + ".MMPS-WFN." + tag + ".";
    const std::string t_info = fr->prefix + ".MPS.INFO." + tag + ".";
    const std::string t_bin  = tag + "-mps_info.bin";
    std::vector<std::filesystem::path> victims;
    for (auto &de : std::filesystem::directory_iterator(fr->mps_dir, ec)) {
        const std::string name = de.path().filename().string();
        if (name.rfind(t_mps, 0) == 0 || name.rfind(t_mmps, 0) == 0 ||
            name.rfind(t_wfn, 0) == 0 || name.rfind(t_info, 0) == 0 || name == t_bin)
            victims.push_back(de.path());
    }
    for (auto &p : victims)
        std::filesystem::remove(p, ec);
}

// save NOPT OpenMP/BLAS thread counts on entering a block2 region and restore them on exit
struct host_threads_guard {
    int omp_saved;
    host_threads_guard() : omp_saved(omp_get_max_threads()) {}
    ~host_threads_guard() {
        omp_set_num_threads(omp_saved);
        openblas_set_num_threads(omp_saved);
    }
    host_threads_guard(const host_threads_guard &) = delete;
    host_threads_guard &operator=(const host_threads_guard &) = delete;
};

// Invariant guard: block2's LIFO int/double stacks must be empty at each macro-iteration boundary.
void assert_stack_clean(const char *where) {
    const size_t du = dalloc_<double>()->used, iu = ialloc_()->used;
    if (du != 0 || iu != 0) {
        fprintf(out_stream, "ERROR: block2 stack leak at %s (dalloc=%zu ialloc=%zu)\n", where, du,
                iu);
        assert(false && "block2 LIFO stack leak");
    }
}

// Build + simplify the conventional quantum-chemistry MPO 
std::shared_ptr<MPO<SU2, double>>
build_qc_mpo(const std::shared_ptr<HamiltonianQC<SU2, double>> &hamil) {
    std::shared_ptr<MPO<SU2, double>> mpo = std::make_shared<MPOQC<SU2, double>>(
        hamil, QCTypes::Conventional, "HQC", hamil->n_sites / 2 / 2 * 2, 1);
    mpo->basis = hamil->basis;
    mpo = std::make_shared<SimplifiedMPO<SU2, double>>(
        mpo, std::make_shared<RuleQC<SU2, double>>(), true, true,
        OpNamesSet({OpNames::R, OpNames::RD}));
    return mpo;
}

// Per-sweep DMRG schedule, expanded from breakpoints.
struct dmrg_schedule {
    std::vector<ubond_t> bond_dims;
    std::vector<double> noises;
    std::vector<double> dav_thrds;
    int n_sweeps = 0;
};

// Mirror of block2's "default" schedule (pyblock2/driver/parser.py::get_schedule):
// ramp startM(=min(250,maxM)) -> maxM over the def_m/def_iter/def_noise/def_tol tables,
// then a final noise-free convergence stage at davidson tol = sweep_tol/10. `user_sweeps`
// (= $DMRG sweeps) is a hard cap on the total; sweep_tol is the solve() early-stop.
dmrg_schedule build_default_schedule(int max_m, int user_sweeps, double sweep_tol) {
    struct point { int start, M; double tol, noise; };
    const int start_m = std::min(250, max_m);
    const std::vector<int> def_m = {50,  100,  250,  500,  1000, 2000, 3000,
                                    4000, 5000, 6000, 7000, 8000, 9000, 10000};
    const std::vector<int> def_it = {8, 8, 8, 8, 8, 4, 4, 4, 4, 4, 4, 4, 4, 4};
    const std::vector<double> def_noise = {1e-3, 1e-3, 1e-3, 1e-4, 1e-4, 5e-5, 5e-5,
                                           5e-5, 5e-5, 5e-5, 5e-5, 5e-5, 5e-5, 5e-5};
    const std::vector<double> def_tol = {1e-4, 1e-4, 1e-4, 1e-5, 1e-5, 5e-6, 5e-6,
                                         5e-6, 5e-6, 5e-6, 5e-6, 5e-6, 5e-6, 5e-6};
    const double final_tol = (sweep_tol <= 0 ? 1e-9 : sweep_tol / 10.0);

    std::vector<point> bp;
    if (start_m == max_m) { // small active space: M fixed at maxM (covers maxM <= 250)
        bp.push_back({0, max_m, 1e-5, 1e-4});
        bp.push_back({8, max_m, 5e-6, 5e-5});
    } else { // ramp 250 -> maxM
        int isweep = 0;
        for (size_t i = 0; i < def_m.size(); i++) {
            if (def_m[i] >= max_m) {
                bp.push_back({isweep, max_m, def_tol[i], def_noise[i]});
                break;
            } else if (def_m[i] >= start_m) {
                bp.push_back({isweep, def_m[i], def_tol[i], def_noise[i]});
                isweep += def_it[i];
            }
        }
    }
    bp.push_back({bp.back().start + 8, max_m, final_tol, 0.0}); // noise-free convergence
    // `sweeps` is the total sweep budget: it both caps a short run and *extends* the noise-free
    // final stage so a hard system can keep sweeping (early-stopping at sweep_tol) past the
    // schedule's nominal length. With sweeps unset it falls back to 4 noise-free sweeps.
    int maxiter = (user_sweeps > 0) ? user_sweeps : (bp.back().start + 4);
    if (maxiter < 1)
        maxiter = 1;

    dmrg_schedule s;
    s.n_sweeps = maxiter;
    s.bond_dims.resize(maxiter);
    s.noises.resize(maxiter);
    s.dav_thrds.resize(maxiter);
    for (int sw = 0; sw < maxiter; sw++) {
        const point *p = &bp.front();
        for (const auto &q : bp)
            if (q.start <= sw)
                p = &q;
            else
                break;
        s.bond_dims[sw] = (ubond_t)p->M;
        s.noises[sw] = p->noise;
        s.dav_thrds[sw] = p->tol;
    }
    return s;
}

// Warm re-solve schedule: no cold ramp (the rotated MPS is already at full M) -- a short noise-free
// run at the target bond dim. The exact rotation gives a near-perfect guess, so no perturbative noise
// is needed to re-expand the bond space. All-zero noise also satisfies block2's convergence rule
// (it declares convergence only once the sweep noise equals the final noise; sweep_algorithm.hpp).
dmrg_schedule build_warm_schedule(int max_m, int warm_sweeps, double sweep_tol) {
    int nsw = warm_sweeps > 0 ? warm_sweeps : 4;
    const double dav_final = (sweep_tol <= 0 ? 1e-9 : sweep_tol / 10.0);
    dmrg_schedule s;
    s.n_sweeps = nsw;
    s.bond_dims.assign(nsw, (ubond_t)max_m);
    s.noises.assign(nsw, 0.0);
    s.dav_thrds.assign(nsw, dav_final);
    return s;
}

} // namespace

// -------------------------- engine: all block2 state ----------------------------------
struct dmrgci_engine {
    dmrg_par cfg;
    int n_act, n_elec, twos, mult, n_s, print_number;
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

    dmrgci_engine(int n_act_, int n_elec_, int twos_, int mult_, int n_s_, int print_number_,
                  const dmrg_par &c)
        : cfg(c), n_act(n_act_), n_elec(n_elec_), twos(twos_), mult(mult_), n_s(n_s_),
          print_number(print_number_), target(n_elec_, twos_, 0), orbsym(n_act_, 0),
          E_states(n_s_, 0.0) {}
};

// ---- P2.0 not-implemented sentinel ---------------------------------------------------
[[noreturn]] static void nyi(const char *what) {
    fprintf(out_stream,
            "ERROR: block2_casci_wrap::%s not implemented yet (P2.0 skeleton)\n", what);
    exit(0);
}

// Compute every state's spatial 2-RDM once per solve and cache them as n_s consecutive block2
// D2[p,q,r,s] blocks. Each root is extracted from the state-averaged MultiMPS to its own single-root
// MPS, then one Expect sweep reads that root's 2-RDM (1-RDM is a partial trace of it, so CASSCF
// needs only this one sweep per root). A single sweep per root means the stale-`forward` landmine
// (consecutive sweeps flipping the center) never arises.
static void ensure_2rdm(dmrgci_engine &e) {
    if (e.d2_valid)
        return;
    host_threads_guard htg;
    const int n = e.n_act;
    const size_t blk = (size_t)n * n * n * n;
    e.d2_cache.assign(blk * e.n_s, 0.0);
    std::vector<double> bt_scratch;
    if (e.localize_on) bt_scratch.resize(blk); // back-transform target (rotate2 forbids aliasing)
    std::vector<double> ro_scratch;            // un-permute target (Fiedler lattice order -> input)
    std::vector<int> iperm;                    // inverse of reorder_perm
    if (!e.reorder_perm.empty()) {
        ro_scratch.resize(blk);
        iperm.resize(n);
        for (int i = 0; i < n; i++) iperm[e.reorder_perm[i]] = i;
    }

    for (int st = 0; st < e.n_s; st++) {
        const std::string xtag = e.mps_info->tag + "-" + std::to_string(st);
        std::shared_ptr<MultiMPS<SU2, double>> imps = e.mps->extract(st, xtag);

        std::shared_ptr<MPO<SU2, double>> p2mpo =
            std::make_shared<PDM2MPOQC<SU2, double>>(e.hamil);
        p2mpo = std::make_shared<SimplifiedMPO<SU2, double>>(
            p2mpo, std::make_shared<RuleQC<SU2, double>>(), true, true,
            OpNamesSet({OpNames::R, OpNames::RD}));
        auto p2me = std::make_shared<MovingEnvironment<SU2, double, double>>(p2mpo, imps, imps,
                                                                             "2PDM");
        p2me->init_environments(false);
        auto ex2 = std::make_shared<Expect<SU2, double, double>>(p2me, (ubond_t)e.cfg.m,
                                                                 (ubond_t)e.cfg.m);
        ex2->iprint = 0; // silence the per-site Expect sweep log
        ex2->solve(true, imps->center == 0);
        std::shared_ptr<GTensor<double>> d2 = ex2->get_2pdm_spatial(); // shape {n,n,n,n}

        // GTensor data is contiguous row-major [p,q,r,s] = the block2 D2 layout we cache; copy the
        // flat buffer straight into this state's block (no per-element accessor, no index leak).
        std::copy(d2->data->begin(), d2->data->end(), e.d2_cache.begin() + (size_t)st * blk);

        if (!e.reorder_perm.empty()) { // map the 2-RDM out of block2's Fiedler lattice order
            double *blkp = e.d2_cache.data() + (size_t)st * blk;
            for (int p = 0; p < n; p++)
                for (int q = 0; q < n; q++)
                    for (int r = 0; r < n; r++)
                        for (int s = 0; s < n; s++)
                            ro_scratch[(((size_t)p * n + q) * n + r) * n + s] =
                                blkp[(((size_t)iperm[p] * n + iperm[q]) * n + iperm[r]) * n +
                                     iperm[s]];
            std::copy(ro_scratch.begin(), ro_scratch.end(), blkp);
        }

        if (e.localize_on) { // rotate this state's 2-RDM back to the delocalized basis
            double *blkptr = e.d2_cache.data() + (size_t)st * blk;
            rotate2(blkptr, e.U_loc.data(), n, bt_scratch.data(), /*forward=*/false);
            std::copy(bt_scratch.begin(), bt_scratch.end(), blkptr);
        }

        p2me->remove_partition_files();
        p2mpo->deallocate();
        remove_tag_files(xtag); // the per-root extract is transient
    }
    e.d2_valid = true;
    assert_stack_clean("2-RDM read"); // the Expect sweeps must leave the LIFO stacks as they found them
}

// Build the full n_s x n_s spin-summed 1-RDM once per solve for the property read-out: the
// state-diagonal blocks and the transition (off-diagonal) blocks, in the delocalized basis. Each
// block is a direct 1-PDM Expect sweep over the extracted root(s) -- one-electron properties read
// one-electron RDMs, never a 2-RDM (that sweep is far heavier and the direct 1-PDM is more accurate
// at finite M). Diagonal i==j is a plain per-state 1-RDM (bra==ket); i<j is the bra!=ket transition
// 1-RDM, with (j,i) recovered as its transpose (exact against the symmetric property integrals).
static void ensure_dm_full(dmrgci_engine &e) {
    if (e.dmfull_valid)
        return;
    host_threads_guard htg;
    const int n = e.n_act;
    const size_t blk = (size_t)n * n;
    e.dmfull_cache.assign(blk * e.n_s * e.n_s, 0.0);

    std::vector<double> bt_scratch;
    if (e.localize_on) bt_scratch.resize(blk); // back-transform target (rotate1 forbids aliasing)
    std::vector<double> ro_scratch;            // un-permute target (Fiedler lattice order -> input)
    std::vector<int> iperm;                    // inverse of reorder_perm
    if (!e.reorder_perm.empty()) {
        ro_scratch.resize(blk);
        iperm.resize(n);
        for (int i = 0; i < n; i++) iperm[e.reorder_perm[i]] = i;
    }

    for (int i = 0; i < e.n_s; i++)
        for (int j = i; j < e.n_s; j++) {
            // Fresh extract of both roots so bra/ket share a canonical center (as in ensure_2rdm).
            const std::string itag = e.mps_info->tag + "-t" + std::to_string(i);
            std::shared_ptr<MultiMPS<SU2, double>> imps = e.mps->extract(i, itag);
            std::string jtag;
            std::shared_ptr<MultiMPS<SU2, double>> jmps = imps;
            if (j != i) {
                jtag = e.mps_info->tag + "-t" + std::to_string(j);
                jmps = e.mps->extract(j, jtag);
            }

            // NoTransposeRule: the transpose-symmetry simplification is invalid when bra != ket
            // (correct-and-uniform for the diagonal too).
            std::shared_ptr<MPO<SU2, double>> p1mpo =
                std::make_shared<PDM1MPOQC<SU2, double>>(e.hamil);
            p1mpo = std::make_shared<SimplifiedMPO<SU2, double>>(
                p1mpo,
                std::make_shared<NoTransposeRule<SU2, double>>(
                    std::make_shared<RuleQC<SU2, double>>()),
                true, true, OpNamesSet({OpNames::R, OpNames::RD}));
            auto p1me = std::make_shared<MovingEnvironment<SU2, double, double>>(p1mpo, imps, jmps,
                                                                                 "1PDM");
            p1me->init_environments(false);
            auto ex1 = std::make_shared<Expect<SU2, double, double>>(p1me, (ubond_t)e.cfg.m,
                                                                     (ubond_t)e.cfg.m);
            ex1->iprint = 0; // silence the per-site Expect sweep log
            ex1->solve(true, jmps->center == 0);
            GMatrix<double> d1 = ex1->get_1pdm_spatial(); // n x n row-major, on the block2 double stack

            double *bij = e.dmfull_cache.data() + (size_t)(i * e.n_s + j) * blk;
            std::copy(d1.data, d1.data + blk, bij);
            d1.deallocate(); // LIFO free of the stack matrix

            if (!e.reorder_perm.empty()) { // map out of block2's Fiedler lattice order
                for (int p = 0; p < n; p++)
                    for (int q = 0; q < n; q++)
                        ro_scratch[(size_t)p * n + q] = bij[(size_t)iperm[p] * n + iperm[q]];
                std::copy(ro_scratch.begin(), ro_scratch.end(), bij);
            }
            if (e.localize_on) { // rotate back to the delocalized basis
                rotate1(bij, e.U_loc.data(), n, bt_scratch.data(), /*forward=*/false);
                std::copy(bt_scratch.begin(), bt_scratch.end(), bij);
            }

            if (j == i) { // symmetrize the state-diagonal block (finite-M mildly breaks it)
                for (int p = 0; p < n; p++)
                    for (int q = p + 1; q < n; q++) {
                        double a = 0.5 * (bij[p * n + q] + bij[q * n + p]);
                        bij[p * n + q] = bij[q * n + p] = a;
                    }
            } else { // (j,i) block is the transpose of (i,j) for symmetric property integrals
                double *bji = e.dmfull_cache.data() + (size_t)(j * e.n_s + i) * blk;
                for (int p = 0; p < n; p++)
                    for (int q = 0; q < n; q++)
                        bji[(size_t)p * n + q] = bij[(size_t)q * n + p];
            }

            p1me->remove_partition_files();
            p1mpo->deallocate();
            remove_tag_files(itag); // the per-root extracts are transient
            if (j != i) remove_tag_files(jtag);
        }
    e.dmfull_valid = true;
    assert_stack_clean("full 1-RDM read");
}

// Reload the retained MultiMPS fresh from disk for a warm restart. After a solve, block2's LIFO
// memory stack is empty and the in-memory MultiMPS is a partially-deallocated shell whose StateInfos
// dangle into freed stack memory -- the on-disk copy is authoritative. load_data/load_mutable
// reallocate every StateInfo/SparseMatrixInfo on the heap (surviving any number of stack resets);
// the info is rebuilt against the CURRENT MPO basis. Mirrors MultiMPS::extract's reload, all roots.
static void reload_retained_mps(dmrgci_engine &e) {
    const std::string tag = e.mps_info->tag;
    auto info = std::make_shared<MultiMPSInfo<SU2>>(e.mpo->n_sites, e.hamil->vacuum,
                                                    std::vector<SU2>{e.target}, e.mpo->basis);
    info->tag = tag;
    info->load_mutable();
    auto mps = std::make_shared<MultiMPS<SU2, double>>(info);
    mps->load_data();
    mps->load_mutable();
    e.mps_info = info;
    e.mps = mps;
}

// Rotate the retained MultiMPS from the previous macro-iteration's active basis into the current one
// (localization-rotation warm start). R (e.R_active, orbital order [a*n+p]) is the overlap of the two
// localized active-MO sets; regularize it to the nearest proper rotation U~, take kappa = log(U~), and
// apply exp(kappa) to the MPS by RK4 time evolution -- block2's orbital-rotation recipe (two-site, so
// noise-free; the bond space rides along the SVD). NOPT already imports the new-basis integrals, so
// only the MPS is rotated (the integrals are not). e.mps is mutated in place on success. Returns false
// (=> caller cold-starts) if the rotation is untrustworthy: too much active-external leakage, an
// improper/complex generator, or a propagated norm that drifts too far from 1.
static bool rotate_retained_mps(dmrgci_engine &e) {
    const int n = e.n_act;
    const size_t nn = (size_t)n * n;

    // The regularization/log temporaries MUST live on the heap, not block2's LIFO stack: at the
    // rotation the stack is empty, so this math and the subsequent MPS-rotation init_environments
    // allocate from the same base region -- and block2 has a latent uninitialized read there, so
    // leftover eigs/log floating-point bytes on the stack cause a data-dependent segfault. (eigs and
    // logarithm already keep their own workspace on the heap; only these caller matrices need it.)
    auto heap = std::make_shared<VectorAllocator<double>>();

    // --- regularize R -> nearest orthogonal U~ = R (R^T R)^{-1/2} (Lowdin) ---
    std::vector<double> S(nn);
    cblas_dgemm(::CblasRowMajor, ::CblasTrans, ::CblasNoTrans, n, n, n, 1.0, e.R_active.data(), n,
                e.R_active.data(), n, 0.0, S.data(), n); // S = R^T R
    MatrixRef Sm(S.data(), n, n);
    DiagonalMatrix w(nullptr, n);
    w.allocate(heap);
    MatrixFunctions::eigs(Sm, w); // Sm rows <- eigenvectors, w <- eigenvalues; S now holds V^T
    double leak = 0.0, wmin = 1e30;
    for (int i = 0; i < n; i++) {
        leak = std::max(leak, std::fabs(w.data[i] - 1.0));
        wmin = std::min(wmin, w.data[i]);
    }
    if (wmin < 1e-6 || leak > 0.2) { // active space moved too much to carry the MPS
        w.deallocate(heap);
        fprintf(out_stream, "  warm-start: rotation leakage %.2e too large -> cold fallback\n", leak);
        return false;
    }
    std::vector<double> tmp(nn), Sinvh(nn), U(nn);
    for (int i = 0; i < n; i++) { // tmp = diag(w^{-1/2}) V^T
        double s = 1.0 / std::sqrt(w.data[i]);
        for (int j = 0; j < n; j++) tmp[(size_t)i * n + j] = s * S[(size_t)i * n + j];
    }
    w.deallocate(heap);
    cblas_dgemm(::CblasRowMajor, ::CblasTrans, ::CblasNoTrans, n, n, n, 1.0, S.data(), n, tmp.data(), n,
                0.0, Sinvh.data(), n); // S^{-1/2} = V diag(w^{-1/2}) V^T
    cblas_dgemm(::CblasRowMajor, ::CblasNoTrans, ::CblasNoTrans, n, n, n, 1.0, e.R_active.data(), n,
                Sinvh.data(), n, 0.0, U.data(), n); // U~ = R S^{-1/2}

    std::vector<double> Udet = U; // det() does an in-place LU, so work on a copy
    MatrixRef Udetm(Udet.data(), n, n);
    if (MatrixFunctions::det(Udetm) < 0.0) { // improper (reflection) -> log is complex
        fprintf(out_stream, "  warm-start: improper rotation (det<0) -> cold fallback\n");
        return false;
    }

    // --- kappa = log(U~) (real antisymmetric generator) ---
    ComplexMatrixRef ck(nullptr, n, n);
    ck.allocate(heap);
    ck.clear();
    ComplexMatrixFunctions::fill_complex(ck, MatrixRef(U.data(), n, n), MatrixRef(nullptr, n, n));
    ComplexMatrixFunctions::logarithm(ck);
    std::vector<double> kre(nn), kim(nn);
    ComplexMatrixFunctions::extract_complex(ck, MatrixRef(kre.data(), n, n),
                                            MatrixRef(kim.data(), n, n));
    ck.deallocate(heap);
    double imnorm = MatrixFunctions::norm(MatrixRef(kim.data(), n, n));
    if (imnorm > 1e-6) { // non-real generator -> not a proper rotation
        fprintf(out_stream, "  warm-start: complex generator (|Im k|=%.2e) -> cold fallback\n", imnorm);
        return false;
    }

    // NC MPO expects kappa^T; reindex into the frozen lattice order (fcidump->reorder(perm) put
    // orbital perm[i] on site i, so the generator on the lattice is kappa^T(perm[i], perm[j])).
    const bool have_perm = !e.reorder_perm.empty();
    std::vector<double> kap(nn);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            int oi = have_perm ? e.reorder_perm[i] : i;
            int oj = have_perm ? e.reorder_perm[j] : j;
            kap[(size_t)i * n + j] = kre[(size_t)oj * n + oi]; // transpose + permute
        }

    // Near convergence the localization-rotation vanishes (exp(-kappa) -> I), so the retained MPS
    // already sits in the current basis to machine precision -- rotating is pointless. It is also
    // unsafe: a numerically-zero one-body generator yields a degenerate rotation MPO whose
    // SimplifiedMPO-pruned operator schema block2's environment build reads past the end of
    // (initialize_tp -> a heap-dependent crash). Skip and reuse the reloaded MPS as-is; this is still
    // a warm restart (the short re-solve follows). The threshold sits far above the generator's noise
    // floor and far below any rotation an m-truncated DMRG could resolve.
    double kmax = 0.0;
    for (size_t t = 0; t < nn; t++) kmax = std::max(kmax, std::fabs(kap[t]));
    if (kmax < 1e-8) {
        return true; // negligible rotation (near convergence): reuse the reloaded MPS as-is
    }

    // --- one-body rotation MPO exp(-kappa t), built from the (lattice-order) generator ---
    auto fd_rot = std::make_shared<FCIDUMP<double>>();
    fd_rot->initialize_h1e(n, e.n_elec, e.twos, /*isym=*/0, 0.0, kap.data(), nn);
    auto hamil_rot = std::make_shared<HamiltonianQC<SU2, double>>(SU2(0), n, e.orbsym, fd_rot);
    std::shared_ptr<MPO<SU2, double>> mpo_rot =
        std::make_shared<MPOQC<SU2, double>>(hamil_rot, QCTypes::NC);
    mpo_rot->basis = hamil_rot->basis;
    mpo_rot = std::make_shared<SimplifiedMPO<SU2, double>>(
        mpo_rot,
        std::make_shared<AntiHermitianRuleQC<SU2, double>>(std::make_shared<RuleQC<SU2, double>>()),
        true);

    // Propagate exp(-kappa) over t in [0,1] as a multi-root TangentSpace TE sweep (block2's own
    // MultiMPS TE is complex-only, ket.size()==2; evolve_sa_multimps drives the per-root apply +
    // shared-basis truncation). -dt convention verified in the tests/block harness.
    const int rot_m = e.cfg.rot_m > 0 ? e.cfg.rot_m : e.cfg.m;
    // exp(-kappa) over t in [0,1] as n_steps TangentSpace TE sweeps (dt = 1/n_steps). One step
    // already reproduces the full rotation for small per-iter rotations (the sweep applies a
    // Krylov matrix-exponential per site, not an Euler step); rot_steps guards larger rotations.
    const int n_steps = e.cfg.rot_steps > 0 ? e.cfg.rot_steps : 10;
    const double dt = 1.0 / n_steps;
    const double norm2 = evolve_sa_multimps(e.mps, mpo_rot, (ubond_t)rot_m, dt, n_steps);
    mpo_rot->deallocate();

    if (std::fabs(norm2 - 1.0) > 1e-3) { // unitary propagation should preserve per-root norms
        fprintf(out_stream, "  warm-start: rotation norm^2=%.6f drifted -> cold fallback\n", norm2);
        return false;
    }
    return true; // rotated MPS in place; benign success is silent (keeps the CASSCF table clean)
}

// ---- ctor / dtor ---------------------------------------------------------------------
block2_casci_wrap::block2_casci_wrap(int n_act, int na, int nb, int mult, int n_s,
                                     int print_number, const dmrg_par &cfg)
    : impl_(std::make_unique<dmrgci_engine>(n_act, na + nb, na - nb, mult, n_s, print_number,
                                            cfg)) {
    host_threads_guard htg;
    int nthr = omp_get_max_threads();
    if (nthr < 1) nthr = 1;
    ensure_block2_runtime(cfg.save_dir, nthr);
}

block2_casci_wrap::~block2_casci_wrap() = default;

// ----------------------- casci_solver contract ----------------------------------------
void block2_casci_wrap::init_state_storage(int n_s, int) {
    impl_->E_states.assign(n_s, 0.0);
    impl_->storage_ready = true;
}
bool block2_casci_wrap::has_coef(int) const { return impl_->storage_ready; }
void block2_casci_wrap::set_act_rep_num(int *) {
    // Symmetry deferred: the DMRG block runs in C1, so orbsym stays all-zero and the
    // per-active-orbital irreps are ignored (mapped in once symmetry lands).
}
void block2_casci_wrap::set_localization(const double *U) {
    // Copy the localizing rotation (caller's buffer is transient); nullptr => no localization.
    dmrgci_engine &e = *impl_;
    if (U == nullptr) { e.localize_on = false; return; }
    e.U_loc.assign(U, U + (size_t)e.n_act * e.n_act);
    e.localize_on = true;
}
void block2_casci_wrap::set_active_rotation(const double *R) {
    // Retain the previous-basis -> current-basis active rotation for the next solve's warm restart.
    dmrgci_engine &e = *impl_;
    if (R == nullptr) { e.have_rotation = false; return; }
    e.R_active.assign(R, R + (size_t)e.n_act * e.n_act);
    e.have_rotation = true;
}
void block2_casci_wrap::set_report_rotation(const double *U) {
    // Retain the active-block canonicalization for the leading-configuration read-out (print_states).
    dmrgci_engine &e = *impl_;
    if (U == nullptr) { e.have_canon = false; return; }
    e.U_canon.assign(U, U + (size_t)e.n_act * e.n_act);
    e.have_canon = true;
}
void block2_casci_wrap::import_integrals(double *aaaa, double *f_act, double e_core) {
    dmrgci_engine &e = *impl_;
    const int n = e.n_act;
    // When localizing, rotate the integrals into the localized basis (E_core is invariant).
    double *h1 = f_act, *h2 = aaaa;
    if (e.localize_on) {
        e.F_loc.resize((size_t)n * n);
        e.g_loc.resize((size_t)n * n * n * n);
        rotate1(f_act, e.U_loc.data(), n, e.F_loc.data(), /*forward=*/true);
        rotate2(aaaa, e.U_loc.data(), n, e.g_loc.data(), /*forward=*/true);
        h1 = e.F_loc.data();
        h2 = e.g_loc.data();
    }
    // In-memory FCIDUMP (classic SU2 path). No rescale(): NOPT passes the embedded 1-e
    // Hamiltonian F_act (frozen core folded in) and chemist (tu|vw) directly
    e.fcidump = std::make_shared<FCIDUMP<double>>();
    e.fcidump->initialize_su2(n, e.n_elec, e.twos, /*isym=*/0, e_core,
                              h1, (size_t)n * n,
                              h2, (size_t)n * n * n * n);
    e.fcidump->set_orb_sym<int>(std::vector<int>(n, 0)); // C1

    // DMRG lattice order. block2's low-level path strings orbitals on the lattice in raw FCIDUMP
    // order; localized orbitals come out scrambled, so order them with block2's own Fiedler on the
    // exchange matrix K_ij = |(ij|ji)| and apply it via FCIDUMP::reorder. RDMs come back in this
    // order and are un-permuted in ensure_2rdm.
    // Warm-start freezes the order: a warm solve reuses the retained MPS's lattice order so the
    // rotated MPS stays consistent across the basis change (Fiedler order is a function of the
    // localized orbitals, so it must be pinned together with the frozen localization). A cold solve
    // recomputes it.
    if (e.have_rotation && !e.reorder_perm.empty()) {
        e.fcidump->reorder(e.reorder_perm); // frozen order (warm restart)
    } else {
        e.reorder_perm.clear();
        if (e.cfg.loc_order == DMRG_LOCORDER_FIEDLER) {
            std::vector<double> kmat((size_t)n * n, 0.0);
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    if (i != j)
                        kmat[(size_t)i * n + j] =
                            std::fabs(h2[(((size_t)i * n + j) * n + j) * n + i]);
            e.reorder_perm = OrbitalOrdering::fiedler((uint16_t)n, kmat);
            e.fcidump->reorder(e.reorder_perm);
        }
    }

    SU2 vacuum(0);
    e.hamil = std::make_shared<HamiltonianQC<SU2, double>>(vacuum, n, e.orbsym, e.fcidump);
    e.hamil->opf->seq->mode = SeqTypes::Tasked;
    e.mpo = build_qc_mpo(e.hamil);
}
int block2_casci_wrap::solve(int, int, bool) {
    dmrgci_engine &e = *impl_;
    host_threads_guard htg;
    assert_stack_clean("solve entry"); // block2 LIFO stacks must be empty between macro-iterations
    e.d2_valid = false; // new wavefunction -> any cached 2-RDM is stale
    e.dmfull_valid = false; // ... and the cached property 1-RDM

    // Warm-start: rotate the retained MPS into the current basis and re-solve from it. Requires a
    // usable retained MPS and a rotation from the host; the rotation itself may decline (return
    // false) if the basis change is too large, in which case we cold-start this iteration.
    bool warm = (e.cfg.warm_start == DMRG_WARM_ON && e.have_rotation && e.mps != nullptr &&
                 e.mps_info != nullptr);
    if (warm) {
        reload_retained_mps(e);        // fresh from disk (the in-memory shell is stale post-solve)
        if (e.cfg.warm_rotate == DMRG_WARM_ON)
            warm = rotate_retained_mps(e); // rotate into the current basis; false => cold fallback
        // else reuse-only: the reloaded MPS is the (unrotated) warm guess; the short re-solve corrects
        // the basis change. Proven crash-free and == cold; the safe fallback if rotation is declined.
    }
    e.have_rotation = false; // consumed; the host supplies a fresh R each warm iteration

    // --- sweep schedule: short warm re-solve vs full cold ramp ---
    dmrg_schedule sch;
    if (warm) {
        sch = build_warm_schedule(e.cfg.m, e.cfg.warm_sweeps, e.cfg.sweep_tol);
    } else if (e.cfg.schedule == DMRG_SCHED_DEFAULT) {
        sch = build_default_schedule(e.cfg.m, e.cfg.sweeps, e.cfg.sweep_tol);
    } else {
        fprintf(out_stream,
                "ERROR: DMRG schedule option not implemented yet (only 'default')\n");
        exit(0);
    }

    if (!warm) {
        // Cold start: drop the previous MPS and build a fresh random one. State-averaged MultiMPS
        // over the single target sector with nroots = n_s (n_s = 1 is just the single-state case).
        // state_specific stays off (block2 default), so the n_s roots share one renormalized basis
        // built from the equal-weight average density matrix -- exactly SA-CASSCF.
        if (e.mps_info != nullptr)
            remove_tag_files(e.mps_info->tag);
        Random::rand_seed(0);
        e.mps_info = std::make_shared<MultiMPSInfo<SU2>>(e.mpo->n_sites, e.hamil->vacuum,
                                                         std::vector<SU2>{e.target}, e.mpo->basis);
        e.mps_info->tag = "work_" + std::to_string(e.solve_count++); // unique per macro-iteration
        // --- initial MPS occupancy (only hf_occ=integral built; others provisioned) ---
        if (e.cfg.hf_occ == DMRG_HF_OCC_INTEGRAL) {
            e.mps_info->set_bond_dimension((ubond_t)e.cfg.m); // full FCI envelope
        } else {
            fprintf(out_stream,
                    "ERROR: DMRG hf_occ option not implemented yet (only 'integral')\n");
            exit(0);
        }
        e.mps = std::make_shared<MultiMPS<SU2, double>>(e.mpo->n_sites, 0, 2, e.n_s);
        e.mps->initialize(e.mps_info);
        e.mps->random_canonicalize();
        e.mps->save_mutable();
        e.mps->deallocate();
        e.mps_info->save_mutable();
        e.mps_info->deallocate_mutable();
    }
    // warm: e.mps is the rotated MultiMPS (in place) under its retained tag; e.mps_info is reused.

    auto me = std::make_shared<MovingEnvironment<SU2, double, double>>(e.mpo, e.mps, e.mps,
                                                                       "DMRG");
    me->delayed_contraction = OpNamesSet::normal_ops();
    me->cached_contraction = true;
    me->init_environments(false);

    auto dmrg = std::make_shared<DMRG<SU2, double, double>>(me, sch.bond_dims, sch.noises);
    dmrg->davidson_conv_thrds = sch.dav_thrds;
    dmrg->noise_type = NoiseTypes::ReducedPerturbativeCollected; // state-averaged perturbative noise
    dmrg->trunc_type = dmrg->trunc_type | TruncationTypes::RealDensityMatrix;
    dmrg->decomp_type = DecompositionTypes::DensityMatrix;
    dmrg->davidson_soft_max_iter = 200;
    dmrg->iprint = 0;
    dmrg->solve(sch.n_sweeps, e.mps->center == 0, e.cfg.sweep_tol);

    // per-root energies (ascending; root 0 = ground state). block2's energy precision is FPLS
    // (long double for FL=double), so bind via auto and narrow to NOPT's double.
    const auto &eng = dmrg->energies.back();
    e.E_states.assign(e.n_s, 0.0);
    for (int s = 0; s < e.n_s && s < (int)eng.size(); s++)
        e.E_states[s] = (double)eng[s];

    // sweeps actually run (block2 early-stops at sweep_tol) and the convergence achieved at the
    // last sweep, measured on the state-averaged energy (block2's own stopping criterion).
    e.last_n_sweeps = (int)dmrg->energies.size();
    e.last_sweep_dE = 0.0;
    if (dmrg->energies.size() >= 2) {
        const auto &e1 = dmrg->energies[dmrg->energies.size() - 1];
        const auto &e0 = dmrg->energies[dmrg->energies.size() - 2];
        double a1 = 0.0, a0 = 0.0;
        const int nr = (int)e1.size();
        for (int r = 0; r < nr; r++) { a1 += (double)e1[r]; a0 += (double)e0[r]; }
        if (nr > 0)
            e.last_sweep_dE = std::fabs(a1 / nr - a0 / nr);
    }
    // Exhausted the sweep budget without meeting block2's sweep_tol early-stop -> possibly
    // under-converged CI vector (flagged in the CAS-SCF table). A run that converges exactly on
    // its last scheduled sweep has last_sweep_dE <= sweep_tol and is not flagged.
    e.last_hit_max = (e.last_n_sweeps >= sch.n_sweeps) && (e.last_sweep_dE > e.cfg.sweep_tol);

    me->remove_partition_files(); // keep mps/mps_info alive for the RDM read-out
    return e.last_n_sweeps;
}
void block2_casci_wrap::calc_DM_diag(double *gamma, int /*i_set*/) {
    dmrgci_engine &e = *impl_;
    const int n = e.n_act;
    if (e.n_elec < 2) { // the 2-RDM trace cannot recover the 1-RDM for <2 electrons
        fprintf(out_stream,
                "ERROR: DMRG 1-RDM via 2-RDM trace needs N_elec>=2 (got %d)\n", e.n_elec);
        exit(0);
    }
    ensure_2rdm(e);
    // Per state: 1-RDM as a partial trace of that state's spatial 2-RDM:
    //   D1[p,s] = 1/(N-1) * sum_k D2[p,k,k,s]   (block2 D2 = <a+_p a+_q a_r a_s>)
    // gamma holds n_s consecutive n_act^2 blocks (CAS_engine averages them with the weights).
    const double inv = 1.0 / (e.n_elec - 1);
    const size_t blk = (size_t)n * n * n * n;
    for (int st = 0; st < e.n_s; st++) {
        const double *d2 = e.d2_cache.data() + (size_t)st * blk;
        double *g = gamma + (size_t)st * n * n;
        for (int p = 0; p < n; p++)
            for (int s = 0; s < n; s++) {
                double v = 0.0;
                for (int k = 0; k < n; k++)
                    v += d2[(((size_t)p * n + k) * n + k) * n + s];
                g[p * n + s] = inv * v;
            }
        // symmetrize (finite-M DMRG mildly breaks it)
        for (int p = 0; p < n; p++)
            for (int q = p + 1; q < n; q++) {
                double a = 0.5 * (g[p * n + q] + g[q * n + p]);
                g[p * n + q] = g[q * n + p] = a;
            }
    }
}
void block2_casci_wrap::G_calc(double *GAMMA) {
    dmrgci_engine &e = *impl_;
    const int n = e.n_act;
    ensure_2rdm(e);
    // Per state, NOPT layout GAMMA[p,q,r,s] = block2 D2[p,r,s,q]
    // GAMMA holds n_s consecutive n_act^4 blocks (CAS_engine averages them with the weights).
    const size_t blk = (size_t)n * n * n * n;
    for (int st = 0; st < e.n_s; st++) {
        const double *d2 = e.d2_cache.data() + (size_t)st * blk;
        double *G = GAMMA + (size_t)st * blk;
        for (int p = 0; p < n; p++)
            for (int q = 0; q < n; q++)
                for (int r = 0; r < n; r++)
                    for (int s = 0; s < n; s++)
                        G[(((size_t)p * n + q) * n + r) * n + s] =
                            d2[(((size_t)p * n + r) * n + s) * n + q];
    }
}
void block2_casci_wrap::calc_DMA(double *gamma, int a, int b) {
    // Properties (dipole/quadrupole/Mulliken/transition dipoles). block2 has a single SA state set,
    // so the host only ever asks for the full n_s x n_s block, gamma[(bra*n_s+ket)*n_act^2].
    if (a != 0 || b != 0) {
        fprintf(out_stream, "ERROR: DMRG backend expects a single state set (got a=%d b=%d)\n", a, b);
        exit(0);
    }
    dmrgci_engine &e = *impl_;
    ensure_dm_full(e);
    // SU2 get_1pdm_spatial is already spin-summed, so calc_DMA carries the full spin-summed 1-RDM
    // and calc_DMB is a no-op add (the host sums DMA+DMB). Add, matching aldet's accumulate.
    const size_t nel = e.dmfull_cache.size();
    for (size_t k = 0; k < nel; k++)
        gamma[k] += e.dmfull_cache[k];
}
void block2_casci_wrap::calc_DMB(double *, int, int) {
    // No-op: the spin-summed 1-RDM is entirely reported by calc_DMA (SU2 is spin-adapted).
}
int    block2_casci_wrap::n_states()      const { return impl_->n_s; }
int    block2_casci_wrap::mult()          const { return impl_->mult; }
double block2_casci_wrap::E_state(int i)  const { return impl_->E_states[i]; }
double block2_casci_wrap::S2_state(int)   const { double S = impl_->twos / 2.0; return S * (S + 1.0); }
double block2_casci_wrap::L2_state(int)   const { return 0.0; } // linear-molecule Lambda: deferred
double block2_casci_wrap::P_state(int)    const { return 0.0; } // parity: deferred
double *block2_casci_wrap::E_states_ptr() const { return impl_->E_states.data(); }
double block2_casci_wrap::last_solve_resid() const { return impl_->last_sweep_dE; }
bool block2_casci_wrap::last_solve_hit_max() const { return impl_->last_hit_max; }
void block2_casci_wrap::gen_ext_ind() { /* aldet determinant index tables; n/a for an MPS backend */ }
// Compress a single MPS to bond dimension target_m by fitting a small-m MPS to it through the identity
// MPO (block2 Linear "Normal" equation). State-preserving projection used only for the read-out: it
// makes the DeterminantTRIE search cheaper (cost ~ m^2) on the (higher-m) rotated canonical MPS, while
// preserving the leading determinants. ctag names the fitted MPS's scratch; the caller removes it. Env
// scratch is cleaned here.
static std::shared_ptr<MPS<SU2, double>>
compress_single_mps(dmrgci_engine &e, const std::shared_ptr<MPS<SU2, double>> &ket,
                    int target_m, const std::string &ctag) {
    std::shared_ptr<MPO<SU2, double>> impo = std::make_shared<IdentityMPO<SU2, double>>(e.hamil);
    impo = std::make_shared<SimplifiedMPO<SU2, double>>(impo, std::make_shared<Rule<SU2, double>>());

    auto binfo = std::make_shared<MPSInfo<SU2>>(e.mpo->n_sites, e.hamil->vacuum, e.target,
                                                e.mpo->basis);
    binfo->set_bond_dimension((ubond_t)target_m);
    binfo->tag = ctag;
    auto bra = std::make_shared<MPS<SU2, double>>(e.mpo->n_sites, ket->center, ket->dot);
    bra->initialize(binfo);
    bra->random_canonicalize();
    bra->save_mutable();
    bra->deallocate();
    binfo->save_mutable();
    binfo->deallocate_mutable();

    auto cme = std::make_shared<MovingEnvironment<SU2, double, double>>(impo, bra, ket, "CPS");
    cme->init_environments(false);
    // ket bond dim from the ket itself (the rotated MPS can exceed the solve m), so the fit does not
    // truncate the ket -- only the bra is compressed to target_m.
    std::vector<ubond_t> bdim{(ubond_t)target_m}, kdim{ket->info->get_max_bond_dimension()};
    std::vector<double> noises{1e-9, 0.0}; // one noisy sweep to seed the fit, then clean
    auto cps = std::make_shared<Linear<SU2, double, double>>(cme, bdim, kdim, noises);
    cps->iprint = 0;
    cps->solve(8, ket->center == 0);

    cme->remove_partition_files();
    impo->deallocate();
    return bra;
}

// The discrete (signed-permutation) part of the solve->canonical basis change, peeled off so only a
// proper rotation reaches the log. g maps canonical orbital -> residual; sigma is the per-canonical-orbital
// sign the rotation cannot carry (see peel_report_permutation).
struct report_peel {
    std::vector<int> g;         // canonical p -> residual index (g = id recovers the full-Rp read-out)
    std::vector<double> sigma;  // canonical p -> +/-1: canon_p = sigma_p * residual_{g[p]}
};

// Peel the discrete part off the solve->canonical basis change Rp (= R^T; column p is canonical orbital p in
// solve coordinates). rotate_multimps_to_canonical evolves under exp(-log W), so a near-permutation Rp
// makes log(Rp) hit the matrix-log branch cut (an eigenvalue near -1 -> complex generator) -- the
// canon-only (non-localized) case, where the canonical orbitals are just the ascending-Fock-eigenvalue
// reordering of the solve orbitals. Factor Rp = Rt * Q with Q a signed permutation and Rt ~ I proper
// (right factor: Q mixes the canonical/column index). Rotating by the residual Rt reads residual orbitals
// (Fact F: site reads a column of the fed matrix). Overwrites Rp in place with Rt and returns Q as
// (g, sigma): canonical orbital p is realized as sigma_p * (residual g[p]), so canon_p = sigma_p *
// residual_{g[p]}. The reflection sigma is a per-orbital orbital-sign gauge -- NOT a free constant: it
// is theory-fixed by Rp and must be reinstated on the read-out (a determinant with p singly occupied
// carries sigma_p; doubly/empty carry sigma_p^2 = +1 / nothing). The det-guard flip (an improper Rp
// realized as a residual column negation) is folded into sigma so nothing is dropped. Greedy
// max-magnitude matching; exact when Rp is a signed permutation.
static report_peel peel_report_permutation(double *Rp, int n) {
    struct Cell { double mag; int i, j; };
    std::vector<Cell> cells;
    cells.reserve((size_t)n * n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            cells.push_back({std::fabs(Rp[(size_t)i * n + j]), i, j});
    std::sort(cells.begin(), cells.end(),
              [](const Cell &a, const Cell &b) { return a.mag > b.mag; });
    std::vector<int> rho(n, -1);        // rho[a] = canonical column matched to solve row a
    std::vector<double> sgn(n, 1.0);    // sgn[b] = aligned-entry sign of residual column b
    std::vector<char> row_used(n, 0), col_used(n, 0);
    for (const Cell &c : cells) {
        if (row_used[c.i] || col_used[c.j]) continue;
        row_used[c.i] = col_used[c.j] = 1;
        rho[c.i] = c.j;
        sgn[c.i] = Rp[(size_t)c.i * n + c.j] < 0.0 ? -1.0 : 1.0;
    }
    // Right factor Rt = Rp * Q: column b of Rt = sgn[b] * (column rho[b] of Rp); diagonal Rt[b][b] > 0.
    std::vector<double> Rt((size_t)n * n);
    for (int a = 0; a < n; a++)
        for (int b = 0; b < n; b++)
            Rt[(size_t)a * n + b] = sgn[b] * Rp[(size_t)a * n + rho[b]];
    // Keep Rt proper (Rp improper is absorbed as an odd Q sign); a negative det would give a -1
    // eigenvalue -> complex log. Flip the least-aligned residual column (smallest diagonal) and fold the
    // flip into sgn so it is reinstated with the rest of Q. exp(-log Rt) still guards a residual that
    // stays complex.
    if (MatrixFunctions::det(MatrixRef(Rt.data(), n, n)) < 0.0) {
        int bmin = 0;
        for (int b = 1; b < n; b++)
            if (std::fabs(Rt[(size_t)b * n + b]) < std::fabs(Rt[(size_t)bmin * n + bmin])) bmin = b;
        for (int a = 0; a < n; a++) Rt[(size_t)a * n + bmin] = -Rt[(size_t)a * n + bmin];
        sgn[bmin] = -sgn[bmin];
    }
    std::copy(Rt.begin(), Rt.end(), Rp);
    report_peel rp;
    rp.g.assign(n, -1);          // report -> residual: g = rho^{-1}
    rp.sigma.assign(n, 1.0);
    for (int a = 0; a < n; a++) rp.g[rho[a]] = a;
    for (int p = 0; p < n; p++) rp.sigma[p] = sgn[rp.g[p]]; // sigma_p = aligned sign of residual g[p]
    return rp;
}

// Rotate a MultiMPS by the residual continuous rotation between the solve and canonical bases. The full
// solve->canonical rotation (localized or canonical; see print_states) is a signed permutation times a
// residual; the caller peels the permutation off (peel_report_permutation) and passes only the residual
// U (n_act x n_act, [a*n+p], ~ I, proper). kappa = log(U) is reindexed into the frozen lattice (Fiedler)
// order and the one-body MPO exp(-kappa) is applied by rot_steps TangentSpace TE sweeps at bond dim
// rot_m -- block2's orbital-rotation recipe, shared with the warm-start (evolve_sa_multimps). Mutates
// mps in place; a no-op when the residual is negligible (an already-canonical or pure-permutation
// solve). Unlike the warm-start there is no cold fallback: the read-out accepts the truncation,
// reported downstream as the captured weight. Returns false (leaving mps untouched) only if the
// generator still comes back complex despite the peel -- the caller then reports in the solve basis
// without the permutation relabel.
static bool rotate_multimps_to_canonical(dmrgci_engine &e,
                                         const std::shared_ptr<MultiMPS<SU2, double>> &mps,
                                         const double *U, int rot_m, int rot_steps) {
    const int n = e.n_act;
    const size_t nn = (size_t)n * n;
    auto heap = std::make_shared<VectorAllocator<double>>(); // log temporaries off block2's LIFO stack

    std::vector<double> Um(U, U + nn);
    ComplexMatrixRef ck(nullptr, n, n);
    ck.allocate(heap);
    ck.clear();
    ComplexMatrixFunctions::fill_complex(ck, MatrixRef(Um.data(), n, n), MatrixRef(nullptr, n, n));
    ComplexMatrixFunctions::logarithm(ck); // kappa = log(U); real antisymmetric for a proper rotation
    std::vector<double> kre(nn), kim(nn);
    ComplexMatrixFunctions::extract_complex(ck, MatrixRef(kre.data(), n, n), MatrixRef(kim.data(), n, n));
    ck.deallocate(heap);
    if (MatrixFunctions::norm(MatrixRef(kim.data(), n, n)) > 1e-6) {
        fprintf(out_stream, "  warning: read-out rotation generator not real -> reporting in solve basis\n");
        return false;
    }

    // NC MPO expects kappa^T; reindex into the frozen lattice order (as the warm-start rotation does).
    const bool have_perm = !e.reorder_perm.empty();
    std::vector<double> kap(nn);
    double kmax = 0.0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            int oi = have_perm ? e.reorder_perm[i] : i;
            int oj = have_perm ? e.reorder_perm[j] : j;
            kap[(size_t)i * n + j] = kre[(size_t)oj * n + oi];
            kmax = std::max(kmax, std::fabs(kap[(size_t)i * n + j]));
        }
    if (kmax < 1e-8) return true; // exp(-kappa) ~ I: the MPS already sits in the reporting basis

    auto fd_rot = std::make_shared<FCIDUMP<double>>();
    fd_rot->initialize_h1e(n, e.n_elec, e.twos, /*isym=*/0, 0.0, kap.data(), nn);
    auto hamil_rot = std::make_shared<HamiltonianQC<SU2, double>>(SU2(0), n, e.orbsym, fd_rot);
    std::shared_ptr<MPO<SU2, double>> mpo_rot =
        std::make_shared<MPOQC<SU2, double>>(hamil_rot, QCTypes::NC);
    mpo_rot->basis = hamil_rot->basis;
    mpo_rot = std::make_shared<SimplifiedMPO<SU2, double>>(
        mpo_rot,
        std::make_shared<AntiHermitianRuleQC<SU2, double>>(std::make_shared<RuleQC<SU2, double>>()),
        true);
    evolve_sa_multimps(mps, mpo_rot, (ubond_t)rot_m, 1.0 / rot_steps, rot_steps); // exp(-kappa t), t in [0,1]
    mpo_rot->deallocate();
    return true;
}

using det_occ = std::vector<uint8_t>; // per spatial orbital: 0 empty, 1 alpha, 2 beta, 3 doubly
using det_expansion = std::map<det_occ, double>;

static void add_det_coeff(det_expansion &dets, const det_occ &occ, double coef) {
    if (std::fabs(coef) < 1e-14) return;
    double &slot = dets[occ];
    slot += coef;
    if (std::fabs(slot) < 1e-14) dets.erase(occ);
}

// Expand a block2 SU2 step-vector CSF into highest-weight Slater determinants in the same lattice
// order. The coefficients are Clebsch-Gordan products for coupling each singly occupied orbital to the
// running spin path (u: S -> S+1/2, d: S -> S-1/2). The determinant phase is still the local
// site-ordered block2 phase; canonicalize_site_det converts it to alpha-string/beta-string order.
static void expand_csf_step_rec(const std::vector<uint8_t> &step, int site, int twos, int twom,
                                det_occ &det, double coef, int target_twos,
                                det_expansion &out) {
    const int n = (int)step.size();
    if (site == n) {
        if (twos == target_twos && twom == target_twos) add_det_coeff(out, det, coef);
        return;
    }

    const uint8_t s = step[site] & 3;
    if (s == 0 || s == 3) {
        det[site] = (s == 3) ? 3 : 0;
        expand_csf_step_rec(step, site + 1, twos, twom, det, coef, target_twos, out);
        det[site] = 0;
        return;
    }

    const double j = 0.5 * twos;
    const double m = 0.5 * twom;
    const double den = 2.0 * j + 1.0;
    if (s == 1) { // spin-path increase
        det[site] = 1; // alpha
        expand_csf_step_rec(step, site + 1, twos + 1, twom + 1, det,
                            coef * std::sqrt(std::max(0.0, (j + m + 1.0) / den)),
                            target_twos, out);
        det[site] = 2; // beta
        expand_csf_step_rec(step, site + 1, twos + 1, twom - 1, det,
                            coef * std::sqrt(std::max(0.0, (j - m + 1.0) / den)),
                            target_twos, out);
    } else if (twos > 0) { // spin-path decrease: block2's left-coupled Condon-Shortley convention
        det[site] = 1; // alpha
        expand_csf_step_rec(step, site + 1, twos - 1, twom + 1, det,
                            -coef * std::sqrt(std::max(0.0, (j - m) / den)),
                            target_twos, out);
        det[site] = 2; // beta
        expand_csf_step_rec(step, site + 1, twos - 1, twom - 1, det,
                            coef * std::sqrt(std::max(0.0, (j + m) / den)),
                            target_twos, out);
    }
    det[site] = 0;
}

static det_expansion expand_csf_step(const std::vector<uint8_t> &step, int target_twos) {
    det_expansion out;
    det_occ det(step.size(), 0);
    expand_csf_step_rec(step, 0, 0, 0, det, 1.0, target_twos, out);
    return out;
}

struct det_image {
    det_occ occ;
    double phase = 1.0;
};

// Convert a determinant from the current MPS lattice/site convention into canonical alpha-string /
// beta-string convention. block2's non-spin-adapted determinant/SCI Fock convention uses interleaved
// spin orbitals (alpha0,beta0,alpha1,beta1,...) for fermion phases, while aldet prints separate
// alpha/beta strings. The inversion count below applies that spin-order phase together with the
// Fiedler/site orbital permutation and the chosen canonical-orbital sign gauge.
static det_image canonicalize_site_det(const det_occ &site_det,
                                       const std::vector<int> &canon_of_lattice,
                                       const std::vector<double> &lattice_sign) {
    const int n = (int)site_det.size();
    det_image img;
    img.occ.assign(n, 0);
    std::vector<int> spin_orbs;
    spin_orbs.reserve((size_t)2 * n);

    for (int k = 0; k < n; k++) {
        const uint8_t occ = site_det[k] & 3;
        const int p = canon_of_lattice[k];
        const double sg = lattice_sign[k] < 0.0 ? -1.0 : 1.0;
        if (occ & 1) {
            img.occ[p] |= 1;
            img.phase *= sg;
            spin_orbs.push_back(p); // alpha block
        }
        if (occ & 2) {
            img.occ[p] |= 2;
            img.phase *= sg;
            spin_orbs.push_back(n + p); // beta block
        }
    }

    int odd = 0;
    for (int i = 0; i < (int)spin_orbs.size(); i++)
        for (int j = i + 1; j < (int)spin_orbs.size(); j++)
            if (spin_orbs[i] > spin_orbs[j]) odd ^= 1;
    if (odd) img.phase = -img.phase;
    return img;
}

static void normalize_report_sign_gauge(const std::vector<int> &canon_of_lattice,
                                        std::vector<double> &lattice_sign) {
    double prod = 1.0;
    for (double s : lattice_sign)
        prod *= (s < 0.0 ? -1.0 : 1.0);
    if (prod > 0.0 || canon_of_lattice.empty())
        return;

    // LAPACK eigenvector signs are arbitrary. After the peel, an odd number of signs in the discrete
    // map is only a printed canonical-MO phase choice; use the equivalent gauge with det(Q)=+1.
    for (int k = 0; k < (int)canon_of_lattice.size(); k++)
        if (canon_of_lattice[k] == 0) {
            lattice_sign[k] = -lattice_sign[k];
            return;
        }
}

static double expansion_weight(const det_expansion &dets) {
    double w = 0.0;
    for (const auto &kv : dets) w += kv.second * kv.second;
    return w;
}

struct su2_readout_expansion {
    det_expansion determinants;
    double weight = 0.0;
};

// Extract the short SU2 step-vector expansion directly from the MPS, then expand the retained block2
// CSFs locally into determinant coefficients in the requested canonical orbital order/sign gauge. This
// avoids the expensive SU2->SZ MPS transform; only the retained read-out subspace is expanded.
static su2_readout_expansion canonical_readout_from_su2_mps(
    const std::shared_ptr<UnfusedMPS<SU2, double>> &su2_umps,
    int n_sites, int target_twos, double cutoff,
    const std::vector<int> &canon_of_lattice,
    const std::vector<double> &lattice_sign) {
    su2_readout_expansion out;

    auto ctrie = std::make_shared<DeterminantTRIE<SU2, double>>(n_sites, true);
    ctrie->evaluate(su2_umps, cutoff);
    for (int t = 0; t < (int)ctrie->size(); t++) {
        const double csf_coef = (double)ctrie->vals[t];
        if (std::fabs(csf_coef) < 1e-14) continue;
        det_expansion site_dets = expand_csf_step((*ctrie)[t], target_twos);
        for (const auto &kv : site_dets) {
            det_image img = canonicalize_site_det(kv.first, canon_of_lattice, lattice_sign);
            add_det_coeff(out.determinants, img.occ, csf_coef * kv.second * img.phase);
        }
    }

    out.weight = expansion_weight(out.determinants);
    return out;
}

static void det_strings(const det_occ &det, std::string &alpha, std::string &beta) {
    alpha.assign(det.size(), '0');
    beta.assign(det.size(), '0');
    for (int i = 0; i < (int)det.size(); i++) {
        if (det[i] & 1) alpha[i] = '1';
        if (det[i] & 2) beta[i] = '1';
    }
}

static void report_leading_determinants(std::FILE *out, int state_idx, double E, double S2,
                                        const det_expansion &dets, int print_number) {
    std::fprintf(out,
                 "State %d  E  = % 18.10f S^2 = %.2f:\n",
                 state_idx, E, S2);

    std::vector<det_expansion::const_iterator> ord;
    ord.reserve(dets.size());
    for (auto it = dets.begin(); it != dets.end(); ++it) ord.push_back(it);
    const int top = (int)std::min<size_t>(print_number, ord.size());
    std::partial_sort(ord.begin(), ord.begin() + top, ord.end(),
                      [](det_expansion::const_iterator a, det_expansion::const_iterator b) {
                          return std::fabs(a->second) > std::fabs(b->second);
                      });
    std::string alpha, beta;
    for (int i = 0; i < top; i++) {
        det_strings(ord[i]->first, alpha, beta);
        std::fprintf(out, "% .10e  | %s | %s\n", ord[i]->second, alpha.c_str(), beta.c_str());
    }
}

static void report_determinant_weights(std::FILE *out, const std::vector<double> &weights) {
    std::fprintf(out, "\nDMRG determinant extraction weights:\n");
    std::fprintf(out, "____________________________\n");
    std::fprintf(out, " State| weight       |\n");
    std::fprintf(out, "______|______________|\n");
    for (int st = 0; st < (int)weights.size(); st++)
        std::fprintf(out, "%5d | %.6f     |\n", st, weights[st]);
    std::fprintf(out, "______|______________|\n");
}

// Report each state's leading determinant expansion in the canonical active-orbital basis. The DMRG
// solve runs in the localized + Fiedler-ordered basis; each state's MPS is rotated by the proper
// (continuous) part of the solve->canonical change. The retained block2 SU2 vector is expanded locally
// into determinants, and the remaining discrete map (lattice orbital -> canonical orbital plus
// orbital order/sign) is applied to those determinants.
void block2_casci_wrap::print_states(int, int, int print) {
    if (!print) return;
    dmrgci_engine &e = *impl_;
    if (e.mps == nullptr) return;
    if (e.cfg.print_dets != DMRG_WARM_ON) { // read-out opted out ($DMRG print_dets = off)
        fprintf(out_stream, "  (leading-configuration read-out disabled: $DMRG print_dets = off)\n");
        return;
    }
    host_threads_guard htg;

    const int n = e.n_act;
    const double S = e.twos / 2.0, S2 = S * (S + 1.0);
    const double cutoff = e.cfg.extract_cutoff;
    const int print_number = e.print_number; // leading determinants per state (p_n, $act_space)
    const int rot_m = e.cfg.det_rot_m, rot_steps = e.cfg.det_rot_steps, cm = e.cfg.extract_m;

    // Fiedler reorder (empty => input order): the DMRG lattice site k holds solve orbital reord[k].
    std::vector<int> reord(e.reorder_perm.begin(), e.reorder_perm.end());
    std::vector<int> iperm;
    if (!reord.empty()) {
        iperm.resize(n);
        for (int i = 0; i < n; i++) iperm[reord[i]] = i;
    }

    // Solve -> canonical basis change: R = <canon|solve> is U_canon*U_loc (localizing), U_canon
    // (canonicalizing a delocalized solve), or U_loc. rotate_multimps_to_canonical consumes Rp = R^T
    // (column p is canonical orbital p in solve coordinates). peel_report_permutation splits Rp = Rt * Q:
    // the proper residual Rt (real, short log) drives the rotation, and the signed permutation Q =
    // (g, sigma) is the discrete map from the extracted (lattice-order) CSFs to the canonical orbitals.
    std::vector<double> R_total, Rp;
    report_peel rp;
    const double *R = nullptr;
    if (e.have_canon && e.localize_on) {
        R_total.assign((size_t)n * n, 0.0);
        cblas_dgemm(::CblasRowMajor, ::CblasNoTrans, ::CblasNoTrans, n, n, n, 1.0, e.U_canon.data(),
                    n, e.U_loc.data(), n, 0.0, R_total.data(), n);
        R = R_total.data();
    } else if (e.have_canon) {
        R = e.U_canon.data();
    } else if (e.localize_on) {
        R = e.U_loc.data();
    }
    if (R != nullptr) {
        Rp.assign((size_t)n * n, 0.0);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) Rp[(size_t)i * n + j] = R[(size_t)j * n + i]; // Rp = R^T
        rp = peel_report_permutation(Rp.data(), n); // Rp <- residual Rt; rp = Q = (g, sigma)
    }

    // Discrete map Q -> per-lattice-site (canonical orbital, sign): canonical orbital p is realized at the
    // lattice site of residual g[p] (Fact F) with sign sigma_p. Without a rotation (R == nullptr) the
    // read-out sits in the solve basis: lattice site k holds solve orbital reord[k], signs +1.
    std::vector<int> canon_of_lattice(n);
    std::vector<double> lattice_sign(n, 1.0);
    if (R != nullptr) {
        for (int p = 0; p < n; p++) {
            const int resid = rp.g[p];
            const int k = reord.empty() ? resid : iperm[resid]; // lattice site of residual resid
            canon_of_lattice[k] = p;
            lattice_sign[k] = rp.sigma[p];
        }
        normalize_report_sign_gauge(canon_of_lattice, lattice_sign);
    } else {
        for (int k = 0; k < n; k++) canon_of_lattice[k] = reord.empty() ? k : reord[k];
    }
    const char *basis_label = e.have_canon ? "canonical orbital" : "active orbital";
    report_dmrg_orbital_map(out_stream, n, canon_of_lattice, lattice_sign, basis_label);

    std::vector<int> solve_of_lattice(n);
    std::vector<double> plus_sign(n, 1.0);
    for (int k = 0; k < n; k++) solve_of_lattice[k] = reord.empty() ? k : reord[k];
    std::vector<double> det_weights((size_t)e.n_s, 0.0);

    for (int st = 0; st < e.n_s; st++) {
        // Work on a copy of this root; e.mps stays in the localized basis for the property read-out.
        const std::string xtag = e.mps_info->tag + "-cf" + std::to_string(st);
        std::shared_ptr<MultiMPS<SU2, double>> imps = e.mps->extract(st, xtag);
        // Rotate by the proper residual. A complex residual (guarded inside) leaves the MPS in the solve
        // basis, where the reported orbital map no longer holds -- flag it.
        bool map_applies = true;
        if (!Rp.empty() && !rotate_multimps_to_canonical(e, imps, Rp.data(), rot_m, rot_steps)) {
            fprintf(out_stream, "  warning: state %d rotation fell back to solve basis; determinant orbital "
                                "map does not apply\n", st);
            map_applies = false;
        }

        const std::string stag = e.mps_info->tag + "-cs" + std::to_string(st);
        std::shared_ptr<MPS<SU2, double>> smps = imps->make_single(stag);

        // Compress the canonical-basis MPS for a cheaper TRIE search (opt-out with extract_m=0); the leading
        // configurations are low-rank-robust, so this preserves them (loss shows up as captured weight).
        std::string ctag;
        std::shared_ptr<MPS<SU2, double>> rmps = smps;
        if (cm > 0) {
            ctag = e.mps_info->tag + "-cc" + std::to_string(st);
            rmps = compress_single_mps(e, smps, cm, ctag);
        }

        std::shared_ptr<UnfusedMPS<SU2, double>> umps =
            std::make_shared<UnfusedMPS<SU2, double>>(rmps);

        const std::vector<int> &det_orb_map = map_applies ? canon_of_lattice : solve_of_lattice;
        const std::vector<double> &det_sign = map_applies ? lattice_sign : plus_sign;
        su2_readout_expansion canonical =
            canonical_readout_from_su2_mps(umps, n, e.twos, cutoff, det_orb_map, det_sign);

        det_weights[st] = canonical.weight;
        report_leading_determinants(out_stream, st, e.E_states[st], S2, canonical.determinants,
                                    print_number);

        remove_tag_files(xtag);
        remove_tag_files(stag);
        if (!ctag.empty()) remove_tag_files(ctag);
    }
    report_determinant_weights(out_stream, det_weights);
    assert_stack_clean("leading-config read");
}
void block2_casci_wrap::write_civec(int, char *) { /* MPS write-out deferred */ }
