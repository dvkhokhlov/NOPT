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
#include <string>
#include <system_error>
#include <unistd.h>     // getpid
#include <omp.h>

#include "common_vars.h"      // out_stream
#include "localized_dmrg.h"   // rotate1/rotate2 (active-space basis transforms)

using namespace block2;

// ---- process-global block2 runtime ---------------------------------------------------
namespace {

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
// always delimited, in every per-tag filename (mps.hpp get_filename):
//   <prefix>.MPS.<tag>.<i>        MPS site tensors
//   <prefix>.MPS.INFO.<tag>.*     MPSInfo bond dims
//   <tag>-mps_info.bin            serialized MPSInfo
// The trailing '.' on the tensor/info prefixes keeps tag "work_1" from matching "work_10".
void remove_tag_files(const std::string &tag) {
    auto fr = frame_<double>();
    if (fr == nullptr)
        return;
    std::error_code ec;
    const std::string t_mps  = fr->prefix + ".MPS." + tag + ".";
    const std::string t_info = fr->prefix + ".MPS.INFO." + tag + ".";
    const std::string t_bin  = tag + "-mps_info.bin";
    std::vector<std::filesystem::path> victims;
    for (auto &de : std::filesystem::directory_iterator(fr->mps_dir, ec)) {
        const std::string name = de.path().filename().string();
        if (name.rfind(t_mps, 0) == 0 || name.rfind(t_info, 0) == 0 || name == t_bin)
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

} // namespace

// -------------------------- engine: all block2 state ----------------------------------
struct dmrgci_engine {
    dmrg_par cfg;
    int n_act, n_elec, twos, mult, n_s;
    SU2 target;                       // active-space symmetry sector (C1: pg = 0)
    std::vector<uint8_t> orbsym;      // per-orbital irrep; C1 -> all 0 (set in set_act_rep_num)
    std::vector<double> E_states;     // per-state energies (returned by E_states_ptr)
    bool storage_ready = false;       // has init_state_storage run?

    // Localization. Empty/off => solve in the given basis.
    std::vector<double> U_loc;        // n_act x n_act, [a*n_act+p]; valid when localize_on
    bool localize_on = false;
    std::vector<double> F_loc, g_loc; // integrals rotated into the localized basis (when localize_on)

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
    std::vector<uint16_t> reorder_perm;         // DMRG lattice order (Fiedler); empty => input order

    dmrgci_engine(int n_act_, int n_elec_, int twos_, int mult_, int n_s_, const dmrg_par &c)
        : cfg(c), n_act(n_act_), n_elec(n_elec_), twos(twos_), mult(mult_), n_s(n_s_),
          target(n_elec_, twos_, 0), orbsym(n_act_, 0), E_states(n_s_, 0.0) {}
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

// ---- ctor / dtor ---------------------------------------------------------------------
block2_casci_wrap::block2_casci_wrap(int n_act, int na, int nb, int mult, int n_s,
                                     const dmrg_par &cfg)
    : impl_(std::make_unique<dmrgci_engine>(n_act, na + nb, na - nb, mult, n_s, cfg)) {
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

    // remove previous step MPS
    if (e.mps_info != nullptr)
        remove_tag_files(e.mps_info->tag);

    // --- sweep schedule (only "default" built; others provisioned) ---
    dmrg_schedule sch;
    if (e.cfg.schedule == DMRG_SCHED_DEFAULT) {
        sch = build_default_schedule(e.cfg.m, e.cfg.sweeps, e.cfg.sweep_tol);
    } else {
        fprintf(out_stream,
                "ERROR: DMRG schedule option not implemented yet (only 'default')\n");
        exit(0);
    }

    // State-averaged MultiMPS over the single target sector with nroots = n_s (n_s = 1 is just the
    // single-state case). state_specific stays off (block2 default), so the n_s roots share one
    // renormalized basis built from the equal-weight average density matrix -- exactly SA-CASSCF.
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
void block2_casci_wrap::gen_ext_ind() { /* aldet determinant index tables; n/a for an MPS backend */ }
void block2_casci_wrap::print_states(int, int, int) { /* MPS CSF/det printout deferred (DeterminantTRIE) */ }
void block2_casci_wrap::write_civec(int, char *) { /* MPS write-out deferred */ }
