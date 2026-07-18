// block2_casci_wrap — casci_solver backend over the external block2 DMRG library.

#include "block2_dmrg_engine.h"   // block2 headers + dmrgci_engine + shared helpers

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cerrno>
#include <csignal>
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
// Shared block2-backend helpers. Named, not anonymous: the split TUs link against the cross-TU ones.
namespace nopt_block2 {


// Remove scratch dirs of runs that died before their destructor ran. Dirs are named <pid>_<unique>;
// one whose pid is alive or unsignalable is left alone. The unique suffix keeps a recycled pid from
// matching a live run's dir. Not covered: a scratch root shared across pid namespaces.
static void reap_orphan_scratch(const std::string &root) {
    std::error_code ec;
    const std::string tag = "nopt_dmrg_";
    for (const auto &de : std::filesystem::directory_iterator(root, ec)) {
        if (!de.is_directory(ec))
            continue;
        const std::string name = de.path().filename().string();
        if (name.rfind(tag, 0) != 0)
            continue;
        const std::string rest = name.substr(tag.size());
        const size_t us = rest.find('_');
        if (us == 0 || us == std::string::npos)
            continue; // no unique suffix: not one of ours to judge
        const std::string digits = rest.substr(0, us);
        if (digits.find_first_not_of("0123456789") != std::string::npos)
            continue;
        const long pid = strtol(digits.c_str(), nullptr, 10);
        if (pid <= 0 || pid == (long)getpid())
            continue;
        if (kill((pid_t)pid, 0) == 0 || errno == EPERM)
            continue; // still running (EPERM: alive, just not ours)
        std::filesystem::remove_all(de.path(), ec);
    }
}

// Process-global block2 runtime: owns the frame_/threading_ globals and the per-process scratch dir,
// and tears them down at program exit.
struct Block2Runtime {
    std::string scratch;

    Block2Runtime(const std::string &save_dir_root, double memory_gb, int n_threads) {
        // Per-run scratch under the configured root (default /dev/shm), so concurrent runs never
        // share block2's renormalized-operator files. Identity is <pid>_<unique>: the pid for
        // reap_orphan_scratch, the mkdtemp suffix so a recycled pid cannot inherit a dead run's dir.
        reap_orphan_scratch(save_dir_root);
        std::string tmpl = save_dir_root + "/nopt_dmrg_" + std::to_string((long)getpid()) + "_XXXXXX";
        std::vector<char> tbuf(tmpl.begin(), tmpl.end());
        tbuf.push_back('\0');
        if (mkdtemp(tbuf.data()) == nullptr) {
            fprintf(out_stream, "ERROR: cannot create DMRG scratch under %s: %s\n",
                    save_dir_root.c_str(), strerror(errno));
            exit(EXIT_FAILURE);
        }
        scratch = tbuf.data();

        Random::rand_seed(0);
        // isize/dsize are BYTE sizes of the integer/double stacks. The double stack holds the
        // renormalized operators and is sized by $DMRG memory (GB); the integer stack is bookkeeping.
        frame_<double>() = std::make_shared<DataFrame<double>>(
            (size_t)1 << 24, (size_t)(memory_gb * (double)((size_t)1 << 30)), scratch);
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
void ensure_block2_runtime(const std::string &save_dir_root, double memory_gb, int n_threads) {
    static Block2Runtime runtime(save_dir_root, memory_gb, n_threads);
    (void)runtime;
}

// Best-effort removal of every scratch file of one MPS tag; keeps the scratch dir from growing
// across macro-iterations. block2 embeds the tag, trailing-'.'-delimited, in
// <prefix>.{MPS,MMPS,MMPS-WFN,MPS.INFO}.<tag>.* and <tag>-mps_info.bin. The trailing '.' keeps tag
// "work_1" from matching "work_10".
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

} // namespace nopt_block2
using namespace nopt_block2;

// struct dmrgci_engine (all block2 state) + host_threads_guard live in block2_dmrg_engine.h.

// Exact energy of the stored MPS for one root: that root's 2-RDM (and its 1-RDM partial trace)
// contracted with the integrals the solver used. block2's sweep energy is the Davidson eigenvalue of
// the untruncated center, not the energy of the stored bond-M MPS. `d2` must still be in the
// solver's lattice-ordered, localized basis, i.e. before ensure_2rdm un-permutes and delocalizes it.
static double rdm_energy(const dmrgci_engine &e, const double *d2) {
    const int n = e.n_act;
    const double inv = 1.0 / (e.n_elec - 1);
    double e1 = 0.0, e2 = 0.0;
    for (int p = 0; p < n; p++)
        for (int s = 0; s < n; s++) {
            double d1 = 0.0; // D1[p,s] = 1/(N-1) sum_k D2[p,k,k,s]
            for (int k = 0; k < n; k++)
                d1 += d2[(((size_t)p * n + k) * n + k) * n + s];
            e1 += e.fcidump->t(p, s) * inv * d1;
        }
    for (int p = 0; p < n; p++)
        for (int q = 0; q < n; q++)
            for (int r = 0; r < n; r++)
                for (int s = 0; s < n; s++)
                    e2 += e.fcidump->v(p, s, q, r) *
                          d2[(((size_t)p * n + q) * n + r) * n + s];
    return (double)e.fcidump->e() + e1 + 0.5 * e2;
}

// One root of the state-averaged MultiMPS as a plain single-root MPS, ready for an Expect sweep.
// Expect dispatches on the canonical form, and its MultiMPS branch resolves Automatic to Normal, a
// serial loop over operators; only a plain MPS reaches the OpenMP-parallel path. The center is then
// put in the one-site form (block2's site_type=1 conversion), which needs it at an end.
static std::shared_ptr<MPS<SU2, double>> extract_root_single(dmrgci_engine &e, int st,
                                                             const std::string &xtag,
                                                             const std::string &stag) {
    std::shared_ptr<MultiMPS<SU2, double>> xmps = e.mps->extract(st, xtag);
    std::shared_ptr<MPS<SU2, double>> imps = xmps->make_single(stag);
    if (imps->dot == 2) {
        imps->dot = 1;
        if (imps->center == imps->n_sites - 2) {
            imps->center = imps->n_sites - 1;
            imps->canonical_form[imps->n_sites - 1] = 'S';
        } else if (imps->center == 0)
            imps->canonical_form[0] = 'K';
        else {
            fprintf(out_stream, "ERROR: RDM read-out expects the MPS center at an end "
                                "(got %d of %d sites)\n", imps->center, imps->n_sites);
            exit(EXIT_FAILURE);
        }
        imps->save_data();
    }
    return imps;
}

// Spatial 2-RDMs, once per solve, one root at a time: one Expect sweep per extracted root yields
// that root's 2-RDM, which gives its energy and its 1-RDM (a partial trace), is folded into the
// running state average, and is dropped. Only the average survives, so the peak is two n_act^4
// tensors and the un-permutation and back-transform run once. A PDM sweep straight on the MultiMPS
// would truncate the shared basis when it moves the center; a single-root MPS moves losslessly.
static void ensure_2rdm(dmrgci_engine &e) {
    if (e.d2_valid)
        return;
    host_threads_guard htg;
    const int n = e.n_act;
    const size_t blk = (size_t)n * n * n * n;
    const size_t blk1 = (size_t)n * n;
    e.d2_av.assign(blk, 0.0);
    e.d1_states.assign(blk1 * e.n_s, 0.0);

    // Weights the host optimizes under; absent a weight vector the roots are equally weighted.
    std::vector<double> w(e.n_s, 1.0);
    if ((int)e.w_state.size() == e.n_s)
        w = e.w_state;
    double wsum = 0.0;
    for (int st = 0; st < e.n_s; st++)
        wsum += w[st];

    // State-independent (a function of e.hamil alone): built once, the environment rebound per root.
    std::shared_ptr<MPO<SU2, double>> p2mpo = std::make_shared<PDM2MPOQC<SU2, double>>(e.hamil);
    p2mpo = std::make_shared<SimplifiedMPO<SU2, double>>(
        p2mpo, std::make_shared<RuleQC<SU2, double>>(), true, true,
        OpNamesSet({OpNames::R, OpNames::RD}));

    const double inv = (e.n_elec >= 2) ? 1.0 / (e.n_elec - 1) : 0.0;

    for (int st = 0; st < e.n_s; st++) {
        const std::string xtag = e.mps_info->tag + "-" + std::to_string(st);
        const std::string stag = xtag + "-s";
        std::shared_ptr<MPS<SU2, double>> imps = extract_root_single(e, st, xtag, stag);

        auto p2me = std::make_shared<MovingEnvironment<SU2, double, double>>(p2mpo, imps, imps,
                                                                             "2PDM");
        p2me->init_environments(false);
        auto ex2 = std::make_shared<Expect<SU2, double, double>>(p2me, (ubond_t)e.cfg.m,
                                                                 (ubond_t)e.cfg.m);
        ex2->iprint = 0; // silence the per-site Expect sweep log
        ex2->solve(true, imps->center == 0);
        std::shared_ptr<GTensor<double>> d2 = ex2->get_2pdm_spatial(); // shape {n,n,n,n}

        // Contiguous row-major [p,q,r,s] = the block2 D2 layout. This root's own buffer, read in place.
        const double *d2p = d2->data->data();

        // This root's true energy, not the solver's pre-truncation sweep value. Taken while the
        // 2-RDM is still in the solver's basis and lattice order, the one e.fcidump is in.
        if (e.n_elec >= 2)
            e.E_states[st] = rdm_energy(e, d2p);

        // This root's 1-RDM, D1[p,s] = 1/(N-1) sum_k D2[p,k,k,s]. The trace commutes with the
        // orthogonal un-permutation and back-transform below, which the n_act^2 matrix carries.
        double *d1 = e.d1_states.data() + (size_t)st * blk1;
        for (int p = 0; p < n; p++)
            for (int s = 0; s < n; s++) {
                double v = 0.0;
                for (int k = 0; k < n; k++)
                    v += d2p[(((size_t)p * n + k) * n + k) * n + s];
                d1[(size_t)p * n + s] = inv * v;
            }

        const double ws = w[st];
#pragma omp parallel for schedule(static)
        for (size_t k = 0; k < blk; k++)
            e.d2_av[k] += ws * d2p[k];

        p2me->remove_partition_files();
        remove_tag_files(xtag); // the per-root extract and its single-MPS copy are transient
        remove_tag_files(stag);
    }
    p2mpo->deallocate();

#pragma omp parallel for schedule(static)
    for (size_t k = 0; k < blk; k++)
        e.d2_av[k] /= wsum;

    if (!e.reorder_perm.empty()) { // map out of block2's Fiedler lattice order
        std::vector<int> iperm(n);                 // inverse of reorder_perm
        for (int i = 0; i < n; i++) iperm[e.reorder_perm[i]] = i;
        std::vector<double> ro(blk);
        for (int p = 0; p < n; p++)
            for (int q = 0; q < n; q++)
                for (int r = 0; r < n; r++)
                    for (int s = 0; s < n; s++)
                        ro[(((size_t)p * n + q) * n + r) * n + s] =
                            e.d2_av[(((size_t)iperm[p] * n + iperm[q]) * n + iperm[r]) * n +
                                    iperm[s]];
        e.d2_av.swap(ro);
        std::vector<double> ro1(blk1);
        for (int st = 0; st < e.n_s; st++) {
            double *d1 = e.d1_states.data() + (size_t)st * blk1;
            for (int p = 0; p < n; p++)
                for (int q = 0; q < n; q++)
                    ro1[(size_t)p * n + q] = d1[(size_t)iperm[p] * n + iperm[q]];
            std::copy(ro1.begin(), ro1.end(), d1);
        }
    }

    if (e.localize_on) { // rotate back to the delocalized basis (rotate1/rotate2 forbid aliasing)
        std::vector<double> bt(blk);
        rotate2(e.d2_av.data(), e.U_loc.data(), n, bt.data(), /*forward=*/false);
        e.d2_av.swap(bt);
        std::vector<double> bt1(blk1);
        for (int st = 0; st < e.n_s; st++) {
            double *d1 = e.d1_states.data() + (size_t)st * blk1;
            rotate1(d1, e.U_loc.data(), n, bt1.data(), /*forward=*/false);
            std::copy(bt1.begin(), bt1.end(), d1);
        }
    }

    for (int st = 0; st < e.n_s; st++) { // symmetrize (finite-M DMRG mildly breaks it)
        double *d1 = e.d1_states.data() + (size_t)st * blk1;
        for (int p = 0; p < n; p++)
            for (int q = p + 1; q < n; q++) {
                double a = 0.5 * (d1[p * n + q] + d1[q * n + p]);
                d1[p * n + q] = d1[q * n + p] = a;
            }
    }

    e.d2_valid = true;
    assert_stack_clean("2-RDM read"); // the Expect sweeps must leave the LIFO stacks as they found them
}

// Full n_s x n_s spin-summed 1-RDM for the property read-out, once per solve, in the delocalized
// basis: one direct 1-PDM Expect sweep per block, never a 2-RDM partial trace. Diagonal i==j is a
// per-state 1-RDM (bra==ket); i<j is the bra!=ket transition 1-RDM, with (j,i) its transpose (exact
// against the symmetric property integrals).
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

    // State-pair-independent: built once, the environment rebound per pair. NoTransposeRule because
    // the transpose-symmetry simplification is invalid when bra != ket.
    std::shared_ptr<MPO<SU2, double>> p1mpo = std::make_shared<PDM1MPOQC<SU2, double>>(e.hamil);
    p1mpo = std::make_shared<SimplifiedMPO<SU2, double>>(
        p1mpo,
        std::make_shared<NoTransposeRule<SU2, double>>(std::make_shared<RuleQC<SU2, double>>()),
        true, true, OpNamesSet({OpNames::R, OpNames::RD}));

    for (int i = 0; i < e.n_s; i++)
        for (int j = i; j < e.n_s; j++) {
            // Fresh extract of both roots so bra/ket share a canonical center. The single-MPS form is
            // all-or-nothing: the effective Hamiltonian asserts bra and ket agree.
            const std::string itag = e.mps_info->tag + "-t" + std::to_string(i);
            const std::string istag = itag + "-s";
            std::shared_ptr<MPS<SU2, double>> imps = extract_root_single(e, i, itag, istag);
            std::string jtag, jstag;
            std::shared_ptr<MPS<SU2, double>> jmps = imps;
            if (j != i) {
                jtag = e.mps_info->tag + "-t" + std::to_string(j);
                jstag = jtag + "-s";
                jmps = extract_root_single(e, j, jtag, jstag);
            }

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
            remove_tag_files(itag); // the per-root extracts and their single-MPS copies are transient
            remove_tag_files(istag);
            if (j != i) {
                remove_tag_files(jtag);
                remove_tag_files(jstag);
            }
        }
    p1mpo->deallocate();
    e.dmfull_valid = true;
    assert_stack_clean("full 1-RDM read");
}

// Reload the retained MultiMPS fresh from disk for a warm restart. Post-solve the in-memory MultiMPS
// is a partially-deallocated shell whose StateInfos dangle into freed stack memory; the on-disk copy
// is authoritative. load_data/load_mutable reallocate every StateInfo/SparseMatrixInfo on the heap,
// rebuilt against the current MPO basis.
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

// Rotate the retained MultiMPS from the previous macro-iteration's active basis into the current one.
// R (e.R_active, order [a*n+p]) is the overlap of the two localized active-MO sets: regularize to the
// nearest proper rotation U~, take kappa = log(U~), apply exp(kappa) by RK4 time evolution. Only the
// MPS is rotated; the new-basis integrals are already imported. e.mps is mutated in place on success.
// Returns false (=> caller cold-starts) on excess leakage, an improper generator, or norm drift.
static bool rotate_retained_mps(dmrgci_engine &e) {
    const int n = e.n_act;
    const size_t nn = (size_t)n * n;

    // These temporaries MUST live on the heap, not block2's LIFO stack: the stack is empty here, so
    // this math and the subsequent init_environments allocate from the same base region, where block2
    // has a latent uninitialized read -- leftover bytes there cause a data-dependent segfault.
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
        fprintf(out_stream, " warm-start: rotation leakage %.2e too large -> cold fallback                |\n", leak);
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

    // exp(-kappa) rotation of the retained MPS (kappa = log U~), shared with the read-out. One step
    // reproduces the full rotation for small per-iter rotations; rot_steps guards larger ones.
    const int rot_m = e.cfg.rot_m > 0 ? e.cfg.rot_m : e.cfg.m;
    const int n_steps = e.cfg.rot_steps > 0 ? e.cfg.rot_steps : 10;
    auto res = apply_orbital_rotation_mps(e.mps, U.data(), n, e.n_elec, e.twos, e.orbsym,
                                          e.reorder_perm, rot_m, n_steps);
    if (res.complex_generator) { // non-real generator -> not a proper rotation
        fprintf(out_stream, "  warm-start: complex generator (|Im k|=%.2e) -> cold fallback\n",
                res.im_norm);
        return false;
    }
    if (res.skipped)
        return true; // negligible rotation (near convergence): reuse the reloaded MPS as-is
    // Unitary propagation should preserve per-root norms. Negated rather than `> 1e-3` so a non-finite
    // norm falls back too: every comparison with NaN is false, so `fabs(NaN) > 1e-3` would accept it.
    if (!(std::fabs(res.norm2 - 1.0) <= 1e-3)) {
        fprintf(out_stream, "  warm-start: rotation norm^2=%.6f drifted -> cold fallback\n", res.norm2);
        return false;
    }
    return true; // rotated MPS in place; benign success is silent (keeps the CASSCF table clean)
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

    // DMRG lattice order: block2 strings orbitals on the lattice in raw FCIDUMP order, so localized
    // orbitals come out scrambled and are ordered by Fiedler on the exchange matrix K_ij = |(ij|ji)|.
    // RDMs come back in this order and are un-permuted in ensure_2rdm. A warm solve reuses the
    // retained MPS's order -- Fiedler order is a function of the localized orbitals, so it is pinned
    // together with the frozen localization. A cold solve recomputes it.
    if (e.have_rotation && !e.reorder_perm.empty()) {
        e.fcidump->reorder(e.reorder_perm); // frozen order (warm restart)
    } else {
        e.reorder_perm.clear();
        if (e.cfg.loc_order == DMRG_LOCORDER_FIEDLER) {
            // block2's metric (pyblock2 parser.py): the exchange graph, with the one-electron
            // coupling as a tie-break so a disconnected or tied graph still orders reproducibly.
            std::vector<double> kmat((size_t)n * n, 0.0);
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    if (i != j)
                        kmat[(size_t)i * n + j] =
                            std::fabs(h2[(((size_t)i * n + j) * n + j) * n + i]) +
                            1e-7 * std::fabs(h1[(size_t)i * n + j]);
            e.reorder_perm = OrbitalOrdering::fiedler((uint16_t)n, kmat);
            e.fcidump->reorder(e.reorder_perm);
        }
    }

    SU2 vacuum(0);
    e.hamil = std::make_shared<HamiltonianQC<SU2, double>>(vacuum, n, e.orbsym, e.fcidump);
    e.hamil->opf->seq->mode = SeqTypes::Tasked;
    e.mpo = build_qc_mpo(e.hamil);
}

// One-site sweeps closing every solve. Mirrors block2's default schedule, which ends two sweeps
// into the noise-free stage in one-site mode (pyblock2/driver/parser.py sets twodot_to_onedot =
// last_iter + 2 unless the input overrides it).
static const int DMRG_ONEDOT_TAIL = 2;

// Put the MPS back into a two-site center after a one-site tail: a one-site sweep ends with a fused
// single-site center ('J'/'T') at a boundary, and every downstream consumer builds its layout from
// canonical_form/center/dot, which must agree. Port of block2's DMRGDriver::adjust_mps(dot=2).
static void adjust_mps_two_dot(dmrgci_engine &e) {
    auto &mps = e.mps;
    auto cg = e.mpo->tf->opf->cg;
    const int n = mps->n_sites;
    auto flip = [&](int c) { // fused-left <-> fused-right; refresh the in-memory copy
        mps->flip_fused_form(c, cg, nullptr);
        mps->save_data();
        mps->load_mutable();
        mps->info->load_mutable();
    };
    if (mps->center == 0 && (mps->canonical_form[0] == 'S' || mps->canonical_form[0] == 'T'))
        flip(0); // S->K, T->J
    mps->dot = 2;
    const char cf = mps->canonical_form[mps->center];
    if (cf == 'L' && mps->center != n - 2)
        mps->center += 1;
    else if ((cf == 'C' || cf == 'M') && mps->center != 0)
        mps->center -= 1;
    if (mps->center == n - 1) {
        if (mps->canonical_form[mps->center] == 'K' || mps->canonical_form[mps->center] == 'J')
            flip(mps->center); // K->S, J->T
        mps->center = n - 2;
    }
    mps->save_data();
}

// Re-order the lattice for a cold fallback: import_integrals commits the order before solve() may
// decline the warm rotation, leaving the FCIDUMP on the previous orbitals' Fiedler order. The
// already-permuted FCIDUMP is re-ordered in place -- Fiedler on a relabelled exchange graph is the
// same ordering up to that relabelling, so the fresh permutation composes with the frozen one.
static void recompute_cold_order(dmrgci_engine &e) {
    if (e.cfg.loc_order != DMRG_LOCORDER_FIEDLER || e.reorder_perm.empty())
        return;
    const int n = e.n_act;
    std::vector<double> kmat((size_t)n * n, 0.0);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (i != j)
                kmat[(size_t)i * n + j] = std::fabs(e.fcidump->v(i, j, j, i)) +
                                          1e-7 * std::fabs(e.fcidump->t(i, j));
    std::vector<uint16_t> p2 = OrbitalOrdering::fiedler((uint16_t)n, kmat);
    bool unchanged = true;
    for (int k = 0; k < n && unchanged; k++)
        unchanged = (p2[k] == (uint16_t)k);
    if (unchanged)
        return; // the frozen order already is the fresh one: nothing to rebuild
    e.fcidump->reorder(p2);
    std::vector<uint16_t> composed((size_t)n); // site k now carries the orbital that sat at site p2[k]
    for (int k = 0; k < n; k++)
        composed[k] = e.reorder_perm[p2[k]];
    e.reorder_perm.swap(composed);

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
    // Same condition import_integrals keyed the frozen lattice order on, before the rotation could decline.
    const bool order_frozen = (e.have_rotation && !e.reorder_perm.empty());
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
    if (!warm && order_frozen)
        recompute_cold_order(e); // the order was pinned to a basis we are no longer in

    // --- sweep schedule: short warm re-solve vs full cold ramp ---
    dmrg_schedule sch;
    if (warm) {
        sch = build_warm_schedule(e.cfg.m, e.cfg.warm_sweeps, e.cfg.sweep_tol);
    } else if (e.cfg.schedule == DMRG_SCHED_DEFAULT) {
        sch = build_default_schedule(e.cfg.m, e.cfg.sweeps, e.cfg.sweep_tol);
    } else {
        fprintf(out_stream,
                "ERROR: DMRG schedule option not implemented yet (only 'default')\n");
        exit(EXIT_FAILURE);
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
            exit(EXIT_FAILURE);
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

    // Convergence of the variational (two-site) phase, taken before the tail appends to the same
    // history. Residual = max over roots of the last two sweeps' |dE|; a mean would cancel
    // opposite-moving roots. A single sweep leaves it unmeasurable (NaN) rather than faking
    // convergence.
    const int n2 = (int)dmrg->energies.size(); // two-site sweeps actually run (may early-stop)
    if (n2 >= 2) {
        const auto &e1 = dmrg->energies[n2 - 1];
        const auto &e0 = dmrg->energies[n2 - 2];
        const int nr = (int)e1.size() < (int)e0.size() ? (int)e1.size() : (int)e0.size();
        double dmax = 0.0;
        for (int r = 0; r < nr; r++) {
            double d = std::fabs((double)e1[r] - (double)e0[r]);
            if (d > dmax) dmax = d;
        }
        e.last_sweep_dE = dmax;
    } else {
        e.last_sweep_dE = std::numeric_limits<double>::quiet_NaN();
    }
    // Exhausted the sweep budget without block2's sweep_tol early-stop -> possibly under-converged
    // (flagged in the CAS-SCF table). An unmeasurable one-sweep residual (NaN) also flags.
    e.last_hit_max = (n2 >= sch.n_sweeps) &&
                     (std::isnan(e.last_sweep_dE) || e.last_sweep_dE > e.cfg.sweep_tol);

    // Close with one-site sweeps at zero noise, tolerance zero. The center move is an exact gauge
    // transformation only for a single root (density-matrix rank <= M); SA roots share a renormalized
    // basis, so for n_s > 1 the tail truncates. Zero noise is a correctness condition, not tuning:
    // block2's perturbative noise raises the density-matrix rank. Tolerance zero because block2
    // compares against the previous history entry, a two-site value measuring a different quantity.
    if (DMRG_ONEDOT_TAIL > 0) {
        const int total = n2 + DMRG_ONEDOT_TAIL;
        dmrg->me->dot = 1;
        dmrg->bond_dims.assign(total, (ubond_t)e.cfg.m);
        dmrg->noises.assign(total, 0.0);
        dmrg->davidson_conv_thrds.assign(total, sch.dav_thrds.back());
        dmrg->solve(total, dmrg->forward, /*tol=*/0.0, n2);
        e.mps->dot = 1; // the sweep switched me->dot; the MPS must say so too, or every consumer
        e.mps->save_data(); // that reads canonical_form/center/dot builds a two-site layout on it
        adjust_mps_two_dot(e); // ... and back to a two-site center for those consumers
    }

    // per-root energies (ascending; root 0 = ground state). block2's energy precision is FPLS
    // (long double for FL=double), so bind via auto and narrow to NOPT's double.
    if (dmrg->energies.empty()) {
        fprintf(out_stream, "ERROR: DMRG produced no sweep energies\n");
        exit(EXIT_FAILURE);
    }
    const auto &eng = dmrg->energies.back();
    if ((int)eng.size() < e.n_s) {
        fprintf(out_stream, "ERROR: DMRG returned %d roots < n_s=%d (requested state count likely"
                            " exceeds the sector dimension)\n", (int)eng.size(), e.n_s);
        exit(EXIT_FAILURE);
    }
    e.E_states.assign(e.n_s, 0.0);
    for (int s = 0; s < e.n_s; s++) {
        if (!std::isfinite((double)eng[s])) {
            fprintf(out_stream, "ERROR: DMRG root %d energy is not finite\n", s);
            exit(EXIT_FAILURE);
        }
        e.E_states[s] = (double)eng[s];
    }

    e.last_n_sweeps = (int)dmrg->energies.size(); // two-site + tail: what this solve actually cost

    me->remove_partition_files(); // keep mps/mps_info alive for the RDM read-out
    return e.last_n_sweeps;
}
void block2_casci_wrap::calc_DM_diag(double *gamma, int /*i_set*/) {
    dmrgci_engine &e = *impl_;
    if (e.n_elec < 2) { // the 2-RDM trace cannot recover the 1-RDM for <2 electrons
        fprintf(out_stream,
                "ERROR: DMRG 1-RDM via 2-RDM trace needs N_elec>=2 (got %d)\n", e.n_elec);
        exit(EXIT_FAILURE);
    }
    ensure_2rdm(e);
    // Per state: the 1-RDM is the partial trace of that state's spatial 2-RDM,
    //   D1[p,s] = 1/(N-1) * sum_k D2[p,k,k,s]   (block2 D2 = <a+_p a+_q a_r a_s>),
    // contracted root by root as the 2-RDMs are read. gamma holds n_s consecutive n_act^2 blocks
    // (CAS_engine averages them with the weights).
    std::copy(e.d1_states.begin(), e.d1_states.end(), gamma);
}
void block2_casci_wrap::G_calc(double *GAMMA) {
    dmrgci_engine &e = *impl_;
    const int n = e.n_act;
    ensure_2rdm(e);
    // The orbital gradient consumes the state-averaged 2-RDM, and that is the only one this backend
    // forms (the per-state tensors never coexist), so GAMMA is a single n_act^4 block:
    // NOPT layout GAMMA[p,q,r,s] = block2 D2[p,r,s,q].
    const double *d2 = e.d2_av.data();
    for (int p = 0; p < n; p++)
        for (int q = 0; q < n; q++)
            for (int r = 0; r < n; r++)
                for (int s = 0; s < n; s++)
                    GAMMA[(((size_t)p * n + q) * n + r) * n + s] =
                        d2[(((size_t)p * n + r) * n + s) * n + q];
}
void block2_casci_wrap::calc_DMA(double *gamma, int a, int b) {
    // Properties (dipole/quadrupole/Mulliken/transition dipoles). block2 has a single SA state set,
    // so the host only ever asks for the full n_s x n_s block, gamma[(bra*n_s+ket)*n_act^2].
    if (a != 0 || b != 0) {
        fprintf(out_stream, "ERROR: DMRG backend expects a single state set (got a=%d b=%d)\n", a, b);
        exit(EXIT_FAILURE);
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
    // No-op: SU2 is spin-adapted, so calc_DMA carries the full spin-summed 1-RDM and the host sums
    // DMA+DMB. WARNING: the A/B split is NOT spin-resolved (all density in A, none in B), so this is
    // valid only for consumers that use the sum. One reading DMA and DMB separately (XMCQDPT2/CDAS PT
    // tensors) needs a real SU2 spin-density (na-nb sector) 1-RDM here first.
}
// Compress a single MPS to bond dimension target_m by fitting a small-m MPS to it through the
// identity MPO (block2 Linear "Normal" equation). State-preserving projection, used only for the
// read-out, where it makes the DeterminantTRIE search cheaper (cost ~ m^2) while preserving the
// leading determinants. ctag names the fitted MPS's scratch; the caller removes it.
std::shared_ptr<MPS<SU2, double>>
nopt_block2::compress_single_mps(dmrgci_engine &e, const std::shared_ptr<MPS<SU2, double>> &ket,
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

