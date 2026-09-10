// block2_casci_wrap density read-outs over the external block2 DMRG library: the general-NPDM
// sweep every reduced density matrix runs through (npdm_lattice), the per-state spin-summed 2-RDM
// (G2_calc_diag) and the 3-body moment (G3_calc_diag). Split from block2_dmrg.cpp, same idiom.

#include "block2_dmrg_engine.h"   // block2 headers + dmrgci_engine + shared helpers
#include "block2_gpu_guard.h" // $DMRG gpu=on: the block2 GPU backend for the sweeps

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <memory>
#include <string>
#include <system_error>
#include <vector>

#include "common_vars.h"      // out_stream
#include "dmrg_log.h"         // DMRG_LOG_IPRINT, dmrg_log_guard
#include "tensor_rotate.h"    // rotate2/rotate3 (active-space basis back-transform)

using namespace block2;
using namespace nopt_block2;

// ---- 2-RDM read-out helpers ----------------------------------------------------------------

// Finish one raw block2 2-RDM block into the NOPT GAMMA convention: un-permute out of the
// Fiedler lattice, rotate back to the delocalized basis, then GAMMA[p,q,r,s] = D2[p,r,s,q].
// perm_scr / rot_scr are caller-owned n_act^4 scratch (rotate2 forbids aliasing).
static void g2full_finish_block(const dmrgci_engine &e, const double *raw,
                                const std::vector<int> &iperm, double *perm_scr,
                                double *rot_scr, double *out) {
    const int n = e.n_act;
    const double *cur = raw;
    if (!iperm.empty()) { // map out of block2's Fiedler lattice order
        for (int p = 0; p < n; p++)
            for (int q = 0; q < n; q++)
                for (int r = 0; r < n; r++)
                    for (int s = 0; s < n; s++)
                        perm_scr[(((size_t)p * n + q) * n + r) * n + s] =
                            raw[(((size_t)iperm[p] * n + iperm[q]) * n + iperm[r]) * n + iperm[s]];
        cur = perm_scr;
    }
    if (e.localize_on) { // rotate back to the delocalized basis
        rotate2(cur, e.U_loc.data(), n, rot_scr, /*forward=*/false);
        cur = rot_scr;
    }
    for (int p = 0; p < n; p++)
        for (int q = 0; q < n; q++)
            for (int r = 0; r < n; r++)
                for (int s = 0; s < n; s++)
                    out[(((size_t)p * n + q) * n + r) * n + s] =
                        cur[(((size_t)p * n + r) * n + s) * n + q];
}

// The lattice -> input orbital map of every read-out: inverse of reorder_perm, empty when the
// solve ran in the input order.
static std::vector<int> lattice_iperm(const dmrgci_engine &e) {
    std::vector<int> iperm;
    if (!e.reorder_perm.empty()) {
        iperm.resize(e.n_act);
        for (int i = 0; i < e.n_act; i++) iperm[e.reorder_perm[i]] = i;
    }
    return iperm;
}

// A converged MPS is the precondition of every read-out below.
static void require_solved(const dmrgci_engine &e, const char *what) {
    if (e.mps == nullptr || e.mps_info == nullptr) {
        fprintf(out_stream, "ERROR: DMRG %s needs a converged MPS (call solve first)\n", what);
        exit(EXIT_FAILURE);
    }
}

// ---- the general NPDM read-out every density matrix runs through ---------------------------

// NOPT_RDM_DUMP=<dir>: the raw lattice-order tensor, one file per root pair and body order.
static void dump_lattice(const GTensor<double> &t, int N, int ket_state, int bra_state) {
    const char *dir = std::getenv("NOPT_RDM_DUMP");
    if (dir == nullptr || *dir == '\0')
        return;
    std::error_code ec;
    std::filesystem::create_directories(dir, ec);
    const std::string path = std::string(dir) + "/N" + std::to_string(N) + "_ket" +
                             std::to_string(ket_state) + "_bra" + std::to_string(bra_state) +
                             ".bin";
    FILE *f = fopen(path.c_str(), "wb");
    if (f == nullptr) {
        fprintf(out_stream, "ERROR: DMRG RDM dump cannot open %s\n", path.c_str());
        exit(EXIT_FAILURE);
    }
    fwrite(t.data->data(), sizeof(double), t.size(), f);
    fclose(f);
}

// One process plays every rank in turn, so block2's cross-rank sum of a pass's NPDM fragment has
// nothing to add: each pass keeps its own share and the caller sums the passes. Every other
// collective keeps the base class's single-rank assertion.
struct npdm_pass_comm : ParallelCommunicator<SU2> {
    npdm_pass_comm(int size, int rank) : ParallelCommunicator<SU2>(size, rank, 0) {}
    void allreduce_sum(double *, size_t) override {}
};

// block2's ParallelRule constructor repoints the distributed scratch prefix at the rank. The
// serial prefix must be restored on every exit from the pass loop, exceptional ones included.
struct prefix_guard {
    std::string prefix;
    bool can_write;
    prefix_guard()
        : prefix(frame_<double>()->prefix_distri),
          can_write(frame_<double>()->prefix_can_write) {}
    ~prefix_guard() {
        frame_<double>()->prefix_distri = prefix;
        frame_<double>()->prefix_can_write = can_write;
    }
    prefix_guard(const prefix_guard &) = delete;
    prefix_guard &operator=(const prefix_guard &) = delete;
};

// SU2 recoupling string of the spin-summed N-body density: "(C+D)0" wrapped N-1 times. block2
// returns the singlet-coupled raw[x0..x_{N-1},y0..y_{N-1}] = sum_spins <a+_x0 .. a_y0> scaled by
// 2^{-N/2}.
static std::string npdm_expr(int N) {
    std::string s = "(C+D)0";
    for (int k = 1; k < N; k++)
        s = "((C+" + s + ")1+D)0";
    return s;
}

// block2 streams the NPDM middle intermediates to <save_dir>/*.NPDM.FRAG.* and leaves them there;
// the Compressed algorithm writes the .fpc spelling, the plain one .npy.
static void remove_npdm_fragments(const MovingEnvironment<SU2, double, double> &me) {
    std::error_code ec;
    for (int i = 0; i < me.n_sites; i++) {
        const std::string base = me.get_npdm_fragment_filename(i);
        std::filesystem::remove(base + ".fpc", ec);
        std::filesystem::remove(base + ".npy", ec);
    }
}

// Canonicalize a one-site extract to the end dir asks for (+1 = center 0 'K', -1 = center
// n_sites-1 'S'); dir 0 keeps the end the solve left. The split runs at the MPS's own maximum
// bond dimension, so it is lossless. The walked site tensors are dropped again: they are on disk.
static void move_to_end(const std::shared_ptr<MPS<SU2, double>> &mps, int dir) {
    if (dir == 0)
        return;
    const int want = dir > 0 ? 0 : mps->n_sites - 1;
    if (mps->center == want)
        return;
    auto cg = std::make_shared<CG<SU2>>();
    mps->info->load_mutable();
    mps->info->bond_dim = std::max(mps->info->bond_dim, mps->info->get_max_bond_dimension());
    while (mps->center != want)
        dir > 0 ? mps->move_left(cg) : mps->move_right(cg);
    mps->save_data();
    mps->info->save_mutable();
    mps->info->deallocate_mutable();
    for (int i = 0; i < mps->n_sites; i++)
        if (mps->tensors[i] != nullptr && mps->tensors[i]->total_memory != 0)
            mps->unload_tensor(i);
    if (mps->center != want || mps->canonical_form[want] != (dir > 0 ? 'K' : 'S')) {
        fprintf(out_stream, "ERROR: DMRG RDM center move ended at center %d form '%c'\n",
                mps->center, mps->canonical_form[want]);
        exit(EXIT_FAILURE);
    }
}

// Sweep direction of every RDM read-out: NOPT_RDM_SWEEP = auto | forward | backward. auto keeps
// the end the solve left, forward sweeps from site 0, backward from the last site.
int nopt_block2::rdm_sweep_dir() {
    const char *v = std::getenv("NOPT_RDM_SWEEP");
    if (v == nullptr)
        return 0;
    return std::strcmp(v, "forward") == 0 ? 1 : (std::strcmp(v, "backward") == 0 ? -1 : 0);
}

// One root pair's spin-summed N-body density in block2's lattice order, from one general-NPDM
// Expect sweep on transient single-root extracts. The result is unscaled (block2's convention);
// callers apply sqrt(2)^N and their own gathers. tag names the environment and the scratch.
std::shared_ptr<GTensor<double>> nopt_block2::npdm_lattice(dmrgci_engine &e, int N, int ket_state,
                                                           int bra_state, int dir,
                                                           const char *tag) {
    require_solved(e, tag);
    const int n = e.n_act;
    size_t nel = 1;
    for (int k = 0; k < 2 * N; k++) nel *= (size_t)n;
    // the divider pays only for the 3-RDM: the lower orders buy no memory with it
    const int passes = (N == 3 && e.cfg.rdm_passes > 1) ? e.cfg.rdm_passes : 1;

    const std::string ktag = e.mps_info->tag + "-" + tag + std::to_string(ket_state);
    const std::string kstag = ktag + "-s";
    const std::string btag = e.mps_info->tag + "-" + tag + "b" + std::to_string(bra_state);
    const std::string bstag = btag + "-s";

    std::shared_ptr<GTensor<double>> acc; // {n}^{2N}, lattice order, summed over the passes
    prefix_guard pfx;

    for (int r = 0; r < passes; r++) {
        // One pass carries one rank's share of the operator set. The rank splits the MPO's left
        // families and its longest right strings; the passes sum to the full density.
        std::shared_ptr<ParallelRuleSimple<SU2, double>> sp_rule;
        if (passes > 1)
            sp_rule = std::make_shared<ParallelRuleSimple<SU2, double>>(
                ParallelSimpleTypes::None,
                std::make_shared<npdm_pass_comm>(passes, r));

        // The GeneralHamiltonian is built fresh for each MPO: its on-site operator tables are
        // populated on first use, and reusing one instance corrupts every operator carrying
        // coincident legs.
        SU2 vacuum(0);
        std::vector<typename SU2::pg_t> gorbsym(n, 0); // C1 site irreps
        auto ghamil = std::make_shared<GeneralHamiltonian<SU2, double>>(vacuum, n, gorbsym);
        const std::string expr = npdm_expr(N);
        auto perm = std::make_shared<SpinPermScheme>(
            SpinPermScheme::initialize_su2(2 * N, expr, /*is_npdm=*/true));
        auto ppmpo = std::make_shared<GeneralNPDMMPO<SU2, double>>(
            ghamil, std::make_shared<NPDMScheme>(perm), /*symbol_free=*/true, 0.0, 0,
            "NPDM" + std::to_string(N));
        ppmpo->delta_quantum = SU2(0, SpinPermRecoupling::get_target_twos(expr), 0);
        ppmpo->iprint = DMRG_LOG_IPRINT >= 2 ? 1 : 0; // per-site operator counts into the sweep log
        ppmpo->parallel_rule = sp_rule; // must be set before build(): it sizes the rank's blocks
        ppmpo->build();
        std::shared_ptr<MPO<SU2, double>> pmpo = std::make_shared<SimplifiedMPO<SU2, double>>(
            ppmpo, std::make_shared<Rule<SU2, double>>(), false, false);
        if (sp_rule != nullptr)
            pmpo = std::make_shared<ParallelMPO<SU2, double>>(pmpo, sp_rule);

        {
            std::shared_ptr<MPS<SU2, double>> ket = extract_root_single(e, ket_state, ktag, kstag);
            std::shared_ptr<MPS<SU2, double>> bra = ket;
            if (bra_state != ket_state)
                bra = extract_root_single(e, bra_state, btag, bstag);
            move_to_end(ket, dir);
            if (bra != ket)
                move_to_end(bra, dir);

            auto me = std::make_shared<MovingEnvironment<SU2, double, double>>(pmpo, bra, ket, tag);
            me->cached_contraction = false; // conflicts with the fused zero-dot contraction
            me->fused_contraction_rotation = true;
            block2_gpu_guard g(e.cfg.gpu, me);
            me->init_environments(DMRG_LOG_IPRINT >= 2);
            auto ex = std::make_shared<Expect<SU2, double, double>>(me, (ubond_t)e.cfg.m,
                                                                    (ubond_t)e.cfg.m);
            ex->algo_type =
                ExpectationAlgorithmTypes::SymbolFree | ExpectationAlgorithmTypes::Compressed;
            ex->zero_dot_algo = true; // extract_root_single leaves the one-dot end-center form
            ex->iprint = DMRG_LOG_IPRINT;
            ex->cutoff = 1e-24;
            ex->solve(true, ket->center == 0);
            std::vector<std::shared_ptr<GTensor<double>>> npdm = ex->get_npdm();
            g.finish();
            remove_npdm_fragments(*me);
            me->remove_partition_files();

            if (npdm.size() != 1 || npdm[0] == nullptr || npdm[0]->size() != nel) {
                fprintf(out_stream, "ERROR: DMRG %d-body npdm shape mismatch (expected one"
                                    " n_act^%d = %zu element tensor)\n", N, 2 * N, nel);
                exit(EXIT_FAILURE);
            }
            if (acc == nullptr)
                acc = npdm[0];
            else {
                double *a = acc->data->data();
                const double *b = npdm[0]->data->data();
#pragma omp parallel for schedule(static)
                for (size_t i = 0; i < nel; i++)
                    a[i] += b[i];
            }
        }
        // No pmpo->deallocate(): the NPDM MPO's numeric legs are heap-owned site operators cached
        // in ghamil and several MPO entries alias the same one, so a tensor-wise deallocate
        // double-frees.
        remove_tag_files(ktag); // the per-root extracts and their single-MPS copies are transient
        remove_tag_files(kstag);
        if (bra_state != ket_state) {
            remove_tag_files(btag);
            remove_tag_files(bstag);
        }
        assert_stack_clean(tag); // the Expect sweep must leave the LIFO stacks as it found them
    }

    dump_lattice(*acc, N, ket_state, bra_state);
    return acc;
}

// ---- per-state 2-RDM (diagonal blocks) -----------------------------------------------------

// The n_s diagonal blocks of the state matrix, GAMMA convention in the delocalized basis, read
// off the per-state 2-RDMs the state-averaged read-out already formed. Overwrites the caller's
// n_s consecutive n_act^4 blocks (aldet's G_calc convention).
void block2_casci_wrap::G2_calc_diag(double *G) {
    dmrgci_engine &e = *impl_;
    require_solved(e, "G2_calc_diag");
    dmrg_log_guard log(false);
    host_threads_guard htg;
    const int n = e.n_act;
    const size_t blk = (size_t)n * n * n * n;
    ensure_2rdm(e);

    const std::vector<int> iperm = lattice_iperm(e);
    std::vector<double> perm_scr, rot_scr; // un-permute / back-transform targets (no aliasing)
    if (!iperm.empty()) perm_scr.resize(blk);
    if (e.localize_on) rot_scr.resize(blk);

    for (int s = 0; s < e.n_s; s++)
        g2full_finish_block(e, e.d2_states.data() + (size_t)s * blk, iperm, perm_scr.data(),
                            rot_scr.data(), G + (size_t)s * blk);
}

// ---- per-state 3-body moment (spin-summed 3-RDM) -------------------------------------------

// Raw lattice-ordered 3-body moment -> NOPT layout G3[p,q,r,i,j,k] = <a+_p a+_q a+_r a_k a_j a_i>
// = scale * raw[p,q,r,k,j,i]: the three annihilation axes reverse. Reversing axis positions and
// relabelling indices (lattice -> input, iperm; empty => lattice is the input order) commute, so
// both happen in this one gather. raw and out must not alias.
static void npdm3_gather(const double *raw, int n, const std::vector<int> &iperm, double scale,
                         double *out) {
    std::vector<int> ix(n);
    for (int a = 0; a < n; a++) ix[a] = iperm.empty() ? a : iperm[a];
    const size_t n2 = (size_t)n * n, n3 = n2 * n;
#pragma omp parallel for schedule(static) collapse(3)
    for (int p = 0; p < n; p++)
        for (int q = 0; q < n; q++)
            for (int r = 0; r < n; r++) {
                const size_t sa = ((size_t)ix[p] * n + ix[q]) * n + ix[r];
                const size_t da = ((size_t)p * n + q) * n + r;
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++) {
                        const double *s = raw + sa * n3 + (size_t)ix[j] * n + ix[i];
                        double *d = out + da * n3 + ((size_t)i * n + j) * n;
                        for (int k = 0; k < n; k++)
                            d[k] = scale * s[(size_t)ix[k] * n2];
                    }
            }
}

// One root's spin-summed 3-body moment in the native (delocalized) active basis, NOPT layout.
// A single NPDM Expect sweep on a transient extract of the root: the stored MultiMPS is only read.
// Overwrites the caller's n_act^6 buffer.
void block2_casci_wrap::G3_calc_diag(double *G3, int state) {
    dmrgci_engine &e = *impl_;
    require_solved(e, "G3_calc_diag");
    if (state < 0 || state >= e.n_s) {
        fprintf(out_stream, "ERROR: DMRG G3_calc_diag root %d out of range (n_s = %d)\n", state,
                e.n_s);
        exit(EXIT_FAILURE);
    }
    dmrg_log_guard log(false);
    host_threads_guard htg;
    const int n = e.n_act;
    const size_t blk6 = (size_t)n * n * n * n * n * n;

    std::shared_ptr<GTensor<double>> raw =
        npdm_lattice(e, 3, state, state, rdm_sweep_dir(), "NPDM3");
    npdm3_gather(raw->data->data(), n, lattice_iperm(e), 2.0 * std::sqrt(2.0), G3);
    raw = nullptr; // drop the raw n_act^6 tensor before the back-transform allocates
    if (e.localize_on) { // rotate back to the delocalized basis
        std::vector<double> rot(blk6);
        rotate3(G3, e.U_loc.data(), n, rot.data(), /*forward=*/false);
        std::copy(rot.begin(), rot.end(), G3);
    }
    if (std::getenv("NOPT_RDM3_STOP") != nullptr) { // the run ends once the 3-RDM is built
        fprintf(out_stream, "[RDM] stop after N=3 (NOPT_RDM3_STOP)\n");
        fflush(out_stream);
        exit(EXIT_SUCCESS);
    }
}

#if 0  // DIRECT lambda3 path (superseded by the explicit lattice-3RDM route; revive for nact >~ 30)

// ---- complementary six-operator overlap (3-RDM-free lambda3) -------------------------------

// Reorder the three active axes of a [np][n^3] tensor onto the Fiedler lattice (out[a,b,c] =
// in[perm[a],perm[b],perm[c]]); the external p axis is untouched. Empty perm => plain copy.
// Same gather direction as block2's FCIDUMP::reorder (site k carries orbital perm[k]).
static void reorder_active_axes(const double *in, int np, int n,
                                const std::vector<uint16_t> &perm, double *out) {
    const size_t na3 = (size_t)n * n * n;
    if (perm.empty()) {
        std::copy(in, in + (size_t)np * na3, out);
        return;
    }
    for (int p = 0; p < np; p++) {
        const double *ip = in + (size_t)p * na3;
        double *op = out + (size_t)p * na3;
        for (int a = 0; a < n; a++)
            for (int b = 0; b < n; b++)
                for (int c = 0; c < n; c++)
                    op[((size_t)a * n + b) * n + c] =
                        ip[((size_t)perm[a] * n + perm[b]) * n + perm[c]];
    }
}

// One external-leg (N-1)-electron operator Tp[w,u,v] -> spin-adapted MPO. The SU2 string
// "((C+D)0+D)1" couples creation w with the singlet-paired annihilation to spin 0, then the free
// annihilation to spin 1/2; the {n^2,1,n} stride triple is the storage->operator-leg transpose.
// Returns nullptr if the whole p-block falls below cutoff. Reproduces block2's get_mpo chain.
// The GeneralHamiltonian is built fresh for each operator: its on-site operator tables are
// populated on first use, and reusing one instance across successive builds corrupts every
// operator carrying coincident (on-site multi-operator) legs after the first.
static std::shared_ptr<MPO<SU2, double>>
build_caa_mpo(const std::vector<typename SU2::pg_t> &gorbsym, const double *Tp,
              int n, const std::vector<int> &orbsym0, double cutoff) {
    auto ghamil = std::make_shared<GeneralHamiltonian<SU2, double>>(SU2(0), n, gorbsym);
    const size_t na3 = (size_t)n * n * n;
    auto gfd = std::make_shared<GeneralFCIDUMP<double>>(ElemOpTypes::SU2);
    gfd->exprs.push_back("((C+D)0+D)1");
    gfd->add_sum_term(Tp, na3, {n, n, n}, {(size_t)n * n, (size_t)1, (size_t)n}, cutoff, 2.0,
                      orbsym0, {}, 0);
    std::shared_ptr<GeneralFCIDUMP<double>> afd = gfd->adjust_order();
    if (afd->exprs.empty())
        return nullptr; // whole p-block below cutoff
    auto gmpo = std::make_shared<GeneralMPO<SU2, double>>(ghamil, afd,
                                                          MPOAlgorithmTypes::FastBipartite,
                                                          0.0, -1, 0);
    gmpo->build();
    std::shared_ptr<MPO<SU2, double>> mpo = std::make_shared<SimplifiedMPO<SU2, double>>(
        gmpo, std::make_shared<Rule<SU2, double>>(), false, false);
    mpo = std::make_shared<IdentityAddedMPO<SU2, double>>(mpo);
    return mpo;
}

// Per root r: omega[r] = sum_p sum_channels <r| x+ y+ w z+ v u |r> contracted with Tbra/Tket.
// The Tket operator is applied to |r> and compressed to an intermediate MPS by a zero-noise
// perturbative-compression fit; the Tbra operator becomes the expectation MPO on that
// intermediate (Forte's argument swap, mathematically symmetric). Compression-approximate: the
// intermediate is truncated to bond dim h2caa_m (auto 2m); the compressed side is Tket (amplitude),
// the same physical side Forte compresses under the cert leg convention. The NOPT MPS is not
// singlet-embedded, so the recoupling channels are the quanta of the operator+target SU2 sum
// (one for a singlet reference, S+-1/2 otherwise). Active legs arrive in the frozen lattice
// basis; only the Fiedler reorder is applied here (the caller owns the U_loc rotation).
void block2_casci_wrap::h2caa_overlap(const double *Tbra, const double *Tket, int np,
                                      double *omega) {
    dmrgci_engine &e = *impl_;
    host_threads_guard htg;
    if (e.mps == nullptr || e.mps_info == nullptr) {
        fprintf(out_stream, "ERROR: DMRG h2caa_overlap needs a converged MPS (call solve first)\n");
        exit(EXIT_FAILURE);
    }
    const int n = e.n_act;
    const size_t na3 = (size_t)n * n * n;

    std::vector<double> Tb((size_t)np * na3), Tk((size_t)np * na3);
    reorder_active_axes(Tbra, np, n, e.reorder_perm, Tb.data());
    reorder_active_axes(Tket, np, n, e.reorder_perm, Tk.data());

    SU2 vacuum(0);
    std::vector<typename SU2::pg_t> gorbsym(n, 0); // C1 site irreps
    std::vector<int> orbsym0(n, 0);                // C1 add_sum_term filter (keeps every term)
    const double op_cutoff = 1e-12;      // operator-build magnitude threshold (Forte default)
    const int cps_m = e.cfg.h2caa_m > 0 ? e.cfg.h2caa_m : 2 * e.cfg.m; // compressed-intermediate bond dim (0 = auto: 2m)
    const int cps_sweeps = 2 * (e.cfg.sweeps > 0 ? e.cfg.sweeps : 8); // <= 2x the solve budget

    for (int r = 0; r < e.n_s; r++) {
        const std::string rtag = e.mps_info->tag + "-o" + std::to_string(r);
        const std::string rstag = rtag + "-s";
        std::shared_ptr<MPS<SU2, double>> psi = extract_root_single(e, r, rtag, rstag);

        double val = 0.0;
        for (int p = 0; p < np; p++) {
            std::shared_ptr<MPO<SU2, double>> kmpo =
                build_caa_mpo(gorbsym, Tk.data() + (size_t)p * na3, n, orbsym0, op_cutoff);
            if (kmpo == nullptr)
                continue; // whole p-block below cutoff
            std::shared_ptr<MPO<SU2, double>> bmpo =
                build_caa_mpo(gorbsym, Tb.data() + (size_t)p * na3, n, orbsym0, op_cutoff);
            if (bmpo == nullptr) {
                kmpo->deallocate();
                continue;
            }

            const SU2 bq = kmpo->op->q_label + psi->info->target; // intermediate target sector(s)
            for (int ch = 0; ch < bq.count(); ch++) {
                const std::string btag =
                    rtag + "-b" + std::to_string(p) + "-" + std::to_string(ch);
                auto binfo = std::make_shared<MPSInfo<SU2>>(n, vacuum, bq[ch], e.mpo->basis);
                binfo->tag = btag;
                // The spin-1/2 operator MPO injects its spin at the left boundary; the
                // intermediate's left vacuum must be the MPO's left vacuum, not the bare vacuum
                // (MovingEnvironment's boundary coupling rejects the default).
                binfo->set_bond_dimension_fci(kmpo->left_vacuum, vacuum);
                binfo->set_bond_dimension((ubond_t)cps_m);
                binfo->bond_dim = (ubond_t)cps_m;
                if (binfo->get_max_bond_dimension() == 0)
                    continue; // empty sector for this recoupling channel

                auto phi = std::make_shared<MPS<SU2, double>>(n, psi->center, psi->dot);
                phi->initialize(binfo);
                phi->random_canonicalize();
                phi->tensors[phi->center]->normalize();
                phi->save_mutable();
                phi->deallocate();
                binfo->save_mutable();
                binfo->deallocate_mutable();

                // phi ~= Tket-op |psi> by a zero-noise perturbative-compression fit. The reference is
                // a fresh copy: the fit mutates its environment side, never the extracted root.
                std::shared_ptr<MPS<SU2, double>> ref = psi->deep_copy(btag + "-ref");
                auto cme = std::make_shared<MovingEnvironment<SU2, double, double>>(kmpo, phi, ref,
                                                                                   "CAA-CPS");
                cme->delayed_contraction = OpNamesSet::normal_ops();
                cme->cached_contraction = true;
                cme->init_environments(true);
                std::vector<ubond_t> bdim{(ubond_t)cps_m},
                    kdim{ref->info->get_max_bond_dimension()};
                std::vector<double> noises{0.0};
                auto cps = std::make_shared<Linear<SU2, double, double>>(cme, bdim, kdim, noises);
                cps->iprint = 0;
                cps->noise_type = NoiseTypes::ReducedPerturbative;
                cps->eq_type = EquationTypes::PerturbativeCompression;
                cps->solve(cps_sweeps, phi->center == 0, 1e-8);
                // Compression is our approximation: surface (never abort) when the fit did not
                // demonstrably converge -- the last sweep's target change exceeds the 1e-8
                // tolerance, or it ran the whole sweep budget. delta_F bounds the overlap error.
                double dF = 0.0;
                if (cps->targets.size() >= 2)
                    dF = std::fabs(cps->targets.back().back() -
                                   cps->targets[cps->targets.size() - 2].back());
                const double max_dw =
                    cps->discarded_weights.empty()
                        ? 0.0
                        : *std::max_element(cps->discarded_weights.begin(),
                                            cps->discarded_weights.end());
                if (dF > 1e-8 || (int)cps->targets.size() >= cps_sweeps)
                    fprintf(out_stream,
                            "NOTE: h2caa compression under-converged (root %d p %d channel %d): "
                            "final |dF| = %.2e (tol 1e-8), max discarded weight = %.2e\n",
                            r, p, ch, dF, max_dw);
                if (phi->center != psi->center)
                    cps->solve(1, psi->center != 0); // align the center for the overlap sweep

                // omega_{p,ch} = <phi| Tbra-op |psi>.
                auto ome = std::make_shared<MovingEnvironment<SU2, double, double>>(bmpo, phi, psi,
                                                                                   "CAA-EXP");
                ome->delayed_contraction = OpNamesSet::normal_ops();
                ome->cached_contraction = true;
                ome->init_environments(false);
                auto ex = std::make_shared<Expect<SU2, double, double>>(ome, (ubond_t)cps_m,
                                                                        (ubond_t)cps_m);
                ex->iprint = 0;
                val += (double)ex->solve(false, psi->center != 0);

                cme->remove_partition_files();
                ome->remove_partition_files();
                remove_tag_files(btag);          // the intermediate MPS and its info are transient
                remove_tag_files(btag + "-ref"); // ... and the compression reference copy
            }
            kmpo->deallocate();
            bmpo->deallocate();
        }
        omega[r] = val; // per-root, unweighted (the caller applies the SA weights)

        remove_tag_files(rtag);
        remove_tag_files(rstag);
    }
    assert_stack_clean("h2caa overlap");
}

#endif
