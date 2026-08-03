// block2_casci_wrap density read-outs over the external block2 DMRG library: the per-state
// spin-summed 2-RDM (G2_calc_diag) and 3-body moment (G3_calc_diag). All port Forte's block2
// primitives onto NOPT's own scaffolds; split from block2_dmrg.cpp, same author idiom.

#include "block2_dmrg_engine.h"   // block2 headers + dmrgci_engine + shared helpers

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <string>
#include <system_error>
#include <vector>

#include "common_vars.h"      // out_stream
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

#if 0  // full transition 2-RDM: no consumer, the driver reads G2_calc_diag. Revive for first-order
       // properties -- the 2-body Mbar needs bra != ket densities between the dressed roots.
// Full n_s x n_s spin-summed 2-RDM, once per solve: one Expect sweep per root pair i<=j of the
// single SA MultiMPS yields the block2 2-RDM <i| a+ a+ a a |j>; the (j,i) block follows by the
// operator adjoint (full index reversal on the raw layout, before the un-permute/back-transform
// chain that all four legs share and therefore commute with). Mirrors ensure_2rdm/ensure_dm_full.
static void ensure_g2full(dmrgci_engine &e) {
    if (e.g2full_valid)
        return;
    host_threads_guard htg;
    const int n = e.n_act;
    const size_t blk = (size_t)n * n * n * n;
    e.dg2full.assign(blk * e.n_s * e.n_s, 0.0);

    std::vector<int> iperm; // inverse of reorder_perm (Fiedler lattice -> input)
    if (!e.reorder_perm.empty()) {
        iperm.resize(n);
        for (int i = 0; i < n; i++) iperm[e.reorder_perm[i]] = i;
    }
    std::vector<double> perm_scr, rot_scr; // un-permute / back-transform targets (no aliasing)
    if (!iperm.empty()) perm_scr.resize(blk);
    if (e.localize_on) rot_scr.resize(blk);

    // State-pair-independent MPO. NoTransposeRule: the transpose-symmetry simplification is
    // invalid when bra != ket (the off-diagonal transition blocks).
    std::shared_ptr<MPO<SU2, double>> p2mpo = std::make_shared<PDM2MPOQC<SU2, double>>(e.hamil);
    p2mpo = std::make_shared<SimplifiedMPO<SU2, double>>(
        p2mpo,
        std::make_shared<NoTransposeRule<SU2, double>>(std::make_shared<RuleQC<SU2, double>>()),
        true, true, OpNamesSet({OpNames::R, OpNames::RD}));

    std::vector<double> raw(blk), rev; // raw block2 tensor + adjoint scratch (i<j only)
    if (e.n_s > 1) rev.resize(blk);

    for (int i = 0; i < e.n_s; i++)
        for (int j = i; j < e.n_s; j++) {
            // Fresh extract of both roots so bra/ket share a canonical center (the single-MPS
            // form is all-or-nothing: the effective Hamiltonian asserts bra and ket agree).
            const std::string itag = e.mps_info->tag + "-g" + std::to_string(i);
            const std::string istag = itag + "-s";
            std::shared_ptr<MPS<SU2, double>> imps = extract_root_single(e, i, itag, istag);
            std::string jtag, jstag;
            std::shared_ptr<MPS<SU2, double>> jmps = imps;
            if (j != i) {
                jtag = e.mps_info->tag + "-g" + std::to_string(j);
                jstag = jtag + "-s";
                jmps = extract_root_single(e, j, jtag, jstag);
            }

            auto p2me = std::make_shared<MovingEnvironment<SU2, double, double>>(p2mpo, imps, jmps,
                                                                                 "T2PDM");
            p2me->init_environments(false);
            auto ex2 = std::make_shared<Expect<SU2, double, double>>(p2me, (ubond_t)e.cfg.m,
                                                                     (ubond_t)e.cfg.m);
            ex2->iprint = 0; // silence the per-site Expect sweep log
            ex2->solve(true, jmps->center == 0);
            std::shared_ptr<GTensor<double>> d2 = ex2->get_2pdm_spatial(); // {n,n,n,n}
            std::copy(d2->data->data(), d2->data->data() + blk, raw.data());

            g2full_finish_block(e, raw.data(), iperm, perm_scr.data(), rot_scr.data(),
                                e.dg2full.data() + (size_t)(i * e.n_s + j) * blk);
            if (j != i) {
                // (j,i) = adjoint of (i,j): d2_ji[p,q,r,s] = d2_ij[s,r,q,p] on the raw layout
                for (int p = 0; p < n; p++)
                    for (int q = 0; q < n; q++)
                        for (int r = 0; r < n; r++)
                            for (int s = 0; s < n; s++)
                                rev[(((size_t)p * n + q) * n + r) * n + s] =
                                    raw[(((size_t)s * n + r) * n + q) * n + p];
                g2full_finish_block(e, rev.data(), iperm, perm_scr.data(), rot_scr.data(),
                                    e.dg2full.data() + (size_t)(j * e.n_s + i) * blk);
            }

            p2me->remove_partition_files();
            remove_tag_files(itag); // the per-root extracts and their single-MPS copies are transient
            remove_tag_files(istag);
            if (j != i) {
                remove_tag_files(jtag);
                remove_tag_files(jstag);
            }
        }
    p2mpo->deallocate();
    e.g2full_valid = true;
    assert_stack_clean("transition 2-RDM read"); // the Expect sweeps must leave the stacks as they found them
}

void block2_casci_wrap::G_calc_full(double *G) {
    dmrgci_engine &e = *impl_;
    ensure_g2full(e);
    // GAMMA convention, delocalized basis; caller zeroes, we accumulate (matches every RDM call).
    const size_t nel = e.dg2full.size();
    for (size_t k = 0; k < nel; k++)
        G[k] += e.dg2full[k];
}
#endif

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

// ---- per-state 2-RDM (diagonal blocks) -----------------------------------------------------

// The n_s diagonal blocks of the state matrix, one Expect sweep per root, GAMMA convention in
// the delocalized basis. Overwrites the caller's n_s consecutive n_act^4 blocks (aldet's
// G_calc convention).
void block2_casci_wrap::G2_calc_diag(double *G) {
    dmrgci_engine &e = *impl_;
    const int n = e.n_act;
    const size_t blk = (size_t)n * n * n * n;
    if (e.g2full_valid) { // the full state matrix is already cached: read its diagonal
        for (int s = 0; s < e.n_s; s++)
            std::copy(e.dg2full.begin() + (size_t)(s * e.n_s + s) * blk,
                      e.dg2full.begin() + (size_t)(s * e.n_s + s + 1) * blk, G + (size_t)s * blk);
        return;
    }
    require_solved(e, "G2_calc_diag");
    host_threads_guard htg;

    const std::vector<int> iperm = lattice_iperm(e);
    std::vector<double> perm_scr, rot_scr; // un-permute / back-transform targets (no aliasing)
    if (!iperm.empty()) perm_scr.resize(blk);
    if (e.localize_on) rot_scr.resize(blk);

    // bra == ket throughout, so RuleQC's transpose-symmetry simplification is valid (the full
    // state matrix has to fall back on NoTransposeRule for its off-diagonal blocks).
    std::shared_ptr<MPO<SU2, double>> p2mpo = std::make_shared<PDM2MPOQC<SU2, double>>(e.hamil);
    p2mpo = std::make_shared<SimplifiedMPO<SU2, double>>(
        p2mpo, std::make_shared<RuleQC<SU2, double>>(), true, true,
        OpNamesSet({OpNames::R, OpNames::RD}));

    for (int s = 0; s < e.n_s; s++) {
        const std::string xtag = e.mps_info->tag + "-d" + std::to_string(s);
        const std::string stag = xtag + "-s";
        std::shared_ptr<MPS<SU2, double>> smps = extract_root_single(e, s, xtag, stag);

        auto p2me = std::make_shared<MovingEnvironment<SU2, double, double>>(p2mpo, smps, smps,
                                                                            "D2PDM");
        p2me->init_environments(false);
        auto ex2 = std::make_shared<Expect<SU2, double, double>>(p2me, (ubond_t)e.cfg.m,
                                                                 (ubond_t)e.cfg.m);
        ex2->iprint = 0; // silence the per-site Expect sweep log
        ex2->solve(true, smps->center == 0);
        std::shared_ptr<GTensor<double>> d2 = ex2->get_2pdm_spatial(); // {n,n,n,n}
        g2full_finish_block(e, d2->data->data(), iperm, perm_scr.data(), rot_scr.data(),
                            G + (size_t)s * blk);

        p2me->remove_partition_files();
        remove_tag_files(xtag); // the per-root extract and its single-MPS copy are transient
        remove_tag_files(stag);
    }
    p2mpo->deallocate();
    assert_stack_clean("per-state 2-RDM read");
}

// ---- per-state 3-body moment (spin-summed 3-RDM) -------------------------------------------

// SU2 recoupling string of the spin-summed 3-body density. block2 returns the singlet-coupled
// raw[x0,x1,x2,y0,y1,y2] = sum_spins <a+_x0 a+_x1 a+_x2 a_y0 a_y1 a_y2> scaled by 2^{-n_cds/4},
// n_cds = 6 operators.
static const char *const npdm3_expr = "((C+((C+(C+D)0)1+D)0)1+D)0";

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
    host_threads_guard htg;
    const int n = e.n_act;
    const size_t blk6 = (size_t)n * n * n * n * n * n;

    // The GeneralHamiltonian is built fresh for each MPO: its on-site operator tables are populated
    // on first use, and reusing one instance corrupts every operator carrying coincident legs.
    SU2 vacuum(0);
    std::vector<typename SU2::pg_t> gorbsym(n, 0); // C1 site irreps
    auto ghamil = std::make_shared<GeneralHamiltonian<SU2, double>>(vacuum, n, gorbsym);

    auto perm = std::make_shared<SpinPermScheme>(
        SpinPermScheme::initialize_su2(6, npdm3_expr, /*is_npdm=*/true));
    auto ppmpo = std::make_shared<GeneralNPDMMPO<SU2, double>>(
        ghamil, std::make_shared<NPDMScheme>(perm), /*symbol_free=*/true, 0.0, 0);
    ppmpo->delta_quantum = SU2(0, SpinPermRecoupling::get_target_twos(npdm3_expr), 0);
    ppmpo->build();
    std::shared_ptr<MPO<SU2, double>> pmpo = std::make_shared<SimplifiedMPO<SU2, double>>(
        ppmpo, std::make_shared<Rule<SU2, double>>(), false, false);

    const std::string xtag = e.mps_info->tag + "-p3" + std::to_string(state);
    const std::string stag = xtag + "-s";

    std::vector<std::shared_ptr<GTensor<double>>> npdm; // {n,n,n,n,n,n}, lattice order
    {
        std::shared_ptr<MPS<SU2, double>> psi = extract_root_single(e, state, xtag, stag);
        auto pme = std::make_shared<MovingEnvironment<SU2, double, double>>(pmpo, psi, psi,
                                                                           "NPDM3");
        pme->cached_contraction = false; // conflicts with the fused zero-dot contraction
        pme->fused_contraction_rotation = true;
        pme->init_environments(false);
        auto ex = std::make_shared<Expect<SU2, double, double>>(pme, (ubond_t)e.cfg.m,
                                                                (ubond_t)e.cfg.m);
        ex->algo_type =
            ExpectationAlgorithmTypes::SymbolFree | ExpectationAlgorithmTypes::Compressed;
        ex->zero_dot_algo = true; // extract_root_single leaves the one-dot end-center form
        ex->iprint = 0;
        ex->cutoff = 1e-24;
        ex->solve(true, psi->center == 0);
        npdm = ex->get_npdm();
        remove_npdm_fragments(*pme);
        pme->remove_partition_files();
    }
    // No pmpo->deallocate(): the NPDM MPO's numeric legs are heap-owned site operators cached in
    // ghamil and several MPO entries alias the same one, so a tensor-wise deallocate double-frees.
    // block2's own npdm driver drops the MPO the same way; assert_stack_clean is the leak check.
    remove_tag_files(xtag); // the per-root extract and its single-MPS copy are transient
    remove_tag_files(stag);
    assert_stack_clean("3-body moment read");

    if (npdm.size() != 1 || npdm[0] == nullptr || npdm[0]->size() != blk6) {
        fprintf(out_stream, "ERROR: DMRG 3-body npdm shape mismatch (expected one n_act^6 = %zu"
                            " element tensor)\n", blk6);
        exit(EXIT_FAILURE);
    }
    npdm3_gather(npdm[0]->data->data(), n, lattice_iperm(e), 2.0 * std::sqrt(2.0), G3);
    npdm[0] = nullptr; // drop the raw n_act^6 tensor before the back-transform allocates
    if (e.localize_on) { // rotate back to the delocalized basis
        std::vector<double> rot(blk6);
        rotate3(G3, e.U_loc.data(), n, rot.data(), /*forward=*/false);
        std::copy(rot.begin(), rot.end(), G3);
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
