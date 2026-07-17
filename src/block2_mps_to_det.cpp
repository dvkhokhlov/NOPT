// block2_mps_to_det — leading-determinant read-out for the block2 DMRG backend.
//
// Each state's MPS is rotated into the canonical active-orbital basis, its Slater determinants are
// enumerated in the M_S = na-nb sector (via the SZ transform), and the leading configurations printed.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>

#include "block2_dmrg_engine.h"   // dmrgci_engine, host_threads_guard, shared helpers
#include "mps_rotation.h"         // evolve_sa_multimps

using namespace block2;
using namespace nopt_block2;

// Print the lattice-site -> canonical-orbital map (with the per-orbital sign gauge) for the read-out.
static void report_dmrg_orbital_map(std::FILE *out, int n_act,
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

// The solve->canonical report basis change. R = <canon|solve> = U_canon*U_loc (localizing) / U_canon
// (canonicalizing a delocalized solve) / U_loc (localized only) / none (report in the solve basis).
// Peel the signed permutation Q off Rp = R^T: a greedy max-magnitude match leaves the proper residual
// Rt (drives the read-out rotation), and Q gives the per-site (reported orbital, sign) map. The sigma
// signs are a theory-fixed per-orbital gauge, reinstated on read-out (canonicalize_site_det).
struct report_basis_change {
    std::vector<double> residual;      // proper residual Rt in solve coords (empty => no rotation)
    std::vector<int> canon_of_lattice; // lattice site k -> reported orbital index
    std::vector<double> lattice_sign;  // per-lattice-site orbital sign gauge
    bool has_rotation = false;         // false => read-out sits in the solve basis
};

static report_basis_change build_report_basis_change(const dmrgci_engine &e, int n,
                                                     const std::vector<int> &reord,
                                                     const std::vector<int> &iperm) {
    report_basis_change out;
    out.canon_of_lattice.assign(n, 0);
    out.lattice_sign.assign(n, 1.0);

    std::vector<double> R_total;
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
    if (R == nullptr) { // no rotation: lattice site k holds solve orbital reord[k], signs +1
        for (int k = 0; k < n; k++) out.canon_of_lattice[k] = reord.empty() ? k : reord[k];
        return out;
    }
    out.has_rotation = true;

    // Rp = R^T; peeled in place to the residual Rt.
    std::vector<double> &Rp = out.residual;
    Rp.assign((size_t)n * n, 0.0);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) Rp[(size_t)i * n + j] = R[(size_t)j * n + i];
    const std::vector<double> Rp0 = Rp; // un-peeled copy: verifies the Rt/Q factorization below

    struct Cell { double mag; int i, j; };
    std::vector<Cell> cells;
    cells.reserve((size_t)n * n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            cells.push_back({std::fabs(Rp[(size_t)i * n + j]), i, j});
    std::sort(cells.begin(), cells.end(),
              [](const Cell &a, const Cell &b) { return a.mag > b.mag; });
    std::vector<int> rho(n, -1);        // rho[a] = canonical column matched to solve row a
    std::vector<double> sgn(n, 1.0);    // aligned-entry sign of each residual column
    std::vector<char> row_used(n, 0), col_used(n, 0);
    for (const Cell &c : cells) {
        if (row_used[c.i] || col_used[c.j]) continue;
        row_used[c.i] = col_used[c.j] = 1;
        rho[c.i] = c.j;
        sgn[c.i] = Rp[(size_t)c.i * n + c.j] < 0.0 ? -1.0 : 1.0;
    }
    std::vector<double> Rt((size_t)n * n); // Rt = Rp * Q: column b = sgn[b] * (column rho[b] of Rp)
    for (int a = 0; a < n; a++)
        for (int b = 0; b < n; b++)
            Rt[(size_t)a * n + b] = sgn[b] * Rp[(size_t)a * n + rho[b]];
    // Keep Rt proper (det<0 -> a -1 eigenvalue -> complex log): flip the least-aligned residual column
    // and fold the flip into sgn so it is reinstated with the rest of Q.
    if (MatrixFunctions::det(MatrixRef(Rt.data(), n, n)) < 0.0) {
        int bmin = 0;
        for (int b = 1; b < n; b++)
            if (std::fabs(Rt[(size_t)b * n + b]) < std::fabs(Rt[(size_t)bmin * n + bmin])) bmin = b;
        for (int a = 0; a < n; a++) Rt[(size_t)a * n + bmin] = -Rt[(size_t)a * n + bmin];
        sgn[bmin] = -sgn[bmin];
    }
    std::copy(Rt.begin(), Rt.end(), Rp.begin());

    // canonical orbital p sits at the lattice site of residual g[p]=rho^{-1}[p], with sign sgn[g[p]].
    std::vector<int> g(n, -1);
    for (int a = 0; a < n; a++) g[rho[a]] = a;
    for (int p = 0; p < n; p++) {
        const int resid = g[p];
        const int k = reord.empty() ? resid : iperm[resid];
        out.canon_of_lattice[k] = p;
        out.lattice_sign[k] = sgn[resid];
    }

    // Rt and Q are a factorization of Rp: Rp[a][j] = Rt[a][g[j]] * sgn[g[j]]. The read-out rotates the
    // MPS by Rt and reinstates Q by relabelling sites and signing single occupations, so a sign carried
    // by only one of the two halves is not a gauge choice -- it negates every determinant with orbital j
    // singly occupied. Q itself may be improper; only Rt must be proper, and it is made so above.
    double dev = 0.0;
    for (int a = 0; a < n; a++)
        for (int j = 0; j < n; j++) {
            const double rec = Rp[(size_t)a * n + g[j]] * sgn[g[j]];
            dev = std::max(dev, std::fabs(rec - Rp0[(size_t)a * n + j]));
        }
    if (dev > 1e-10) {
        fprintf(out_stream, "ERROR: DMRG read-out basis change does not factorize into a rotation and a "
                            "signed permutation (max deviation %.3e)\n", dev);
        exit(EXIT_FAILURE);
    }

    return out;
}

// Rotate a MultiMPS from the solve basis into the canonical reporting basis by the residual proper
// rotation U (n_act x n_act, [a*n+p]) the caller peeled off the signed permutation. Unlike the
// warm-start there is no norm-drift rejection -- the read-out accepts the truncation (captured weight).
// Returns false (leaving mps in the solve basis) only if the generator still comes back complex.
static bool rotate_multimps_to_canonical(dmrgci_engine &e,
                                         const std::shared_ptr<MultiMPS<SU2, double>> &mps,
                                         const double *U, int rot_m, int rot_steps) {
    auto res = apply_orbital_rotation_mps(mps, U, e.n_act, e.n_elec, e.twos, e.orbsym,
                                          e.reorder_perm, rot_m, rot_steps);
    if (res.complex_generator) {
        fprintf(out_stream, "  warning: read-out rotation generator not real -> reporting in solve basis\n");
        return false;
    }
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

struct det_image {
    det_occ occ;
    double phase = 1.0;
};

// Convert a determinant from the MPS lattice/site convention into canonical alpha-string/beta-string
// convention. block2's determinant Fock convention interleaves spin orbitals (a0,b0,a1,b1,...) for
// fermion phases; aldet prints separate alpha/beta strings. The inversion count applies that
// spin-order phase together with the Fiedler permutation and the canonical-orbital sign gauge.
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

static double expansion_weight(const det_expansion &dets) {
    double w = 0.0;
    for (const auto &kv : dets) w += kv.second * kv.second;
    return w;
}

struct su2_readout_expansion {
    det_expansion determinants;
    double weight = 0.0;
};

// Put the MPS in the one-site read-out form: left-canonical throughout with a single right-fused
// center at the last site. For a non-singlet target the SU2 -> SZ unfused transform rebuilds only
// LEFT-canonical site tensors, and a left-end center makes backward_right_fused assert on
// disagreeing SZ bond dims. The end is otherwise uncontrolled (the compression fit early-stops, so
// the sweep count varies with the state), so pinning it makes the read-out deterministic.
static void pin_right_end(const std::shared_ptr<MPS<SU2, double>> &mps,
                          const std::shared_ptr<CG<SU2>> &cg) {
    mps->load_mutable();
    mps->info->load_mutable();
    const int n = mps->n_sites;
    std::string &cf = mps->canonical_form;
    if (mps->center == 0 && cf[0] == 'C' && n > 1 && cf[1] == 'R')
        cf[0] = 'K';                                     // one-site left-fused center at 0
    else if (n > 1 && cf[n - 1] == 'C' && cf[n - 2] == 'L') {
        cf[n - 1] = 'S';                                 // one-site right-fused center at n-1
        mps->center = n - 1;
    } else if (mps->center == n - 2 && n > 1 && cf[n - 2] == 'L')
        mps->center = n - 1;
    while (mps->center < n - 1) mps->move_right(cg, nullptr);
    mps->save_data();
}

// Put the single MPS into a one-site canonical center at site 0 (the form to_singlet_embedding_wfn
// needs), then singlet-embed it. The compressed read-out MPS lands in one of block2's post-sweep
// canonical forms; the branches below relabel that boundary center as a one-site fused center exactly
// as block2main's trans_mps_to_singlet_embedding does, then sweep any interior center to 0.
static void singlet_embed_su2_mps(const std::shared_ptr<MPS<SU2, double>> &mps,
                                  const std::shared_ptr<CG<SU2>> &cg) {
    mps->load_mutable();
    mps->info->load_mutable();
    const int n = mps->n_sites;
    std::string &cf = mps->canonical_form;
    if (mps->center == 0 && cf[0] == 'C' && n > 1 && cf[1] == 'R')
        cf[0] = 'K';                                     // one-site left-fused center at 0
    else if (n > 1 && cf[n - 1] == 'C' && cf[n - 2] == 'L') {
        cf[n - 1] = 'S';                                 // one-site right-fused center at n-1
        mps->center = n - 1;
    } else if (mps->center == n - 2 && n > 1 && cf[n - 2] == 'L')
        mps->center = n - 1;
    while (mps->center > 0) mps->move_left(cg, nullptr); // sweep the center down to site 0
    mps->to_singlet_embedding_wfn(cg);
}

// Enumerate the retained MPS's leading Slater determinants (SZ DeterminantTRIE) in the requested
// canonical orbital order/sign gauge, routed through the SZ M_S = twosz/2 sector: block2 cannot
// enumerate from a spin-adapted MPS of a non-singlet state. twosz == 0 takes the plain SU2->SZ
// transform. twosz != 0 needs singlet-embed -> M_S = 0 -> resolve_singlet_embedding(twosz), whose
// resolved sector carries norm^2 = 1/(2S+1), so coefficients are rescaled by sqrt(2S+1).
static su2_readout_expansion canonical_readout_from_su2_mps(
    const std::shared_ptr<MPS<SU2, double>> &su2_mps,
    int n_sites, int n_elec, int twos, int twosz, double cutoff,
    const std::vector<int> &canon_of_lattice,
    const std::vector<double> &lattice_sign,
    const std::shared_ptr<CG<SU2>> &cg, const std::string &sztag) {
    su2_readout_expansion out;

    std::shared_ptr<UnfusedMPS<SZ, double>> sz_umps;
    double scale = 1.0;
    if (twosz == 0) {
        pin_right_end(su2_mps, cg);  // non-singlet targets only transform from the left-canonical form
        auto su2_umps = std::make_shared<UnfusedMPS<SU2, double>>(su2_mps);
        SZ targetz(n_elec, 0, 0);  // M_S = 0: the native (aldet) projection sector
        sz_umps = TransUnfusedMPS<SU2, SZ, double>::forward(su2_umps, sztag, cg, targetz);
    } else {
        singlet_embed_su2_mps(su2_mps, cg);  // target twos -> 0, n -> n_elec + twos (the M_S = 0 embedding)
        auto su2_umps = std::make_shared<UnfusedMPS<SU2, double>>(su2_mps);
        SZ targetz(su2_mps->info->target.n(), 0, 0);
        sz_umps = TransUnfusedMPS<SU2, SZ, double>::forward(su2_umps, sztag, cg, targetz);
        sz_umps->resolve_singlet_embedding(twosz);  // project onto M_S = na-nb, restore the physical n
        scale = std::sqrt((double)(twos + 1));       // undo the 1/sqrt(2S+1) singlet-embedding norm
    }
    // Round-trip through a finalized MPS: the raw unfused transform has no clean single-target final
    // state, so evaluate would trip its final-target assert; finalize rebuilds a proper one.
    sz_umps = std::make_shared<UnfusedMPS<SZ, double>>(sz_umps->finalize());
    auto dtrie = std::make_shared<DeterminantTRIE<SZ, double>>(n_sites, true);
    // The singlet embedding leaves this sector at norm^2 = 1/(2S+1), and the TRIE prunes on the
    // coefficient magnitude in the MPS's own normalization -- so it must be given the cutoff in that
    // same normalization. Passing the physical one would make it act as cutoff*sqrt(2S+1): stricter
    // for every non-singlet, dropping determinants that were asked for and biasing the weight low.
    dtrie->evaluate(sz_umps, cutoff / scale);
    for (int t = 0; t < (int)dtrie->size(); t++) {
        const double c = (double)dtrie->vals[t] * scale;
        if (std::fabs(c) < 1e-14) continue;
        // canonicalize_site_det already applies the interleaved->alpha/beta spin-order parity together
        // with the orbital-order parity, so do NOT call dtrie->convert_phase: it would apply the same
        // two parities a second time and flip signs.
        det_image img = canonicalize_site_det((*dtrie)[t], canon_of_lattice, lattice_sign);
        add_det_coeff(out.determinants, img.occ, c * img.phase);
    }
    remove_tag_files(sztag);

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

    if (dets.empty()) {  // read-out captured zero weight (SZ transform empty in this M_S sector)
        std::fprintf(out, "  (no determinants captured -- the canonical-basis open-shell rotation left "
                          "nothing to enumerate here; use localize=pm for this state's read-out)\n");
        return;
    }

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

// `reported` marks the roots that actually produced an expansion; a root whose read-out was declined
// has no weight to quote and must not be read as one that captured nothing.
static void report_determinant_weights(std::FILE *out, const std::vector<double> &weights,
                                       const std::vector<char> &reported) {
    std::fprintf(out, "\nDMRG determinant extraction weights:\n");
    std::fprintf(out, "____________________________\n");
    std::fprintf(out, " State| weight       |\n");
    std::fprintf(out, "______|______________|\n");
    for (int st = 0; st < (int)weights.size(); st++)
        if (reported[st])
            std::fprintf(out, "%5d | %.6f     |\n", st, weights[st]);
        else
            std::fprintf(out, "%5d | n/a          |\n", st);
    std::fprintf(out, "______|______________|\n");

    // The printed configurations are a qualitative picture, not a quantitative norm; say so when the
    // extraction keeps too little of the vector for even that to be trusted.
    const double note_below = 0.9;
    double lo = 1.0;
    bool any = false;
    for (int st = 0; st < (int)weights.size(); st++)
        if (reported[st]) { lo = std::min(lo, weights[st]); any = true; }
    if (any && lo < note_below)
        std::fprintf(out, "  note: only %.3f of the CI weight is captured -- the leading-configuration\n"
                          "        picture is incomplete (lower $DMRG extract_cutoff, and raise or drop\n"
                          "        extract_m if the extraction MPS is being compressed)\n", lo);
}

// Report each state's leading determinant expansion in the canonical active-orbital basis. The solve
// runs in the localized + Fiedler-ordered basis, so each MPS is rotated by the proper (continuous)
// part of the solve->canonical change, expanded into determinants, and the remaining discrete map
// (lattice orbital -> canonical orbital, order and sign) applied to those determinants.
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

    report_basis_change rbc = build_report_basis_change(e, n, reord, iperm);
    const char *basis_label = e.have_canon ? "canonical orbital" : "active orbital";
    report_dmrg_orbital_map(out_stream, n, rbc.canon_of_lattice, rbc.lattice_sign, basis_label);

    std::vector<double> det_weights((size_t)e.n_s, 0.0);
    std::vector<char> det_reported((size_t)e.n_s, 1);

    // CG factors for the SU2->SZ read-out transform (live off the Hamiltonian; a fresh one is sized in
    // its constructor should the Hamiltonian ever be released before the print).
    std::shared_ptr<CG<SU2>> cg = e.mpo->tf->opf->cg;
    if (cg == nullptr) cg = std::make_shared<CG<SU2>>();

    for (int st = 0; st < e.n_s; st++) {
        // Work on a copy of this root; e.mps stays in the localized basis for the property read-out.
        const std::string xtag = e.mps_info->tag + "-cf" + std::to_string(st);
        std::shared_ptr<MultiMPS<SU2, double>> imps = e.mps->extract(st, xtag);
        // Rotate by the proper residual. A complex residual (guarded inside) leaves the MPS in the solve
        // basis, while the orbitals this run writes carry U_canon: printing the solve-basis determinants
        // would hand out coefficients and orbitals that refer to different bases. Report nothing for this
        // root rather than something incompatible with the orbital file shipped alongside it.
        if (rbc.has_rotation && !rotate_multimps_to_canonical(e, imps, rbc.residual.data(), rot_m, rot_steps)) {
            fprintf(out_stream, "State %d  E  = % 18.10f S^2 = %.2f:\n", st, e.E_states[st], S2);
            fprintf(out_stream, "  (no determinant read-out: the solve -> canonical rotation generator is "
                                "not real, so this root cannot be expressed in the orbital basis the CAS\n"
                                "   orbitals are written in)\n");
            det_reported[st] = 0;
            remove_tag_files(xtag);
            continue;
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

        const std::string sztag = e.mps_info->tag + "-sz" + std::to_string(st);
        su2_readout_expansion canonical =
            canonical_readout_from_su2_mps(rmps, n, e.n_elec, e.twos, e.twosz, cutoff,
                                           rbc.canon_of_lattice, rbc.lattice_sign, cg, sztag);

        det_weights[st] = canonical.weight;
        report_leading_determinants(out_stream, st, e.E_states[st], S2, canonical.determinants,
                                    print_number);

        remove_tag_files(xtag);
        remove_tag_files(stag);
        if (!ctag.empty()) remove_tag_files(ctag);
    }
    report_determinant_weights(out_stream, det_weights, det_reported);
    assert_stack_clean("leading-config read");
}
void block2_casci_wrap::write_civec(int, char *) { /* MPS write-out deferred */ }
