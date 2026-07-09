// block2_mps_to_det — leading-determinant read-out for the block2 DMRG backend.
//
// Each state's MPS is rotated into the canonical active-orbital basis, its determinants are
// extracted (SU2 step-vector CSFs expanded locally), and the leading configurations printed.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
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
