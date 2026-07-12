#include "localizer.h"

#include "blas_link.h"    // cblas_dgemm
#include "molecule.h"     // molecule
#include "common_vars.h"  // out_stream

#include <cmath>
#include <cstdio>
#include <cstring>
#include <algorithm>

namespace {

inline FILE* loc_os() { return out_stream ? out_stream : stderr; }

inline void set_identity(double* U, int n) {
    std::memset(U, 0, (size_t)n * n * sizeof(double));
    for (int i = 0; i < n; ++i) U[i * n + i] = 1.0;
}

// PM functional L = sum_A sum_i (Q^A_ii)^2.
double pm_functional(const double* Qpops, int n_atoms, int n_blk) {
    double L = 0.0;
    for (int a = 0; a < n_atoms; ++a) {
        const double* Q = Qpops + (size_t)a * n_blk * n_blk;
        for (int i = 0; i < n_blk; ++i) {
            double q = Q[i * n_blk + i];
            L += q * q;
        }
    }
    return L;
}

} // namespace

loc_result pm_jacobi_sweep(double* Qpops, int n_atoms, int n_blk,
                           double* U, const loc_options& opt) {
    const double angle_eps = 1e-14;   // skip rotations below this
    const double L0 = pm_functional(Qpops, n_atoms, n_blk);

    loc_result res;
    res.L = L0;

    double L = L0;
    double last_grad = 0.0;
    for (int sweep = 0; sweep < opt.max_sweeps; ++sweep) {
        double max_grad = 0.0;
        for (int s = 0; s < n_blk; ++s)
        for (int t = s + 1; t < n_blk; ++t) {
            // pair coefficients A_st, B_st accumulated over atoms
            double Acoef = 0.0, Bcoef = 0.0;
            for (int a = 0; a < n_atoms; ++a) {
                const double* Q = Qpops + (size_t)a * n_blk * n_blk;
                double d    = Q[s * n_blk + t];
                double diff = Q[s * n_blk + s] - Q[t * n_blk + t];
                Acoef += d * d - 0.25 * diff * diff;
                Bcoef += d * diff;
            }
            max_grad = std::max(max_grad, std::fabs(Bcoef));

            // Optimal angle. The atan2 form self-escapes the A>0,B~0 symmetric-bond saddle by
            // returning alpha = pi/4 (a plain 0.5*atan(B/A) would stall there). See localizer.h.
            double alpha = 0.25 * std::atan2(Bcoef, -Acoef);
            double c  = std::cos(alpha);
            double sn = std::sin(alpha);
            if (std::fabs(sn) < angle_eps) continue;

            // rotate every Q^A: Q' = R^T Q R, R = [[c,-sn],[sn,c]] in the (s,t) subspace.
            for (int a = 0; a < n_atoms; ++a) {
                double* Q = Qpops + (size_t)a * n_blk * n_blk;
                for (int r = 0; r < n_blk; ++r) {
                    if (r == s || r == t) continue;
                    double qsr = Q[s * n_blk + r];
                    double qtr = Q[t * n_blk + r];
                    double nsr =  c * qsr + sn * qtr;
                    double ntr = -sn * qsr + c * qtr;
                    Q[s * n_blk + r] = nsr; Q[r * n_blk + s] = nsr;
                    Q[t * n_blk + r] = ntr; Q[r * n_blk + t] = ntr;
                }
                // 2x2 block (read a_/b_/d_ before writing)
                double a_ = Q[s * n_blk + s];
                double b_ = Q[t * n_blk + t];
                double d_ = Q[s * n_blk + t];
                Q[s * n_blk + s] = c * c * a_ + 2.0 * c * sn * d_ + sn * sn * b_;
                Q[t * n_blk + t] = sn * sn * a_ - 2.0 * c * sn * d_ + c * c * b_;
                double nst = (c * c - sn * sn) * d_ + c * sn * (b_ - a_);
                Q[s * n_blk + t] = nst; Q[t * n_blk + s] = nst;
            }
            // accumulate the rotation into U:  U <- U * R  (columns s,t)
            for (int row = 0; row < n_blk; ++row) {
                double us = U[row * n_blk + s];
                double ut = U[row * n_blk + t];
                U[row * n_blk + s] =  c * us + sn * ut;
                U[row * n_blk + t] = -sn * us + c * ut;
            }
        }

        double L_new = pm_functional(Qpops, n_atoms, n_blk);
        double dL = L_new - L;          // monotone non-decreasing (each 2x2 rotation maximizes its pair)
        L = L_new;
        last_grad = max_grad;
        if (opt.verbose)
            fprintf(loc_os(), "localizer: sweep %4d  L=%.12f  max|grad|=%.3e\n",
                    sweep, L_new, max_grad);

        if (max_grad < opt.grad_tol || dL < opt.ftol) {
            res.converged  = true;
            res.n_sweeps   = sweep + 1;
            res.final_grad = max_grad;
            res.L          = L;
            return res;
        }
    }

    // not converged within max_sweeps: fall back to delocalized (U=I), observable via converged=false.
    set_identity(U, n_blk);
    if (opt.verbose)
        fprintf(loc_os(), "localizer: NOT converged after %d sweeps (max|grad|=%.3e) -> U=I fallback\n",
                opt.max_sweeps, last_grad);
    res.converged  = false;
    res.n_sweeps   = opt.max_sweeps;
    res.final_grad = last_grad;
    res.L          = L0;               // report the delocalized functional, matching the U=I we return
    return res;
}

loc_result orbital_localizer::localize(const double* C, int n_ao, int n_sub,
                                       const int* orb_irrep, double* U,
                                       const loc_options& opt) {
    if (opt.by_irrep) {
        // provisioned but not implemented this step — error loudly, no silent fallback.
        fprintf(loc_os(),
                "localizer: irrep-separated localization (by_irrep=true) not implemented yet\n");
        set_identity(U, n_sub);
        loc_result r;          // converged=false, n_sweeps=0
        return r;
    }

    // seed U (identity, or the optional warm-start U0)
    if (opt.U0) std::memcpy(U, opt.U0, (size_t)n_sub * n_sub * sizeof(double));
    else        set_identity(U, n_sub);

    (void)orb_irrep;           // reserved until by_irrep is implemented
    return localize_block(C, n_ao, n_sub, U, opt);
}

pm_localizer::pm_localizer(molecule& mol) {
    n_ao_    = mol.n_ao;
    n_atoms_ = mol.n_atoms;
    if (mol.S_AO == nullptr) {
        fprintf(loc_os(),
                "pm_localizer: molecule::S_AO is null (call calc_S_AO() first); localizer disabled\n");
        return;                // S_AO_ stays null; localize_block guards on it
    }

    // AO -> atom from the shell centers (AOs of a shell are contiguous and sit on its center).
    ao_atom_.assign((size_t)n_ao_, -1);
    int ao_num = 0;
    for (size_t sh = 0; sh < mol.s.size(); ++sh) {
        int sh_s = mol.s[sh].contr[0].size();
        for (int i = 0; i < sh_s && ao_num + i < n_ao_; ++i)
            ao_atom_[ao_num + i] = mol.shell_center[sh];
        ao_num += sh_s;
    }
    if (ao_num != n_ao_) {
        fprintf(loc_os(),
                "pm_localizer: shells span %d AOs, expected %d; localizer disabled\n", ao_num, n_ao_);
        ao_atom_.clear();
        return;
    }
    S_AO_ = mol.S_AO;
}

loc_result pm_localizer::localize_block(const double* C_blk, int n_ao, int n_g,
                                        double* Ug, const loc_options& opt) {
    if (S_AO_ == nullptr || n_ao != n_ao_) {
        if (n_ao != n_ao_)
            fprintf(loc_os(), "pm_localizer: n_ao mismatch (%d vs cached %d) -> U=I fallback\n",
                    n_ao, n_ao_);
        set_identity(Ug, n_g);
        loc_result r;          // converged=false
        return r;
    }

    // Seed-consistent orbitals: build Q^A from C*Ug so a non-identity Ug (warm start) is honored.
    // For the default Ug=I this is just a copy of C_blk.
    std::vector<double> Ceff((size_t)n_ao_ * n_g, 0.0);
    for (int ao = 0; ao < n_ao_; ++ao)
        for (int p = 0; p < n_g; ++p) {
            double acc = 0.0;
            for (int o = 0; o < n_g; ++o)
                acc += C_blk[(size_t)ao * n_g + o] * Ug[o * n_g + p];
            Ceff[(size_t)ao * n_g + p] = acc;
        }

    // SC = S_AO * Ceff (n_ao x n_g), then each atom's Q^A (n_g x n_g, symmetric) straight from its
    // AO rows: Q^A_pq = 1/2 sum_{mu in A} ( Ceff[mu,p] SC[mu,q] + SC[mu,p] Ceff[mu,q] ).
    std::vector<double> SC((size_t)n_ao_ * n_g, 0.0);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n_ao_, n_g, n_ao_, 1.0,
                S_AO_, n_ao_,
                Ceff.data(), n_g, 0.0,
                SC.data(), n_g);

    std::vector<double> Qpops((size_t)n_atoms_ * n_g * n_g, 0.0);
    for (int ao = 0; ao < n_ao_; ++ao) {
        double* Q = Qpops.data() + (size_t)ao_atom_[ao] * n_g * n_g;
        const double* ce = Ceff.data() + (size_t)ao * n_g;
        const double* sc = SC.data()   + (size_t)ao * n_g;
        for (int p = 0; p < n_g; ++p)
            for (int q = 0; q < n_g; ++q)
                Q[p * n_g + q] += 0.5 * (ce[p] * sc[q] + sc[p] * ce[q]);
    }

    if (opt.verbose) {
        // self-check: sum_A Q^A should be the identity for S-orthonormal C.
        double diagerr = 0.0, offmax = 0.0;
        for (int i = 0; i < n_g; ++i)
            for (int j = 0; j < n_g; ++j) {
                double s = 0.0;
                for (int a = 0; a < n_atoms_; ++a)
                    s += Qpops[(size_t)a * n_g * n_g + i * n_g + j];
                if (i == j) diagerr = std::max(diagerr, std::fabs(s - 1.0));
                else        offmax  = std::max(offmax, std::fabs(s));
            }
        fprintf(loc_os(), "pm_localizer: sum_A Q^A - I : max diag err=%.3e, max offdiag=%.3e\n",
                diagerr, offmax);
    }

    return pm_jacobi_sweep(Qpops.data(), n_atoms_, n_g, Ug, opt);
}
