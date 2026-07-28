#include "aldet_casci_wrap.h"

#include "blas_link.h"      // nopt_par_dgemm, CBLAS enums, num_threads
#include "matr.h"           // set_zero_matr

#include <vector>
#include <functional>

// Out-of-line key function: anchors aldet_casci_wrap's vtable in this single TU.
// The adapter owns no resources (ci_/par_ are non-owning), so the default suffices.
aldet_casci_wrap::~aldet_casci_wrap() = default;

// Full n_s x n_s spin-summed 2-RDM state matrix in the NOPT GAMMA convention
// (diagonal = per-state, off-diagonal = transition). Mirrors aldet_data::G_calc's
// OpenMP reduction but with the non-_diag workers. Accumulates (+=): the caller zeroes G,
// here only the private per-thread buffers are zeroed.
void aldet_casci_wrap::G_calc_full(double* G) {
    aldet_data* ci = ci_;
    const int  n_s   = ci->n_states[0];
    const int  N     = ci->n_act;
    const long block = (long)n_s * n_s * N * N * N * N;

    std::vector<double*> Gth(num_threads);
    std::vector<std::vector<double>> priv(num_threads > 1 ? num_threads - 1 : 0);
    Gth[0] = G;
    for (int i = 1; i < num_threads; i++) { priv[i - 1].assign(block, 0.0); Gth[i] = priv[i - 1].data(); }

    #pragma omp parallel num_threads(num_threads)
    {
        int nt = omp_get_thread_num();
        double* B = Gth[nt];
        aldet_calc_DM_2body_AA(B, n_s, n_s, ci->coef[0],     N, ci->na, ci->Na, ci->Nb, ci->fa, ci->vec_a, nt, num_threads);
        aldet_calc_DM_2body_AA(B, n_s, n_s, ci->coef_bas[0], N, ci->nb, ci->Nb, ci->Na, ci->fb, ci->vec_b, nt, num_threads);
        aldet_calc_DM_2body_AB(B, n_s, n_s, ci->coef[0],     N, ci->na, ci->nb, ci->Na, ci->Nb, ci->fa, ci->fb, ci->vec_a, ci->vec_b, nt, num_threads);
    }

    for (int i = 1; i < num_threads; i++)
        for (long j = 0; j < block; j++) Gth[0][j] += Gth[i][j];
}

namespace {

// One (Tbra, Tket, np) pair of the complementary overlap plus its per-root output.
struct h2caa_leg { const double* Tbra; const double* Tket; int np; double* omega; };

// Complementary six-operator overlap (3-RDM-free lambda3), per root, overwritten.
// omega[r] = sum_p sum_{spins} Tbra[p,w,x,y] Tket[p,z,u,v] <r| x+ y+ w z+ v u |r>
// = sum{ B2 . SF_G2 } - sum{ Tbra . D3 . Tket }, the certified normal-ordered moment.
// SF_G2 = per-root 2-RDM (GAMMA[x,u,y,v]); SF_D3 is spin-summed from the AAA/AAB workers.
void h2caa_kernel(aldet_data* ci, const h2caa_leg* leg, int nleg) {
    const int  n_s = ci->n_states[0];
    const int  no  = ci->n_act;
    const int  na2 = no * no;
    const int  na3 = no * no * no;
    const long na4 = (long)na3 * no;
    const long na6 = (long)na3 * na3;

    for (int l = 0; l < nleg; l++)
        for (int r = 0; r < n_s; r++) leg[l].omega[r] = 0.0;

    // Root-independent 2-body T-fold B2[(x,y),(u,v)] = sum_{p,w} Tbra[p,w,x,y] Tket[p,w,u,v],
    // viewing Tbra/Tket as [np*no x na2] (rows (p,w)); B2 = Tbra^T Tket.
    std::vector<std::vector<double>> B2(nleg);
    for (int l = 0; l < nleg; l++) {
        B2[l].assign((long)na2 * na2, 0.0);
        nopt_par_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, na2, na2, leg[l].np * no,
                       1.0, leg[l].Tbra, na2, leg[l].Tket, na2, 0.0, B2[l].data(), na2);
    }

    // Per-root 2-RDM diagonal blocks (n_s blocks of na4) in the NOPT GAMMA convention.
    std::vector<double> G2diag((long)n_s * na4, 0.0);
    ci->G_calc(G2diag.data());

    // Per-thread private-buffer reduction accumulating ONE worker's diagonal (r,r) block into dst
    // (size na6); dst is NOT zeroed here -- the caller zeroes it once, then sums the two spin copies.
    auto reduce_diag = [&](std::vector<double>& dst, const std::function<void(double*, int, int)>& w) {
        std::vector<std::vector<double>> priv(num_threads > 1 ? num_threads - 1 : 0);
        std::vector<double*> th(num_threads);
        th[0] = dst.data();
        for (int i = 1; i < num_threads; i++) { priv[i - 1].assign(na6, 0.0); th[i] = priv[i - 1].data(); }
        #pragma omp parallel num_threads(num_threads)
        {
            int nt = omp_get_thread_num();
            w(th[nt], nt, num_threads);
        }
        for (int i = 1; i < num_threads; i++)
            for (long j = 0; j < na6; j++) th[0][j] += th[i][j];
    };

    // 6-index flat address inside one na6 block.
    auto id6 = [no](int a, int b, int c, int d, int e, int f) -> long {
        return (((((long)a * no + b) * no + c) * no + d) * no + e) * no + f; };

    // Only diagonal (r,r) blocks are needed, so the roots loop OUTSIDE: each worker is invoked with
    // n_s=1, ld=n_states[0], base = coef[0]+r, which selects state r as element 0 of every
    // determinant's coefficient stride and writes just that block. Peak per-thread storage is na6,
    // not n_s^2*na6 (~9x smaller at n_s=3): one na6 pair (AAAsum = AAA_a+AAA_b, AABsum = AAB+BBA)
    // is reused across roots.
    std::vector<double> AAAsum(na6), AABsum(na6);
    std::vector<double> D3M(na6);
    std::vector<std::vector<double>> Z(nleg);
    for (int l = 0; l < nleg; l++) Z[l].assign((long)leg[l].np * na3, 0.0);

    for (int r = 0; r < n_s; r++) {
        set_zero_matr(AAAsum.data(), na6);
        reduce_diag(AAAsum, [&](double* o, int it, int nt) {
            aldet_calc_DM_3body_AAA(o, 1, n_s, ci->coef[0]     + r, no, ci->na, ci->Na, ci->Nb, ci->fa, ci->vec_a, it, nt); });
        reduce_diag(AAAsum, [&](double* o, int it, int nt) {
            aldet_calc_DM_3body_AAA(o, 1, n_s, ci->coef_bas[0] + r, no, ci->nb, ci->Nb, ci->Na, ci->fb, ci->vec_b, it, nt); });
        set_zero_matr(AABsum.data(), na6);
        reduce_diag(AABsum, [&](double* o, int it, int nt) {
            aldet_calc_DM_3body_AAB(o, 1, n_s, ci->coef[0]     + r, no, ci->na, ci->nb, ci->Na, ci->Nb, ci->fa, ci->fb, ci->vec_a, ci->vec_b, it, nt); });
        reduce_diag(AABsum, [&](double* o, int it, int nt) {
            aldet_calc_DM_3body_AAB(o, 1, n_s, ci->coef_bas[0] + r, no, ci->nb, ci->na, ci->Nb, ci->Na, ci->fb, ci->fa, ci->vec_b, ci->vec_a, it, nt); });

        const double* Gr = G2diag.data() + (long)r * na4;

        // Moment tensor D3M[(w,x,y),(z,u,v)] = SF_D3_r[x,y,z,u,v,w]. The two mixed-spin
        // families reordered into the (a,a,b) worker slots each cross one physical fermion
        // anticommutation (aab, aba), so they enter with a minus; aaa and the baa family do not.
        const double* A3 = AAAsum.data();
        const double* B3 = AABsum.data();
        for (int w = 0; w < no; w++)
        for (int x = 0; x < no; x++)
        for (int y = 0; y < no; y++)
        for (int z = 0; z < no; z++)
        for (int u = 0; u < no; u++)
        for (int v = 0; v < no; v++)
            D3M[(long)((w * no + x) * no + y) * na3 + ((z * no + u) * no + v)] =
                  A3[id6(x, y, z, w, v, u)]
                - B3[id6(x, y, z, w, u, v)]
                - B3[id6(x, z, y, v, u, w)]
                + B3[id6(y, z, x, w, v, u)];

        for (int l = 0; l < nleg; l++) {
            const double* Tbra = leg[l].Tbra;
            const double* Tket = leg[l].Tket;
            const int     np   = leg[l].np;

            // 2-body completion: omega_2 = sum B2[(x,y),(u,v)] * GAMMA_r[x,u,y,v].
            double o2 = 0.0;
            for (int x = 0; x < no; x++)
            for (int y = 0; y < no; y++)
            for (int u = 0; u < no; u++)
            for (int v = 0; v < no; v++)
                o2 += B2[l][(long)(x * no + y) * na2 + (u * no + v)]
                    * Gr[(((long)x * no + u) * no + y) * no + v];

            // 3-body contraction: Z = Tket . D3M^T, then omega_3 = -sum Tbra . Z.
            nopt_par_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, np, na3, na3,
                           1.0, Tket, na3, D3M.data(), na3, 0.0, Z[l].data(), na3);
            double o3 = 0.0;
            for (long k = 0; k < (long)np * na3; k++) o3 += Tbra[k] * Z[l][k];

            leg[l].omega[r] = o2 - o3;
        }
    }
}

}   // namespace

void aldet_casci_wrap::h2caa_overlap(const double* Tbra, const double* Tket, int np, double* omega) {
    const h2caa_leg leg = { Tbra, Tket, np, omega };
    h2caa_kernel(ci_, &leg, 1);
}

void aldet_casci_wrap::h2caa_overlap2(const double* Tbra1, const double* Tket1, int np1, double* omega1,
                                      const double* Tbra2, const double* Tket2, int np2, double* omega2) {
    const h2caa_leg leg[2] = { { Tbra1, Tket1, np1, omega1 }, { Tbra2, Tket2, np2, omega2 } };
    h2caa_kernel(ci_, leg, 2);
}
