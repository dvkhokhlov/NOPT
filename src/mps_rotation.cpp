// mps_rotation — externally-driven multi-root TangentSpace TE rotation of a state-averaged MultiMPS.
// See include/mps_rotation.h. This is a faithful reduction of block2's
// TimeEvolution::update_multi_two_dot (mode = TangentSpace, serial, noise = 0), with the
// size==2-locked effective apply replaced by a per-root real Krylov expo loop -- no block2 change.

#include "mps_rotation.h"

#include <cmath>
#include <vector>

using namespace block2;

namespace {

// exp(t*(H_eff+const_e)) applied independently to each root's center wavefunction. Equivalent to
// what block2's EffectiveFunctions::expo_apply does for size==2 (a real generator decouples the
// re/im channels), generalized to arbitrary nroots. Returns the total (summed-over-roots) norm^2.
double expo_apply_multi_roots(
    const std::shared_ptr<EffectiveHamiltonian<SU2, double, MultiMPS<SU2, double>>> &h_eff,
    double t, double const_e, double conv_thrd, int krylov_size) {
    h_eff->precompute();
    double anorm =
        MatrixFunctions::norm(MatrixRef(h_eff->diag->data, (MKL_INT)h_eff->diag->total_memory, 1));
    const bool tasked = (h_eff->tf->opf->seq->mode == SeqTypes::Auto) ||
                        (h_eff->tf->opf->seq->mode & SeqTypes::Tasked);
    auto g = [&h_eff, tasked](const GMatrix<double> &a, const GMatrix<double> &b) {
        if (tasked)
            h_eff->tf->operator()(a, b, (double)1.0);
        else
            (*h_eff)(a, b, 0, (double)1.0);
    };
    double nsq = 0.0;
    for (int r = 0; r < (int)h_eff->ket.size(); r++) {
        GMatrix<double> v(h_eff->ket[r]->data, (MKL_INT)h_eff->ket[r]->total_memory, 1);
        IterativeMatrixFunctions<double>::expo_apply(
            g, t, anorm, v, const_e, /*symmetric=*/false, /*iprint=*/false,
            (std::shared_ptr<ParallelCommunicator<SU2>>)nullptr, conv_thrd, krylov_size);
        double nr = MatrixFunctions::norm(v);
        nsq += nr * nr;
    }
    h_eff->post_precompute();
    return nsq;
}

// One two-site TangentSpace TE update for a MultiMPS of arbitrary nroots. Returns the total norm^2
// of the (forward) evolved center before truncation.
double multi_two_dot_ext(const std::shared_ptr<MovingEnvironment<SU2, double, double>> &me, int i,
                         bool forward, double beta, ubond_t bond_dim, double cutoff,
                         double conv_thrd, int krylov_size) {
    auto mket = std::dynamic_pointer_cast<MultiMPS<SU2, double>>(me->ket);
    frame_<double>()->activate(0);
    if (mket->tensors[i] != nullptr || mket->tensors[i + 1] != nullptr)
        MovingEnvironment<SU2, double, double>::contract_multi_two_dot(i, mket);
    else {
        mket->load_tensor(i);
        mket->tensors[i] = mket->tensors[i + 1] = nullptr;
    }
    std::vector<std::shared_ptr<SparseMatrixGroup<SU2, double>>> old_wfns = mket->wfns;
    auto h_eff = me->multi_eff_ham(FuseTypes::FuseLR, forward, true);
    double nsq = expo_apply_multi_roots(h_eff, -beta, me->mpo->const_e, conv_thrd, krylov_size);
    h_eff->deallocate();

    std::vector<double> wfn_spectra;
    std::shared_ptr<SparseMatrix<SU2, double>> dm =
        MovingEnvironment<SU2, double, double>::density_matrix_with_multi_target(
            mket->info->vacuum, old_wfns, mket->weights, forward, 0.0, NoiseTypes::None);
    MovingEnvironment<SU2, double, double>::multi_split_density_matrix(
        dm, old_wfns, (int)bond_dim, forward, false, mket->wfns,
        forward ? mket->tensors[i] : mket->tensors[i + 1], cutoff, false, wfn_spectra,
        TruncationTypes::Physical);

    std::shared_ptr<StateInfo<SU2>> info = nullptr;
    if (forward) {
        info = mket->tensors[i]->info->extract_state_info(forward);
        mket->info->bond_dim = std::max(mket->info->bond_dim, (ubond_t)info->n_states_total);
        mket->info->left_dims[i + 1] = info;
        mket->info->save_left_dims(i + 1);
        mket->canonical_form[i] = 'L';
        mket->canonical_form[i + 1] = 'M';
    } else {
        info = mket->tensors[i + 1]->info->extract_state_info(forward);
        mket->info->bond_dim = std::max(mket->info->bond_dim, (ubond_t)info->n_states_total);
        mket->info->right_dims[i + 1] = info;
        mket->info->save_right_dims(i + 1);
        mket->canonical_form[i] = 'M';
        mket->canonical_form[i + 1] = 'R';
    }
    info->deallocate();
    if (forward) {
        mket->save_wavefunction(i + 1);
        mket->save_tensor(i);
        mket->unload_wavefunction(i + 1);
        mket->unload_tensor(i);
    } else {
        mket->save_tensor(i + 1);
        mket->save_wavefunction(i);
        mket->unload_tensor(i + 1);
        mket->unload_wavefunction(i);
    }
    dm->info->deallocate();
    dm->deallocate();
    for (int k = mket->nroots - 1; k >= 0; k--)
        old_wfns[k]->deallocate();
    old_wfns[0]->deallocate_infos();

    // TangentSpace back-propagation: exp(+beta*H_eff) on the moved single-site center.
    if (forward && i + 1 != me->n_sites - 1) {
        me->move_to(i + 1, true);
        mket->load_wavefunction(i + 1);
        auto k_eff = me->multi_eff_ham(FuseTypes::FuseR, forward, true);
        expo_apply_multi_roots(k_eff, beta, me->mpo->const_e, conv_thrd, krylov_size);
        k_eff->deallocate();
        mket->save_wavefunction(i + 1);
        mket->unload_wavefunction(i + 1);
    } else if (!forward && i != 0) {
        me->move_to(i - 1, true);
        mket->load_wavefunction(i);
        auto k_eff = me->multi_eff_ham(FuseTypes::FuseL, forward, true);
        expo_apply_multi_roots(k_eff, beta, me->mpo->const_e, conv_thrd, krylov_size);
        k_eff->deallocate();
        mket->save_wavefunction(i);
        mket->unload_wavefunction(i);
    }
    MovingEnvironment<SU2, double, double>::propagate_multi_wfn(i, 0, me->n_sites, mket, forward,
                                                                me->mpo->tf->opf->cg);
    mket->save_data();
    return nsq;
}

// One full TangentSpace TE sweep for a MultiMPS of arbitrary nroots. Returns the last update's
// total norm^2 (a unitary-propagation accuracy gauge).
double multi_te_sweep(const std::shared_ptr<MovingEnvironment<SU2, double, double>> &me,
                      bool forward, double beta, ubond_t bond_dim, double cutoff, double conv_thrd,
                      int krylov_size) {
    me->prepare();
    std::vector<int> sweep_range;
    if (forward)
        for (int it = me->center; it < me->n_sites - me->dot + 1; it++)
            sweep_range.push_back(it);
    else
        for (int it = me->center; it >= 0; it--)
            sweep_range.push_back(it);
    double nsq = 0.0;
    for (int i : sweep_range) {
        me->move_to(i);
        nsq = multi_two_dot_ext(me, i, forward, beta, bond_dim, cutoff, conv_thrd, krylov_size);
    }
    return nsq;
}

} // namespace

double evolve_sa_multimps(const std::shared_ptr<MultiMPS<SU2, double>> &mps,
                          const std::shared_ptr<MPO<SU2, double>> &mpo_rot, ubond_t rot_m, double dt,
                          int n_steps) {
    auto saved_seq = threading_()->seq_type;
    threading_()->seq_type = SeqTypes::Simple; // Tasked corrupts time evolution
    auto rme = std::make_shared<MovingEnvironment<SU2, double, double>>(mpo_rot, mps, mps, "ROT");
    rme->init_environments(false);
    double nsq = mps->nroots; // if n_steps == 0, treat as unrotated (mean per-root norm^2 = 1)
    for (int s = 0; s < n_steps; s++) {
        bool forward = (mps->center == 0);
        nsq = multi_te_sweep(rme, forward, -dt, rot_m, 1e-20, 1e-12, 40); // -dt: toward the new basis
    }
    rme->remove_partition_files();
    threading_()->seq_type = saved_seq;
    return mps->nroots > 0 ? nsq / mps->nroots : nsq; // mean per-root norm^2
}

mps_rotation_result apply_orbital_rotation_mps(
    const std::shared_ptr<MultiMPS<SU2, double>> &mps, const double *U, int n, int n_elec, int twos,
    const std::vector<uint8_t> &orbsym, const std::vector<uint16_t> &reorder_perm,
    int rot_m, int rot_steps) {
    const size_t nn = (size_t)n * n;
    mps_rotation_result res;

    // kappa = log(U) (real antisymmetric for a proper rotation). Log/complex temporaries live on the
    // heap, not block2's LIFO stack, which the subsequent environment build reads from.
    auto heap = std::make_shared<VectorAllocator<double>>();
    std::vector<double> Um(U, U + nn);
    ComplexMatrixRef ck(nullptr, n, n);
    ck.allocate(heap);
    ck.clear();
    ComplexMatrixFunctions::fill_complex(ck, MatrixRef(Um.data(), n, n), MatrixRef(nullptr, n, n));
    ComplexMatrixFunctions::logarithm(ck);
    std::vector<double> kre(nn), kim(nn);
    ComplexMatrixFunctions::extract_complex(ck, MatrixRef(kre.data(), n, n), MatrixRef(kim.data(), n, n));
    ck.deallocate(heap);
    res.im_norm = MatrixFunctions::norm(MatrixRef(kim.data(), n, n));
    if (res.im_norm > 1e-6) { // non-real generator -> not a proper rotation
        res.complex_generator = true;
        return res;
    }

    // NC MPO expects kappa^T; reindex into the frozen lattice order (site i holds orbital perm[i]).
    const bool have_perm = !reorder_perm.empty();
    std::vector<double> kap(nn);
    double kmax = 0.0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            int oi = have_perm ? reorder_perm[i] : i;
            int oj = have_perm ? reorder_perm[j] : j;
            kap[(size_t)i * n + j] = kre[(size_t)oj * n + oi];
            kmax = std::max(kmax, std::fabs(kap[(size_t)i * n + j]));
        }
    // Below threshold exp(-kappa) ~ I: the MPS already sits in the target basis, and a numerically-zero
    // generator makes a degenerate rotation MPO that block2's environment build reads past. Skip.
    if (kmax < 1e-8) {
        res.skipped = true;
        return res;
    }

    // one-body rotation MPO exp(-kappa t), built from the (lattice-order) generator
    auto fd_rot = std::make_shared<FCIDUMP<double>>();
    fd_rot->initialize_h1e(n, n_elec, twos, /*isym=*/0, 0.0, kap.data(), nn);
    auto hamil_rot = std::make_shared<HamiltonianQC<SU2, double>>(SU2(0), n, orbsym, fd_rot);
    std::shared_ptr<MPO<SU2, double>> mpo_rot =
        std::make_shared<MPOQC<SU2, double>>(hamil_rot, QCTypes::NC);
    mpo_rot->basis = hamil_rot->basis;
    mpo_rot = std::make_shared<SimplifiedMPO<SU2, double>>(
        mpo_rot,
        std::make_shared<AntiHermitianRuleQC<SU2, double>>(std::make_shared<RuleQC<SU2, double>>()),
        true);
    res.norm2 = evolve_sa_multimps(mps, mpo_rot, (ubond_t)rot_m, 1.0 / rot_steps, rot_steps);
    mpo_rot->deallocate();
    return res;
}
