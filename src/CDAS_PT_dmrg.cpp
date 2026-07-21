// CDAS_PT_dmrg — DMRG-backed CDAS-PT2 driver (EE family). Owns a fresh CAS_engine, runs the bare
// SA-DMRG solve to freeze the localized lattice, builds the dressed operator in that lattice basis
// (spin-free EE tensor engine), then imports it and re-solves cold on the SAME engine. Talks only to
// the casci_solver virtual surface, so it compiles in every build config (dressed import is DMRG-only).

# include "blas_link.h"
# include "matr.h"
# include "timer.h"
# include "RI.h"
# include "inp_out.h"
# include "CAS.h"
# include "cdas_sf_tensors.h"
# include "localized_dmrg.h"
# include "CDAS_PT.h"

# include <vector>
# include <cmath>
# include <algorithm>
# include <cstdio>
# include <cstdlib>

// SF tensor dump writer (src/cdas_sf_spincase.cpp); the frozen header is untouched.
void cdas_sf_write_dump(const char* path_prefix, const char* scheme,
                        const char* basis, double eps_A,
                        const cdas_sf_tensors& t);

// Split an MO-row matrix V[orb*d_all+ao] into per-space column vectors [ao*n+orb]
// (core|active|virtual by row offset). Mirrors CDAS_PT.cpp's copy_MO_to_CVEC.
static void mo_rows_to_cvec(const double* V, int d1, int d2, int d3, int d_all,
                            double* O1, double* O2, double* O3){
    for(int j=0;j<d_all;j++) for(int i=0;i<d1;i++) O1[j*d1+i]=V[i*d_all+j];
    for(int i=0;i<d2;i++)    for(int j=0;j<d_all;j++) O2[j*d2+i]=V[(i+d1)*d_all+j];
    for(int j=0;j<d_all;j++) for(int i=0;i<d3;i++) O3[j*d3+i]=V[(i+d1+d2)*d_all+j];
}

int CDAS_PT_dmrg(molecule * M, cdas_par * cdas, char * job_name){

    if(RI==0){
        fprintf(out_stream,"WARNING: RI=0 is not supported for any PT\n");
        exit(0);
    }

    // --- 1. validate the input scheme FIRST (before any RI/engine work) ---
    // EE family only, and only the forms whose active energy is a single scalar. The active block
    // of M->orb_energy is zero after a DMRG-CASSCF, so any Fock-derived active energy is undefined.
    if(cdas->IPEA || cdas->MPPT){
        fprintf(out_stream,"ERROR: IPEA/MPPT CDAS-PT is not implemented for the DMRG backend"
                           " (EE family only: HOMO or constant ENERGY)\n");
        exit(EXIT_FAILURE);
    }
    if(cdas->actual || cdas->orb_e || cdas->mult_e || cdas->fit_e){
        fprintf(out_stream,"ERROR: DMRG-CDAS supports the EE family only (HOMO or constant ENERGY);"
                           " ACTUAL/ENERGIES/USE_ORB_FOR_ENERGY/USE_FIREFLY_FIT_ENERGY need active"
                           " Fock energies that are undefined under a DMRG CAS\n");
        exit(EXIT_FAILURE);
    }
    if(cdas->sing_e && (cdas->have_eps || cdas->eps.size()!=1)){
        fprintf(out_stream,"ERROR: DMRG-CDAS ENERGY scheme needs a single constant value\n");
        exit(EXIT_FAILURE);
    }
    if(!(cdas->HOMO || cdas->sing_e)){
        fprintf(out_stream,"ERROR: DMRG-CDAS needs an EE scheme (HOMO or constant ENERGY=<value>)\n");
        exit(EXIT_FAILURE);
    }
    if(fabs(cdas->edshift)>1E-8){
        fprintf(out_stream,"ERROR: DMRG-CDAS requires edshift=0 (bare resolvent)\n");
        exit(EXIT_FAILURE);
    }
    if(cdas->pt1_d)
        fprintf(out_stream,"NOTE: PT1 dipole correction not computed (DMRG-CDAS driver)\n");

    gen_RI_AA(M);
    cdas->write_info(M->n_act_el_alp[0], M->n_act_el_bet[0], M->n_act_orb[0], M->CI[0].mult);

    // --- 2. space sizes / handles (D5-aware) ---
    int n_ao  = M->n_ao;
    int n_act = M->n_act_orb[0];
    int n_cor = M->n_cor_orb;
    int n_mo  = n_ao;
    if(D5){ M->calc_d5_n_ao(); n_mo = M->n_ao_d5; }
    int n_virt = n_mo - n_cor - n_act;
    RI_core_realloc(n_cor+n_act, n_ao);

    // --- 3. fresh engine + bare SA-DMRG solve (freezes U_loc + Fiedler lattice) ---
    CAS_engine CAS;
    CAS.init(cdas->cas, M);
    CAS.calc_DM_C();
    CAS.CI_calc(1,0,0);

    // Snapshot the bare energies before any RDM call: the first RDM read overwrites the wrapper's
    // E_states cache with rdm_energy over the bare integrals (block2_dmrg.cpp).
    std::vector<double> E_bare(CAS.n_s);
    for(int i=0;i<CAS.n_s;i++) E_bare[i]=CAS.CI->E_state(i);

    // --- 4. eps_A (single EE scalar) ---
    double eps_A;
    if(cdas->sing_e){
        eps_A = cdas->eps[0];
    }
    else{
        // HOMO: eps_A = active-Fock eigenvalue #(n_alpha-1). Build the SA generalized Fock in the
        // driver (perform_diag=0 leaves M->orb_energy/MO_VEC pristine), then diagonalize its active
        // block for the eigenvalues only (no report rotation). SCF_alloc provides F_tot/GEN_CVEC/
        // MO_BUF -- init alone does not (the ENERGY scheme never needs them).
        CAS.SCF_alloc();
        CAS.av_DM_and_F_calc(0);
        std::vector<double> Fblk((size_t)CAS.n_act*CAS.n_act), ev(CAS.n_act);
        for(int i=0;i<CAS.n_act;i++)
        for(int j=0;j<CAS.n_act;j++)
            Fblk[i*CAS.n_act+j]=CAS.F_tot[(i+CAS.n_core)*CAS.n_ao+(j+CAS.n_core)];
        lapack_diag(Fblk.data(), ev.data(), CAS.n_act);
        eps_A = ev[M->n_act_el_alp[0]-1];
    }

    // --- 5. interlacing + intruder triad (hard contract for this driver) ---
    double * eps   = M->orb_energy;              // core|active(zero)|virtual, post-CASSCF
    double * eps_e = eps + n_cor + n_act;
    const bool hasc=(n_cor>0), hasv=(n_virt>0), has2c=(n_cor>=2), has2v=(n_virt>=2);
    double c1=eps_A, c2=eps_A, v1=eps_A, v2=eps_A;
    if(hasc){
        c1=eps[0]; c2=-1e300;
        for(int i=1;i<n_cor;i++){ double e=eps[i]; if(e>c1){ c2=c1; c1=e; } else if(e>c2) c2=e; }
    }
    if(hasv){
        v1=eps_e[0]; v2=1e300;
        for(int i=1;i<n_virt;i++){ double e=eps_e[i]; if(e<v1){ v2=v1; v1=e; } else if(e<v2) v2=e; }
    }
    fprintf(out_stream,"\nDMRG-CDAS EE intruder diagnostics:\n");
    fprintf(out_stream,"  max core eps    = ");   hasc? fprintf(out_stream,"%.6f\n",c1):fprintf(out_stream,"(none)\n");
    fprintf(out_stream,"  eps_A           = %.6f\n", eps_A);
    fprintf(out_stream,"  min virtual eps = ");   hasv? fprintf(out_stream,"%.6f\n",v1):fprintf(out_stream,"(none)\n");
    fprintf(out_stream,"  per-class min|D| (drive-to-zero labels):\n");
    auto pcl=[&](const char* nm, bool ok, double D){
        if(ok) fprintf(out_stream,"    %-5s = %.6f\n", nm, fabs(D));
        else   fprintf(out_stream,"    %-5s = (class empty)\n", nm); };
    pcl("CCVV", has2c&&has2v, c1+c2-v1-v2);
    pcl("CAVV", hasc &&has2v, c1+eps_A-v1-v2);
    pcl("AAVV", has2v,        2*eps_A-v1-v2);
    pcl("CCAV", has2c&&hasv,  c1+c2-eps_A-v1);
    pcl("CCAA", has2c,        c1+c2-2*eps_A);
    pcl("CV",   hasc &&hasv,  c1-v1);
    pcl("AV",   hasv,         eps_A-v1);
    pcl("CA",   hasc,         c1-eps_A);
    if((hasc && !(c1<eps_A)) || (hasv && !(eps_A<v1))){
        fprintf(out_stream,"\nERROR: EE active energy eps_A=%.6f violates interlacing"
                " (max core eps=%.6f, min virtual eps=%.6f); the EE resolvent denominators change"
                " sign -- CDAS-PT2 is not defined here\n", eps_A, hasc?c1:eps_A, hasv?v1:eps_A);
        exit(EXIT_FAILURE);
    }

    // --- 6. lattice-basis dressing inputs ---
    const bool has_loc = (bool)CAS.localizer_;

    // Private MO_VEC copy; active rows replaced by ACT_loc = U_loc^T * ACT (U_loc[a*n_act+p]).
    std::vector<double> MO_VEC_loc((size_t)n_ao*n_ao);
    std::copy(M->MO_VEC, M->MO_VEC+(size_t)n_ao*n_ao, MO_VEC_loc.begin());
    if(has_loc){
        double* ACT = M->MO_VEC + (size_t)n_cor*n_ao;   // (n_act x n_ao), read-only
        std::vector<double> ACT_loc((size_t)n_act*n_ao);
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                    n_act, n_ao, n_act, 1.0,
                    CAS.U_loc.data(), n_act,
                    ACT, n_ao,
                    0.0, ACT_loc.data(), n_ao);
        std::copy(ACT_loc.begin(), ACT_loc.end(), MO_VEC_loc.begin()+(size_t)n_cor*n_ao);
    }

    // Fresh RI_data for the lattice B-tensors (never re-MO_calc one RI_data -- it leaks).
    RI_data R_loc;
    R_loc.set_par(M, n_cor, n_act, n_virt);
    R_loc.MO_calc(MO_VEC_loc.data(), n_ao);

    // Core Fock F^c = H_AO + 2J - K from the core density; active-rotation-invariant, projected to
    // blocks with the localized active column vectors.
    std::vector<double> J((size_t)n_ao*n_ao), K((size_t)n_ao*n_ao);
    std::vector<double> act_scr((size_t)n_act*n_act*n_act*n_act);
    set_zero_matr(J.data(),(long)n_ao*n_ao);
    set_zero_matr(K.data(),(long)n_ao*n_ao);
    set_zero_matr(act_scr.data(),(long)n_act*n_act*n_act*n_act);
    calc_2el_MO_INTS_RI(n_ao, CAS.DM_C, J.data(), K.data(), act_scr.data(),
                        MO_VEC_loc.data(), MO_VEC_loc.data()+(size_t)n_ao*n_cor,
                        MO_VEC_loc.data()+(size_t)n_ao*n_cor, n_cor, n_act, 0);
    std::vector<double> H_core((size_t)n_ao*n_ao);
    for(size_t i=0;i<(size_t)n_ao*n_ao;i++) H_core[i]=M->H_AO[i]+2*J[i]-K[i];

    std::vector<double> COR_VEC((size_t)n_ao*n_cor), ACT_VEC((size_t)n_ao*n_act),
                        VIRT_VEC((size_t)n_ao*n_virt);
    mo_rows_to_cvec(MO_VEC_loc.data(), n_cor, n_act, n_virt, n_ao,
                    COR_VEC.data(), ACT_VEC.data(), VIRT_VEC.data());

    std::vector<double> H_AV((size_t)n_act*n_virt), H_CA((size_t)n_cor*n_act), H_CV((size_t)n_cor*n_virt);
    transform_from_col_MO(H_AV.data(), H_core.data(), n_ao, ACT_VEC.data(), n_act, VIRT_VEC.data(), n_virt);
    transform_from_col_MO(H_CA.data(), H_core.data(), n_ao, COR_VEC.data(), n_cor, ACT_VEC.data(), n_act);
    transform_from_col_MO(H_CV.data(), H_core.data(), n_ao, COR_VEC.data(), n_cor, VIRT_VEC.data(), n_virt);

    // --- 7. spin-free EE tensor build (gauge/projection/Hermitization all internal) ---
    cdas_sf_kernel Ksf;
    Ksf.e_IP.assign(n_act, eps_A);
    Ksf.e_EA.assign(n_act, eps_A);
    Ksf.deriv = 0;
    cdas_sf_tensors out;
    cdas_sf_build(R_loc, eps, n_cor, n_act, n_virt, H_AV.data(), H_CA.data(), H_CV.data(), Ksf, out);
    if(cdas->DUMP_TENSORS)
        cdas_sf_write_dump(job_name, "EE", "lattice", eps_A, out);
    printf_timer("DMRG-CDAS SF tensors");
    fflush(out_stream);

    // --- 8. total dressed operator (bare integrals rotated into the lattice basis + dressing) ---
    std::vector<double> h1((size_t)n_act*n_act);
    std::vector<double> h2((size_t)n_act*n_act*n_act*n_act);
    if(has_loc){
        rotate1(CAS.F_act,    CAS.U_loc.data(), n_act, h1.data(), true);
        rotate2(CAS.aaaa_ints, CAS.U_loc.data(), n_act, h2.data(), true);
    }
    else{
        std::copy(CAS.F_act,    CAS.F_act+(size_t)n_act*n_act, h1.begin());
        std::copy(CAS.aaaa_ints, CAS.aaaa_ints+(size_t)n_act*n_act*n_act*n_act, h2.begin());
    }
    for(size_t i=0;i<(size_t)n_act*n_act;i++)                  h1[i]+=out.g1[i];
    for(size_t i=0;i<(size_t)n_act*n_act*n_act*n_act;i++)      h2[i]+=out.g2[i];
    double const_total = CAS.E_core + out.E0;

    // --- 9. dressed import + cold dressed solve on the SAME engine (never a second CI_calc) ---
    if(!CAS.CI->supports_dressed_import()){
        fprintf(out_stream,"ERROR: the active CI backend does not support dressed-operator import"
                           " (need cisolver=dmrg)\n");
        exit(EXIT_FAILURE);
    }
    CAS.CI->import_dressed_operator(h1.data(), h2.data(), out.g3.data(), const_total);
    CAS.CI->solve(0,0,false);

    std::vector<double> E_dressed(CAS.n_s);
    for(int i=0;i<CAS.n_s;i++) E_dressed[i]=CAS.CI->E_state(i);

    const bool   hit_max = CAS.CI->last_solve_hit_max();
    const double resid   = CAS.CI->last_solve_resid();
    fprintf(out_stream,"\nBare SA-DMRG state energies:\n");
    for(int i=0;i<CAS.n_s;i++) fprintf(out_stream,"  bare[%3d] = %.10f\n", i, E_bare[i]);
    fprintf(out_stream,"Dressed DMRG solve: residual |dE| = %.3e, sweep-max reached: %s\n",
            resid, hit_max? "YES (under-converged -- raise $DMRG sweeps or m)":"no");

    // --- 10. output contract: the single dressed-energy summary, printed from the local copy ---
    fprintf(out_stream,"\n\nCDAS-PT2 Energy summary:\n");
    PrintEnergy(E_dressed.data(), CAS.n_s, 1);

    printf_timer("CDAS-PT2 (DMRG)");
    fflush(out_stream);
    return 0;
}
