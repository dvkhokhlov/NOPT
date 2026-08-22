# include "blas_link.h"

# include "avas.h"
# include "matr.h"
# include "libint_link.h"
# include "basis_lib_read.h"
# include "common_vars.h"
# include "defaults.h"

# include <vector>
# include <cstdio>
# include <cmath>

//------------------------------------------------------------------------------------------------------------------------
// sigma spectra beyond the selection boundary are printed only this far
#define AVAS_PRINT_TAIL 20
// smallest eigenvalue of the target overlap that still defines a target space
#define AVAS_S22_MIN 1.0e-6
// an unselected occupied sigma at least this large, summing to one with a selected virtual one,
// is one target direction shared between the two blocks
#define AVAS_SPLIT_MIN 0.05
#define AVAS_SPLIT_TOL 0.02

static const char * avas_l_labels = "spdfghik";

// Reference shells of the requested atoms carrying a requested nl label. Within one atom the
// k-th reference shell of angular momentum l is the principal number n = k+l+1. `below` collects
// the non-requested shells of the same atoms up to the highest requested n.
static std::vector<Shell> avas_ref_shells(molecule * M, const avas_par & A,
                                          const std::vector<Shell> & all,
                                          const std::vector<int> & center,
                                          std::vector<Shell> & below)
{
    std::vector<Shell> out;
    int n_kw = int(A.shell_n.size());

    int n_top = 0;
    for(int k=0;k<n_kw;k++)if(A.shell_n[k]>n_top)n_top=A.shell_n[k];

    for(int i_sel=0; i_sel<int(A.atoms.size()); i_sel++){

        int i_a = A.atoms[i_sel]-1;
        if(i_a>=M->n_atoms){
            fprintf(out_stream,"ERROR: $AVAS atom %d is out of range (the molecule has %d atoms)\n",
                                A.atoms[i_sel],M->n_atoms);
            exit(EXIT_FAILURE);
        }

        std::vector<int> found(n_kw,0);
        int n_l[8]={0,0,0,0,0,0,0,0};

        for(int i=0;i<int(all.size());i++){
            if(center[i]!=i_a)continue;
            int l = all[i].contr[0].l;
            if(l>7)continue;
            n_l[l]++;
            int n = n_l[l]+l;
            int is_target = 0;
            for(int k=0;k<n_kw;k++)
                if((A.shell_l[k]==l)&&(A.shell_n[k]==n)){
                    out.push_back(all[i]);
                    found[k]=1;
                    is_target = 1;
                }
            if((is_target==0)&&(n<=n_top))below.push_back(all[i]);
        }

        for(int k=0;k<n_kw;k++)
            if(found[k]==0){
                fprintf(out_stream,"ERROR: $AVAS reference shell %d%c is absent for atom %d (%s)\n",
                                    A.shell_n[k],avas_l_labels[A.shell_l[k]],i_a+1,M->atom_names[i_a]);
                fprintf(out_stream,"       basis %s, library %s\n",A.ref_basis.c_str(),NOPT_LIB);
                exit(EXIT_FAILURE);
            }
    }

    return out;
}

// A = (C S12) S22^-1 (C S12)^T for one MO row block; U <- eigenvectors in rows, sigma ascending.
static void avas_block_projector(const double * C, int n_b, int n_ao, int n_ref,
                                 const double * S12, const double * S22i,
                                 double * U, double * sigma)
{
    if(n_b<=0)return;

    std::vector<double> Mb(size_t(n_b)*n_ref);
    std::vector<double> P (size_t(n_b)*n_ref);

    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   n_b,n_ref,n_ao,1.0,
                   C,n_ao,
                   S12,n_ref,0.0,
                   Mb.data(),n_ref);

    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                   n_b,n_ref,n_ref,1.0,
                   Mb.data(),n_ref,
                   S22i,n_ref,0.0,
                   P.data(),n_ref);

    nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                   n_b,n_b,n_ref,1.0,
                   P.data(),n_ref,
                   Mb.data(),n_ref,0.0,
                   U,n_b);

    symmetrization(U,n_b);
    lapack_diag(U,sigma,n_b);
}

// index of the largest drop in a descending spectrum
static int avas_max_gap(const double * sigma, int n)
{
    int best=-1;
    double g=-1.0;

    for(int i=0;i+1<n;i++)
        if(sigma[i]-sigma[i+1]>g){
            g=sigma[i]-sigma[i+1];
            best=i;
        }

    return best;
}

static void avas_print_sigma(const char * title, const double * sigma, int n, int n_sel)
{
    int n_p = n_sel+AVAS_PRINT_TAIL;
    if(n_p>n)n_p=n;

    fprintf(out_stream,"%s (%d of %d selected)\n",title,n_sel,n);
    for(int i=0,c=0;i<n_p;i++,c++){
        if((i==n_sel)&&(i)){ fprintf(out_stream,"\n  ---- selection boundary ----\n"); c=0; }
        else if((i)&&(c==6)){ fprintf(out_stream,"\n"); c=0; }
        fprintf(out_stream," %11.8f",sigma[i]);
    }
    if(n_p<n)fprintf(out_stream,"\n  ... %d smaller values not printed",n-n_p);
    fprintf(out_stream,"\n\n");
}

//------------------------------------------------------------------------------------------------------------------------
int avas_steer(molecule * M, const avas_par & A, char * job_name)
{
    fprintf(out_stream,"\n\n\n");
    fprintf(out_stream,"_____________________Starting_AVAS_orbital_steering____________________\n\n");

    A.write_info();

    int n_ao  = M->n_ao;
    int n_occ = M->n_el_calc/2;
    int n_vir = M->n_mo-n_occ;
    int k_o   = n_occ-M->n_cor_orb;
    int k_v   = M->n_act_orb[0]-k_o;

    if((k_o<0)||(k_v<0)||(k_v>n_vir)){
        fprintf(out_stream,"ERROR: AVAS cannot fill the active window: %d occupied and %d virtual orbitals\n",k_o,k_v);
        fprintf(out_stream,"       are needed, %d and %d are available\n",n_occ,n_vir);
        exit(EXIT_FAILURE);
    }
    // several fragments interleave several active blocks, which is not one window to fill
    if(M->n_frag>1){
        fprintf(out_stream,"ERROR: AVAS supports a single active space (got %d fragments; N_MOL must be 1)\n",M->n_frag);
        exit(EXIT_FAILURE);
    }

    for(int i=0;i<M->n_mo;i++)
        for(int j=0;j<n_ao;j++)
            if(!std::isfinite(M->MO_VEC[size_t(i)*n_ao+j])){
                fprintf(out_stream,"ERROR: AVAS input orbitals are not finite (first bad MO %d)\n",i+1);
                fprintf(out_stream,"       RHF=0 needs _MO INP= orbitals covering core+active; supply them or set RHF=1\n");
                exit(EXIT_FAILURE);
            }

    std::vector<int>   ref_center;
    std::vector<Shell> ref_all = basis_lib_read_gbs(M,A.ref_basis.c_str(),1,0,
                                                    nullptr,&ref_center,true,nullptr,nullptr);
    std::vector<Shell> low_s;
    std::vector<Shell> ref_s   = avas_ref_shells(M,A,ref_all,ref_center,low_s);

    int n_ref=0;
    for(auto &sh: ref_s)n_ref+=sh.size();
    int n_low=0;
    for(auto &sh: low_s)n_low+=sh.size();

    fprintf(out_stream,"Reference functions:              %d\n",n_ref);
    fprintf(out_stream,"Orthogonalised against:           %d lower reference functions\n",n_low);
    fprintf(out_stream,"Active window:                    %d occupied + %d virtual\n\n",k_o,k_v);

    // AO_1el_from_2shells accumulates, so both overlaps start at zero
    std::vector<double> S22 (size_t(n_ref)*n_ref,0.0);
    std::vector<double> S12 (size_t(n_ao )*n_ref,0.0);
    AO_1el_from_2shells(S22.data(),ref_s,ref_s,n_ref,n_ref,'s',0);
    AO_1el_from_2shells(S12.data(),M->s ,ref_s,n_ao ,n_ref,'s',0);

    // In a segmented basis a requested shell is not orthogonal to the shells below it (def2-SVP Cr:
    // <3s|4s>=0.40, <3p|4p>=0.32), so span{target} holds semicore the virtual block cannot reach.
    // Removed over all requested atoms at once - per-atom removal strands more, not less.
    if(n_low){
        std::vector<double> SLL(size_t(n_low)*n_low,0.0);
        std::vector<double> SLT(size_t(n_low)*n_ref,0.0);
        std::vector<double> S1L(size_t(n_ao )*n_low,0.0);
        AO_1el_from_2shells(SLL.data(),low_s,low_s,n_low,n_low,'s',0);
        AO_1el_from_2shells(SLT.data(),low_s,ref_s,n_low,n_ref,'s',0);
        AO_1el_from_2shells(S1L.data(),M->s ,low_s,n_ao ,n_low,'s',0);

        inv_matr_constr(SLL.data(),n_low);

        std::vector<double> W(size_t(n_low)*n_ref);
        nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                       n_low,n_ref,n_low, 1.0,
                       SLL.data(),n_low,
                       SLT.data(),n_ref,0.0,
                       W  .data(),n_ref);

        nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                       n_ao,n_ref,n_low,-1.0,
                       S1L.data(),n_low,
                       W  .data(),n_ref,1.0,
                       S12.data(),n_ref);

        nopt_par_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                       n_ref,n_ref,n_low,-1.0,
                       SLT.data(),n_ref,
                       W  .data(),n_ref,1.0,
                       S22.data(),n_ref);
    }

    // a target left inside the shells below it has no direction of its own to project on
    std::vector<double> S22e(S22), s22_eig(n_ref);
    symmetrization(S22e.data(),n_ref);
    lapack_diag(S22e.data(),s22_eig.data(),n_ref);
    if(s22_eig[0]<AVAS_S22_MIN){
        fprintf(out_stream,"ERROR: $AVAS target overlap is singular, smallest eigenvalue %.3e\n",s22_eig[0]);
        fprintf(out_stream,"       a requested shell of %s is spanned by the shells below it\n",A.ref_basis.c_str());
        exit(EXIT_FAILURE);
    }

    std::vector<double> S22i(S22);
    inv_matr_constr(S22i.data(),n_ref);

    std::vector<double> U_o(size_t(n_occ)*n_occ), sig_o(n_occ);
    std::vector<double> U_v(size_t(n_vir)*n_vir), sig_v(n_vir);

    avas_block_projector(M->MO_VEC                    ,n_occ,n_ao,n_ref,S12.data(),S22i.data(),U_o.data(),sig_o.data());
    avas_block_projector(M->MO_VEC+size_t(n_occ)*n_ao ,n_vir,n_ao,n_ref,S12.data(),S22i.data(),U_v.data(),sig_v.data());

    // rotate the two row blocks in place; ascending sigma leaves the k_o selected occupieds as
    // the last occupied rows, and the reversal below leaves the k_v selected virtuals first
    std::vector<double> B(size_t(n_occ>n_vir?n_occ:n_vir)*n_ao);

    if(n_occ>0){
        cblas_dcopy(n_occ*n_ao,M->MO_VEC,1,B.data(),1);
        nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                       n_occ,n_ao,n_occ,1.0,
                       U_o.data(),n_occ,
                       B.data(),n_ao,0.0,
                       M->MO_VEC,n_ao);
    }

    if(n_vir>0){
        cblas_dcopy(n_vir*n_ao,M->MO_VEC+size_t(n_occ)*n_ao,1,B.data(),1);
        nopt_par_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                       n_vir,n_ao,n_vir,1.0,
                       U_v.data(),n_vir,
                       B.data(),n_ao,0.0,
                       M->MO_VEC+size_t(n_occ)*n_ao,n_ao);
    }

    for(int i=0;i<n_vir/2;i++){
        cblas_dswap(n_ao,M->MO_VEC+size_t(n_occ+i)*n_ao,1,M->MO_VEC+size_t(n_occ+n_vir-1-i)*n_ao,1);
        double t=sig_v[i]; sig_v[i]=sig_v[n_vir-1-i]; sig_v[n_vir-1-i]=t;
    }

    for(int i=0;i<n_occ;i++)M->orb_energy[i]       = sig_o[i];
    for(int a=0;a<n_vir;a++)M->orb_energy[n_occ+a] = sig_v[a];

    std::vector<double> sig_o_d(n_occ);
    for(int i=0;i<n_occ;i++)sig_o_d[i]=sig_o[n_occ-1-i];

    avas_print_sigma("AVAS occupied-block sigma (descending)",sig_o_d.data(),n_occ,k_o);
    avas_print_sigma("AVAS virtual-block  sigma (descending)",sig_v  .data(),n_vir,k_v);

    // sigma sums to the target dimension over both blocks, so what the window misses is the
    // target weight the active space cannot carry
    double w_sel=0.0, w_out_o=0.0, w_out_v=0.0;
    for(int i=0    ;i<k_o  ;i++)w_sel  +=sig_o_d[i];
    for(int a=0    ;a<k_v  ;a++)w_sel  +=sig_v  [a];
    for(int i=k_o  ;i<n_occ;i++)w_out_o+=sig_o_d[i];
    for(int a=k_v  ;a<n_vir;a++)w_out_v+=sig_v  [a];

    fprintf(out_stream,"Target weight in the active space: %.3f of %d\n",w_sel,n_ref);
    fprintf(out_stream,"          left in occupied orbitals: %.3f\n",w_out_o);
    fprintf(out_stream,"          left in virtual  orbitals: %.3f\n\n",w_out_v);

    // one target direction shared by an occupied and a virtual orbital sums to one over the pair,
    // and the window fixed by n_alp/n_bet can hold only the virtual half of such a pair
    int    n_split=0;
    double w_split=0.0;
    for(int i=k_o;i<n_occ;i++){
        if(sig_o_d[i]<AVAS_SPLIT_MIN)break;
        for(int a=0;a<k_v;a++)
            if(fabs(sig_o_d[i]+sig_v[a]-1.0)<AVAS_SPLIT_TOL){
                n_split++;
                if(sig_o_d[i]>w_split)w_split=sig_o_d[i];
                break;
            }
    }
    if(n_split)
        fprintf(out_stream,"NOTE: %d target directions are split between the two blocks; the largest leaves\n"
                           "      %.3f outside the occupied window that n_alp and n_bet fix\n",n_split,w_split);

    if((k_o>0)&&(k_o<n_occ)&&(avas_max_gap(sig_o_d.data(),n_occ)!=k_o-1))
        fprintf(out_stream,"NOTE: the occupied boundary after %d orbitals is not the largest gap of its spectrum\n",k_o);
    if((k_v>0)&&(k_v<n_vir)&&(avas_max_gap(sig_v.data(),n_vir)!=k_v-1))
        fprintf(out_stream,"NOTE: the virtual boundary after %d orbitals is not the largest gap of its spectrum\n",k_v);

    // a selected orbital with sigma ~ 0 lies in the projector kernel: the reference span is
    // exhausted and its identity is arbitrary within the block
    int nz_o=0,nz_v=0;
    for(int i=0;i<k_o;i++)if(sig_o_d[i]>1e-8)nz_o++;
    for(int a=0;a<k_v;a++)if(sig_v  [a]>1e-8)nz_v++;
    if(nz_o<k_o)
        fprintf(out_stream,"NOTE: only %d of the %d selected occupied orbitals overlap the reference "
                           "(%d functions); the rest are arbitrary within the occupied space\n",nz_o,k_o,n_ref);
    if(nz_v<k_v)
        fprintf(out_stream,"NOTE: only %d of the %d selected virtual orbitals overlap the reference "
                           "(%d functions); the rest are arbitrary within the virtual space\n",nz_v,k_v,n_ref);

    M->MO_gamess_format();

    char name[BUF_LINE_LENGTH];

    fprintf(out_stream,"\n");
    fprintf(out_stream,"Writing AVAS orbitals:\n");

    sprintf(name,"%s_AVAS.out",job_name);
    M->GAMESS_type_out_print(name,-1);
    fprintf(out_stream,"visualization file: %s\n",name);

    sprintf(name,"%s_AVAS.orb",job_name);
    M->MO_print(name);
    fprintf(out_stream,"data file         : %s\n",name);

    fprintf(out_stream,"\n");
    fprintf(out_stream,"_______________________________________________________________________\n\n\n");

    return 0;
}
