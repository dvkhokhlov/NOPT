// Spin-free g -> aldet RF spin-case map + the tensor dump writer. Rank <= 2
// tables come from the final g1/g2; rank-3 spin cases come from the RAW
// both-double tables (raw_av+raw_ca) via the certified read-off images
// (conventions.md §7 (ii)/(iii)), NEVER from g3. Outputs are aldet-native
// packings in the exact state PT2_import_data expects.

#include "cdas_sf_tensors.h"

#include <vector>
#include <cstdio>
#include <cstddef>

namespace {

// ½(X+Xᵀ) on an m×m row-major matrix — the stock E2_calc_EE finalization step.
void symmetrize_mm(double* M, int m){
    for(int i=0;i<m;i++)
    for(int j=i+1;j<m;j++){
        double s = 0.5*(M[(size_t)i*m+j] + M[(size_t)j*m+i]);
        M[(size_t)i*m+j] = M[(size_t)j*m+i] = s;
    }
}

} // namespace

// Builds the six aldet RF arrays from the final tensors + raw rank-3 tables.
// Caller owns storage: RF_PH n_a², RF_PV_* n_a⁴, RF_P3_* n_a⁶, RF_PS scalar.
int cdas_sf_to_rf(const cdas_sf_tensors& t,
                  double* RF_PS, double* RF_PH,
                  double* RF_PV_JK, double* RF_PV_AB,
                  double* RF_P3_JK, double* RF_P3_AB){
    const int n = t.n_a;
    const size_t n2=(size_t)n*n, n4=n2*n2, n6=n4*n2;
    auto i4 = [n](int t_,int u_,int v_,int w_){
        return (((size_t)t_*n+u_)*n+v_)*n+w_; };
    auto i6 = [n](int a,int b,int c,int d,int e,int f){
        return (((((size_t)a*n+b)*n+c)*n+d)*n+e)*n+f; };

    // --- rank 0/1: scalar and 1-body are identities (gate A: expected 1.0). ---
    *RF_PS = t.E0;
    for(size_t k=0;k<n2;k++) RF_PH[k] = t.g1[k];

    // --- rank 2: g2 -> aldet ABAB pairing. RF_PV_AB is the middle-index swap
    // (opposite-spin read-off, coeff g2_{tu,vw}); RF_PV_JK is the same-spin
    // antisymmetrized image in that pairing (gate A pins both factors). ---
    for(int a=0;a<n;a++) for(int b=0;b<n;b++)
    for(int c=0;c<n;c++) for(int d=0;d<n;d++){
        RF_PV_AB[i4(a,b,c,d)] = t.g2[i4(a,c,b,d)];
        RF_PV_JK[i4(a,b,c,d)] = t.g2[i4(a,d,b,c)] - t.g2[i4(a,c,b,d)];
    }

    // --- rank 3: raw both-double tables (empty class contributes zero). ---
    const double* rav = t.raw_av.empty() ? nullptr : t.raw_av.data();
    const double* rca = t.raw_ca.empty() ? nullptr : t.raw_ca.data();
    auto raw = [&](int a,int b,int c,int d,int e,int f)->double{
        double x=0.0; size_t k=i6(a,b,c,d,e,f);
        if(rav) x+=rav[k]; if(rca) x+=rca[k]; return x; };
    // read-off (ii): mixed AAB coeff = raw_{tu,vw,xy} − raw_{tw,vu,xy} (u↔w).
    auto mixed = [&](int t_,int u_,int v_,int w_,int x_,int y_)->double{
        return raw(t_,u_,v_,w_,x_,y_) - raw(t_,w_,v_,u_,x_,y_); };

    // RF_P3_AB: the mixed image is the coeff of the aldet AAB string in the
    // POST-re-index order [t,u,v,w,x,y]. Replicate the stock finalization frame
    // exactly — Hermitize (½(X+Xᵀ), n³) in the PRE-re-index native order
    // [t,u,x,v,w,y], then AABAAB re-index — so the table matches PT2_import_data.
    {
        // native[a,b,c,d,e,f] = mixed(a,d,b,e,c,f) (arrangement pinned against
        // the stock AV/CA builders, gate B), then Hermitize + AABAAB re-index.
        std::vector<double> nat(n6);
        for(int a=0;a<n;a++) for(int b=0;b<n;b++)
        for(int c=0;c<n;c++) for(int d=0;d<n;d++)
        for(int e=0;e<n;e++) for(int f=0;f<n;f++)
            nat[i6(a,b,c,d,e,f)] = mixed(a,d,b,e,c,f);
        symmetrize_mm(nat.data(), (int)(n2*n));       // n³×n³ Hermitize
        for(int t_=0;t_<n;t_++) for(int u_=0;u_<n;u_++)
        for(int v_=0;v_<n;v_++) for(int w_=0;w_<n;w_++)
        for(int x_=0;x_<n;x_++) for(int y_=0;y_<n;y_++)
            RF_P3_AB[i6(t_,u_,v_,w_,x_,y_)] = nat[i6(t_,u_,x_,v_,w_,y_)];
    }

    // RF_P3_JK (same-spin ααα): (1/3)·mixed — aldet's 9-term e3 sum supplies the
    // ×3. The stock table itself is a different gauge representative and is NOT
    // reproducible from raw_av; only the e3 weight is invariant. aldet also
    // self-contracts this table into T2_AA, and RF_PV_JK already carries the
    // full 2-body from g2 — that shadow is subtracted below (no double count).
    for(int a=0;a<n;a++) for(int b=0;b<n;b++)
    for(int c=0;c<n;c++) for(int d=0;d<n;d++)
    for(int e=0;e<n;e++) for(int f=0;f<n;f++)
        RF_P3_JK[i6(a,b,c,d,e,f)] = (1.0/3.0)*mixed(a,c,b,d,e,f);
    for(int i=0;i<n;i++) for(int j=0;j<n;j++)
    for(int k=0;k<n;k++) for(int l=0;l<n;l++){
        double sh=0.0;
        for(int m=0;m<n;m++) sh += RF_P3_JK[i6(j,i,m,l,k,m)];
        RF_PV_JK[i4(i,j,k,l)] -= sh;
    }

    return 0;
}

// Tensor dump (conventions.md §10 verbatim). External linkage; the driver
// reaches it through a local extern declaration (the frozen header is untouched).
// Binary little-endian, ASCII header lines, raw g1/g2/g3 in the §4 layouts.
void cdas_sf_write_dump(const char* path_prefix, const char* scheme,
                        const char* basis, double eps_A,
                        const cdas_sf_tensors& t){
    char name[4096];
    snprintf(name, sizeof name, "%s_sf_tensors.dump", path_prefix);
    FILE* fp = fopen(name, "wb");
    if(!fp){ fprintf(stderr,"cdas_sf_write_dump: cannot open %s\n", name); return; }
    const int n = t.n_a;
    const size_t n2=(size_t)n*n, n4=n2*n2, n6=n4*n2;
    fprintf(fp, "CDAS_SF_DUMP v1\n");
    fprintf(fp, "n_act %d  scheme %s  basis %s\n", n, scheme, basis);
    fprintf(fp, "eps_A %.17g\n", eps_A);
    fprintf(fp, "E0 %.17g\n", t.E0);
    fwrite(t.g1.data(), sizeof(double), n2, fp);
    fwrite(t.g2.data(), sizeof(double), n4, fp);
    fwrite(t.g3.data(), sizeof(double), n6, fp);
    fclose(fp);
    fprintf(stdout, "CDAS SF tensor dump written: %s\n", name);
}
