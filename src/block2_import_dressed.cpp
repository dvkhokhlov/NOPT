// block2_import_dressed — encode a TOTAL dressed active-space operator (0/1/2/3-body) as one
// spin-adapted GeneralFCIDUMP -> GeneralMPO and swap it into the engine for a cold dressed solve.
// The bare fcidump/hamil stay intact so the RDM read-out (block2_dmrg.cpp) keeps working.

#include "block2_dmrg_engine.h"   // dmrgci_engine, block2 API, shared helpers

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "common_vars.h"          // out_stream
#include "matr.h"
#include "timer.h"

using namespace block2;
using namespace nopt_block2;

namespace {

// Gather one axis-permuted tensor: out[i,j,...] = in[perm[i], perm[j], ...] over `rank` axes of
// length n. This is block2's FCIDUMP::reorder direction (TInt/V1Int::reorder, integral.hpp:107/176:
// site i takes input orbital perm[i]), so the dressed tensors land on the bare solve's lattice.
static std::vector<double> reorder_gather(const double *in, int n, int rank,
                                          const std::vector<uint16_t> &perm) {
    size_t len = 1;
    for (int r = 0; r < rank; r++) len *= (size_t)n;
    std::vector<double> out(len);
    if (perm.empty()) { // no frozen order => identity copy
        std::copy(in, in + len, out.begin());
        return out;
    }
    std::vector<int> idx(rank);
    for (size_t f = 0; f < len; f++) {
        size_t rem = f;
        for (int r = rank - 1; r >= 0; r--) { idx[r] = (int)(rem % n); rem /= n; }
        size_t src = 0;
        for (int r = 0; r < rank; r++) src = src * n + perm[idx[r]];
        out[f] = in[src];
    }
    return out;
}

// Append one dense-tensor group (cutoff 0, axis->slot rperm). Mirrors the classic add path:
// push the coupling string, then add_sum_term over the full dense tensor.
static void add_group(const std::shared_ptr<GeneralFCIDUMP<double>> &gfd, const std::string &expr,
                      const std::vector<double> &T, int n, int rank, double factor,
                      std::vector<uint16_t> rperm) {
    std::vector<int> shape(rank, n);
    std::vector<size_t> strides(rank);
    size_t s = 1;
    for (int r = rank - 1; r >= 0; r--) { strides[r] = s; s *= (size_t)n; }
    gfd->exprs.push_back(expr);
    gfd->add_sum_term(T.data(), s, shape, strides, 0.0, factor, {}, rperm);
}

// Null or all-zero 3-body tensor => skip the 3-body group.
static bool all_zero(const double *T, size_t len) {
    if (T == nullptr) return true;
    for (size_t i = 0; i < len; i++)
        if (T[i] != 0.0) return false;
    return true;
}

} // namespace

void block2_casci_wrap::import_dressed_operator(const double *h1_total, const double *h2_total,
                                                const double *h3_total, double const_total) {
    dmrgci_engine &e = *impl_;
    const int n = e.n_act;

    // A Fiedler lattice must already be frozen by the bare import; without it the dressed tensors
    // have no site->orbital map to land on.
    if (e.cfg.loc_order == DMRG_LOCORDER_FIEDLER && e.reorder_perm.empty()) {
        fprintf(out_stream, "ERROR: dressed import before bare import (no frozen Fiedler order)\n");
        exit(EXIT_FAILURE);
    }

    // Reorder every axis of each tensor onto the frozen lattice (identity if reorder_perm empty).
    std::vector<double> h1 = reorder_gather(h1_total, n, 2, e.reorder_perm);
    std::vector<double> h2 = reorder_gather(h2_total, n, 4, e.reorder_perm);
    const size_t len3 = (size_t)n * n * n * n * n * n;
    const bool have_3body = !all_zero(h3_total, len3);
    std::vector<double> h3;
    if (have_3body)
        h3 = reorder_gather(h3_total, n, 6, e.reorder_perm);

    // One spin-adapted GeneralFCIDUMP; factors 2^{k/2}/k! complete the 1/k! generator sums.
    auto gfd = std::make_shared<GeneralFCIDUMP<double>>(ElemOpTypes::SU2);
    gfd->const_e = const_total;
    add_group(gfd, "(C+D)0", h1, n, 2, std::sqrt(2.0), {});
    add_group(gfd, "((C+(C+D)0)1+D)0", h2, n, 4, 1.0, {0, 3, 1, 2});
    if (have_3body)
        add_group(gfd, "((C+((C+(C+D)0)1+D)0)1+D)0", h3, n, 6, std::sqrt(2.0) / 3.0,
                  {0, 5, 1, 4, 2, 3});
    std::shared_ptr<GeneralFCIDUMP<double>> afd = gfd->adjust_order();

    // General site Hamiltonian (same vacuum/n/orbsym as the bare path) -> exact FastBipartite MPO.
    SU2 vacuum(0);
    std::vector<typename SU2::pg_t> orbsym(n, 0); // C1
    auto hamil = std::make_shared<GeneralHamiltonian<SU2, double>>(vacuum, n, orbsym);
    auto gmpo = std::make_shared<GeneralMPO<SU2, double>>(
        hamil, afd, MPOAlgorithmTypes::FastBipartite, 0.0, -1, 0);
    gmpo->build();
    std::shared_ptr<MPO<SU2, double>> mpo = std::make_shared<SimplifiedMPO<SU2, double>>(
        gmpo, std::make_shared<Rule<SU2, double>>(), false, false);
    // const_e reaches the sweep once (added at sweep_algorithm.hpp:717 after the const-free eigs);
    // IdentityAddedMPO copies const_e unchanged and injects only a coeff-1 identity operator, so it
    // adds no second constant in a sweep, while making the operator usable for expectations.
    mpo = std::make_shared<IdentityAddedMPO<SU2, double>>(mpo);

    e.mpo = mpo; // the existing solve() runs the dressed MPO cold (have_rotation stays false)
}

void block2_casci_wrap::PT2_import_data(double * ext_T3,
                                       double * ext_T3_AB,
                                       double * ext_T2,
                                       double * ext_T2_AB,
                                       double * ext_T1,
                                       double   ext_T0){
    
    int error=0;
    if(g1.size()!=n_act_*n_act_            ) error=1;
    if(g2.size()!=n_act_*n_act_*n_act_*n_act_) error=1;
    if(error){
        fprintf(out_stream,"ERROR: block2_casci_wrap::PT2_import_data found inconsistent n_act\n");
        fprintf(out_stream,"       check your code\n");
        exit(0);
    }
    // std::vector<double> h1((size_t)n_act_*n_act_);
    // std::vector<double> h2((size_t)n_act_*n_act_*n_act_*n_act_);
    // std::copy(H_AA,    H_AA+(size_t)n_act_*n_act_, h1.begin());
    // std::copy(act_INTS, act_INTS+(size_t)n_act_*n_act_*n_act_*n_act_, h2.begin());
        
    for(size_t i=0;i<(size_t)n_act_*n_act_;i++)g1[i]+=ext_T1[i];
    
    for(int a=0;a<n_act_;a++) 
    for(int c=0;c<n_act_;c++) 
    for(int b=0;b<n_act_;b++)
    for(int d=0;d<n_act_;d++){
        g2[((a*n_act_+c)*n_act_+b)*n_act_+d]+=ext_T2_AB[((a*n_act_+b)*n_act_+c)*n_act_+d];
    }
    
    // #pragma omp parallel for collapse(2) schedule(static)
    g3.resize(n_act_*n_act_*n_act_*n_act_*n_act_*n_act_);
    for(int t=0;t<n_act_;t++) 
    for(int u=0;u<n_act_;u++)
    for(int v=0;v<n_act_;v++) 
    for(int w=0;w<n_act_;w++)
    for(int x=0;x<n_act_;x++) 
    for(int y=0;y<n_act_;y++){
       g3[((((t*n_act_+u)*n_act_+v)*n_act_+w)*n_act_+x)*n_act_+y] = 
              ( 2.0*ext_T3_AB[((((t*n_act_+v)*n_act_+u)*n_act_+w)*n_act_+x)*n_act_+y] 
              + 2.0*ext_T3_AB[((((t*n_act_+x)*n_act_+u)*n_act_+y)*n_act_+v)*n_act_+w]
              + 2.0*ext_T3_AB[((((v*n_act_+t)*n_act_+w)*n_act_+u)*n_act_+x)*n_act_+y] 
              + 2.0*ext_T3_AB[((((v*n_act_+x)*n_act_+w)*n_act_+y)*n_act_+t)*n_act_+u]
              + 2.0*ext_T3_AB[((((x*n_act_+t)*n_act_+y)*n_act_+u)*n_act_+v)*n_act_+w] 
              + 2.0*ext_T3_AB[((((x*n_act_+v)*n_act_+y)*n_act_+w)*n_act_+t)*n_act_+u]
              -     ext_T3_AB[((((t*n_act_+v)*n_act_+u)*n_act_+w)*n_act_+x)*n_act_+y] 
              +     ext_T3_AB[((((v*n_act_+t)*n_act_+u)*n_act_+w)*n_act_+x)*n_act_+y]
              +     ext_T3_AB[((((x*n_act_+v)*n_act_+u)*n_act_+w)*n_act_+t)*n_act_+y] 
              +     ext_T3_AB[((((t*n_act_+x)*n_act_+u)*n_act_+w)*n_act_+v)*n_act_+y]
              -     ext_T3_AB[((((v*n_act_+x)*n_act_+u)*n_act_+w)*n_act_+t)*n_act_+y] 
              -     ext_T3_AB[((((x*n_act_+t)*n_act_+u)*n_act_+w)*n_act_+v)*n_act_+y] ) / 12.0;
    }
    int n=n_act_;
    auto ix = [n](int t,int u,int v,int w,int x,int y)->size_t {
    return (((((size_t)t*n+u)*n+v)*n+w)*n+x)*n+y; };
    #pragma omp parallel for collapse(2) schedule(static)
        for(int t=0;t<n;t++) for(int u=0;u<n_act_;u++)
        for(int v=0;v<n_act_;v++) for(int w=0;w<n_act_;w++)
        for(int x=0;x<n_act_;x++) for(int y=0;y<n_act_;y++){
            const size_t i = ix(t,u,v,w,x,y), id = ix(u,t,w,v,y,x);
            if(i < id){ double s = 0.5*(g3[i]+g3[id]); g3[i] = g3[id] = s; }
        }
    
    
    
    double const_total = g0 + ext_T0;
    
    import_dressed_operator(g1.data(), g2.data(), g3.data(), const_total);
    
    
    return;
}

int average_DM_aldet_diag(double * G_out, double * G, std::vector<double> avecoe,int na_p, int n_s);

int block2_casci_wrap::calc_IPEA_single(double * U_IP, double * H_IP, 
                                        double * U_EA, double * H_EA,
                                        int a, std::vector<double> avecoe) {
     
    int n_s = n_states();
    
    //U_IP
    double * gamma = new double[n_s*n_act_*n_act_];
    set_zero_matr(gamma,n_act_*n_act_*n_s);
    calc_DM_diag(gamma,a);
    for(int i=0;i<n_s*n_act_*n_act_;i++)gamma[i]=gamma[i]*0.5;
    // memcpy(U_IP,gamma,sizeof(double)*n_act_*n_act_);
    average_DM_aldet_diag(U_IP,gamma, avecoe,n_act_*n_act_,n_s);
    
    
    //U_EA
    for(int i=0;i<n_act_*n_act_;i++)U_EA[i]=-U_IP[i];
    for(int i=0;i<n_act_;i++)U_EA[i*(n_act_+1)]=1.0+U_EA[i*(n_act_+1)];
    
    
    double * GAMMA = new double[n_s*n_act_*n_act_*n_act_*n_act_];
    set_zero_matr(GAMMA,n_s*n_act_*n_act_*n_act_*n_act_);
    G_calc(GAMMA);
    
    //H_IP
    set_zero_matr(H_IP,n_act_*n_act_);

    for(int t=0;t<n_act_;t++)
    for(int u=0;u<n_act_;u++)
    for(int w=0;w<n_act_;w++)
        H_IP[t*n_act_+u]+= U_IP[t*n_act_+w]*g1[u*n_act_+w];//restricted variant
    
    
    for(int t=0;t<n_act_;t++)
    for(int u=0;u<n_act_;u++)
    for(int v=0;v<n_act_;v++)
    for(int x=0;x<n_act_;x++)
    for(int y=0;y<n_act_;y++)
        H_IP[t*n_act_+u]+=GAMMA[((t*n_act_+y)*n_act_+v)*n_act_+x]*0.5*
                      (g2[((u*n_act_+y)*n_act_+v)*n_act_+x]);//restricted variant
    
    
    // fprintf(out_stream,"H:\n");
    // PrintMatr(H_IP,n_act_,n_act_,0);
    
    //H_EA
    set_zero_matr(H_EA,n_act_*n_act_);

    for(int t=0;t<n_act_;t++)
    for(int u=0;u<n_act_;u++)
    for(int v=0;v<n_act_;v++)
        H_EA[t*n_act_+u]+= U_EA[t*n_act_+v]*g1[u*n_act_+v];//restricted variant
    
    for(int t=0;t<n_act_;t++)
    for(int u=0;u<n_act_;u++)
    for(int v=0;v<n_act_;v++)
    for(int x=0;x<n_act_;x++){
        H_EA[t*n_act_+u]+=U_IP[x*n_act_+v]*(2*g2[((t*n_act_+u)*n_act_+v)*n_act_+x]-g2[((t*n_act_+x)*n_act_+v)*n_act_+u]);
    }
    
    for(int t=0;t<n_act_;t++)
    for(int u=0;u<n_act_;u++)
    for(int v=0;v<n_act_;v++)
    for(int w=0;w<n_act_;w++)
    for(int x=0;x<n_act_;x++)
        H_EA[t*n_act_+u]-=GAMMA[((v*n_act_+t)*n_act_+x)*n_act_+w]*0.5*
                      (      g2[((v*n_act_+u)*n_act_+w)*n_act_+x]);//restricted variant
    
   
    // fprintf(out_stream,"H_EA:\n");
    // PrintMatr(H_EA,n_act_,n_act_,0);
    // exit(0);
    
    printf_timer("calculation of IPEA matrices");
    delete[] GAMMA;
    delete[] gamma;
    return 0;
}
           


