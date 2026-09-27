# include "blas_link.h"

# include <algorithm>

# include "sad_guess.h"
# include "matr.h"
# include "libint_link.h"
# include "basis_lib_read.h"
# include "common_vars.h"
# include "RI.h"

#define SAD_BASIS  "cc-pvtz-minao"
#define SAD_Z_MAX  36 //last element of SAD_BASIS
#define SAD_L_MAX   3
#define SAD_N_MAX   8

//electrons of a neutral atom over the (n,l) subshells, Madelung order (n+l, then n)
static int atomic_config(int Z, int occ[SAD_L_MAX+1][SAD_N_MAX+1]){

    for(int l=0;l<=SAD_L_MAX;l++)
    for(int n=0;n<=SAD_N_MAX;n++)
        occ[l][n]=0;

    int e=Z;
    for(int s=1; (s<=2*SAD_N_MAX)&&(e>0); s++)
    for(int l=(s-1)/2; (l>=0)&&(e>0); l--){
        int n=s-l;
        if((l>SAD_L_MAX)||(n>SAD_N_MAX))continue;
        occ[l][n]=std::min(e,2*(2*l+1));
        e-=occ[l][n];
    }

    if(e)return 1;

    //aufbau exceptions below SAD_Z_MAX
    if((Z==24)||(Z==29)){
        occ[0][4]-=1;
        occ[2][3]+=1;
    }

    return 0;
}

int sad_guess(molecule * M){

    for(int i_a=0;i_a<M->n_atoms;i_a++)
        if(M->nucl_charges_full[i_a]>SAD_Z_MAX+0.5){
            fprintf(out_stream,"NOTE: SAD guess covers H-Kr only, using the Huckel guess\n");
            return 0;
        }

    if(M->PP.size()){
        fprintf(out_stream,"NOTE: SAD guess does not support ECP, using the Huckel guess\n");
        return 0;
    }

    fprintf(out_stream,"                      using SAD orbitals\n");
    sprintf(M->read_basis_name,"%s",SAD_BASIS);

    std::vector<int> shell_center_r;

    M->n_ro=0;
    M->read_s = basis_lib_read_gbs(M,SAD_BASIS,1,0,nullptr,&shell_center_r, true, nullptr,nullptr);
    for(auto &s: M->read_s)M->n_ro+=s.size();

    int n_ro = M->n_ro;

    M->MO_VEC_R = new double[n_ro*n_ro];

    double * S_22 = new double[n_ro*n_ro];
    double * DM   = new double[n_ro*n_ro];
    double * NO   = new double[n_ro*n_ro];
    double * BUF  = new double[n_ro*n_ro];
    double * occ  = new double[n_ro];

    set_zero_matr(S_22, n_ro*n_ro);
    set_zero_matr(DM  , n_ro*n_ro);
    AO_1el_from_2shells(S_22,M->read_s,M->read_s,n_ro,n_ro,'s',0);

    int i_s=0;
    int b_0=0;
    for(int i_a=0;i_a<M->n_atoms;i_a++){

        if(M->nucl_charges_full[i_a]==0)continue;

        int Z=int(M->nucl_charges_full[i_a]);
        int config[SAD_L_MAX+1][SAD_N_MAX+1];
        if(atomic_config(Z,config)){
            fprintf(out_stream,"ERROR: SAD guess could not place %d electrons of atom %d\n",Z,i_a+1);
            exit(1);
        }

        //k-th shell of a given l in a minimal basis is the (k+l)-th principal number
        int n_bf=0;
        int n_l[SAD_L_MAX+1]={0,0,0,0};
        int placed=0;
        std::vector<double> n_occ;
        for(int i=i_s; (i<M->read_s.size())&&(shell_center_r[i]==i_a); i++){
            int l=M->read_s[i].contr[0].l;
            if((l>SAD_L_MAX)||(n_l[l]+l>=SAD_N_MAX)){
                fprintf(out_stream,"ERROR: SAD guess got an unexpected shell (l=%d) in %s\n",l,SAD_BASIS);
                exit(1);
            }
            n_l[l]++;
            int e=config[l][n_l[l]+l];
            for(int j=0;j<M->read_s[i].size();j++)
                n_occ.push_back(0.5*double(e)/double(2*l+1));
            n_bf+=M->read_s[i].size();
            placed+=e;
        }
        i_s+=n_l[0]+n_l[1]+n_l[2]+n_l[3];

        if(placed!=Z){
            fprintf(out_stream,"ERROR: SAD guess placed %d of %d electrons of atom %d in %s\n",
                                placed,Z,i_a+1,SAD_BASIS);
            exit(1);
        }

        //atomic block: D = S^-05 n S^-05, so that n are its occupation numbers
        double * S_AA = new double[n_bf*n_bf];
        double * X    = new double[n_bf*n_bf];
        double * nX   = new double[n_bf*n_bf];

        for(int i=0;i<n_bf;i++)
        for(int j=0;j<n_bf;j++)
            S_AA[i*n_bf+j]=S_22[(b_0+i)*n_ro+b_0+j];

        S05_calc(S_AA,X,nX,n_bf);

        for(int i=0;i<n_bf;i++)
        for(int j=0;j<n_bf;j++)
            nX[i*n_bf+j]=n_occ[i]*X[i*n_bf+j];

        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                    n_bf,n_bf,n_bf,1.0,
                    X ,n_bf,
                    nX,n_bf,0.0,
                    DM+b_0*n_ro+b_0,n_ro);

        delete[] S_AA;
        delete[] X   ;
        delete[] nX  ;

        b_0+=n_bf;
    }

    //natural orbitals of the superposition, most occupied first
    HSC_CE(DM, S_22, NO, BUF, occ, n_ro);

    M->guess_occ.resize(n_ro);
    for(int i=0;i<n_ro;i++){
        memcpy(M->MO_VEC_R+i*n_ro, NO+(n_ro-1-i)*n_ro, n_ro*sizeof(double));
        M->guess_occ[i]=occ[n_ro-1-i];
    }

    delete[] S_22;
    delete[] DM  ;
    delete[] NO  ;
    delete[] BUF ;
    delete[] occ ;

    return 1;
}

int sad_guess_fock(molecule * M, double * F){

    int n_ao=M->n_ao;

    //the atomic occupations are fractional, so the density needs its own vectors:
    //the ones the SCF loop starts from carry n_cor of them, filled to 1
    int n_g=0;
    while((n_g<M->guess_occ.size())&&(n_g<n_ao)&&(M->guess_occ[n_g]>0.0))n_g++;

    double * C  = new double[n_g*n_ao];
    double * DM = new double[n_ao*n_ao];
    double * J  = new double[n_ao*n_ao];
    double * K  = new double[n_ao*n_ao];

    for(int i=0;i<n_g;i++)
    for(int j=0;j<n_ao;j++)
        C[i*n_ao+j]=sqrt(M->guess_occ[i])*M->MO_VEC[i*n_ao+j];

    nopt_par_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                   n_ao,n_ao,n_g,1.0,
                   C,n_ao,
                   C,n_ao,0.0,
                   DM,n_ao);

    set_zero_matr(J,n_ao*n_ao);
    set_zero_matr(K,n_ao*n_ao);

    if(RI==0){
        double * J_p[1]={ J};
        double * K_p[1]={ K};
        double * D_p[1]={DM};
        DM_to_F_transform(J_p,K_p,D_p,1,D_p,1,M->s,n_ao);
    }
    else{
        RI_core_realloc(n_g,n_ao);
        DM_to_F_transform_RI(J,K,DM,C,n_g,n_ao);
        RI_core_realloc(M->n_el_calc/2,n_ao);
    }

    for(int i=0;i<n_ao*n_ao;i++)F[i]=M->H_AO[i]+2*J[i]-K[i];

    M->diag_X_AO_in_MO(F);
    M->guess_occ.clear();

    delete[] C ;
    delete[] DM;
    delete[] J ;
    delete[] K ;

    return 0;
}
