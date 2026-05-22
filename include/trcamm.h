#ifndef __trcamm
#define __trcamm

#include <vector>
#include "libint2/shell.h"
using libint2::Shell;


int gen_atomic_DM(double * A_DM, double * DM, int atom, std::vector<Shell> * s, std::vector<int> * shell_center, int n_ao);


// int gen_atomic_DM_anti(double * A_DM, double * DM, int atom, int * ao_center, int n_ao);

// int calc_TrCAMM_by_DM(molecule * A, double ** DM, double * S_WF, int n_s);
// 
// int calc_TrCAMM_by_DM_2(molecule * A, double ** DM, double * S_WF, int n_s1, int n_s2);


// inline double calc_atom_multipole(double * A_DM, double * DM, double * PM, int atom, int * ao_center, int n_ao, int anti){
//     
//     double m;
//     
//     if(anti==0)gen_atomic_DM     (A_DM, PM, atom, ao_center, n_ao);
//     if(anti==1)gen_atomic_DM_anti(A_DM, PM, atom, ao_center, n_ao);
//     m= E_1el_calc(A_DM, DM, n_ao, n_ao);
// //     m=-E_1el_calc(A_DM, DM, n_ao, n_ao);
//     
//     return m;
// }

#endif




 
