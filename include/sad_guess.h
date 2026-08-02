#ifndef __sad_guess
#define __sad_guess

#include "molecule.h"

// Superposition of atomic densities in a minimal basis. Fills read_s/n_ro/MO_VEC_R
// with the natural orbitals of the superposition, sorted by occupation, and leaves
// the projection to the calculation basis to the caller.
// Returns 0 - the molecule is out of the guess coverage, nothing was changed.
int sad_guess(molecule * M);

// Fock of the atomic density and the orbitals it defines; consumes guess_occ.
int sad_guess_fock(molecule * M, double * F);

#endif
