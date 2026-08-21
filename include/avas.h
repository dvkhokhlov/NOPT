#ifndef __avas
#define __avas

#include "molecule.h"
#include "inp_par_read.h"

// Atomic Valence Active Space. Rotates the occupied and the virtual row block of MO_VEC so
// that the orbitals overlapping the requested atomic reference shells fill the active window
// [n_cor_orb, n_cor_orb+n_act_orb); the window sizes come from $ACT_SPACE and are not changed.
// C1 and all-electron only. Placement is by construction - the orbitals must not be re-sorted.
int avas_steer(molecule * M, const avas_par & avas, char * job_name);

#endif
