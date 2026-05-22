#ifndef MOL_LINK_H
#define MOL_LINK_H
//user
# include "molecule.h"
int STATES_link(molecule * AB, molecule * A, molecule * B, int nAB, int nA, int nB);

int MO_link(double * AB, double * A, double * B, int dimA, int dimB, int corA, int corB, int actA, int actB);

int MO_link_2(double * AB, double * A, double * B, int dimA, int dimB, int corA, int corB, int actA, int actB);

int geom_link(molecule * AB, molecule * A, molecule * B);

int mol_link(molecule * AB,molecule * A, molecule * B, int gen_states);

int mol_array_link(molecule * O, molecule * IN, int n, int gen_states);

int mol_cpy(molecule *  O, molecule * IN , int gen_states);
#endif
