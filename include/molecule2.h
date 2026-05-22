
#ifndef __molecule2
#define __molecule2
#include "molecule.h"
#include "doCI_data.h"


int calc_2el_CI(molecule* cD, doCI_data * D);

int calc_2el_CI_1st_order(molecule* cD, doCI_data * D, int nshA, int nshB);

#endif
