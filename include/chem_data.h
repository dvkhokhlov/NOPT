#include <stdio.h>

#ifndef __CHEM_DATA_H 
#define __CHEM_DATA_H

double covalentR(int);
int charge(char * atom, FILE * data);
double mass(char * atom, FILE * data);
int ghost(char * atom, FILE * data);
int atom_name_by_charge(int,char*);
double mass_by_charge(int);

#endif
