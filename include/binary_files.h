#ifndef _XDR_LINK_H
#define _XDR_LINK_H

int get_full_name(char ** f, char * p, const char * n);

int binary_array_read(double * M, const char * name, char * prefix, int n);

int binary_array_write(double * M, const char * name, char * prefix, int n, int ask);

int binary_VEC_write(const char * out_name, double * MO_VEC, double * orb_energy,int n_mo, int n_ao, double * MO_VEC_B, double * orb_energy_B);

int binary_VEC_read(const char * inp_name, double * MO_VEC, double * orb_energy, int n_mo, int n_ao, double * MO_VEC_B, double * orb_energy_B);

#endif
