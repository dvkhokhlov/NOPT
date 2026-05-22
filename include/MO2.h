#ifndef __MO2
#define __MO2

int VEC_read(char *inp_name, int start_line, double* vec, int n_ao, int n_mo, char ab);

int energy_read(char *inp_name, int start_line, double* en, int n_mo, char ab);

#endif
