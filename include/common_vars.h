#include "inp_par_read.h"

extern int D5;
extern int RI;
extern int RI_flag;
extern int SO;
extern int nQM;

extern int engines_initialized;
extern int num_threads;
extern int IS_SYM;
extern int LINEAR;
extern int write_orbs;
extern int write_ci;
extern int print_dipole;
extern int print_dispersion;
extern int print_quadrupole;
extern int print_mulliken;

extern FILE * out_stream;
int set_out_stream_by_name(const char * name,const char *mode);
#define nopt_printf(...) fprintf(out_stream, __VA_ARGS__)
int set_out_stream_by_FILE(FILE * file, const char * name);
int set_out_stream_stdout();
int set_out_stream_stderr();

extern double * RI_AO_M;

extern backup_par backup;

extern int * nopt_parallel_index;

int nopt_parallel_init(int n);

int nopt_parallel_finalize();
