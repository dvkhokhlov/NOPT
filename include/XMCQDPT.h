#ifndef NOPT_H
#define NOPT_H

int three_block_diagonalization(double * M, double * V, double * e1,
                                int d1, int d2, int d3, int d_all,
                                double * O1, double * O2, double* O3, double * O4);

int QDPT2(molecule * A, xmc_par * xmc, char * job_name);
//xmc =0 for MCQDPT; xmc =1 for XMCQDPT


#endif
