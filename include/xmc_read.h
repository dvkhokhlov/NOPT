////ALL XMC reading files once will be here
#ifndef XMC_RAED_H
#define XMC_RAED_H

std::vector<double> xmc_ifitd_read(char * xmc_name);

int H0PP_read(double * Hpp_XMC, int n_s, char * xmc_name);

int read_H02_with_H0PP(double * H2, double * H0PP, int n_s, char * xmc_name);

#endif
