#ifndef __bas_lib_read
#define __bas_lib_read

#include "molecule.h"

std::vector<Shell> basis_lib_read_gbs(molecule * M,  const char * lib_name,
                                      int in_lib, int start, 
                                      std::vector<std::vector<double>> * lib_coef, 
                                      std::vector<int> * shell_center, bool pure,
                                      std::vector<double> * energy,
                                      std::vector<int> * is_core);

std::vector<Shell> basis_lib_read_exp(molecule * M,  const char * lib_name, 
                                      int in_lib, int start, 
                                      std::vector<std::vector<double>> * lib_coef, 
                                      std::vector<int> * shell_center, bool pure);


int ecp_lib_read_gbs (molecule * M, int in_lib, int start, const char * lib_name, const char * suffix);

int ecp_lib_read_exp (molecule * M, int in_lib, int start, const char * lib_name, const char * suffix);

int ecp_lib_read_grpp(molecule * M, int in_lib, int start, const char * lib_name, const char * suffix);

#endif
