#ifndef GEOM_H
#define GEOM_H
# include "molecule.h"
# include "inp_out.h"

int orbital_orthogonalize(molecule * A);

int shift_molecule(molecule * A, double * dir, double l);


class coord_list;

class coord_list{
    
    public:
        std::vector<double> value;
        std::vector<double> value_delta;
        std::vector<std::vector<int>> refs;
        std::vector<std::vector<char>> name;
        std::vector<int> type;//0 - dist, 1 - angle
        
        
        coord_list();
        ~coord_list();
    
};

int read_coord(char * l, double * v, int p, int t);

int read_cl(recursive_file * in);

int dist_scale(double c);

int coord_sync(double * v);

int print_coord_list();

int fprint_coord_list(const char * name);

int set_std_delta();

#endif
