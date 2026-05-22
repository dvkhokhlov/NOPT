# ifndef ETC_H
# define ETC_H
# include <vector>

double * new_double_w_check(long long A, const char * stage);

double dist(double x1, double y1, double z1,
            double x2, double y2, double z2,
            int control);

double angle(double x1, double y1, double z1,
             double x2, double y2, double z2,
             double x3, double y3, double z3,
             int control);

double dihedral(double x1, double y1, double z1,
                double x2, double y2, double z2,
                double x3, double y3, double z3,
                double x4, double y4, double z4,
                int control);

// int write_file_name(char * source, int b, char * target, int e_l, int p); 

int compare_int_ar(int * A, int * B, int n);

int D_change(char * line);

int cut_zero(std::vector<double>* c);

class multi_iterator
{
    public:
        int * number;
        int set(int * max_number_inp, int dim_inp);
        int zero();
        int next();
        int not_ended();
        ~multi_iterator();
        
        
    private:
        int dim;
        int * max_number;
        int not_finished;
        int step_for_next(int n);
};


#endif
