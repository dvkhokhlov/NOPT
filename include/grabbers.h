extern int grab_E;
extern std::vector<double> E_grabbed;

extern int grab_D;
extern int grab_HS_data;
extern int grabber_n_states;
extern double *  H_grabbed;
extern double *  S_grabbed;
extern double * Dx_grabbed;
extern double * Dy_grabbed;
extern double * Dz_grabbed;

extern FILE * e_scan;
extern FILE * d_scan;

int alloc_HS_grabber(int n, int m);
int alloc_D_grabber(int n, int m);
