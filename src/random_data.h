// random_data.h

double *brownian ( int m, int n, int *seed );
double d_epsilon ( void );
double d_max ( double x, double y );
double d_min ( double x, double y );
int d_nint ( double x );
double d_normal_01 ( int *seed );
double d_pi ( void );
double d_random ( double rlo, double rhi, int *seed );
double d_uniform_01 ( int *seed );
double *dge_mxv ( int m, int n, double a[], double x[] );
void direction_random_nd ( int m, int *seed, double w[] );
void dmat_print ( int m, int n, double a[], char *title );
double *dpo_fa ( int n, double a[] );
double *dut_mxv ( int m, int n, double a[], double x[] );
void dvec_normal_01 ( int n, int *seed, double x[] );
void dvec_print ( int n, double a[], char *title );
void dvec_uniform_01 ( int n, int *seed, double r[] );
unsigned long get_seed ( void );
double *grid_in_cube01 ( int m, int n, int center, int *seed );
int grid_side ( int m, int n );
bool halton_base_check ( int ndim, int base[] );
double *halton_in_circle01_accept ( int m, int n, int *seed );
double *halton_in_circle01_map ( int m, int n, int *seed );
double *halton_in_cube01 ( int m, int n, int *seed );
bool halton_leap_check ( int ndim, int leap[] );
bool halton_n_check ( int n );
bool halton_ndim_check ( int ndim );
bool halton_seed_check ( int ndim, int seed[] );
bool halton_step_check ( int step );
bool hammersley_base_check ( int ndim, int base[] );
double *hammersley_in_cube01 ( int m, int n, int *seed );
bool hammersley_leap_check ( int ndim, int leap[] );
bool hammersley_n_check ( int n );
bool hammersley_ndim_check ( int ndim );
bool hammersley_seed_check ( int ndim, int seed[] );
bool hammersley_step_check ( int step );
int i_factorial ( int n );
int i_max ( int i1, int i2 );
int i_min ( int i1, int i2 );
int i_modp ( int i, int j );
int i_random ( int ilo, int ihi, int *seed );
void i_to_halton ( int ndim, int step, int seed[], int leap[], int base[], 
  double r[] );
void i_to_halton_sequence ( int ndim, int n, int step, int seed[], int leap[],
  int base[], double r[] );
void i_to_hammersley_sequence ( int ndim, int n, int step, int seed[], int leap[],
  int base[], double r[] );
void ivec_transpose_print ( int n, int a[], char *title );
void ksub_random2 ( int n, int k, int *seed, int a[] );
double *polygon_centroid_2d ( int n, double v[] );
int prime ( int n );
unsigned long random_initialize ( unsigned long seed );
int s_len_trim ( char *s );
void scale_from_simplex01 ( int m, int n, double t[], double x[] );
void scale_to_ball01 ( int m, int n, double x[] );
void scale_to_block01 ( int m, int n, double x[] );
void scale_to_cube01 ( int m, int n, double x[] );
void timestamp ( void );
char *timestring ( void );
double triangle_area_2d ( double v1[2], double v2[2], double v3[2] );
void tuple_next_fast ( int m, int n, int rank, int x[] );
double *normal ( int m, int n, double r[], double mu[], int *seed );
double *normal_circular ( int m, int n, int *seed );
double *normal_multivariate ( int m, int n, double r[], double mu[], 
  int *seed );
double *normal_simple ( int m, int n, int *seed );
double *uniform_in_annulus01_accept ( int m, int n, double r, int *seed );
double *uniform_in_circle01_map ( int m, int n, int *seed );
double *uniform_in_cube01 ( int m, int n, int *seed );
double *uniform_in_ellipsoid_map ( int m, int n, double a[], double r, int *seed );
double *uniform_in_parallelogram_map ( double v1[2], double v2[2], 
  double v3[2], int n, int *seed );
double *uniform_in_polygon_map ( int nv, double v[], int n, int *seed );
double *uniform_in_sector_map ( double r1, double r2, double t1, 
  double t2, int n, int *seed );
double *uniform_in_simplex01_map ( int m, int n, int *seed );
double *uniform_in_sphere01_map ( int m, int n, int *seed );
double *uniform_in_triangle_map1 ( double v1[2], double v2[2], double v3[2],
  int n, int *seed );
double *uniform_in_triangle_map2 ( double v1[2], double v2[2], double v3[2],
  int n, int *seed );
double *uniform_in_triangle01_map ( int n, int *seed );
double *uniform_on_ellipsoid_map ( int m, int n, double a[], 
  double r, int *seed );
double *uniform_on_simplex01_map ( int m, int n, int *seed );
double *uniform_on_sphere01_map ( int m, int n, int *seed );
double *uniform_walk ( int m, int n, int *seed );
void write_data ( int ndim, int n, double r[], char *file_out_name, 
  char *title );
