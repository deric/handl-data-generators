
#define PI 3.141592653589793

void usage();
void gen_data0();
void gen_data1(int num_small, int num_ellipse, int num_out);
void gen_data2(int num_noise);
void llrand();
double scale(double value, double fromRangeMin, double fromRangeMax, double toRangeMin, double toRangeMax);
//check if coordinates does not collide with existing cluster
bool inside_big(double* vec, double delta);
bool inside_elly(double* vec, double xoffset, double delta);
bool inside_small(double* vec, double yoffset, double delta);
bool on_line(double* vec, double lval, double delta);