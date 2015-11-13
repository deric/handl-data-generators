#define PI 3.141592653589793

void usage();
void gen_data(int num_quad, int num_noise, double r, double gap);
double scale(double value, double fromRangeMin, double fromRangeMax, double toRangeMin, double toRangeMax);

void draw_circle(int num_pts, double amin, double amax, double* half, int label, double cr);
void draw_spiral(int num_pts, double cx, double cy, int label, double cr, double a, double b, double c, int p, int q );
