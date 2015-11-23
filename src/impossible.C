
using namespace std;
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <iostream>
#include "random_data.h"
#include "random.h"
#include "impossible.h"

/* Random number stuff */
#define RN ran0(&seed)
#define GN gaussian(&seed)

/* Random number generator */
/* Copyright Numerical Recipes in C */
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-14
#define RNMX (1.0-EPS)
#define MASK 123459876

double gaussian(long *idum);

int npoints;
long int seed;
int dim;
//output file
ofstream * out;

double xmin;
double xmax;

void usage() {
    fprintf(stderr, "===== impossible data generator\n");
    fprintf(stderr, "Generates set of data that is according to Jain 2010 impossible to cluster\n");
    fprintf(stderr, "by default data are in interval [-10, 10]\n");
    fprintf(stderr, "\nThe command line parameters for this generator are:\n\n");
    fprintf(stderr, "$ ./impossible -n <num points> [-s <seed>] [-d <dimension>] [-l <x/y min>] [-m <x/y max>] [-b ratio noise]\n\n");
}

int main(int argc, char **argv) {
    int i;
    int type = 0;
    dim = 2;
    xmin = -10.0;
    xmax = 10.0;
    double r = 0.5;
    double b = 0.1;
    double gap = 0.1;
    //std::random_device rd;
    //std::mt19937 gen(rd());
    //std::uniform_int_distribution<unsigned long long> dis(lowerBorder, upperBorder);
    //seed = dis(gen)

    srand(time(NULL));

    seed = rand();

    if (argc <= 2) {
        usage();
        exit(1);
    }

    if ((argc > 2)&&(argc % 2 == 1)) {
      for (i = 1; i < argc - 1; i += 2) {
        if (strcmp("-n", argv[i]) == 0) {
            npoints = atoi(argv[i + 1]);
        } else if (strcmp("-s", argv[i]) == 0) {
            seed = atol(argv[i + 1]);
        } else if (strcmp("-t", argv[i]) == 0) {
            type = atoi(argv[i + 1]);
        } else if (strcmp("-d", argv[i]) == 0) {
            dim = atoi(argv[i + 1]);
        } else if (strcmp("-l", argv[i]) == 0) {
            xmin = atof(argv[i + 1]);
        } else if (strcmp("-m", argv[i]) == 0) {
            xmax = atof(argv[i + 1]);
        } else if (strcmp("-r", argv[i]) == 0) {
            r = atof(argv[i + 1]);
        } else if (strcmp("-b", argv[i]) == 0) {
            b = atof(argv[i + 1]);
        } else if (strcmp("-g", argv[i]) == 0) {
            gap = atof(argv[i + 1]);
        } else {
            fprintf(stderr, "Unrecognized parameter sequence.\n$ ./impossible -h\nfor help.\n");
            exit(1);
        }
      }
    } else {
        fprintf(stderr, "Unrecognized parameter sequence.\n$ ./impossible -h\nfor help.\n");
        exit(1);
    }
    cout << "using seed: " << seed << endl;


    char name[50];
    sprintf(name, "impossible-%dn.dat", npoints);
    cerr << "Data written to " << name << endl;
    out = new ofstream(name);
    int num_quad;
    int num_noise;
    switch (type) {
        case 0:
            num_quad = (int) ((npoints - npoints * b) / 4);
            num_noise = npoints - 4 * num_quad;
            gen_data(num_quad, num_noise, r, gap);
            break;
        default:
            cout << "Error, type " << type << " is not supported " << endl;
            break;
    }
    cout << "generated impossible type " << type << " dataset with " << npoints << " data items." << endl;
    out->close();

    exit(0);
}

/**
 * Generate impossible set
 */
void gen_data(int num_quad, int num_noise, double r, double gap) {
  int label = 0;
  double * points;
  double val;
  int s = (int) seed;
  double* half = new double[dim];
  for(int i = 0; i < dim; i++){
     half[i] = (xmin + xmax) / 2;
  }
  //top-left gaussian circle
  points = normal_circular(dim, num_quad, &s);
  //find min,max for generated data
  double* gmin = new double[dim];
  double* gmax = new double[dim];
  for(int j = 0; j < num_quad; j++){
    for (int k = 0; k < dim; k++) {
      val = points[j * dim + k];
      if(val < gmin[k]){
        gmin[k] = val;
      }
      if(val > gmax[k]){
        gmax[k] = val;
      }
    }
  }

  //write scaled gaussian to upper left corner
  for (int j = 0; j < num_quad; j++) {
    for (int k = 0; k < dim; k++) {
      if(k == 1){
        val = scale(points[j * dim + k], gmin[k], gmax[k], half[k], xmax);
      }else {
        val = scale(points[j * dim + k], gmin[k], gmax[k], xmin, half[k]);
      }
      *out << val << " ";
    }
    *out << label << endl;
  }
  delete[] points;

  //3 circles inside each other
  int num_circles = 2*num_quad - num_quad / 2;
  int c1 = (int) num_circles * 0.5;
  int c2 = (int) num_circles * 0.3;
  double cr = (half[0] - xmin) / 2;
  double cx = (half[0] - xmin) / 2;
  double cy = (half[1] - xmax) / 2;

  draw_circle(c1, cx, cy, half, 1, 0.8 * cr);
  draw_circle(c2, cx, cy, half, 2, 0.5 * cr);
  draw_circle(num_circles - c1 - c2, cx, cy, half, 3, 0.2 * cr);

  //spirals
  cx =  (half[0] + xmax) / 2;
  cy = (half[1] + xmax) / 2;
  int num_spiral = num_quad / 2;

  // spirals that are very close
  //draw_spiral(num_spiral, cx, cy, 4, 0.9 * cr, 0.0, 0.2, 0.07, 1, 1);
  //draw_spiral(num_spiral, cx-0.35, cy-0.35, 5, 0.9 * cr, -1.9, 0.9, 0.0698, 0, -1);

  draw_spiral(num_spiral, cx, cy, 4, 0.8 * cr, 0.0, 0.2, 0.03, 1, 1);
  draw_spiral(num_spiral, cx-0.5, cy-0.6, 5, 0.8 * cr, -0.5, 0.4, 0.025, 0, -1);
  //draw_spiral(num_spiral, cx - 0.3, cy - 0.4, 5, 0.8 * cr, 0.2, 0.5, -0.07, 0);

  cx =  (xmin + half[0]) / 2;
  cy = (xmin + half[1]) / 2;
 // draw_circle(num_quad, cx, cy, half, 6, cr);
  //move half of poits to "circles cluster"
  draw_elly(num_quad / 2, cx, cy, 6, cr, xmin, half);

  delete[] half;
}

void draw_elly(int length, double cx, double cy, int label, double cr, double xmin, double* half){
  int i, j;
  double sum;
  double *cart;
  double *gauss;
  double *origin;
  double *sph;
  cart = (double *) malloc(dim * sizeof (double));
  gauss = (double *) malloc(dim * sizeof (double));
  origin = (double *) malloc(dim * sizeof (double));
  sph = (double *) malloc(dim * sizeof (double));
  double* foc = (double *) malloc(dim * sizeof (double));
  // 1 - sumdist, where sumdist is the maximum allowed sum of the two distances between a point and the two foci. (This relates to the eccentricity of the cluster, hence the name)
  double ecc = 0.05 + (RN * 0.1);
  double* points = new double[length * dim];
  double* gmin = new double[dim];
  double* gmax = new double[dim];

  // first generate the second on a unit hypersphere around the origin to get a random orientation
  sum = 0.0;
  for (j = 0; j < dim; j++) {
    origin[j] = 0.0;
    gauss[j] = GN;
    if (RN < 0.5){
      gauss[j] *= -1;
    }
    sum += gauss[j] * gauss[j];
  }
  for (j = 0; j < dim; j++){
    foc[j] = (1.0 / sqrt(sum)) * gauss[j] * length;
  }


  // then, generate the points a Gaussian distributed distance in a random direction away from a
  // uniformly random point somewhere along the major axis of the ellipsoid between the foci. Accept
  // the point if it is within the ellipsoid.
  i = 0;
  bool success = true;
  while (i < length) {
      double ax = RN;
      double rad = GN;
      sum = 0.0;

      if (success) {
          for (j = 0; j < dim; j++) {
              gauss[j] = GN;
              sum += gauss[j] * gauss[j];
          }
          for (j = 0; j < dim; j++) {
              sph[j] = (1.0 / sqrt(sum)) * gauss[j];
              //  printf("%g ", sph[j]);
          }
      }
      // printf("\n");

      for (j = 0; j < dim; j++) {
          cart[j] = ax * foc[j] + sph[j] * rad * ecc * length * 2; // the Cartesian coordinates of the point
          // printf("%g ", cart[j]);
      }

      double diff;

      // is the point inside the ellipsoid ?
      if ((diff = ((Eucdist(cart, origin) + Eucdist(cart, foc)) - (1.0 + ecc) * length)) < 0) {
          //  if(RN< diff/(ecc*(-length)))
          //    {
          for (j = 0; j < dim; j++) {
              // printf("%g ", cart[j]);
              //datitem[i][j] = cart[j];
              //*out << cart[j] << " ";
              points[j + i * dim] = cart[j];
              if(cart[j] < gmin[j]){
                gmin[j] = cart[j];
              }
              if(cart[j] > gmax[j]){
                gmax[j] = cart[j];
              }
          }
          //*out << label << endl;
          //  printf("\n");
          i++;
          success = true;
          //  }
      } else {
          success = false;
      }
  }

  double val;
  for (int j = 0; j < length; j++) {
    for (int k = 0; k < dim; k++) {
      val = scale(points[j * dim + k], gmin[k], gmax[k], xmin, 0.8*half[k]);
      *out << val << " ";
    }
    *out << label << endl;
  }
  delete[] points;
}

void draw_spiral(int num_spiral, double cx, double cy, int label, double cr, double a, double b, double c, int p, int q){
  double angle;
  double smin = 0;
  double smax = 0;
  double* points = new double[num_spiral * dim];
  double val;
  for(int i = 0; i < num_spiral; i++){
    for (int k = 0; k < dim; k++) {
      angle = c * i *q;
      if(k == p){
        val = (a + b * angle)*sin(angle)/3;
        //val = scale(val, -230, 230, cy - cr, cy + cr);
      }else {
        val = (a + b * angle)*cos(angle)/3;
        //val = scale(val, -230, 230, cx - cr, cx + cr);
      }
      if(val < smin){
        smin = val;
      }
      if(val > smax){
        smax = val;
      }
      points[k + i * dim] = val;
    }
  }
  cout << "min = "<< smin << ", max = " << smax << endl;

  for (int j = 0; j < num_spiral; j++) {
    for (int k = 0; k < dim; k++) {
      if(k == 1){
        val = scale(points[j * dim + k], smin, smax, cy - cr, cy + cr);
      }else {
        val = scale(points[j * dim + k], smin, smax, cx - cr, cx + cr);
      }
      *out << val << " ";
    }
    *out << label << endl;
  }
  delete [] points;
}


void draw_circle(int num_pts, double cx, double cy, double* half, int label, double r){
  double * points;
  double val;
  int s = (int) seed;
  points = uniform_in_annulus01_accept(dim, num_pts, 0.9, &s);
  for (int j = 0; j < num_pts; j++) {
    for (int k = 0; k < dim; k++) {
      if(k == 1){
        val = scale(points[j * dim + k], -1.0, 1.0, cy - r, cy + r);
      }else {
        val = scale(points[j * dim + k], -1.0, 1.0, cx - r, cx + r);
      }
      *out << val << " ";
    }
    *out << label << endl;
  }
  delete [] points;
}

double scale(double value, double fromRangeMin, double fromRangeMax, double toRangeMin, double toRangeMax) {
  return ((value - fromRangeMin) * (toRangeMax - toRangeMin) / (fromRangeMax - fromRangeMin) + toRangeMin);
}

double Eucdist(double *a, double *b) {
    double sum = 0.0;
    for (int i = 0; i < dim; i++) {
        sum += (a[i] - b[i])*(a[i] - b[i]);
    }
    return (sqrt(sum));
}

/* Generate a N(0,1) r.v. */
double gaussian(long *idum) {
    static int iset = 0;
    static double gset;
    double fac, r, v1, v2;
    if (iset == 0) {
        do {
            v1 = 2.0 * RN - 1.0;
            v2 = 2.0 * RN - 1.0;
            r = v1 * v1 + v2*v2;
        } while (r >= 1.0 || r == 0.0);
        fac = sqrt(-2.0 * log(r) / r);
        gset = v1*fac;
        iset = 1;
        return v2*fac;
    } else {
        iset = 0;
        return gset;
    }
}
