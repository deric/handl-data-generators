
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

  for (int k = 0; k < dim; k++) {
    cout << "gaussian dim "<< k <<": [" << gmin[k] << ", " << gmax[k] << "] " << endl;
    if(k == 1){
      cout << "scaling to: " << half[k] << ", " << xmax << endl;
    }else {
      cout << "scaling to: " << xmin << ", " << half[k] << endl;
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

  label = 1;
  double amin = xmin * 0.8;
  double amax = xmax * 0.8;

  //3 circles inside each other
  points = uniform_in_annulus01_accept(dim, num_quad, 0.9, &s);
  for (int j = 0; j < num_quad; j++) {
    for (int k = 0; k < dim; k++) {
      if(k == 1){
        val = scale(points[j * dim + k], -1.0, 1.0, half[k], amax);
      }else {
        val = scale(points[j * dim + k], -1.0, 1.0, amin, half[k]);
      }
      *out << val << " ";
    }
    *out << label << endl;
  }
  delete [] points;




  delete[] half;
}

double scale(double value, double fromRangeMin, double fromRangeMax, double toRangeMin, double toRangeMax) {
  return ((value - fromRangeMin) * (toRangeMax - toRangeMin) / (fromRangeMax - fromRangeMin) + toRangeMin);
}