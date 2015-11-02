

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
#include "disk.h"


int npoints;
long int seed;
int dim;
//output file
ofstream * out;

double xmin;
double xmax;

void usage() {
    fprintf(stderr, "===== disk-in-disk data generator\n");
    fprintf(stderr, "Generates two circles, one inside another\n");
    fprintf(stderr, "by default data are in interval [-10, 10]\n");
    fprintf(stderr, "\nThe command line parameters for this generator are:\n\n");
    fprintf(stderr, "$ ./disk -n <num points> [-s <seed>] [-d <dimension>] [-l <x/y min>] [-m <x/y max>] [-t type of data] [-r inner radius]\n\n");
}

int main(int argc, char **argv) {
    int i;
    int type = 0;
    dim = 2;
    xmin = -10.0;
    xmax = 10.0;
    double r = 0.5;
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
        } else {
            fprintf(stderr, "Unrecognized parameter sequence.\n$ ./disk -h\nfor help.\n");
            exit(1);
        }
      }
    } else {
        fprintf(stderr, "Unrecognized parameter sequence.\n$ ./disk -h\nfor help.\n");
        exit(1);
    }
    cout << "using seed: " << seed << endl;


    char name[50];
    sprintf(name, "disk-t%d-%dn.dat", type, npoints);
    cerr << "Data written to " << name << endl;
    out = new ofstream(name);
    int num_inner;
    int num_outer;
    switch (type) {
        case 0:
            num_inner = (int) npoints * 0.5;
            num_outer = npoints - num_inner;
            gen_data0(num_inner, num_outer, r);
            break;
        default:
            cout << "Error, type " << type << " is not supported " << endl;
            break;
    }
    cout << "generated disk-in-disk type " << type << " dataset with " << npoints << " data items." << endl;
    out->close();

    exit(0);
}

/**
 * Generate disk-in-disk
 * r - inner radius of an annulus (outer radius is 1.0)
 */
void gen_data0(int num_inner, int num_outer, double r) {
  int label = 0;
  double * points;
  double val;
  int s = (int) seed;

  points = uniform_in_annulus01_accept(dim, num_outer, r, &s);
  for (int j = 0; j < num_outer; j++) {
    for (int k = 0; k < dim; k++) {
      val = scale(points[j * dim + k], -1.0, 1.0, xmin, xmax);
      *out << val << " ";
    }
    *out << label << endl;
  }
  delete [] points;
}

double scale(double value, double fromRangeMin, double fromRangeMax, double toRangeMin, double toRangeMax) {
  return ((value - fromRangeMin) * (toRangeMax - toRangeMin) / (fromRangeMax - fromRangeMin) + toRangeMin);
}
