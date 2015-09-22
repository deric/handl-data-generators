/**************************

Generate datasets similar to CURE's datasets.

Guha, Sudipto, Rajeev Rastogi, and Kyuseok Shim.
"CURE: an efficient clustering algorithm for
large databases." ACM SIGMOD Record. Vol. 27. No. 2. ACM, 1998.


This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

 */

using namespace std;
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <fstream>
#include <iostream>
#include <iostream>
#include <math.h>
#include "random_data.h"
#include "random.h"
#include "cure.h"


int npoints;
long int seed;
int dim;
double xmin;
double xmax;
//output file
ofstream * out;
double * points;

void usage() {
    fprintf(stderr, "===== CURE data generator\n");
    fprintf(stderr, "Generates square dataset with several geometic shapes\n");
    fprintf(stderr, "by default data are in interval [-1, 1]\n");
    fprintf(stderr, "\nThe command line parameters for this generator are:\n\n");
    fprintf(stderr, "$ ./cure -n <num points> [-s <seed>] [-l <x/y min>] [-m <x/y max>] [-t type of data]\n\n");
}

int main(int argc, char **argv) {
    int i;
    int type = 1;
    dim = 2;
    xmin = -1.0;
    xmax = 1.0;
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
            } else {
                fprintf(stderr, "Unrecognized parameter sequence.\n$ ./cure -h\nfor help.\n");
                exit(1);
            }
        }
    } else {
        fprintf(stderr, "Unrecognized parameter sequence.\n$ ./cure -h\nfor help.\n");
        exit(1);
    }
    cout << "using seed: " << seed << endl;


    char name[50];
    sprintf(name, "cure-t%d-%dn-%dD.dat", type, npoints, dim);
    cerr << "Data written to " << name << endl;
    out = new ofstream(name);
    switch (type) {
        case 1:
            gen_data1();
            break;
        default:
            cout << "Error, type " << type << " is not supported " << endl;
            break;
    }
    cout << "generated CURE type " << type << " dataset with " << npoints << " data items." << endl;
    out->close();

    exit(0);
}

/**
 * Generate 3 circles (a large one and two small ones)
 */
void gen_data1() {
    int label = 0;
    int s = (int) seed;
    double val;


    //use 10% of point to form small circle
    int num_small = (int) npoints * 0.1;
    int num_big = npoints - 2 * num_small;

    //big circle - 80% of data points
    points = uniform_in_circle01_map(dim, num_big, &s);
    double b_max = xmin + (xmax - xmin) * 0.7;
    for (int j = 0; j < num_big; j++) {
        for (int k = 0; k < dim; k++) {
            if (k == 0) {
                val = scale(points[j * dim + k], xmin, xmax, xmin, b_max);
            } else {
                val = points[j * dim + k];
            }

            *out << val << " ";
        }
        *out << label << endl;
    }

    //upper left smaller circle
    label = 1;
    //double gap = (xmax - xmin) * 0.15;
    points = uniform_in_circle01_map(dim, num_small, &s);
    double quad = (xmax - xmin) * 0.12;
    double half = xmin + (xmax - xmin) / 2;
    for (int j = 0; j < num_small; j++) {
        for (int k = 0; k < dim; k++) {
            if (k == 0) { //x dimension
                val = scale(points[j * dim + k], xmin, xmax, xmax - quad, xmax);
            } else {
                val = scale(points[j * dim + k], xmin, xmax, half + quad, half + 2 * quad);
            }

            *out << val << " ";
        }
        *out << label << endl;
    }

    //bottom left smaller circle
    label = 2;
    points = uniform_in_circle01_map(dim, num_small, &s);
    for (int j = 0; j < num_small; j++) {
        for (int k = 0; k < dim; k++) {
            if (k == 1) { //y dimension
                val = scale(points[j * dim + k], xmin, xmax, half - 2 * quad, half - quad);
            } else {
                val = scale(points[j * dim + k], xmin, xmax, xmax - quad, xmax);
            }

            *out << val << " ";
        }
        *out << label << endl;
    }
}

double scale(double value, double fromRangeMin, double fromRangeMax, double toRangeMin, double toRangeMax) {
    return ((value - fromRangeMin) * (toRangeMax - toRangeMin) / (fromRangeMax - fromRangeMin) + toRangeMin);
}
