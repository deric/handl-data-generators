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
//output file
ofstream * out;
double * points;

void usage() {
    fprintf(stderr, "===== CURE data generator\n");
    fprintf(stderr, "\nThe command line parameters for this generator are:\n\n");
    fprintf(stderr, "$ ./cure [-n <num points>] [-s <seed>]\n\n");
}

int main(int argc, char **argv) {
    int i;
    int type = 1;
    dim = 2;
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
            } else {
                fprintf(stderr, "Unrecognised parameter sequence.\n$ ./cure -h\nfor help.\n");
                exit(1);
            }
        }
    } else {
        fprintf(stderr, "Unrecognised parameter sequence.\n$ ./cure -h\nfor help.\n");
        exit(1);
    }
    cout << "using seed: " << seed << endl;


    char name[50];
    sprintf(name, "cure-t%d-%dn-%dD.dat", type, npoints, dim);
    cerr << "Data written to " << name << endl;
    out = new ofstream(name);
      switch ( type ) {
        case 1:
            gen_data1();
        break;
        default:
            cout<<"Error, type " << type << " is not supported "<< endl;
        break;
    }
    cout << "generated CURE type " << type << " dataset with " << npoints << " data items." << endl;
    out->close();

    exit(0);
}

void gen_data1() {
    int label = 0;
    int s = (int) seed;
    points = uniform_in_circle01_map(dim, npoints, &s);
    // Print data to file
    for (int j=0; j<npoints; j++) {
      for (int k=0; k<dim; k++) {
        *out << points[j*dim+k] << " ";
      }
      *out << label << endl;
    }
}
