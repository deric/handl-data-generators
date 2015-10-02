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

//===============
// derived values
double xmin;
double xmax;
//uppper bounds of big circle
double b_max;
double quad;
double half;
double gap;
double e_max;

void usage() {
    fprintf(stderr, "===== CURE data generator\n");
    fprintf(stderr, "Generates square dataset with several geometic shapes\n");
    fprintf(stderr, "by default data are in interval [-1, 1]\n");
    fprintf(stderr, "\nThe command line parameters for this generator are:\n\n");
    fprintf(stderr, "$ ./cure -n <num points> [-s <seed>] [-d <dimension>] [-l <x/y min>] [-m <x/y max>] [-t type of data]\n\n");
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
    int num_small;
    int num_ellipse;
    int num_out;
    int num_noise;
    switch (type) {
        case 0:
            gen_data0();
            break;
        case 1:
            //use 10% of point to form small circle
            num_small = (int) npoints * 0.1;
            num_ellipse = (int) npoints * 0.15;
            num_out = (int) npoints * 0.01;
            //big circle will contain 49% of data points
            gen_data1(num_small, num_ellipse, num_out);
            break;
        case 2:
            //use 10% of point to form small circle
            num_small = (int) npoints * 0.1;
            num_ellipse = (int) npoints * 0.15;
            // 5% noise
            num_noise = npoints * 0.05;
            num_out = (int) (npoints * 0.01 + num_noise);
            //big circle will contain 44% of data points
            gen_data1(num_small, num_ellipse, num_out);
            gen_data2(num_noise);
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
 * Generate 3 circular clusters
 */
void gen_data0() {
    int label = 0;
    int s = (int) seed;
    double val;

    //use 10% of point to form small circle
    int num_small = (int) npoints * 0.1;
    int num_big = npoints - 2 * num_small;

    //big circle - 80% of data points
    points = uniform_in_circle01_map(dim, num_big, &s);
    b_max = xmin + (xmax - xmin) * 0.7;
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
    delete [] points;

    //upper right smaller circle
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
    delete [] points;

    //bottom right smaller circle
    label = 2;
    points = uniform_in_circle01_map(dim, num_small, &s);
    for (int j = 0; j < num_small; j++) {
        for (int k = 0; k < dim; k++) {
            if (k == 0) { //x dimension
                val = scale(points[j * dim + k], xmin, xmax, xmax - quad, xmax);
            } else {
                val = scale(points[j * dim + k], xmin, xmax, half - 2 * quad, half - quad);
            }

            *out << val << " ";
        }
        *out << label << endl;
    }
    delete [] points;
}

/**
 * Generate 5 high density clusters (one big circle, 2 small circles and two ellipsoids)
 */
void gen_data1(int num_small, int num_ellipse, int num_out) {
    int label = 0;
    int s = (int) seed;
    double val;

    int num_big = npoints - 2 * num_small - 2 * num_ellipse - num_out;

    points = uniform_in_circle01_map(dim, num_big, &s);
    b_max = xmin + (xmax - xmin) * 0.7;
    for (int j = 0; j < num_big; j++) {
        for (int k = 0; k < dim; k++) {
            val = scale(points[j * dim + k], xmin, xmax, xmin, b_max);
            *out << val << " ";
        }
        *out << label << endl;
    }
    delete [] points;

    //upper right smaller circle
    label = 1;
    points = uniform_in_circle01_map(dim, num_small, &s);
    quad = (b_max - xmin) * 0.12;
    half = xmin + (b_max - xmin) / 2;
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
    delete [] points;

    //bottom right smaller circle
    label = 2;
    points = uniform_in_circle01_map(dim, num_small, &s);
    for (int j = 0; j < num_small; j++) {
        for (int k = 0; k < dim; k++) {
            if (k == 0) { //x dimension
                val = scale(points[j * dim + k], xmin, xmax, xmax - quad, xmax);
            } else {
                val = scale(points[j * dim + k], xmin, xmax, half - 2 * quad, half - quad);
            }

            *out << val << " ";
        }
        *out << label << endl;
    }
    delete [] points;

    //ellipsoid
    gap = (xmax - xmin) * 0.05;
    half = xmin + (xmax - xmin) / 2;
    e_max = (xmax - xmin) * 0.2;
    label = 3;
    points = uniform_in_circle01_map(dim, num_ellipse, &s);
    for (int j = 0; j < num_ellipse; j++) {
        for (int k = 0; k < dim; k++) {
            if (k == 0) { //x dimension
                val = scale(points[j * dim + k], xmin, xmax, xmin, half - gap);
            } else {
                val = scale(points[j * dim + k], xmin, xmax, xmax - e_max, xmax);
            }

            *out << val << " ";
        }
        *out << label << endl;
    }
    delete [] points;

    //second ellipsoid
    label = 4;
    points = uniform_in_circle01_map(dim, num_ellipse, &s);
    for (int j = 0; j < num_ellipse; j++) {
        for (int k = 0; k < dim; k++) {
            if (k == 0) { //x dimension
                val = scale(points[j * dim + k], xmin, xmax, half + gap, xmax);
            } else {
                val = scale(points[j * dim + k], xmin, xmax, xmax - e_max, xmax);
            }

            *out << val << " ";
        }
        *out << label << endl;
    }
    delete [] points;

    //outliers connecting both ellipsoids
        //second ellipsoid
    label = 5;
    for (int j = 0; j < num_out; j++) {
        for (int k = 0; k < dim; k++) {
            if (k == 0) { //x dimension
                //val = d_uniform_01(&s);
                val = d_random(half - gap, half + gap, &s);
                //val = scale(val, xmin, xmax, half - gap, half + gap);
            } else {
                val = xmax - e_max / 2;
            }
            *out << val << " ";
        }
        *out << label << endl;
    }
}

/**
* Add sixth cluster to Type 1 dataset containing only noise
*/
void gen_data2(int num_noise) {
    int label = 6;
    int s = (int) seed;
    double* vec;
    int i = 0;
    vec = new double[dim];
    double delta, e1xoffset, e2xoffset, abshalf, a, sm1_off, sm2_off, hf;
    double sm_delta = (xmax - xmin) * 0.005;
    //neighbourhood around cluster without points
    delta = (xmax - xmin) * 0.02;

    abshalf = (xmax - xmin) / 2.0;
    a = (abshalf - gap) / 2.0;
    e1xoffset = half - gap - a;
    e2xoffset = half + gap + a;
    hf = xmin + (b_max - xmin) / 2;
    sm1_off = hf + quad;
    sm2_off = hf - 2 * quad;
    double lval = xmax - e_max / 2;
    while(i < num_noise){
        dvec_uniform_01(dim, &s, vec);
        //scale the vector
        for (int k = 0; k < dim; k++) {
            vec[k] = scale(vec[k], 0.0, 1.0, xmin, xmax);
        }

        if(inside_big(vec, delta)){
            //skip the point
        }else if (inside_elly(vec, e1xoffset, delta)){
            //skip the point
        }else if (inside_elly(vec, e2xoffset, delta)){
            //skip the point
        }else if (inside_small(vec, sm1_off, sm_delta)){
            //skip the point
        }else if (inside_small(vec, sm2_off, sm_delta)){
            //skip the point
        }else if (on_line(vec, lval, delta)){
            //skip the point
        }else{
            //create outlier
            for (int k = 0; k < dim; k++) {
                *out << vec[k] << " ";
            }
            *out << label << endl;
            i++;
        }

    }
    delete [] vec;
}

double scale(double value, double fromRangeMin, double fromRangeMax, double toRangeMin, double toRangeMax) {
    return ((value - fromRangeMin) * (toRangeMax - toRangeMin) / (fromRangeMax - fromRangeMin) + toRangeMin);
}

/**
* check if points coordinates are inside big circle
*/
bool inside_big(double* vec, double delta){
    double total, tmp, r, a, rsq;
    int i;

    total = 0.0;
    r = (b_max - xmin) / 2.0;
    rsq = r * r + delta;
    a = b_max - r;
    for (i = 0; i < dim; i++) {
        tmp = (vec[i] - a);
        total += tmp * tmp;
    }
    if (total <= rsq) {
        return true;
    }
    return false;
}

/**
* check if points coordinates are inside ellipsoid
*/
bool inside_elly(double* vec, double xoffset, double delta){
    double total, tmp, a, b, yoffset, abshalf;
    int i;

    //x: [xmin, half - gap]
    //y: [xmax - e_max, xmax]
    total = 0.0;
    //half size of are we paint
    abshalf = (xmax - xmin) / 2.0;
    a = (abshalf - gap) / 2.0;
    //TODO: this might not work with shift
    b = e_max / 2.0;
    yoffset = abshalf - b;
    a += delta;
    b += delta;
    for (i = 0; i < dim; i++) {
        if(i == 0){
            tmp = (vec[i] - xoffset) / a;
        }else {
            tmp = (vec[i] - yoffset) / b;
        }
        total += tmp * tmp;
    }
    if (total <= 1.0) {
        return true;
    }
    return false;
}

bool inside_small(double* vec, double yoffset, double delta){
    double total, tmp, r, rsq, xoffset;
    int i;
    //x: [xmax - quad, xmax]
    //y: [half + quad, half + 2 * quad]
    total = 0.0;
    r = quad / 2.0;
    rsq = r * r + delta;
    xoffset = xmax - r;
    for (i = 0; i < dim; i++) {
        if(i == 0){
            tmp = (vec[i] - xoffset);
        }else {
            tmp = (vec[i] - yoffset -r);
        }
        total += tmp * tmp;
    }
    if (total <= rsq) {
        return true;
    }
    return false;
}

bool on_line(double* vec, double yval, double delta){
    if(vec[0] > (half - gap) && vec[0] < (half + gap)){
        if(vec[1] > (yval - delta) && vec[1] < (yval + delta)){
            return true;
        }
    }
    return false;
}