/*  Synthetic test data for clustering
    Copyright (C) 2004 Julia Handl
    Email: Julia.Handl@gmx.de

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/


using namespace std;
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include "random.h"
#include "gasdev.h"
#include "random_data.h"
#include <cmath>


#define DIM 2        // dimensionality of the data
#define NUM 40      // number of clusters

#define MAXMU 10    // mean in each dimension is in range [0,MAXMU]
#define MINMU -10
#define MINSIGMA 0
#define MAXSIGMA 20*sqrt(DIM) // standard deviation (to be added on top
// of row sum in each dimension is in range [0,MAXSIGMA]
#define MAXSIZE 100  // size of each cluster is in range [MINSIZE,MAXSIZE]
#define MINSIZE 10
#define RUNS 10      // number of data sets to be generated


double mean[DIM];
double sigma[DIM*DIM];
int size[NUM];
double * points;
double pp[DIM*NUM*(MAXSIZE+MINSIZE)];
int pctr;
int number;
long int idum;
int idum2;

double meanmem[NUM][DIM];
double sigmamem[NUM][DIM];

double myabs(double x);
double max(double x, double y);
double min(double x, double y);
double square(double x);

void generate_config();
ofstream * out;


int main(int argc, char ** argv) {
  idum = 2717239;
  idum2 = 3417624;

  for (int j=0; j<RUNS; j++) {
    char name[50];
    sprintf(name, "%dd-%dc-no%d.dat", DIM, NUM, j);
    cerr << "Data written to " << name << endl;
    out = new ofstream(name);
    generate_config();
    cout << "Run " << j << " with " << number << " data items." << endl;
    out->close();
  }
}


void generate_config() {

  number = 0;
  pctr = 0;



  for (int i=0; i<NUM; i++) {
   	size[i] = MINSIZE + (int)ceil(ran0(&idum)*MAXSIZE);

	// Generate mean values
	for (int j=0; j<DIM; j++) {
	  mean[j] = -MAXMU+ran0(&idum)*(MAXMU-MINMU);
	  meanmem[i][j] = mean[j];


	  // Generate off-diagonal entries
	  for (int l=0; l<j; l++) {
	    if (ran0(&idum) < 0.5) {
	      sigma[j*DIM+l] = square(ran0(&idum));
	    }
	    else {
	      sigma[j*DIM+l] = -square(ran0(&idum));
	    }
	    sigma[l*DIM+j] = sigma[j*DIM+l];
	    }
	}

	// Generate diagonal entries
	for (int j=0; j<DIM; j++) {
	  double a=0.0;
	  for (int l=0; l<DIM; l++) {
	    if (l==j) continue;
	    else a += myabs(sigma[j*DIM+l]);
	  }

	  sigma[j*DIM+j] = a+MINSIGMA+square(ran0(&idum))*(MAXSIGMA-MINSIGMA);

	  sigmamem[i][j] = sigma[j*DIM+j];
	}

	// Construct and sample from multivariate normal distribution
	double * r = dpo_fa(DIM, sigma);
	if ( !r )
	  {
	    cout << "\n";
	    cout << "TEST04 - Fatal error!\n";
	    cout << "  Variance-covariance matrix factorization failed.\n";
	    exit ( 1 );
	  }
	points = normal_multivariate(DIM,size[i],r,mean, &idum2);

	// Check for violations and reject cluster if violation occurs
	//int violation = 0;
	if (i > 0) {
	  int found = 0;
	  for (int k=0; k<size[i]; k++) {
	    double mind = 1e10;
	    for (int l=0; l<size[i]; l++) {
	      double d = 0.0;
	      if (l==k) continue;
	      for (int j=0;j<DIM;j++) {
		d += square(points[k*DIM+j]-points[l*DIM+j]);
	      }
	      d = sqrt(d);
	      mind = min(d, mind);
	    }
	    double mind2 = 1e10;
	    for (int l=0; l<number; l++) {
	      double d = 0.0;
	      for (int j=0;j<DIM;j++) {
		d += square(points[k*DIM+j]-pp[l*DIM+j]);
	      }
	      d = sqrt(d);
	      mind2 = min(d, mind2);
	    }
	    // Accept cluster only if the nearest neighbour of each
	    // data item belongs to the same cluster
	    if (mind2 < mind) {
	      found = 1;
	      break;
	    }
	  }

	  if (found == 1) {
	    cerr << "Cluster " << i << " invalid" << endl;
	    i--;
	    continue;
	  }
	}

	number += size[i];

	// Print data to file
	for (int j=0; j<size[i]; j++) {
	  for (int k=0; k<DIM; k++) {
	    *out << points[j*DIM+k] << " ";
	  }
	  *out << i << endl;
	}
	for (int j=0; j<size[i]; j++) {
	  for (int k=0; k<DIM; k++) {
	    pp[pctr++] = points[j*DIM+k];
	  }
	}
    }
}


double max(double x, double y) {
  if (x < y) return y;
  else return x;
}

double min(double x, double y) {
  if (x > y) return y;
  else return x;
}


double myabs(double x) {
  if (x < 0) return -x;
  else return x;
}


double square(double x) {
  return x*x;
}
