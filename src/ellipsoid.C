
/***************************************************

ellipsoid.cc   Copyright (C) 2005 Joshua Knowles

Create elongated ellipsoidal clusters, compactly
arranged, in a high dimensional space


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



** Compile: **

   $ g++ ellipsoid.cc -o elly -lm -Wall -pedantic -O3


** Run: **

   $  ./elly -h
   for instructions on command line parameters


** Output: **

   The output to stdout is a file with d+1 columns where
   d is the number of dimensions of the data set. Each
   row is a data point, with the last column indicating the
   class label, i.e. the cluster to which the point belongs.


   In addition, there is some output to stderr.
   The overall deviation value is what the EA tries to
   minimize, so that the clusters are arranged to nearly
   overlap. The violation values indicate if the clusters
   overlap each other or not - no violation indicates that
   they are separated. At the end of the EA run, the lowest
   violation solution in the population is selected. The
   violation of this solution is the last line of output to
   stderr, and is labelled "Lowest_violation".
   By checking this value, one can see if the clusters overlap too
   much. Small violation values ( < 0.01 ) are fine
   and indicate that the clusters are close but do not overlap
   very much.

******************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#define PI 3.141592654

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

double ran0(long *idum);
double gaussian(long *idum);
/* End Numerical recipes */


FILE *fp;

typedef struct popmem
{
  double **shift;
  double deviation;
  double violation;
}P;

P *pop;


typedef struct clust
{
  double **datitem;
  double *focus;
  double *shift; // vector which shifts the whole away thing from the origin
  double length;
  double ecc;
  int npoints;
}C;

C *clu;

int npoints;
int nclust;
int nobjs; // dimension of hyperspace
long int seed;


double Eucdist(double *a, double *b);
void bounding_vector();
void cluster_at_origin(double length, double ecc, int npoints, double **datitem, double *foc);
void usage();
double overall_dev();
double feasible();


int main(int argc, char **argv)
{
  int i;
  nobjs=2;
  nclust=5;

  if(argc<=2)
    {
      usage();
      exit(1);
    }

  if((argc>2)&&(argc%2==1))
    {
      for(i=1;i<argc-1;i+=2)
	{
	  if(strcmp("-k", argv[i])==0)
            {
              nclust = atoi(argv[i+1]);
            }
	  else if(strcmp("-s", argv[i])==0)
	    {
	      seed = atol(argv[i+1]);
	    }
	  else if(strcmp("-d", argv[i])==0)
	    {
	      nobjs = atoi(argv[i+1]);
	    }
	  else if(strcmp("-n", argv[i])==0)
	    {
	      npoints = atoi(argv[i+1]);
	    }
	  else
	    {
	      fprintf(stderr, "Unrecognised parameter sequence.\n$ ./elly -h\nfor help.\n");
	      exit(1);
	    }

	}
    }
  else
    {
      fprintf(stderr,"Unrecognised parameter sequence.\n$ ./elly -h\nfor help.\n");
      exit(1);
    }
  ran0(&seed);


  double longlength=0;
  clu = (C *)malloc(nclust*sizeof(C));
  for(int c=0;c < nclust; c++)
    {
      double len=1.0+(RN*2.0);  // the length of the major axis of the ellipsoid
      double ecc=0.05+(RN*0.1); // 1 - sumdist, where sumdist is the maximum allowed sum of the two distances between a point and the two foci. (This relates to the eccentricity of the cluster, hence the name)
      clu[c].shift=(double *)malloc(nobjs*sizeof(double));
      for(int j=0;j<nobjs; j++)
	clu[c].shift[j]=c*0.2;
      clu[c].focus=(double *)malloc(nobjs*sizeof(double));
      if(nclust<20)
	clu[c].npoints=50+int(RN*(501-50));  // cluster size range when the number of clusters is smaller
      else
	clu[c].npoints=10+int(RN*(101-10)); // cluster size range when the number of clusters is larger
      clu[c].datitem=(double **)malloc(clu[c].npoints*sizeof(double *));
      for(int i=0;i<clu[c].npoints;i++)
	clu[c].datitem[i]=(double *)malloc(nobjs*sizeof(double));

      if(len>longlength)
	longlength=len;
      clu[c].length=len;
      clu[c].ecc=ecc;
      fprintf(stderr, "generating cluster %d:\n", c+1);
      cluster_at_origin(clu[c].length,clu[c].ecc,clu[c].npoints, clu[c].datitem, clu[c].focus);
    }

  // begin the steady-state EA to arrange the clusters
  fprintf(stderr, "Now arranging cluster origins...\n");

  // allocate and initialize the population of shift vectors
  int popsize=50;
  pop = (P *)malloc(popsize*sizeof(P));
  for(int m=0;m<popsize;m++)
    {
      pop[m].shift = (double **)malloc(nclust*sizeof(double *));
      for (int c=0;c<nclust;c++)
	{
	  pop[m].shift[c]=(double *)malloc(nobjs*sizeof(double));
	}

      for(int c=0;c<nclust;c++)
	for(int j=0;j<nobjs;j++)
	  pop[m].shift[c][j]=GN*sqrt((0.15*0.15)/double(nobjs))*longlength*0.67;
    }

  int best=-1;
  int worst=-1;
  double min_dev=1.0e100;
  double worst_dev=0;
  for(int m=0;m<popsize;m++)
    {
      for(int c=0;c<nclust;c++)
	for(int j=0;j<nobjs;j++)
	  clu[c].shift[j]=pop[m].shift[c][j];


      pop[m].deviation = overall_dev();
      pop[m].violation = feasible();

      pop[m].deviation += pop[m].violation*50*pop[m].deviation;

      if(pop[m].deviation < min_dev)
	{
	  min_dev = pop[m].deviation;
	  best=m;
	}
      if(pop[m].deviation > worst_dev)
	{
	  worst_dev = pop[m].deviation;
	  worst=m;
	}
    }

  double **newshift;
  newshift = (double **)malloc(nclust*sizeof(double *));
  for (int c=0;c<nclust;c++)
    {
      newshift[c]=(double *)malloc(nobjs*sizeof(double));
    }

  int iter=0;
  int parent_a;
  double viol=1.0;
  while(iter<500) // total number of EA iterations
    {
      fprintf(stderr,"Iteration= %d\n", iter);
      int x=int(RN*popsize);
      int y=int(RN*popsize);
      // binary tournament selection
      if(pop[x].deviation < pop[y].deviation)
	parent_a=x;
      else
	parent_a=y;

      for(int c=0;c<nclust;c++)
	for(int j=0;j<nobjs;j++)
	  {
	    if(RN < 1.0/(nclust*nobjs)) // mutation
	      {
		double gaussdev;
		gaussdev=GN*longlength*sqrt((0.15*0.15)/double(nobjs))*(double(nclust)/60.0);
		newshift[c][j]=gaussdev;
	      }
	    else
	      newshift[c][j]=pop[parent_a].shift[c][j];
	  }
      for(int c=0;c<nclust;c++)
	for(int j=0;j<nobjs;j++)
	  clu[c].shift[j]=newshift[c][j];


      double newdev = overall_dev();
      viol=feasible();
      newdev+= viol*50*newdev;

      fprintf(stderr,"Violation= %g\n", viol);
      fprintf(stderr,"New_deviation= %g Worst_deviation= %g\n",newdev, worst_dev);
      if(newdev < worst_dev)
	{
	  for(int c=0;c<nclust;c++)
	    for(int j=0;j<nobjs;j++)
	      pop[worst].shift[c][j]=newshift[c][j];
	  pop[worst].deviation=newdev;
	}

      if(newdev < min_dev)
	{
	  min_dev = newdev;
	  best=worst;
	}
      iter++;

      worst_dev=0.0;
      for(int m=0;m<popsize;m++)
	if(pop[m].deviation > worst_dev)
	  {
	    worst_dev = pop[m].deviation;
	    worst=m;
	  }
      fprintf(stderr, "Worst_overall_deviation= %g\n",worst_dev);
    }

  fprintf(stderr, "Minimum_overall_deviation= %g\n",min_dev);

   for(int c=0;c<nclust;c++)
     for(int j=0;j<nobjs;j++)
       clu[c].shift[j]=pop[best].shift[c][j];


   double lowest=1.0;
   for(int m=0;m<popsize;m++)
     {
       for(int c=0;c<nclust;c++)
	 for(int j=0;j<nobjs;j++)
	   clu[c].shift[j]=pop[m].shift[c][j];
       viol=feasible();
       if(viol < lowest)
	 {
	   lowest=viol;
	   best=m;
	 }
     }

   fprintf(stderr, "Lowest_violation= %g\n", lowest);

   for(int c=0;c<nclust;c++)
	for(int j=0;j<nobjs;j++)
	  clu[c].shift[j]=pop[best].shift[c][j];

  // print the clusters
  for(int c=0;c < nclust; c++)
    {
      for(int n=0; n < clu[c].npoints; n++)
	{
	  for(int j=0;j<nobjs;j++)
	    {
	      printf("%g ", clu[c].datitem[n][j]+clu[c].shift[j]);
	      // fprintf(stderr, "c=%d shift=%g ", c, clu[c].shift[j]);
	    }
	  //  fprintf(stderr, "\n");
	  printf("%d\n", c);
	}
      //  printf("\n\n");
      //  fprintf(stderr, "\n\n");
    }

  exit(0);

}

double feasible()
{
  double *shiftitem_c;
  double *focus_d;
  shiftitem_c = (double *)malloc(nobjs*sizeof(double));
  focus_d = (double *)malloc(nobjs*sizeof(double));

  int violation=0;
  int ntests=0;
  for (int c=0; c<nclust;c++)
    {
      for(int d=0;d<nclust;d++)
	{
	  if(c!=d)
	    {
	      for (int n=0; n < clu[c].npoints; n++)
		{
		  for(int j=0;j<nobjs;j++)
		    {
		      focus_d[j]=clu[d].focus[j]+clu[d].shift[j];
		      shiftitem_c[j] = clu[c].datitem[n][j]+clu[c].shift[j];
		    }

		  if(Eucdist(shiftitem_c,clu[d].shift)+Eucdist(shiftitem_c, focus_d) < (1.0+clu[d].ecc)*clu[d].length)
		    {
		      violation++;
		    }
		  ntests++;
		}
	    }
	}
    }
  return(double(violation)/double(ntests));
}

double overall_dev()
{
  //measure the overall deviation of the set of clusters
  double *mean;
  int total=0;
  mean = (double *)malloc(nobjs*sizeof(double));

  for(int j=0;j<nobjs;j++)
    mean[j]=0.0;

  for(int c=0; c< nclust; c++)
    {
      for (int n=0; n < clu[c].npoints; n++)
	{
	  total++;
	  for(int j=0;j<nobjs;j++)
	    {
	      mean[j]+=clu[c].datitem[n][j]+clu[c].shift[j];
	    }
	}
    }
  for(int j=0;j<nobjs;j++)
    mean[j]/=(double)total;

  double dev=0.0;

  double *sp; //shifted point;
  sp = (double *)malloc(nobjs*sizeof(double));


  for(int c=0; c< nclust; c++)
    {
      for (int n=0; n < clu[c].npoints; n++)
	{
	  for(int j=0;j<nobjs;j++)
	    sp[j]=clu[c].datitem[n][j]+clu[c].shift[j];
	  dev+=Eucdist(sp,mean);
	}
    }
  fprintf(stderr, "Overall_deviation= %g\n", dev);
  return(dev);
}



void usage()
{
  fprintf(stderr, "\nYou have requested help. The command line parameters for this generator are:\n\n");
  fprintf(stderr, "$ ./elly [-k <nclust>] [-d <dimension>] [-s <seed>]\n\n");
  fprintf(stderr, "where all parameters are optional and:\n  <nclust> is a positive int >= 2\n  <dimension> is a positive int >= 2\n  <seed> is a long int.\n\n");
  fprintf(stderr, "Use the following line as an example\n\n$ ./elly -k 10 -d 20 -s 30298343 > 20d10c.dat\n");
  fprintf(stderr, "...and change options as required.\n\n");
  fprintf(stderr, "Note 1: the ranges of size of clusters, the ranges of lengths of major axes, and the ranges of eccentricities are all hard-coded. However, they are commented in the source code to facilitate adapting them.\n\n");
  fprintf(stderr, "Note 2: the useable range of the dimension and k parameters is currently d in [10,100] and k in [2,40]. Outside these ranges, adjustment of the initialization and mutation operators might be needed in order to generate clusters that are close together but not overlapping too much.\n\n");
}

double Eucdist(double *a, double *b)
{
  double sum=0.0;
  for(int i=0;i<nobjs;i++)
    {
      sum+=(a[i]-b[i])*(a[i]-b[i]);
    }
  return(sqrt(sum));
}

void cluster_at_origin(double length, double ecc, int npoints, double **datitem, double *foc)
{
  // length controls the distance between the two foci of the cluster;
  // ecc controls the eccentricity of the cluster: it is the constraint
  //   specifying the maximum allowed sum for the two distances of a point from the two
  //   foci. It is usually in the range [0,1] but can go from [0,infty] too.

  int i,j;
  double sum;
  double *cart;
  double *gauss;
  double *origin;
  double *sph;
  cart = (double *)malloc(nobjs*sizeof(double));
  gauss = (double *)malloc(nobjs*sizeof(double));
  origin = (double *)malloc(nobjs*sizeof(double));
  sph = (double *)malloc(nobjs*sizeof(double));

  // first generate the second on a unit hypersphere around the origin to get a random orientation
  sum=0.0;
  for(j=0;j<nobjs;j++)
    {
      origin[j]=0.0;
      gauss[j]=GN;
      if(RN<0.5)
     	gauss[j]*=-1;
      sum+=gauss[j]*gauss[j];
    }
  for(j=0;j<nobjs;j++)
    foc[j]=(1.0/sqrt(sum))*gauss[j]*length;


  // then, generate the points a Gaussian distributed distance in a random direction away from a
  // uniformly random point somewhere along the major axis of the ellipsoid between the foci. Accept
  // the point if it is within the ellipsoid.
  i=0;
  bool success=true;
  while(i<npoints)
    {
      double ax=RN;
      double rad=GN;
      sum=0.0;

      if(success)
	{
	  for(j=0;j<nobjs;j++)
	    {
	      gauss[j]=GN;
	      sum+=gauss[j]*gauss[j];
	    }
	  for(j=0;j<nobjs;j++)
	    {
	      sph[j]=(1.0/sqrt(sum))*gauss[j];
	      //  printf("%g ", sph[j]);
	    }
	}
      // printf("\n");

      for(j=0;j<nobjs;j++)
	{
	  cart[j]=ax*foc[j]+sph[j]*rad*ecc*length*2; // the Cartesian coordinates of the point
	  // printf("%g ", cart[j]);
	}

      double diff;

      // is the point inside the ellipsoid ?
      if( (diff = ((Eucdist(cart,origin)+Eucdist(cart,foc)) - (1.0+ecc)*length)) < 0)
	{
	  //  if(RN< diff/(ecc*(-length)))
	  //    {
	      for(j=0;j<nobjs;j++)
		{
		  // printf("%g ", cart[j]);
		  datitem[i][j]=cart[j];
		}
	      //  printf("\n");
	      i++;
	      fprintf(stderr,".");
	      success=true;
	      //  }
	}
      else
	success=false;

    }
  //  printf("\n\n");
  fprintf(stderr,"\n\n");

}


/* Copyright Numerical Recipes in C */

double ran0(long *idum)
{
	long k;
	double ans;

	*idum ^= MASK;
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	ans=AM*(*idum);
	*idum ^= MASK;
	return ans;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef MASK


/* Generate a N(0,1) r.v. */
double gaussian(long *idum)
{
    static int iset=0;
    static double gset;
    double fac,r,v1,v2;
    if  (iset == 0) {
        do {
            v1=2.0*RN-1.0;
            v2=2.0*RN-1.0;
            r=v1*v1+v2*v2;
        } while (r >= 1.0 || r == 0.0);
        fac=sqrt(-2.0*log(r)/r);
        gset=v1*fac;
        iset=1;
        return v2*fac;
    } else {
        iset=0;
        return gset;
    }
}
