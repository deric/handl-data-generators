# J. Handl: data generators

Syntetic data generators for creating datasets with Gaussian distribution. The code was taken from [official website](http://personalpages.manchester.ac.uk/mbs/Julia.Handl/generators.html) and slightly modified for modern compilers.

## Requirements

You'll need `g++` compiler.

  * Debian: `apt-get install build-essential` should be enough

and run

```
$ make
```

## mult_generator

Currently does not take any parameters, all settings is hard-coded in constants:

```c
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
```

simply run:

```
$ ./mult_generator
```

## elly

Ellipsoid generator

```
$ ./elly [-k <nclust>] [-d <dimension>] [-s <seed>]
```
where all parameters are optional and:
  * `<nclust>` is a positive int >= 2
  * `<dimension>` is a positive int >= 2
  * `<seed>` is a long int.


## Authors

  * Julia Handl
  * Joshua Knowles
  * John Burkardt