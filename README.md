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

![mult_generator](https://raw.githubusercontent.com/deric/handl-data-generators/screens/img/2d-4c-no9.png)

## elly

Ellipsoid generator

```
$ ./elly [-k <nclust>] [-d <dimension>] [-s <seed>]
```
where all parameters are optional and:
  * `<nclust>` is a positive int >= 2
  * `<dimension>` is a positive int >= 2
  * `<seed>` is a long int.


![elly example](https://raw.githubusercontent.com/deric/handl-data-generators/screens/img/elly-2d10c13s.png)

## cure

CURE data sets generator. See Guha, Sudipto, Rajeev Rastogi, and Kyuseok Shim. "CURE: an efficient clustering algorithm for
large databases." ACM SIGMOD Record. Vol. 27. No. 2. ACM, 1998. for more details.

The distribution of data points is just approximated

```
$ ./cure -n <npoints> [-d <dimension>] [-s <seed>] [-l <x/y min>] [-m <x/y max>] [-t type of data]
```
where:
  * `-l` minimal x/y value
  * `-m` maximal x/y value
  * `-t` type of dataset, currently supports values 0-1

![cure t0](https://raw.githubusercontent.com/deric/handl-data-generators/screens/img/cure-t0-2k-2d.png)

![cure t1](https://raw.githubusercontent.com/deric/handl-data-generators/screens/img/cure-t1-2k-2d.png)

## Authors

  * Julia Handl
  * Joshua Knowles
  * John Burkardt
  * Tomas Barton
