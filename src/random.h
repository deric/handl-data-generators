/* Random number generator */
/* Copyright Numerical Recipes in C */
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876  
#define RN ran0(&seed)

double ran0(long *idum);
/* End copyright Numerical Recipes in C */


/*
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef MASK
*/

/* End copyright Numerical Recipes in C */
