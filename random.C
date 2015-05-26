#include "random.h"

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
