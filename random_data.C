# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "random_data.H"

//******************************************************************************

double *brownian ( int m, int n, int *seed )

//******************************************************************************
//
//  Purpose:
//
//    BROWNIAN creates Brownian motion points.
//
//  Discussion:
//
//    A starting point is generated at the origin.  The next point
//    is generated at a uniformly random angle and a (0,1) normally
//    distributed distance from the previous point.
//
//    It is up to the user to rescale the data, if desired.
//
//  Modified:
//
//    19 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the dimension of the space.
//
//    Input, int N, the number of points.
//
//    Input, int *SEED, a seed for the random number generator.
//
//    Output, double BROWNIAN[M*N], the Brownian motion points.
//
{
  double *direction;
  int i;
  int j;
  double r;
  double *x;

  direction = new double[m];
  x = new double[m*n];
//
//  Initial point.
//
  j = 0;
  for ( i = 0; i < m; i++ )
  {
    x[i+j*m] = 0.0;
  }
//
//  Generate angles and steps.
//
  for ( j = 1; j < n; j++ )
  {
    r = d_normal_01 ( seed );
    r = fabs ( r );

    direction_random_nd ( m, seed, direction );

    for ( i = 0; i < m; i++ )
    {
      x[i+j*m] = x[i+(j-1)*m] + r * direction[i];
    }

  }

  delete [] direction;

  return x;
}
//******************************************************************************

double d_epsilon ( void )

//******************************************************************************
//
//  Purpose:
//
//    D_EPSILON returns the round off unit for double precision arithmetic.
//
//  Discussion:
//
//    D_EPSILON is a number R which is a power of 2 with the property that,
//    to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but 
//      1 = ( 1 + R / 2 )
//
//  Modified:
//
//    01 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double D_EPSILON, the double precision round-off unit.
//
{
  double r;

  r = 1.0E+00;

  while ( 1.0E+00 < ( double ) ( 1.0E+00 + r )  )
  {
    r = r / 2.0E+00;
  }

  return ( 2.0E+00 * r );
}
//*********************************************************************

double d_max ( double x, double y )

//*********************************************************************
//
//  Purpose:
//
//    D_MAX returns the maximum of two double precision values.
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output double D_MAX, the maximum of X and Y.
//
{
  if ( y < x )
  {
    return x;
  } 
  else
  {
    return y;
  }
}
//*********************************************************************

double d_min ( double x, double y )

//*********************************************************************
//
//  Purpose:
//
//    D_MIN returns the minimum of two double precision values.
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output double D_MIN, the minimum of X and Y.
//
{
  if ( x < y )
  {
    return x;
  } 
  else
  {
    return y;
  }
}
//******************************************************************************

int d_nint ( double x )

//******************************************************************************
//
//  Purpose:
//
//    D_NINT returns the nearest integer to a double precision real value.
//
//  Examples:
//
//        X         D_NINT
//
//      1.3         1
//      1.4         1
//      1.5         1 or 2
//      1.6         2
//      0.0         0
//     -0.7        -1
//     -1.1        -1
//     -1.6        -2
//
//  Modified:
//
//    26 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the real value.
//
//    Output, int D_NINT, the nearest integer to X.
//
{
  int s;

  if ( x < 0.0 )
  {
    s = -1;
  }
  else
  {
    s = 1;
  }

  return ( s * ( int ) ( fabs ( x ) + 0.5 ) );
}
//******************************************************************************

double d_normal_01 ( int *seed )

//******************************************************************************
//
//  Purpose:
//
//    D_NORMAL_01 samples the standard normal probability distribution.
//
//  Discussion:
//
//    The standard normal probability distribution function (PDF) has
//    mean 0 and standard deviation 1.
//
//  Method:
//
//    The Box-Muller method is used, which is efficient, but
//    generates two values at a time.
//
//  Modified:
//
//    16 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double D_NORMAL_01, a sample of the standard normal PDF.
//
{
# define PI 3.141592653589793

  double r1;
  double r2;
  static int used = 0;
  double x;
  static double y = 0.0E+00;
//
//  If we've used an even number of values so far, generate two more,
//  return one and save one.
//
  if ( ( used % 2 ) == 0 )
  {
    for ( ; ; )
    {
      r1 = d_uniform_01 ( seed );

      if ( r1 != 0.0E+00 )
      {
        break;
      }
    }

    r2 = d_uniform_01 ( seed );

    x = sqrt ( -2.0 * log ( r1 ) ) * cos ( 2.0 * PI * r2 );
    y = sqrt ( -2.0 * log ( r1 ) ) * sin ( 2.0 * PI * r2 );
  }
//
//  Otherwise, return the second, saved, value.
//
  else
  {
    x = y;
  }

  used = used + 1;

  return x;
# undef PI
}
//******************************************************************************

double d_pi ( void )

//******************************************************************************
//
//  Purpose:
//
//    D_PI returns the value of PI to 16 digits.
//
//  Modified:
//
//    17 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double D_PI, the value of PI.
//
{
  return 3.141592653589793;
}
//****************************************************************************

double d_random ( double rlo, double rhi, int *seed )

//****************************************************************************
//
//  Purpose:
//
//    D_RANDOM returns a random double in a given range.
//
//  Modified:
//
//    16 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double RLO, RHI, the minimum and maximum values.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double D_RANDOM, the randomly chosen value.
//
{
  double t;

  t = d_uniform_01 ( seed );

  t = ( 1.0E+00 - t ) * rlo + t * rhi;

  return t;
}
//******************************************************************************

double d_uniform_01 ( int *seed )

//******************************************************************************
//
//  Purpose:
//
//    D_UNIFORM_01 is a portable pseudorandom number generator.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2**31 - 1 )
//      unif = seed / ( 2**31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Modified:
//
//    11 August 2004
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, L E Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//  Parameters:
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double D_UNIFORM_01, a new pseudorandom variate, strictly between
//    0 and 1.
//
{
  int k;
  double r;

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }
//
//  Although SEED can be represented exactly as a 32 bit integer,
//  it generally cannot be represented exactly as a 32 bit real number!
//
  r = ( double ) ( *seed ) * 4.656612875E-10;

  return r;
}
//******************************************************************************

double *dge_mxv ( int m, int n, double a[], double x[] )

//******************************************************************************
//
//  Purpose:
//
//    DGE_MXV multiplies a DGE matrix times a vector.
//
//  Discussion:
//
//    The DGE storage format is used for a general M by N matrix.  A physical storage 
//    space is made for each logical entry.  The two dimensional logical
//    array is mapped to a vector, in which storage is by columns.
//
//  Modified:
//
//    11 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//    M must be positive.
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, double A[M*N], the SGE matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double DGE_MXV[M], the product A * x.
//
{
  double *b;
  int i;
  int j;

  b = new double[m];

  for ( i = 0; i < m; i++ )
  {
    b[i] = 0.0E+00;
    for ( j = 0; j < n; j++ )
    {
      b[i] = b[i] + a[i+j*m] * x[j];
    }
  }

  return b;
}
//******************************************************************************

void direction_random_nd ( int m, int *seed, double w[] )

//******************************************************************************
//
//  Purpose:
//
//    DIRECTION_RANDOM_ND generates a random direction vector in ND.
//
//  Discussion:
//
//    This is actually simply a random point on the unit sphere in ND.
//
//  Modified:
//
//    19 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the dimension of the space.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double W[M], a random direction vector, with unit norm.
//
{
  int i;
  double norm;
//
//  Sample the standard normal distribution.
//
  dvec_normal_01 ( m, seed, w );
//
//  Compute the length of the vector.
//
  norm = 0.0;
  for ( i = 0; i < m; i++ )
  {
    norm = norm + w[i] * w[i];
  }
  norm = sqrt ( norm );
//
//  Normalize the vector.
//
  for ( i = 0; i < m; i++ )
  {
    w[i] = w[i] / norm;
  }

  return;
}
//******************************************************************************

void dmat_print ( int m, int n, double a[], char *title )

//******************************************************************************
//
//  Purpose:
//
//    DMAT_PRINT prints a double precision matrix, with an optional title.
//
//  Discussion:
//
//    The doubly dimensioned array A is treated as a one dimensional vector,
//    stored by COLUMNS.  Entry A(I,J) is stored as A[I+J*M]
//
//  Modified:
//
//    29 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, double A[M*N], the M by N matrix.
//
//    Input, char *TITLE, a title to be printed.
//
{
  int i;
  int j;
  int jhi;
  int jlo;
//
  if ( 0 < s_len_trim ( title ) ) 
  {
    cout << "\n";
    cout << title << "\n";
  }

  for ( jlo = 0; jlo < n; jlo = jlo + 6 )
  {
    jhi = jlo + 6;
    if ( n < jhi ) 
    {
      jhi = n;
    }
    cout << "\n";
    cout << "   Col  ";
    for ( j = jlo; j < jhi; j++ )
    {
      cout << setw(5) << j + 1 << "       ";
    }
    cout << "\n";
    cout << "   Row\n";
    cout << "\n";

    for ( i = 0; i < m; i++ )
    {
      cout << setw(6) << i + 1 << "  ";
      for ( j = jlo; j < jhi; j++ )
      {
        cout << setw(10) << a[j*m+i] << "  ";
      }
      cout << "\n";
    }
  }

  return;
}
//******************************************************************************

double *dpo_fa ( int n, double a[] )

//******************************************************************************
//
//  Purpose:
//
//    DPO_FA factors a DPO matrix.
//
//  Discussion:
//
//    The DPO storage format is appropriate for a symmetric positive definite 
//    matrix and its inverse.  (The Cholesky factor of a DPO matrix is an
//    upper triangular matrix, so it will be in DGE storage format.)
//
//    Only the diagonal and upper triangle of the square array are used.
//    This same storage format is used when the matrix is factored by
//    DPO_FA, or inverted by DPO_INVERSE.  For clarity, the lower triangle
//    is set to zero.
//
//    The positive definite symmetric matrix A has a Cholesky factorization
//    of the form:
//
//      A = R' * R
//
//    where R is an upper triangular matrix with positive elements on
//    its diagonal.  
//
//  Modified:
//
//    04 February 2004
//
//  Reference:
//
//    Dongarra, Bunch, Moler, Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[N*N], the matrix in DPO storage.
//
//    Output, double DPO_FA[N*N], the Cholesky factor in DGE
//    storage, or NULL if there was an error.
//
{
  double *b;
  int i;
  int j;
  int k;
  double s;

  b = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      b[i+j*n] = a[i+j*n];
    }
  }

  for ( j = 0; j < n; j++ )
  {
    for ( k = 0; k <= j-1; k++ )
    {
      for ( i = 0; i <= k-1; i++ )
      {
        b[k+j*n] = b[k+j*n] - b[i+k*n] * b[i+j*n];
      }
      b[k+j*n] = b[k+j*n] / b[k+k*n];
    }

    s = b[j+j*n];
    for ( i = 0; i <= j-1; i++ )
    {
      s = s - b[i+j*n] * b[i+j*n];
    }

    if ( s <= 0.0E+00 )
    {
      delete [] b;
      return NULL;
    }

    b[j+j*n] = sqrt ( s );
  }
//
//  Since the Cholesky factor is in DGE format, zero out the lower triangle.
//
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < i; j++ )
    {
      b[i+j*n] = 0.0E+00;
    }
  }

  return b;
}
//******************************************************************************

double *dut_mxv ( int m, int n, double a[], double x[] )

//******************************************************************************
//
//  Purpose:
//
//    DUT_MXV multiplies an DUT matrix times a vector.
//
//  Discussion:
//
//    The DUT storage format is used for an M by N upper triangular matrix,
//    and allocates space even for the zero entries.
//
//  Modified:
//
//    28 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//    M must be positive.
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, double A[M*N], the DUT matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double DUT_MXV[M], the product A * x.
//
{
  double *b;
  int i;
  int j;

  b = new double[m];

  for ( i = 0; i < m; i++ )
  {
    b[i] = 0.0E+00;
    for ( j = i; j < n; j++ )
    {
      b[i] = b[i] + a[i+j*m] * x[j];
    }
  }

  return b;
}
//*****************************************************************************

void dvec_normal_01 ( int n, int *seed, double x[] )

//*****************************************************************************
//
//  Purpose:
//
//    DVEC_NORMAL_01 samples the standard normal probability distribution.
//
//  Discussion:
//
//    The standard normal probability distribution function (PDF) has
//    mean 0 and standard deviation 1.
//
//    This routine can generate a vector of values on one call.  It
//    has the feature that it should provide the same results
//    in the same order no matter how we break up the task.
//
//    Before calling this routine, the user may call RANDOM_SEED
//    in order to set the seed of the random number generator.
//
//  Method:
//
//    The Box-Muller method is used, which is efficient, but
//    generates an even number of values each time.  On any call
//    to this routine, an even number of new values are generated.
//    Depending on the situation, one value may be left over.
//    In that case, it is saved for the next call.
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of values desired.  If N is negative,
//    then the code will flush its internal memory; in particular,
//    if there is a saved value to be used on the next call, it is
//    instead discarded.  This is useful if the user has reset the
//    random number seed, for instance.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double X[N], a sample of the standard normal PDF.
//
{
  int i;
  int m;
  double pi = 3.141592653589793;
  double *r;
  static int made = 0;
  static int saved = 0;
  int xhi;
  int xlo;
  static double y = 0.0;
//
//  I'd like to allow the user to reset the random number seed.
//  But this won't work properly if we have a saved value Y.
//  I'm making a crock option that allows the user to signal
//  explicitly that any internal memory should be flushed,
//  by passing in a negative value for N.
//
  if ( n < 0 )
  {
    made = 0;
    saved = 0;
    y = 0.0;
    return;
  }
  else if ( n == 0 )
  {
    return;
  }
//
//  Record the range of X we need to fill in.
//
  xlo = 1;
  xhi = n;
//
//  Use up the old value, if we have it.
//
  if ( saved == 1 )
  {
    x[0] = y;
    saved = 0;
    xlo = 2;
  }
//
//  If we don't need any more values, return.
//
  if ( xhi - xlo + 1 == 0 )
  {
    return;
  }

  r = new double[n+1];
//
//  If we need just one new value, do that here to avoid null arrays.
//
  if ( xhi - xlo + 1 == 1 )
  {
    dvec_uniform_01 ( 2, seed, r );

    x[xhi-1] = sqrt ( -2.0 * log ( r[0] ) ) * cos ( 2.0 * pi * r[1] );
    y =        sqrt ( -2.0 * log ( r[0] ) ) * sin ( 2.0 * pi * r[1] );

    saved = 1;

    made = made + 2;
  }
//
//  If we require an even number of values, that's easy.
//
  else if ( ( ( xhi-xlo+1) % 2 ) == 0 )
  {
    m = ( xhi-xlo+1 ) / 2;

    dvec_uniform_01 ( 2*m, seed, r );

    for ( i = 0; i < 2*m; i = i + 2 )
    {
      x[xlo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * pi * r[i+1] );
      x[xlo+i]   = sqrt ( -2.0 * log ( r[i] ) ) * sin ( 2.0 * pi * r[i+1] );
    }

    made = made + xhi - xlo + 1;
  }
//
//  If we require an odd number of values, we generate an even number,
//  and handle the last pair specially, storing one in X(N), and
//  saving the other for later.
//
  else
  {
    xhi = xhi - 1;

    m = ( xhi-xlo+1 ) / 2 + 1;

    dvec_uniform_01 ( 2*m, seed, r );

    for ( i = 0; i < 2*m-2; i = i + 2 )
    {
      x[xlo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * pi * r[i+1] );
      x[xlo+i  ] = sqrt ( -2.0 * log ( r[i] ) ) * sin ( 2.0 * pi * r[i+1] );
    }

    x[n-1] = sqrt ( -2.0 * log ( r[2*m-2] ) ) * cos ( 2.0 * pi * r[2*m-1] );
    y =      sqrt ( -2.0 * log ( r[2*m-2] ) ) * sin ( 2.0 * pi * r[2*m-1] );

    saved = 1;

    made = made + xhi - xlo + 2;

  }

  delete [] r;

  return;
}
//********************************************************************

void dvec_print ( int n, double a[], char *title )

//********************************************************************
//
//  Purpose:
//
//    DVEC_PRINT prints a real vector.
//
//  Modified:
//
//    16 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, double A[N], the vector to be printed.
//
//    Input, char *TITLE, a title to be printed first.
//    TITLE may be blank.
//
{
  int i;

  if ( 0 < s_len_trim ( title ) )
  {
    cout << "\n";
    cout << title << "\n";
  }

  cout << "\n";
  for ( i = 0; i <= n-1; i++ ) 
  {
    cout << setw(6)  << i + 1 << "  " 
         << setw(14) << a[i]  << "\n";
  }

  return;
}
//******************************************************************************

void dvec_uniform_01 ( int n, int *seed, double r[] )

//******************************************************************************
//
//  Purpose:
//
//    DVEC_UNIFORM_01 fills a double precision vector with pseudorandom values.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2**31 - 1 )
//      unif = seed / ( 2**31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Modified:
//
//    19 August 2004
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, L E Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R[N], the vector of pseudorandom values.
//
{
  int i;
  int k;

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + 2147483647;
    }

    r[i] = ( double ) ( *seed ) * 4.656612875E-10;
  }

  return;
}
//********************************************************************

unsigned long get_seed ( void )

//********************************************************************
//
//  Purpose:
//
//    GET_SEED returns a random seed for the random number generator.
//
//  Modified:
//
//    17 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, unsigned long GET_SEED, a random seed value.
//
{
# define UNSIGNED_LONG_MAX 4294967295UL
  time_t clock;
  int i;
  int hours;
  int minutes;
  int seconds;
  struct tm *lt;
  static unsigned long seed = 0;
  time_t tloc;
//
//  If the internal seed is 0, generate a value based on the time.
//
  if ( seed == 0 )
  {
    clock = time ( &tloc );
    lt = localtime ( &clock );
//
//  Extract HOURS.
//
    hours = lt->tm_hour;
//
//  In case of 24 hour clocks, shift so that HOURS is between 1 and 12.
//
    if ( 12 < hours )
    {
      hours = hours - 12;
    }
//
//  Move HOURS to 0, 1, ..., 11
//
    hours = hours - 1;

    minutes = lt->tm_min;

    seconds = lt->tm_sec;

    seed = seconds + 60 * ( minutes + 60 * hours );
//
//  We want values in [1,43200], not [0,43199].
//
    seed = seed + 1;
//
//  Remap SEED from [1,43200] to [1,UNSIGNED_LONG_MAX].
//
    seed = ( unsigned long ) 
      ( ( ( double ) seed )
      * ( ( double ) UNSIGNED_LONG_MAX ) / ( 60.0E+00 * 60.0E+00 * 12.0E+00 ) );
  }
//
//  Never use a seed of 0.
//
  if ( seed == 0 )
  {
    seed = 1;
  }

  return seed;

#undef UNSIGNED_LONG_MAX
}
//******************************************************************************

double *grid_in_cube01 ( int m, int n, int center, int *seed )

//******************************************************************************
//
//  Purpose:
//
//    GRID_IN_CUBE01 generates a grid dataset in the unit hypercube.
//
//  Discussion:
//
//    N points are needed in an M dimensional space.
//
//    The points are to lie on a uniform grid of side N_SIDE.
//
//    Unless the N = N_SIDE**M for some N_SIDE, we can't use all the
//    points on a grid.  What we do is find the smallest N_SIDE
//    that's big enough, and randomly omit some points.
//
//    If N_SIDE is 4, then the choices in 1D are:
//
//    A: 0,   1/3, 2/3, 1
//    B: 1/5, 2/5, 3/5, 4/5
//    C: 0,   1/4, 2/4, 3/4
//    D: 1/4, 2/4, 3/4, 1
//    E: 1/8, 3/8, 5/8, 7/8
//
//  Modified:
//
//    17 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, int CENTER, specifies the 1D grid centering:
//    1: first point is 0.0, last point is 1.0;
//    2: first point is 1/(N+1), last point is N/(N+1);
//    3: first point is 0, last point is (N-1)/N;
//    4: first point is 1/N, last point is 1;
//    5: first point is 1/(2*N), last point is (2*N-1)/(2*N);
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double GRID_IN_CUBE01[M*N], the points.
//
{
  int i;
  int j;
  int n_grid;
  int n_side;
  double *r;
  int rank;
  int *rank_list;
  int *tuple;
//
//  Find the dimension of the smallest grid with N points.
//
  n_side = grid_side ( m, n );
//
//  We need to select N points out of N_SIDE**M set.
//
  n_grid = ( int ) pow ( ( double ) n_side, m );
//
//  Generate a random subset of N items from a set of size N_GRID.
//
  rank_list = new int[n];

  ksub_random2 ( n_grid, n, seed, rank_list );
//
//  Must make one dummy call to TUPLE_NEXT_FAST with RANK = -1.
//
  rank = -1;
  tuple = new int[m];
  tuple_next_fast ( n_side, m, rank, tuple );
//
//  Now generate the appropriate indices, and "center" them.
//
  r = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    rank = rank_list[j] - 1;

    tuple_next_fast ( n_side, m, rank, tuple );

    if ( center == 1 )
    {
      for ( i = 0; i < m; i++ )
      {
        r[i+j*m] = ( double ) ( tuple[i] - 1 ) / ( double ) ( n_side - 1 );
      }
    }
    else if ( center == 2 )
    {
      for ( i = 0; i < m; i++ )
      {
        r[i+j*m] = ( double ) ( tuple[i] ) / ( double ) ( n_side + 1 );
      }
    }
    else if ( center == 3 )
    {
      for ( i = 0; i < m; i++ )
      {
        r[i+j*m] = ( double ) ( tuple[i] - 1 ) / ( double ) ( n_side );
      }
    }
    else if ( center == 4 )
    {
      for ( i = 0; i < m; i++ )
      {
        r[i+j*m] = ( double ) ( tuple[i] ) / ( double ) ( n_side );
      }
    }
    else if ( center == 5 )
    {
      for ( i = 0; i < m; i++ )
      {
        r[i+j*m] = ( double ) ( 2 * tuple[i] - 1 ) / ( double ) ( 2 * n_side );
      }
    }
  }

  delete [] rank_list;
  delete [] tuple;

  return r;
} 
//******************************************************************************

int grid_side ( int m, int n )

//******************************************************************************
//
//  Purpose:
//
//    GRID_SIDE finds the smallest M-D grid containing at least N points.
//
//  Discussion:
//
//    Each coordinate of the grid will have N_SIDE distinct values.
//    Thus the total number of points in the grid is N_SIDE**M.
//    This routine seeks the smallest N_SIDE such that N <= N_SIDE**M.
//
//  Modified:
//
//    17 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Output, int GRID_SIDE, the length of one side of the smallest 
//    grid in M dimensions that contains at least N points.
//
{
  double exponent;
  int n_grid;
  int n_side;

  if ( n <= 0 )
  {
    n_side = 0;
    return n_side;
  }

  if ( m <= 0 )
  {
    n_side = -1;
    return n_side;
  }

  exponent = 1.0E+00 / ( double ) m;

  n_side = ( int ) pow ( n, exponent );

  if ( pow ( ( double ) n_side, ( double ) m ) < n )
  {
    n_side = n_side + 1;
  }

  return n_side;
}
//**********************************************************************

bool halton_base_check ( int ndim, int base[] )

//**********************************************************************
//
//  Purpose:
//
//    HALTON_BASE_CHECK is TRUE if BASE is legal.
//
//  Modified:
//
//    16 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int BASE[HALTON_NDIM], the Halton bases.
//    Each base must be greater than 1.
//
//    Output, bool HALTON_BASE_CHECK.
//
{
  int i;
  bool value;

  value = true;

  for ( i = 0; i < ndim; i++ )
  {
    if ( base[i] <= 1 ) 
    {
      cout << "\n";
      cout << "HALTON_BASE_CHECK - Fatal error!\n";
      cout << "  Bases must be greater than 1.\n";
      cout << "  base[" << i << "] = " << base[i] << "\n";
      value = false;
      break;
    }
  }

  return value;
}
//******************************************************************************

double *halton_in_circle01_accept ( int m, int n, int *seed )

//******************************************************************************
//
//  Purpose:
//
//    HALTON_IN_CIRCLE01_ACCEPT accepts Halton points in a unit circle.
//
//  Discussion:
//
//    The acceptance/rejection method is used.
//
//    The unit circle is centered at the origin and has radius 1.
//
//  Modified:
//
//    25 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the dimension of the space.
//
//    Input, int N, the number of points.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double HALTON_IN_CIRCLE01_ACCEPT[M*N], the points.
//
{
# define M 2

  int base[M];
  int have;
  int i;
  int j;
  int leap[M];
  int seed_vec[M];
  int step;
  double total;
  double u[M];
  double *x;
//
  x = new double[m*n];

  have = 0;

  for ( i = 0; i < M; i++ )
  {
    seed_vec[i] = 0;
  }
  for ( i = 0; i < M; i++ )
  {
    leap[i] = 1;
  }
  for ( i = 0; i < M; i++ )
  {
    base[i] = prime ( i + 1 );
  }

  while ( have < n )
  {
    step = *seed;

    i_to_halton ( m, step, seed_vec, leap, base, u );

    *seed = *seed + 1;

    total = 0.0;
    for ( i = 0; i < M; i++ )
    {
      u[i] = 2.0 * u[i] - 1.0;
      total = total + u[i] * u[i];
    }

    if ( total <= 1.0 )
    {
      for ( i = 0; i < m; i++ )
      {
        x[i+have*m] = u[i];
      }
      have = have + 1;
    }
  }

  return x;
}
//******************************************************************************

double *halton_in_circle01_map ( int m, int n, int *seed )

//******************************************************************************
//
//  Purpose:
//
//    HALTON_IN_CIRCLE01_MAP maps Halton points into a unit circle.
//
//  Discussion:
//
//    The unit circle is centered at the origin and has radius 1.
//
//  Modified:
//
//    25 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the dimension of the space.
//
//    Input, int N, the number of points.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double HALTON_IN_CIRCLE01_MAP[M*N], the points.
//
{
# define PI 3.141592653589793

  int base[1];
  int j;
  int leap[1];
  double *r;
  double rval;
  int step;
  int seed_vec[1];
  double *t;
  double *x;
//
  r = new double[n];
  t = new double[n];
  x = new double[m*n];

  step = 0;
  seed_vec[0] = *seed;
  leap[0] = 1;
  base[0] = prime ( 1 );

  i_to_halton_sequence ( 1, n, step, seed_vec, leap, base, r );

  for ( j = 0; j < n; j++ )
  {
    r[j] = sqrt ( r[j] );
  }

  step = 0;
  seed_vec[0] = *seed;
  leap[0] = 1;
  base[0] = prime ( 2 );

  i_to_halton_sequence ( 1, n, step, seed_vec, leap, base, t );

  for ( j = 0; j < n; j++ )
  {
    t[j] = 2.0 * PI * t[j];
  }

  for ( j = 0; j < n; j++ )
  {
    x[0+j*m] = r[j] * cos ( t[j] );
    x[1+j*m] = r[j] * sin ( t[j] );
  }

  *seed = *seed + n;

  delete [] r;
  delete [] t;

  return x;
# undef PI
}
//******************************************************************************

double *halton_in_cube01 ( int m, int n, int *seed )

//******************************************************************************
//
//  Purpose:
//
//    HALTON_IN_CUBE01 generates Halton points in the unit hypercube.
//
//  Modified:
//
//    17 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of elements.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double HALTON_IN_CUBE01[M*N], the points
//
{
  int *base;
  int i;
  int *leap;
  int *seed_vec;
  int step;
  double *x;
//
  base = new int[m];
  leap = new int[m];
  seed_vec = new int[m];
  x = new double[m*n];

  step = *seed;
  for ( i = 0; i < m; i++ )
  {
    seed_vec[i] = 0;
  }
  for ( i = 0; i < m; i++ )
  {
    leap[i] = 1;
  }
  for ( i = 0; i < m; i++ )
  {
    base[i] = prime ( i + 1 );
  }

  i_to_halton_sequence ( m, n, step, seed_vec, leap, base, x );

  *seed = *seed + n;

  delete [] base;
  delete [] leap;
  delete [] seed_vec;

  return x;
}
//**********************************************************************

bool halton_leap_check ( int ndim, int leap[] )

//**********************************************************************
//
//  Purpose:
//
//    HALTON_LEAP_CHECK is TRUE if LEAP is legal.
//
//  Modified:
//
//    16 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int LEAP[HALTON_NDIM], the successive jumps in the Halton sequence.
//    Each entry must be greater than 0.
//
//    Output, bool HALTON_LEAP_CHECK.
//
{
  int i;
  bool value;

  value = true;

  for ( i = 0; i < ndim; i++ )
  {
    if ( leap[i] < 1 ) 
    {
      cout << "\n";
      cout << "HALTON_LEAP_CHECK - Fatal error!\n";
      cout << "  Leap entries must be greater than 0.\n";
      cout << "  leap[" << i << "] = " << leap[i] << "\n";
      value = false;
      break;
    }
  }

  return value;
}
//**********************************************************************

bool halton_n_check ( int n )

//**********************************************************************
//
//  Purpose:
//
//    HALTON_N_CHECK is TRUE if N is legal.
//
//  Modified:
//
//    16 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of the subsequence.
//    N must be positive.
//
//    Output, bool HALTON_N_CHECK.
//
{
  bool value;

  if ( n < 1 ) 
  {
    cout << "\n";
    cout << "HALTON_N_CHECK - Fatal error!\n";
    cout << "  N < 0.";
    cout << "  N = " << n << "\n";
    value = false;
  }
  else
  {
    value = true;
  }

  return value;
}
//**********************************************************************

bool halton_ndim_check ( int ndim )

//**********************************************************************
//
//  Purpose:
//
//    HALTON_NDIM_CHECK is TRUE if NDIM is legal.
//
//  Modified:
//
//    16 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NDIM, the spatial dimension.
//    NDIM must be positive.
//
//    Output, bool HALTON_NDIM_CHECK.
//
{
  bool value;

  if ( ndim < 1 ) 
  {
    cout << "\n";
    cout << "HALTON_NDIM_CHECK - Fatal error!\n";
    cout << "  NDIM < 0.";
    cout << "  NDIM = " << ndim << "\n";
    value = false;
  }
  else
  {
    value = true;
  }

  return value;
}
//**********************************************************************

bool halton_seed_check ( int ndim, int seed[] )

//**********************************************************************
//
//  Purpose:
//
//    HALTON_SEED_CHECK is TRUE if SEED is legal.
//
//  Modified:
//
//    16 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int SEED[NDIM], the Halton sequence index
//    corresponding to STEP = 0.  Each entry must be 0 or greater.
//
//    Output, bool HALTON_SEED_CHECK.
{
  int i;
  bool value;

  value = true;

  for ( i = 0; i < ndim; i++ )
  {
    if ( seed[i] < 0 ) 
    {
      cout << "\n";
      cout << "HALTON_SEED_CHECK - Fatal error!\n";
      cout << "  SEED entries must be nonnegative.\n";
      cout << "  seed[" << i << "] = " << seed[i] << "\n";
      value = false;
      break;
    }
  }

  return value;
}
//**********************************************************************

bool halton_step_check ( int step )

//**********************************************************************
//
//  Purpose:
//
//    HALTON_STEP_CHECK is TRUE if STEP is legal.
//
//  Modified:
//
//    16 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int STEP, the index of the subsequence element.
//    STEP must be 1 or greater.
//
//    Output, bool HALTON_STEP_CHECK.
{
  int i;
  bool value;

  if ( step < 0 ) 
  {
    cout << "\n";
    cout << "HALTON_STEP_CHECK - Fatal error!\n";
    cout << "  STEP < 0.";
    cout << "  STEP = " << step << "\n";
    value = false;
  }
  else
  {
    value = true;
  }

  return value;
}
//**********************************************************************

bool hammersley_base_check ( int ndim, int base[] )

//**********************************************************************
//
//  Purpose:
//
//    HAMMERSLEY_BASE_CHECK is TRUE if BASE is legal.
//
//  Modified:
//
//    16 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NDIM, the spatial dimension.
//
//    Input, int BASE[NDIM], the bases.
//
//    Output, bool HAMMERSLEY_BASE_CHECK.
{
  int i;
  bool value;

  value = true;

  for ( i = 0; i < ndim; i++ )
  {
    if ( base[i] == 0 || base[i] == 1 ) 
    {
      cout << "\n";
      cout << "HAMMERSLEY_BASE_CHECK - Fatal error!\n";
      cout << "  Bases may not be 0 or 1.\n";
      ivec_transpose_print ( ndim, base, "BASE:  " );
      value = false;
      break;
    }
  }

  return value;
}
//******************************************************************************

double *hammersley_in_cube01 ( int m, int n, int *seed )

//******************************************************************************
//
//  Purpose:
//
//    HAMMERSLEY_IN_CUBE01 computes Hammersley points in the unit hypercube.
//
//  Modified:
//
//    25 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of elements.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double HAMMERSLEY_IN_CUBE01[M*N], the points.
//
{
  int *base;
  int i;
  int *leap;
  int *seed_vec;
  int step;
  double *x;
//
  base = new int[m];
  leap = new int[m];
  seed_vec = new int[m];
  x = new double[m*n];

  step = *seed;
  for ( i = 0; i < m; i++ )
  {
    seed_vec[i] = 0;
  }
  for ( i = 0; i < m; i++ )
  {
    leap[i] = 1;
  }
  base[0] = -n;
  for ( i = 1; i < m; i++ )
  {
    base[i] = prime ( i );
  }

  i_to_hammersley_sequence ( m, n, step, seed_vec, leap, base, x );

  *seed = *seed + n;

  delete [] base;
  delete [] leap;
  delete [] seed_vec;

  return x;
}
//**********************************************************************

bool hammersley_leap_check ( int ndim, int leap[] )

//**********************************************************************
//
//  Purpose:
//
//    HAMMERSLEY_LEAP_CHECK is TRUE if LEAP is legal.
//
//  Modified:
//
//    16 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NDIM, the spatial dimension.
//
//    Input, int LEAP[NDIM], the successive jumps in the Hammersley sequence.
//    Each entry must be greater than 0.
//
//    Output, bool HAMMERSLEY_LEAP_CHECK.
{
  int i;
  bool value;

  value = true;

  for ( i = 0; i < ndim; i++ )
  {
    if ( leap[i] < 1 ) 
    {
      cout << "\n";
      cout << "HAMMERSLEY_LEAP_CHECK - Fatal error!\n";
      cout << "  Leap entries must be greater than 0.\n";
      ivec_transpose_print ( ndim, leap, "LEAP:  " );
      value = false;
      break;
    }
  }

  return value;
}
//**********************************************************************

bool hammersley_n_check ( int n )

//**********************************************************************
//
//  Purpose:
//
//    HAMMERSLEY_N_CHECK is TRUE if N is legal.
//
//  Modified:
//
//    16 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Output, bool HAMMERSLEY_N_CHECK.
{
  bool value;

  if ( n < 1 ) 
  {
    cout << "\n";
    cout << "HAMMERSLEY_N_CHECK - Fatal error!\n";
    cout << "  N < 0.";
    cout << "  N = " << n << "\n";
    value = false;
  }
  else
  {
    value = true;
  }

  return value;
}
//**********************************************************************

bool hammersley_ndim_check ( int ndim )

//**********************************************************************
//
//  Purpose:
//
//    HAMMERSLEY_NDIM_CHECK is TRUE if NDIM is legal.
//
//  Modified:
//
//    16 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NDIM, the spatial dimension.
//    NDIM must be positive.
//
//    Output, bool HAMMERSLEY_NDIM_CHECK.
{
  bool value;

  if ( ndim < 1 ) 
  {
    cout << "\n";
    cout << "HAMMERSLEY_NDIM_CHECK - Fatal error!\n";
    cout << "  NDIM < 0.";
    cout << "  NDIM = " << ndim << "\n";
    value = false;
  }
  else
  {
    value = true;
  }

  return value;
}
//**********************************************************************

bool hammersley_seed_check ( int ndim, int seed[] )

//**********************************************************************
//
//  Purpose:
//
//    HAMMERSLEY_SEED_CHECK is TRUE if SEED is legal.
//
//  Modified:
//
//    16 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NDIM, the spatial dimension.
//
//    Input, int SEED[NDIM], the Hammersley sequence index
//    corresponding to STEP = 0.  Each entry must be 0 or greater.
//
//    Output, bool HAMMERSLEY_SEED_CHECK.
//
{
  int i;
  bool value;

  value = true;

  for ( i = 0; i < ndim; i++ )
  {
    if ( seed[i] < 0 ) 
    {
      cout << "\n";
      cout << "HAMMERSLEY_SEED_CHECK - Fatal error!\n";
      cout << "  SEED entries must be nonnegative.\n";
      ivec_transpose_print ( ndim, seed, "SEED:  " );
      value = false;
      break;
    }
  }

  return value;
}
//**********************************************************************

bool hammersley_step_check ( int step )

//**********************************************************************
//
//  Purpose:
//
//    HAMMERSLEY_STEP_CHECK is TRUE is STEP is legal.
//
//  Modified:
//
//    16 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int STEP, the index of the subsequence element.
//    STEP must be 1 or greater.
//
//    Output, bool HAMMERSLEY_STEP_CHECK.
{
  bool value;

  if ( step < 0 ) 
  {
    cout << "\n";
    cout << "HAMMERSLEY_STEP_CHECK - Fatal error!\n";
    cout << "  STEP < 0.";
    cout << "  STEP = " << step << "\n";
    value = false;
  }
  else
  {
    value = true;
  }
  return value;
}
//**********************************************************************

int i_factorial ( int n )

//**********************************************************************
//
//  Purpose:
//
//    I_FACTORIAL returns N!.
//
//  Definition:
//
//    N! = Product ( 1 <= I <= N ) I
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the argument of the factorial function.
//    0 <= N.
//
//    Output, int I_FACTORIAL, the factorial of N.
//
{
  int fact;
  int i;
//
//  Check.
//
  if ( n < 0 )
  {
    cout << "\n";
    cout << "I_FACTORIAL - Fatal error!\n";
    cout << "  N < 0.\n";
    return 0;
  }

  fact = 1;

  for ( i = 2; i <= n; i++ )
  {
    fact = fact * i;
  }

  return fact;
}
//****************************************************************************

int i_max ( int i1, int i2 )

//****************************************************************************
//
//  Purpose:
//
//    I_MAX returns the maximum of two integers.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, are two integers to be compared.
//
//    Output, int I_MAX, the larger of I1 and I2.
//
//
{
  if ( i2 < i1 )
  {
    return i1;
  }
  else
  {
    return i2;
  }

}
//****************************************************************************

int i_min ( int i1, int i2 )

//****************************************************************************
//
//  Purpose:
//
//    I_MIN returns the smaller of two integers.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I_MIN, the smaller of I1 and I2.
//
//
{
  if ( i1 < i2 )
  {
    return i1;
  }
  else
  {
    return i2;
  }

}
//*********************************************************************

int i_modp ( int i, int j )

//*********************************************************************
//
//  Purpose:
//
//    I_MODP returns the nonnegative remainder of integer division.
//
//  Formula:
//
//    If 
//      NREM = I_MODP ( I, J ) 
//      NMULT = ( I - NREM ) / J
//    then
//      I = J * NMULT + NREM
//    where NREM is always nonnegative.
//
//  Comments:
//
//    The MOD function computes a result with the same sign as the
//    quantity being divided.  Thus, suppose you had an angle A,
//    and you wanted to ensure that it was between 0 and 360.
//    Then mod(A,360) would do, if A was positive, but if A
//    was negative, your result would be between -360 and 0.
//
//    On the other hand, I_MODP(A,360) is between 0 and 360, always.
//
//  Examples:
//
//        I         J     MOD  I_MODP   I_MODP Factorization
// 
//      107        50       7       7    107 =  2 *  50 + 7
//      107       -50       7       7    107 = -2 * -50 + 7
//     -107        50      -7      43   -107 = -3 *  50 + 43
//     -107       -50      -7      43   -107 =  3 * -50 + 43
//
//  Modified:
//
//    26 May 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the number to be divided.
//
//    Input, int J, the number that divides I.
//
//    Output, int I_MODP, the nonnegative remainder when I is 
//    divided by J.
//
{
  int value;

  if ( j == 0 )
  {
    cout << "\n";
    cout << "I_MODP - Fatal error!\n";
    cout << "  I_MODP ( I, J ) called with J = " << j << "\n";
    exit ( 1 );
  }

  value = i % j;

  if ( value < 0 )
  {
    value = value + abs ( j );
  }

  return value;
}
//********************************************************************

int i_random ( int ilo, int ihi, int *seed )

//********************************************************************
//
//  Purpose:
//
//    I_RANDOM returns a random integer in a given range.
//
//  Modified:
//
//    19 May 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ILO, IHI, the minimum and maximum values acceptable
//    for I.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, int I_RANDOM, the randomly chosen integer.
//
{
  int i;
  double r;
  double rhi;
  double rlo;
//
//  Pick a random number in (0,1).
//
  r = d_uniform_01 ( seed );
//
//  Set an interval [RLO,RHI] which contains the integers [ILO,IHI],
//  each with a "neighborhood" of width 1.
//
  rlo = ( ( double ) ilo ) - 0.5E+00;
  rhi = ( ( double ) ihi ) + 0.5E+00;
//
//  Set I to the integer that is nearest the scaled value of R.
//
  r = ( ( 1.0E+00 - r ) * rlo + r * rhi );

  if ( r < 0.0E+00 )
  {
    r = r - 0.5E+00;
  }
  else
  {
    r = r + 0.5E+00;
  }

  i = ( int ) r;
//
//  In case of oddball events at the boundary, enforce the limits.
//
  if ( i < ilo )
  {
    i = ilo;
  }

  if ( ihi < i )
  {
    i = ihi;
  }

  return i;
}
//**********************************************************************

void i_to_halton ( int ndim, int step, int seed[], int leap[], int base[], 
  double r[] )

//******************************************************************************
//
//  Purpose:
//
//    I_TO_HALTON computes one element of a leaped Halton subsequence.
//
//  Modified:
//
//    16 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    J H Halton,
//    On the efficiency of certain quasi-random sequences of points
//    in evaluating multi-dimensional integrals,
//    Numerische Mathematik,
//    Volume 2, 1960, pages 84-90.
//
//    J H Halton and G B Smith,
//    Algorithm 247: Radical-Inverse Quasi-Random Point Sequence,
//    Communications of the ACM,
//    Volume 7, 1964, pages 701-702.
//
//    Ladislav Kocis and William Whiten,
//    Computational Investigations of Low-Discrepancy Sequences,
//    ACM Transactions on Mathematical Software,
//    Volume 23, Number 2, 1997, pages 266-294.
//
//  Parameters:
//
//    Input, int NDIM, the spatial dimension.
//    1 <= NDIM is required.
//
//    Input, int STEP, the index of the subsequence element.
//    0 <= STEP is required.
//
//    Input, int SEED[NDIM], the Halton sequence index corresponding 
//    to STEP = 0.
//    0 <= SEED(1:NDIM) is required.
//
//    Input, int LEAP[NDIM], the successive jumps in the Halton sequence.
//    1 <= LEAP(1:NDIM) is required.
//
//    Input, int BASE[NDIM], the Halton bases.
//    1 < BASE(1:NDIM) is required.
//
//    Output, double R[NDIM], the STEP-th element of the leaped 
//    Halton subsequence.
//
{
  double base_inv;
  int digit;
  int i;
  int seed2;
//
//  Check the input.
//
  if ( !halton_ndim_check ( ndim ) )
  {
    exit ( 1 );
  }

  if ( !halton_step_check ( step ) )
  {
    exit ( 1 );
  }

  if ( !halton_seed_check ( ndim, seed ) )
  {
    exit ( 1 );
  }

  if ( !halton_leap_check ( ndim, leap ) )
  {
    exit ( 1 );
  }

  if ( !halton_base_check ( ndim, base ) )
  {
    exit ( 1 );
  }
//
//  Calculate the data.
//
  for ( i = 0; i < ndim; i++ )
  {
    seed2 = seed[i] + step * leap[i];

    r[i] = 0.0E+00;

    base_inv = 1.0E+00 / ( ( double ) base[i] );

    while ( seed2 != 0 )
    {
      digit = seed2 % base[i];
      r[i] = r[i] + ( ( double ) digit ) * base_inv;
      base_inv = base_inv / ( ( double ) base[i] );
      seed2 = seed2 / base[i];
    }
  }

  return;
}
//**********************************************************************

void i_to_halton_sequence ( int ndim, int n, int step, int seed[], int leap[],
  int base[], double r[] )

//******************************************************************************
//
//  Purpose:
//
//    I_TO_HALTON_SEQUENCE computes N elements of a leaped Halton subsequence.
//
//  Discussion:
//
//    The NDIM-dimensional Halton sequence is really NDIM separate
//    sequences, each generated by a particular base.
//
//    This routine selects elements of a "leaped" subsequence of the
//    Halton sequence.  The subsequence elements are indexed by a
//    quantity called STEP, which starts at 0.  The STEP-th subsequence
//    element is simply element
//
//      SEED(1:NDIM) + STEP * LEAP(1:NDIM)
//
//    of the original Halton sequence.
//
//
//    The data to be computed has two dimensions.
//
//    The number of data items is NDIM * N, where NDIM is the spatial dimension
//    of each element of the sequence, and N is the number of elements of the sequence.
//
//    The data is stored in a one dimensional array R.  The first element of the
//    sequence is stored in the first NDIM entries of R, followed by the NDIM entries
//    of the second element, and so on.
//
//    In particular, the J-th element of the sequence is stored in entries
//    0+(J-1)*NDIM through (NDIM-1) + (J-1)*NDIM.
//
//  Modified:
//
//    16 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    J H Halton,
//    On the efficiency of certain quasi-random sequences of points
//    in evaluating multi-dimensional integrals,
//    Numerische Mathematik,
//    Volume 2, 1960, pages 84-90.
//
//    J H Halton and G B Smith,
//    Algorithm 247: Radical-Inverse Quasi-Random Point Sequence,
//    Communications of the ACM,
//    Volume 7, 1964, pages 701-702.
//
//    Ladislav Kocis and William Whiten,
//    Computational Investigations of Low-Discrepancy Sequences,
//    ACM Transactions on Mathematical Software,
//    Volume 23, Number 2, 1997, pages 266-294.
//
//  Parameters:
//
//    Input, int NDIM, the spatial dimension.
//
//    Input, int N, the number of elements.
//
//    Input, int STEP, the index of the subsequence element.
//    0 <= STEP is required
//
//    Input, int SEED[NDIM], the Halton sequence index corresponding
//    to STEP = 0. 
//
//    Input, int LEAP[NDIM], the succesive jumps in the Halton sequence.
//
//    Input, int BASE[NDIM], the Halton bases.
//
//    Output, double R[NDIM*N], the next N elements of the
//    leaped Halton subsequence, beginning with element STEP.
//
{
  double base_inv;
  int digit;
  int i;
  int j;
  int *seed2;
//
//  Check the input.
//
  if ( !halton_ndim_check ( ndim ) )
  {
    exit ( 1 );
  }

  if ( !halton_n_check ( n ) )
  {
    exit ( 1 );
  }

  if ( !halton_step_check ( step ) )
  {
    exit ( 1 );
  }

  if ( !halton_seed_check ( ndim, seed ) )
  {
    exit ( 1 );
  }

  if ( !halton_leap_check ( ndim, leap ) )
  {
    exit ( 1 );
  }

  if ( !halton_base_check ( ndim, base ) )
  {
    exit ( 1 );
  }
//
//  Calculate the data.
//
  seed2 = new int[n];

  for ( i = 0; i < ndim; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      seed2[j] = seed[i] + ( step + j ) * leap[i];
    }

    for ( j = 0; j < n; j++ )
    {
      r[i+j*ndim] = 0.0E+00;
    }

    for ( j = 0; j < n; j++ )
    {
      base_inv = 1.0E+00 / ( ( double ) base[i] );

      while ( seed2[j] != 0 )
      {
        digit = seed2[j] % base[i];
        r[i+j*ndim] = r[i+j*ndim] + ( ( double ) digit ) * base_inv;
        base_inv = base_inv / ( ( double ) base[i] );
        seed2[j] = seed2[j] / base[i];
      }
    }
  }

  delete [] seed2;

  return;
}
//**********************************************************************

void i_to_hammersley_sequence ( int ndim, int n, int step, int seed[], int leap[],
  int base[], double r[] )

//******************************************************************************
//
//  Purpose:
//
//    I_TO_HAMMERSLEY_SEQUENCE computes N elements of a leaped Hammersley subsequence.
//
//  Discussion:
//
//    The NDIM-dimensional Hammersley sequence is really NDIM separate
//    sequences, each generated by a particular base.  If the base is 
//    greater than 1, a standard 1-dimensional
//    van der Corput sequence is generated.  But if the base is 
//    negative, this is a signal that the much simpler sequence J/(-BASE) 
//    is to be generated.  For the standard Hammersley sequence, the
//    first spatial coordinate uses a base of (-N), and subsequent
//    coordinates use bases of successive primes (2, 3, 5, 7, 11, ...).
//    This program allows the user to specify any combination of bases,
//    included nonprimes and repeated values.
//
//    This routine selects elements of a "leaped" subsequence of the
//    Hammersley sequence.  The subsequence elements are indexed by a
//    quantity called STEP, which starts at 0.  The STEP-th subsequence
//    element is simply element
//
//      SEED(1:NDIM) + STEP * LEAP(1:NDIM)
//
//    of the original Hammersley sequence.
//
//
//    The data to be computed has two dimensions.
//
//    The number of data items is NDIM * N, where NDIM is the spatial dimension
//    of each element of the sequence, and N is the number of elements of the sequence.
//
//    The data is stored in a one dimensional array R.  The first element of the
//    sequence is stored in the first NDIM entries of R, followed by the NDIM entries
//    of the second element, and so on.
//
//    In particular, the J-th element of the sequence is stored in entries
//    0+(J-1)*NDIM through (NDIM-1) + (J-1)*NDIM.
//
//  Modified:
//
//    20 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    J M Hammersley,
//    Monte Carlo methods for solving multivariable problems,
//    Proceedings of the New York Academy of Science,
//    Volume 86, 1960, pages 844-874.
//
//    Ladislav Kocis and William Whiten,
//    Computational Investigations of Low-Discrepancy Sequences,
//    ACM Transactions on Mathematical Software,
//    Volume 23, Number 2, 1997, pages 266-294.
//
//  Parameters:
//
//    Input, int NDIM, the spatial dimension.
//
//    Input, int N, the number of elements.
//
//    Input, int STEP, the index of the subsequence element.
//    0 <= STEP is required
//
//    Input, int SEED[NDIM], the Hammersley sequence index corresponding
//    to STEP = 0. 
//
//    Input, int LEAP[NDIM], the succesive jumps in the Hammersley sequence.
//
//    Input, int BASE[NDIM], the Hammersley bases.
//
//    Output, double R[NDIM*N], the next N elements of the
//    leaped Hammersley subsequence, beginning with element STEP.
//
{
  double base_inv;
  int digit;
  int i;
  int j;
  int *seed2;
  int temp;
//
//  Check the input.
//
  if ( !hammersley_ndim_check ( ndim ) )
  {
    exit ( 1 );
  }

  if ( !hammersley_n_check ( n ) )
  {
    exit ( 1 );
  }

  if ( !hammersley_step_check ( step ) )
  {
    exit ( 1 );
  }

  if ( !hammersley_seed_check ( ndim, seed ) )
  {
    exit ( 1 );
  }

  if ( !hammersley_leap_check ( ndim, leap ) )
  {
    exit ( 1 );
  }

  if ( !hammersley_base_check ( ndim, base ) )
  {
    exit ( 1 );
  }
//
//  Calculate the data.
//
  seed2 = new int[n];

  for ( i = 0; i < ndim; i++ )
  {
    if ( 1 < base[i] )
    {
      for ( j = 0; j < n; j++ )
      {
        seed2[j] = seed[i] + ( step + j ) * leap[i];
      }

      for ( j = 0; j < n; j++ )
      {
        r[i+j*ndim] = 0.0E+00;
      }

      for ( j = 0; j < n; j++ )
      {
        base_inv = 1.0E+00 / ( ( double ) base[i] );
  
        while ( seed2[j] != 0 )
        {
          digit = seed2[j] % base[i];
          r[i+j*ndim] = r[i+j*ndim] + ( ( double ) digit ) * base_inv;
          base_inv = base_inv / ( ( double ) base[i] );
          seed2[j] = seed2[j] / base[i];
        }
      }
    }
    else
    {
      for ( j = 0; j < n; j++ )
      {
        temp = ( seed[i] + ( step + j ) * leap[i] ) % ( -base[i] );

        r[i+j*ndim] = ( double ) ( temp ) 
                    / ( double ) ( -base[i] );
      }
    }
  }

  delete [] seed2;

  return;
}
//******************************************************************************

void ivec_transpose_print ( int n, int a[], char *title )

//******************************************************************************
//
//  Purpose:
//
//    IVEC_TRANSPOSE_PRINT prints an integer vector "transposed".
//
//  Example:
//
//    A = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 }
//    TITLE = "My vector:  "
//
//    My vector:      1    2    3    4    5
//                    6    7    8    9   10
//                   11
//
//  Modified:
//
//    03 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, int A[N], the vector to be printed.
//
//    Input, char *TITLE, a title to be printed first.
//    TITLE may be blank or NULL.
//
{
  int i;
  int ihi;
  int ilo;
  int title_len;

  if ( 0 < s_len_trim ( title ) )
  {
    title_len = strlen ( title );

    for ( ilo = 1; ilo <= n; ilo = ilo + 5 )
    {
      ihi = i_min ( ilo + 5 - 1, n );
      if ( ilo == 1 )
      {
        cout << title;
      }
      else
      {
        for ( i = 1; i <= title_len; i++ )
        {
          cout << " ";
        }
      }
      for ( i = ilo; i <= ihi; i++ )
      {
        cout << setw(12) << a[i-1];
      }
      cout << "\n";
    }
  }
  else
  {
    for ( ilo = 1; ilo <= n; ilo = ilo + 5 )
    {
      ihi = i_min ( ilo + 5 - 1, n );
      for ( i = ilo; i <= ihi; i++ )
      {
        cout << setw(12) << a[i-1];
      }
      cout << "\n";
    }
  }

  return;
}
//******************************************************************************

void ksub_random2 ( int n, int k, int *seed, int a[] )

//******************************************************************************
//
//  Purpose:
//
//    KSUB_RANDOM2 selects a random subset of size K from a set of size N.
//
//  Modified:
//
//    17 May 2003
//
//  Reference:
//
//    A Nijenhuis and H Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, the size of the set from which subsets are drawn.
//
//    Input, int K, number of elements in desired subsets.  K must
//    be between 0 and N.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, int A[K].  A(I) is the I-th element of the
//    output set.  The elements of A are in order.
//
{
  int available;
  int candidate;
  int have;
  int need;
  double r;

  if ( k < 0 || n < k )
  {
    cout << "\n";
    cout << "KSUB_RANDOM2 - Fatal error!\n";
    cout << "  N = " << n << "\n";
    cout << "  K = " << k << "\n";
    cout << "  but 0 <= K <= N is required!\n";
    exit ( 1 );
  }

  if ( k == 0 )
  {
    return;
  }

  need = k;
  have = 0;
  available = n;
  candidate = 0;

  for ( ; ; )
  {
    candidate = candidate + 1;

    r = d_uniform_01 ( seed );

    if ( r * ( double ) available <= ( double ) need )
    {
      need = need - 1;
      a[have] = candidate;
      have = have + 1;

      if ( need <= 0 )
      {
        break;
      }

    }

    available = available - 1;

  }

  return;
}
//******************************************************************************

double *normal ( int m, int n, double r[], double mu[], int *seed )

//******************************************************************************
//
//  Purpose:
//
//    NORMAL creates normally distributed points in M space.
//
//  Discussion:
//
//    The multivariate normal distribution for the M dimensional vector X
//    has the form:
//
//      pdf(X) = (2*pi*det(V))**(-M/2) * exp(-0.5*(X-MU)'*inverse(V)*(X-MU))
//
//    where MU is the mean vector, and V is a positive definite symmetric
//    matrix called the variance-covariance matrix.
//
//    This routine requires that the user supply the upper triangular
//    Cholesky factor R, which has the property that
//
//      V = R' * R
//
//    This factorization always exists if V is actually symmetric and
//    positive definite.  This factorization can be computed by the
//    routine SPO_FA.
//
//    The user also supplies the mean vector MU.
//
//  Modified:
//
//    19 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the dimension of the space.
//
//    Input, int N, the number of points.
//
//    Input, double R[M*M], the upper triangular Cholesky factor
//    of the variance-covariance matrix.
//
//    Input, double MU[M], the mean vector.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double NORMAL[M*N], the random points.
//
{
  int i;
  int j;
  int k;
  double *v;
  double *x;

  v = new double[m];
  x = new double[m*n];
//
//  Get a matrix V of normal data.
//  Compute X = MU + R' * V.
//  We actually carry out this computation in the equivalent form X' * R.
//
  for ( j = 0; j < n; j++ )
  {
    dvec_normal_01 ( m, seed, v );

    for ( i = 0; i < m; i++ )
    {
      x[i+j*m] = mu[i];
      for ( k = 0; k <= i; k++ )
      {
        x[i+j*m] = x[i+j*m] + v[k] * r[k+i*m];
      }
    }
  }

  delete [] v;

  return x;
}
//******************************************************************************

double *normal_circular ( int m, int n, int *seed )

//******************************************************************************
//
//  Purpose:
//
//    NORMAL_CIRCULAR creates circularly normal points in 2 space.
//
//  Modified:
//
//    20 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964, page 936.
//
//  Parameters:
//
//    Input, int M, the dimension of the space, which must be 2.
//
//    Input, int N, the number of points.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double NORMAL_CIRULAR[M*N], the random points.
//
{
# define PI 3.141592653589793

  int j;
  double *r;
  double *t;
  double *x;

  r = new double[n];
  t = new double[n];
  x = new double[m*n];
//
//  The angle varies uniformly from 0 to 2 pi.
//
  dvec_uniform_01 ( n, seed, t );

  for ( j = 0; j < n; j++ )
  {
    t[j] = 2.0 * PI * t[j];
  }
//
//  The radius is normally distributed.
//
  dvec_normal_01 ( n, seed, r );

  for ( j = 0; j < n; j++ )
  {
    x[0+j*m] = r[j] * cos ( t[j] );
    x[1+j*m] = r[j] * sin ( t[j] );
  }

  delete [] r;
  delete [] t;

  return x;
# undef PI
}
//******************************************************************************

double *normal_multivariate ( int m, int n, double r[], double mu[], 
  int *seed )

//******************************************************************************
//
//  Purpose:
//
//    NORMAL_MULTIVARIATE samples a multivariate normal distribution.
//
//  Discussion:
//
//    The multivariate normal distribution for the M dimensional vector X
//    has the form:
//
//      pdf(X) = (2*pi*det(V))**(-M/2) * exp(-0.5*(X-MU)'*inverse(V)*(X-MU))
//
//    where MU is the mean vector, and V is a positive definite symmetric
//    matrix called the variance-covariance matrix.
//
//    This routine samples points associated with the M-dimensional
//    normal distribution with mean MU and covariance matrix V.
//
//    This routine requires that the user supply the upper triangular
//    Cholesky factor R of V, which has the property that
//
//      V = R' * R
//
//    This factorization always exists if V is actually symmetric and
//    positive definite.  This factorization can be computed by the
//    routine DPO_FA.
//
//    The user also supplies the mean vector MU.
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Russell Cheng,
//    Random Variate Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998, pages 167-168.
//
//  Parameters:
//
//    Input, int M, the dimension of the space.
//
//    Input, int N, the number of points.
//
//    Input, double R[M*M], the upper triangular Cholesky factor
//    of the variance-covariance matrix.
//
//    Input, double MU[M], the mean vector.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double NORMAL_MULTIVARIATE[M*N], corresponding 
//    points associated with the multivariate normal distribution.
//
{
  int i;
  int j;
  int k;
  double *v;
  double *x;

  v = new double[m];
  x = new double[m*n];
//
//  Compute X = MU + R' * V.
//  We actually carry out this computation in the equivalent form MU + V' * R.
//
  for ( j = 0; j < n; j++ )
  {
    dvec_normal_01 ( m, seed, v );

    for ( i = 0; i < m; i++ )
    {
      x[i+j*m] = mu[i];
      for ( k = 0; k <= i; k++ )
      {
        x[i+j*m] = x[i+j*m] + v[k] * r[k+i*m];
      }
    }
  }

  delete [] v;

  return x;
}
//******************************************************************************

double *normal_simple ( int m, int n, int *seed )

//******************************************************************************
//
//  Purpose:
//
//    NORMAL_SIMPLE creates normally distributed points in M space.
//
//  Discussion:
//
//    The multivariate normal distribution has the form:
//
//      f(x) = (2*pi*det(V))**(-n/2) * exp(-0.5*(x-mu)'*inverse(V)*(x-mu))
//
//    where mu is the mean vector, and V is a positive definite symmetric
//    matrix called the variance-covariance matrix.
//
//    This routine implements the simplest version of a multivariate
//    normal distribution.  The variance-covariance matrix is the identity,
//    and the mean vector is entirely zero.  Thus, a sample on N points
//    is simply M*N scalar values generated under the univariate
//    normal distribution with zero mean and unit variance.
//
//  Modified:
//
//    19 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the dimension of the space.
//
//    Input, int N, the number of points.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double NORMAL_SIMPLE[M*N], the random points.
//
{
  double *x;

  x = new double[m*n];

  dvec_normal_01 ( m*n, seed, x );

  return x;
}
//******************************************************************************

double *polygon_centroid_2d ( int n, double v[] )

//******************************************************************************
//
//  Purpose:
//
//    POLYGON_CENTROID_2D computes the centroid of a polygon in 2D.
//
//  Formula:
//
//    Denoting the centroid coordinates by CENTROID, then
//
//      CENTROID(1) = Integral ( Polygon interior ) x dx dy / Area ( Polygon )
//      CENTROID(2) = Integral ( Polygon interior ) y dx dy / Area ( Polygon ).
//
//    Green's theorem states that
//
//      Integral ( Polygon boundary ) ( M dx + N dy ) =
//      Integral ( Polygon interior ) ( dN/dx - dM/dy ) dx dy.
//
//    Using M = 0 and N = x * x / 2, we get:
//
//      CENTROID(1) = 0.5 * Integral ( Polygon boundary ) x * x dy,
//
//    which becomes
//
//      CENTROID(1) = 1/6 Sum ( 1 <= I <= N )
//        ( X(I+1) + X(I) ) * ( X(I) * Y(I+1) - X(I+1) * Y(I))
//
//    where, when I = N, the index "I+1" is replaced by 1.
//
//    A similar calculation gives us a formula for CENTROID(2).
//
//  Modified:
//
//    17 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Gerard Bashein and Paul Detmer,
//    Centroid of a Polygon,
//    Graphics Gems IV, edited by Paul Heckbert,
//    AP Professional, 1994.
//
//  Parameters:
//
//    Input, int N, the number of sides of the polygonal shape.
//
//    Input, double V[2*N], the coordinates of the vertices
//    of the shape.
//
//    Output, double POLYGON_CENTROID_2D[2], the coordinates of the
//    centroid of the shape.
//
{
  double area;
  double *centroid;
  int i;
  int ip1;
  double temp;
//
  area = 0.0;
  centroid = new double[2];
  centroid[0] = 0.0;
  centroid[1] = 0.0;

  for ( i = 0; i < n; i++ )
  {
    if ( i < n-1 )
    {
      ip1 = i + 1;
    }
    else
    {
      ip1 = 0;
    }

    temp = ( v[0+i*2] * v[1+ip1*2] - v[0+ip1*2] * v[1+i*2] );

    area = area + temp;

    centroid[0] = centroid[0] + ( v[0+ip1*2] + v[0+i*2] ) * temp;
    centroid[1] = centroid[1] + ( v[1+ip1*2] + v[1+i*2] ) * temp;

  }

  area = area / 2.0;

  centroid[0] = centroid[0] / ( 6.0E+00 * area );
  centroid[1] = centroid[1] / ( 6.0E+00 * area );

  return centroid;
}
//******************************************************************************

int prime ( int n )

//******************************************************************************
//
//  Purpose:
//
//    PRIME returns any of the first PRIME_MAX prime numbers.
//
//  Discussion:
//
//    PRIME_MAX is 1600, and the largest prime stored is 13499.
//
//    Thanks to Bart Vandewoestyne for pointing out a typo, 18 February 2005.
//
//  Modified:
//
//    18 February 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964, pages 870-873.
//
//    Daniel Zwillinger,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996, pages 95-98.
//
//  Parameters:
//
//    Input, int N, the index of the desired prime number.
//    In general, is should be true that 0 <= N <= PRIME_MAX.
//    N = -1 returns PRIME_MAX, the index of the largest prime available.
//    N = 0 is legal, returning PRIME = 1.
//
//    Output, int PRIME, the N-th prime.  If N is out of range, PRIME
//    is returned as -1.
//
{
# define PRIME_MAX 1600

  int npvec[PRIME_MAX] = {
        2,    3,    5,    7,   11,   13,   17,   19,   23,   29,
       31,   37,   41,   43,   47,   53,   59,   61,   67,   71,
       73,   79,   83,   89,   97,  101,  103,  107,  109,  113,
      127,  131,  137,  139,  149,  151,  157,  163,  167,  173,
      179,  181,  191,  193,  197,  199,  211,  223,  227,  229,
      233,  239,  241,  251,  257,  263,  269,  271,  277,  281,
      283,  293,  307,  311,  313,  317,  331,  337,  347,  349,
      353,  359,  367,  373,  379,  383,  389,  397,  401,  409,
      419,  421,  431,  433,  439,  443,  449,  457,  461,  463,
      467,  479,  487,  491,  499,  503,  509,  521,  523,  541,
      547,  557,  563,  569,  571,  577,  587,  593,  599,  601,
      607,  613,  617,  619,  631,  641,  643,  647,  653,  659,
      661,  673,  677,  683,  691,  701,  709,  719,  727,  733,
      739,  743,  751,  757,  761,  769,  773,  787,  797,  809,
      811,  821,  823,  827,  829,  839,  853,  857,  859,  863,
      877,  881,  883,  887,  907,  911,  919,  929,  937,  941,
      947,  953,  967,  971,  977,  983,  991,  997, 1009, 1013,
     1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069,
     1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151,
     1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223,
     1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 
     1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 
     1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 
     1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 
     1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 
     1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, 
     1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 
     1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, 
     1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889,
     1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987,
     1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053,
     2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129,
     2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213,
     2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287,
     2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, 
     2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, 
     2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, 
     2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, 
     2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687, 
     2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741,
     2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, 
     2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903, 
     2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, 
     3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, 
     3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181, 
     3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, 
     3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, 
     3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413, 
     3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, 
     3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571,
     3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643,
     3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727,
     3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821,
     3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907,
     3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989,
     4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057,
     4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139,
     4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231,
     4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297,
     4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409,
     4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493, 
     4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583, 
     4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657, 
     4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751, 
     4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, 
     4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937, 
     4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003, 
     5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087, 
     5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179, 
     5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279,
     5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387,
     5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443,
     5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521,
     5527, 5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623, 5639,
     5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693,
     5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791,
     5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857,
     5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939,
     5953, 5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053,
     6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133,
     6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221,
     6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301,
     6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367, 
     6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473, 
     6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571, 
     6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673, 
     6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761, 
     6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833, 
     6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917, 
     6947, 6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997,
     7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103, 
     7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207, 
     7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297, 
     7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411, 
     7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499, 
     7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561, 
     7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643, 
     7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723, 
     7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829, 
     7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919,
     7927, 7933, 7937, 7949, 7951, 7963, 7993, 8009, 8011, 8017,
     8039, 8053, 8059, 8069, 8081, 8087, 8089, 8093, 8101, 8111,
     8117, 8123, 8147, 8161, 8167, 8171, 8179, 8191, 8209, 8219,
     8221, 8231, 8233, 8237, 8243, 8263, 8269, 8273, 8287, 8291,
     8293, 8297, 8311, 8317, 8329, 8353, 8363, 8369, 8377, 8387,
     8389, 8419, 8423, 8429, 8431, 8443, 8447, 8461, 8467, 8501,
     8513, 8521, 8527, 8537, 8539, 8543, 8563, 8573, 8581, 8597,
     8599, 8609, 8623, 8627, 8629, 8641, 8647, 8663, 8669, 8677,
     8681, 8689, 8693, 8699, 8707, 8713, 8719, 8731, 8737, 8741, 
     8747, 8753, 8761, 8779, 8783, 8803, 8807, 8819, 8821, 8831,
     8837, 8839, 8849, 8861, 8863, 8867, 8887, 8893, 8923, 8929,
     8933, 8941, 8951, 8963, 8969, 8971, 8999, 9001, 9007, 9011,
     9013, 9029, 9041, 9043, 9049, 9059, 9067, 9091, 9103, 9109,
     9127, 9133, 9137, 9151, 9157, 9161, 9173, 9181, 9187, 9199,
     9203, 9209, 9221, 9227, 9239, 9241, 9257, 9277, 9281, 9283,
     9293, 9311, 9319, 9323, 9337, 9341, 9343, 9349, 9371, 9377,
     9391, 9397, 9403, 9413, 9419, 9421, 9431, 9433, 9437, 9439,
     9461, 9463, 9467, 9473, 9479, 9491, 9497, 9511, 9521, 9533,
     9539, 9547, 9551, 9587, 9601, 9613, 9619, 9623, 9629, 9631,
     9643, 9649, 9661, 9677, 9679, 9689, 9697, 9719, 9721, 9733,
     9739, 9743, 9749, 9767, 9769, 9781, 9787, 9791, 9803, 9811,
     9817, 9829, 9833, 9839, 9851, 9857, 9859, 9871, 9883, 9887,
     9901, 9907, 9923, 9929, 9931, 9941, 9949, 9967, 9973,10007,
    10009,10037,10039,10061,10067,10069,10079,10091,10093,10099,
    10103,10111,10133,10139,10141,10151,10159,10163,10169,10177,
    10181,10193,10211,10223,10243,10247,10253,10259,10267,10271,
    10273,10289,10301,10303,10313,10321,10331,10333,10337,10343,
    10357,10369,10391,10399,10427,10429,10433,10453,10457,10459,
    10463,10477,10487,10499,10501,10513,10529,10531,10559,10567,
    10589,10597,10601,10607,10613,10627,10631,10639,10651,10657,
    10663,10667,10687,10691,10709,10711,10723,10729,10733,10739,
    10753,10771,10781,10789,10799,10831,10837,10847,10853,10859,
    10861,10867,10883,10889,10891,10903,10909,10937,10939,10949,
    10957,10973,10979,10987,10993,11003,11027,11047,11057,11059,
    11069,11071,11083,11087,11093,11113,11117,11119,11131,11149,
    11159,11161,11171,11173,11177,11197,11213,11239,11243,11251,
    11257,11261,11273,11279,11287,11299,11311,11317,11321,11329,
    11351,11353,11369,11383,11393,11399,11411,11423,11437,11443,
    11447,11467,11471,11483,11489,11491,11497,11503,11519,11527,
    11549,11551,11579,11587,11593,11597,11617,11621,11633,11657,
    11677,11681,11689,11699,11701,11717,11719,11731,11743,11777,
    11779,11783,11789,11801,11807,11813,11821,11827,11831,11833,
    11839,11863,11867,11887,11897,11903,11909,11923,11927,11933,
    11939,11941,11953,11959,11969,11971,11981,11987,12007,12011,
    12037,12041,12043,12049,12071,12073,12097,12101,12107,12109,
    12113,12119,12143,12149,12157,12161,12163,12197,12203,12211,
    12227,12239,12241,12251,12253,12263,12269,12277,12281,12289,
    12301,12323,12329,12343,12347,12373,12377,12379,12391,12401,
    12409,12413,12421,12433,12437,12451,12457,12473,12479,12487,
    12491,12497,12503,12511,12517,12527,12539,12541,12547,12553,
    12569,12577,12583,12589,12601,12611,12613,12619,12637,12641, 
    12647,12653,12659,12671,12689,12697,12703,12713,12721,12739, 
    12743,12757,12763,12781,12791,12799,12809,12821,12823,12829, 
    12841,12853,12889,12893,12899,12907,12911,12917,12919,12923, 
    12941,12953,12959,12967,12973,12979,12983,13001,13003,13007, 
    13009,13033,13037,13043,13049,13063,13093,13099,13103,13109, 
    13121,13127,13147,13151,13159,13163,13171,13177,13183,13187, 
    13217,13219,13229,13241,13249,13259,13267,13291,13297,13309, 
    13313,13327,13331,13337,13339,13367,13381,13397,13399,13411, 
    13417,13421,13441,13451,13457,13463,13469,13477,13487,13499 };

  if ( n == -1 )
  {
    return PRIME_MAX;
  }
  else if ( n == 0 )
  {
    return 1;
  }
  else if ( n <= PRIME_MAX )
  {
    return npvec[n-1];
  }
  else
  {
    cout << "\n";
    cout << "PRIME - Fatal error!\n";
    cout << "  Unexpected input value of n = " << n << "\n";
    exit ( 1 );
  }

  return 0;
# undef PRIME_MAX
}
//******************************************************************************

unsigned long random_initialize ( unsigned long seed )

//******************************************************************************
//
//  Purpose:
//
//    RANDOM_INITIALIZE initializes the RANDOM random number generator.
//
//  Discussion:
//
//    If you don't initialize RANDOM, the random number generator, 
//    it will behave as though it were seeded with value 1.  
//    This routine will either take a user-specified seed, or
//    (if the user passes a 0) make up a "random" one.  In either
//    case, the seed is passed to SRANDOM (the appropriate routine 
//    to call when setting the seed for RANDOM).  The seed is also
//    returned to the user as the value of the function.
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, unsigned long SEED, is either 0, which means that the user
//    wants this routine to come up with a seed, or nonzero, in which
//    case the user has supplied the seed.
//
//    Output, unsigned long RANDOM_INITIALIZE, is the value of the seed
//    passed to SRANDOM, which is either the user's input value, or if
//    that was zero, the value selected by this routine.
//
{
# define DEBUG 0

  if ( seed != 0 )
  {
    if ( DEBUG )
    {
      cout << "\n";
      cout << "RANDOM_INITIALIZE\n";
      cout << "  Initialize RANDOM with user SEED = " << seed << "\n";
    }
  }
  else
  {
    seed = get_seed ( );
    if ( DEBUG )
    {
      cout << "\n";
      cout << "RANDOM_INITIALIZE\n";
      cout << "  Initialize RANDOM with arbitrary SEED = " << seed << "\n";
    }
  }
//
//  Now set the seed.
//
  srandom ( seed );

  return seed;
# undef DEBUG
}
//******************************************************************************

int s_len_trim ( char *s )

//******************************************************************************
//
//  Purpose:
//
//    S_LEN_TRIM returns the length of a string to the last nonblank.
//
//  Modified:
//
//    26 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, a pointer to a string.
//
//    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
//    If S_LEN_TRIM is 0, then the string is entirely blank.
//
{
  int n;
  char* t;

  n = strlen ( s );
  t = s + strlen ( s ) - 1;

  while ( 0 < n ) 
  {
    if ( *t != ' ' )
    {
      return n;
    }
    t--;
    n--;
  }

  return n;
}
//******************************************************************************

void scale_from_simplex01 ( int m, int n, double t[], double x[] )

//******************************************************************************
//
//  Purpose:
//
//    SCALE_FROM_SIMPLEX01 rescales data from a unit to non-unit simplex.
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Reuven Rubinstein,
//    Monte Carlo Optimization, Simulation, and Sensitivity 
//      of Queueing Networks,
//    Wiley, 1986.
//
//  Parameters:
//
//    Input, int M, the dimension of the space.
//
//    Input, int N, the number of points.
//
//    Input, double T[M*(M+1)], the coordinates of the M+1 points that
//    define the simplex.  T[0:M-1,0] corresponds to the origin,
//    and T[0:M-1,J] will be the image of the J-th unit coordinate vector.
//
//    Input/output, double X[M*N], the data to be modified.
//
{
  double *a;
  int i;
  int j;
  double *v;

  a = new double[m*m];
  v = new double[m];

  for ( j = 0; j < m; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = t[i+j*(m+1)] - t[i+0*(m+1)];
    }
  }

  for ( j = 0; j < n; j++ )
  {

    for ( i = 0; i < m; i++ )
    {
      v[i] = x[i+j*m];
    }

    for ( i = 0; i < m; i++ )
    {
      x[i+j*m] = t[i+0*(m+1)];
      for ( j = 0; j < n; j++ )
      {
        x[i+j*m] = x[i+j*m] + a[i+j*m] * v[j];
      }
    }

  }

  delete [] a;
  delete [] v;

  return;
}
//******************************************************************************

void scale_to_ball01 ( int m, int n, double x[] )

//******************************************************************************
//
//  Purpose:
//
//    SCALE_TO_BALL01 translates and rescales data to fit within the unit ball.
//
//  Discussion:
//
//    Completely arbitrary input data is given.
//
//    The average of the data is computed, and taken as the coordinates
//    of the center C of a sphere.  The radius R of that sphere is the
//    distance from the center to the furthest point in the data set.
//
//    Then each point is transformed to the ball of center 0 and radius
//    1 by subtracting C and dividing by R:
//
//      X(1:M,J) -> ( X(1:M,J) - C(1:M) ) / R
//
//  Modified:
//
//    17 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the dimension of the space.
//
//    Input, int N, the number of points.
//
//    Input/output, double X[M*N], the data to be modified.
//
{
  int i;
  int j;
  double r;
  double scale;
  double *xave;
//
//  Determine the center.
//
  xave = new double[m];

  for ( i = 0; i < m; i++ )
  {
    xave[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      xave[i] = xave[i] + x[i+j*m];
    }
    xave[i] = xave[i] / ( double ) n;
  }
//
//  Determine the maximum distance of any point from the center.
//
  for ( j = 0; j < n; j++ )
  {
    r = 0.0;
    for ( i = 0; i < m; i++ )
    {
      r = r + pow ( x[i+j*m] - xave[i], 2 );
    }
    if ( scale < r )
    {
      scale = r ;
    }
  }

  scale = sqrt ( scale );


  if ( 0.0 < scale )
  {
    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < m; i++)
      {
        x[i+j*m] = ( x[i+j*m] - xave[i] ) / scale;
      }
    }
  }
  else
  {
    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < m; i++)
      {
        x[i+j*m] = 0.0;
      }
    }
  }

  delete [] xave;

  return;
}
//******************************************************************************

void scale_to_block01 ( int m, int n, double x[] )

//******************************************************************************
//
//  Purpose:
//
//    SCALE_TO_BLOCK01 translates and rescales data to fit in the unit block.
//
//  Discussion:
//
//    The minimum and maximum coordinate values M1(I) and M2(I) are
//    determined, and the maximum of M2(I) - M1(I) is used to scale
//    all the coordinates by the same factor.
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the dimension of the space.
//
//    Input, int N, the number of points.
//
//    Input/output, double X[M*N], the data to be modified.
//
{
  int i;
  int j;
  double *xmax;
  double *xmin;
  double xrange;
  double xrange2;

  xmax = new double[m];
  xmin = new double[m];
//
//  Determine the extremes in each dimension.
//
  xrange = 0.0;
  for ( i = 0; i < m; i++ )
  {
    xmin[i] = x[i+0*m];
    xmax[i] = x[i+0*m];
    for ( j = 1; j < n; j++ )
    {
      xmin[i] = d_min ( xmin[i], x[i+j*m] );
      xmax[i] = d_max ( xmax[i], x[i+j*m] );
    }
    xrange = d_max ( xrange, xmax[i] - xmin[i] );
  }
//
//  Extend all the extremes so that the range is the same in each dimension.
//
  for ( i = 0; i < m; i++ )
  {
    xrange2 = xrange - ( xmax[i] - xmin[i] );
    xmax[i] = xmax[i] + 0.5 * xrange2;
    xmin[i] = xmin[i] - 0.5 * xrange2;
  }
//
//  Now map the data to [0,1], using a single dilation factor 
//  for all dimensions.
//
  if ( 0.0 == xrange )
  {
    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        x[i+j*m] = 0.5;
      }
    }
  }
  else
  {
    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        x[i+j*m] = ( x[i+j*m] - xmin[i] ) / xrange;
      }
    }
  }

  delete [] xmax;
  delete [] xmin;

  return;
}
//******************************************************************************

void scale_to_cube01 ( int m, int n, double x[] )

//******************************************************************************
//
//  Purpose:
//
//    SCALE_TO_CUBE01 translates and rescales data to the unit hypercube.
//
//  Discussion:
//
//    In each coordinate dimension I, the minimum and maximum coordinate
//    values M1(I) and M2(I) are determined.
//
//    Then, in each coordinate, the points are rescaled as
//
//      X(I) -> ( X(I) - M1(I) ) / ( M2(I) - M1(I) ).
//
//  Modified:
//
//    17 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the dimension of the space.
//
//    Input, int N, the number of points.
//
//    Input/output, double X[M*N], the data to be modified.
//
{
  int i;
  int j;
  double xmax;
  double xmin;

  for ( i = 0; i < m; i++ )
  {
    xmin = x[i+0*m];
    xmax = x[i+0*m];
    for ( j = 1; j < n; j++ )
    {
      if ( x[i+j*m] < xmin )
      {
        xmin = x[i+j*m];
      }
      if ( xmax < x[i+j*m] )
      {
        xmax = x[i+j*m];
      }
    }

    if ( 0.0 < xmax - xmin )
    {
      for ( j = 0; j < n; j++ )
      {
        x[i+j*m] = ( x[i+j*m] - xmin ) / ( xmax - xmin );
      }
    }
    else
    {
      for ( j = 0; j < n; j++ )
      {
        x[i+j*m] = 0.0;
      }
    }
  }

  return;
}
//**********************************************************************

void timestamp ( void )

//**********************************************************************
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    May 31 2001 09:45:54 AM
//
//  Modified:
//
//    24 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
#define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
#undef TIME_SIZE
}
//**********************************************************************

char *timestring ( void )

//**********************************************************************
//
//  Purpose:
//
//    TIMESTRING returns the current YMDHMS date as a string.
//
//  Example:
//
//    May 31 2001 09:45:54 AM
//
//  Modified:
//
//    24 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, char *TIMESTRING, a string containing the current YMDHMS date.
//
{
#define TIME_SIZE 40

  const struct tm *tm;
  size_t len;
  time_t now;
  char *s;

  now = time ( NULL );
  tm = localtime ( &now );

  s = new char[TIME_SIZE];

  len = strftime ( s, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  return s;
#undef TIME_SIZE
}
//******************************************************************************

double triangle_area_2d ( double v1[2], double v2[2], double v3[2] )

//******************************************************************************
//
//  Purpose:
//
//    TRIANGLE_AREA_2D computes the area of a triangle in 2D.
//
//  Modified:
//
//    17 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double V1[2], V2[2], V3[2], the triangle vertices.
//
//    Output, double TRIANGLE_AREA_2D, the absolute area of the triangle.
//
{
  double area;

  area = 0.5 * fabs ( 
    ( v1[0] * ( v2[1] - v3[1] ) 
    + v2[0] * ( v3[1] - v1[1] ) 
    + v3[0] * ( v1[1] - v2[1] ) ) );

  return area;
}
//******************************************************************************

void tuple_next_fast ( int m, int n, int rank, int x[] )

//******************************************************************************
//
//  Purpose:
//
//    TUPLE_NEXT_FAST computes the next element of a tuple space, "fast".
//
//  Discussion:
//
//    The elements are N vectors.  Each entry is constrained to lie
//    between 1 and M.  The elements are produced one at a time.
//    The first element is
//      (1,1,...,1)
//    and the last element is
//      (M,M,...,M)
//    Intermediate elements are produced in lexicographic order.
//
//  Example:
//
//    N = 2,
//    M = 3
//
//    INPUT        OUTPUT
//    -------      -------
//    Rank          X
//    ----          ----
//   -1            -1 -1
//
//    0             1  1
//    1             1  2
//    2             1  3
//    3             2  1
//    4             2  2
//    5             2  3
//    6             3  1
//    7             3  2
//    8             3  3
//    9             1  1
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the maximum entry in each component.
//    M must be greater than 0.
//
//    Input, int N, the number of components.
//    N must be greater than 0.
//
//    Input, integer RANK, indicates the rank of the tuples.
//    Typically, 0 <= RANK < N**M; values larger than this are legal
//    and meaningful, and are equivalent to the corresponding value
//    MOD N**M.  If RANK < 0, this indicates that this is the first
//    call for the given values of (M,N).  Initialization is done,
//    and X is set to a dummy value.
//
//    Output, int X[N], the next tuple, or a dummy value if initialization
//    is being done.
//
{
  static int *base = NULL;
  int i;
//
  if ( rank < 0 )
  {
    if ( m <= 0 )
    {
      cout << "\n";
      cout << "TUPLE_NEXT_FAST - Fatal error!\n";
      cout << "  The value M <= 0 is not legal.\n";
      cout << "  M = " << m << "\n";
      exit ( 1 );
    }
    if ( n <= 0 )
    {
      cout << "\n";
      cout << "TUPLE_NEXT_FAST - Fatal error!\n";
      cout << "  The value N <= 0 is not legal.\n";
      cout << "  N = " << n << "\n";
      exit ( 1 );
    }

    if ( base )
    {
      delete [] base;
    }
    base = new int[n];

    base[n-1] = 1;
    for ( i = n-2; 0 <= i; i-- )
    {
      base[i] = base[i+1] * m;
    }
    for ( i = 0; i < n; i++ )
    {
      x[i] = -1;
    }
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      x[i] = ( ( rank / base[i] ) % m ) + 1;
    }
  }
  return;
}
//******************************************************************************

double *uniform_in_annulus01_accept ( int m, int n, double r, int *seed )

//******************************************************************************
//
//  Purpose:
//
//    UNIFORM_IN_ANNULUS01_ACCEPT accepts uniform points in a unit annulus.
//
//  Discussion:
//
//    The acceptance/rejection method is used.
//
//    The annulus is the area between two concentric spheres.  The
//    inner sphere has radius R, and the outer sphere has radius 1.
//
//  Modified:
//
//    17 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the dimension of the space.
//
//    Input, int N, the number of points.
//
//    Input, double R, the inner radius of the annulus, which must
//    be less than 1.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double UNIFORM_IN_ANNULUS01_ACCEPT[M*N], the points.
//
{
  int i;
  int j;
  double norm;
  double *x;
//
  if ( 1.0 <= r * r )
  {
    cout << "\n";
    cout << "UNIFORM_IN_ANNULUS01_ACCEPT - Fatal error!\n";
    cout << "  1 <= R * R.\n";
    cout << "  R = " << r << "\n";
    return NULL;
  }

  x = new double[m*n];
//
//  Generate points in the square.
//  Accept points that are inside the unit circle and outside the circle of radius R.
//
  for ( j = 0; j < n; j++ )
  {
    for ( ; ; )
    {
      norm = 0.0;
      for ( i = 0; i < m; i++ )
      {
        x[i+j*m] = d_uniform_01 ( seed );
        x[i+j*m] = 2.0 * x[i+j*m] - 1.0;
        norm = norm + pow ( x[i+j*m], 2 );
      }
      if ( r * r <= norm && norm <= 1.0E+00 )
      {
        break;
      }
    }
  }

  return x;
}
//******************************************************************************

double *uniform_in_circle01_map ( int m, int n, int *seed )

//******************************************************************************
//
//  Purpose:
//
//    UNIFORM_IN_CIRCLE01_MAP maps uniform points into the unit circle.
//
//  Discussion:
//
//    The unit circle is centered at the origin and has radius 1.
//
//  Modified:
//
//    24 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the dimension of the space.
//
//    Input, int N, the number of points.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double UNIFORM_IN_CIRCLE01_MAP[M*N], the points.
//
{
# define PI 3.141592653589793

  int j;
  double *r;
  double *t;
  double *x;
//
  r = new double[n];
  t = new double[n];
  x = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    r[j] = d_uniform_01 ( seed );
    r[j] = sqrt ( r[j] );
  }

  for ( j = 0; j < n; j++ )
  {
    t[j] = 2.0 * PI * d_uniform_01 ( seed );
  }

  for ( j = 0; j < n; j++ )
  {
    x[0+j*m] = r[j] * cos ( t[j] );
    x[1+j*m] = r[j] * sin ( t[j] );
  }

  delete [] r;
  delete [] t;

  return x;
# undef PI
}
//******************************************************************************

double *uniform_in_cube01 ( int m, int n, int *seed )

//******************************************************************************
//
//  Purpose:
//
//    UNIFORM_IN_CUBE01 creates uniform points in the unit hypercube.
//
//  Discussion:
//
//    The unit hypercube is defined as points all of whose components are between
//    0 and 1.
//
//  Modified:
//
//    17 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the dimension of the space.
//
//    Input, int N, the number of points.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double UNIFORM_IN_CUBE01[M*N], the points.
//
{
  int i;
  int j;
  double *x;

  x = new double[m*n];

  dvec_uniform_01 ( m*n, seed, x );

  return x;
}
//******************************************************************************

double *uniform_in_ellipsoid_map ( int m, int n, double a[], double r, int *seed )

//******************************************************************************
//
//  Purpose:
//
//    UNIFORM_IN_ELLIPSOID_MAP maps uniform points into an ellipsoid.
//
//  Discussion:
//
//    The points X in the ellipsoid are described by an M by M positive
//    definite symmetric matrix A, and a "radius" R, such that
//
//      X' * A * X <= R * R
//
//  Modified:
//
//    20 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Reuven Rubinstein,
//    Monte Carlo Optimization, Simulation, and Sensitivity
//      of Queueing Networks,
//    Wiley, 1986.
//
//  Parameters:
//
//    Input, int M, the dimension of the space.
//
//    Input, int N, the number of points.
//
//    Input, double A[M*M], the matrix that describes the ellipsoid.
//
//    Input, double R, the right hand side of the ellipsoid equation.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double UNIFORM_IN_ELLIPSOID_MAP[M*N], the points.
//
{
  int i;
  int j;
  int k;
  double *u;
  double *x;
  double *y;
//
//  Get the upper triangular Cholesky factor U of A.
//
  u = dpo_fa ( m, a );

  if ( !u )
  {
    cout << "\n";
    cout << "UNIFORM_IN_ELLIPSOID_MAP - Fatal error!\n";
    cout << "  SPO_FA reports that the matrix A\n";
    cout << "  is not positive definite symmetric.\n";
    exit ( 1 );
  }

  y = uniform_in_sphere01_map ( m, n, seed );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      y[i+j*m] = r * y[i+j*m];
    }
  }

  x = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      x[i+j*m] = 0.0;
      for ( k = 0; k <= i; k++ )
      {
        x[i+j*m] = x[i+j*m] + u[k+i*m] * y[k+j*m];
      }
    }
  }

  delete [] u;
  delete [] y;

  return x;
}
//******************************************************************************

double *uniform_in_parallelogram_map ( double v1[2], double v2[2], 
  double v3[2], int n, int *seed )

//******************************************************************************
//
//  Purpose:
//
//    UNIFORM_IN_PARALLELOGRAM_MAP maps uniform points into a parallelogram.
//
//  Discussion:
//
//    The parallelogram is defined by three vertices, V1, V2 and V3.
//    The missing vertex V4 is equal to V2+V3-V1.
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Greg Turk,
//    Generating Random Points in a Triangle,
//    in Graphics Gems,
//    edited by Andrew Glassner,
//    AP Professional, 1990, pages 24-28.
//
//  Parameters:
//
//    Input, double V1[2], V2[2], V3[2], the vertices.
//
//    Input, int N, the number of points.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double UNIFORM_IN_PARALLELOGRAM_MAP[2*N], the points.
//
{
# define M 2

  int i;
  int j;
  double r;
  double s;
  double *x;

  x = new double[M*n];

  for ( j = 0; j < n; j++ )
  {
    r = d_uniform_01 ( seed );
    s = d_uniform_01 ( seed );

    for ( i = 0; i < M; i++ )
    {
      x[i+j*M] = ( 1.0 - r - s ) * v1[i] 
                       + r       * v2[i] 
                           + s   * v3[i];
    }
  }

  return x;
# undef M
}
//******************************************************************************

double *uniform_in_polygon_map ( int nv, double v[], int n, int *seed )

//******************************************************************************
//
//  Purpose:
//
//    UNIFORM_IN_POLYGON_MAP maps uniform points into a polygon.
//
//  Discussion:
//
//    If the polygon is regular, or convex, or at least star-shaped,
//    this routine will work.
//
//    This routine assumes that all points between the centroid and
//    any point on the boundary lie within the polygon.
//
//  Modified:
//
//    19 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NV, the number of vertices.
//
//    Input, double V[2*NV], the vertices of the polygon, listed in
//    clockwise or counterclockwise order.
//
//    Input, int N, the number of points.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double UNIFORM_IN_POLYGON_MAP[2*N], the points.
//
{
# define M 2

  double *area;
  double *centroid;
  int i;
  int ip1;
  int j;
  int k;
  double r[M];
  int t;
  int tp1;
  double total;
  double u;
  double *x;

  area = new double[nv];
  x = new double[M*n];
//
//  Find the centroid.
//
  centroid = polygon_centroid_2d ( nv, v );
//
//  Determine the areas of each triangle.
//
  total = 0.0;
  for ( i = 0; i < nv; i++ )
  {
    if ( i < nv-1 )
    {
      ip1 = i + 1;
    }
    else
    {
      ip1 = 0;
    }
    
    area[i] = triangle_area_2d ( v+2*i, v+2*ip1, centroid );
    total = total + area[i];
  }
//
//  Normalize the areas.
//
  for ( i = 0; i < nv; i++ )
  {
    area[i] = area[i] / total;
  }
//
//  Replace each area by the sum of itself and all previous ones.
//
  for ( i = 1; i < nv; i++ )
  {
    area[i] = area[i] + area[i-1];
  }

  for ( j = 0; j < n; j++ )
  {
//
//  Choose a triangle T at random, based on areas.
//
    u = d_uniform_01 ( seed );

    for ( k = 0; k < nv; k++ )
    {
      if ( u <= area[k] )
      {
        t = k;
        break;
      }
    }
//
//  Now choose a point at random in the triangle.
//
    if ( t < nv-1 )
    {
      tp1 = t + 1;
    }
    else
    {
      tp1 = 0;
    }

    dvec_uniform_01 ( M, seed, r );

    if ( 1.0 < r[0] + r[1] )
    {
      r[0] = 1.0 - r[0];
      r[1] = 1.0 - r[1];
    }

    for ( i = 0; i < M; i++ )
    {
      x[i+j*M] = ( 1.0 - r[0] - r[1] ) * v[i+M*t] 
                 +       r[0]          * v[i+M*tp1] 
                 +              r[1]   * centroid[i];
    }

  }

  delete [] area;
  delete [] centroid;

  return x;
# undef M
}
//******************************************************************************

double *uniform_in_sector_map ( double r1, double r2, double t1, 
  double t2, int n, int *seed )

//******************************************************************************
//
//  Purpose:
//
//    UNIFORM_IN_SECTOR_MAP maps uniform points into a circular sector.
//
//  Discussion:
//
//    The sector lies between circles with center at 0 and radius R1 and R2,
//    and between rays from the center at the angles T1 and T2.
//
//  Modified:
//
//    20 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Shirley,
//    Nonuniform Random Point Sets Via Warping,
//    Graphics Gems III,
//    Edited by David Kirk,
//    Academic Press, 1992.
//
//  Parameters:
//
//    Input, double R1, R2, the two radii.
//
//    Input, double T1, T2, the two angles, which should
//    be measured in radians, with T1 < T2.
//
//    Input, int N, the number of points.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double UNIFORM_IN_SECTOR_MAP[2*N], the points.
//
{
  int j;
  double *r;
  double *t;
  double *u;
  double *v;
  double *x;

  x = new double[2*n];
  r = new double[n];
  t = new double[n];
  u = new double[n];
  v = new double[n];

  dvec_uniform_01 ( n, seed, u );
  dvec_uniform_01 ( n, seed, v );

  for ( j = 0; j < n; j++ )
  {
    t[j] = ( 1.0 - u[j] ) * t1 + u[j] * t2;
    r[j] = sqrt ( ( 1.0 - v[j] ) * r1 * r1 + v[j] * r2 * r2 );

    x[0+j*2] = r[j] * cos ( t[j] );
    x[1+j*2] = r[j] * sin ( t[j] );
  }

  delete [] r;
  delete [] t;
  delete [] u;
  delete [] v;

  return x;
}
//******************************************************************************

double *uniform_in_simplex01_map ( int m, int n, int *seed )

//******************************************************************************
//
//  Purpose:
//
//    UNIFORM_IN_SIMPLEX01 maps uniform points into the unit simplex.
//
//  Discussion:
//
//    The interior of the unit M-dimensional simplex is the set of points X(1:M)
//    such that each X(I) is nonnegative, and sum(X(1:M)) <= 1.
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Reuven Rubinstein,
//    Monte Carlo Optimization, Simulation, and Sensitivity 
//      of Queueing Networks,
//    Wiley, 1986.
//
//  Parameters:
//
//    Input, int M, the dimension of the space.
//
//    Input, int N, the number of points.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double UNIFORM_IN_SIMPLEX01_MAP[M*N], the points.
//
{
  double *e;
  int i;
  int j;
  double total;
  double *x;
//
//  The construction begins by sampling M+1 points from the
//  exponential distribution with parameter 1.
//
  e = new double[m+1];
  x = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    dvec_uniform_01 ( m+1, seed, e );

    for ( i = 0; i <= m; i++ )
    {
      e[i] = -log ( e[i] );
    }

    total = 0.0;
    for ( i = 0; i <= m; i++ )
    {
      total = total + e[i];
    }

    for ( i = 0; i < m; i++ )
    {
      x[i+m*j] = e[i] / total;
    }

  }

  delete [] e;

  return x;
}
//******************************************************************************

double *uniform_in_sphere01_map ( int m, int n, int *seed )

//******************************************************************************
//
//  Purpose:
//
//    UNIFORM_IN_SPHERE01_MAP maps uniform points into the unit sphere.
//
//  Discussion:
//
//    The sphere has center 0 and radius 1.
//
//    We first generate a point ON the sphere, and then distribute it
//    IN the sphere.
//
//  Modified:
//
//    17 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Russell Cheng,
//    Random Variate Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998, pages 168.
//
//    Reuven Rubinstein,
//    Monte Carlo Optimization, Simulation, and 
//    Sensitivity of Queueing Networks,
//    Wiley, 1986, page 232.
//
//  Parameters:
//
//    Input, int M, the dimension of the space.
//
//    Input, int N, the number of points.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double X[M*N], the points.
//
{
  double exponent;
  int i;
  int j;
  double norm;
  double r;
  double *v;
  double *x;
//
  exponent = 1.0E+00 / ( double ) ( m );

  v = new double[m];
  x = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
//
//  Fill a vector with normally distributed values.
//
    dvec_normal_01 ( m, seed, v );
//
//  Compute the length of the vector.
//
    norm = 0.0;
    for ( i = 0; i < m; i++ )
    {
      norm = norm + pow ( v[i], 2 );
    }
    norm = sqrt ( norm );
//
//  Normalize the vector.
//
    for ( i = 0; i < m; i++ )
    {
      v[i] = v[i] / norm;
    }
//
//  Now compute a value to map the point ON the sphere INTO the sphere.
//
    r = d_uniform_01 ( seed );
    r = pow ( r, exponent );

    for ( i = 0; i < m; i++ )
    {
      x[i+j*m] = r * v[i];
    }
  }

  delete [] v;

  return x;
}
//******************************************************************************

double *uniform_in_triangle_map1 ( double v1[2], double v2[2], double v3[2],
  int n, int *seed )

//******************************************************************************
//
//  Purpose:
//
//    UNIFORM_IN_TRIANGLE_MAP1 maps uniform points into a triangle.
//
//  Discussion:
//
//    This routine uses Turk's rule 1.
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Greg Turk,
//    Generating Random Points in a Triangle,
//    in Graphics Gems,
//    edited by Andrew Glassner,
//    AP Professional, 1990, pages 24-28.
//
//  Parameters:
//
//    Input, double V1[2], V2[2], V3[2], the vertices.
//
//    Input, int N, the number of points.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double UNIFORM_IN_TRIANGLE_MAP1[2*N], the points.
//
{
# define M 2

  double a;
  double b;
  double c;
  int i;
  int j;
  double r[M];
  double total;
  double *x;

  x = new double[M*n];

  for ( j = 0; j < n; j++ )
  {
    dvec_uniform_01 ( M, seed, r );

    r[1] = sqrt ( r[1] );

    a = 1.0 - r[1];
    b = ( 1.0E+00 - r[0] ) * r[1];
    c = r[0] * r[1];

    for ( i = 0; i < M; i++ )
    {
      x[i+j*M] = a * v1[i] + b * v2[i] + c * v3[i];
    }
  }

  return x;
# undef M
}
//******************************************************************************

double *uniform_in_triangle_map2 ( double v1[2], double v2[2], double v3[2],
  int n, int *seed )

//******************************************************************************
//
//  Purpose:
//
//    UNIFORM_IN_TRIANGLE_MAP2 maps uniform points into a triangle.
//
//  Discussion:
//
//    The triangle is defined by three vertices.
//
//    This routine uses Turk's rule 2.
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Greg Turk,
//    Generating Random Points in a Triangle,
//    in Graphics Gems,
//    edited by Andrew Glassner,
//    AP Professional, 1990, pages 24-28.
//
//  Parameters:
//
//    Input, double V1[2], V2[2], V3[2], the vertices.
//
//    Input, int N, the number of points.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double UNIFORM_IN_TRIANGLE_MAP2[2*N], the points.
//
{
# define M 2

  int i;
  int j;
  double r[M];
  double total;
  double *x;
//
  x = new double[M*n];

  for ( j = 0; j < n; j++ )
  {
    dvec_uniform_01 ( M, seed, r );

    if ( 1.0 < r[0] + r[1] )
    {
      r[0] = 1.0 - r[0];
      r[1] = 1.0 - r[1];
    }

    for ( i = 0; i < M; i++ )
    {
      x[i+j*M] = ( 1.0E+00 - r[0] - r[1] ) * v1[i] 
               +             r[0]          * v2[i] 
               +                    r[1]   * v3[i];
    }
  }

  return x;
# undef M
}
//******************************************************************************

double *uniform_in_triangle01_map ( int n, int *seed )

//******************************************************************************
//
//  Purpose:
//
//    UNIFORM_IN_TRIANGLE01_MAP maps uniform points into the unit triangle.
//
//  Discussion:
//
//    The triangle is defined by the three vertices (1,0), (0,1) and (0,0).
//    Because this is a right triangle, it is easy to generate sample points.
//    In the case of a general triangle, more care must be taken.
//
//  Modified:
//
//    19 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input/output, int *SEED, a seed for the random number geneator.
//
//    Output, double UNIFORM_IN_TRIANGLE01_MAP[2*N], the points.
//
{
# define M 2

  int i;
  int j;
  double r[M];
  double total;
  double *x;
//
  x = new double[M*n];

  for ( j = 0; j < n; j++ )
  {
    dvec_uniform_01 ( M, seed, r );

    total = 0.0;
    for ( i = 0; i < M; i++ )
    {
      total = total + r[i];
    }
    if ( 1.0 < total )
    {
      for ( i = 0; i < M; i++ )
      {
        r[i] = 1.0 - r[i];
      }
    }

    for ( i = 0; i < M; i++ )
    {
      x[i+j*M] = r[i];
    }
  }

  return x;
# undef M
}
//******************************************************************************

double *uniform_on_ellipsoid_map ( int m, int n, double a[], 
  double r, int *seed )

//******************************************************************************
//
//  Purpose:
//
//    UNIFORM_ON_ELLIPSOID_MAP maps uniform points onto an ellipsoid.
//
//  Discussion:
//
//    The points X on the ellipsoid are described by an M by M positive
//    definite symmetric matrix A, and a "radius" R, such that
//
//      X' * A * X = R * R
//
//  Modified:
//
//    20 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Reuven Rubinstein,
//    Monte Carlo Optimization, Simulation, and Sensitivity 
//    of Queueing Networks,
//    Wiley, 1986.
//
//  Parameters:
//
//    Input, int M, the dimension of the space.
//
//    Input, int N, the number of points.
//
//    Input, double A[M*M], the matrix that describes the ellipsoid.
//
//    Input, double R, the right hand side of the ellipsoid equation.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double UNIFORM_ON_ELLIPSOID_MAP[M*N], the points.
//
{
  int i;
  int j;
  int k;
  double *u;
  double *v;
  double *x;
//
//  Get the factor U.
//
  u = dpo_fa ( m, a );

  if ( !u )
  {
    cout << "\n";
    cout << "UNIFORM_ON_ELLIPSOID_MAP - Fatal error!\n";
    cout << "  DPO_FA reports that the matrix A \n";
    cout << "  is not positive definite symmetric.\n";
    exit ( 1 );
  }

  v = new double[m];
  x = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    direction_random_nd ( m, seed, v );

    for ( i = 0; i < m; i++ )
    {
      v[i] = r * v[i];
    }
    for ( i = 0; i < m; i++ )
    {
      x[i+j*m] = 0.0;
      for ( k = 0; k <= i; k++ )
      {
        x[i+j*m] = x[i+j*m] + v[k] * u[k+i*m];
      }
    }
  }
  delete [] u;
  delete [] v;

  return x;
}
//******************************************************************************

double *uniform_on_simplex01_map ( int m, int n, int *seed )

//******************************************************************************
//
//  Purpose:
//
//    UNIFORM_ON_SIMPLEX01_MAP maps uniform points onto the unit simplex.
//
//  Discussion:
//
//    The surface of the unit M-dimensional simplex is the set of points
//    X(1:M) such that each X(I) is nonnegative,
//    every X(I) is no greater than 1, and
//
//    ( X(I) = 0 for some I, or sum ( X(1:M) ) = 1. )
//
//    In M dimensions, there are M sides, and one main face.
//    This code picks a point uniformly with respect to "area".
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Reuven Rubinstein,
//    Monte Carlo Optimization, Simulation, and Sensitivity 
//      of Queueing Networks,
//    Wiley, 1986.
//
//  Parameters:
//
//    Input, int M, the dimension of the space.
//
//    Input, int N, the number of points.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double UNIFORM_ON_SIMPLEX01_MAP[M*N], the points.
//
{
  double area1;
  double area2;
  double *e;
  int i;
  int j;
  double r;
  double total;
  double u;
  double *x;
//
//  The construction begins by sampling M points from the
//  exponential distribution with parameter 1.
//
  e = new double[m];
  x = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    dvec_uniform_01 ( m, seed, e );

    for ( i = 0; i < m; i++ )
    {
      e[i] = -log ( e[i] );
    }

    total = 0.0;
    for ( i = 0; i < m; i++ )
    {
      total = total + e[i];
    }
//
//  Based on their relative areas, choose a side of the simplex,
//  or the main face.
//
    for ( i = 0; i < m; i++ )
    {
      x[i+j*m] = e[i] / total;
    }

    area1 = sqrt ( ( double ) m );
    area2 = ( double ) m;

    r = d_uniform_01 ( seed );

    if ( area1 / ( area1 + area2 ) < r )
    {
      i = i_random ( 0, m-1, seed );
      x[i+j*m] = 0.0;
    }

  }
  delete [] e;

  return x;
}
//******************************************************************************

double *uniform_on_sphere01_map ( int m, int n, int *seed )

//******************************************************************************
//
//  Purpose:
//
//    UNIFORM_ON_SPHERE01_MAP maps uniform points onto the unit sphere.
//
//  Discussion:
//
//    The sphere has center 0 and radius 1.
//
//  Modified:
//
//    19 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Russell Cheng,
//    Random Variate Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998, pages 168.
//
//    Reuven Rubinstein,
//    Monte Carlo Optimization, Simulation, and Sensitivity 
//    of Queueing Networks,
//    Wiley, 1986, page 234.
//
//  Parameters:
//
//    Input, int M, the dimension of the space.
//
//    Input, int N, the number of points.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double UNIFORM_ON_SPHERE01_MAP[M*N], the points.
//
{
  int i;
  int j;
  double norm;
  double *u;
  double *x;

  u = new double[m];
  x = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
//
//  Fill a vector with normally distributed values.
//
    dvec_normal_01 ( m, seed, u );
//
//  Compute the length of the vector.
//
    norm = 0.0;
    for ( i = 0; i < m; i++ )
    {
      norm = norm + u[i] * u[i];
    }
    norm = sqrt ( norm );
//
//  Normalize the vector.
//
    for ( i = 0; i < m; i++ )
    {
      x[i+j*m] = u[i] / norm;
    }

  }

  delete [] u;
  return x;
}
//******************************************************************************

double *uniform_walk ( int m, int n, int *seed )

//******************************************************************************
//
//  Purpose:
//
//    UNIFORM_WALK generates points on a uniform random walk.
//
//  Discussion:
//
//    The first point is at the origin.  Uniform random numbers are
//    generated to determine the direction of the next step, which
//    is always of length 1, and in coordinate direction.
//
//  Modified:
//
//    26 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the dimension of the space.
//
//    Input, int N, the number of points.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double UNIFORM_WALK[M*N], the points.
//
{
  double arg;
  double dir;
  int i;
  int j;
  double *x;

  x = new double[m*n];

  j = 0;
  for ( i = 0; i < m; i++ )
  {
    x[i+j*m] = 0.0;
  }

  for ( j = 1; j < n; j++ )
  {
    dir = d_uniform_01 ( seed );
    dir = ( double ) ( 2 * m ) * ( dir - 0.5 );

    for ( i = 0; i < m; i++ )
    {
      x[i+j*m] = x[i+(j-1)*m];
    }
    arg = fabs ( dir ) + 0.5;
    i = d_nint ( arg );
    i = i_min ( i, m );
    i = i_max ( i, 1 );
    i = i - 1;

    if ( dir < 0.0 )
    {
      x[i+j*m] = x[i+j*m] - 1.0;
    }
    else
    {
      x[i+j*m] = x[i+j*m] + 1.0;
    }

  }

  return x;
}
//******************************************************************************

void write_data ( int ndim, int n, double r[], char *file_out_name, 
  char *title )

//******************************************************************************
//
//  Purpose:
//
//    WRITE_DATA writes out a table of X, Y data.
//
//  Modified:
//
//    16 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NDIM, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double R[NDIM*N], the points.
//
//    Input, char *FILE_OUT_NAME, the name of the output file.
//
//    Input, character ( len = * ) TITLE, a title for the data.
//
{
  ofstream file_out;
  int i;
  int j;
  int mhi;
  int mlo;
  char *s;

  file_out.open ( file_out_name );

  if ( !file_out )
  {
    cout << "\n";
    cout << "WRITE_DATA - Fatal error!\n";
    cout << "  Could not open the output file.\n";
    exit ( 1 );
  }

  s = timestring ( );

  file_out << "#  " << file_out_name << "\n";
  file_out << "#  created by routine WRITE_DATA.CC" << "\n";
  file_out << "#  at " << s << "\n";
  file_out << "#\n";
  file_out << "#  Title: " << title << "\n";
  file_out << "#\n";
  file_out << "#  NDIM =   " << setw(12) << ndim   << "\n";
  file_out << "#  N =      " << setw(12) << n      << "\n";
  file_out << "#  EPSILON (unit roundoff) = " << d_epsilon ( ) << "\n";
  file_out << "#\n";

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < ndim; i++ )
    {
      file_out << setw(10) << r[i+j*ndim] << "  ";
    }
    file_out << "\n";
  }

  file_out.close ( );

  return;
}
