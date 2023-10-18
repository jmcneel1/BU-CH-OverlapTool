#ifndef _bueht_util
#define _bueht_util

/*
  James McNeely

  Some small utility functions used throughout the program.

  So far, this included a function to calculate the factorial
  and a binomial coefficient.
  
*/

namespace BUEHT
{

  int factorial ( int n )
  {
    int result = 1;
    for ( int j = 2; j <= n; j++ )
    {
      result = result * j;
    }
    return result;
  }

  double binomial_coeff ( int n, int k )
  {
    if ( ( 0 <= k ) && ( k < n ) )
    {
      return factorial(n)/(factorial(k)*factorial(n-k));
    }
    else
    {
      return 0;
    }
  }

}

#endif