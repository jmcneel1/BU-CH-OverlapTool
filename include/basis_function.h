#ifndef _bueht_basis_function
#define _bueht_basis_function

#include <iostream>
#include <stdlib.h>
#include <bueht_parameters.h>
#include <math.h>
#include <bueht_util.h>
#include <vector>

/*
  This BasisFunction class stores relevant information for a shell.
*/

namespace BUEHT
{

class BasisFunction
{
  public:

    BasisFunction () = default;

    BasisFunction ( int n, int l, int num_primitive,
                   const std::vector<double> & coeffs, 
                   const std::vector<double> & zeta )
    {
      if ( n > 0 && n < 5 ) // We are only going to Kr, with n=4 as the max
      {
        myN = n;
        if ( l < n )
        {
          myL = l;
          if ( (num_primitive == coeffs.size()) && (num_primitive == zeta.size()) )
          {
            this->myNumPrimitives = num_primitive;
            this->myZetas = zeta;
            this->myCoefficients = coeffs;
          }
          else
          {
            std::cout << "Incompatible primitize count and vector size! Exiting...\n";
            std::exit(1);
          }
        }
        else
        {
          std::cout << "Invalid Angular Momentum Quantum Number (l)! Exiting...\n";
          std::exit(1);
        }
      }
      else
      {
        std::cout<< "Invalid Principle Quantum Number (n)! Exiting...\n"; 
        std::exit(1);
      }
    }

    BasisFunction ( const int & n, const int & l, const int & atomic_num ) 
    {
      bool found = false;
      int count = 0;
      while ( ! found )
      {
        if ( (*(bueht_param_ptr[atomic_num-1]+count)).an == atomic_num )
        {
          if ( ((*(bueht_param_ptr[atomic_num-1]+count)).n == n) &&
               ((*(bueht_param_ptr[atomic_num-1]+count)).l == l) )
          {
            myN = n;
            myL = l;
            if ( std::abs((*(bueht_param_ptr[atomic_num-1]+count)).c2) < 1e-8 )
            {
              myNumPrimitives = 1;
              myZetas.push_back((*(bueht_param_ptr[atomic_num-1]+count)).z1);
              myCoefficients.push_back((*(bueht_param_ptr[atomic_num-1]+count)).c1);
            }
            else
            {
              myNumPrimitives = 2;
              myZetas.push_back((*(bueht_param_ptr[atomic_num-1]+count)).z1);
              myZetas.push_back((*(bueht_param_ptr[atomic_num-1]+count)).z2);
              myCoefficients.push_back((*(bueht_param_ptr[atomic_num-1]+count)).c1);
              myCoefficients.push_back((*(bueht_param_ptr[atomic_num-1]+count)).c2);
            }
            found = true;
          }
        }
        else
        {
          std::cout << "Unknown N/L for Atomic Number " << atomic_num << ". Exiting...\n";
          std::exit(1);
        }
        count++;
      }
    }

    double GetExponent ( const int & index ) const
    {
      if  ( (index >= 0) && (index < myNumPrimitives) )
      {
        return myZetas[index];
      }
      else
      {
        std::cout << "Invalid index for AO coefficient vector! Exiting...\n";
        std::exit(1);
      }
    }

    double GetCoefficient ( const int & index ) const
    {
      if  ( (index >= 0) && (index < myNumPrimitives) )
      {
        return myCoefficients[index];
      }
      else
      {
        std::cout << "Invalid index for AO coefficient vector! Exiting...\n";
        std::exit(1);
      }
    }

    std::vector<double> GetCoefficients () const
    {
      return myCoefficients;
    }

    std::vector<double> GetExponents () const
    {
      return myZetas;
    }

    void SetN (int n)
    {
      if ( n < 5 )
      {
        if ( std::abs(myL) < n )
        {
          myN = n;
        }
        else
        {
          std::cout << "Invalid Angular Momentum Quantum Number (l)! Exiting...\n";
          std::exit(1);
        }
      }
      else
      {
        std::cout<< "Invalid Principle Quantum Number (n)! Exiting...\n";
        std::exit(1);
      }
    }

    int GetN () const { return myN; }

    int GetL () const { return myL; }

    unsigned int GetNumPrimitives () const { return myNumPrimitives; }

    friend std::ostream & operator << ( std::ostream & buffer, 
                                   const BUEHT::BasisFunction & basis_function )
    {
      buffer << "N: " << basis_function.myN << "\n";
      buffer << "L: " << basis_function.myL << "\n";
      buffer << "Number of Primitives: " << basis_function.myNumPrimitives << "\n";
      for ( unsigned int i = 0; i < basis_function.myNumPrimitives; i++ )
      {
        buffer << "{" << basis_function.myCoefficients[i] << ", ";
        buffer << basis_function.myZetas[i] << "}";
        if ( i == ( basis_function.myNumPrimitives - 1 ) )
        {
          buffer << "\n";
        }
        else
        {
          buffer << ", ";
        }
      }
      return buffer;
    }

  private:

    int myN;
    int myL;
    unsigned int myNumPrimitives;
    std::vector<double> myCoefficients;
    std::vector<double> myZetas;

};

}

#endif
