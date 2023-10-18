#ifndef _bueht_atom
#define _bueht_atom

#include <stdlib.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <atomic_properties.h>
#include <bueht_constants.h>
#include <bueht_parameters.h>

namespace BUEHT
{

/*

  James McNeely

  This contains the Atom class

  The main private variables of an atom object are as follows:

  1: The atomic number.
  2: The atomic mass.
  3: The number of valence shells for the atom.
  4. The coordinates for the atom

  This class also contains features that are irrelevent for this simple tool, but
  were left in the class for transferrability sake.

  The constructors are listed below in the code, and are the easiest way to 
  construct an Atom object.

*/

class Atom
{
  public:

    Atom() = default; 

    Atom( const Atom & old_atom ) = default;

    /*
      There are no coordinates in the following constructor. Useful as a 
      placeholder.
    */
    
    Atom( int anum )                     
    {
      if ( anum >= 0 && anum < (NUM_ELEMENTS+1) )
      {
        myAtomicNumber = anum;
        myAtomicMass = bueht_atomic_masses[anum-1];
        myNumShells = (bueht_param_ptr[anum] - bueht_param_ptr[anum-1]);
      }
      else
      {
        std::cout << "Invalid Atomic Number ... Exiting.\n";
        exit(1);
      }
    }

    /*
      This constructor provides an object fully capable of being used
      to calculate the overlap integral
    */
  
    Atom( int anum, const std::vector<double> & coords ) : myX{coords[0]},
                                                           myY{coords[1]},
                                                           myZ{coords[2]}
    {
      Atom tmp(anum);
      this->myAtomicNumber = tmp.myAtomicNumber;
      this->myAtomicMass = tmp.myAtomicMass;
      this->myNumShells = tmp.myNumShells;
    }
  
    /*
      This constructor provides an object fully capable of being used
      to calculate the overlap integral
    */

    Atom( int anum, double x, double y, double z) : myX{x},
                                                    myY{y},
                                                    myZ{z}
    {
      Atom tmp(anum);
      this->myAtomicNumber = tmp.myAtomicNumber;
      this->myAtomicMass = tmp.myAtomicMass;
      this->myNumShells = tmp.myNumShells;
    }

    /*
      The following is useful for printing
    */

    std::string GetSymbol() const
    {
      return bueht_atomic_symbols[myAtomicNumber-1];
    }

    /*
      The following is useful for printing, although
      it doesn't really belong to the object. This might be 
      moved to a utility header in the future.
    */

    std::string GetSymbol( int anum )
    {
      Atom tat(anum);
      return tat.GetSymbol();
    }

    double GetX() const { return myX; }

    double GetY() const { return myY; }

    double GetZ() const { return myZ; }

    int GetNumShells () const { return myNumShells; }

    std::vector<double> GetCoordinates() const
    {
      std::vector<double> coords{myX,myY,myZ};
      return coords;
    }

    int GetAtomicNumber () const
    {
      return myAtomicNumber;
    }

    double GetAtomicMass () const
    {
      return myAtomicMass;
    }

    /*
      This is not functional in the current tool.
    */

    double GetMulliken () const
    {
      return myMulliken;
    }

    Atom & operator = ( const Atom & old_atom ) = default;

    bool operator == ( const Atom & rhs ) const
    {
      return myAtomicNumber == rhs.myAtomicNumber &&
             std::abs(myX-rhs.myX) < 1e-6 &&
             std::abs(myY-rhs.myY) < 1e-6 &&
             std::abs(myZ-rhs.myZ) < 1e-6 &&
             std::abs(myAtomicMass-rhs.myAtomicMass) < 1e-6;
    }

    bool operator != ( const Atom & rhs ) const
    {
      return myAtomicNumber != rhs.myAtomicNumber ||
             std::abs(myX-rhs.myX) >= 1e-6 ||
             std::abs(myY-rhs.myY) >= 1e-6 ||
             std::abs(myZ-rhs.myZ) >= 1e-6 ||
             std::abs(myAtomicMass-rhs.myAtomicMass) >= 1e-6;
    }

    friend std::ostream & operator << ( std::ostream & buffer, 
                                   const Atom & atom )
    {
      buffer << std::setw(3) << atom.GetSymbol() << std::setw(19) 
         << std::setprecision(14) << atom.myX*bueht_bohr_to_angstroem 
         << std::setw(19) << std::setprecision(14) << atom.myY*bueht_bohr_to_angstroem 
         << std::setw(19) << std::setprecision(14) << atom.myZ*bueht_bohr_to_angstroem;
      return buffer;
    }

    friend std::istream & operator >> ( std::istream & buffer, 
                                        Atom & atom )
    {
      std::string temp;
      double x, y, z;
      buffer >> temp;
      if ( temp.find_first_not_of("0123456789") == std::string::npos )
      {
        std::stringstream ss;
        ss >> temp;
        ss << atom.myAtomicNumber;
        atom.myAtomicMass = bueht_atomic_masses[atom.myAtomicNumber-1];
        atom.myNumShells = (bueht_param_ptr[atom.myAtomicNumber]-
                            bueht_param_ptr[atom.myAtomicNumber-1]);
      }
      else
      {
        atom.ChangeAtomType(temp);
      }
      buffer >> x >> y >> z;
      atom.myX = x/bueht_bohr_to_angstroem;
      atom.myY = y/bueht_bohr_to_angstroem;
      atom.myZ = z/bueht_bohr_to_angstroem;
      return buffer;
    }

    void ChangeAtomType ( const std::string & symbol )
    {
      bool found = false;
      int index = 0;
      while ( !found && index < NUM_ELEMENTS )
      {
        if ( symbol == bueht_atomic_symbols[index] )
        {
          found = true;
        }
        else
        {
          index++;
        }
      }
      if (found)
      {
        myAtomicNumber = index + 1;
        myAtomicMass = bueht_atomic_masses[index];
        myNumShells = bueht_param_ptr[myAtomicNumber] - bueht_param_ptr[myAtomicNumber-1];
      }
      else
      {
        std::cout << symbol << " is not an element! Exiting...\n";
        std::exit(1);
      }
    }

    void ChangeAtomType ( int anum )
    {
      if ( anum >= 0 && anum < (NUM_ELEMENTS+1) )
      {
        myAtomicNumber = anum;
        myAtomicMass = bueht_atomic_masses[anum-1];
        myNumShells = bueht_param_ptr[myAtomicNumber] - bueht_param_ptr[myAtomicNumber-1];
      }
    }

    void ChangeAtomicMass ( double new_mass )
    {
      myAtomicMass = new_mass;
    }

    void ChangeCoordinates( double x, double y, double z ) 
    {
      myX = x; myY = y; myZ = z;
    }

    void ChangeCoordinates( const std::vector<double> & new_coords )
    {
      if ( new_coords.size() != 3 )
      {
        std::cout << "Invalid coordinates! Exiting...";
        std::exit(1);
      }
      else
      {
        myX = new_coords[0];
        myY = new_coords[1];
        myZ = new_coords[2];
      }
    }

    void SetMulliken ( double mull_charge )
    {
      myMulliken = mull_charge;
    }

    ~Atom() = default;

  private:

    double myX;
    double myY;
    double myZ;
    double myAtomicMass;
    double myMulliken;
    int myAtomicNumber;
    int myNumShells;
    int myNumCore;

};

}

#endif
