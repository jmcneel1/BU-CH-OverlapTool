#ifndef _bueht_overlap
#define _bueht_overlap

#include "basis_function.h"
#include "bueht_constants.h"
#include "atom.h"
#include "rotate.h"
#include <math.h>
#include <vector>

namespace BUEHT
{
    /*

      James McNeely

      Here we use basis functions specific for EHT.
      Namely each BF is a contraction of MAX two STO primitives.
      Thus, there are 4 pairs to evaluate, but in many cases this 
      is drastically reduced due to there being a single STO primitive

      We also have the norms for each of the primitives stored already

      This is the kernel of the program.

      Here we define SS, SP, SD, PP, PD, and DD overlaps.

      Each "type" starts with prolate spheroidal coordinates
      to calculate the overlap (for each pair of primitives), 
      and proceeds to rotate the overlaps accordingly to
      get back to the input orientation.

      For convenience here, I provide a derivation of the PrimitiveOverlap
      function, where terms in the derivation are used in comments that follow

      ---------------------------------------------

      We start with the definition of a normalized spherical harmonic:
   
         

    */

    /*
      This defines the C coeffi
    */

    inline double Overlap_C_lmj ( const int & l, const int & m,
                           const int & j )
    {
      return std::pow(double(factorial(l-m))/factorial(l+m),0.5) *
             factorial(2.0*(l-j))*std::pow(-1.0,j)/
             (std::pow(2.0,l)*factorial(j)*factorial(l-j)*
             factorial(l-m-2*j));
    }

    double Overlap_A ( const int & n, const double & a )
    {
      double sum = 0.0;
      for ( unsigned int i = 0; i <= n; i++ )
      {
        
        sum+=std::pow(a,i)/double(factorial(i));
      }
      return std::exp(-a)*double(factorial(n))/std::pow(a,n+1)*sum;
    }

    double Overlap_B ( const int & n, const double & a )
    {
      double sum1 = 0.0;
      double sum2 = 0.0;
      if ( std::abs(a) < 1e-8 ) 
      {
        if ( (n%2) == 0 )
        {
          return 2.0/((double)(n+1));
        }
        else
        {
          return 0.0;
        }
      }
      for ( unsigned int i = 0; i <= n; i++ )
      {
        sum1+=std::pow(a,i)/factorial(i);
        sum2+=std::pow(-a,i)/factorial(i);
      }
      return factorial(n)/std::pow(a,n+1)*(-std::exp(-a)*sum1+std::exp(a)*sum2);
    }

    double PrimitiveOverlap (int na, int la, int ma, int nb, int lb, int mb,
                             double zeta1, double zeta2,
                             const Atom & atom1, const Atom & atom2 )
    {
      double distance = std::pow(std::pow(atom1.GetX()-atom2.GetX(),2)+
                                 std::pow(atom1.GetY()-atom2.GetY(),2)+
                                 std::pow(atom1.GetZ()-atom2.GetZ(),2),0.5)/bueht_bohr_to_angstroem;
      double sum1 = 0., termja, termjb, termka, termkb, termpa, termpb, rhoa, rhob;
      double termqa, termqb, termpre;
      if (ma != mb) return 0.;
      int ja_end, jb_end;
      rhoa = 1./2.*(zeta1+zeta2)*distance;
      rhob = 1./2.*(zeta1-zeta2)*distance;
      ja_end = (la-ma)/2;
      jb_end = (lb-mb)/2;
      termpre = std::pow((double)factorial(ma),2);
      for ( unsigned int ja = 0; ja <= ja_end; ja++ )
      {
        int pa_end = na-la+2*ja;
        int qa_end = la-ma-2*ja;
        termja = termpre;
        termja *= Overlap_C_lmj(la,ma,ja);
        termja *= (double)factorial(la-ma-2*ja);
        termja *= (double)factorial(na-la+2*ja);
        for ( unsigned int jb = 0; jb <= jb_end; jb++ )
        {
          int pb_end = nb-lb+2*jb;
          int qb_end = lb-mb-2*jb;
          termjb = termja;
          termjb *= Overlap_C_lmj(lb,mb,jb);
          termjb *= (double)factorial(lb-mb-2*jb);
          termjb *= (double)factorial(nb-lb+2*jb);
          for ( unsigned int ka = 0; ka <= ma; ka++ )
          {
            termka = termjb;
            termka /= (double)factorial(ma-ka);
            termka /= (double)factorial(ka);
            for ( unsigned int kb = 0; kb <= mb; kb++ )
            {
              termkb = termka;
              termkb /= (double)factorial(kb);
              termkb /= (double)factorial(mb-kb);
              for ( unsigned int pa = 0; pa <= pa_end; pa++ )
              {
                termpa = termkb;
                termpa /= (double)factorial(na-la+2*ja-pa);
                termpa /= (double)factorial(pa);
                for ( unsigned int pb = 0; pb <= pb_end; pb++ )
                {
                  termpb = termpa;
                  termpb /= (double)factorial(nb-lb+2*jb-pb);
                  termpb /= (double)factorial(pb);
                  for ( unsigned int qa = 0; qa <= qa_end; qa++ )
                  {
                    termqa = termpb;
                    termqa /= (double)factorial(la-ma-2*ja-qa);
                    termqa /= (double)factorial(qa);
                    for ( unsigned int qb = 0; qb <= qb_end; qb++ )
                    {
                      termqb = termqa;
                      termqb /= (double)factorial(lb-mb-2*jb-qb);
                      termqb /= (double)factorial(qb);
                      sum1+= termqb*
                           //std::pow(-1.0,-ka-kb+nb-pb)*
                           std::pow(-1.0,ka+kb+ma+pb+qb)*
                           //Overlap_B(2*kb+qa+qb+na-la+2*ja-pa+nb-lb+2*jb-pb,rhob)*
                           //Overlap_A(2.0*ka+qa+qb+pa+pb,rhoa); 
                           Overlap_B(2*ka+pa+pb+qa+qb,rhob)*
                           Overlap_A(2*kb+na-la+2*ja+nb-lb+2*jb-pa-pb+qa+qb,rhoa);
                    }
                  }
                }
              }
            }
          }
        }
      }
      return 0.5*std::pow(zeta1,na+0.5)*std::pow(zeta2,nb+0.5)*
             std::pow(1.0*(2*la+1)*(2*lb+1)/((double)factorial(2*na)*(double)factorial(2*nb)),0.5)*
             std::pow(distance,na+nb+1)*
             sum1;
    }

    double OverlapSS (const BasisFunction & bf1,
                    const BasisFunction & bf2,
                    const Atom & atom1,
                    const Atom & atom2 )
    {
      /*
        We assume here that all S basis functions are single-zeta.
      */
      if ( ( bf1.GetL() == 0 ) && ( bf2.GetL() == 0 ) )
      {
        return PrimitiveOverlap(bf1.GetN(),0,0,bf2.GetN(),0,0,
                                bf1.GetExponent(0), bf2.GetExponent(0),
                                atom1, atom2);
      }
      else
      {
        std::cout<<"ERROR! One of the BFs for SS Overlap isn't an S-Orbital!!!!\n";
        std::exit(1);
      }
    }   

    std::vector<double> OverlapSP ( const BasisFunction & bf1,
                       const BasisFunction & bf2,
                       const Atom & atom1,
                       const Atom & atom2 )
    {
      // We assume here that all S & P BFs are single-zeta
      if ( ( ( bf1.GetL() == 0 ) && ( bf2.GetL() == 1 ) ) || 
           ( ( bf1.GetL() == 1 ) && ( bf2.GetL() == 0 ) ) )
      {
        double distance, spz_overlap;
        double xvec, yvec, zvec;
        std::vector<double> ovlp {0.,0.,0.};
        std::vector<double> tcoor1 {0.,0.,0.};
        xvec = atom1.GetX() - atom2.GetX();
        yvec = atom1.GetY() - atom2.GetY();
        zvec = atom1.GetZ() - atom2.GetZ();
        distance = std::pow(xvec*xvec+yvec*yvec+zvec*zvec,0.5);
        std::vector<double> tcoor2 {0.,0.,distance};
        Atom tatom1(atom1.GetAtomicNumber(),tcoor1);
        Atom tatom2(atom2.GetAtomicNumber(),tcoor2);
        spz_overlap = PrimitiveOverlap(bf1.GetN(),bf1.GetL(),0,
                                       bf2.GetN(),bf2.GetL(),0,
                                       bf1.GetExponent(0),bf2.GetExponent(0),
                                       tatom1, tatom2 );
        ovlp[2] = ( bf2.GetL() == 0 ) ? -spz_overlap : spz_overlap;
        RotateSP(ovlp,atom1,atom2,distance);
        return ovlp;
      }
      else
      {
        std::cout<<"ERROR! One of the BFs for SP Overlap isn't correct!!!!\n";
        std::exit(1);
      }
    }   

    std::vector<double> OverlapPP ( const BasisFunction & bf1,
                                    const BasisFunction & bf2,
                                    const Atom & atom1,
                                    const Atom & atom2 )
    {
      // We assume here all P are single zeta
      // We return here a vector of length 9
      // The elements of the vector are as follows:
      // pxpx, pxpy, pxpz, pypx, pypy, pypz, pzpx, pzpy, pzpz
      // Thus ... Row Major

      // First we check to make sure both BFs are P-Shells
      if ( ( bf1.GetL() == 1 ) && ( bf2.GetL() == 1 ) )
      {
        double distance, pzpz_overlap, pxpx_overlap;
        double xvec, yvec, zvec;
        // Declare the empty overlap vector
        std::vector<double> ovlp {0.,0.,0.,0.,0.,0.,0.,0.,0.};
        // We'll translate the dimer to the origin...
        std::vector<double> tcoor1 {0.,0.,0.};
        // Calculate the relative orientation of the dimer
        xvec = atom1.GetX()-atom2.GetX();
        yvec = atom1.GetY() - atom2.GetY();
        zvec = atom1.GetZ() - atom2.GetZ();
        // Calculate the distance
        distance = std::pow(xvec*xvec+yvec*yvec+zvec*zvec,0.5);
        // Place the second atom on the Z-Axis and then 
        // declare atoms in the new coordinate frame
        std::vector<double> tcoor2 {0.,0.,distance};
        Atom tatom1(atom1.GetAtomicNumber(),tcoor1);
        Atom tatom2(atom2.GetAtomicNumber(),tcoor2);
        // Calculate the sigma and pi overlaps and place them in
        // The overlap matrix on the diagonal
        // RotatePP will use these values...
        pzpz_overlap = PrimitiveOverlap(bf1.GetN(),bf1.GetL(),0,
                                       bf2.GetN(),bf2.GetL(),0,
                                       bf1.GetExponent(0),bf2.GetExponent(0),
                                       tatom1, tatom2 );
        pxpx_overlap = PrimitiveOverlap(bf1.GetN(),bf1.GetL(),1,
                                       bf2.GetN(),bf2.GetL(),1,
                                       bf1.GetExponent(0),bf2.GetExponent(0),
                                       tatom1, tatom2 );
        ovlp[8] = -pzpz_overlap;
        ovlp[0] = ovlp[3] = pxpx_overlap;
        // Now perform the rotation and return the results
        RotatePP(ovlp,atom1,atom2,distance);
        return ovlp;
      }
      else
      {
        std::cout << "ERROR! One of the BFs for PP Overlap isn't correct!!!\n";
        std::exit(1);
      }
    }    

    std::vector<double> OverlapSD ( const BasisFunction & bf1,
                       const BasisFunction & bf2,
                       const Atom & atom1,
                       const Atom & atom2 )
    {
      // Starting with D-Orbitals, we have to check whether the BF
      // is single or double-zeta
      if ( ( ( bf1.GetL() == 0 ) && ( bf2.GetL() == 2 ) ) || 
           ( ( bf1.GetL() == 2 ) && ( bf2.GetL() == 0 ) ) )
      {
        double distance, sz2_overlap;
        double xvec, yvec, zvec;
        std::vector<double> ovlp {0.,0.,0.,0.,0.};
        std::vector<double> tcoor1 {0.,0.,0.};
        xvec = atom1.GetX() - atom2.GetX();
        yvec = atom1.GetY() - atom2.GetY();
        zvec = atom1.GetZ() - atom2.GetZ();
        distance = std::pow(xvec*xvec+yvec*yvec+zvec*zvec,0.5);
        std::vector<double> tcoor2 {0.,0.,distance};
        Atom tatom1(atom1.GetAtomicNumber(),tcoor1);
        Atom tatom2(atom2.GetAtomicNumber(),tcoor2);
        if ( bf1.GetL() == 2 )
        {
          if ( bf1.GetNumPrimitives() == 1 )
          {
            sz2_overlap = PrimitiveOverlap(bf1.GetN(),bf1.GetL(),0,
                                           bf2.GetN(),bf2.GetL(),0,
                                           bf1.GetExponent(0),bf2.GetExponent(0),
                                           tatom1, tatom2 );
          }
          else
          {
            sz2_overlap = bf1.GetCoefficient(0)*PrimitiveOverlap(bf1.GetN(),bf1.GetL(),0,
                                                                 bf2.GetN(),bf2.GetL(),0,
                                                                 bf1.GetExponent(0),
                                                                 bf2.GetExponent(0),
                                                                 tatom1, tatom2 ) +
                          bf1.GetCoefficient(1)*PrimitiveOverlap(bf1.GetN(),bf1.GetL(),0,
                                                                 bf2.GetN(),bf2.GetL(),0,
                                                                 bf1.GetExponent(1),
                                                                 bf2.GetExponent(0),
                                                                 tatom1, tatom2 );
          } 
        }                         
        else
        {
          if ( bf2.GetNumPrimitives() == 1 )
          {
            sz2_overlap = PrimitiveOverlap(bf1.GetN(),bf1.GetL(),0,
                                           bf2.GetN(),bf2.GetL(),0,
                                           bf1.GetExponent(0),bf2.GetExponent(0),
                                           tatom1, tatom2 );
          }
          else
          {
            sz2_overlap = bf2.GetCoefficient(0)*PrimitiveOverlap(bf1.GetN(),bf1.GetL(),0,
                                                                 bf2.GetN(),bf2.GetL(),0,
                                                                 bf1.GetExponent(0),
                                                                 bf2.GetExponent(0),
                                                                 tatom1, tatom2 ) +
                          bf2.GetCoefficient(1)*PrimitiveOverlap(bf1.GetN(),bf1.GetL(),0,
                                                                 bf2.GetN(),bf2.GetL(),0,
                                                                 bf1.GetExponent(0),
                                                                 bf2.GetExponent(1),
                                                                 tatom1, tatom2 );
          }
        }
        ovlp[2] = sz2_overlap;
        RotateSD(ovlp,atom1,atom2,distance);
        return ovlp;
      }
      else
      {
        std::cout<<"ERROR! One of the BFs for SD Overlap isn't correct!!!!\n";
        std::exit(1);
      }
    } 

    std::vector<double> OverlapPD ( const BasisFunction & bf1,
                       const BasisFunction & bf2,
                       const Atom & atom1,
                       const Atom & atom2 )
    {
      // Starting with D-Orbitals, we have to check whether the BF
      // is single or double-zeta
      if ( ( ( bf1.GetL() == 1 ) && ( bf2.GetL() == 2 ) ) || 
           ( ( bf1.GetL() == 2 ) && ( bf2.GetL() == 1 ) ) )
      {
        double distance, pzz2_overlap, pxxz_overlap;
        double xvec, yvec, zvec;
        std::vector<double> ovlp (15,0.0);
        std::vector<double> tcoor1 {0.,0.,0.};
        xvec = atom1.GetX() - atom2.GetX();
        yvec = atom1.GetY() - atom2.GetY();
        zvec = atom1.GetZ() - atom2.GetZ();
        distance = std::pow(xvec*xvec+yvec*yvec+zvec*zvec,0.5);
        std::vector<double> tcoor2 {0.,0.,distance};
        Atom tatom1(atom1.GetAtomicNumber(),tcoor1);
        Atom tatom2(atom2.GetAtomicNumber(),tcoor2);
        if ( bf1.GetL() == 2 )
        {
          if ( bf1.GetNumPrimitives() == 1 )
          {
            pzz2_overlap = PrimitiveOverlap(bf1.GetN(),bf1.GetL(),0,
                                           bf2.GetN(),bf2.GetL(),0,
                                           bf1.GetExponent(0),bf2.GetExponent(0),
                                           tatom1, tatom2 );
            pxxz_overlap = PrimitiveOverlap(bf1.GetN(),bf1.GetL(),1,
                                            bf2.GetN(),bf2.GetL(),1,
                                            bf1.GetExponent(0),bf1.GetExponent(0),
                                            tatom1,tatom2 );
          }
          else
          {
            pzz2_overlap = bf1.GetCoefficient(0)*PrimitiveOverlap(bf1.GetN(),bf1.GetL(),0,
                                                                 bf2.GetN(),bf2.GetL(),0,
                                                                 bf1.GetExponent(0),
                                                                 bf2.GetExponent(0),
                                                                 tatom1, tatom2 ) +
                          bf1.GetCoefficient(1)*PrimitiveOverlap(bf1.GetN(),bf1.GetL(),0,
                                                                 bf2.GetN(),bf2.GetL(),0,
                                                                 bf1.GetExponent(1),
                                                                 bf2.GetExponent(0),
                                                                 tatom1, tatom2 );
            pxxz_overlap = bf1.GetCoefficient(0)*PrimitiveOverlap(bf1.GetN(),bf1.GetL(),1,
                                                                  bf2.GetN(),bf2.GetL(),1,
                                                                  bf1.GetExponent(0),
                                                                  bf2.GetExponent(0),
                                                                  tatom1,tatom2 )+
                           bf1.GetCoefficient(1)*PrimitiveOverlap(bf1.GetN(),bf1.GetL(),1,
                                                                  bf2.GetN(),bf2.GetL(),1,
                                                                  bf1.GetExponent(1),
                                                                  bf2.GetExponent(0),
                                                                  tatom1,tatom2 );
          } 
        }                         
        else
        {
          if ( bf2.GetNumPrimitives() == 1 )
          {
            pzz2_overlap = PrimitiveOverlap(bf1.GetN(),bf1.GetL(),0,
                                           bf2.GetN(),bf2.GetL(),0,
                                           bf1.GetExponent(0),bf2.GetExponent(0),
                                           tatom1, tatom2 );
            pxxz_overlap = PrimitiveOverlap(bf1.GetN(),bf1.GetL(),1,
                                           bf2.GetN(),bf2.GetL(),1,
                                           bf1.GetExponent(0),bf2.GetExponent(0),
                                           tatom1, tatom2 );
          }
          else
          {
            pzz2_overlap = bf2.GetCoefficient(0)*PrimitiveOverlap(bf1.GetN(),bf1.GetL(),0,
                                                                 bf2.GetN(),bf2.GetL(),0,
                                                                 bf1.GetExponent(0),
                                                                 bf2.GetExponent(0),
                                                                 tatom1, tatom2 ) +
                          bf2.GetCoefficient(1)*PrimitiveOverlap(bf1.GetN(),bf1.GetL(),0,
                                                                 bf2.GetN(),bf2.GetL(),0,
                                                                 bf1.GetExponent(0),
                                                                 bf2.GetExponent(1),
                                                                 tatom1, tatom2 );
            pxxz_overlap = bf2.GetCoefficient(0)*PrimitiveOverlap(bf1.GetN(),bf1.GetL(),1,
                                                                 bf2.GetN(),bf2.GetL(),1,
                                                                 bf1.GetExponent(0),
                                                                 bf2.GetExponent(0),
                                                                 tatom1, tatom2 ) +
                          bf2.GetCoefficient(1)*PrimitiveOverlap(bf1.GetN(),bf1.GetL(),1,
                                                                 bf2.GetN(),bf2.GetL(),1,
                                                                 bf1.GetExponent(0),
                                                                 bf2.GetExponent(1),
                                                                 tatom1, tatom2 );
          }
        }
        ovlp[12] = (bf1.GetL() == 1) ? -pzz2_overlap : pzz2_overlap;
        ovlp[3] = (bf1.GetL() == 1) ? pxxz_overlap : -pxxz_overlap;
        RotatePD(ovlp,atom1,atom2,distance);
        return ovlp;
      }
      else
      {
        std::cout<<"ERROR! One of the BFs for PD Overlap isn't correct!!!!\n";
        std::exit(1);
      }
    } 

    std::vector<double> OverlapDD ( const BasisFunction & bf1,
                       const BasisFunction & bf2,
                       const Atom & atom1,
                       const Atom & atom2 )
    {
      // Starting with D-Orbitals, we have to check whether the BF
      // is single or double-zeta
      if ( ( ( bf1.GetL() == 2 ) && ( bf2.GetL() == 2 ) ) || 
           ( ( bf1.GetL() == 2 ) && ( bf2.GetL() == 2 ) ) )
      {
        double distance, z2z2_overlap, xzxz_overlap, xyxy_overlap;
        double xvec, yvec, zvec;
        std::vector<double> ovlp (25,0.0);
        std::vector<double> tcoor1 {0.,0.,0.};
        xvec = atom1.GetX() - atom2.GetX();
        yvec = atom1.GetY() - atom2.GetY();
        zvec = atom1.GetZ() - atom2.GetZ();
        distance = std::pow(xvec*xvec+yvec*yvec+zvec*zvec,0.5);
        std::vector<double> tcoor2 {0.,0.,distance};
        Atom tatom1(atom1.GetAtomicNumber(),tcoor1);
        Atom tatom2(atom2.GetAtomicNumber(),tcoor2);
        if ( bf1.GetNumPrimitives() == 1 )
        {
          if ( bf2.GetNumPrimitives() == 1 )
          {
            z2z2_overlap = PrimitiveOverlap(bf1.GetN(),bf1.GetL(),0,
                                          bf2.GetN(),bf2.GetL(),0,
                                          bf1.GetExponent(0),bf2.GetExponent(0),
                                          tatom1, tatom2 );
            xzxz_overlap = PrimitiveOverlap(bf1.GetN(),bf1.GetL(),1,
                                           bf2.GetN(),bf2.GetL(),1,
                                           bf1.GetExponent(0),bf1.GetExponent(0),
                                           tatom1,tatom2 );
            xyxy_overlap = PrimitiveOverlap(bf1.GetN(),bf1.GetL(),2,
                                           bf2.GetN(),bf2.GetL(),2,
                                           bf1.GetExponent(0),bf1.GetExponent(0),
                                           tatom1,tatom2 );
          }
          else
          {
            z2z2_overlap = bf2.GetCoefficient(0)*PrimitiveOverlap(bf1.GetN(),bf1.GetL(),0,
                                                                  bf2.GetN(),bf2.GetL(),0,
                                                                  bf1.GetExponent(0),bf2.GetExponent(0),
                                                                  tatom1, tatom2 ) +
                           bf2.GetCoefficient(1)*PrimitiveOverlap(bf1.GetN(),bf1.GetL(),0,
                                                                  bf2.GetN(),bf2.GetL(),0,
                                                                  bf1.GetExponent(0),bf2.GetExponent(1),
                                                                  tatom1, tatom2 );
            xzxz_overlap = bf2.GetCoefficient(0)*PrimitiveOverlap(bf1.GetN(),bf1.GetL(),1,
                                                                  bf2.GetN(),bf2.GetL(),1,
                                                                  bf1.GetExponent(0),bf2.GetExponent(0),
                                                                  tatom1, tatom2 ) +
                           bf2.GetCoefficient(1)*PrimitiveOverlap(bf1.GetN(),bf1.GetL(),1,
                                                                  bf2.GetN(),bf2.GetL(),1,
                                                                  bf1.GetExponent(0),bf2.GetExponent(1),
                                                                  tatom1, tatom2 );
            xyxy_overlap = bf2.GetCoefficient(0)*PrimitiveOverlap(bf1.GetN(),bf1.GetL(),2,
                                                                  bf2.GetN(),bf2.GetL(),2,
                                                                  bf1.GetExponent(0),bf2.GetExponent(0),
                                                                  tatom1, tatom2 ) +
                           bf2.GetCoefficient(1)*PrimitiveOverlap(bf1.GetN(),bf1.GetL(),2,
                                                                  bf2.GetN(),bf2.GetL(),2,
                                                                  bf1.GetExponent(0),bf2.GetExponent(1),
                                                                  tatom1, tatom2 );
          }
        }
        else
        {
          if ( bf2.GetNumPrimitives() == 1 )
          {
            z2z2_overlap = bf1.GetCoefficient(0)*PrimitiveOverlap(bf1.GetN(),bf1.GetL(),0,
                                                                  bf2.GetN(),bf2.GetL(),0,
                                                                  bf1.GetExponent(0),
                                                                  bf2.GetExponent(0),
                                                                  tatom1, tatom2 ) +
                           bf1.GetCoefficient(1)*PrimitiveOverlap(bf1.GetN(),bf1.GetL(),0,
                                                                  bf2.GetN(),bf2.GetL(),0,
                                                                  bf1.GetExponent(1),
                                                                  bf2.GetExponent(0),
                                                                  tatom1, tatom2 );
            xzxz_overlap = bf1.GetCoefficient(0)*PrimitiveOverlap(bf1.GetN(),bf1.GetL(),1,
                                                                  bf2.GetN(),bf2.GetL(),1,
                                                                  bf1.GetExponent(0),
                                                                  bf2.GetExponent(0),
                                                                  tatom1,tatom2 )+
                           bf1.GetCoefficient(1)*PrimitiveOverlap(bf1.GetN(),bf1.GetL(),1,
                                                                  bf2.GetN(),bf2.GetL(),1,
                                                                  bf1.GetExponent(1),
                                                                  bf2.GetExponent(0),
                                                                  tatom1,tatom2 );
            xyxy_overlap = bf1.GetCoefficient(0)*PrimitiveOverlap(bf1.GetN(),bf1.GetL(),2,
                                                                  bf2.GetN(),bf2.GetL(),2,
                                                                  bf1.GetExponent(0),
                                                                  bf2.GetExponent(0),
                                                                  tatom1,tatom2 )+
                           bf1.GetCoefficient(1)*PrimitiveOverlap(bf1.GetN(),bf1.GetL(),2,
                                                                  bf2.GetN(),bf2.GetL(),2,
                                                                  bf1.GetExponent(1),
                                                                  bf2.GetExponent(0),
                                                                  tatom1,tatom2 );
          }
          else
          {
            z2z2_overlap = bf1.GetCoefficient(0)*bf2.GetCoefficient(0)* 
                           PrimitiveOverlap(bf1.GetN(),bf1.GetL(),0,bf2.GetN(),bf2.GetL(),0,
                                            bf1.GetExponent(0),bf2.GetExponent(0),
                                            tatom1, tatom2 ) +
                           bf1.GetCoefficient(1)*bf2.GetCoefficient(0)* 
                           PrimitiveOverlap(bf1.GetN(),bf1.GetL(),0,bf2.GetN(),bf2.GetL(),0,
                                            bf1.GetExponent(1),bf2.GetExponent(0),
                                            tatom1, tatom2 ) +
                           bf1.GetCoefficient(0)*bf2.GetCoefficient(1)* 
                           PrimitiveOverlap(bf1.GetN(),bf1.GetL(),0,bf2.GetN(),bf2.GetL(),0,
                                            bf1.GetExponent(0),bf2.GetExponent(1),
                                            tatom1, tatom2 ) +
                           bf1.GetCoefficient(1)*bf2.GetCoefficient(1)* 
                           PrimitiveOverlap(bf1.GetN(),bf1.GetL(),0,bf2.GetN(),bf2.GetL(),0,
                                            bf1.GetExponent(1),bf2.GetExponent(1),
                                            tatom1, tatom2 );
            xzxz_overlap = bf1.GetCoefficient(0)*bf2.GetCoefficient(0)* 
                           PrimitiveOverlap(bf1.GetN(),bf1.GetL(),1,bf2.GetN(),bf2.GetL(),1,
                                            bf1.GetExponent(0),bf2.GetExponent(0),
                                            tatom1, tatom2 ) +
                           bf1.GetCoefficient(1)*bf2.GetCoefficient(0)* 
                           PrimitiveOverlap(bf1.GetN(),bf1.GetL(),1,bf2.GetN(),bf2.GetL(),1,
                                            bf1.GetExponent(1),bf2.GetExponent(0),
                                            tatom1, tatom2 ) +
                           bf1.GetCoefficient(0)*bf2.GetCoefficient(1)* 
                           PrimitiveOverlap(bf1.GetN(),bf1.GetL(),1,bf2.GetN(),bf2.GetL(),1,
                                            bf1.GetExponent(0),bf2.GetExponent(1),
                                            tatom1, tatom2 ) +
                           bf1.GetCoefficient(1)*bf2.GetCoefficient(1)* 
                           PrimitiveOverlap(bf1.GetN(),bf1.GetL(),1,bf2.GetN(),bf2.GetL(),1,
                                            bf1.GetExponent(1),bf2.GetExponent(1),
                                            tatom1, tatom2 );
            xyxy_overlap = bf1.GetCoefficient(0)*bf2.GetCoefficient(0)* 
                           PrimitiveOverlap(bf1.GetN(),bf1.GetL(),2,bf2.GetN(),bf2.GetL(),2,
                                            bf1.GetExponent(0),bf2.GetExponent(0),
                                            tatom1, tatom2 ) +
                           bf1.GetCoefficient(1)*bf2.GetCoefficient(0)* 
                           PrimitiveOverlap(bf1.GetN(),bf1.GetL(),2,bf2.GetN(),bf2.GetL(),2,
                                            bf1.GetExponent(1),bf2.GetExponent(0),
                                            tatom1, tatom2 ) +
                           bf1.GetCoefficient(0)*bf2.GetCoefficient(1)* 
                           PrimitiveOverlap(bf1.GetN(),bf1.GetL(),2,bf2.GetN(),bf2.GetL(),2,
                                            bf1.GetExponent(0),bf2.GetExponent(1),
                                            tatom1, tatom2 ) +
                           bf1.GetCoefficient(1)*bf2.GetCoefficient(1)* 
                           PrimitiveOverlap(bf1.GetN(),bf1.GetL(),2,bf2.GetN(),bf2.GetL(),2,
                                            bf1.GetExponent(1),bf2.GetExponent(1),
                                            tatom1, tatom2 );
          }
        } 
        ovlp[12] = z2z2_overlap;
        ovlp[18] = -xzxz_overlap;
        ovlp[0] = xyxy_overlap;
        RotateDD(ovlp,atom1,atom2,distance);
        return ovlp;
      }
      else
      {
        std::cout<<"ERROR! One of the BFs for DD Overlap isn't correct!!!!\n";
        std::exit(1);
      }
    } 
}

#endif
