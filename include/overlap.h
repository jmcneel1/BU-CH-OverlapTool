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
   
      φ = (2ξ)^(n+0.5)/sqrt((2n)!)*r^(n-1)e^(-ξr)Y(θ,φ)

      We'll stick with the complex spherical harmonics, because the real 
      spherical harmonics, as given below:

      RYlm = (-1)^m/sqrt(2)*(Ylm + Ylm*) m > 0
      RYlm = Yl0 m == 0
      RYlm = (-1)^m/(i*sqrt(2))*(Yl|m| - Yl|m|*) m < 0

      As such, any reccurence relations we derive from the complex SH
      can be subsequently applied to real SH. This will be done later 
      to very equivalence...

      So ... the complex SH can be expressed as

      Y = (-1)^m*sqrt((2l+1)(l-m)!/(4π(l+m)!))*P(cosθ)e^(imφ)
        = (-1)^m*e^(imφ)*sqrt((2l+1)(l-m)!/(4π(l+m)!))*
          ((sinθ^m * d^(l+m)(cosθ^2-1)^l)/(2^l*l!*(d cosθ)^(l+m)))
      
      If we do a binomial expansion of (cosθ^2-1)^l, we have

      Y = (-1)^m*e^(imφ)*sqrt((2l+1)(l-m)!/(4π(l+m)!))*
          ((sinθ^m * d^(l+m) sum([k=0..l](l!/(k!(l-k)!))*
          cosθ^(2k)(-1)^(l-k)))/(2^l*l!*(d cosθ)^(l+m)))
      d^(l+m)/(d cosθ)^(l+m) (l!/(k!(l-k)!)) cosθ^(2k)^(-1)^(l-k) =
        (-1)^(l-k)(2k)!*(cosθ)^(2k-l-m)*l!/(k!*(l-k)!*(2k-l-m)!)

      So that

      Y = (-1)^m*e^(imφ)*sqrt((2l+1)(l-m)!/(4π(l+m)!))*
          ((sinθ^m * sum([k=0..l](l!/(k!(l-k)!))*
          (2k)!*cosθ^(2k-l-m)(-1)^(l-k)/(2k-l-m)!))/(2^l*l!))

      We now note that the factorials are only defined for 2k >= l+m, so

      Y = (-1)^m*e^(imφ)*sqrt((2l+1)(l-m)!/(4π(l+m)!))*
          ((sinθ^m * sum([k=1/2(l+m)..l](l!/(k!(l-k)!))*
          (2k)!*cosθ^(2k-l-m)(-1)^(l-k)/(2k-l-m)!))/(2^l*l!))

      Making the substitution j = l-k to standardize the summation

        = (-1)^m*e^(imφ)*sqrt((2l+1)(l-m)!/(4π(l+m)!))*
          ((sinθ^m * sum([j=0..1/2(l-m)](l!/((l-j)!(j)!))*
          (2(l-j))!*cosθ^(2(l-j)-l-m)(-1)^(j)/(2(l-j)-l-m)!))/(2^l*l!))

        = (-1)^m * e^(imφ) * sqrt((2l+1)(l-m)!/(4π(l+m)!)) *
          ( ( sinθ^m / (2^l*l!) * sum( [j=0..1/2(l-m)] (l!/((l-j)!(j)!)) *
          ( 2(l-j))! * cosθ^(l-m-2j) * (-1)^(j) / 
          ( l-m-2j )! )) )

      Now let's collect some of the terms above into a coefficient, which is 
      Overlap_C_lmj defined below

        Clmj = sqrt((l-m)!/(l+m)!) * ( (2(l-j))! * (-1)^j )/(2^l * j! * (l-j)! * (l-m-2j)! )

      This make the STO:

        φ = (2ξ)^(n+0.5)/sqrt((2n)!) * r^(n-1) * e^(-ξr) * (-1)^m * e^(imφ) * sinθ^m * 
          sqrt(2l+1/(4π)) * sum( [j=0..1/2(l-m)] Clmj * cosθ^(l-m-2j) ) 

      If we now look forward a bit, and realize that the Kronecker delta 
      δ(l1,l2) doesn't hold but δ(m1,m2) DOES hold moving from spherical to prolate 
      spheroidal coordinates, we can express the overlap integral as follows (also notice that
      the phase factor disappears because m is required to be an integer..., so
      (-1)^(2m) is 1 ):

        <φa|φb> = (2ξa)^(na+0.5) * (2ξb)^(nb+0.5) * sqrt( (2la+1) * (2lb+1) /
                  ((2*na)!*(2*nb)!)) / 2 *
                  Int[ sinθa^m * sinθb^m * ra^(na-1) * rb^(nb-1) * 
                       e^(-ξara) * e^(-ξbrb) * 
                       sum( [ja=0..1/2(la-m)] Clamja * cosθa^(la-m-2ja) ) *
                       sum( [jb=0..1/2(lb-m)] Clbmjb * cosθb^(lb-m-2jb) )
                  ,{θ,0,π},{r,0,inf}]

      This integral is still not feasible to solve, because of the disparate coordinate
      systems for the two atomic centers... now we move to prolate spheroidal coordinates.

      We first note that in our system contains a symmetry axis which is the vector
      connecting the two atoms whose overlap integral is being calculated. 

      Here, we assume that the vector joining the atomic centers is placed upon the z-axis.
      Without loss of generality, we can assume that φ is the same in both the PSC and spherical
      coordinate systems. Furthermore, the interatomic distance, R, is a constant to use in the 
      PSC.

      In PSC, for a single value of φ, we use the choice of PSC where z = 1/2 * R * mu * nu,
      x = 1/2 * R * sqrt( mu^2 - 1) * sqrt( 1 - nu^2 ) * cos(φ)
      y = 1/2 * R * sqrt( mu^2 - 1) * sqrt( 1 - nu^2 ) * sin(φ)

      The variables mu and nu thus have physical representations, in that when we hold
      mu constant, nu is a prolate spheroid (an ellipse rotated around the z-axis by phi)
      an if nu is held constant, we get hyperbola. Thus, any point is specified as a point
      intersecting an ellipse hyperbola rotated along the symmetry axis by φ.

      A nice property that can be divined from the above equations is that for any value of φ,
      the distance from any point to atom2 plus the distance to atom1 is a constant = R*mu. If
      the distances are subtracted, the constant is R*nu. 

      Let's first choose φ = 0, then we are in the xz plane. The distance between any point and
      atom2 is sqrt((1/2*R*mu*nu-1/2*R)^2+(1/2*R*sqrt(mu^2-1)*sqrt(1-nu^2))^2)
      = sqrt(1/4*R^2*(mu^2*nu^2-2*mu*nu+1+mu^2-mu^2*nu^2-1+nu^2))
      = sqrt(1/4*R^2*(-2*mu*nu+mu^2+nu^2)) = 1/2 * R * (mu-nu)
      For atom 1
      = sqrt((1/2*R*mu*nu+1/2*R)^2+1/4*R^2*(mu^2-1)*(1-nu^2))
      = sqrt(1/4*R^2*(mu^2*nu^2+2*mu*nu+1+mu^2-mu^2*nu^2-1+nu^2))
      = sqrt(1/4*R^2*(2*mu*nu+mu^2+nu^2)) = 1/2 * R * (mu+nu)

      This can be shown to hold for any value of φ, because for the foci, x=y=0, and the 
      x and y coordinates when squared will be 1/2 * R * (mu+nu) (cosφ^2+sinφ^2) ...

      With the above results defining rA = 1/2 * R * (mu+nu) and rB = 1/2 * R * (mu-nu),
      we now want to derive expressions for cosθa, cosθb, sinθa, and sinθb... so we'll
      just use the law of cosines because we have all three lengths now...

      cosθa = ( R^2 + rA^2 - rB^2 ) / ( 2 * R * rA )
            = ( R^2 + (1/2*R*(mu+nu))^2 - (1/2*R*(mu-nu))^2 ) / ( 2 * R * 1/2 * R * (mu+nu) )
            = ( 1/4 * R^2 * (4 + (mu+nu))^2 - (mu-nu)^2 ) / ( R^2 * (mu+nu) )
            = (1 + 1/4*(mu^2 + nu^2 + 2*mu*nu - mu^2 - nu^2 + 2 * mu * nu ))/(mu+nu)
            = (1+mu*nu)/(mu+nu)
      cosθb =  ( R^2 + rB^2 - rA^2 ) / ( 2 * R * rB )
            = ( R^2 + (1/2*R*(mu-nu))^2 - (1/2*R*(mu+nu))^2 ) / ( 2 * R * 1/2 * R * (mu-nu) )
            = ( 1/4 * R^2 * (4 + (mu-nu))^2 - (mu+nu)^2 ) / ( R^2 * (mu-nu) )
            = (1 + 1/4*(mu^2 + nu^2 - 2*mu*nu - mu^2 - nu^2 - 2 * mu * nu ))/(mu-nu)
            = (1-mu*nu)/(mu-nu)
      sinθa = sqrt(1-cosθa^2)
            = sqrt(1-(1+mu*nu)^2/(mu+nu)^2)
            = sqrt(((mu+nu)^2-(1+mu*nu)^2)/((mu+nu)^2))
            = sqrt((mu^2+nu^2+2*mu*nu-1-2*mu*nu-mu^2*nu^2)/(mu+nu)^2)
            = sqrt((mu^2+nu^2-1-mu^2*nu^2)/(mu+nu)^2)
            = sqrt((mu^2+1)*(1-nu^2))/(mu+nu)
      sinθb = sqrt(1-cosθb^2)
            = sqrt(1-(1-mu*nu)^2/(mu-nu)^2)
            = sqrt(((mu-nu)^2-(1-mu*nu)^2)/((mu-nu)^2))
            = sqrt((mu^2+nu^2-2*mu*nu-1+2*mu*nu+mu^2*nu^2)/(mu-nu)^2)
            = sqrt((mu^2+nu^2-1-mu^2*nu^2)/(mu-nu)^2)
            = sqrt((mu^2+1)*(1-nu^2))/(mu-nu)

      The integration volume element r^2*sinθ*dr*dθ*dφ also changes...

      The matrix of derivatives is as follows:

                          dx                          dy                    dz
            ________________________________________________________________________
           | - R * nu * sqrt(mu^2-1)*      - R * nu * sqrt(mu^2-1)*      1/2*R*mu
      dnu  | Cos[phi]/(2*sqrt(1-nu^2))     Sin[phi]/(2*sqrt(1-nu^2))
           |           | 
           | R * mu * sqrt(1-nu^2)*        R * mu * sqrt(1-nu^2)*        1/2*R*nu
      dmu  | Cos[phi]/(2*sqrt(mu^2-1))     Sin[phi]/(2*sqrt(mu^2-1))
           |
      dphi | - 1/2 * R * sqrt(mu^2-1) *     1/2 * R * sqrt(mu^2-1) *        0
           | sqrt(1-nu^2)*Sin[phi]          sqrt(1-nu^2)*Sin[phi]

      The determinant of this matrix is:

      R^3/8 * (mu-nu) * (mu+nu)

      Now we plug back in:

      <φa|φb> = (2ξa)^(na+0.5) * (2ξb)^(nb+0.5) * sqrt( (2la+1) * (2lb+1) / ((2na)!*(2nb)!)) / 2 *
                  Int[ (sqrt((mu^2+1)*(1-nu^2))/(mu+nu))^m * (sqrt((mu^2+1)*(1-nu^2))/(mu-nu))^m * 
                       (1/2*R*(mu+nu))^(na-1) * (1/2*R*(mu-nu))^(nb-1) * 
                       e^(-ξa(1/2*R*(mu+nu))) * e^(-ξb(1/2*R*(mu-nu))) * 
                       sum( [ja=0..1/2(la-m)] Clamja * ((1+mu*nu)/(mu+nu))^(la-m-2ja) ) *
                       sum( [jb=0..1/2(lb-m)] Clbmjb * ((1-mu*nu)/(mu-nu))^(lb-m-2jb) ) *
                       R^3/8*(nu-mu)*(mu+nu)
                  ,{mu,1,inf},{nu,-1,1}]

              = (2ξa)^(na+0.5)*(2ξb)^(nb+0.5)*sqrt((2la+1)*(2lb+1)/((2na)!*(2nb)!)) / 2 *
                  Int[ (R^(3+na-1+nb-1)*(mu^2+1)^m*(1-nu^2)^m *
                       (mu+nu)^(-m+na-1) * (mu-nu)^(-m+nb-1) * 
                       (1/2)^(na-1)*(1/2)^(nb-1) * 
                       e^(-ξa(1/2*R*(mu+nu))) * e^(-ξb(1/2*R*(mu-nu))) * 
                       sum( [ja=0..1/2(la-m)] Clamja * ((1+mu*nu)/(mu+nu))^(la-m-2ja) ) *
                       sum( [jb=0..1/2(lb-m)] Clbmjb * ((1-mu*nu)/(mu-nu))^(lb-m-2jb) ) *
                       1/8*(nu-mu)*(mu+nu)
                  ,{mu,1,inf},{nu,-1,1}]
              = (ξa)^(na+0.5)*(ξb)^(nb+0.5)*sqrt((2la+1)*(2lb+1)/((2na)!*(2nb)!)) / 2 *
                  Int[ (R^(na+nb+1)*(mu^2+1)^m*(1-nu^2)^m *
                       e^(-1/2*(ξa+ξb)*R*mu) * e^(-1/2*(ξa-ξb)*R*nu) * 
                       sum( [ja=0..1/2(la-m)]Clamja*(mu+nu)^(na-la+2ja)*((1+mu*nu))^(la-m-2ja))*
                       sum( [jb=0..1/2(lb-m)]Clbmjb*(mu-nu)^(nb-lb+2jb)*((1-mu*nu))^(lb-m-2jb))
                  ,{mu,1,inf},{nu,-1,1}]

            We can expand the terms as a binomial expansion:

            (mu^2+1)^m = sum([k=0..m] m!/(k!(m-k)!) * mu^2k)
            (1-nu^2)^m = sum([k=0..m] m!/(k!(m-k)!) * (-1)^k * nu^2k) 
            (mu+nu)^(na-la+2ja) = sum([k=0..na-la+2ja] (na-la+2ja)!/(k!(na-la+2ja-k)!) 
                                                       * mu^(na-la+2*ja-k) * nu^k )
            (mu-nu)^(nb-lb+2jb) = sum([k=0..nb-lb+2jb] (nb-lb+2jb)!/(k!(nb-lb+2jb-k)!) 
                                                       * (-1)^k * mu^(nb-lb+2*jb-k) * nu^k )
            (1+mu*nu)^(la-m-2ja) = sum([k=0..la-m-2ja] (la-m-2ja)!/(k!(la-m-2ja-k)!)
                                                       * mu^k * nu^k )
            (1-mu*nu)^(lb-m-2jb) = sum([k=0..lb-m-2jb] (lb-m-2jb)!/(k!(lb-m-2jb-k)!)
                                                       * (-1)^k * mu^k * nu^k )

            We now make variable susbstitutions ... where in order they will be ka, kb, 
            pa, pb, qa, and qb:

            (mu^2+1)^m = sum([ka=0..m] m!/(ka!(m-ka)!) * mu^2ka)
            (1-nu^2)^m = sum([kb=0..m] m!/(kb!(m-kb)!) * (-1)^kb * nu^2kb)
            (mu+nu)^(na-la+2ja) = sum([pa=0..na-la+2ja] (na-la+2ja)!/(pa!(na-la+2ja-pa)!) 
                                                        * mu^(na-la+2*ja-pa) * nu^pa )
            (mu-nu)^(nb-lb+2jb) = sum([pb=0..nb-lb+2jb] (nb-lb+2jb)!/(pb!(na-la+2ja-pb)!) 
                                                        * (-1)^pb * mu^(nb-lb+2*jb-pb) * nu^pb )
            (1+mu*nu)^(la-m-2ja) = sum([qa=0..la-m-2ja] (la-m-2ja)!/(qa!(la-m-2ja-qa)!)
                                                        * mu^qa * nu^qa )
            (1-mu*nu)^(lb-m-2jb) = sum([qb=0..lb-m-2jb] (lb-m-2jb)!/(qb!(lb-m-2jb-qb)!)
                                                        * (-1)^qb * mu^qb * nu^qb )

            Now we combine terms ...

            mu^(2*ka+na-la+2*ja-pa+nb-lb+2*jb-pb+qa+qb)
            nu^(2*kb+pa+pb+qa+qb)
            (-1)^(kb+pb+qb)
    */

    /*
      This defines the C coefficients discussed above...
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
