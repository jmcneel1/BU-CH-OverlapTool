#ifndef _bueht_rotate
#define _bueht_rotate

#include "atom.h"

/*

  James McNeely

  This header provides all the necessary functions to rotate Overlaps in the
  prolate spheroidal coordinate system to the molecular system.

  Derivations are explicitly given below

*/

namespace BUEHT
{

  /*
    Throughout the following, we use the following rotation matrix:

    | cosBcosA  cosBsinA  -sinB |
    |    -sinA      cosA      0 |
    | sinBcosA  sinBsinA   cosB |

    Thus, the sigma component (which is in the Z-direction in the dimer
    system), is given by the bottom row of the above...
  */
  
  void RotateSP ( std::vector<double> & ovlp,
                  const BUEHT::Atom & atom1,
                  const BUEHT::Atom & atom2,
                  const double & distance )
  {
    /*
      ovlp[2] contains the sigma overlap in the dimer system

      There is no pi overlap for S-P interactions, and thus
      we simply need to rotate from the Z direction to x,y,z
      which involves the last row in the rotation matrix shown
      above

      R(Z) = sinBcosA x + sinBsinA y + cosB z
    */
    double xy_norm, cosbeta, sinbeta, cosalpha, sinalpha;
    double sigma = ovlp[2];
    std::vector<double> coord1 = atom1.GetCoordinates();
    std::vector<double> coord2 = atom2.GetCoordinates();
    coord2[0] -= coord1[0]; coord2[1] -= coord1[1]; coord2[2] -= coord1[2];
    xy_norm = std::pow(std::pow(distance,2.e0)-coord2[2]*coord2[2],0.5);
    if ( xy_norm < 1.e-5 )
    {
      cosalpha = 1.0;
      sinalpha = 0.0;
      sinbeta = 0.0;
    }
    else
    {
      cosalpha = coord2[0]/xy_norm;
      sinalpha = coord2[1]/xy_norm;
      sinbeta = xy_norm/distance;
    }
    cosbeta = coord2[2]/distance;
    ovlp[0] = sigma*sinbeta*cosalpha;
    ovlp[1] = sigma*sinbeta*sinalpha;
    ovlp[2] = sigma*cosbeta;
  }

  void RotatePP ( std::vector<double> & ovlp,
                  const BUEHT::Atom & atom1,
                  const BUEHT::Atom & atom2,
                  const double & distance )
  {
    double xy_norm, cosbeta, sinbeta, cosalpha, sinalpha;
    double sigma = ovlp[8];
    double pi = ovlp[0];
    std::vector<double> coord1 = atom1.GetCoordinates();
    std::vector<double> coord2 = atom2.GetCoordinates();
    coord2[0] -= coord1[0]; coord2[1] -= coord1[1]; coord2[2] -= coord1[2];
    xy_norm = std::pow(std::pow(distance,2.e0)-coord2[2]*coord2[2],0.5);
    if ( xy_norm < 1.e-5 )
    {
      cosalpha = 1.0;
      sinalpha = 0.0;
      sinbeta = 0.0;
    }
    else
    {
      cosalpha = coord2[0]/xy_norm;
      sinalpha = coord2[1]/xy_norm;
      sinbeta = xy_norm/distance;
    }
    cosbeta = coord2[2]/distance;
    ovlp[0] = sigma*sinbeta*cosalpha*sinbeta*cosalpha + // xx
              ( cosbeta*cosalpha*cosbeta*cosalpha +
                sinalpha*sinalpha )*pi;
    ovlp[1] = sigma*sinbeta*cosalpha*sinbeta*sinalpha + // xy
              ( cosbeta*cosalpha*cosbeta*sinalpha -
                sinalpha*cosalpha)*pi;
    ovlp[2] = sigma*sinbeta*cosalpha*cosbeta - // xz
              cosbeta*cosalpha*sinbeta*pi;
    ovlp[3] = sigma*sinbeta*sinalpha*sinbeta*cosalpha + // yx
              ( cosbeta*sinalpha*cosbeta*cosalpha -
                cosalpha*sinalpha ) * pi;
    ovlp[4] = sigma*sinbeta*sinalpha*sinbeta*sinalpha + // yy
              ( cosbeta*sinalpha*cosbeta*sinalpha +
                cosalpha*cosalpha ) * pi;
    ovlp[5] = sigma*sinbeta*sinalpha*cosbeta - // yz
              cosbeta*sinalpha*sinbeta*pi;
    ovlp[6] = sigma*cosbeta*sinbeta*cosalpha - // zx
              sinbeta*cosbeta*cosalpha*pi;
    ovlp[7] = sigma*cosbeta*sinbeta*sinalpha - 
              sinbeta*cosbeta*sinalpha*pi; // zy
    ovlp[8] = sigma*cosbeta*cosbeta + sinbeta*sinbeta*pi; // zz
  }

  void RotateSD ( std::vector<double> & ovlp,
                  const BUEHT::Atom & atom1,
                  const BUEHT::Atom & atom2,
                  const double & distance )
  {
    /*
      R(2Z2-X2-Y2) = 2(sBcA x + sBsA y + cB z)^2 - (cBcA x + cBsA y - sB z)^2-
                      (-sA x + cA y )^2
                   = x^2 ( 2s2Bc2A - c2Bc2A - s2A ) +
                     y^2 ( 2s2Bs2A - c2Bs2A - c2A ) +
                     z^2 ( 2c2B - s2B ) +
                     xy  ( 4sBcAsBsA - 2cBcAcBsA + 2sAcA ) +
                     xz  ( 4sBcAcB + 2cBcAsB ) + 
                     yz  ( 4sBsAcB + 2cBsAsB )
                   = x^2 ( 2s2Bc2A - c2Bc2A - s2A ) +
                     y^2 ( 2s2Bs2A - c2Bs2A - c2A ) +
                     2z^2(1 - 3/2 s2B) +
                     xy ( 2sAcA (2s2B - c2B + 1)) + 
                     xz ( 6sBcAcB ) +
                     yz ( 6sBsAcB )
                   = -x^2 (1 - 3/2 s2B - 2s2Bc2A + c2Bc2A + s2A - 1 + 3/2 s2B) +
                     -y^2 (1 - 3/2 s2B - 2s2Bs2A + c2Bs2A + c2A - 1 + 3/2 s2B)
                     +2z^2(1 - 3/2 s2B) +
                     xy ( 2sAcA (2s2B - 1 + s2B + 1)) +
                     xz ( 6sBcacB ) +
                     yz ( 6sBsAcB ) 
                   = 2z^2-x^2-y^2 (1 - 3/2 s2B ) +
                     -x^2 ( -2s2Bc2A + c2Bc2A + s2A - 1 + 3/2 s2B ) +
                     -y^2 ( -2s2Bs2A + c2Bs2A + c2A - 1 + 3/2 s2B ) +
                     xy ( 6 sAcAs2B ) +
                     xz ( 6 sBcAcB ) + 
                     yz ( 6 sBsAcB )
                   = 2z^2-x^2-y^2 (1 - 3/2 s2B ) +
                     xy ( 6 sAcAs2B ) +
                     xz ( 6 sBcAcB ) + 
                     yz ( 6 sBsAcB ) +
                     x^2 ( 2s2Bc2A - c2Bc2A - s2A + 1 - 3/2 s2B ) +
                     -y^2 ( -2s2Bs2A + c2Bs2A + 1 - s2A - 1 + 3/2 s2B ) 
                   = 2z^2-x^2-y^2 (1 - 3/2 s2B ) +
                     xy ( 6 sAcAs2B ) +
                     xz ( 6 sBcAcB ) + 
                     yz ( 6 sBsAcB ) +
                     x^2 ( 2s2B(1-s2A) - c2B(1-s2A) - s2A + 1 - 3/2 s2B ) +
                     -y^2 ( -s2B(2s2A-3/2) + c2Bs2A - s2A )
                   = 2z^2-x^2-y^2 (1 - 3/2 s2B ) +
                     xy ( 6 sAcAs2B ) +
                     xz ( 6 sBcAcB ) + 
                     yz ( 6 sBsAcB ) +
                     x^2 ( 2s2B - 2s2Bs2A - c2B + c2Bs2A - s2A + 1 - 3/2 s2B ) +
                     -y^2 ( -2s2Bs2A + c2Bs2A - s2A + 3/2 s2B )
                   = 2z^2-x^2-y^2 (1 - 3/2 s2B ) +
                     xy ( 6 sAcAs2B ) +
                     xz ( 6 sBcAcB ) + 
                     yz ( 6 sBsAcB ) +
                     x^2 ( 2s2B - 2s2Bs2A + s2B + c2Bs2A - s2A - 3/2 s2B ) +
                     -y^2 ( -2s2Bs2A + c2Bs2A - s2A + 3/2 s2B )
                   = 2z^2-x^2-y^2 (1 - 3/2 s2B ) +
                     xy ( 6 sAcAs2B ) +
                     xz ( 6 sBcAcB ) + 
                     yz ( 6 sBsAcB ) +
                     x^2 ( -2s2Bs2A + c2Bs2A - s2A + 3/2 s2B ) +
                     -y^2 ( -2s2Bs2A + c2Bs2A - s2A + 3/2 s2B )
                   = 2z^2-x^2-y^2 (1 - 3/2 s2B ) +
                     xy ( 6 sAcAs2B ) +
                     xz ( 6 sBcAcB ) + 
                     yz ( 6 sBsAcB ) +
                     x^2 ( s2B(-2s2A+3/2) + (1-s2B)s2A - s2A ) +
                     -y^2 ( s2B(-2s2A+3/2) + (1-s2B)s2A - s2A  )
                   = 2z^2-x^2-y^2 (1 - 3/2 s2B ) +
                     xy ( 6 sAcAs2B ) +
                     xz ( 6 sBcAcB ) + 
                     yz ( 6 sBsAcB ) +
                     x^2 ( s2B(-3s2A+3/2)) +
                     -y^2 ( s2B(-3s2A+3/2))
                   = 2z^2-x^2-y^2 (1 - 3/2 s2B ) +
                     xy ( 6 sAcAs2B ) +
                     xz ( 6 sBcAcB ) + 
                     yz ( 6 sBsAcB ) +
                     x^2 ( 3/2 s2B(-2s2A+1)) +
                     -y^2 ( 3/2 s2B(-2s2A+1))
                   = 2z^2-x^2-y^2 (1 - 3/2 s2B ) +
                     xy ( 6 sAcAs2B ) +
                     xz ( 6 sBcAcB ) + 
                     yz ( 6 sBsAcB ) +
                     x^2 ( 3/2 s2B(c2A-s2A)) +
                     -y^2 ( 3/2 s2B(c2A-s2A))
                N(2z2-x2-y2)/N(xy) = 1/sqrt(12)
                N(2z2-x2-y2)/N(xz) = 1/sqrt(12)
                N(2z2-x2-y2)/N(yz) = 1/sqrt(12)
                N(2z2-x2-y2)/N(x2-y2) = 1/sqrt(3)
    */
    double xy_norm, cosbeta, sinbeta, cosalpha, sinalpha, sqrt3;
    sqrt3 = std::sqrt(3);
    // We remember the order here is m=-2,-1,0,1,2, or 
    // xy, yz, z2, xz, x2-y2
    double sigma = ovlp[2];
    std::vector<double> coord1 = atom1.GetCoordinates();
    std::vector<double> coord2 = atom2.GetCoordinates();
    coord2[0] -= coord1[0]; coord2[1] -= coord1[1]; coord2[2] -= coord1[2];
    xy_norm = std::pow(std::pow(distance,2.e0)-coord2[2]*coord2[2],0.5);
    if ( xy_norm < 1.e-5 )
    {
      cosalpha = 1.0;
      sinalpha = 0.0;
      sinbeta = 0.0;
    }
    else
    {
      cosalpha = coord2[0]/xy_norm;
      sinalpha = coord2[1]/xy_norm;
      sinbeta = xy_norm/distance;
    }
    cosbeta = coord2[2]/distance;
    ovlp[0] = sigma*sqrt3*cosalpha*sinalpha*sinbeta*sinbeta;
    ovlp[1] = sigma*sqrt3*cosbeta*sinbeta*sinalpha;
    ovlp[2] = sigma*(1.0-1.5*sinbeta*sinbeta);
    ovlp[3] = sigma*sqrt3*cosbeta*sinbeta*cosalpha;
    ovlp[4] = sigma*sqrt3*0.5*sinbeta*sinbeta*(cosalpha*cosalpha-sinalpha*sinalpha);
  }

  void RotatePD ( std::vector<double> & ovlp,
                  const BUEHT::Atom & atom1,
                  const BUEHT::Atom & atom2,
                  const double & distance )
  /*
    the pi contribution contains the xz and yz orbital being rotated...
    thus
    Correcting for the norm, we know Nxz = Nyz = Nxy
    xz = (cBcA x + cBsA y - sB z)*(sBcA x + sBsA y + cB z)
       =  cBcAsBcA x^2 + cBsAsBsA y^2 - sBcB z^2
        + 2cBcAsBsA xy + ( cBcAcB - sBsBcA ) xz 
        + ( cBsAcB - sBsBsA ) yz
       = cBcAsBcA x^2 + cBsAsBsA y^2 - sBcB z^2
        + 2cBcAsBsA xy 
        + cA ( cBcB - sBsB ) xz 
        + ( cBsAcB - sBsBsA ) yz
       = cBcAsBcA x^2 + cBsAsBsA y^2 - sBcB z^2
        + 2cBcAsBsA xy 
        + cA ( cBcB - sBsB ) xz 
        + sA ( cBcB - sBsB ) yz
       = 2z2(-1/2 sBcB)
        - x2 ( -1/2 sBcB ) + ( cBcAsBcA - 1/2 sBcB ) x^2 + cBsAsBsA y^2 
        + 2cBcAsBsA xy 
        + cA ( cBcB - sBsB ) xz 
        + sA ( cBcB - sBsB ) yz
       = 2z2(-1/2 sBcB)
        - x2 ( -1/2 sBcB ) + ( cBcAsBcA - 1/2 sBcB ) x^2 
        - y2 ( -1/2 sBcB ) - ( -cBsAsBsA + 1/2 sBcB ) y^2 
        + 2cBcAsBsA xy 
        + cA ( cBcB - sBsB ) xz 
        + sA ( cBcB - sBsB ) yz
       = 2z2(-1/2 sBcB) - x2 (-1/2 sBcB) - y2 (-1/2 sBcB)
        + ( cBcAsBcA - 1/2 sBcB ) x^2 
        - ( -cBsAsBsA + 1/2 sBcB ) y^2 
        + 2cBcAsBsA xy 
        + cA ( cBcB - sBsB ) xz 
        + sA ( cBcB - sBsB ) yz
       = 2z2(-1/2 sBcB) - x2 (-1/2 sBcB) - y2 (-1/2 sBcB)
        + ( c2AcBsB - 1/2 sBcB ) x^2 
        - ( -s2AcBsB + 1/2 sBcB ) y^2 
        + 2cBcAsBsA xy 
        + cA ( cBcB - sBsB ) xz 
        + sA ( cBcB - sBsB ) yz
       = 2z2(-1/2 sBcB) - x2 (-1/2 sBcB) - y2 (-1/2 sBcB)
        + ( c2AcBsB - 1/2 sBcB ) x^2 
        - ( -(1-c2A)cBsB + 1/2 sBcB ) y^2 
        + 2cBcAsBsA xy 
        + cA ( cBcB - sBsB ) xz 
        + sA ( cBcB - sBsB ) yz
       = 2z2(-1/2 sBcB) - x2 (-1/2 sBcB) - y2 (-1/2 sBcB)
        + ( c2AcBsB - 1/2 sBcB ) x^2 
        - ( c2AcBsB - 1/2 sBcB ) y^2 
        + 2cBcAsBsA xy 
        + cA ( cBcB - sBsB ) xz 
        + sA ( cBcB - sBsB ) yz
       = -sqrt(12)/2 * (2z2(sBcB) - x2 (sBcB) - y2 (sBcB)
        + 2 ( c2AcBsB - 1/2 sBcB ) ( x^2 - y^2  )
        + 2cBcAsBsA xy 
        + cA ( cBcB - sBsB ) xz 
        + sA ( cBcB - sBsB ) yz
       = -sqrt(3) * (2z2(sBcB) - x2 (sBcB) - y2 (sBcB)
        + ( 2c2AcBsB - sBcB ) ( x^2 - y^2  )
        + 2cBcAsBsA xy 
        + cA ( cBcB - sBsB ) xz 
        + sA ( cBcB - sBsB ) yz
       = -sqrt(3) * (2z2(sBcB) - x2 (sBcB) - y2 (sBcB)
        + ( cBsB(2c2A - 1) ) ( x^2 - y^2  )
        + 2cBcAsBsA xy 
        + cA ( cBcB - sBsB ) xz 
        + sA ( cBcB - sBsB ) yz
       = -sqrt(3) * (2z2(sBcB) - x2 (sBcB) - y2 (sBcB)
        + ( cBsB(c2A-s2A) ) ( x^2 - y^2  )
        + 2cBcAsBsA xy 
        + cA ( cBcB - sBsB ) xz 
        + sA ( cBcB - sBsB ) yz

    yz = (-sA x + cA y )*(sBcA x + sBsA y + cB z)
       = -sAsBcA x2 + cAsBsA y2 + (-sAsBsA + cAsBcA ) xy
         + ( -sAcB ) xz + cAcB yz
       = 2*((-sAsBcA) x2 - (-sAsBcA) y2 )
         + (-sAsBsA + cAsBcA ) xy
         + (-sAcB) xz 
         + cAcB yz
       = -2*((sAsBcA) (x2-y2))
         + (-sAsBsA + cAsBcA ) xy
         + (-sAcB) xz 
         + cAcB yz
       = -((sB2sAcA) (x2-y2))
         + (-sAsBsA + cAsBcA ) xy
         + (-sAcB) xz 
         + cAcB yz
       = -((sB2sAcA) (x2-y2))
         + (sB(-s2A+c2A)) xy
         + (-sAcB) xz 
         + cAcB yz
  */
  {
    double xy_norm, cosbeta, sinbeta, cosalpha, sinalpha, sqrt3;
    sqrt3 = std::sqrt(3);
    // We remember the order here is m=-2,-1,0,1,2, or 
    // xy, yz, z2, xz, x2-y2
    // p orbitals is x,y,z
    double sigma = ovlp[12];
    double pi = ovlp[3];
    std::vector<double> coord1 = atom1.GetCoordinates();
    std::vector<double> coord2 = atom2.GetCoordinates();
    coord2[0] -= coord1[0]; coord2[1] -= coord1[1]; coord2[2] -= coord1[2];
    xy_norm = std::pow(std::pow(distance,2.e0)-coord2[2]*coord2[2],0.5);
    if ( xy_norm < 1.e-5 )
    {
      cosalpha = 1.0;
      sinalpha = 0.0;
      sinbeta = 0.0;
    }
    else
    {
      cosalpha = coord2[0]/xy_norm;
      sinalpha = coord2[1]/xy_norm;
      sinbeta = xy_norm/distance;
    }
    cosbeta = coord2[2]/distance;
    ovlp[0] = sigma*sqrt3*cosalpha*sinalpha*sinbeta*sinbeta*sinbeta*cosalpha + // PZ-DZ2 -> px - dxy
              pi * ( cosalpha*cosbeta*2*cosbeta*cosalpha*sinbeta*sinalpha - 
                   sinalpha*sinbeta*(cosalpha*cosalpha-sinalpha*sinalpha) ); // PX-DXZ + PY-DYZ -> px-dxy 
    ovlp[1] = sigma*sqrt3*cosbeta*sinbeta*sinalpha*sinbeta*cosalpha + // PZ-DZ2 -> px-dyz
              pi * ( cosbeta*cosalpha*sinalpha*(cosbeta*cosbeta - sinbeta*sinbeta ) -
                     sinalpha*cosalpha*cosbeta ); // PX-DXZ + PY-DYZ -> px-dyz
    ovlp[2] = sigma*(1.0-1.5*sinbeta*sinbeta)*sinbeta*cosalpha + // PZ-DZ2 -> px - dz2
              pi * ( -cosbeta*cosalpha*sqrt3*cosbeta*sinbeta ); // PX-DXZ + PY-DYZ -> px-dz2
    ovlp[3] = sigma*sinbeta*cosalpha*sqrt3*cosbeta*sinbeta*cosalpha + // PZ-DZ2 -> px-dxz
              pi * ( cosbeta*cosalpha*cosalpha*(cosbeta*cosbeta - sinbeta*sinbeta) +
                     sinalpha*sinalpha*cosbeta ); // PX-DXZ + PY-DYZ -> px-dxz
    ovlp[4] = sigma*sinbeta*cosalpha*sqrt3*0.5*sinbeta*sinbeta*(cosalpha*cosalpha-sinalpha*sinalpha) +
              pi * ( cosbeta*cosalpha*cosbeta*sinbeta*(cosalpha*cosalpha-sinalpha*sinalpha) +
                     2*sinalpha*sinbeta*sinalpha*cosalpha ); // -> px-x2-y2
    ovlp[5] = sigma*sinbeta*sinalpha*sqrt3*cosalpha*sinalpha*sinbeta*sinbeta + 
              pi * ( cosbeta*sinalpha*2*cosbeta*cosalpha*sinbeta*sinalpha +
                     cosalpha*sinbeta*(cosalpha*cosalpha-sinalpha*sinalpha) ); // py-dxy
    ovlp[6] = sigma*sinbeta*sinalpha*sqrt3*sinbeta*sinalpha*cosbeta + 
              pi * ( cosbeta*sinalpha*sinalpha*(cosbeta*cosbeta-sinbeta*sinbeta) +
                     cosalpha*cosalpha*cosbeta ); // py-dyz
    ovlp[7] = sigma*sinbeta*sinalpha*(1-1.5*sinbeta*sinbeta) +
              pi * ( -cosbeta*sinalpha*sqrt3*sinbeta*cosbeta ); // py-dz2
    ovlp[8] = sigma*sinbeta*sinalpha*sqrt3*sinbeta*cosalpha*cosbeta +
              pi * ( cosbeta*sinalpha*cosalpha*(cosbeta*cosbeta-sinbeta*sinbeta) -
                     cosalpha*sinalpha*cosbeta ); // py-dxz
    ovlp[9] = sigma*sinbeta*sinalpha*sqrt(3)*0.5*sinbeta*sinbeta*(cosalpha*cosalpha-sinalpha*sinalpha) +
              pi * ( cosbeta*sinalpha*cosbeta*sinbeta*(cosalpha*cosalpha-sinalpha*sinalpha) -
                     cosalpha*2*sinbeta*sinalpha*cosalpha ); // py-dx2-y2
    ovlp[10] = sigma*cosbeta*sqrt3*sinalpha*cosalpha*sinbeta*sinbeta +
               pi * ( -sinbeta*2*cosbeta*cosalpha*sinbeta*sinalpha ); // pz-dxy
    ovlp[11] = sigma*cosbeta*sqrt3*sinbeta*sinalpha*cosbeta +
               pi * ( -sinbeta*sinalpha*(cosbeta*cosbeta-sinbeta*sinbeta)); // pz-dyz
    ovlp[12] = sigma*cosbeta*(1-1.5*sinbeta*sinbeta) + 
               pi * ( sinbeta*sqrt3*sinbeta*cosbeta ); // pz-dz2
    ovlp[13] = sigma*cosbeta*sqrt3*sinbeta*cosalpha*cosbeta +
               pi * ( -sinbeta*cosalpha*(cosbeta*cosbeta-sinbeta*sinbeta) ); // pz-dxz
    // ovlp[14] = sigma*cosbeta*sqrt3*0.5*sinbeta*sinbeta*(cosalpha*cosalpha-sinalpha*sinalpha) +
    //           pi * ( sinbeta*sqrt3*sinbeta*cosbeta ); // pz-dx2-y2
    ovlp[14] = sigma*cosbeta*sqrt3*0.5*sinbeta*sinbeta*(cosalpha*cosalpha-sinalpha*sinalpha) -
               pi * sinbeta*cosbeta*sinbeta*( cosalpha*cosalpha-sinalpha*sinalpha ); // pz-dx2-y2
  }

  void RotateDD ( std::vector<double> & ovlp,
                  const BUEHT::Atom & atom1,
                  const BUEHT::Atom & atom2,
                  const double & distance )
  /*
    We also need to rotate xy and x2-y2 here
    xy = (cBcA x + cBsA y - sB z)*(-sA x + cA y)
       = -cBcAsA x2 + cBsAcA y2 + (cBcAcA-sAcBsA) xy
         + sAsB xz - cAsB yz
       = -cBsAcA(x2-y2) +
         + cB(c2A-s2A) xy
         + sAsB xz
         + -cAsB yz
       =-2cBsAcA(x2-y2) +
         + cB(c2A-s2A) xy
         + sAsB xz
         + -cAsB yz
    x2-y2 = (cBcA x + cBsA y - sB z)2 - (-sA x + cA y)2
          = c2Bc2A x2 + c2Bs2A y2 + s2B z2 + 2cBcAcBsA xy -
            2cBcAsB xz - 2cBsAsB yz - s2A x2 - c2A y2 + 
            2 sAcA xy
          = (c2Bc2A-s2A) x2 + (c2Bs2A-c2A) y2 + s2B z2
            - 2cBcAsB xz 
            - 2cBsAsB yz
            + 2(cBcAcBsA+sAcA) xy
          = 2z2 (1/2s2B) - x2 (1/2s2B) - y2 (1/2 s2B)
            + (c2Bc2A-s2A+1/2s2B) x2 - (-c2Bs2A+c2A-1/2s2B) y2
            - 2cBcAsB xz 
            - 2cBsAsB yz
            + 2(cBcAcBsA+sAcA) xy
          = 2z2 (1/2s2B) - x2 (1/2s2B) - y2 (1/2 s2B)
            + (c2Bc2A-s2A+1/2s2B) x2 - (-c2B(1-c2A)+c2A-1/2s2B) y2
            - 2cBcAsB xz 
            - 2cBsAsB yz
            + 2(cBcAcBsA+sAcA) xy
          = 2z2 (1/2s2B) - x2 (1/2s2B) - y2 (1/2 s2B)
            + (c2Bc2A-s2A+1/2s2B) x2 - (-c2B+c2Bc2A+c2A-1/2s2B) y2
            - 2cBcAsB xz 
            - 2cBsAsB yz
            + 2(cBcAcBsA+sAcA) xy
          = 2z2 (1/2s2B) - x2 (1/2s2B) - y2 (1/2 s2B)
            + (c2Bc2A-s2A+1/2s2B) x2 - (-(1-s2B)+c2Bc2A+c2A-1/2s2B) y2
            - 2cBcAsB xz 
            - 2cBsAsB yz
            + 2(cBcAcBsA+sAcA) xy
          = 2z2 (1/2s2B) - x2 (1/2s2B) - y2 (1/2 s2B)
            + (c2Bc2A-s2A+1/2s2B) x2 - (-1+c2Bc2A+c2A+1/2s2B) y2
            - 2cBcAsB xz 
            - 2cBsAsB yz
            + 2(cBcAcBsA+sAcA) xy
          = 2z2 (1/2s2B) - x2 (1/2s2B) - y2 (1/2 s2B)
            + (c2Bc2A-s2A+1/2s2B) x2 - (-1+c2Bc2A+(1-s2A)+1/2s2B) y2
            - 2cBcAsB xz 
            - 2cBsAsB yz
            + 2(cBcAcBsA+sAcA) xy
          = 2z2 (1/2s2B) - x2 (1/2s2B) - y2 (1/2 s2B)
            + (c2Bc2A-s2A+1/2s2B) x2 - (c2Bc2A-s2A+1/2s2B) y2
            - 2cBcAsB xz 
            - 2cBsAsB yz
            + 2(cBcAcBsA+sAcA) xy  --> now include N corrections....
          = sqrt(3)/2*s2B*(2z2 - x2 - y2)
            + ((1-s2B)c2A-s2A+1/2s2B) x2 - ((1-s2B)c2A-s2A+1/2s2B) y2
            - cBcAsB xz 
            - cBsAsB yz
            + (sAcA(cBcB+1)) xy
          = sqrt(3)/2*s2B*(2z2 - x2 - y2)
            + (c2A+(1-c2B)(1/2-c2A)-(1-c2A)) x2 - (c2A+(1-cos2B)(1/2-c2A)-(1-c2A)) y2
            - cBcAsB xz 
            - cBsAsB yz
            + (sAcA(cBcB+1)) xy
          = sqrt(3)/2*s2B*(2z2 - x2 - y2)
            + (c2A+(-1/2c2B+c2Bc2A)-1/2) x2 - (c2A+(-1/2c2B+c2Bc2A)-1/2) y2
            - cBcAsB xz 
            - cBsAsB yz
            + (sAcA(cBcB+1)) xy
          = sqrt(3)/2*s2B*(2z2 - x2 - y2)
            + (1/2(2c2A-c2B+2c2Bc2A-1)) x2 - (1/2(2c2A-c2B+2c2Bc2A-1)) y2
            - cBcAsB xz 
            - cBsAsB yz
            + (sAcA(cBcB+1)) xy
          = sqrt(3)/2*s2B*(2z2 - x2 - y2)
            + (1/2(c2A-s2A-c2B+2c2Bc2A)) x2 - (1/2(c2A-s2A-c2B+2c2Bc2A-1)) y2
            - cBcAsB xz 
            - cBsAsB yz
            + (sAcA(cBcB+1)) xy
          = sqrt(3)/2*s2B*(2z2 - x2 - y2)
            + (1/2(c2A-s2A+c2B(2c2A-1)) x2 - (1/2(c2A-s2A+c2B(2c2A-1)) y2
            - cBcAsB xz 
            - cBsAsB yz
            + (sAcA(cBcB+1)) xy
          = sqrt(3)/2*s2B*(2z2 - x2 - y2)
            + (1/2(c2A-s2A+c2B(c2A-s2A)) x2 - (1/2(c2A-s2A+c2B(c2A-s2A)) y2
            - cBcAsB xz 
            - cBsAsB yz
            + (sAcA(cBcB+1)) xy
          = sqrt(3)/2*s2B*(2z2 - x2 - y2)
            + (1/2((c2A-s2A)(1+c2B)) x2 - (1/2((c2A-s2A)(1+c2B)) y2
            - cBcAsB xz 
            - cBsAsB yz
            + (sAcA(cBcB+1)) xy
  */
  {
    double xy_norm, cosbeta, sinbeta, cosalpha, sinalpha, sqrt3;
    sqrt3 = std::sqrt(3);
    // We remember the order here is m=-2,-1,0,1,2, or 
    // xy, yz, z2, xz, x2-y2
    double sigma = ovlp[12];
    double pi = ovlp[18];
    double delta = ovlp[0];
    std::vector<double> coord1 = atom1.GetCoordinates();
    std::vector<double> coord2 = atom2.GetCoordinates();
    coord2[0] -= coord1[0]; coord2[1] -= coord1[1]; coord2[2] -= coord1[2];
    xy_norm = std::pow(std::pow(distance,2.e0)-coord2[2]*coord2[2],0.5);
    if ( xy_norm < 1.e-5 )
    {
      cosalpha = 1.0;
      sinalpha = 0.0;
      sinbeta = 0.0;
    }
    else
    {
      cosalpha = coord2[0]/xy_norm;
      sinalpha = coord2[1]/xy_norm;
      sinbeta = xy_norm/distance;
    }
    cosbeta = coord2[2]/distance;
    std::vector<double> rotmatrix (25,0.0);
    rotmatrix[0] = cosbeta*(cosalpha*cosalpha-sinalpha*sinalpha); // xy-xy
    rotmatrix[1] = -cosalpha*sinbeta; // xy-yz
    rotmatrix[3] = sinalpha*sinbeta; // xy-xz
    rotmatrix[4] = -2.0*cosbeta*sinalpha*cosalpha; //xy-x2y2
    rotmatrix[5] = sinbeta*(cosalpha*cosalpha-sinalpha*sinalpha);  //yz-xy
    rotmatrix[6] = cosalpha*cosbeta; // yz-yz
    rotmatrix[8] = -sinalpha*cosbeta; // yz-xz
    rotmatrix[9] = -2.0*sinbeta*sinalpha*cosalpha; // yz-x2y2
    rotmatrix[10] = sqrt3*sinalpha*cosalpha*sinbeta*sinbeta; //z2-xy
    rotmatrix[11] = sqrt3*sinbeta*sinalpha*cosbeta; //z2-yz
    rotmatrix[12] = 1.0-1.5*sinbeta*sinbeta; //z2-z2
    rotmatrix[13] = sqrt3*sinbeta*cosalpha*cosbeta; //z2-xz
    rotmatrix[14] = sqrt3*0.5*sinbeta*sinbeta*(cosalpha*cosalpha-sinalpha*sinalpha); //z2-x2y2
    rotmatrix[15] = 2.0*cosbeta*cosalpha*sinbeta*sinalpha; // xz-xy
    rotmatrix[16] = sinalpha*(cosbeta*cosbeta-sinbeta*sinbeta); // xz-yz
    rotmatrix[17] = -sqrt3*cosbeta*sinbeta; // xz-z2
    rotmatrix[18] = cosalpha*(cosbeta*cosbeta-sinbeta*sinbeta); // xz-xz
    rotmatrix[19] = sinbeta*cosbeta*(cosalpha*cosalpha-sinalpha*sinalpha); // xz-x2y2
    rotmatrix[20] = sinalpha*cosalpha*(1.0+cosbeta*cosbeta); // x2y2-xy
    rotmatrix[21] = -cosbeta*sinalpha*sinbeta; // x2y2-yz
    rotmatrix[22] = sqrt3*0.5*sinbeta*sinbeta; // x2y2-z2
    rotmatrix[23] = -cosbeta*cosalpha*sinbeta; // x2y2-xz
    rotmatrix[24] = 0.5*(cosalpha*cosalpha-sinalpha*sinalpha)*(1.0+cosbeta*cosbeta); // x2y2-x2y2
    for ( unsigned int i = 0; i < 5; i++ )
    {
      for ( unsigned int j = i; j < 5; j++ )
      {
        ovlp[5*i+j] = sigma*rotmatrix[10+i]*rotmatrix[10+j] +
                     pi*rotmatrix[5+i]*rotmatrix[5+j] +
                     pi*rotmatrix[15+i]*rotmatrix[15+j] +
                     delta*rotmatrix[i]*rotmatrix[j] +
                     delta*rotmatrix[20+i]*rotmatrix[20+j];
        if ( i != j ) ovlp[5*j+i] = ovlp[5*i+j];
      }
    }          
  }

}

#endif