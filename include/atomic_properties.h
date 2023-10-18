#ifndef _bueht_atomic_properties
#define _bueht_atomic_properties

/*

  James McNeely

  This file contains some useful arrays that
  store information about the atoms modeled 
  with BUEHT

  The individual arrays are described below
*/

namespace BUEHT
{

int bueht_core_electrons [] = {
  0,  // H
  0,  // He
  2,  // Li 1s
  2,  // Be 1s
  2,  // B  1s
  2,  // C  1s
  2,  // N  1s
  2,  // O  1s
  2,  // F  1s
  2,  // Ne 1s
  10, // Na 1s2s2p
  10, // Mg 1s2s2p
  10, // Al 1s2s2p
  10, // Si 1s2s2p
  10, // P  1s2s2p
  10, // S  1s2s2p
  10, // Cl 1s2s2p
  10, // Ar 1s2s2p
  18, // K  1s2s2p3s3p
  18, // Ca 1s2s2p3s3p
  18, // Sc 1s2s2p3s3p
  18, // Ti 1s2s2p3s3p
  18, // V  1s2s2p3s3p
  18, // Cr 1s2s2p3s3p
  18, // Mn 1s2s2p3s3p
  18, // Fe 1s2s2p3s3p
  18, // Co 1s2s2p3s3p
  18, // Ni 1s2s2p3s3p
  18, // Cu 1s2s2p3s3p
  18, // Zn 1s2s2p3s3p
  28, // Ga 1s2s2p3s3p3d
  28, // Ge 1s2s2p3s3p3d
  28, // As 1s2s2p3s3p3d
  28, // Se 1s2s2p3s3p3d
  28, // Br 1s2s2p3s3p3d
  28  // Kr 1s2s2p3s3p3d
};

float bueht_atomic_masses [] = {
  1.008, // H
  4.002602, // He
  6.94, // Li
  9.0121831, // Be
  10.81, // B
  12.011, // C
  14.007, // N
  15.999, // O
  18.998403163, // F
  20.1797, // Ne
  22.98976928, // Na
  24.305, // Mg
  26.9815385, // Al
  28.085, // Si
  30.973761998, // P
  32.06, // S
  35.45, // Cl
  39.948, // Ar
  39.0983, // K
  40.078, // Ca
  44.955908, // Sc
  47.867, // Ti
  50.9415, // V
  51.9961, // Cr
  54.938044, //Mn
  55.845, // Fe
  58.933194, // Co
  58.6934, // Ni
  63.546, // Cu
  65.38, // Zn
  69.723, // Ga
  72.630, // Ge
  74.921595, // As
  78.971, // Se
  79.904, // Br
  83.798 // Kr
};

std::string bueht_atomic_symbols[] = {
      "H",
      "He",
      "Li",
      "Be",
      "B",
      "C",
      "N",
      "O",
      "F",
      "Ne",
      "Na",
      "Mg",
      "Al",
      "Si",
      "P",
      "S",
      "Cl",
      "Ar",
      "K",
      "Ca",
      "Sc",
      "Ti",
      "V",
      "Cr",
      "Mn",
      "Fe",
      "Co",
      "Ni",
      "Cu",
      "Zn",
      "Ga",
      "Ge",
      "As",
      "Se",
      "Br",
      "Kr"
    };

} // End of namespace BUEHT

#endif
