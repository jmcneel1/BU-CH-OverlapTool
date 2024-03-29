#ifndef _bueht_parameters
#define _bueht_parameters

#include <bueht_constants.h>

namespace BUEHT
{

/*
  James McNeely

  The following definitions store parameters used by the program
  including basis function information.
  
*/

/*
  A paramater object is a "row" of the parameters
  listed below
*/

struct Parameter
{
  short an;
  short n;
  short l;
  double hii;
  double z1, z2;
  double c1, c2;

  bool operator < ( const Parameter & rhs )
  {
    if ( rhs.an  == an )
    {
      if ( rhs.n > n )
      {
        return true;
      }
      else
      {
        if ( rhs.n == n )
        {
          if ( rhs.l > l )
          {
            return true;
          }
        }
      }
    }
    else
    {
      if ( rhs.an > an )
      {
        return true;
      }
    }
    return false;
  }

};  // End of struct

/*
{Atomic Number, n, l, Hii (hartree), Exp1, Exp2, Coeff1, Coeff2}
These are taken from Greg Landrum's YAeHMOP
(https://github.com/greglandrum/yaehmop.git)
*/

Parameter bueht_eht_parameters [] = {
{1,1,0,-0.499791,1.3000,0.0000,1.0000,0.0000}, // H 1s Updated yaehmop 
{2,1,0,-0.859934,1.6880,0.0000,1.0000,0.0000}, // He 1s Updated yaehmop
{3,2,0,-0.198446,0.6500,0.0000,1.0000,0.0000}, // Li 2s Updated yaehmop
{3,2,1,-0.128623,0.6500,0.0000,1.0000,0.0000}, // Li 2p Updated yaehmop
{4,2,0,-0.367493,0.9750,0.0000,1.0000,0.0000}, // Be 2s Updated yaehmop
{4,2,1,-0.220496,0.9750,0.0000,1.0000,0.0000}, // Be 2p Updated yaehmop
{5,2,0,-0.558590,1.3000,0.0000,1.0000,0.0000}, // B 2s Updated yaehmop
{5,2,1,-0.312369,1.3000,0.0000,1.0000,0.0000}, // B 2p Updated yaehmop
{6,2,0,-0.786435,1.6250,0.0000,1.0000,0.0000}, // C 2s Updated yaehmop
{6,2,1,-0.418942,1.6250,0.0000,1.0000,0.0000}, // C 2p Updated yaehmop
{7,2,0,-0.955482,1.9500,0.0000,1.0000,0.0000}, // N 2s Updated yaehmop
{7,2,1,-0.492441,1.9500,0.0000,1.0000,0.0000}, // N 2p Updated yaehmop
{8,2,0,-1.187003,2.2750,0.0000,1.0000,0.0000}, // O 2s Updated yaehmop
{8,2,1,-0.543890,2.2750,0.0000,1.0000,0.0000}, // O 2p Updated yaehmop
{9,2,0,-1.469973,2.4250,0.0000,1.0000,0.0000}, // F 2s Updated yaehmop
{9,2,1,-0.665163,2.4250,0.0000,1.0000,0.0000}, // F 2p Updated yaehmop
{10,2,0,-1.587571,2.8790,0.0000,1.0000,0.0000}, // Ne 2s Updated yaehmop
{10,2,1,-0.734986,2.8790,0.0000,1.0000,0.0000}, // Ne 2p Updated yaehmop
{11,3,0,-0.187422,0.7330,0.0000,1.0000,0.0000}, // Na 3s Updated yaehmop
{11,3,1,-0.110248,0.7330,0.0000,1.0000,0.0000}, // Na 3p Updated yaehmop
{12,3,0,-0.330744,1.1000,0.0000,1.0000,0.0000}, // Mg 3s Updated yaehmop
{12,3,1,-0.165372,1.1000,0.0000,1.0000,0.0000}, // Mg 3p Updated yaehmop
{13,3,0,-0.452017,1.1670,0.0000,1.0000,0.0000}, // Al 3s Updated yaehmop
{13,3,1,-0.238871,1.1670,0.0000,1.0000,0.0000}, // Al 3p Updated yaehmop
{13,3,2,-0.000000,1.1670,0.0000,1.0000,0.0000}, // Al 3d Updated yaehmop
{14,3,0,-0.635763,1.3830,0.0000,1.0000,0.0000}, // Si 3s Updated yaehmop
{14,3,1,-0.338094,1.3830,0.0000,1.0000,0.0000}, // Si 3p Updated yaehmop
{15,3,0,-0.683537,1.7500,0.0000,1.0000,0.0000}, // P 3s Updated yaehmop
{15,3,1,-0.514491,1.3000,0.0000,1.0000,0.0000}, // P 3p Updated yaehmop
{16,3,0,-0.734986,2.1220,0.0000,1.0000,0.0000}, // S 3s Updated yaehmop
{16,3,1,-0.404243,1.8270,0.0000,1.0000,0.0000}, // S 3p Updated yaehmop
{17,3,0,-0.966507,2.1830,0.0000,1.0000,0.0000}, // Cl 3s Updated yaehmop
{17,3,1,-0.521840,1.7330,0.0000,1.0000,0.0000}, // Cl 3p Updated yaehmop
{18,3,0,-0.000000,0.0000,0.0000,1.0000,0.0000}, // Ar 3s Updated yaehmop
{18,3,1,-0.000000,0.0000,0.0000,1.0000,0.0000}, // Ar 3p Updated yaehmop
{19,4,0,-0.159492,0.8740,0.0000,1.0000,0.0000}, // K 4s Updated yaehmop
{19,4,1,-0.100326,0.8740,0.0000,1.0000,0.0000}, // K 4p Updated yaehmop
{20,3,2,-0.000000,0.0000,0.0000,1.0000,0.0000}, // Ca 3d Updated yaehmop
{20,4,0,-0.257245,1.2000,0.0000,1.0000,0.0000}, // Ca 4s Updated yaehmop
{20,4,1,-0.146997,1.2000,0.0000,1.0000,0.0000}, // Ca 4p Updated yaehmop
{21,3,2,-0.312737,4.3500,1.7000,0.4228,0.7276}, // Sc 3d Updated yaehmop
{21,4,0,-0.325966,1.3000,0.0000,1.0000,0.0000}, // Sc 4s Updated yaehmop
{21,4,1,-0.101061,1.3000,0.0000,1.0000,0.0000}, // Sc 4p Updated yaehmop
{22,3,2,-0.397260,4.5500,1.4000,0.4206,0.7839}, // Ti 3d Updated yaehmop
{22,4,0,-0.329641,1.0750,0.0000,1.0000,0.0000}, // Ti 4s Updated yaehmop
{22,4,1,-0.199916,1.0750,0.0000,1.0000,0.0000}, // Ti 4p Updated yaehmop
{23,3,2,-0.404243,4.7500,1.7000,0.4755,0.7052}, // V 3d Updated yaehmop
{23,4,0,-0.323762,1.3000,0.0000,1.0000,0.0000}, // V 4s Updated yaehmop
{23,4,1,-0.202856,1.3000,0.0000,1.0000,0.0000}, // V 4p Updated yaehmop
{24,3,2,-0.412327,4.9500,1.8000,0.5060,0.6750}, // Cr 3d Updated yaehmop
{24,4,0,-0.318249,1.7000,0.0000,1.0000,0.0000}, // Cr 4s Updated yaehmop
{24,4,1,-0.192566,1.7000,0.0000,1.0000,0.0000}, // Cr 4p Updated yaehmop
{25,3,2,-0.428865,5.1500,1.7000,0.5139,0.6929}, // Mn 3d Updated yaehmop
{25,4,0,-0.358306,0.9700,0.0000,1.0000,0.0000}, // Mn 4s Updated yaehmop
{25,4,1,-0.216454,0.9700,0.0000,1.0000,0.0000}, // Mn 4p Updated yaehmop
{26,3,2,-0.463041,5.3500,2.0000,0.5505,0.6260}, // Fe 3d Updated yaehmop
{26,4,0,-0.334419,1.9000,0.0000,1.0000,0.0000}, // Fe 4s Updated yaehmop
{26,4,1,-0.195506,1.9000,0.0000,1.0000,0.0000}, // Fe 4p Updated yaehmop
{27,3,2,-0.484356,5.5500,2.1000,0.5680,0.6060}, // Co 3d Updated yaehmop
{27,4,0,-0.338461,2.0000,0.0000,1.0000,0.0000}, // Co 4s Updated yaehmop
{27,4,1,-0.194404,2.0000,0.0000,1.0000,0.0000}, // Co 4p Updated yaehmop
{28,3,2,-0.521840,5.7500,2.3000,0.5683,0.6292}, // Ni 3d Updated yaehmop
{28,4,0,-0.402405,2.1000,0.0000,1.0000,0.0000}, // Ni 4s Updated yaehmop
{28,4,1,-0.230418,2.1000,0.0000,1.0000,0.0000}, // Ni 4p Updated yaehmop
{29,3,2,-0.514491,5.9500,2.3000,0.5933,0.5744}, // Cu 3d Updated yaehmop
{29,4,0,-0.418942,2.2000,0.0000,1.0000,0.0000}, // Cu 4s Updated yaehmop
{29,4,1,-0.222701,2.2000,0.0000,1.0000,0.0000}, // Cu 4p Updated yaehmop
{30,4,0,-0.456059,2.0100,0.0000,1.0000,0.0000}, // Zn 4s Updated yaehmop
{30,4,1,-0.239973,1.7000,0.0000,1.0000,0.0000}, // Zn 4p Updated yaehmop
{31,4,0,-0.424602,1.8080,0.0000,1.0000,0.0000}, // Ga 4s
{31,4,1,-0.208516,1.3140,0.0000,1.0000,0.0000}, // Ga 4p
{32,4,0,-0.557046,2.0240,0.0000,1.0000,0.0000}, // Ge 4s
{32,4,1,-0.269336,1.5500,0.0000,1.0000,0.0000}, // Ge 4p
{33,4,0,-0.695150,2.2220,0.0000,1.0000,0.0000}, // As 4s
{33,4,1,-0.330156,1.7570,0.0000,1.0000,0.0000}, // As 4p
{34,4,0,-0.840163,2.4090,0.0000,1.0000,0.0000}, // Se 4s
{34,4,1,-0.392520,1.9490,0.0000,1.0000,0.0000}, // Se 4p
{35,4,0,-0.992709,2.5880,0.0000,1.0000,0.0000}, // Br 4s
{35,4,1,-0.457088,2.1310,0.0000,1.0000,0.0000}, // Br 4p
{36,4,0,-1.152936,2.7620,0.0000,1.0000,0.0000}, // Kr 4s
{36,4,1,-0.524192,2.3060,0.0000,1.0000,0.0000}  // Kr 4p
};

Parameter *bueht_param_ptr[] = {
  &bueht_eht_parameters[0],  // H
  &bueht_eht_parameters[1],  // He
  &bueht_eht_parameters[2],  // Li
  &bueht_eht_parameters[4],  // Be
  &bueht_eht_parameters[6],  // B
  &bueht_eht_parameters[8],  // C
  &bueht_eht_parameters[10], // Nd
  &bueht_eht_parameters[12], // O
  &bueht_eht_parameters[14], // F
  &bueht_eht_parameters[16], // Ne
  &bueht_eht_parameters[18], // Na
  &bueht_eht_parameters[20], // Mg
  &bueht_eht_parameters[22], // Al
  &bueht_eht_parameters[24], // Si
  &bueht_eht_parameters[26], // P
  &bueht_eht_parameters[28], // S
  &bueht_eht_parameters[30], // Cl
  &bueht_eht_parameters[32], // Ar
  &bueht_eht_parameters[34], // K
  &bueht_eht_parameters[37], // Ca
  &bueht_eht_parameters[40], // Sc
  &bueht_eht_parameters[43], // Ti
  &bueht_eht_parameters[46], // V
  &bueht_eht_parameters[49], // Cr
  &bueht_eht_parameters[52], // Mn
  &bueht_eht_parameters[55], // Fe
  &bueht_eht_parameters[58], // Co
  &bueht_eht_parameters[61], // Ni
  &bueht_eht_parameters[64], // Cu
  &bueht_eht_parameters[67], // Zn
  &bueht_eht_parameters[70], // Ga
  &bueht_eht_parameters[72], // Ge
  &bueht_eht_parameters[74], // As
  &bueht_eht_parameters[76], // Se
  &bueht_eht_parameters[78], // Br
  &bueht_eht_parameters[80]  // Kr
  
}; 

} // End of namespace BUEHT

#endif
