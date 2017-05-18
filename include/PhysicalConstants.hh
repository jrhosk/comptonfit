#ifndef PhysicalConstants_h
#define PhysicalConstants_h

#include <iostream>
#include <cstdlib>

namespace pchep { 

  static const double Avogadro = 6.02214179e+23;    // per mole

  static const double c_light = 2.99792458e+8;
  static const double c_light_lol = 0.299792458;
  static const double c_squared = c_light * c_light;

  static const double h_Planck = 6.62606896e-34;
  static const double hbar_Planck = h_Planck/2*3.14159265358979323846;
  static const double hbarc = hbar_Planck * c_light;
  static const double hbarc_squared = hbarc * hbarc;
  static const double fine_struct = 0.007297352;


  static const double electron_charge = -1; 
  static const double e_SI = 1.602176487e-19; 
  static const double coulomb = electron_charge/e_SI;
  static const double classic_e_radius = 2.8179403227e-15;

  static const double electron_mass     = 0.510998910;       // MeV
  static const double electron_mass_GeV = 0.000510998910;    // GeV
  static const double   proton_mass     = 938.272013;        // MeV
  static const double  neutron_mass     = 939.56536;         // MeV
  static const double           amu_c2  = 931.494028;        // MeV
  static const double           amu     = amu_c2/c_squared;


};

#endif
