
#ifndef CONSTANT_h
#define CONSTANT_h

/*----------------------------------------------------------------------------*/

    /*######################################################################
     
     revision log:
        10 Sept. 2021.
     
     ######################################################################*/

/*----------------------------------------------------------------------------*/

//Ratio of circumference to diameter
//3.1415926535897932384626433832795028841971
#define C_Pi 3.141592653589793

#define C_Pi2 1.5707963267948966

#define C_sqrtpi 1.7724538509055159

//Unified atomic mass unit kg^-1
#define C_m0 1.66053904e-27

#define C_me 9.109382e-31

#define Electron_charge 1.602176462e-19

//Larmor_frequency C_Nul*B in s^-1 (B in Gauss)
#define C_Nul 1.3996e6

//Speed of light
#define C_c 2.99792458E+08

//Boltzmann constant
#define C_Kb 1.3806e-23

//Planck constant
#define C_h 6.62606896e-34

//solar radius (cm)
#define C_Solarradus 6.9626e10

//Thomson scattering cross section m^2
#define C_sigmaT 6.65246e-29

//Hydrogen abundance in the sun (mass)
#define C_H_Ratio 0.7346

//Helium abundance in the sun (mass)
#define C_He_Ratio 0.2485

//Hydrogen
#define C_mH 1.008

//Helium
#define C_mHe 4.0026

//Ferrum
#define C_mFe 55.845

//C_H2E = (H_Density+He_Density)/H_Density, 
//    where H_Density = C_H_Ratio/C_mH He_Density = C_He_Ratio/C_mHe*2

// ratio between hydron and eletron (from Chianti guide)
#define C_H2E 0.83

// 1/C_H2E
#define C_E2H 1.2048192771084338

//log10(C_H2E)
//#define C_log10H2E 0.06832764757541908

//squrt(2.0)
#define C_sqrt2 1.4142135623730951

//squrt(3.0)
#define C_sqrt3 1.7320508075688772

/*----------------------------------------------------------------------------*/

#endif /* CONSTANT_h */
