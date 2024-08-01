
#ifndef READ_INPUT_h
#define READ_INPUT_h

/*----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <complex.h>
#include "STR.h"
#include "ALLOCATION.h"
#include "SFB_TRANSFORM.h"
#include "MPI_CONTROL.h"
#include "ERROR.h"
#include "CONSTANT.h"
#include "LL04.h"

/*----------------------------------------------------------------------------*/

enum keywordtype {KEYWORD_REQUIRED, KEYWORD_DEFAULT, KEYWORD_OPTIONAL};

// keywords
typedef struct Struct_Keywords{

    char keyword[Key_Length];
    char line[Max_Line_Length];
    bool Set, Required;

}STRUCT_KEYS;

/*----------------------------------------------------------------------------*/

// lines
typedef struct Struct_Stokes_grid{

    double Stk[4], Stksav[4];

}STRUCT_STOKES_GRID;

// Thomson scattering
typedef struct Struct_Thomson_Grid{

    //double wavelength;
    double Kr, Kt;
    double val[2], valsav[2];

}STRUCT_THOMSON_GRID;

// incidenct radiation tensor
typedef struct Struct_Incident_grid{

    //J^0_0 J^2_q
    double *J00, **J2q;

}STRUCT_INCIDENT_GRID;

// magnetic field structure (local vertical coordinate)
typedef struct Struct_Field{

    // strength, inclination, azimuth 
    double B, ThetaB, PhiB;

    // three components of the magnetic field vector.
    double Br, Bt, Bp;

}STRUCT_MAG;

// model parameter strucure
typedef struct Struct_Parameter{

    // temperature, hydrogen density, electron density, 
    // velocity along x direction.
    double T, Tsav, ne, nesav, nH, nHsav, Vx, Bx;

    // magnetic field
    STRUCT_MAG Mag, Magsav;

    double *Ion, *Ionsav;

}STRUCT_PARA;

typedef struct Struct_Model_Grid{

    STRUCT_STOKES_GRID *Line;
    STRUCT_THOMSON_GRID *Thom;
    STRUCT_INCIDENT_GRID *Jkq;
    STRUCT_PARA Para;

    // coordinate
    double R, Theta, Phi, Rsq0;

    // spherical Bessel function;
    double **Legendre, **LegendreD, **sBessel, **sBesselD;
    complex double *Phiarray;

    complex double **T2Q;

    // pixel index
    int ipixel, ipspec;

}STRUCT_GRID;

/*----------------------------------------------------------------------------*/

// synthesis
typedef struct Struct_synthesis{

    // grids
    STRUCT_GRID *Grids;

    // numbers of grids and pixels, and the pixel number for spectra
    int ngrids, npixels, npspec;

}STRUCT_SYN;

/*----------------------------------------------------------------------------*/

// SFB coefficients
typedef struct Struct_SFB_Coefficients{

    complex double ***Coeff;
    int N, L, Lsquare, NLsquare, threshold_n;
    double sparsity;


    // invert the parameter or not.
    bool invt, in;

}STRUCT_SFB_COEFF;

// Thomson scattering information
typedef struct Struct_Thomson{

    //lambda, intensity
    double lambda, Intensity;

    //limb darkening coefficients;
    double u1, u2;

    //a constant
    double Const;

}STRUCT_THOMSON;

// 
typedef struct Struct_Spectrum{

    // number of grids
    int Nl;

    // Spectrum range and wavelength
    double range[2], *Lambda;

}STRUCT_SPEC;

// Input
typedef struct Struct_Input{

    //SFB coefficients for 4 parameters;
    STRUCT_SFB_COEFF Para[4];

    //zeros
    double **Zeros;

    // max L and N of the zeros;
    int Lzero, Nzero;

    //boundition condition
    enum boundary_condition BC;

    //Thomson scattering
    STRUCT_THOMSON *Thom;

    //verbose 
    int Verbose;

    //INVERSION or synthesis
    int Mode;

    //Path to the files
    char Path_Zero[Max_Line_Length], Path_Output[Max_Line_Length], \
      Path_Observation[Max_Line_Length], Path_Atmo[Max_Line_Length], \
      Path_coeff[Max_Line_Length];

    char **Path_Atom; 

    //number of atoms;
    int Natom;

    //N and K of the SFB coefficients
    int Nmax, Lmax;

    //max k value for the tensor, and delte k;
    int Kmax, Kdelta;
    
    //max J, L and S of atoms
    int rJmax, rLmax, rSmax;

    // number of the coefficients
    int ncoeff;
    
    //number of lines, Thomson, transitions, max atom transition number
    int Nline, NThom, Ntrans, Maxntran, nJKQ;

    // radius boundary for the synthesis
    double rmin, rmax, rsqmin, rsqmax; 
    
    // radius boundary for integration
    double rint, rsqint;

    //field of view, field of view (spectrum);
    double FOV[2][2], FOVSPEC[2][2];
    
    //delta x, y, z
    double dx, dy, dz;

    // size in solar y direction, and z dirction.
    int ny, nz, nys, nzs;
    
    //pertub values of each coefficients
    double perturb[4];

    //satured hanle effect
    bool symmetry, Saturated;

    //weights 
    int nweight;
    double *weight;
    
    //max step size
    double step_size;

    //coordinate for single grid synthesis
    double Ypos, Zpos, Xpos;
  
    // Tkq tensor
    complex double ***Tkq;

    // rotation matrix (from vertical to magnetic field coordinate)
    complex double **Dkmn[3];

    // max k value for T and J;
    int TQmax, JQmax;

    double Per;

    int inner_flag;

    int Ncut, Lcut, Niter;

    bool Bpotential;

    double Rlim;

    // number of spectrum;
    int Nspec, Nl;
    STRUCT_SPEC *Spec;

    // output V image
    bool OutputV;

    int Nstk;

}STRUCT_INPUT;

/*----------------------------------------------------------------------------*/

// level structure
typedef struct Struct_Level{

    // energy of the level
    double Energy;

    //J L S
    double J, L, S;

    //degeneracy, square root of degeneracy
    double deg, sqrt_deg;

    //lande factor
    double g;

    //max value of K, number
    int Kmax, nKQ;

    //rho^K_Q
    complex double **Rho;

}STRUCT_LV;

// transition structure
typedef struct Struct_Transition{

    // indexes of upper and lower levels
    int au, al;

    // Einstein coefficients
    double Aul, Bul, Blu;

    // frequency, and lambda, intensity, effective lande factor, 
    // and the delta factor in Egidio's book (Eq.13.28)
    double nu, lambda, Intens, geff, delta;

    // two constants to compute the Plank function;
    double cplank1, cplank2;

    // coefficients ll04
    double w[3];

    //complex double **J_KQ, **J_KQ0;
    complex double **JKQ;

    // limb darkening coefficients
    double u1, u2;

    // M1 transition or not, output or not
    bool M1;

    // Rates (LL04)
    double ***TA, *TE, ***TS, ***RA, ***RS;

    // const1: LL04 Chapter 10. 10.50 (\sqrt((2*Jl+1)/(2*Ju+1))*Blu/Aul
    // const2: h\nu/4/\pi*\sqrt(2*Ju+1)*Aul*dx*Rsun.
    double const1, const2;

}STRUCT_TRANS;

// collision structure
typedef struct Struct_Collision{

    // size of temperature;
    int NT;

    // temperature array;
    double *T;

    // collisional strength
    double ***Strength;

    // collisional rates for a given temperature
    double **Rates;

}STRUCT_COL;

// collision structure
typedef struct Struct_Ion{

    // size of temperature;
    int NT;

    // temperature array;
    double **frac;

}STRUCT_ION;

// Atom structure
typedef struct Struct_Atom{

    // level structure
    STRUCT_LV *LV;

    // transition structure
    STRUCT_TRANS *TR;

    // collision structure
    STRUCT_COL *col;

    // ion structure
    STRUCT_ION Ion;

    // collision data
    bool collision;

    // element
    char Element[10];

    // mass, abundance
    double Mass, Abund;

    // equation index of the level
    int *eqindx;

    // number of levels, transition, output lines, dimension of equations
    // number of J^K_Q
    int Nlevel, Ntrans, Nline, nEq;

    // two level
    bool TwoLv;

    int *iout;

    double cDopp;

}STRUCT_ATOM;

/*----------------------------------------------------------------------------*/

typedef struct Struct_Output{

    // synthesis and perturb
    double **synloc, **syntot;

    // radius, norm
    double *R, *norm;

    // Spectrum
    double **specloc, **spectot;

}STRUCT_OUT;

/*----------------------------------------------------------------------------*/

typedef struct Struct_Obervation{

    // pixel conuts for inversion, 
    int npixels, num;

    // Observation data[npixels, num]
    double **Data;

}STRUCT_OBSERVATION;

/*----------------------------------------------------------------------------*/

extern int Keywords_Conversion(STRUCT_KEYS Keywords[], \
    STRUCT_INPUT *Input, STRUCT_MPI *Mpi);

extern int RDINPUT(char Filename[], STRUCT_INPUT *Input, \
    STRUCT_MPI *Mpi);

extern int FREE_INPUT(STRUCT_INPUT *Input);

/*----------------------------------------------------------------------------*/

#endif /* READ_INPUT_h */
