
#ifndef READ_ATMO_h
#define READ_ATMO_h

/*----------------------------------------------------------------------------*/

#include <stdio.h>
#include "CONSTANT.h"
#include "READ_INPUT.h"
#include "LL04.h"
#include "STR.h"

/*----------------------------------------------------------------------------*/

// atmosphere model structure
typedef struct Struct_atmosphere{

    // type of the density (0: log10 nh; 1: nh; 2: log10 ne; 3: ne),
    // type of the magnetic field (0: Br Bt Bp; 1: B thetaB phiB), 
    // and type of the velocity (0: Vr Vt Vp; 1: V thetaV phiV) 
    int rhotype, btype, vtype;

    // sizes of R, Theta, Phi
    int nR, nTheta, nPhi;

    // radius, theta, and phi arrays
    double *R, *Theta, *Phi;

    // temperature and density
    double ***T, ***rho;

    // Br, Btheta, and Bphi (btype=0) or B, ThetaB, and PhiB (btype=1)
    double ***B1, ***B2, ***B3;

    // Vr, Vtheta, and Vphi (vtype=0) or V, ThetaV, and PhiV (vtype=1)
    double ***V1, ***V2, ***V3;

}STRUCT_ATMO;

/*----------------------------------------------------------------------------*/

extern void READ_ATMO(char *Filename, STRUCT_ATMO *Atmo);

extern void WRITE_ATMO(char *Filename, STRUCT_ATMO *Atmo);

extern void FREE_MODEL(STRUCT_ATMO *Atmo);

/*----------------------------------------------------------------------------*/

#endif /* READ_ATMO_h */
