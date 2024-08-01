
#ifndef INIT_h
#define INIT_h

/*----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ALLOCATION.h"
#include "SPECIAL_FUNCTIONS.h"
#include "READ_INPUT.h"
#include "MPI_CONTROL.h"


#include <ctype.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <mpi.h>
#include "MPI_CONTROL.h"
#include "READ_ATMO.h"
#include "IO.h"

#include "FCTSG.h"
#include "WIGNER.h"

/*----------------------------------------------------------------------------*/

extern int Randomseeds(STRUCT_MPI *Mpi);

extern double Interpol_Linear3D(double ***Data, int i, int j, int k, \
    double ir, double jr, double kr);

extern int Init(STRUCT_INPUT *Input, STRUCT_ATOM *Atom, \
    STRUCT_SYN *Syn, STR_FCTSG *fctsg, STRUCT_OUT *Output, \
    STRUCT_MPI *Mpi);
    
extern int FREE_GRIDS(STRUCT_INPUT *Input, STRUCT_SYN *Syn, \
    STRUCT_MPI *Mpi);

extern void Get_Para(STRUCT_ATMO *Atmo, STRUCT_ATOM *Atom, \
    STRUCT_SYN *Syn, STRUCT_INPUT *Input, STRUCT_MPI *Mpi);
    
extern int Ion_fraction(STRUCT_ATOM *Atom, STRUCT_SYN *Syn, \
    STRUCT_INPUT *Input, STRUCT_MPI *Mpi);
    
extern int INIT_COEFF(STRUCT_SFB_COEFF *Parain, STRUCT_INPUT *Input);

extern double Baumbach(double radius);

extern void FREE_OUTPUT(STRUCT_OUT *Output, int Mode);

/*----------------------------------------------------------------------------*/

#endif /* INIT_h */
