
#ifndef INIT_h
#define INIT_h

/*----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ALLOCATION.h"
#include "SPECIAL_FUNCTIONS.h"
#include "READ_INPUT.h"



#include <ctype.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <mpi.h>
#include "MPI_CTRL.h"
#include "READ_ATMO.h"
#include "IO.h"
#include "FCTSG.h"
#include "WIGNER.h"

/*----------------------------------------------------------------------------*/

extern int Init(STRUCT_INPUT *Input, STRUCT_ATOM *Atom, \
    STRUCT_SYN *Syn, STR_FCTSG *fctsg, STRUCT_OUT *Output, \
    STRUCT_MPI *Mpi);
    
extern int INIT_COEFF(STRUCT_SFB_COEFF *Parain, STRUCT_INPUT *Input);

extern double Baumbach(double radius);

/*----------------------------------------------------------------------------*/

#endif /* INIT_h */
