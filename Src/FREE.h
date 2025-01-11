
#ifndef FREE_h
#define FREE_h

/*----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <complex.h>
#include "READ_ATMO.h"
#include "READ_ATOM.h"
#include "READ_INPUT.h"

/*----------------------------------------------------------------------------*/

extern int FREE_INPUT(STRUCT_INPUT *Input);

extern int FREE_COEFF(STRUCT_SFB_COEFF *Para);

extern void FREE_MODEL(STRUCT_ATMO *Atmo);

extern void FREE_ATOM(STRUCT_ATOM *Atom, int Natom);
    
extern int FREE_GRIDS(STRUCT_INPUT *Input, STRUCT_SYN *Syn, \
    STRUCT_MPI *Mpi);

extern void FREE_OUTPUT(STRUCT_OUT *Output, STRUCT_INPUT *Input);

/*----------------------------------------------------------------------------*/

#endif /* FREE_h */
