
#ifndef FISTA_h
#define FISTA_h

#include <stdio.h>
#include "FORWARD.h"
#include "SFB_TRANSFORM.h"
#include "INIT.h"
#include "READ_ATOM.h"
#include "IO.h"

/*----------------------------------------------------------------------------*/

extern void FISTA(STRUCT_INPUT *Input, STRUCT_ATOM *Atom, STRUCT_SYN *Syn, \
    STRUCT_OUT *Output, STRUCT_OBSERVATION *Obs, STR_FCTSG *fctsg, \
    STRUCT_MPI *Mpi);

extern double LOSS_FUNCTION(STRUCT_OBSERVATION *Obs, STRUCT_OUT *Output, \
    STRUCT_INPUT *Input);

extern void array2coeff(double *qk,  STRUCT_INPUT *Input);

extern void coeff2array(double *qk,  STRUCT_INPUT *Input);

extern void Hard_thresholding(double *qk,  STRUCT_INPUT *Input);

/*----------------------------------------------------------------------------*/

#endif /* FISTA_h */
