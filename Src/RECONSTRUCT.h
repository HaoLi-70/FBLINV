
#ifndef RECONSTRUCT_h
#define RECONSTRUCT_h

/*----------------------------------------------------------------------------*/

#include <complex.h>

#include "COEFF.h"
#include "READ_ATMO.h"
#include "READ_ATOM.h"
#include "READ_INPUT.h"

/*----------------------------------------------------------------------------*/

extern int SFB_RECONSTRUTION(STRUCT_INPUT *Input, STRUCT_SYN *Syn, \
    STRUCT_MPI *Mpi);

extern int SFB_RECONSTR_DELTA(STRUCT_INPUT *Input, STRUCT_SYN *Syn, \
    int icoeff, STRUCT_MPI *Mpi);

/*----------------------------------------------------------------------------*/

#endif /* RECONSTRUCT_h */
