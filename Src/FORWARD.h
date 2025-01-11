
#ifndef FORWARD_h
#define FORWARD_h

/*----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdbool.h>
#include <complex.h>


#include "LL04.h"
#include "ALLOCATION.h"
#include "CONSTANT.h"
#include "COLLISION.h"
#include "RANDOM_NUM.h"
#include "READ_INPUT.h"
#include "READ_ATOM.h"
#include "READ_ATMO.h"
#include "ERROR.h"
#include "MPI_CTRL.h"
#include "SEE.h"
#include "FCTSG.h"
#include "INIT.h"
#include "LU.h"

/*----------------------------------------------------------------------------*/

extern int FORWARD(STRUCT_SYN *Syn, STRUCT_ATOM *Atom, \
    STRUCT_INPUT *Input, STR_FCTSG *fctsg, STRUCT_MPI *Mpi, bool sparse);

extern int Copy_Par(STRUCT_SYN *Syn, STRUCT_INPUT *Input, STRUCT_MPI *Mpi);

extern int Grid2Pixel(STRUCT_SYN *Syn, STRUCT_ATOM *Atom, \
    STRUCT_INPUT *Input, STRUCT_OUT *Output, STRUCT_MPI *Mpi);

/*----------------------------------------------------------------------------*/

#endif /* FORWARD_h */
