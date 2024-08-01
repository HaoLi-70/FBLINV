
#ifndef FORWARD_h
#define FORWARD_h

/*----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdbool.h>
#include <complex.h>


#include "LL04.h"
#include "ALLOCATION.h"
#include "CONSTANT.h"
#include "RANDOM_NUMBER.h"
#include "READ_INPUT.h"
#include "READ_ATOM.h"
#include "READ_ATMO.h"
#include "ERROR.h"
#include "MPI_CONTROL.h"
#include "SEE.h"
#include "FCTSG.h"
#include "INIT.h"
#include "LU.h"

/*----------------------------------------------------------------------------*/

extern int FORWARD(STRUCT_SYN *Syn, STRUCT_ATOM *Atom, \
    STRUCT_INPUT *Input, STR_FCTSG *fctsg, STRUCT_MPI *Mpi, bool sparse);

extern int Copy_Par(STRUCT_SYN *Syn, STRUCT_INPUT *Input, STRUCT_MPI *Mpi);

extern int SFB_RECONSTRUTION(STRUCT_INPUT *Input, STRUCT_SYN *Syn, \
    STRUCT_MPI *Mpi);

extern int SFB_RECONSTR_DELTA(STRUCT_INPUT *Input, STRUCT_SYN *Syn, \
    int icoeff, STRUCT_MPI *Mpi);

extern int Grid2Pixel(STRUCT_SYN *Syn, STRUCT_ATOM *Atom, \
    STRUCT_INPUT *Input, STRUCT_OUT *Output, STRUCT_MPI *Mpi);

extern void Coeff2NLM(STRUCT_INPUT *Input, int icoeff, int *ipara, \
    int *in, int *il, int *im, int *real);

extern void NLM2Coeff(STRUCT_INPUT *Input, int ipara, int in, \
    int il, int im, int real, int *icoeff);

/*----------------------------------------------------------------------------*/

#endif /* FORWARD_h */
