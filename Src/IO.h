
#ifndef IO_h
#define IO_h

#include <stdio.h>
#include <stdbool.h>
#include <string.h>

#include "ALLOCATION.h"
#include "READ_INPUT.h"
#include "STR.h"

/*----------------------------------------------------------------------------*/

extern int READ_ZEROS(char Filename[], STRUCT_INPUT *Input, STRUCT_MPI *Mpi);

extern int WRITE_ZEROS(char *Filename, int Zero_L, int Zero_N, \
    enum boundary_condition type);

extern int READ_OBSERVATION(STRUCT_INPUT *Input, STRUCT_OBSERVATION *Obs);

extern int WRITE_SYNTHESIS(char *Filename, STRUCT_INPUT *Input, \
    STRUCT_ATOM *Atom, STRUCT_OUT *Output, STRUCT_SYN *Syn);

extern void READ_COEFF(char *Filename, STRUCT_SFB_COEFF *Para, double *Rlim);

extern void WRITE_COEFF(char *Filename, STRUCT_SFB_COEFF *Para, int btype, \
    int rhotype, double Rlim);

extern int FREE_COEFF(STRUCT_SFB_COEFF *Para);

/*----------------------------------------------------------------------------*/


#endif /* IO_h */
