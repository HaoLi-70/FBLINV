#ifndef INTERPOL_h
#define INTERPOL_h

/*----------------------------------------------------------------------------*/

#include "READ_ATMO.h"
#include "READ_ATOM.h"
#include "READ_INPUT.h"

/*----------------------------------------------------------------------------*/

extern double Interpol_Linear3D(double ***Data, int i, int j, int k, \
    double ir, double jr, double kr);

extern void Get_Para(STRUCT_ATMO *Atmo, STRUCT_ATOM *Atom, \
    STRUCT_SYN *Syn, STRUCT_INPUT *Input, STRUCT_MPI *Mpi);
    
extern int Ion_fraction(STRUCT_ATOM *Atom, STRUCT_SYN *Syn, \
    STRUCT_INPUT *Input, STRUCT_MPI *Mpi);
    
/*----------------------------------------------------------------------------*/

#endif /* INTERPOL_h */