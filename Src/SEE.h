
#ifndef SEE_h
#define SEE_h

/*----------------------------------------------------------------------------*/

#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include "CONSTANT.h"
#include "ALLOCATION.h"
#include "LL04.h"
#include "FCTSG.h"
#include "WIGNER.h"
#include "READ_ATOM.h"



/*----------------------------------------------------------------------------*/

#define tiny 1e-10

/*----------------------------------------------------------------------------*/

extern int Transition_Rates(STRUCT_ATOM *Atom, STR_FCTSG *fctsg, \
    STRUCT_INPUT *Input);
    
extern complex double **SEE(STRUCT_ATOM *Atom, STR_FCTSG *fctsg, \
    STRUCT_PARA *Para, STRUCT_INPUT *Input);

/*----------------------------------------------------------------------------*/

#endif /* LL04_h */
