
#ifndef TRATE_h
#define TRATE_h

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
#include "READ_INPUT.h"

/*----------------------------------------------------------------------------*/

extern int TRates(STRUCT_ATOM *Atom, STR_FCTSG *fctsg, \
    STRUCT_INPUT *Input);
    
/*----------------------------------------------------------------------------*/

#endif /* TRATE_h */