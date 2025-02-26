
#ifndef READ_ATOM_h
#define READ_ATOM_h

/*----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "STR.h"
#include "CONSTANT.h"
#include "ALLOCATION.h"
#include "LL04.h"
#include "ERROR.h"
#include "READ_INPUT.h"
#include "WIGNER.h"

/*----------------------------------------------------------------------------*/

extern void RDATOM(char filename[], STRUCT_INPUT *Input, STRUCT_ATOM *Atom);

/*----------------------------------------------------------------------------*/

#endif /* READ_ATOM_h */
