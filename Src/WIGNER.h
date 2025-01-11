
#ifndef WIGNER_h
#define WIGNER_h

/*----------------------------------------------------------------------------*/

#include <math.h>
#include "FCTSG.h"
#include "MEMOR.h"

/*----------------------------------------------------------------------------*/

extern double WIGNER_3J(double J1, double J2, double J3, double M1, \
    double M2, double M3, STR_FCTSG *fctsg);

extern double WIGNER_6J(double J1, double J2, double J3, double J4, \
    double J5, double J6, STR_FCTSG *fctsg);

extern double WIGNER_9J(double J1, double J2, double J3, double J4, \
    double J5, double J6, double J7, double J8, double J9, \
    STR_FCTSG *fctsg);

/*----------------------------------------------------------------------------*/

#endif /* WIGNER_h */
