
#ifndef COEFF_h
#define COEFF_h

/*----------------------------------------------------------------------------*/

#include "READ_INPUT.h"

/*----------------------------------------------------------------------------*/

extern void Coeff2NLM(STRUCT_INPUT *Input, int icoeff, int *ipara, \
    int *in, int *il, int *im, int *real);

extern void NLM2Coeff(STRUCT_INPUT *Input, int ipara, int in, \
    int il, int im, int real, int *icoeff);

extern void array2coeff(double *qk,  STRUCT_INPUT *Input);

extern void coeff2array(double *qk,  STRUCT_INPUT *Input);

/*----------------------------------------------------------------------------*/

#endif /* COEFF_h */
