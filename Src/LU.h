
#ifndef LU_H
#define LU_H

/*----------------------------------------------------------------------------*/

#include <complex.h>
#include <math.h>
#include "ALLOCATION.h"

/*----------------------------------------------------------------------------*/

//A small number.
#define TINY_IN_LU 1.0e-20;

/*----------------------------------------------------------------------------*/

extern int ludcmp_dbl(double **a, int n, int *indx);

extern void lubksb_dbl(double **a, int n, int *indx, double *b);

extern void improve_dbl(double **a, double **alud, int n, int *indx, \
    double *b, double *x);

extern int ludcmp_flt(float **a, int n, int *indx);

extern void lubksb_flt(float **a, int n, int *indx, float *b);

extern void improve_flt(float **a, float **alud, int n, int *indx, \
    float *b, float *x);

extern int ludcmp_cplx(complex double **a, int n, int *indx);

extern void lubksb_cplx(complex double **a, int n, int *indx, \
    complex double *b);

extern void improve_cplx(complex double **a, complex double **alud, \
    int n, int *indx, complex double *b, complex double *x);

/*----------------------------------------------------------------------------*/

#endif /* LU_H */
