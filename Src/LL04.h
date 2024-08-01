
#ifndef LL04_h
#define LL04_h

/*----------------------------------------------------------------------------*/

#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include "ALLOCATION.h"
#include "CONSTANT.h"
#include "FCTSG.h"

/*----------------------------------------------------------------------------*/

#define KQ_index(k,q,dk) (dk==2?(k*(k+1)/2+q):(k*(k-1)+q)) 

#define cPlank(c1,c2,T) (c1/(exp(c2/T)-1))

#define Wien_limit(c1,c2,T) (c1*exp(-c2/T))

/*----------------------------------------------------------------------------*/

extern double Geffect(double Gu, double Gl, double Ju, double Jl);

extern double Gfactor(double J, double L, double S);

extern complex double Djmn(double Alpha, double Beta, double Gamma, \
    double J, double M, double N, STR_FCTSG *fctsg);

extern int Rotmat(double Alpha, double Beta, double Gamma, \
    complex double ***Djmn, int Rank);

extern int TKP(complex double ***T);

extern void TKQ(complex double ***T, double Chi, double Theta, \
    double Gamma);

extern void TKQ90(complex double ***T);

extern double Planck(double nu, double T);

extern int Jkq_off(double u1, double u2, double r, double *Jkq);

extern double J00_off(double u1, double u2, double r);

extern double J20_off(double u1, double u2, double r);

extern void Limb_Darkening(double Lambda, double *u1, double *u2);

extern void Thom_Scat_van(double r, double *PAB, double q);

extern void Thom_Scat(double r, double *PAB, double u1, double u2);

/*----------------------------------------------------------------------------*/

#endif /* LL04_h */
