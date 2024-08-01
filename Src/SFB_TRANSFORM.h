
#ifndef SFB_TRANSFORM_h
#define SFB_TRANSFORM_h

/*----------------------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "SORT.h"
#include "CONSTANT.h"
#include "ALLOCATION.h"
#include "SPECIAL_FUNCTIONS.h"

/*----------------------------------------------------------------------------*/

typedef struct Struct_Model{
    double ***Data;
    int nR, nTheta, nPhi;
    double *R;
    double *Theta;
    double *Phi;
}STRUCT_MODEL;

typedef struct Struct_Order{
    int L_Max;
    int K_Max;
    int R_Max;
}STRUCT_ORDER;

typedef struct Struct_Sphere_Filename{
    char *file_R;
    char *file_THETA;
    char *file_PHI;
    char *file_Data;
}Str_Sph_File;

typedef struct {
    int X;
    int Y;
    int Z;
    double Ratio_X;
    double Ratio_Y;
    double Ratio_Z;
}Position_Struct;

/*----------------------------------------------------------------------------*/

extern void SB_DECOMP_SINGLE(double *Field_Array, double *Radius, \
    int Num_R, double **Zeros, int Kmax, int il, int Int_Flag, \
    enum boundary_condition BC_Flag, double *Coef_SB);

extern void SB_RECONSTRUCT_SINGLE(double *Coef_SB, int Kmax, int il, \
    double **Zeros, double R_Max, double *Radius, int Num_R, \
    double *field);

extern void SHB_DECOMP(STRUCT_MODEL *Model, STRUCT_ORDER *Order, \
    double **Zeros, enum boundary_condition BC_Flag, \
    complex double ***Coeff_klm);




extern void Sparse_Set(STRUCT_ORDER *Order, double Percent, \
    complex double ***Coeff_klm);

extern void SH_DECOMP(STRUCT_MODEL *Model, STRUCT_ORDER *Order, \
    complex double ***Coeff_rlm);

extern void SB_DECOMP(complex double ***Coeff_rlm, STRUCT_ORDER *Order, \
    double **Zeros, double *Radius, int Int_Flag, \
    enum boundary_condition BC_Flag, complex double ***Coeff_klm);

extern void SFB_DECOMP(STRUCT_MODEL *Model, STRUCT_ORDER *Order, \
    double **Zeros, int Int_Flag, enum boundary_condition BC_Flag, \
    complex double ***Coeff_klm);

extern void SB_RECONSTRUCT(STRUCT_ORDER *Order, double **Zeros, \
    double *Radius, complex double ***Coeff_klm, complex double ***Coeff_rlm);

extern void SH_RECONSTRUCT(STRUCT_ORDER *Order, complex double ***Coeff_rlm, \
    STRUCT_MODEL *Model);

extern void SFB_RECONSTRUCT(complex double ***Coeff_klm, STRUCT_ORDER *Order, \
    double **Zeros, STRUCT_MODEL *Model);




extern double Xvalue(double x, int N);

extern double Inl(double x, int n, int l, double *J);

extern double IntgBessel(double x0, double x1, int l, double knl);


extern void SFBde(STRUCT_MODEL *Model, STRUCT_ORDER *Order, \
    double **Zeros, int Int_Flag, enum boundary_condition BC_Flag, \
    complex double ***Coeff_klm);

extern void SFBre(STRUCT_ORDER *Order, complex double ***Coeff_rlm, 
    STRUCT_MODEL *Model, double R_Max, double **Zeros, \
    enum boundary_condition BC_Flag);

/*----------------------------------------------------------------------------*/

#endif /* SFB_TRANSFORM_h */
