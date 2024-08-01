
#ifndef MEMOR_h
#define MEMOR_h

/*----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------*/

typedef struct Struct_Scalar{
    double *dat;
}STR_SCALAR;

typedef struct Struct_1D{
    int bounds[2];
    STR_SCALAR *dat;
}STR_1D;

typedef struct Struct_2D{
    int bounds[2];
    STR_1D *dat;
}STR_2D;

typedef struct Struct_3D{
    int bounds[2];
    STR_2D *dat;
}STR_3D;

typedef struct Struct_4D{
    int bounds[2];
    STR_3D *dat;
}STR_4D;

typedef struct Struct_5D{
    int bounds[2];
    STR_4D *dat;
}STR_5D;

typedef struct Struct_6D{
    int bounds[2];
    STR_5D *dat;
}STR_6D;

typedef struct Struct_7D{
    int bounds[2];
    STR_6D *dat;
}STR_7D;

typedef struct Struct_8D{
    int bounds[2];
    STR_7D *dat;
}STR_8D;

typedef struct Struct_9D{
  int bounds[2];
    STR_8D *dat;
}STR_9D;

/*----------------------------------------------------------------------------*/

extern int inster_6d(STR_6D *ptr, double val, int indx1, int indx2, int indx3, \
    int indx4, int indx5, int indx6);

extern int inster_9d(STR_9D *ptr, double val, int indx1, int indx2, int indx3, \
    int indx4, int indx5, int indx6, int indx7, int indx8, int indx9);

extern double *ele6d(STR_6D *ptr, int indx1, int indx2, int indx3, int indx4, 
    int indx5, int indx6);

extern double *ele9d(STR_9D *ptr, int indx1, int indx2, int indx3, int indx4, 
    int indx5, int indx6, int indx7, int indx8, int indx9);

/*----------------------------------------------------------------------------*/

#endif /* ERROR_h */
