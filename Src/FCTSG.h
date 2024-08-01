
#ifndef FCTSG_h
#define FCTSG_h

/*----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "ALLOCATION.h"
#include "ERROR.h"
#include "MEMOR.h"


typedef struct Struct_Factor_Sign{

    bool memo;

    STR_6D *J3;
    STR_6D *J6;
    STR_9D *J9;

    int nmax;
    int *sg;
    double *fct, *fct2;

}STR_FCTSG;

/*----------------------------------------------------------------------------*/

extern int INIT_FCTSG(STR_FCTSG *fctsg);

extern int FREE_FCTSG(STR_FCTSG *fctsg);

/*----------------------------------------------------------------------------*/

#endif /* ERROR_h */
