
#ifndef MPI_CONTROL_h
#define MPI_CONTROL_h

/*----------------------------------------------------------------------------*/

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdbool.h>
#include <mpi.h>

/*----------------------------------------------------------------------------*/

extern int pid;

/*----------------------------------------------------------------------------*/

typedef struct Struct_Mpi{

    // number of processor, number of grids, grid indexes, and the rank
    int nprocs, ngrids, gindx[2], rank;

    // random seeds
    long *idum;

}STRUCT_MPI;

extern void CONTROL(void);

/*----------------------------------------------------------------------------*/

#endif /* MPI_CONTROL_h */
