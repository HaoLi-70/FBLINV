
#ifndef MPI_CTRL_h
#define MPI_CTRL_h

/*----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <time.h>
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

extern void ABORT(void);

extern int Randomseeds(STRUCT_MPI *Mpi);

/*----------------------------------------------------------------------------*/

#endif /* MPI_CTRL_h */
