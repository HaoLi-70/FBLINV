
#include "MPI_CTRL.h"

/*----------------------------------------------------------------------------*/

    /*######################################################################
     
     revision log:

      30 Otc. 2024
          --- update: reduced the RAM cost in Randomseeds.
     
    ######################################################################*/

/*----------------------------------------------------------------------------*/

int pid;

/*----------------------------------------------------------------------------*/

extern void CONTROL(void){

    /*######################################################################
      Purpose:
        controls if any CPU has crashed and stops if needed..
      Record of revisions:
        8 Sept. 2021
      Input parameters:
        .
    ######################################################################*/
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    int a1 = 0, a2 = 0;
    
    MPI_Allreduce(&a1, &a2, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
    
    if(!a2) return;
    
    MPI_Finalize();
    
    return;
}

/*----------------------------------------------------------------------------*/

extern void ABORT(void){

    /*######################################################################
      Purpose:
        stop the code.
      Record of revisions:
        8 Sept. 2021
      Input parameters:
        .
      Return:
        return the current conts.
    ######################################################################*/
    
    MPI_Finalize();

    exit(0);
}

/*----------------------------------------------------------------------------*/

extern int Randomseeds(STRUCT_MPI *Mpi){

    /*######################################################################
      Purpose:
        initialize the mpi structure.
      Record of revisions:
        30 Otc. 2024
      Input parameters:
        Mpi, a structure save the information for mpi.
      Output parameters:
        Mpi, the structure save the information for mpi.
    ######################################################################*/
    
    int ipid;

    Mpi->idum = (long *)malloc(sizeof(long));
    long *tmp = (long *)malloc(sizeof(long)*Mpi->nprocs);

    if (Mpi->rank == 0){
      srand((unsigned int)time(NULL));
      srand(rand());  
      for(ipid=0; ipid<Mpi->nprocs; ipid++){  
        tmp[ipid] = -rand();
      }
    }
    MPI_Bcast(tmp, Mpi->nprocs, MPI_LONG, 0, MPI_COMM_WORLD);

    *(Mpi->idum) = tmp[Mpi->rank];
    free(tmp);

    return 0;
}

/*----------------------------------------------------------------------------*/