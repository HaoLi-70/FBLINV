#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>

#include <mpi.h>


#include <complex.h>

#include "ALLOCATION.h"
#include "READ_ATMO.h"
#include "READ_ATOM.h"
#include "READ_INPUT.h"
#include "TRATE.h"


#include "FCTSG.h"

#include "TIME_PRINT.h"
#include "INIT.h"
#include "SEE.h"
#include "FORWARD.h"
#include "IO.h"
#include "ERROR.h"
#include "FISTA.h"
#include "FREE.h"
#include "INTERPOL.h"

#define Path_Input "./input.inv"

/*----------------------------------------------------------------------------*/

int main(int argc, char *argv[]){
  
    /*######################################################################
      Purpose:
        Read the input file for the forbidden line calculation.
      Record of revisions:
        10 Sept. 2021.
      Input parameters:
        Filename[], the input file.
      Output parameters:
        Input, a structure saved the input information.
        Input, a structure saved the inversion information.
        Output, a structure saved the output information.
        Mpi, a structure saved the Mpi information.

        X point to the observer
        Z vertical
        Y horizontal

     ######################################################################*/
    

    /*---------- Initialize the MPI Execution Environment ----------*/
    
    MPI_Init(&argc,&argv);
    
    STRUCT_MPI *Mpi;
    Mpi = (STRUCT_MPI *)malloc(sizeof(STRUCT_MPI));
    MPI_Comm_rank(MPI_COMM_WORLD, &(Mpi->rank));
    MPI_Comm_size(MPI_COMM_WORLD, &(Mpi->nprocs));
    pid = Mpi->rank;
    
    /*---------- Begin to calculate the Running Time ----------*/
    
    if(Mpi->rank == 0) Time_Print();

    STRUCT_INPUT *Input = (STRUCT_INPUT *)malloc(sizeof(STRUCT_INPUT));

    char *Path = Path_Input;
    if(argc>=2){
      Path = argv[1];
    }

    if(!FILE_EXIST(Path)){
      Error(enum_error, "main", "input doesn't exist!");
    }

    RDINPUT(Path, Input, Mpi);
    int iatom, ipara;

    if(Input->Mode<4){

      STRUCT_ATOM *Atom = (STRUCT_ATOM *)malloc(Input->Natom \
          *sizeof(STRUCT_ATOM));
      STR_FCTSG *fctsg = (STR_FCTSG *)malloc(sizeof(STR_FCTSG));
      STRUCT_SYN *Syn = (STRUCT_SYN *)malloc(sizeof(STRUCT_SYN));
      STRUCT_OUT *Output = (STRUCT_OUT *)malloc(sizeof(STRUCT_OUT));
      STRUCT_ATMO *Atmo = (STRUCT_ATMO *)malloc(sizeof(STRUCT_ATMO));
      READ_ATMO(Input->Path_Atmo, Atmo);

      MPI_Barrier(MPI_COMM_WORLD);
      if(Mpi->rank == 0) fprintf(stderr, " Atmosphere model read \n");

      for(iatom = 0; iatom < Input->Natom; iatom++){
        RDATOM(Input->Path_Atom[iatom], Input, Atom+iatom);
      }

      if(Mpi->rank == 0) fprintf(stderr, " Atoms read \n");

      if(Input->Mode==0){

        STRUCT_OBSERVATION *Obs = (STRUCT_OBSERVATION *)malloc( \
            sizeof(STRUCT_OBSERVATION));
        READ_OBSERVATION(Input, Obs);
        if(Mpi->rank == 0) fprintf(stderr, " obervation read \n");

        Init(Input, Atom, Syn, fctsg, Output, Mpi);
        if(Mpi->rank == 0) fprintf(stderr, " Initalized \n");

        if(Obs->npixels!=Syn->npixels){
          Error(enum_error, "main", "pixel number mismatch");
        }

        for(iatom = 0; iatom < Input->Natom; iatom++){
          if(Atom[iatom].TwoLv) continue;
          if(Mpi->rank == 0) fprintf(stderr, " computing rates for "\
              "atom %d \n", iatom);
          TRates(Atom+iatom, fctsg, Input);
        }

        Get_Para(Atmo, Atom, Syn, Input, Mpi);
        if(Mpi->rank == 0) fprintf(stderr, " parameter got \n");
        FREE_MODEL(Atmo);


        STRUCT_SFB_COEFF *Para = (STRUCT_SFB_COEFF *)malloc(4*sizeof \
            (STRUCT_SFB_COEFF));

        READ_COEFF(Input->Path_coeff, Para, &(Input->Rlim));

        INIT_COEFF(Para, Input);

        FISTA(Input, Atom, Syn, Output, Obs, fctsg, Mpi);

      }else{

        Init(Input, Atom, Syn, fctsg, Output, Mpi);

        

        if(Mpi->rank == 0) fprintf(stderr, " Initalized \n");

        for(iatom = 0; iatom < Input->Natom; iatom++){
          if(Atom[iatom].TwoLv) continue;
          if(Mpi->rank == 0) fprintf(stderr, " computing rates for "\
              "atom %d \n", iatom);
          TRates(Atom+iatom, fctsg, Input);
        }

        Get_Para(Atmo, Atom, Syn, Input, Mpi);
        if(Mpi->rank == 0) fprintf(stderr, " parameter got \n");

        FREE_MODEL(Atmo);


//    MPI_Finalize();
//    return 0;
        if(Mpi->rank == 0) fprintf(stderr, " model freed \n");

        FORWARD(Syn, Atom, Input, fctsg, Mpi, false);

        if(Mpi->rank == 0) fprintf(stderr, " synthesis finished \n");


        Grid2Pixel(Syn, Atom, Input, Output, Mpi);

        if(Mpi->rank == 0) WRITE_SYNTHESIS(Input->Path_Output, Input, \
            Atom, Output, Syn);

        if(Mpi->rank == 0) fprintf(stderr," synthesis written \n");


        if(Input->OutputPara) WRITE_PARA(Input, Syn, Atom, Mpi);

      }

      FREE_FCTSG(fctsg);
      FREE_GRIDS(Input, Syn, Mpi);
      FREE_OUTPUT(Output, Input);
      FREE_ATOM(Atom, Input->Natom);

    }else{

      if(Mpi->rank == 0){

        STRUCT_ATMO *Atmo = (STRUCT_ATMO *)malloc(sizeof(STRUCT_ATMO));
        STRUCT_ORDER *Order = (STRUCT_ORDER *)malloc(sizeof(STRUCT_ORDER));
        STRUCT_MODEL *Model = (STRUCT_MODEL *)malloc(sizeof(STRUCT_MODEL));

        fprintf(stderr, " reading Atmosphere model\n");
        READ_ATMO(Input->Path_Atmo, Atmo);
        fprintf(stderr, " reading zeros\n");
        READ_ZEROS(Input->Path_Zero, Input, Mpi);

        Order->R_Max = Atmo->nR;
        Model->R = Atmo->R;
        Model->Theta = Atmo->Theta;
        Model->Phi = Atmo->Phi;
        Model->nR = Atmo->nR;
        Model->nTheta = Atmo->nTheta;
        Model->nPhi = Atmo->nPhi;



        if(Input->Mode==4){

          fprintf(stderr, " Decomposition \n");
          fprintf(stderr,"btype = %d rho type = %d \n", Atmo->btype, \
              Atmo->rhotype);

          for(ipara=0;ipara<4;ipara++){
            fprintf(stderr,"decom %d %d\n",ipara,Input->Para[ipara].invt);

            if(!Input->Para[ipara].invt) continue;
            Order->L_Max = Input->Para[ipara].L;
            Order->K_Max = Input->Para[ipara].N;
            Input->Para[ipara].Coeff = (complex double ***) \
                TENSOR_RHO_CPLX(Order->K_Max, Order->L_Max, true);
        
            switch(ipara){
              case 0:
                Model->Data = Atmo->T;
                break;

              case 1:
                Model->Data = Atmo->rho;
                break;

              case 2:
                Model->Data = Atmo->B2;
                break;

              case 3:
                Model->Data = Atmo->B3;
                break;

              default:
                break;
            }

            SFBde(Model, Order, Input->Zeros, Input->inner_flag, \
                Input->BC, Input->Para[ipara].Coeff);
  
            //SFB_DECOMP(Model, Order, Input->Zeros, Input->inner_flag, 
            //    Input->BC, Input->Para[ipara].Coeff);

          }

          WRITE_COEFF(Input->Path_Output, Input->Para, Atmo->btype, \
              Atmo->rhotype, Model->R[Model->nR-1]);

          fprintf(stderr,"Finish writing coefficients \n");
        
        }else{
          READ_COEFF(Input->Path_coeff, Input->Para, &(Input->Rlim));

          for(ipara=0;ipara<4;ipara++){
            fprintf(stderr,"recon %d %d\n",ipara,Input->Para[ipara].invt);
            if(!Input->Para[ipara].invt) continue;
            Order->L_Max = Input->Para[ipara].L;
            Order->K_Max = Input->Para[ipara].N;

            switch(ipara){
              case 0:
                Model->Data = Atmo->T;
                break;

              case 1:
                Model->Data = Atmo->rho;
                break;

              case 2:
                Model->Data = Atmo->B2;
                break;

              case 3:
                Model->Data = Atmo->B3;
                break;

              default:
                break;
            }

            Sparse_Set(Order, Input->Para[ipara].sparsity, \
                Input->Para[ipara].Coeff);

            SFB_RECONSTRUCT(Input->Para[ipara].Coeff, Order, Input->Zeros, \
                Model);
SFBre(Order, Input->Para[ipara].Coeff, 
    Model, Input->Rlim, Input->Zeros, \
    Input->BC);
          }  
          
          WRITE_ATMO(Input->Path_Output, Atmo);

          fprintf(stderr,"Finish writing the atmosphere \n");
  
        }

        FREE_MODEL(Atmo);
        free(Model);
        free(Order);

      }
    }

    FREE_INPUT(Input);

    if(Mpi->rank == 0) Time_Print();

    free(Mpi);

    MPI_Finalize();
    return 0;

    /*----------                       END                      ----------*/
    
}

/*----------------------------------------------------------------------------*/
