
#include "READ_ATMO.h"

/*----------------------------------------------------------------------------*/

    /*######################################################################
    
      revision log:

        30 Otc. 2024
          --- update: moved subroutine FREE_MODEL to FREE.c  
                      moved subroutine WRITE_ATMO to IO.c  
    
    ######################################################################*/

/*----------------------------------------------------------------------------*/

extern void READ_ATMO(char *Filename, STRUCT_ATMO *Atmo){
  
    /*######################################################################
      Purpose:
        read the coronal model.
      Record of revisions:
        16 Jan. 2024.
      Input parameters:
        Filename, path to the mode atmosphere.
      Output parameters:
        Atmo, a structure with the coronal model.
    ######################################################################*/

    const char *routine_name = "READ_ATMO";

    FILE *fa = fopen(Filename, "rb");

    size_t nsize;

    char tmp[5];
    nsize = fread(tmp,sizeof(char),4,fa);
    tmp[4] = '\0';

    if(strcmp(tmp, "fbmd") != 0){
      Error(enum_error, routine_name, "not a model file");
    }

    nsize = fread(&(Atmo->nR),sizeof(int),1,fa);
    if(Atmo->nR<0||Atmo->nR>1e4){
      Error(enum_error, routine_name, "nR error! \n");
    }

    nsize = fread(&(Atmo->nTheta),sizeof(int),1,fa);
    if(Atmo->nTheta<0||Atmo->nTheta>1e4){
      Error(enum_error, routine_name, "nR error! \n");
    }

    nsize = fread(&(Atmo->nPhi),sizeof(int),1,fa);
    if(Atmo->nPhi<0||Atmo->nPhi>1e4){
      Error(enum_error, routine_name, "nR error! \n");
    }

    Atmo->R = (double *)malloc(Atmo->nR*sizeof(double));
    Atmo->Theta = (double *)malloc(Atmo->nTheta*sizeof(double));
    Atmo->Phi = (double *)malloc(Atmo->nPhi*sizeof(double));

    nsize = fread(Atmo->R,sizeof(double),Atmo->nR,fa);
    nsize = fread(Atmo->Theta,sizeof(double),Atmo->nTheta,fa);
    nsize = fread(Atmo->Phi,sizeof(double),Atmo->nPhi,fa);

    Atmo->T = TENSOR_DBL(0, Atmo->nR-1, 0, Atmo->nTheta-1, 0, \
                          Atmo->nPhi-1, false);

    Atmo->rho = TENSOR_DBL(0, Atmo->nR-1, 0, Atmo->nTheta-1, 0, \
                                Atmo->nPhi-1, false);

    Atmo->B1 = TENSOR_DBL(0, Atmo->nR-1, 0, Atmo->nTheta-1, 0, \
                          Atmo->nPhi-1, false);
    Atmo->B2 = TENSOR_DBL(0, Atmo->nR-1, 0, Atmo->nTheta-1, 0, \
                              Atmo->nPhi-1, false);                      

    Atmo->B3 = TENSOR_DBL(0, Atmo->nR-1, 0, Atmo->nTheta-1, 0, \
                            Atmo->nPhi-1, false);

    int num = Atmo->nR*Atmo->nTheta*Atmo->nPhi;

    nsize = fread(Atmo->T[0][0],sizeof(double),num,fa);

    nsize = fread(&(Atmo->rhotype),sizeof(int),1,fa);
    if(Atmo->rhotype<0||Atmo->rhotype>3){
      Error(enum_error, routine_name, "rho type error! \n");
    }
    nsize = fread(Atmo->rho[0][0],sizeof(double),num,fa);

    nsize = fread(&(Atmo->btype),sizeof(int),1,fa);
    if(Atmo->btype<0||Atmo->btype>1){
      Error(enum_error, routine_name, "B type error! \n");
    }
    nsize = fread(Atmo->B1[0][0],sizeof(double),num,fa);
    nsize = fread(Atmo->B2[0][0],sizeof(double),num,fa);
    nsize = fread(Atmo->B3[0][0],sizeof(double),num,fa);

    if(fread(&(Atmo->vtype),sizeof(int),1,fa)>0){
      if(Atmo->vtype>=0){
        Atmo->V1 = TENSOR_DBL(0, Atmo->nR-1, 0, Atmo->nTheta-1, 0, \
                              Atmo->nPhi-1, false);
        Atmo->V2 = TENSOR_DBL(0, Atmo->nR-1, 0, Atmo->nTheta-1, 0, \
                                  Atmo->nPhi-1, false);                      
        Atmo->V3 = TENSOR_DBL(0, Atmo->nR-1, 0, Atmo->nTheta-1, 0, \
                                Atmo->nPhi-1, false);

        nsize = fread(Atmo->V1[0][0],sizeof(double),num,fa);
        nsize = fread(Atmo->V2[0][0],sizeof(double),num,fa);
        nsize = fread(Atmo->V3[0][0],sizeof(double),num,fa);

      }
    }else{
      Atmo->vtype = -1;
    }

    fclose(fa);
    
    return;
}

/*----------------------------------------------------------------------------*/
