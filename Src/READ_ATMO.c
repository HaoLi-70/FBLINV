
#include "READ_ATMO.h"

/*----------------------------------------------------------------------------*/

    /*######################################################################
    
      revision log:
        11 Aug. 2023.
    
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

    char tmp[5];
    fread(tmp,sizeof(char),4,fa);
    tmp[4] = '\0';

    if(strcmp(tmp, "fbmd") != 0){
      Error(enum_error, routine_name, "not a model file");
    }

    fread(&(Atmo->nR),sizeof(int),1,fa);
    if(Atmo->nR<0||Atmo->nR>1e4){
      Error(enum_error, routine_name, "nR error! \n");
    }

    fread(&(Atmo->nTheta),sizeof(int),1,fa);
    if(Atmo->nTheta<0||Atmo->nTheta>1e4){
      Error(enum_error, routine_name, "nR error! \n");
    }

    fread(&(Atmo->nPhi),sizeof(int),1,fa);
    if(Atmo->nPhi<0||Atmo->nPhi>1e4){
      Error(enum_error, routine_name, "nR error! \n");
    }

    Atmo->R = (double *)malloc(Atmo->nR*sizeof(double));
    Atmo->Theta = (double *)malloc(Atmo->nTheta*sizeof(double));
    Atmo->Phi = (double *)malloc(Atmo->nPhi*sizeof(double));

    fread(Atmo->R,sizeof(double),Atmo->nR,fa);
    fread(Atmo->Theta,sizeof(double),Atmo->nTheta,fa);
    fread(Atmo->Phi,sizeof(double),Atmo->nPhi,fa);

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

    fread(Atmo->T[0][0],sizeof(double),num,fa);

    fread(&(Atmo->rhotype),sizeof(int),1,fa);
    if(Atmo->rhotype<0||Atmo->rhotype>3){
      Error(enum_error, routine_name, "rho type error! \n");
    }
    fread(Atmo->rho[0][0],sizeof(double),num,fa);

    fread(&(Atmo->btype),sizeof(int),1,fa);
    if(Atmo->btype<0||Atmo->btype>1){
      Error(enum_error, routine_name, "B type error! \n");
    }
    fread(Atmo->B1[0][0],sizeof(double),num,fa);
    fread(Atmo->B2[0][0],sizeof(double),num,fa);
    fread(Atmo->B3[0][0],sizeof(double),num,fa);

    if(fread(&(Atmo->vtype),sizeof(int),1,fa)>0){
      if(Atmo->vtype>=0){
        Atmo->V1 = TENSOR_DBL(0, Atmo->nR-1, 0, Atmo->nTheta-1, 0, \
                              Atmo->nPhi-1, false);
        Atmo->V2 = TENSOR_DBL(0, Atmo->nR-1, 0, Atmo->nTheta-1, 0, \
                                  Atmo->nPhi-1, false);                      
        Atmo->V3 = TENSOR_DBL(0, Atmo->nR-1, 0, Atmo->nTheta-1, 0, \
                                Atmo->nPhi-1, false);

        fread(Atmo->V1[0][0],sizeof(double),num,fa);
        fread(Atmo->V2[0][0],sizeof(double),num,fa);
        fread(Atmo->V3[0][0],sizeof(double),num,fa);

      }
    }else{
      Atmo->vtype = -1;
    }

    fclose(fa);
    
    return;
}

/*----------------------------------------------------------------------------*/

extern void WRITE_ATMO(char *Filename, STRUCT_ATMO *Atmo){
  
    /*######################################################################
      Purpose:
        write the coronal model.
      Record of revisions:
        16 Jan. 2024.
      Input parameters:
        Filename, path to the mode atmosphere.
        Atmo, a structure with the coronal model.
    ######################################################################*/

    FILE *fa = fopen(Filename, "wb");

    char tmp[5] = "fbmd";
    fwrite(tmp, sizeof(char), 4, fa);

    fwrite(&(Atmo->nR),sizeof(int),1,fa);
    fwrite(&(Atmo->nTheta),sizeof(int),1,fa);
    fwrite(&(Atmo->nPhi),sizeof(int),1,fa);

    fwrite(Atmo->R,sizeof(double),Atmo->nR,fa);
    fwrite(Atmo->Theta,sizeof(double),Atmo->nTheta,fa);
    fwrite(Atmo->Phi,sizeof(double),Atmo->nPhi,fa);

    int num = Atmo->nR*Atmo->nTheta*Atmo->nPhi;

    fwrite(Atmo->T[0][0],sizeof(double),num,fa);

    fwrite(&(Atmo->rhotype),sizeof(int),1,fa);
    fwrite(Atmo->rho[0][0],sizeof(double),num,fa);

    fwrite(&(Atmo->btype),sizeof(int),1,fa);
    fwrite(Atmo->B1[0][0],sizeof(double),num,fa);
    fwrite(Atmo->B2[0][0],sizeof(double),num,fa);
    fwrite(Atmo->B3[0][0],sizeof(double),num,fa);

    if(Atmo->vtype>=0){
      fwrite(&(Atmo->vtype),sizeof(int),1,fa);
      fwrite(Atmo->V1[0][0],sizeof(double),num,fa);
      fwrite(Atmo->V2[0][0],sizeof(double),num,fa);
      fwrite(Atmo->V3[0][0],sizeof(double),num,fa);
    }

    fclose(fa);
    
    return;
}

/*----------------------------------------------------------------------------*/

extern void FREE_MODEL(STRUCT_ATMO *Atmo){
  
    /*######################################################################
      Purpose:
        free the memory of the coronal model.
      Record of revisions:
        16 Jan. 2024.
      Input parameters:
        Atmo, a structure with the coronal model.
    ######################################################################*/
  
    free(Atmo->R);
    free(Atmo->Theta);
    free(Atmo->Phi);
    FREE_TENSOR_DBL(Atmo->B1, 0, 0, 0);
    FREE_TENSOR_DBL(Atmo->B2, 0, 0, 0);
    FREE_TENSOR_DBL(Atmo->B3, 0, 0, 0);
    FREE_TENSOR_DBL(Atmo->rho, 0, 0, 0);
    FREE_TENSOR_DBL(Atmo->T, 0, 0, 0);
    if(Atmo->vtype>=0){
      FREE_TENSOR_DBL(Atmo->V1, 0, 0, 0);
      FREE_TENSOR_DBL(Atmo->V2, 0, 0, 0);
      FREE_TENSOR_DBL(Atmo->V3, 0, 0, 0);
    }
    free(Atmo);
        
    return;
}

/*----------------------------------------------------------------------------*/
