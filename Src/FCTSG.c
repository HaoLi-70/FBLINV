
#include "FCTSG.h"

/*----------------------------------------------------------------------------*/

    /*######################################################################
     
      revision log:

        8 Aug. 2023.
     
    ######################################################################*/

/*----------------------------------------------------------------------------*/

extern int INIT_FCTSG(STR_FCTSG *fctsg){

    /*######################################################################
      Purpose:
        compute the facorial and the sign.
      Record of revisions:
        8 Aug. 2023
      Input parameters:
        fctsg, the structure.
      Output parameters:
        fctsg, the structure.
    ######################################################################*/

    if(fctsg->nmax>171){
      fctsg->nmax = 171;
      Error(enum_warning, "INIT_FCTSG", "nmax is large than 171.");
    }else if(fctsg->nmax<2){
      fctsg->nmax = 2;
      Error(enum_warning, "INIT_FCTSG", "nmax is smaller than 2.");
    }

    fctsg->J3 = (STR_6D *)malloc(sizeof(STR_6D));
    fctsg->J6 = (STR_6D *)malloc(sizeof(STR_6D));
    fctsg->J9 = (STR_9D *)malloc(sizeof(STR_9D));
    fctsg->J3->dat = NULL;
    fctsg->J6->dat = NULL;
    fctsg->J9->dat = NULL;

    fctsg->fct = (double *)malloc((fctsg->nmax+1)*sizeof(double));
    fctsg->fct2 = (double *)malloc((fctsg->nmax+1)*sizeof(double));
    fctsg->sg = (int *)VECTOR(-fctsg->nmax,fctsg->nmax,enum_int,false);

    fctsg->sg[0] = 1;
    fctsg->sg[1] = -1;
    fctsg->sg[-1] = -1;
    fctsg->fct[0] = 1.0;
    fctsg->fct[1] = 1.0;
    fctsg->fct2[0] = 1.0;
    fctsg->fct2[1] = 1.0;

    int i;
    for(i=2;i<=fctsg->nmax;i++){
      fctsg->fct[i] = fctsg->fct[i-1]*i;
      fctsg->sg[i] = -fctsg->sg[i-1];
      fctsg->sg[-i] = fctsg->sg[i];
      fctsg->fct2[i] = fctsg->fct2[i-2]*i;
    }

    return 0;
}

/*----------------------------------------------------------------------------*/

extern int FREE_FCTSG(STR_FCTSG *fctsg){

    /*######################################################################
      Purpose:
        free the memory in sctsg.
      Record of revisions:
        8 Aug. 2023
      Input parameters:
        fctsg, the structure.
    ######################################################################*/

    if(fctsg->memo){

      /////TBC
      free(fctsg->J3);
      free(fctsg->J6);
      free(fctsg->J9);
    }else{
      free(fctsg->J3);
      free(fctsg->J6);
      free(fctsg->J9);
    }

    free(fctsg->fct);
    free(fctsg->fct2);
    FREE_VECTOR(fctsg->sg, -fctsg->nmax, enum_int);

    free(fctsg);

    return 0;

}

/*----------------------------------------------------------------------------*/
