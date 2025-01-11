
#include "FISTA.h"

/*----------------------------------------------------------------------------*/

    /*######################################################################
 
     revision log:

      30 Otc. 2024
          --- update: moved subroutine array2coeff and coeff2array to 
                      COEFF.c
 
    ######################################################################*/

/*----------------------------------------------------------------------------*/

extern void FISTA(STRUCT_INPUT *Input, STRUCT_ATOM *Atom, STRUCT_SYN *Syn, \
    STRUCT_OUT *Output, STRUCT_OBSERVATION *Obs, STR_FCTSG *fctsg, \
    STRUCT_MPI *Mpi){

    /*######################################################################
      Purpose:
        FISTA.
      Record of revisions:
        1 Sept. 2021.
      Input parameters:
        Input, a structure with the input information.
        Atom, a structure with the atomic information.
        Syn, a structure with forward synthesis.
        fctsg, a structure with factoral, signs.
        Mpi, a structure with the Mpi information.
      Output parameters:
    ######################################################################*/
    
    double *qk, *q, *yk, *yk1, *Derivative;
    double chisq, chisq_new;
    int ipara, i, k = 1, k_restart = 0;
    q = (double *)malloc(sizeof(double)*Input->ncoeff);
    qk = (double *)malloc(sizeof(double)*Input->ncoeff);
    yk = (double *)malloc(sizeof(double)*Input->ncoeff);
    yk1 = (double *)malloc(sizeof(double)*Input->ncoeff);
    Derivative = (double *)malloc(sizeof(double)*Input->ncoeff);

    coeff2array(q, Input);
    Hard_thresholding(q, Input);
    for(i=0; i<Input->ncoeff; i++){
      yk[i] = q[i];
    }

    array2coeff(yk,  Input);
    
    SFB_RECONSTRUTION(Input, Syn, Mpi);

    FORWARD(Syn, Atom, Input, fctsg, Mpi, false);

    Copy_Par(Syn, Input, Mpi);

    Grid2Pixel(Syn, Atom, Input, Output, Mpi);

    if(Mpi->rank == 0) chisq = LOSS_FUNCTION(Obs, Output, Input);

    if(Mpi->rank == 0) WRITE_SYNTHESIS("./rsyn_init", Input, \
        Atom, Output, Syn);

    if(Mpi->rank == 0) WRITE_COEFF("./rcoeff_init", Input->Para, 1, \
        0, Input->rmax);

    MPI_Bcast(&chisq, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if(Mpi->rank == 0) fprintf(stderr, "chisq = %e \n",chisq);

    if(Input->Niter<=0) goto jump;
    /*
    MPI_Bcast(q, Input->ncoeff, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(qk, Input->ncoeff, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(yk, Input->ncoeff, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(yk1, Input->ncoeff, MPI_DOUBLE, 0, MPI_COMM_WORLD);
*/

    bool converge = false;
    double tk = 1, tk1, sum = 0;


    /*
    for(i=0; i<Input->ncoeff; i++){
yk[i] = i;
    }
array2coef(yk,  Input);
return;
       */
    do{  
      if(Mpi->rank == 0) fprintf(stderr, "iter = %d\n", k);
        
      for(i=0; i<Input->ncoeff; i++){

        ipara = SFB_RECONSTR_DELTA(Input, Syn, i, Mpi);

        if(ipara ==0) Ion_fraction(Atom, Syn, Input, Mpi);

        FORWARD(Syn, Atom, Input, fctsg, Mpi, true);

        Grid2Pixel(Syn, Atom, Input, Output, Mpi);

        if(Mpi->rank == 0) chisq_new = LOSS_FUNCTION(Obs, Output, Input);
        MPI_Bcast(&chisq_new, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        Derivative[i] = (chisq_new - chisq)/Input->perturb[ipara];
        qk[i] = yk[i]-Input->step_size*Derivative[i];

        if(Mpi->rank == 0 && i%100 == 99) \
            fprintf(stderr, "%d %e %e\n",i,yk[i],qk[i]);

      }
        
      if(Mpi->rank == 0) fprintf(stderr, "check point 1 \n");

      Hard_thresholding(qk, Input);
      tk1 = 0.5+0.5*sqrt(1.0+4.0*tk*tk);
      for(i=0; i<Input->ncoeff; i++){
        yk1[i] = qk[i]+(tk-1.0)/tk1*(qk[i]-q[i]);
	    }

      sum = 0;
      for(i=0; i<Input->ncoeff; i++){
        sum += Derivative[i]*(qk[i]-q[i]);
      }
        
      for(i=0; i<Input->ncoeff; i++){
        q[i] = qk[i];
      }
	    //if(Mpi->rank == 0) fprintf(stderr,"sum = %e \n",sum);
      if(sum > 0){
        tk = 1.0;
        for(i=0; i<Input->ncoeff; i++){
          yk[i] = qk[i];
        }
	      k_restart = 0;
      }else{
        tk = tk1;
        for(i=0; i<Input->ncoeff; i++){
          yk[i] = yk1[i];
        }
	      k_restart++;
      }
      MPI_Bcast(q, Input->ncoeff, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(yk, Input->ncoeff, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      if(Mpi->rank == 0) fprintf(stderr, "sum %e \n", sum);

      array2coeff(yk,  Input);

      SFB_RECONSTRUTION(Input, Syn, Mpi);

      FORWARD(Syn, Atom, Input, fctsg, Mpi, false);

      Copy_Par(Syn, Input, Mpi);

      Grid2Pixel(Syn, Atom, Input, Output, Mpi);

      if(Mpi->rank == 0) chisq_new = LOSS_FUNCTION(Obs, Output, Input);

      MPI_Bcast(&chisq_new, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      if(Mpi->rank == 0) fprintf(stderr, "chisq = %e %e \n", chisq_new, chisq);
	    if((chisq-chisq_new)/chisq < 3e-3 || chisq_new > chisq) converge = true;
      chisq = chisq_new;
      k++;
      
    } while (!converge && k< Input->Niter);

jump:
   if(Mpi->rank == 0) WRITE_SYNTHESIS("./rsyn", Input, \
        Atom, Output, Syn);

    if(Mpi->rank == 0) WRITE_COEFF("./rcoeff", Input->Para, 1, \
        0, Input->rmax);

    free(q);
    free(qk);
    free(yk);
    free(yk1);
    free(Derivative);

    return;
}

/*----------------------------------------------------------------------------*/

extern double LOSS_FUNCTION(STRUCT_OBSERVATION *Obs, STRUCT_OUT *Output, \
    STRUCT_INPUT *Input){
  
    /*######################################################################
      Purpose:
        compute the loss function.
      Record of revisions:
        21 Aug. 2023.
      Input parameters:
        Obs, structure with obbservation data.
        Output, structure with synthesis.
        Input, a structure with the inversion information.
      Return:
        return the chisq.
    ######################################################################*/
    

    double Chisq = 0;
    /*
    int i, j;
        
    for(i=0; i<Obs->npixels; i++){
      for(j=0; j<Obs->num; j++){
        Chisq += (Output->syntot[i][j]-Obs->Data[i][j]) \
            *(Output->syntot[i][j]-Obs->Data[i][j]) \
            *Output->norm[i]*Input->weight[j]*Input->weight[j];
      }
    }
    */
    return Chisq*1e14;
}

/*----------------------------------------------------------------------------*/

extern void Hard_thresholding(double *qk,  STRUCT_INPUT *Input){
    
    /*######################################################################
      Purpose:
        nonlinear hard thresholding operator.
      Record of revisions:
        9 Sept. 2021.
      Input parameters:
        qk, the coefficient array.
        Input, a structure with the inversion information.
      Output parameters:
        qk, the coefficient array..
    ######################################################################*/
    
    int icoeff, ipara, start=0;
    double *qk_abs = (double *)malloc(sizeof(double)*Input->ncoeff);
    int *indx = (int *)malloc(sizeof(int)*Input->ncoeff);

    for(icoeff=0; icoeff<Input->ncoeff; icoeff++){
      qk_abs[icoeff] = fabs(qk[icoeff]);
    }
        
    for(ipara=0; ipara<4; ipara++){
      if(!Input->Para[ipara].invt) continue;
      qsort_index(start, start+Input->Para[ipara].NLsquare-1, qk_abs, indx);
      for(icoeff=start; icoeff<start+Input->Para[ipara].NLsquare \
          -Input->Para[ipara].threshold_n; icoeff++){
        qk[indx[icoeff]]=0;
      }
      start += Input->Para[ipara].NLsquare;
    }

    free(qk_abs);
    free(indx);
    
    return;
}

/*----------------------------------------------------------------------------*/

