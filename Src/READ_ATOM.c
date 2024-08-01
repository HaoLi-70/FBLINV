
#include "READ_ATOM.h"

/*----------------------------------------------------------------------------*/

    /*######################################################################
     
      revision log:
        22 Aug. 2023.
     
    ######################################################################*/   

/*----------------------------------------------------------------------------*/

extern void RDATOM(char filename[], STRUCT_INPUT *Input, STRUCT_ATOM *Atom){
    
    /*######################################################################
      Purpose:
        read the atomic information from a file used for the forbidden
            line calculation.
      Record of revisions:
        10 Aug. 2023.
      Input parameters:
        Filename[], the input file.
        Input, a structure with the input information.
      Output parameters:
        Input, a structure with the input information.
        Atom, a structure with the atomic information.
    ######################################################################*/
    
    FILE *Fa = fopen(filename, "r");
    
    char lines[Max_Line_Length], parameter[Key_Length];
    char tmp1[Key_Length], tmp2[Key_Length];
    int indx, i, j, au, al, ilevel, itrans, iline, nread, \
        outarry[100];
    double Arates, T = 5800;
 
    Read_line(lines, Fa);

    memmove(Atom->Element, lines, 10);
    //if (Mpi_pid == 0) fprintf(stderr, "%s \n", Atom->Element);
          
    Read_line(lines, Fa);
    Atom->Mass = atof(lines);
   // if (Mpi_pid == 0) fprintf(stderr, "%e \n", Atom->Mass);
  
    Atom->cDopp = sqrt(2.*C_Kb/Atom->Mass/C_m0)/C_c;

    Read_line(lines, Fa);
    Atom->Abund = atof(lines);
    //if (Mpi_pid == 0) fprintf(stderr, "%e \n", Atom->Abundance);
 
    Read_line(lines, Fa);
    sscanf(lines, "%d  %d",&Atom->Nlevel, &Atom->Ntrans);

    Atom->TwoLv = Atom->Nlevel==2?true:false;
        
    Atom->LV = (STRUCT_LV *)malloc(Atom->Nlevel*sizeof(STRUCT_LV));
    Atom->TR = (STRUCT_TRANS *)malloc(Atom->Ntrans*sizeof(STRUCT_TRANS));
    Atom->eqindx = (int *)malloc(Atom->Nlevel*sizeof(int));
    Atom->nEq = -1;

    for(ilevel=0; ilevel<Atom->Nlevel; ilevel++){
      Read_line(lines, Fa);
      nread = sscanf(lines,"%d, %lf, %lf, %lf, %lf", &indx, \
          &Atom->LV[ilevel].Energy, &Atom->LV[ilevel].J, \
          &Atom->LV[ilevel].L, &Atom->LV[ilevel].S);

      if(nread!=5){
        Error(enum_error, "READ_ATOM", "error in reading levels \n");
      }
      if(indx!=ilevel){
        Error(enum_warning, "READ_ATOM", "level index warning \n");
      }

      Atom->LV[ilevel].deg = 2*Atom->LV[ilevel].J+1;
      Atom->LV[ilevel].sqrt_deg = sqrt(Atom->LV[ilevel].deg);

      if(Input->rJmax<Atom->LV[ilevel].J) Input->rJmax=Atom->LV[ilevel].J;
      if(Input->rLmax<Atom->LV[ilevel].L) Input->rLmax=Atom->LV[ilevel].L;
      if(Input->rSmax<Atom->LV[ilevel].S) Input->rSmax=Atom->LV[ilevel].S;

      Atom->LV[ilevel].Energy *= 100;

      Atom->LV[ilevel].Kmax = lround(Atom->LV[ilevel].J*2);
      if(Input->Kdelta==2) Atom->LV[ilevel].Kmax=Atom->LV[ilevel].Kmax/2*2;

      if(Atom->LV[ilevel].Kmax>Input->Kmax) Atom->LV[ilevel].Kmax = \
          Input->Kmax;

      if(Input->Saturated){
        Atom->LV[ilevel].nKQ = 1+Atom->LV[ilevel].Kmax/2;
      }else{
        if(Input->Kdelta==2){
          Atom->LV[ilevel].nKQ = (Atom->LV[ilevel].Kmax/2+1) \
              *(Atom->LV[ilevel].Kmax/2*2+1);
        }else{
          Atom->LV[ilevel].nKQ = (Atom->LV[ilevel].Kmax+1) \
              *(Atom->LV[ilevel].Kmax+1);
        }
      }
      Atom->nEq += Atom->LV[ilevel].nKQ;

      Atom->LV[ilevel].Rho = (complex double **)MATRIX_RHO_CPLX( \
          Atom->LV[ilevel].Kmax, true);

      Atom->LV[ilevel].g = Gfactor(Atom->LV[ilevel].J, \
          Atom->LV[ilevel].L, Atom->LV[ilevel].S);
     
      /*
      if (Mpi_pid == 0) fprintf(stderr, "%d, %lf, %lf, %lf, %lf\n", \
          indx, Atom->LV[ilevel].Energy, Atom->LV[ilevel].J, \
          Atom->LV[ilevel].L, Atom->LV[ilevel].S);
      */
    }

    Atom->eqindx[0] = 0;
    for(ilevel=0; ilevel<Atom->Nlevel-1; ilevel++){
      Atom->eqindx[ilevel+1] = Atom->eqindx[ilevel]+Atom->LV[ilevel].nKQ;
    }

    Atom->Nline = 0;
    if(Input->Maxntran<Atom->Ntrans) Input->Maxntran=Atom->Ntrans;
    for(itrans=0; itrans<Atom->Ntrans; itrans++){
      Read_line(lines, Fa);
      nread = sscanf(lines,"%d %d %lf %s %s", &au, &al, \
          &Arates, tmp1, tmp2);

      if(nread>5||nread<4){
        Error(enum_error, "READ_ATOM", "error in reading transiton \n");
      }

      if(au>=Atom->Nlevel||al<0||au<al){
        Error(enum_error, "READ_ATOM", "index error for the "\
            "transitons \n");
      }

      Atom->TR[itrans].Aul = Arates;
      Atom->TR[itrans].au = au;
      Atom->TR[itrans].al = al;
      Trim(tmp1, 3);
      if(strcmp(tmp1, "E1") == 0 || strcmp(tmp1, "E1") == 0){
        Atom->TR[itrans].M1 = false;
      }else{
        Atom->TR[itrans].M1 = true;
      }

      Atom->TR[itrans].lambda = 1/(Atom->LV[au].Energy \
          -Atom->LV[al].Energy);
      Atom->TR[itrans].nu = C_c/Atom->TR[itrans].lambda;
      Atom->TR[itrans].cplank1 = 2*C_h*Atom->TR[itrans].nu \
          *Atom->TR[itrans].nu*Atom->TR[itrans].nu/C_c/C_c;
      Atom->TR[itrans].cplank2 = C_h*Atom->TR[itrans].nu/C_Kb;
      Atom->TR[itrans].Intens = cPlank(Atom->TR[itrans].cplank1, \
          Atom->TR[itrans].cplank2, T);

      Atom->TR[itrans].Bul = 1/(2*C_h*C_c)*Atom->TR[itrans].Aul \
          *Atom->TR[itrans].lambda*Atom->TR[itrans].lambda \
          *Atom->TR[itrans].lambda;
        //fprintf(stderr, "%e %e %e \n", Atom->TR[itrans].lambda,C_h,C_c);
      Atom->TR[itrans].Blu = Atom->TR[itrans].Bul*(2*Atom->LV[au].J+1) \
          /(2*Atom->LV[al].J+1);

      Limb_Darkening(Atom->TR[itrans].lambda, &Atom->TR[itrans].u1, \
          &Atom->TR[itrans].u2);

      Atom->TR[itrans].JKQ = (complex double **)MATRIX_RHO_CPLX(2, true);

      /*
      if(Mpi_pid == 0) fprintf(stderr, "%s %d %d %e %d %d %d %d\n", \
          Atom->Element, j, k, Atom->A[au][al], Atom->TR[itrans].M2, \
          Atom->TR[itrans].Output, Atom->TwoLv, true);
      */

      if(nread==5){
        Trim(tmp2, 3);
        String_to_Upper(tmp2);

        if(strcmp(tmp2, "YES") == 0 || strcmp(tmp2, "Y") == 0 ){
          outarry[Atom->Nline] = itrans;
          Atom->TR[itrans].geff = Geffect(Atom->LV[au].g, \
              Atom->LV[al].g,Atom->LV[au].J,Atom->LV[al].J);

          Atom->TR[itrans].const1 = Atom->LV[al].sqrt_deg \
              /Atom->LV[au].sqrt_deg*Atom->TR[itrans].Blu \
              /Atom->TR[itrans].Aul;
          Atom->TR[itrans].const2 = C_h*Atom->TR[itrans].nu/4./C_Pi \
              *Atom->LV[au].sqrt_deg*Atom->TR[itrans].Aul \
              *Input->dx*C_Solarradus;
          Atom->Nline++;
        }else{
          Atom->TR[itrans].geff = -100.;
          Atom->TR[itrans].delta = -100.;
          Atom->TR[itrans].w[0] = -100.;
          Atom->TR[itrans].w[1] = -100.;
          Atom->TR[itrans].w[2] = -100.;
        }

      }
    }

    Input->Nline += Atom->Nline;
    Input->Ntrans += Atom->Ntrans;

    Atom->iout = (int *)malloc(sizeof(int)*Atom->Nline);
    for(iline=0; iline<Atom->Nline; iline++){
      Atom->iout[iline] = outarry[iline];
    }

    Read_line(lines, Fa);
    Atom->Ion.NT = atof(lines);

    Atom->Ion.frac = (double **)MATRIX(0, 1, 0, Atom->Ion.NT-1, \
        enum_dbl, false);
    for(indx=0; indx<Atom->Ion.NT; indx++){
      Read_line(lines, Fa);
      nread = sscanf(lines,"%lf, %lf", &Atom->Ion.frac[0][indx], \
          &Atom->Ion.frac[1][indx]);
      Atom->Ion.frac[0][indx] = log10(Atom->Ion.frac[0][indx]);
    }

    if(Read_line(lines, Fa)>0){
      Atom->collision = true;
      Atom->col = (STRUCT_COL *)malloc(sizeof(STRUCT_COL));
      Atom->col->NT = String_elements(lines);
      Atom->col->T = (double *)VECTOR(0, Atom->col->NT-1, enum_dbl, false);
      for(indx=0; indx<Atom->col->NT; indx++){
        String_Split(parameter, lines);
        Atom->col->T[indx] = atof(parameter);
      }

      Atom->col->Strength = (double ***)TENSOR_DBL(0, Atom->Nlevel-1, \
          0, Atom->Nlevel-1, 0, Atom->col->NT-1, true);
      Atom->col->Rates = (double **)MATRIX(0, Atom->Nlevel-1, 0, \
          Atom->Nlevel-1, enum_dbl, false);
      
      for(al=0; al<Atom->Nlevel-1; al++){
        for(au=al+1; au<Atom->Nlevel; au++){
          Read_line(lines, Fa);

          String_Split(parameter, lines);
          i = atoi(parameter);
          String_Split(parameter, lines);
          j = atoi(parameter);

          if(al==i && au==j){
            for(indx=0; indx<Atom->col->NT; indx++){
              String_Split(parameter, lines);
              Atom->col->Strength[al][au][indx] = atof(parameter);
            }
          }else{
            
            Error(enum_error, "READ_ATOM", "index err of collision");
          }
        }
      }
    }else{
      Atom->collision = false;
    }
     
    fclose(Fa);
    return;
}

/*----------------------------------------------------------------------------*/

extern void FREE_ATOM(STRUCT_ATOM *Atom, int Natom){
    
    /*######################################################################
      Purpose:
        free the memory of atomic structure and collisional rates.
      Record of revisions:
        30 Nov. 2019.
      Input parameters:
        Atom, a structure saved the atomic information.
        Natom, number of atoms;
    ######################################################################*/

    int iatom, itrans, ilevel;

    for(iatom=0;iatom<Natom;iatom++){

      // free transitions
      for (itrans=0; itrans<Atom[iatom].Ntrans; itrans++) {
        FREE_MATRIX_RHO_CPLX(Atom->TR[itrans].JKQ);
        if(Atom[iatom].TwoLv) continue;
        FREE_TENSOR_DBL(Atom->TR[itrans].TA, 0, 0, 0);
        FREE_TENSOR_DBL(Atom->TR[itrans].TS, 0, 0, 0);
        FREE_TENSOR_DBL(Atom->TR[itrans].RA, 0, 0, 0);
        FREE_TENSOR_DBL(Atom->TR[itrans].RS, 0, 0, 0);
        FREE_VECTOR(Atom->TR[itrans].TE, 0, enum_dbl);
      }
      free(Atom[iatom].TR);
     
      // free levels
      for (ilevel=0; ilevel<Atom[iatom].Nlevel; ilevel++) {
        FREE_MATRIX_RHO_CPLX(Atom[iatom].LV[ilevel].Rho);
      }
      free(Atom[iatom].LV);
      
      // free ionization
      FREE_MATRIX(Atom[iatom].Ion.frac, 0, 0, enum_dbl);

      // free collision
      if(Atom[iatom].collision){
        FREE_VECTOR(Atom[iatom].col->T, 0, enum_dbl);
        FREE_TENSOR_DBL(Atom[iatom].col->Strength, 0, 0, 0);
        FREE_MATRIX(Atom[iatom].col->Rates, 0, 0, enum_dbl);


        free(Atom->col);
      }

      free(Atom[iatom].iout);
    }

    free(Atom->eqindx);
    free(Atom);
    
    return;
}

/*----------------------------------------------------------------------------*/

extern void Collisional_Rates(STRUCT_ATOM *Atom, int Natom, double T){
    
    /*######################################################################
      Purpose:
        interpolate a collisional rates matrix for a temperature T.
      Record of revisions:
        10 Aug. 2023.
      Input parameters:
        Atom, a structure saved the atomic information.
        Natom, number of atoms;
        T, the temperature.
      Output parameters:
        Atom[iatom].col->Rates[][], a matrix saved the collisional rates.
      Reference:
        Aggarwal 2014 MNRAS
    ######################################################################*/

    int iatom, iT, al, au, indx = 0;
    double DT1 = 0, DT2 = 0;
    double Strength_tmp, tmp;
    double T_log = log10(T);

    for(iatom=0;iatom<Natom;iatom++){
      if(Atom[iatom].collision){
        if((T_log < Atom[iatom].col->T[0])||(T_log > \
            Atom[iatom].col->T[Atom[iatom].col->NT-1])){
          for(au=0; au<Atom[iatom].Nlevel; au++){
            for(al=0; al<Atom[iatom].Nlevel; al++){
              Atom[iatom].col->Rates[au][al] = 0;
            }
          }
        }else{
          indx = -1;
          for(iT=0; iT<Atom[iatom].col->NT-2; iT++){
            if(T_log>=Atom[iatom].col->T[iT] && \
                T_log<Atom[iatom].col->T[iT+1]){
              indx = iT;
              DT1 = (T_log-Atom[iatom].col->T[iT]) \
                  /(Atom[iatom].col->T[iT+1]-Atom[iatom].col->T[iT]);
              DT2 = 1-DT1;
              break;
            }
          }
            
          for(al=0; al<Atom[iatom].Nlevel; al++){
            for(au=al+1; au<Atom[iatom].Nlevel; au++){
              Strength_tmp = Atom[iatom].col->Strength[al][au][indx]*DT2 \
                  +Atom[iatom].col->Strength[al][au][indx+1]*DT1;
              tmp = 8.63e-6*Strength_tmp/sqrt(T);
              Atom[iatom].col->Rates[au][al] = \
                  tmp/Atom[iatom].LV[au].deg;
              Atom[iatom].col->Rates[al][au] = \
                  tmp/Atom[iatom].LV[al].deg \
                  *exp((Atom[iatom].LV[al].Energy \
                  -Atom[iatom].LV[au].Energy)*C_h*C_c/T/C_Kb);
            }
          }
        }
      }
    }
    return;
}

/*----------------------------------------------------------------------------*/
