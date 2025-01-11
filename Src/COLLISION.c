
#include "COLLISION.h"

/*----------------------------------------------------------------------------*/

    /*######################################################################
    
      revision log:

        30 Otc. 2024
          --- initial comment 
    
    ######################################################################*/

/*----------------------------------------------------------------------------*/

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
