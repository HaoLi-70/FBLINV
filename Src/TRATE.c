
#include "TRATE.h"

/*----------------------------------------------------------------------------*/

    /*######################################################################
    
      revision log:

        30 Otc. 2024
          --- initial comment 
    
    ######################################################################*/

/*----------------------------------------------------------------------------*/

extern int TRates(STRUCT_ATOM *Atom, STR_FCTSG *fctsg, \
    STRUCT_INPUT *Input){

    /*######################################################################
      Purpose:
        Precompute the transition rates TA TE TS RA RE RS (JKQ is not 
            included).
      Record of revisions:
        17 Aug. 2023.
      Input parameters:
        Atom, a structure with the atoms.
        fctsg, a structure with factoral, signs.
        Input, a structure with the inputs.
      Output parameters:
        Atom, a structure with the atoms.
      Reference:
        Precompute coefficients in Eqs 7.14 without JKQ
        TA[KQ,KlQl,KrQr]
        TE = Aul               TE[KQ]    min(KuQu,KlQl)
        TS[KQ,KuQu,KrQr]
        RA[KQ,KpQp,KrQr]       ->KlQl
        RE/Aul
        RS[KQ,KpQp,KrQr]       ->KuQu
    ######################################################################*/

    int indx, itrans, au, al, K, Ku, Kl, Kr, Kp, Qu, Ql, Qr, Qp;
    int indxu, indxl, indxr, indx1, indx2;
    int sizeup, sizelow, Kmin;
    double Ju, Jl, tmp, tmp1, tmp2;

    const double c_tiny = 1e-15;

    for(itrans = 0; itrans<Atom->Ntrans; itrans++){

      au = Atom->TR[itrans].au;
      al = Atom->TR[itrans].al;
      Ju = Atom->LV[au].J;
      Jl = Atom->LV[al].J;

      sizeup = Atom->LV[au].nKQ-1;
      sizelow = Atom->LV[al].nKQ-1;

      Kmin = Atom->LV[au].Kmax>Atom->LV[al].Kmax?Atom->LV[al].Kmax: \
        Atom->LV[au].Kmax;

      Atom->TR[itrans].TA = (double ***)TENSOR_DBL(0, sizeup, \
          0, sizelow, 0, Input->nJKQ-1, true);
      Atom->TR[itrans].TS = (double ***)TENSOR_DBL(0, sizelow, \
          0, sizeup, 0, Input->nJKQ-1, true);

      Atom->TR[itrans].RA = (double ***)TENSOR_DBL(0, sizelow, \
          0, sizelow, 0, Input->nJKQ-1, true);
      Atom->TR[itrans].RS = (double ***)TENSOR_DBL(0, sizeup, \
          0, sizeup, 0, Input->nJKQ-1, true);

      if(Input->Saturated){
        Atom->TR[itrans].TE = (double *)VECTOR(0, Kmin/2,enum_dbl, true);

      }else{
        Atom->TR[itrans].TE = (double *)VECTOR(0, Kmin,enum_dbl, true);
      }

      if(Input->Saturated){
        for(Kr=0; Kr<=2; Kr+=Input->Kdelta){
          indxr = Kr/2;
          for(Ku=0; Ku<=Atom->LV[au].Kmax; Ku+=Input->Kdelta){
            indxu = Ku/2;
            for(Kl=0; Kl<=Atom->LV[al].Kmax; Kl+=Input->Kdelta){
              indxl = Kl/2;
              tmp = sqrt(3.*(2.*Kl+1.)*(2.*Ku+1.)*(2.*Kr+1.));
              Atom->TR[itrans].TA[indxu][indxl][indxr] = tmp \
                  *Atom->LV[al].deg*Atom->TR[itrans].Blu*fctsg->sg[Kl] \
                  *WIGNER_9J(Ju,Jl,1,Ju,Jl,1,Ku,Kl,Kr,fctsg) \
                  *WIGNER_3J(Ku,Kl,Kr,0.,0.,0.,fctsg);
              Atom->TR[itrans].TS[indxl][indxu][indxr] = tmp \
                  *Atom->LV[au].deg*Atom->TR[itrans].Bul*fctsg->sg[Kr+Ku] \
                  *WIGNER_9J(Jl,Ju,1,Jl,Ju,1,Kl,Ku,Kr,fctsg) \
                  *WIGNER_3J(Kl,Ku,Kr,0.,0.,0.,fctsg);
            }

            for(Kp=0; Kp<=Atom->LV[au].Kmax; Kp+=Input->Kdelta){
              indx1 = Kp/2;
              if ((Ku+Kp+Kr)%2==0){
                Atom->TR[itrans].RS[indxu][indx1][indxr] = \
                    sqrt(3.*(2.*Ku+1.)*(2.*Kp+1.)*(2.*Kr+1.)) \
                    *Atom->LV[au].deg*Atom->TR[itrans].Bul \
                    *WIGNER_6J(Ku,Kp,Kr,Ju,Ju,Ju,fctsg) \
                    *WIGNER_6J(1,1,Kr,Ju,Ju,Jl,fctsg) \
                    *fctsg->sg[lround(Jl-Ju)+1] \
                    *WIGNER_3J(Ku,Kp,Kr,0.,0.,0.,fctsg); 
              }
            }
          }
     
          for(Kl=0; Kl<=Atom->LV[al].Kmax; Kl+=Input->Kdelta){
            indx1 = Kl/2;
            for(Kp=0; Kp<=Atom->LV[al].Kmax; Kp+=Input->Kdelta){
              indx2 = Kp/2;
              if ((Kl+Kp+Kr)%2==0){
                Atom->TR[itrans].RA[indx1][indx2][indxr] = \
                    sqrt(3.*(2.*Kl+1.)*(2.*Kp+1.)*(2.*Kr+1.)) \
                    *Atom->LV[al].deg*Atom->TR[itrans].Blu \
                    *WIGNER_6J(Kl,Kp,Kr,Jl,Jl,Jl,fctsg) \
                    *WIGNER_6J(1,1,Kr,Jl,Jl,Ju,fctsg) \
                    *fctsg->sg[lround(Ju-Jl)+Kr+1] \
                    *WIGNER_3J(Kl,Kp,Kr,0.,0.,0.,fctsg) ;
              }
            }
          }
        }

        indx = 0;
        for(K=0; K<=Kmin; K+=Input->Kdelta){
          Atom->TR[itrans].TE[indx] = Atom->LV[au].deg \
              *Atom->TR[itrans].Aul*fctsg->sg[lround(Jl+Ju)+1+K] \
              *WIGNER_6J(Ju,Ju,K,Jl,Jl,1,fctsg);
          indx++;
        }

      }else{
        for(Kr=0; Kr<=2; Kr+=Input->Kdelta){
          for(Ku=0; Ku<=Atom->LV[au].Kmax; Ku+=Input->Kdelta){
            for(Kl=0; Kl<=Atom->LV[al].Kmax; Kl+=Input->Kdelta){
              tmp2 = sqrt(3.*(2.*Kl+1.)*(2.*Ku+1.)*(2.*Kr+1.));
              tmp1 = tmp2*WIGNER_9J(Ju,Jl,1,Ju,Jl,1,Ku,Kl,Kr,fctsg) \
                  *Atom->TR[itrans].Blu*Atom->LV[al].deg;
              tmp2 *= WIGNER_9J(Jl,Ju,1,Jl,Ju,1,Kl,Ku,Kr,fctsg) \
                  *Atom->TR[itrans].Bul*Atom->LV[au].deg;

              if(fabs(tmp1)<=c_tiny && fabs(tmp2)<=c_tiny) continue;

              for(Qu=-Ku; Qu<=Ku; Qu++){
                indxu = KQ_index(Ku,Qu,Input->Kdelta);
                for(Ql=-Kl; Ql<=Kl; Ql++){
                  indxl = KQ_index(Kl,Ql,Input->Kdelta);
                  for(Qr=-Kr; Qr<=Kr; Qr++){
                    indxr = KQ_index(Kr,Qr,Input->Kdelta);

                    Atom->TR[itrans].TA[indxu][indxl][indxr] = tmp1 \
                        *fctsg->sg[Kl+Ql] \
                        *WIGNER_3J(Ku,Kl,Kr,-Qu,Ql,-Qr,fctsg);
                    Atom->TR[itrans].TS[indxl][indxu][indxr] = tmp2 \
                        *fctsg->sg[Kr+Ku+Qu] \
                        *WIGNER_3J(Kl,Ku,Kr,-Ql,Qu,-Qr,fctsg);

                  }
                }
              }
            }

            for(Kp=0; Kp<=Atom->LV[au].Kmax; Kp+=Input->Kdelta){
              if ((Ku+Kp+Kr)%2==0){
                tmp1 = sqrt(3.*(2.*Ku+1.)*(2.*Kp+1.)*(2.*Kr+1.)) \
                    *Atom->LV[au].deg*Atom->TR[itrans].Bul \
                    *WIGNER_6J(Ku,Kp,Kr,Ju,Ju,Ju,fctsg) \
                    *WIGNER_6J(1,1,Kr,Ju,Ju,Jl,fctsg);

                if(fabs(tmp1)<=c_tiny) continue;

                for(Qu=-Ku; Qu<=Ku; Qu++){
                  indx1 = KQ_index(Ku,Qu,Input->Kdelta);
                  for(Qp=-Kp; Qp<=Kp; Qp++){
                    indx2 = KQ_index(Kp,Qp,Input->Kdelta);
                    for(Qr=-Kr; Qr<=Kr; Qr++){
                      indxr = KQ_index(Kr,Qr,Input->Kdelta);

                      Atom->TR[itrans].RS[indx1][indx2][indxr] = tmp1 \
                          *fctsg->sg[lround(Jl-Ju)+Qp+1] \
                          *WIGNER_3J(Ku,Kp,Kr,Qu,-Qp,Qr,fctsg);

                    }
                  }
                }
              }
            }
          }

          for(Kl=0; Kl<=Atom->LV[al].Kmax; Kl+=Input->Kdelta){
            for(Kp=0; Kp<=Atom->LV[al].Kmax; Kp+=Input->Kdelta){
              if ((Kl+Kp+Kr)%2==0){
                tmp1 = sqrt(3.*(2.*Kl+1.)*(2.*Kp+1.)*(2.*Kr+1.)) \
                    *Atom->LV[al].deg*Atom->TR[itrans].Blu \
                    *WIGNER_6J(Kl,Kp,Kr,Jl,Jl,Jl,fctsg) \
                    *WIGNER_6J(1,1,Kr,Jl,Jl,Ju,fctsg);

                if(fabs(tmp1)<=c_tiny) continue;

                for(Ql=-Kl; Ql<=Kl; Ql++){
                  indx1 = KQ_index(Kl,Ql,Input->Kdelta);
                  for(Qp=-Kp; Qp<=Kp; Qp++){
                    indx2 = KQ_index(Kp,Qp,Input->Kdelta);
                    for(Qr=-Kr; Qr<=Kr; Qr++){
                      indxr = KQ_index(Kr,Qr,Input->Kdelta);

                      Atom->TR[itrans].RA[indx1][indx2][indxr] = tmp1 \
                          *fctsg->sg[lround(Ju-Jl)+Kr+Qp+1] \
                          *WIGNER_3J(Kl,Kp,Kr,Ql,-Qp,Qr,fctsg);

                    }
                  }
                }
              }
            }
          }
        }
          
        indx = 0;
        for(K=0; K<=Kmin; K+=Input->Kdelta){
          Atom->TR[itrans].TE[indx] = Atom->LV[au].deg \
              *Atom->TR[itrans].Aul*fctsg->sg[lround(Jl+Ju)+1+K] \
              *WIGNER_6J(Ju,Ju,K,Jl,Jl,1,fctsg);
          indx++;
        }
      }
    }
    
    return 0;
}

/*----------------------------------------------------------------------------*/