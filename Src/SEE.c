
#include "SEE.h"

/*----------------------------------------------------------------------------*/

    /*######################################################################
     
      revision log:

        19 Sep. 2024
          --- update: Parameter structure.

        30 Otc. 2024
          --- update: Simplification in the calculation of some 
                      coefficients.
                      moved subroutine Transition_Rates to TRATE.c
          --- bugfix: fix an invalid address reading.
     
    ######################################################################*/

/*----------------------------------------------------------------------------*/

extern complex double **SEE(STRUCT_ATOM *Atom, STR_FCTSG *fctsg, \
    STRUCT_PARA *Para, STRUCT_INPUT *Input){

    /*######################################################################
      Purpose:
        Build the statistic equilibrium equation.
      Record of revisions:
        19 Sep. 2024.
      Input parameters:
        Atom, a structure with the atoms.
        fctsg, a structure with factoral, signs.
        Para, a structure with physical parameter.
        Input, a structure with the inputs.
      Return:
        the matrix.
    ######################################################################*/
    
    const double c_tiny = 1e-15;

    complex double **matrix = (complex double **)MATRIX(0,Atom->nEq, \
        0, Atom->nEq, enum_cplx, true);

    int au, al, K, Kr, Ku, Kl, Kp, Q, Qr, Qu, Ql, Qp, Kmin;
    int indx, itrans, indxl, indxu, indxp, indxr, iu, il, ip;
    double tmp, tmp1, tmp2, Ju, Jl;

    if(Input->Saturated){

      for(itrans = 0; itrans<Atom->Ntrans; itrans++){
        au = Atom->TR[itrans].au;
        al = Atom->TR[itrans].al;
        Ju = Atom->LV[au].J;
        Jl = Atom->LV[al].J;

        for(Kr=0;Kr<=2;Kr+=Input->Kdelta){
          indxr = Kr/2;
          for(Ku=0; Ku<=Atom->LV[au].Kmax; Ku+=Input->Kdelta){
            indxu = Ku/2;
            iu = indxu+Atom->eqindx[au];
            for(Kl=0; Kl<=Atom->LV[al].Kmax; Kl+=Input->Kdelta){
              indxl = Kl/2;
              il = indxl+Atom->eqindx[al];
              matrix[iu][il] += Atom->TR[itrans].JKQ[Kr][0] \
                  *Atom->TR[itrans].TA[indxu][indxl][indxr];
              matrix[il][iu] += Atom->TR[itrans].JKQ[Kr][0] \
                  *Atom->TR[itrans].TS[indxl][indxu][indxr];
                  
            }

            for(Kp=0; Kp<=Atom->LV[au].Kmax; Kp+=Input->Kdelta){
              indxp = Kp/2;
              ip = indxp+Atom->eqindx[au];
              matrix[iu][ip] -= Atom->TR[itrans].JKQ[Kr][0] \
                  *Atom->TR[itrans].RS[indxu][indxp][indxr];
                  
            }                
          }

          for(Kl=0; Kl<=Atom->LV[al].Kmax; Kl+=Input->Kdelta){
            indxl = Kl/2;
            il = indxl+Atom->eqindx[al];
            for(Kp=0; Kp<=Atom->LV[al].Kmax; Kp+=Input->Kdelta){
              indxp = Kp/2;
              ip = indxp+Atom->eqindx[al];

              matrix[il][ip] -= Atom->TR[itrans].JKQ[Kr][0] \
                  *Atom->TR[itrans].RA[indxl][indxp][indxr];
                     
            }
          }
        }

        Kmin = Atom->LV[au].Kmax<Atom->LV[al].Kmax? \
            Atom->LV[au].Kmax:Atom->LV[al].Kmax;

        for(K=0;K<=Kmin;K+=Input->Kdelta){
          indx = K/2;
          iu = indx+Atom->eqindx[au];
          il = indx+Atom->eqindx[al];
          matrix[il][iu] +=  Atom->TR[itrans].TE[indx];
          
        }

        for(K=0;K<=Atom->LV[au].Kmax;K+=Input->Kdelta){  
          indx = K/2;
          iu = indx+Atom->eqindx[au];
          //RE 
          matrix[iu][iu] -= Atom->TR[itrans].Aul;
        }    
       
      }

      if(!Atom->collision) return matrix;

      for(al=0; al<Atom->Nlevel-1; al++){
        for(au=al+1; au<Atom->Nlevel; au++){
          if(Atom->col->Rates[au][al]<=0) continue;
            
          tmp = Atom->LV[au].sqrt_deg/Atom->LV[al].sqrt_deg;
          Kmin = Atom->LV[au].Kmax<Atom->LV[al].Kmax? \
              Atom->LV[au].Kmax:Atom->LV[al].Kmax;
        
          Ju = Atom->LV[au].J;
          Jl = Atom->LV[al].J;

          for(K=0; K<=Kmin; K+=Input->Kdelta){

            if(K==0){
              tmp1 = 1.0;
              tmp2 = 1.0;
            }else{
              //inelastic collision Clu lower to upper level
              tmp1 = WIGNER_6J(Ju,Ju,0.0,Jl,Jl,1.0,fctsg);
              //super elastic collision Cul upper to lower level
              tmp2 = WIGNER_6J(Jl,Jl,0.0,Ju,Ju,1.0,fctsg);
              if(fabs(tmp1)<c_tiny){
                tmp1 = 0.0;
              }else{
                tmp1 = fctsg->sg[K]*WIGNER_6J(Ju,Ju,K, \
                    Jl,Jl,1.0,fctsg)/tmp1;               
              }
              if(fabs(tmp2)<c_tiny){
                tmp2 = 0.0;
              }else{
                tmp2 = fctsg->sg[K]*WIGNER_6J(Jl,Jl,K, \
                    Ju,Ju,1.0,fctsg)/tmp2;             
              }
            }
              
            indx = K/2;
            iu = indx+Atom->eqindx[au];
            il = indx+Atom->eqindx[al];
            matrix[iu][il] += tmp1*Atom->col->Rates[al][au]*Para->ne/tmp;
            matrix[il][iu] += tmp*tmp2*Atom->col->Rates[au][al]*Para->ne;

            matrix[iu][iu] -= Atom->col->Rates[au][al]*Para->ne;
            matrix[il][il] -= Atom->col->Rates[al][au]*Para->ne;

          }    
        }
      }

    }else{

      //if Hanle regime
      for(al=0; al<Atom->Nlevel; al++){
        if(Atom->LV[al].Kmax<2) continue;
        for(K=0;K<=Atom->LV[al].Kmax;K+=Input->Kdelta){
          for(Q=-K; Q<=K; Q++){
            indxu = KQ_index(K,Q,Input->Kdelta);
            iu = indxu+Atom->eqindx[al];
            matrix[iu][iu] -= 0.879e7*Atom->LV[al].g*Para->B*Q*I;
          }
        }
      }

      for(itrans = 0; itrans<Atom->Ntrans; itrans++){
        au = Atom->TR[itrans].au;
        al = Atom->TR[itrans].al;
        Ju = Atom->LV[au].J;
        Jl = Atom->LV[al].J;

        for(Kr=0;Kr<=2;Kr+=Input->Kdelta){
          for(Qr=-Kr;Qr<=Kr;Qr++){
            indxr = KQ_index(Kr,Qr,Input->Kdelta);
            for(Ku=0; Ku<=Atom->LV[au].Kmax; Ku+=Input->Kdelta){
              for(Qu=-Ku; Qu<=Ku; Qu++){
                indxu = KQ_index(Ku,Qu,Input->Kdelta);
                iu = indxu+Atom->eqindx[au];
                for(Kl=0; Kl<=Atom->LV[al].Kmax; Kl+=Input->Kdelta){
                  for(Ql=-Kl; Ql<=Kl; Ql++){
                    indxl = KQ_index(Kl,Ql,Input->Kdelta);
                    il = indxl+Atom->eqindx[al];
                    matrix[iu][il] += Atom->TR[itrans].JKQ[Kr][Qr] \
                        *Atom->TR[itrans].TA[indxu][indxl][indxr];
                    matrix[il][iu] += Atom->TR[itrans].JKQ[Kr][Qr] \
                        *Atom->TR[itrans].TS[indxl][indxu][indxr];
                  }
                }

                for(Kp=0; Kp<=Atom->LV[au].Kmax; Kp+=Input->Kdelta){
                  for(Qp=-Kp; Qp<=Kp; Qp++){
                    indxp = KQ_index(Kp,Qp,Input->Kdelta);
                    ip = indxp+Atom->eqindx[au];
                    matrix[iu][ip] -= Atom->TR[itrans].JKQ[Kr][Qr] \
                        *Atom->TR[itrans].RS[indxu][indxp][indxr];
                  }
                }
              }            
            }

            for(Kl=0; Kl<=Atom->LV[al].Kmax; Kl+=Input->Kdelta){
              for(Ql=-Kl; Ql<=Kl; Ql++){
                indxl = KQ_index(Kl,Ql,Input->Kdelta);
                il = indxl+Atom->eqindx[al];
                for(Kp=0; Kp<=Atom->LV[al].Kmax; Kp+=Input->Kdelta){
                  for(Qp=-Kp; Qp<=Kp; Qp++){
                    indxp = KQ_index(Kp,Qp,Input->Kdelta);
                    ip = indxp+Atom->eqindx[al];
                    matrix[il][ip] -= Atom->TR[itrans].JKQ[Kr][Qr] \
                        *Atom->TR[itrans].RA[indxl][indxp][indxr];
                    
                  }
                }
              }
            }
          }
        }


        Kmin = Atom->LV[au].Kmax<Atom->LV[al].Kmax? \
            Atom->LV[au].Kmax:Atom->LV[al].Kmax;

        for(K=0;K<=Kmin;K+=Input->Kdelta){
          for(Q=-K;Q<=K;Q++){
            indx = KQ_index(K,Q,Input->Kdelta);
            iu = indx+Atom->eqindx[au];
            il = indx+Atom->eqindx[al];
            matrix[il][iu] +=  Atom->TR[itrans].TE[K];
          }
        }

        for(K=0;K<=Atom->LV[au].Kmax;K+=Input->Kdelta){
          for(Q=-K;Q<=K;Q++){
            indx = KQ_index(K,Q,Input->Kdelta);
            iu = indx+Atom->eqindx[au];
            //RE 
            matrix[iu][iu] -= Atom->TR[itrans].Aul;
          }            
        }
      }

      if(!Atom->collision) return matrix;

      // collision with electron
      for(al=0; al<Atom->Nlevel-1; al++){
        for(au=al+1; au<Atom->Nlevel; au++){
          if(Atom->col->Rates[au][al]<=0) continue;
            
          tmp = Atom->LV[au].sqrt_deg/Atom->LV[al].sqrt_deg;
          Kmin = Atom->LV[au].Kmax<Atom->LV[al].Kmax? \
              Atom->LV[au].Kmax:Atom->LV[al].Kmax;

          Ju = Atom->LV[au].J;
          Jl = Atom->LV[al].J;
          for(K=0; K<=Kmin; K+=Input->Kdelta){

            if(K==0){
              tmp1 = 1.0;
              tmp2 = 1.0;
            }else{
              //inelastic collision Clu lower to upper level
              tmp1 = WIGNER_6J(Ju,Ju,0.0,Jl,Jl,1.0,fctsg);
              //super elastic collision Cul upper to lower level
              tmp2 = WIGNER_6J(Jl,Jl,0.0,Ju,Ju,1.0,fctsg);
              if(fabs(tmp1)<c_tiny){
                tmp1 = 0.0;
              }else{
                tmp1 = fctsg->sg[K]*WIGNER_6J(Ju,Ju,K, \
                    Jl,Jl,1.0,fctsg)/tmp1;               
              }
              if(fabs(tmp2)<c_tiny){
                tmp2 = 0.0;
              }else{
                tmp2 = fctsg->sg[K]*WIGNER_6J(Jl,Jl,K, \
                    Ju,Ju,1.0,fctsg)/tmp2;             
              }
            }
              
            for(Q=-K; Q<=K; Q++){
              indx = KQ_index(K,Q,Input->Kdelta);

              iu = indx+Atom->eqindx[au];
              il = indx+Atom->eqindx[al];
    
              matrix[iu][il] += tmp1*Atom->col->Rates[al][au]*Para->ne/tmp;
              matrix[il][iu] += tmp*tmp2*Atom->col->Rates[au][al]*Para->ne;

              matrix[iu][iu] -= Atom->col->Rates[au][al]*Para->ne;
              matrix[il][il] -= Atom->col->Rates[al][au]*Para->ne;
            }
          }    
        }
      }      
    }

    return matrix;
}

/*----------------------------------------------------------------------------*/


