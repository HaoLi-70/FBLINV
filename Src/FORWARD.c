
#include "FORWARD.h"

/*----------------------------------------------------------------------------*/

    /*######################################################################
    
      revision log:

        30 Otc. 2024
          --- update: moved subroutines Coeff2NLM and NLM2Coeff to 
                      COEFF.c
                      moved subroutines FB_RECONSTRUTION and 
                      SFB_RECONSTR_DELTA to RECONSTRUCT.c

        28 Otc. 2024
          --- bugfix: rho2q was used even for the transition J=1/2.
          
        19 Sep. 2024
          --- update: synthesis for multi LOS.

    ######################################################################*/

/*----------------------------------------------------------------------------*/

extern int FORWARD(STRUCT_SYN *Syn, STRUCT_ATOM *Atom, \
    STRUCT_INPUT *Input, STR_FCTSG *fctsg, STRUCT_MPI *Mpi, bool sparse){
  
    /*######################################################################
      Purpose:
        forward synthesis of the polarization.
      Record of revisions:
        28 Otc. 2024
      Input parameters:
        Syn, a structure with forward synthesis.
        Atom, a structure with the atomic information.
        Input, a structure with the input information.
        fctsg, a structure with factoral, signs.
        Mpi, a structure with the Mpi information.
        sparse, a flag used in the inversion
      Output parameters:
        Syn, a structure with forward synthesis.
    ######################################################################*/
  
    STRUCT_GRID *pgrid;
    STRUCT_LOS *plos;
    double cont2, epsilon, Bp, Hu, norm;
    int Si, ieq, ikq, Q, Qp, K, indx, igrid = 0, au, ilos;
    int itrans, iline, ilines, iatom, ilevel, iThom, shift;

    complex double **matrix_tmp,  *b;
    int *LUindx;

    for(igrid=0; igrid<Mpi->ngrids; igrid++){

      pgrid = Syn->Grids+igrid;

      for(ilos=0;ilos<Input->nlos;ilos++){
        plos = pgrid->los+ilos;

        if(sparse){
          if(Ran1(Mpi->idum) < Input->Per){ 
            
            if(Input->NThom > 0){
              for(iThom=0;iThom<Input->NThom;iThom++){
                plos->Thom[iThom].val[0] = \
                    plos->Thom_unperturb[iThom].val[0];
                plos->Thom[iThom].val[1] = \
                    plos->Thom_unperturb[iThom].val[1];
              }
            }
            if(Input->Nline > 0){
              for(iline=0;iline<Input->Nline;iline++){
                plos->Line[iline].Stk[0] = \
                    plos->Line_unperturb[iline].Stk[0];
                plos->Line[iline].Stk[1] = \
                    plos->Line_unperturb[iline].Stk[1];
                plos->Line[iline].Stk[2] = \
                    plos->Line_unperturb[iline].Stk[2];
              }
            }
            continue;
          }
        }

        Collisional_Rates(Atom, Input->Natom, plos->Para.T);

        Rotmat(plos->Para.PhiB, plos->Para.ThetaB, 0, Input->Dkmn, 2);

        Input->Tkq[0][0][0] = 1;
        for(K=0; K<3; K++){
          for(Qp=-Input->TQmax; Qp<=Input->TQmax; Qp++){
            Input->Tkq[K][2][Qp] = 0;
            for(Q=-2; Q<=2; Q++){
              Input->Tkq[K][2][Qp] += pgrid->T2Q[K][Q] \
                  *Input->Dkmn[2][Q][Qp];
            }
          }
          if(Input->Verbose >= 2  && (Mpi->rank) == 0){
            fprintf(stderr,"%d Input->Tkq new == %e %e %e %e %e \n", K, \
                creal(Input->Tkq[K][2][-2]), creal(Input->Tkq[K][2][-1]), \
                creal(Input->Tkq[K][2][0]), creal(Input->Tkq[K][2][1]), \
                creal(Input->Tkq[K][2][2]));
          }
        }


        shift = 0;
        for(iatom=0; iatom<Input->Natom; iatom++){

          if(Input->Verbose >= 2  && (Mpi->rank) == 0){
            fprintf(stderr,"\n ion = %e, ne = %e, T = %e\n", \
                plos->Para.Ion[iatom], plos->Para.ne,plos->Para.T);

            fprintf(stderr,"collision rates = %e \n", \
                Atom[iatom].col->Rates[1][0]);

            fprintf(stderr,"bolzman %e", \
                exp(-C_h*Atom[iatom].TR[0].nu/C_Kb/plos->Para.T) \
                *Atom[iatom].LV[1].deg/Atom[iatom].LV[0].deg);
            fprintf(stderr,"nline = %d %d\n", Atom[iatom].Nline, \
                Input->Natom);
          }

          if(Atom[iatom].Nline<=0) continue;

          if(Input->Verbose >= 2  && (Mpi->rank) == 0){

            fprintf(stderr,"jkq vertical == %e %e  \n", \
                pgrid->Jkq[iatom].J00[0], pgrid->Jkq[iatom].J2q[0][0]);

          }
          // rotate Jkq to the coordinate of magnetic field.
          // LL04 Chapter 7, page 331, Eq. 7.81
          // J^K_Qnew = sum_Qp(D^K_QpQ*J^K_Qp)
          // if only J^2_0 (local vertical) is non-zero
          for(itrans=0; itrans<Atom[iatom].Ntrans; itrans++){
            Atom[iatom].TR[itrans].JKQ[0][0] = \
                pgrid->Jkq[iatom].J00[itrans];
            for(Qp=-Input->TQmax; Qp<=Input->TQmax; Qp++){
              Atom[iatom].TR[itrans].JKQ[2][Qp] = 0;
              for(Q=-Input->JQmax; Q<=Input->JQmax; Q++){
                Atom[iatom].TR[itrans].JKQ[2][Qp] += \
                    pgrid->Jkq[iatom].J2q[itrans][Q]*Input->Dkmn[2][Q][Qp];
              }
            }
          }

          if(Input->Verbose >= 2  && (Mpi->rank) == 0){

            fprintf(stderr,"j00 new == %e  \n", \
                creal(Atom[iatom].TR[0].JKQ[0][0]));
            fprintf(stderr,"j2 new == %e %e %e %e %e \n", \
                creal(Atom[iatom].TR[0].JKQ[2][-2]),\
                creal(Atom[iatom].TR[0].JKQ[2][-1]),\
                creal(Atom[iatom].TR[0].JKQ[2][-0]),\
                creal(Atom[iatom].TR[0].JKQ[2][1]),\
                creal(Atom[iatom].TR[0].JKQ[2][2]));
          }

          //0.5*(3.0*Cr*Cr-1);
          // Chapter 10, Page 532, Eq. 10.50
          if(Atom[iatom].TwoLv){
            epsilon = Atom[iatom].col->Rates[1][0] \
                *plos->Para.ne/Atom[iatom].TR[0].Aul;

            // notic Bp in Eq. 10.50 is not a Planck function, it is
            // actually the coefficient due to the conversions between
            // Blu to Aul and CI to CS, finally turned to the 
            // Wien limit
            Bp = Wien_limit(Atom[iatom].TR[0].cplank1, \
                Atom[iatom].TR[0].cplank2, plos->Para.T);

            Atom[iatom].LV[1].Rho[0][0] = Atom[iatom].TR[0].const1 \
                *(epsilon*Bp+Atom[iatom].TR[0].JKQ[0][0])/(1+epsilon);

            if(Atom[iatom].LV[1].Kmax>=2){
              if(Input->Saturated){

                Atom[iatom].LV[1].Rho[2][0] = Atom[iatom].TR[0].const1 \
                    *Atom[iatom].TR[0].w[2]*Atom[iatom].TR[0].JKQ[2][0] \
                    /(1+epsilon);  

              }else{

                Hu = 0.879e7*Atom[iatom].LV[1].g*plos->Para.B \
                    /Atom[iatom].TR[0].Aul;

                for(Q=-2; Q<=2; Q++){
                  Atom[iatom].LV[1].Rho[2][Q] = Atom[iatom].TR[0].const1 \
                      *(Atom[iatom].TR[0].w[2] \
                      *Atom[iatom].TR[0].JKQ[2][-Q] \
                      *fctsg->sg[Q])/(1+epsilon+Q*Hu*I);
                }  
              }
            }
      

          }else{

            matrix_tmp = SEE(Atom+iatom, fctsg, &(plos->Para), Input);

            b = (complex double *)VECTOR(1, Atom[iatom].nEq, enum_cplx, \
                false);
            LUindx = (int *)VECTOR(1, Atom[iatom].nEq, enum_int, false);

            for (ieq=1; ieq<=Atom[iatom].nEq; ieq++) {
              b[ieq] = -matrix_tmp[ieq][0];
            }

            //do the LU decomposition
            ludcmp_cplx(matrix_tmp, Atom[iatom].nEq, LUindx);
            
            //sove the equations output the solution in Rho_M
            lubksb_cplx(matrix_tmp, Atom[iatom].nEq, LUindx, b);

            norm = Atom[iatom].LV[0].sqrt_deg;

            for(ilevel=1;ilevel<Atom[iatom].Nlevel;ilevel++){
              norm += Atom[iatom].LV[ilevel].sqrt_deg*creal(b[ilevel]);
            }  

            
            for(K=2;K<=Atom[iatom].LV[0].Kmax;K+=Input->Kdelta){
              if(Input->Saturated){
                indx = Atom[iatom].eqindx[0]+K/2;
                Atom[iatom].LV[0].Rho[K][0] = b[indx];
              }else{
                for(Q=-K;Q<=K;Q++){
                  ikq = KQ_index(K,Q,Input->Kdelta);
                  indx = Atom[iatom].eqindx[0]+ikq;
                  Atom[iatom].LV[0].Rho[K][Q] = b[indx];
                }
              }
            }

            for(ilevel=1;ilevel<Atom[iatom].Nlevel;ilevel++){
              for(K=0;K<=Atom[iatom].LV[ilevel].Kmax;K+=Input->Kdelta){
                if(Input->Saturated){
                  indx = Atom[iatom].eqindx[ilevel]+K/2;
                  Atom[iatom].LV[ilevel].Rho[K][0] = b[indx];

                }else{
                  for(Q=-K;Q<=K;Q++){
                    ikq = KQ_index(K,Q,Input->Kdelta);
                    indx = Atom[iatom].eqindx[ilevel]+ikq;
                    Atom[iatom].LV[ilevel].Rho[K][Q] = b[indx];

                  }
                }
              }
            }
                      
            FREE_MATRIX(matrix_tmp, 0, 0, enum_cplx);
            FREE_VECTOR(b, 1, enum_cplx);
            FREE_VECTOR(LUindx, 1, enum_int);
          }

          if(Input->Verbose >= 2  && (Mpi->rank) == 0){
            fprintf(stderr,"rho00 %e \n", \
                creal(Atom[iatom].LV[1].Rho[0][0]));

            if(Atom[iatom].LV[1].Kmax>=2){
              fprintf(stderr,"rho20 %e \n", \
                  creal(Atom[iatom].LV[1].Rho[2][0]));
            }
          }
          
          norm = Atom[iatom].LV[0].sqrt_deg;

          for(ilevel=1;ilevel<Atom[iatom].Nlevel;ilevel++){
            norm += Atom[iatom].LV[ilevel].sqrt_deg \
                *Atom[iatom].LV[ilevel].Rho[0][0];
          } 
          

          for(iline=0;iline<Atom[iatom].Nline;iline++){
            ilines = iline+shift;
            plos->Line[ilines].Stk[0] = 0;
            plos->Line[ilines].Stk[1] = 0;
            plos->Line[ilines].Stk[2] = 0;
            itrans = Atom[iatom].iout[iline];
            au = Atom[iatom].TR[Atom[iatom].iout[iline]].au;
            if(Input->Saturated){
              for(K=0;K<=Atom[iatom].LV[au].Kmax;K+=Input->Kdelta){
                for(Si=0;Si<3;Si++){
                  plos->Line[ilines].Stk[Si] += \
                      Atom[iatom].TR[itrans].const2 \
                      *Atom[iatom].TR[iline].w[K] \
                      *creal(Input->Tkq[Si][K][0] \
                      *Atom[iatom].LV[au].Rho[K][0]);
                }
              }

            }else{
              for(K=0;K<=Atom[iatom].LV[au].Kmax;K+=Input->Kdelta){
                for(Q=-K;Q<=K;Q++){
                  for(Si=0;Si<3;Si++){
                    plos->Line[ilines].Stk[Si] += \
                        Atom[iatom].TR[itrans].const2 \
                        *Atom[iatom].TR[iline].w[K] \
                        *creal(Input->Tkq[Si][K][Q] \
                        *Atom[iatom].LV[au].Rho[K][Q]);

                  }
                }
              }
            }
            if(Input->OutputV){
              if(Atom[iatom].LV[au].Kmax>=2){
                plos->Line[ilines].Stk[3] = Atom[iatom].TR[itrans].const2 \
                    *plos->Blos*C_Nul*Atom[iatom].TR[itrans].lambda \
                    *Atom[iatom].TR[itrans].lambda/C_c \
                    *(Atom[iatom].TR[itrans].geff \
                    *Atom[iatom].LV[au].Rho[0][0] \
                    +Atom[iatom].TR[itrans].delta \
                    *Atom[iatom].LV[au].Rho[2][0]) \
                    *plos->Para.Ion[iatom]*plos->Para.nH/norm;
              }else{
                plos->Line[ilines].Stk[3] = Atom[iatom].TR[itrans].const2 \
                    *plos->Blos*C_Nul*Atom[iatom].TR[itrans].lambda \
                    *Atom[iatom].TR[itrans].lambda/C_c \
                    *(Atom[iatom].TR[itrans].geff \
                    *Atom[iatom].LV[au].Rho[0][0]) \
                    *plos->Para.Ion[iatom]*plos->Para.nH/norm;                
              }
              if(Input->Verbose >= 2  && (Mpi->rank) == 0){
                fprintf(stderr, "V %e %e %e %e %e \n", \
                    plos->Line[ilines].Stk[3], \
                    Atom[iatom].TR[itrans].const2, \
                    plos->Blos*C_Nul*Atom[iatom].TR[itrans].lambda \
                    *Atom[iatom].TR[itrans].lambda/C_c, 
                    creal(Atom[iatom].TR[itrans].geff \
                    *Atom[iatom].LV[au].Rho[0][0] \
                    +Atom[iatom].TR[itrans].delta \
                    *Atom[iatom].LV[au].Rho[2][0]), \
                    plos->Para.Ion[iatom]*plos->Para.nH/norm);
                fprintf(stderr, "geff %e delta %e au =%d\n", \
                    Atom[iatom].TR[itrans].geff, \
                    Atom[iatom].TR[itrans].delta, au);
                fprintf(stderr, "natom %e  %e norm =%e\n", \
                    plos->Para.Ion[iatom], \
                    plos->Para.nH, norm);

              }

            }

            plos->Line[ilines].Stk[0] *= plos->Para.Ion[iatom] \
                *plos->Para.nH/norm;
            if(Atom[iatom].TR[Atom[iatom].iout[iline]].M1){
              plos->Line[ilines].Stk[1] *= -plos->Para.Ion[iatom] \
                  *plos->Para.nH/norm;
              plos->Line[ilines].Stk[2] *= -plos->Para.Ion[iatom] \
                  *plos->Para.nH/norm;
            }else{
              plos->Line[ilines].Stk[1] *= plos->Para.Ion[iatom] \
                  *plos->Para.nH/norm;
              plos->Line[ilines].Stk[2] *= plos->Para.Ion[iatom] \
                  *plos->Para.nH/norm;
            }
            
            if(Input->Verbose >= 2  && (Mpi->rank) == 0){
              fprintf(stderr, "iline %d %d %d \n",iline,ilines,iatom);
              fprintf(stderr, "ion %e  \n", plos->Para.Ion[iatom] \
                  *plos->Para.nH/norm);
              fprintf(stderr, "IQU %e %e %e \n", \
                  plos->Line[ilines].Stk[0], \
                  plos->Line[ilines].Stk[1], plos->Line[ilines].Stk[2]);
              if(Input->OutputV){
                fprintf(stderr, "V %e \n",plos->Line[ilines].Stk[3]);
              }

            }
        
          }

          shift += Atom[iatom].Nline;
        }  


        for(iThom=0; iThom<Input->NThom; iThom++){
          cont2 = plos->Para.ne*Input->dx*Input->Thom[iThom].Const;
          plos->Thom[iThom].val[0] = (pgrid->Kr[iThom]+pgrid->Kt[iThom]) \
              *cont2;
          plos->Thom[iThom].val[1] = (pgrid->Kt[iThom]-pgrid->Kr[iThom]) \
              *cont2;
        }

      }

    }

    return 0;
}

/*----------------------------------------------------------------------------*/

extern int Copy_Par(STRUCT_SYN *Syn, STRUCT_INPUT *Input, STRUCT_MPI *Mpi){

    /*######################################################################
      Purpose:
        Save the unperturbed parameters.
      Record of revisions:
        19 Sep. 2024
      Input parameters:
        Syn, a structure with forward synthesis.
        Input, a structure with the input information.
        Mpi, a structure with the Mpi information.
      Output parameters:
        Syn, a structure with forward synthesis.
    ######################################################################*/

    STRUCT_GRID *pgrid;
    STRUCT_LOS *plos;

    int igrid, iatom, iThom, iline, ilos;

    for(igrid=0; igrid<Mpi->ngrids; igrid++){
      pgrid = Syn->Grids+igrid;

      for(ilos=0; ilos<Input->nlos; ilos++){
        plos = pgrid->los+ilos;
        plos->Para_unperturb->T = plos->Para.T;
        plos->Para_unperturb->nH = plos->Para.nH;
        plos->Para_unperturb->ne = plos->Para.ne;
        plos->Para_unperturb->B = plos->Para.B;
        plos->Para_unperturb->ThetaB = plos->Para.ThetaB;
        plos->Para_unperturb->PhiB = plos->Para.PhiB;
        plos->Para_unperturb->Br = plos->Para.Br;
        plos->Para_unperturb->Bt = plos->Para.Bt;
        plos->Para_unperturb->Bp = plos->Para.Bp;

        for(iatom=0; iatom<Input->Natom; iatom++){
          plos->Para_unperturb->Ion[iatom] = plos->Para.Ion[iatom];
  
        }

        if(Input->NThom > 0){
          for(iThom=0;iThom<Input->NThom;iThom++){
            plos->Thom_unperturb[iThom].val[0] = plos->Thom[iThom].val[0];
            plos->Thom_unperturb[iThom].val[1] = plos->Thom[iThom].val[1];
          }
        }
        if(Input->Nline > 0){
          for(iline=0;iline<Input->Nline;iline++){
            plos->Line_unperturb[iline].Stk[0] = plos->Line[iline].Stk[0];
            plos->Line_unperturb[iline].Stk[1] = plos->Line[iline].Stk[1];
            plos->Line_unperturb[iline].Stk[2] = plos->Line[iline].Stk[2];
          }
        }

      }

      
    }

    return 0;
}

/*----------------------------------------------------------------------------*/

extern int Grid2Pixel(STRUCT_SYN *Syn, STRUCT_ATOM *Atom, \
    STRUCT_INPUT *Input, STRUCT_OUT *Output, STRUCT_MPI *Mpi){
  
    /*######################################################################
      Purpose:
        convert the synthetic data to the observation style.
      Record of revisions:
        19 Sep. 2024
      Input parameters:
        Syn, a structure with the grids.
        Atom, a structure with the atoms.
        Input, a structure with the inputs.
        Mpi, a structure with the Mpi information.
        perturb, compute the synthesis or a difference
      Output parameters:
        Output, a structure with the output matrix
     ######################################################################*/
  
    int igrid, istk, iline, ilines, iThom, iatom, itrans, ispec, \
        ilambda, ilos, indx, shift;
    double clambda = 25.0, losshift, lcenter;
    double sqrtT, Dopp, prof, *dptmp, dlambda;

    STRUCT_GRID *pgrid;
    STRUCT_LOS *plos;

    int itmp1 = Input->Nstk*Input->Nline;
    int itmp2 = Syn->npspec*Input->Nl*Input->Nstk;
    int itmp3 = Syn->npixels*(itmp1+2*Input->NThom);


    for(ilos=0; ilos<Input->nlos; ilos++){
      // initialize to 0;
      memset(Output->synloc[0], 0, itmp3*sizeof(double));

      if(Input->Nspec>0){      
        memset(Output->specloc[0], 0, itmp2*sizeof(double));
      }


      for(igrid=0; igrid<Mpi->ngrids; igrid++){
        pgrid = Syn->Grids+igrid;

        plos = pgrid->los+ilos;

  //fprintf(stderr,"grid %d \n",igrid);
        for(iline=0; iline<Input->Nline; iline++){
          for(istk=0; istk<Input->Nstk; istk++){
            Output->synloc[pgrid->ipixel][istk+Input->Nstk*iline] += \
                plos->Line[iline].Stk[istk];            
          }
        }

        for(iThom=0; iThom<Input->NThom; iThom++){
          dptmp = Output->synloc[pgrid->ipixel]+2*iThom+itmp1;
          *dptmp += plos->Thom[iThom].val[0];
          *(dptmp+1) += plos->Thom[iThom].val[1];
        }

        // compute the spectrum
        if(Input->Nspec>0){
          if(pgrid->ipspec>=0){
            shift = 0;
            losshift = plos->Vlos*1e3/C_c;
  //fprintf(stderr,"vlos = %e \n",plos->Para.Vx);
  //losshift = 0.;
            sqrtT = sqrt(plos->Para.T);

            for(iatom=0; iatom<Input->Natom; iatom++){

              for(iline=0;iline<Atom[iatom].Nline;iline++){
                ilines = iline+shift;
                itrans = Atom[iatom].iout[iline];

                lcenter = 1e10*Atom[iatom].TR[itrans].lambda \
                    *(1.-losshift);
                Dopp = Atom[iatom].cDopp*sqrtT \
                    *Atom[iatom].TR[itrans].lambda*1e10;

                indx = 0;
  //fprintf(stderr,"atom = %d line = %d %e\n",iatom,iline,lcenter);

                for(ispec=0;ispec<Input->Nspec;ispec++){

                  if(lcenter<Input->Spec[ispec].range[0]-clambda \
                      || lcenter>Input->Spec[ispec].range[1]+clambda){
                    indx += Input->Spec[ispec].Nl*4;
                    continue;
                  }
  //fprintf(stderr,"atom = %d line = %d %e\n",iatom,iline,lcenter);

                  for(ilambda=0;ilambda<Input->Spec[ispec].Nl;ilambda++){
                    dlambda = (Input->Spec[ispec].Lambda[ilambda] \
                        -lcenter)/Dopp;
                    prof = exp(-dlambda*dlambda)/C_sqrtpi/Dopp;
                    dptmp = Output->specloc[pgrid->ipspec]+indx+ilambda;
                    *dptmp += prof*plos->Line[ilines].Stk[0];
                    *(dptmp+Input->Spec[ispec].Nl) += \
                        prof*plos->Line[ilines].Stk[1];
                    *(dptmp+Input->Spec[ispec].Nl*2) += \
                        prof*plos->Line[ilines].Stk[2];

                    if(Input->Nstk==4){
                      *(dptmp+Input->Spec[ispec].Nl*3) += \
                          plos->Line[ilines].Stk[3] \
                          *prof*(2.*dlambda/Dopp)*1e10;

                      if(Input->Verbose >= 2  && (Mpi->rank) == 0){
                        fprintf(stderr,"%d %e %e %e %e %e \n", ilambda, \
                            dlambda, *dptmp, *(dptmp+Input->Spec[ispec].Nl), \
                            *(dptmp+Input->Spec[ispec].Nl*2), \
                            *(dptmp+Input->Spec[ispec].Nl*3));
                      }
                    }


                  }
                  indx += Input->Spec[ispec].Nl*4;
                }
              }
              shift += Atom[iatom].Nline;
            }
          }
        }
      }


      MPI_Reduce(Output->synloc[0], Output->los[ilos].syn[0], itmp3, MPI_DOUBLE, \
          MPI_SUM, 0, MPI_COMM_WORLD);

      if(Input->Nspec>0){
        MPI_Reduce(Output->specloc[0], Output->los[ilos].spec[0], \
            itmp2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      }


    }

    return 0;
}

extern double Gaussian(double Shift, double lambdaD){

    return 1./C_sqrtpi/lambdaD*exp(-Shift*Shift/lambdaD/lambdaD);

}

/*----------------------------------------------------------------------------*/

