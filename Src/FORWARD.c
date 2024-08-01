
#include "FORWARD.h"

/*----------------------------------------------------------------------------*/

    /*######################################################################
    
      revision log:
        20 Aug. 2023.
    
    ######################################################################*/

/*----------------------------------------------------------------------------*/

extern int FORWARD(STRUCT_SYN *Syn, STRUCT_ATOM *Atom, \
    STRUCT_INPUT *Input, STR_FCTSG *fctsg, STRUCT_MPI *Mpi, bool sparse){
  
    /*######################################################################
      Purpose:
        forward synthesis of the polarization.
      Record of revisions:
        25 Oct. 2019
      Input parameters:
        Syn, a structure with forward synthesis.
        Atom, a structure with the atomic information.
        Input, a structure with the input information.
        fctsg, a structure with factoral, signs.
        Mpi, a structure with the Mpi information.
      Output parameters:
        Syn, a structure with forward synthesis.
    ######################################################################*/
  
    STRUCT_GRID *pgrid;

    double cont2, epsilon, Bp, Hu, norm;
    int i, j, Q, Qp, K, indx, igrid = 0, au;
    int itrans, iline, ilines, iatom, ilevel, iThom, shift;

    complex double **matrix_tmp,  *b;
    int *LUindx;

    for(igrid=0; igrid<Mpi->ngrids; igrid++){

      pgrid = Syn->Grids+igrid;

      if(sparse){
        if(Ran1(Mpi->idum+Mpi->rank) < Input->Per){ 
          
          if(Input->NThom > 0){
            for(iThom=0;iThom<Input->NThom;iThom++){
              pgrid->Thom[iThom].val[0] = pgrid->Thom[iThom].valsav[0];
              pgrid->Thom[iThom].val[1] = pgrid->Thom[iThom].valsav[1];
            }
          }
          if(Input->Nline > 0){
            for(iline=0;iline<Input->Nline;iline++){
              pgrid->Line[iline].Stk[0] = pgrid->Line[iline].Stksav[0];
              pgrid->Line[iline].Stk[1] = pgrid->Line[iline].Stksav[1];
              pgrid->Line[iline].Stk[2] = pgrid->Line[iline].Stksav[2];
            }
          }
          continue;
        }
      }
/*
pgrid->Para.Mag.PhiB = 0.;
pgrid->Para.Mag.ThetaB = 0.;
pgrid->Para.Bx = 10.;
*/
      Collisional_Rates(Atom, Input->Natom, pgrid->Para.T);

      Rotmat(pgrid->Para.Mag.PhiB, pgrid->Para.Mag.ThetaB, 0, Input->Dkmn, 2);

      Input->Tkq[0][0][0] = 1;
      for(i=0; i<3; i++){
        for(Qp=-Input->TQmax; Qp<=Input->TQmax; Qp++){
          Input->Tkq[i][2][Qp] = 0;
          for(Q=-2; Q<=2; Q++){
            Input->Tkq[i][2][Qp] += pgrid->T2Q[i][Q]*Input->Dkmn[2][Q][Qp];
          }
        }
        if(Input->Verbose >= 2  && (Mpi->rank) == 0){
          fprintf(stderr,"%d Input->Tkq new == %e %e %e %e %e \n", i, \
              creal(Input->Tkq[i][2][-2]), creal(Input->Tkq[i][2][-1]), \
              creal(Input->Tkq[i][2][0]), creal(Input->Tkq[i][2][1]), \
              creal(Input->Tkq[i][2][2]));
        }
      }

      shift = 0;
      for(iatom=0; iatom<Input->Natom; iatom++){

        if(Input->Verbose >= 2  && (Mpi->rank) == 0){
          fprintf(stderr,"\n ion = %e, ne = %e, T = %e\n",pgrid->Para.Ion[iatom], \
              pgrid->Para.ne,pgrid->Para.T);

          fprintf(stderr,"collision rates = %e \n",Atom[iatom].col->Rates[1][0]);

          fprintf(stderr,"bolzman %e", \
              exp(-C_h*Atom[iatom].TR[0].nu/C_Kb/pgrid->Para.T) \
              *Atom[iatom].LV[1].deg/Atom[iatom].LV[0].deg);
          fprintf(stderr,"nline = %d %d\n", Atom[iatom].Nline, Input->Natom);
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
          Atom[iatom].TR[itrans].JKQ[0][0] = pgrid->Jkq[iatom].J00[itrans];
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
              *pgrid->Para.ne/Atom[iatom].TR[0].Aul;

          // notic Bp in Eq. 10.50 is not a Planck function, it is
          // actually the coefficient due to the conversions between
          // Blu to Aul and CI to CS, finally turned to the 
          // Wien limit
          Bp = Wien_limit(Atom[iatom].TR[0].cplank1, \
              Atom[iatom].TR[0].cplank2, pgrid->Para.T);

          Atom[iatom].LV[1].Rho[0][0] = Atom[iatom].TR[0].const1 \
              *(epsilon*Bp+Atom[iatom].TR[0].JKQ[0][0])/(1+epsilon);

          if(Atom[iatom].LV[1].Kmax>=2){
            if(Input->Saturated){

              Atom[iatom].LV[1].Rho[2][0] = Atom[iatom].TR[0].const1 \
                  *Atom[iatom].TR[0].w[2]*Atom[iatom].TR[0].JKQ[2][0] \
                  /(1+epsilon);  

            }else{

              Hu = 0.879e7*Atom[iatom].LV[1].g*pgrid->Para.Mag.B \
                  /Atom[iatom].TR[0].Aul;

              for(Q=-2; Q<=2; Q++){
                Atom[iatom].LV[1].Rho[2][Q] = Atom[iatom].TR[0].const1 \
                    *(Atom[iatom].TR[0].w[2]*Atom[iatom].TR[0].JKQ[2][-Q] \
                    *fctsg->sg[Q])/(1+epsilon+Q*Hu*I);
              }  
            }
          }
     

        }else{

          matrix_tmp = SEE(Atom+iatom, fctsg, &(pgrid->Para), Input);

          b = (complex double *)VECTOR(1, Atom[iatom].nEq, enum_cplx, \
              false);
          LUindx = (int *)VECTOR(1, Atom[iatom].nEq, enum_int, false);

          for (i=1; i<=Atom[iatom].nEq; i++) {
            b[i] = -matrix_tmp[i][0];
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
                i = KQ_index(K,Q,Input->Kdelta);
                indx = Atom[iatom].eqindx[0]+i;
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
                  i = KQ_index(K,Q,Input->Kdelta);
                  indx = Atom[iatom].eqindx[ilevel]+i;
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
          norm += Atom[iatom].LV[ilevel].sqrt_deg*Atom[iatom].LV[ilevel].Rho[0][0];
        } 
        

        for(iline=0;iline<Atom[iatom].Nline;iline++){
          ilines = iline+shift;
          pgrid->Line[ilines].Stk[0] = 0;
          pgrid->Line[ilines].Stk[1] = 0;
          pgrid->Line[ilines].Stk[2] = 0;
          itrans = Atom[iatom].iout[iline];
          au = Atom[iatom].TR[Atom[iatom].iout[iline]].au;
          if(Input->Saturated){
            for(K=0;K<=Atom[iatom].LV[au].Kmax;K+=Input->Kdelta){
              for(i=0;i<3;i++){
                pgrid->Line[ilines].Stk[i] += Atom[iatom].TR[itrans].const2 \
                    *Atom[iatom].TR[iline].w[K] \
                    *creal(Input->Tkq[i][K][0]*Atom[iatom].LV[au].Rho[K][0]);
              }
            }

          }else{
            for(K=0;K<=Atom[iatom].LV[au].Kmax;K+=Input->Kdelta){
              for(Q=-K;Q<=K;Q++){
                for(i=0;i<3;i++){
                  pgrid->Line[ilines].Stk[i] += Atom[iatom].TR[itrans].const2 \
                      *Atom[iatom].TR[iline].w[K] \
                      *creal(Input->Tkq[i][K][Q]*Atom[iatom].LV[au].Rho[K][Q]);

                }
              }
            }
          }
          if(Input->OutputV){
            pgrid->Line[ilines].Stk[3] = Atom[iatom].TR[itrans].const2 \
                *pgrid->Para.Bx*C_Nul*Atom[iatom].TR[itrans].lambda \
                *Atom[iatom].TR[itrans].lambda/C_c \
                *(Atom[iatom].TR[itrans].geff \
                *Atom[iatom].LV[au].Rho[0][0] \
                +Atom[iatom].TR[itrans].delta*Atom[iatom].LV[au].Rho[2][0]) \
                *pgrid->Para.Ion[iatom]*pgrid->Para.nH/norm;
            if(Input->Verbose >= 2  && (Mpi->rank) == 0){
              fprintf(stderr, "V %e %e %e %e %e \n", pgrid->Line[ilines].Stk[3], \
                  Atom[iatom].TR[itrans].const2, \
                  pgrid->Para.Bx*C_Nul*Atom[iatom].TR[itrans].lambda \
                  *Atom[iatom].TR[itrans].lambda/C_c, 
                  creal(Atom[iatom].TR[itrans].geff \
                  *Atom[iatom].LV[au].Rho[0][0] \
                  +Atom[iatom].TR[itrans].delta*Atom[iatom].LV[au].Rho[2][0]), \
                  pgrid->Para.Ion[iatom]*pgrid->Para.nH/norm);
              fprintf(stderr, "geff %e delta %e au =%d\n", \
                  Atom[iatom].TR[itrans].geff, \
                  Atom[iatom].TR[itrans].delta, au);
              fprintf(stderr, "natom %e  %e norm =%e\n", \
                  pgrid->Para.Ion[iatom], \
                  pgrid->Para.nH, norm);

            }

          }

          pgrid->Line[ilines].Stk[0] *= pgrid->Para.Ion[iatom] \
              *pgrid->Para.nH/norm;
          if(Atom[iatom].TR[Atom[iatom].iout[iline]].M1){
            pgrid->Line[ilines].Stk[1] *= -pgrid->Para.Ion[iatom] \
                *pgrid->Para.nH/norm;
            pgrid->Line[ilines].Stk[2] *= -pgrid->Para.Ion[iatom] \
                *pgrid->Para.nH/norm;
          }else{
            pgrid->Line[ilines].Stk[1] *= pgrid->Para.Ion[iatom] \
                *pgrid->Para.nH/norm;
            pgrid->Line[ilines].Stk[2] *= pgrid->Para.Ion[iatom] \
                *pgrid->Para.nH/norm;
          }
          
          if(Input->Verbose >= 2  && (Mpi->rank) == 0){
            fprintf(stderr, "iline %d %d %d \n",iline,ilines,iatom);
            fprintf(stderr, "ion %e  \n", pgrid->Para.Ion[iatom] \
                *pgrid->Para.nH/norm);
            fprintf(stderr, "IQU %e %e %e \n", pgrid->Line[ilines].Stk[0], \
                pgrid->Line[ilines].Stk[1], pgrid->Line[ilines].Stk[2]);
            if(Input->OutputV){
              fprintf(stderr, "V %e \n",pgrid->Line[ilines].Stk[3]);
            }

          }
       
        }

        shift += Atom[iatom].Nline;
      }  

      for(j=0; j<Input->NThom; j++){
        cont2 = pgrid->Para.ne*Input->dx*Input->Thom[j].Const;
        pgrid->Thom[j].val[0] = (pgrid->Thom[j].Kr+pgrid->Thom[j].Kt) \
            *cont2;
        pgrid->Thom[j].val[1] = (pgrid->Thom[j].Kt-pgrid->Thom[j].Kr) \
            *cont2;
      }
    }

    return 0;
}

/*----------------------------------------------------------------------------*/

extern int Copy_Par(STRUCT_SYN *Syn, STRUCT_INPUT *Input, STRUCT_MPI *Mpi){


    STRUCT_GRID *pgrid;

    int igrid, iatom, iThom, iline;

    for(igrid=0; igrid<Mpi->ngrids; igrid++){
      pgrid = Syn->Grids+igrid;

      pgrid->Para.Tsav = pgrid->Para.T;
      pgrid->Para.nHsav = pgrid->Para.nH;
      pgrid->Para.nesav = pgrid->Para.ne;
      pgrid->Para.Magsav.B = pgrid->Para.Mag.B;
      pgrid->Para.Magsav.ThetaB = pgrid->Para.Mag.ThetaB;
      pgrid->Para.Magsav.PhiB = pgrid->Para.Mag.PhiB;
      pgrid->Para.Magsav.Br = pgrid->Para.Mag.Br;
      pgrid->Para.Magsav.Bt = pgrid->Para.Mag.Bt;
      pgrid->Para.Magsav.Bp = pgrid->Para.Mag.Bp;

      for(iatom=0; iatom<Input->Natom; iatom++){
        pgrid->Para.Ionsav[iatom] = pgrid->Para.Ion[iatom];
      }

      if(Input->NThom > 0){
        for(iThom=0;iThom<Input->NThom;iThom++){
          pgrid->Thom[iThom].valsav[0] = pgrid->Thom[iThom].val[0];
          pgrid->Thom[iThom].valsav[1] = pgrid->Thom[iThom].val[1];
        }
      }
      if(Input->Nline > 0){
        for(iline=0;iline<Input->Nline;iline++){
          pgrid->Line[iline].Stksav[0] = pgrid->Line[iline].Stk[0];
          pgrid->Line[iline].Stksav[1] = pgrid->Line[iline].Stk[1];
          pgrid->Line[iline].Stksav[2] = pgrid->Line[iline].Stk[2];
        }
      }
    }

    return 0;
}

/*----------------------------------------------------------------------------*/

extern int SFB_RECONSTRUTION(STRUCT_INPUT *Input, STRUCT_SYN *Syn, \
    STRUCT_MPI *Mpi){
  
    /*######################################################################
      Purpose:
        recconstruction the data from the SFB coefficients.
      Record of revisions:
        25 Oct. 2019
      Input parameters:
        Input, a structure with the inversion information.
        Mpi, a structure with the Mpi information.
      Output parameters:
        Input, a structure with the inversion information.
    ######################################################################*/
  
    int igrid = 0, in, il, im, ipara;
    complex double tmp[4];
    STRUCT_GRID *pgrid;

    for(igrid=0; igrid<Mpi->ngrids; igrid++){
      pgrid = Syn->Grids+igrid;
      for(ipara=0; ipara<4; ipara++){
        if(Input->Bpotential && ipara==2) break;
        if(!Input->Para[ipara].invt) continue;
        tmp[ipara] = 0;
        for(in=1; in<=Input->Para[ipara].N; in++){
          for(il=0; il<=Input->Para[ipara].L; il++){
            for(im=0; im<=il; im++){
              if(im == 0){
                tmp[ipara] += Input->Para[ipara].Coeff[in][il][im] \
                    *pgrid->sBessel[in][il]*pgrid->Legendre[il][im] \
                    *pgrid->Phiarray[im];
              }else{
                tmp[ipara] += Input->Para[ipara].Coeff[in][il][im] \
                    *pgrid->sBessel[in][il]*pgrid->Legendre[il][im] \
                    *pgrid->Phiarray[im]*2;
              }
            }
          }
        }
      }
      if(Input->Bpotential){
        tmp[2] = 0;
        tmp[3] = 0;

        ipara = 2;
        for(in=1; in<=Input->Para[ipara].N; in++){
          for(il=0; il<=Input->Para[ipara].L; il++){
            for(im=0; im<=il; im++){
              if(im == 0){
                tmp[2] += Input->Para[ipara].Coeff[in][il][im] \
                    *pgrid->sBesselD[in][il]*pgrid->Legendre[il][im] \
                    *pgrid->Phiarray[im];
              }else{
                tmp[2] += Input->Para[ipara].Coeff[in][il][im] \
                    *pgrid->sBesselD[in][il]*pgrid->Legendre[il][im] \
                    *pgrid->Phiarray[im]*2;
                tmp[3] += Input->Para[ipara].Coeff[in][il][im] \
                    *pgrid->sBessel[in][il]*pgrid->Legendre[il][im] \
                    *pgrid->Phiarray[im]*2*im*I;
              }
            }
          }
        }

        pgrid->Para.Mag.Br = -creal(tmp[3])/pgrid->R/sin(pgrid->Theta);
        pgrid->Para.Mag.Bp = creal(tmp[2]);

        tmp[2] = 0;
        tmp[3] = 0;

        ipara = 3;
        for(in=1; in<=Input->Para[ipara].N; in++){
          for(il=0; il<=Input->Para[ipara].L; il++){
            for(im=0; im<=il; im++){
              if(im == 0){
                tmp[2] += Input->Para[ipara].Coeff[in][il][im] \
                    *pgrid->sBesselD[in][il]*pgrid->Legendre[il][im] \
                    *pgrid->Phiarray[im];
                tmp[3] += Input->Para[ipara].Coeff[in][il][im] \
                    *pgrid->sBessel[in][il]*pgrid->LegendreD[il][im] \
                    *pgrid->Phiarray[im]*2;
              }else{
                tmp[2] += Input->Para[ipara].Coeff[in][il][im] \
                    *pgrid->sBesselD[in][il]*pgrid->Legendre[il][im] \
                    *pgrid->Phiarray[im]*2;
                tmp[3] += Input->Para[ipara].Coeff[in][il][im] \
                    *pgrid->sBessel[in][il]*pgrid->LegendreD[il][im] \
                    *pgrid->Phiarray[im]*2;
              }
            }
          }
        }

        pgrid->Para.Mag.Br += creal(tmp[3])/pgrid->R;
        pgrid->Para.Mag.Bt = creal(tmp[2]);

        pgrid->Para.Mag.B = sqrt(pgrid->Para.Mag.Br \
            *pgrid->Para.Mag.Br+pgrid->Para.Mag.Bt*pgrid->Para.Mag.Bt \
            +pgrid->Para.Mag.Bp*pgrid->Para.Mag.Bp);
        pgrid->Para.Mag.ThetaB = acos(pgrid->Para.Mag.Br \
            /pgrid->Para.Mag.B);
        pgrid->Para.Mag.PhiB = atan2(pgrid->Para.Mag.Bp, \
            pgrid->Para.Mag.Bt);

      }else{
        if(Input->Para[2].invt) pgrid->Para.Mag.ThetaB = \
            creal(tmp[2]);
        if(Input->Para[3].invt) pgrid->Para.Mag.PhiB = \
            creal(tmp[3]);
      }
      
      if(Input->Para[0].invt) pgrid->Para.T = creal(tmp[0]);

      if(Input->Para[1].invt){
        pgrid->Para.nH = pow(10,creal(tmp[1]));
        pgrid->Para.ne = pgrid->Para.nH/C_H2E;
      }
    
    }
    
    return 0;
}

/*----------------------------------------------------------------------------*/

extern int SFB_RECONSTR_DELTA(STRUCT_INPUT *Input, STRUCT_SYN *Syn, \
    int icoeff, STRUCT_MPI *Mpi){

    /*######################################################################
      Purpose:
        recconstruction the delta data from the SFB coefficients.
      Record of revisions:
        25 Oct. 2019
      Input parameters:
        Input, a structure with the inversion information.
        Mpi, a structure with the Mpi information.
      Output parameters:
        Input, a structure with the inversion information.
      return:
        the parameter index
    ######################################################################*/
  
    int igrid = 0, in, il, im, real, ipara;
    complex double tmp = 0, tmp1;
    double tmp2, tmp3;
    STRUCT_GRID *pgrid;

    Coeff2NLM(Input, icoeff, &ipara, &in, &il, &im, &real);
    
    for(igrid=0; igrid<Mpi->ngrids; igrid++){
      pgrid = Syn->Grids+igrid;
      tmp = 0;
        
      switch (ipara){
        case 0:
          if(im == 0){
            tmp = Input->perturb[ipara]*pgrid->sBessel[in][il] \
                *pgrid->Legendre[il][im]*pgrid->Phiarray[im];
          }else{
            tmp = Input->perturb[ipara]*pgrid->sBessel[in][il] \
                *pgrid->Legendre[il][im]*pgrid->Phiarray[im]*2;
          }

          if(real == 1){
            tmp2 = creal(tmp);
          }else{
            tmp2 = -cimag(tmp);
          }
          pgrid->Para.T = pgrid->Para.Tsav+tmp2;
          pgrid->Para.nH = pgrid->Para.nHsav;
          pgrid->Para.ne = pgrid->Para.nesav;
          pgrid->Para.Mag.B = pgrid->Para.Magsav.B;
          pgrid->Para.Mag.ThetaB = pgrid->Para.Magsav.ThetaB;
          pgrid->Para.Mag.PhiB = pgrid->Para.Magsav.PhiB;
          break;
            
        case 1:
          if(im == 0){
            tmp = Input->perturb[ipara]*pgrid->sBessel[in][il] \
                *pgrid->Legendre[il][im]*pgrid->Phiarray[im];
          }else{
            tmp = Input->perturb[ipara]*pgrid->sBessel[in][il] \
                *pgrid->Legendre[il][im]*pgrid->Phiarray[im]*2;
          }

          if(real == 1){
            tmp2 = creal(tmp);
          }else{
            tmp2 = -cimag(tmp);
          }

          pgrid->Para.nH = pgrid->Para.nHsav*pow(10,tmp2);
          pgrid->Para.ne = pgrid->Para.nHsav/C_H2E;
          pgrid->Para.T = pgrid->Para.Tsav;
          pgrid->Para.Mag.B = pgrid->Para.Magsav.B;
          pgrid->Para.Mag.ThetaB = pgrid->Para.Magsav.ThetaB;
          pgrid->Para.Mag.PhiB = pgrid->Para.Magsav.PhiB;        
          break;
            
        case 2:

          pgrid->Para.T = pgrid->Para.Tsav;
          pgrid->Para.nH = pgrid->Para.nHsav;
          pgrid->Para.ne = pgrid->Para.nesav;

          if(Input->Bpotential){

            if(im == 0){
              tmp = Input->Para[ipara].Coeff[in][il][im] \
                  *pgrid->sBesselD[in][il]*pgrid->Legendre[il][im] \
                  *pgrid->Phiarray[im];
              tmp1 = 0.0;

            }else{
              tmp = Input->Para[ipara].Coeff[in][il][im] \
                  *pgrid->sBesselD[in][il]*pgrid->Legendre[il][im] \
                  *pgrid->Phiarray[im]*2;
              tmp1 = Input->Para[ipara].Coeff[in][il][im] \
                  *pgrid->sBessel[in][il]*pgrid->Legendre[il][im] \
                  *pgrid->Phiarray[im]*2*im*I;
            }

            if(real == 1){
              tmp2 = creal(tmp);
              tmp3 = creal(tmp1);
            }else{
              tmp2 = -cimag(tmp);
              tmp3 = -cimag(tmp1);
            }

            pgrid->Para.Mag.Br = pgrid->Para.Magsav.Br \
                -tmp3/pgrid->R/sin(pgrid->Theta);
            pgrid->Para.Mag.Bp = pgrid->Para.Magsav.Bp+tmp2;

            pgrid->Para.Mag.B = sqrt(pgrid->Para.Mag.Br \
                *pgrid->Para.Mag.Br+pgrid->Para.Mag.Bt*pgrid->Para.Mag.Bt \
                +pgrid->Para.Mag.Bp*pgrid->Para.Mag.Bp);
            pgrid->Para.Mag.ThetaB = acos(pgrid->Para.Mag.Br \
                /pgrid->Para.Mag.B);
            pgrid->Para.Mag.PhiB = atan2(pgrid->Para.Mag.Bp, \
                pgrid->Para.Mag.Bt);

          }else{
            if(im == 0){
              tmp = Input->perturb[ipara]*pgrid->sBessel[in][il] \
                  *pgrid->Legendre[il][im]*pgrid->Phiarray[im];
            }else{
              tmp = Input->perturb[ipara]*pgrid->sBessel[in][il] \
                  *pgrid->Legendre[il][im]*pgrid->Phiarray[im]*2;
            }

            if(real == 1){
              tmp2 = creal(tmp);
            }else{
              tmp2 = -cimag(tmp);
            }
            pgrid->Para.Mag.ThetaB = pgrid->Para.Magsav.ThetaB+tmp2;
            pgrid->Para.Mag.B = pgrid->Para.Magsav.B;
            pgrid->Para.Mag.PhiB = pgrid->Para.Magsav.PhiB;
          }


          break;
            
        case 3:

          pgrid->Para.T = pgrid->Para.Tsav;
          pgrid->Para.nH = pgrid->Para.nHsav;
          pgrid->Para.ne = pgrid->Para.nesav;

          if(Input->Bpotential){
            if(im == 0){
              tmp += Input->perturb[ipara] \
                  *pgrid->sBesselD[in][il]*pgrid->Legendre[il][im] \
                  *pgrid->Phiarray[im];
              tmp1 += Input->perturb[ipara] \
                  *pgrid->sBessel[in][il]*pgrid->LegendreD[il][im] \
                  *pgrid->Phiarray[im]*2;
            }else{
              tmp += Input->perturb[ipara] \
                  *pgrid->sBesselD[in][il]*pgrid->Legendre[il][im] \
                  *pgrid->Phiarray[im]*2;
              tmp1 += Input->perturb[ipara] \
                  *pgrid->sBessel[in][il]*pgrid->LegendreD[il][im] \
                  *pgrid->Phiarray[im]*2;
            }

            if(real == 1){
              tmp2 = creal(tmp);
              tmp3 = creal(tmp1);
            }else{
              tmp2 = -cimag(tmp);
              tmp3 = -cimag(tmp1);
            }

            pgrid->Para.Mag.Br += tmp3/pgrid->R;
            pgrid->Para.Mag.Bt = tmp2;

            pgrid->Para.Mag.B = sqrt(pgrid->Para.Mag.Br \
                *pgrid->Para.Mag.Br+pgrid->Para.Mag.Bt*pgrid->Para.Mag.Bt \
                +pgrid->Para.Mag.Bp*pgrid->Para.Mag.Bp);
            pgrid->Para.Mag.ThetaB = acos(pgrid->Para.Mag.Br \
                /pgrid->Para.Mag.B);
            pgrid->Para.Mag.PhiB = atan2(pgrid->Para.Mag.Bp, \
                pgrid->Para.Mag.Bt);


          }else{
            if(im == 0){
              tmp = Input->perturb[ipara]*pgrid->sBessel[in][il] \
                  *pgrid->Legendre[il][im]*pgrid->Phiarray[im];
            }else{
              tmp = Input->perturb[ipara]*pgrid->sBessel[in][il] \
                  *pgrid->Legendre[il][im]*pgrid->Phiarray[im]*2;
            }

            if(real == 1){
              tmp2 = creal(tmp);
            }else{
              tmp2 = -cimag(tmp);
            }

            pgrid->Para.Mag.PhiB = pgrid->Para.Magsav.PhiB+tmp2;
            pgrid->Para.Mag.B = pgrid->Para.Magsav.B;
            pgrid->Para.Mag.ThetaB = pgrid->Para.Magsav.ThetaB;

          }

          break;
            
        default:
          break;
      }
    }

    return ipara;
}

/*----------------------------------------------------------------------------*/

extern int Grid2Pixel(STRUCT_SYN *Syn, STRUCT_ATOM *Atom, \
    STRUCT_INPUT *Input, STRUCT_OUT *Output, STRUCT_MPI *Mpi){
  
    /*######################################################################
      Purpose:
        convert the synthetic data to the observation style.
      Record of revisions:
        23 Aug. 2023
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
        ilambda, shift;
    STRUCT_GRID *pgrid;

    int itmp1 = Input->Nstk*Input->Nline;
    int itmp2 = Syn->npspec*Input->Nl*Input->Nstk;
    int itmp3 = Syn->npixels*(itmp1+2*Input->NThom);

    // initialize to 0;
    memset(Output->synloc[0], 0, itmp3*sizeof(double));

    if(Input->Nspec>0){      
      memset(Output->specloc[0], 0, itmp2*sizeof(double));
    }

    double clambda = 25.0;
    double losshift, lcenter;
    double sqrtT, Dopp, prof;
    int indx;
    double *dptmp, dlambda;

    for(igrid=0; igrid<Mpi->ngrids; igrid++){
      pgrid = Syn->Grids+igrid;
//fprintf(stderr,"grid %d \n",igrid);
      for(iline=0; iline<Input->Nline; iline++){
        for(istk=0; istk<Input->Nstk; istk++){
          Output->synloc[pgrid->ipixel][istk+Input->Nstk*iline] += \
              pgrid->Line[iline].Stk[istk];            
        }
      }

      for(iThom=0; iThom<Input->NThom; iThom++){
        dptmp = Output->synloc[pgrid->ipixel]+2*iThom+itmp1;
        *dptmp += pgrid->Thom[iThom].val[0];
        *(dptmp+1) += pgrid->Thom[iThom].val[1];
      }

      // compute the spectrum
      if(Input->Nspec>0){
        if(pgrid->ipspec>=0){
          shift = 0;
          losshift = pgrid->Para.Vx*1e3/C_c;
//fprintf(stderr,"vlos = %e \n",pgrid->Para.Vx);
//losshift = 0.;
          sqrtT = sqrt(pgrid->Para.T);

          for(iatom=0; iatom<Input->Natom; iatom++){

            for(iline=0;iline<Atom[iatom].Nline;iline++){
              ilines = iline+shift;
              itrans = Atom[iatom].iout[iline];

              lcenter = 1e10*Atom[iatom].TR[itrans].lambda \
                  *(1.-losshift);
              Dopp = Atom[iatom].cDopp*sqrtT*Atom[iatom].TR[itrans].lambda*1e10;

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
                  dlambda = (Input->Spec[ispec].Lambda[ilambda]-lcenter)/Dopp;
                  prof = exp(-dlambda*dlambda)/C_sqrtpi/Dopp;
                  dptmp = Output->specloc[pgrid->ipspec]+indx+ilambda;
                  *dptmp += prof*pgrid->Line[ilines].Stk[0];
                  *(dptmp+Input->Spec[ispec].Nl) += \
                      prof*pgrid->Line[ilines].Stk[1];
                  *(dptmp+Input->Spec[ispec].Nl*2) += \
                      prof*pgrid->Line[ilines].Stk[2];

                  if(Input->Nstk==4){
                    *(dptmp+Input->Spec[ispec].Nl*3) += \
                        pgrid->Line[ilines].Stk[3] \
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


    MPI_Reduce(Output->synloc[0], Output->syntot[0], itmp3, MPI_DOUBLE, \
        MPI_SUM, 0, MPI_COMM_WORLD);

    if(Input->Nspec>0){
      MPI_Reduce(Output->specloc[0], Output->spectot[0], \
          itmp2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    return 0;
}

extern double Gaussian(double Shift, double lambdaD){

    return 1./C_sqrtpi/lambdaD*exp(-Shift*Shift/lambdaD/lambdaD);
}

/*----------------------------------------------------------------------------*/

extern void Coeff2NLM(STRUCT_INPUT *Input, int icoeff, int *ipara, \
    int *in, int *il, int *im, int *real){
  
    /*######################################################################
      Purpose:
        convert SFB coefficients index to the orders N, L, M.
      Record of revisions:
        10 Sept. 2021.
      Input parameters:
        Input, a structure with the inversion information.
        icoeff, the index of the coefficient.
      Output parameters:
        ipara, in, il, im, real (1:real; 0:imagenary), the orders.
    ######################################################################*/
  
    int tmp_coef = icoeff, indxpara;
    for(indxpara=0;indxpara<4;indxpara++){
      if(!Input->Para[indxpara].invt) continue;
      if (tmp_coef>=Input->Para[indxpara].NLsquare) {
        tmp_coef -= Input->Para[indxpara].NLsquare;
      }else{
        *ipara = indxpara;
        break;
      }      
    }
    
    *in = tmp_coef/Input->Para[*ipara].Lsquare+1;
    *im = tmp_coef%Input->Para[*ipara].Lsquare;
    *il = (int)sqrt((double)(*im));
    *im -= (*il)*(*il);
    
    if(*im == 0){
      *real = 1;
    }else{
      *real = (*im)%2;
      *im = (*im+1)/2;
    }
    
    return;
}

/*----------------------------------------------------------------------------*/

extern void NLM2Coeff(STRUCT_INPUT *Input, int ipara, int in, \
    int il, int im, int real, int *icoeff){
  
    /*######################################################################
      Purpose:
        convert the orders N, L, M to SFB coefficients index.
      Record of revisions:
        10 Sept. 2021.
      Input parameters:
        Input, a structure with the inversion information.
        ipara, in, il, im, real, the orders.
      Output parameters:
        icoeff, the index of the coefficient.
    ######################################################################*/
  
    int i;
    *icoeff = 0;
    for(i=0; i<ipara; i++){
      *icoeff += Input->Para[i].NLsquare;
    }
    
    *icoeff += (in-1)*Input->Para[ipara].Lsquare+il*il+im*2;
    
    if(real > 0 && im > 0){
      *icoeff = *icoeff-1;
    }
    
    return;
}

/*----------------------------------------------------------------------------*/
