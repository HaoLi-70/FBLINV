

#include "RECONSTRUCT.h"

/*----------------------------------------------------------------------------*/

    /*######################################################################
    
      revision log:

        30 Otc. 2024
          --- initial comment
          --- update: totally broken, will be updated later
    
    ######################################################################*/

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

/*  
    int igrid = 0, in, il, im, ipara, ilos;
    complex double tmp[4];
    STRUCT_GRID *pgrid;
    STRUCT_LOS *plos;

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

        plos->Para.Br = -creal(tmp[3])/pgrid->R/sin(pgrid->Theta);
        plos->Para.Bp = creal(tmp[2]);

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

        plos->Para.Br += creal(tmp[3])/pgrid->R;
        plos->Para.Bt = creal(tmp[2]);

        plos->Para.B = sqrt(plos->Para.Br \
            *plos->Para.Br+plos->Para.Bt*plos->Para.Bt \
            +plos->Para.Bp*plos->Para.Bp);
        plos->Para.ThetaB = acos(plos->Para.Br \
            /plos->Para.B);
        plos->Para.PhiB = atan2(plos->Para.Bp, \
            plos->Para.Bt);

      }else{
        if(Input->Para[2].invt) plos->Para.ThetaB = \
            creal(tmp[2]);
        if(Input->Para[3].invt) plos->Para.PhiB = \
            creal(tmp[3]);
      }
      
      if(Input->Para[0].invt) plos->Para.T = creal(tmp[0]);

      if(Input->Para[1].invt){
        plos->Para.nH = pow(10,creal(tmp[1]));
        plos->Para.ne = plos->Para.nH/C_H2E;
      }
    
    }
*/    
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
  
/*
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
          plos->Para.T = plos->Para_unperturb->T+tmp2;
          plos->Para.nH = plos->Para_unperturb->nH;
          plos->Para.ne = plos->Para_unperturb->ne;
          plos->Para.B = plos->Para_unperturb->B;
          plos->Para.ThetaB = plos->Para_unperturb->ThetaB;
          plos->Para.PhiB = plos->Para_unperturb->PhiB;
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

          plos->Para.nH = plos->Para_unperturb->nH*pow(10,tmp2);
          plos->Para.ne = plos->Para_unperturb->nH/C_H2E;
          plos->Para.T = plos->Para_unperturb->T;
          plos->Para.B = plos->Para_unperturb->B;
          plos->Para.ThetaB = plos->Para_unperturb->ThetaB;
          plos->Para.PhiB = plos->Para_unperturb->PhiB;        
          break;
            
        case 2:

          plos->Para.T = plos->Para_unperturb->T;
          plos->Para.nH = plos->Para_unperturb->nH;
          plos->Para.ne = plos->Para_unperturb->ne;

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

            plos->Para.Br = plos->Para_unperturb->Br \
                -tmp3/pgrid->R/sin(pgrid->Theta);
            plos->Para.Bp = plos->Para_unperturb->Bp+tmp2;

            plos->Para.B = sqrt(plos->Para.Br \
                *plos->Para.Br+plos->Para.Bt*plos->Para.Bt \
                +plos->Para.Bp*plos->Para.Bp);
            plos->Para.ThetaB = acos(plos->Para.Br \
                /plos->Para.B);
            plos->Para.PhiB = atan2(plos->Para.Bp, \
                plos->Para.Bt);

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
            plos->Para.ThetaB = plos->Para_unperturb->ThetaB+tmp2;
            plos->Para.B = plos->Para_unperturb->B;
            plos->Para.PhiB = plos->Para_unperturb->PhiB;
          }


          break;
            
        case 3:

          plos->Para.T = plos->Para_unperturb->T;
          plos->Para.nH = plos->Para_unperturb->nH;
          plos->Para.ne = plos->Para_unperturb->ne;

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

            plos->Para.Br += tmp3/pgrid->R;
            plos->Para.Bt = tmp2;

            plos->Para.B = sqrt(plos->Para.Br \
                *plos->Para.Br+plos->Para.Bt*plos->Para.Bt \
                +plos->Para.Bp*plos->Para.Bp);
            plos->Para.ThetaB = acos(plos->Para.Br \
                /plos->Para.B);
            plos->Para.PhiB = atan2(plos->Para.Bp, \
                plos->Para.Bt);


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

            plos->Para.PhiB = plos->Para_unperturb->PhiB+tmp2;
            plos->Para.B = plos->Para_unperturb->B;
            plos->Para.ThetaB = plos->Para_unperturb->ThetaB;

          }

          break;
            
        default:
          break;
      }
    }

    return ipara;
*/
    return 0;  
}

/*----------------------------------------------------------------------------*/