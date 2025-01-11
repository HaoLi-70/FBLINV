
#include "INTERPOL.h"

/*----------------------------------------------------------------------------*/

    /*######################################################################
    
      revision log:

        30 Otc. 2024
          --- initial comment 

        27 Dec. 2024
          --- update: Capable of inputting the magnetic field for the 
                      grid mode.  
    
    ######################################################################*/

/*----------------------------------------------------------------------------*/

extern double Interpol_Linear3D(double ***Data, int i, int j, int k, \
    double ir, double jr, double kr){
  
    /*######################################################################
      Purpose:
        3D linear Interpolation.
      Record of revisions:
        26 Sep. 2019.
      Input parameters:
        Data[][][], 3D data cube.
        i, j, k, the grid point close the interpolation position.
        ir, jr, kr, the ratios.
      Return:
        return the interpolated value.
    ######################################################################*/
  
    double tmp1, tmp2, tmp3, tmp4;
    double irp = 1-ir;
    
    tmp1 = irp*Data[i][j][k]+ir*Data[i+1][j][k];
    tmp2 = irp*Data[i][j+1][k]+ir*Data[i+1][j+1][k];
    tmp3 = (1-jr)*tmp1+jr*tmp2;
    tmp1 = irp*Data[i][j][k+1]+ir*Data[i+1][j][k+1];
    tmp2 = irp*Data[i][j+1][k+1]+ir*Data[i+1][j+1][k+1];
    tmp4 = (1-jr)*tmp1+jr*tmp2;
    
    return (1-kr)*tmp3+kr*tmp4;
}

/*----------------------------------------------------------------------------*/

extern void Get_Para(STRUCT_ATMO *Atmo, STRUCT_ATOM *Atom, \
    STRUCT_SYN *Syn, STRUCT_INPUT *Input, STRUCT_MPI *Mpi){
  
    /*######################################################################
      Purpose:
        get model parameters of a grid.
      Record of revisions:
        27 Dec. 2024.
      Input parameters:
        Atmo, a structure with the coronal model.
        Position, the position.
      Output parameters:
        Para, a structure with the model paramter.
    ######################################################################*/

    // constant
    //const double H_Density = C_H_Ratio/C_mH;
    //const double He_Density = C_He_Ratio/C_mHe*2;
    //const double ratio = (H_Density+He_Density)/H_Density;
    //ratio = C_E2H
    //fprintf(stderr,"aa = %e %e \n",ratio,C_E2H);

    STRUCT_GRID *pgrid;
    STRUCT_LOS *plos;
    int indx, igrid, ir = -1, itheta = -1, iphi = -1, ilos, iatom;
    double Ratio_R = 0.0, Ratio_Theta = 0.0, Ratio_Phi = 0.0, lg10T, DT;
    double Density, Br = 0.0, Btheta = 0.0, Bphi = 0.0, Phi;
    double Vr = 0.0, Vtheta = 0.0, Vphi = 0.0, V, ThetaV, PhiV;

    for(igrid=0; igrid<Mpi->ngrids; igrid++){
      pgrid = Syn->Grids+igrid;
     
      if((pgrid->R<Atmo->R[0]) || \
          (pgrid->R>Atmo->R[Atmo->nR-1])){

        const char *routine_name = "Get_Para";  
        fprintf(stderr,"r = %e %e %e \n", pgrid->R, Atmo->R[0], \
            Atmo->R[Atmo->nR-1]);
        Error(enum_error, routine_name, "R is out of the model range! \n");
      }else{
        for(indx=0; indx<Atmo->nR-1; indx++){
          if((pgrid->R>=Atmo->R[indx]) && \
              (pgrid->R<=Atmo->R[indx+1])){

            ir = indx;
            Ratio_R = (pgrid->R-Atmo->R[ir]) \
                /(Atmo->R[ir+1]-Atmo->R[ir]);
            break;
          }
        }
      }
 
      if(pgrid->Theta<Atmo->Theta[0]){
        itheta = 0;
        Ratio_Theta = 0.0;
      }else if(pgrid->Theta>Atmo->Theta[Atmo->nTheta-1]){
        itheta = Atmo->nTheta-2;
        Ratio_Theta = 1.0;
      }else{
        for(indx=0; indx<Atmo->nTheta-1; indx++){
          if((pgrid->Theta>=Atmo->Theta[indx]) &&\
              (pgrid->Theta<=Atmo->Theta[indx+1])){
            itheta = indx;
            Ratio_Theta = (pgrid->Theta-Atmo->Theta[itheta])\
                /(Atmo->Theta[itheta+1]-Atmo->Theta[itheta]);
            break;
          }
        }
      }

      for(ilos=0;ilos<Input->nlos;ilos++){
        plos = pgrid->los+ilos;

        Phi = pgrid->Phi+Input->los[ilos];
        if(Phi<0) Phi += C_2Pi;
        if(Phi>C_2Pi) Phi -= C_2Pi;

        if((Phi<Atmo->Phi[0])){
          iphi = 0;
          Ratio_Phi = 0.0;
        }else if(Phi>Atmo->Phi[Atmo->nPhi-1]){
          iphi = Atmo->nPhi-2;
          Ratio_Phi = 1.0;
        }else{
              
          for(indx=0; indx<Atmo->nPhi-1; indx++){
            if(Phi >= Atmo->Phi[indx] && Phi <= Atmo->Phi[indx+1]){
              iphi = indx;
              Ratio_Phi = (Phi-Atmo->Phi[iphi]) \
                      /(Atmo->Phi[iphi+1]-Atmo->Phi[iphi]);
              break;
            }
          }
        }

        Density = Interpol_Linear3D(Atmo->rho, ir, \
            itheta, iphi, Ratio_R, Ratio_Theta, Ratio_Phi);

        switch (Atmo->rhotype){
          case 0:
            plos->Para.nH = pow(10,Density);
            plos->Para.ne = plos->Para.nH/C_H2E;
            break;

          case 1:
            plos->Para.nH = Density;
            plos->Para.ne = plos->Para.nH/C_H2E;
            break;

          case 2:
            plos->Para.ne = pow(10,Density);
            plos->Para.nH = plos->Para.ne*C_H2E;
            break;

          case 3:
            plos->Para.ne = Density;
            plos->Para.nH = plos->Para.ne*C_H2E;
            break;
        
          default:
            break;
        }

        plos->Para.T = Interpol_Linear3D(Atmo->T, ir, itheta, \
            iphi, Ratio_R, Ratio_Theta, Ratio_Phi);

        lg10T = log10(plos->Para.T);

        switch (Atmo->btype){
          case 0:
            Br = Interpol_Linear3D(Atmo->B1, ir, itheta, \
                iphi, Ratio_R, Ratio_Theta, Ratio_Phi);
            Btheta = Interpol_Linear3D(Atmo->B2, ir, \
                itheta, iphi, Ratio_R, Ratio_Theta, Ratio_Phi);
            Bphi = Interpol_Linear3D(Atmo->B3, ir, itheta, \
                iphi, Ratio_R, Ratio_Theta, Ratio_Phi);
            plos->Para.B = sqrt(Br*Br+Btheta*Btheta+Bphi*Bphi);
            plos->Para.ThetaB = acos(Br/plos->Para.B);

            plos->Para.PhiB = atan2(Bphi, Btheta); //atan2(y,x)

            if(plos->Para.PhiB<0) plos->Para.PhiB += C_2Pi;
  
            break;

          case 1:
            plos->Para.B = \
                Interpol_Linear3D(Atmo->B1, ir, itheta, \
                iphi, Ratio_R, Ratio_Theta, Ratio_Phi);
            plos->Para.ThetaB = \
                Interpol_Linear3D(Atmo->B2, ir, \
                itheta, iphi, Ratio_R, Ratio_Theta, Ratio_Phi);
            plos->Para.PhiB = \
                Interpol_Linear3D(Atmo->B3, ir, itheta, \
                iphi, Ratio_R, Ratio_Theta, Ratio_Phi);

            if(plos->Para.PhiB<0) plos->Para.PhiB += C_2Pi;

            Br = plos->Para.B*cos(plos->Para.ThetaB);
            Btheta = plos->Para.B*sin(plos->Para.ThetaB) \
                *cos(plos->Para.PhiB);
            Bphi = plos->Para.B*sin(plos->Para.ThetaB) \
                *sin(plos->Para.PhiB);   

            break;
        
          default:
            break;
        }

        if(Input->Nspec>0){

          switch (Atmo->vtype){
            case 0:
              Vr = Interpol_Linear3D(Atmo->V1, ir, itheta, \
                  iphi, Ratio_R, Ratio_Theta, Ratio_Phi);
              Vtheta = Interpol_Linear3D(Atmo->V2, ir, \
                  itheta, iphi, Ratio_R, Ratio_Theta, Ratio_Phi);
              Vphi = Interpol_Linear3D(Atmo->V3, ir, itheta, \
                  iphi, Ratio_R, Ratio_Theta, Ratio_Phi);
              V = sqrt(Vr*Vr+Vtheta*Vtheta+Vphi*Vphi);
              ThetaV = acos(Vr/V);
              PhiV = atan2(Vphi, Vtheta);

              break;

            case 1:
              V = Interpol_Linear3D(Atmo->V1, ir, itheta, \
                  iphi, Ratio_R, Ratio_Theta, Ratio_Phi);
              ThetaV = Interpol_Linear3D(Atmo->V2, ir, \
                  itheta, iphi, Ratio_R, Ratio_Theta, Ratio_Phi);
              PhiV = Interpol_Linear3D(Atmo->V3, ir, itheta, \
                  iphi, Ratio_R, Ratio_Theta, Ratio_Phi);

              if(PhiV<0) plos->Para.PhiB += C_2Pi;

              Vr = V*cos(ThetaV);
              Vtheta = V*sin(ThetaV)*cos(PhiV);
              Vphi = V*sin(ThetaV)*sin(PhiV);   

              break;

            default:
              break;
          }

          plos->Blos = Br*sin(pgrid->Theta)*cos(pgrid->Phi) \
              +Btheta*sin(pgrid->Theta+C_Pi2)*cos(pgrid->Phi) \
              +Bphi*sin(pgrid->Phi);

          plos->Vlos = Vr*sin(pgrid->Theta)*cos(pgrid->Phi) \
              +Vtheta*sin(pgrid->Theta+C_Pi2)*cos(pgrid->Phi) \
              +Vphi*sin(pgrid->Phi);

  //Vr, Vtheta, Vphi, V, ThetaV, PhiV;

        }

        for(iatom=0; iatom<Input->Natom; iatom++){
          if((lg10T<Atom[iatom].Ion.frac[0][0]) || \
              (lg10T>Atom[iatom].Ion.frac[0][Atom[iatom].Ion.NT-1])){

            plos->Para.Ion[iatom] = 0;

          }else{     
            for(indx=0; indx<Atom[iatom].Ion.NT-1; indx++){
              if((lg10T>=Atom[iatom].Ion.frac[0][indx]) \
                  && (lg10T<=Atom[iatom].Ion.frac[0][indx+1])){
                DT = (lg10T-Atom[iatom].Ion.frac[0][indx]) \
                    /(Atom[iatom].Ion.frac[0][indx+1] \
                    -Atom[iatom].Ion.frac[0][indx]);
                plos->Para.Ion[iatom] = ((1.0-DT) \
                  *Atom[iatom].Ion.frac[1][indx] \
                  +DT*Atom[iatom].Ion.frac[1][indx+1])*Atom[iatom].Abund;
                break;
              }
            }
          }
        }
      }
    }

    if(Input->Mode == 2){
      for(igrid=0; igrid<Mpi->ngrids; igrid++){
        pgrid = Syn->Grids+igrid;

        plos->Para.B = Input->Bvec[0];
        plos->Para.ThetaB = Input->Bvec[1];
        plos->Para.PhiB = Input->Bvec[2];

        Br = plos->Para.B*cos(plos->Para.ThetaB);
        Btheta = plos->Para.B*sin(plos->Para.ThetaB) \
            *cos(plos->Para.PhiB);
        Bphi = plos->Para.B*sin(plos->Para.ThetaB) \
            *sin(plos->Para.PhiB);   

        plos->Blos = Br*sin(pgrid->Theta)*cos(pgrid->Phi) \
            +Btheta*sin(pgrid->Theta+C_Pi2)*cos(pgrid->Phi) \
            +Bphi*sin(pgrid->Phi);

      }
    }

    return;
}

/*----------------------------------------------------------------------------*/

extern int Ion_fraction(STRUCT_ATOM *Atom, STRUCT_SYN *Syn, \
    STRUCT_INPUT *Input, STRUCT_MPI *Mpi){

    /*######################################################################
      Purpose:
        get ion fraction.
      Record of revisions:
        12 Aug. 2023.
      Input parameters:
        Atom, a structure with the atomic information.
        Syn, a structure with forward synthesis.
        Input, a structure with the input information.
        Mpi, a structure with the Mpi information.
      Output parameters:
    ######################################################################*/

    STRUCT_GRID *pgrid;
    STRUCT_LOS *plos;
    int indx, igrid, iatom, ilos;
    double lg10T, DT;

    for(igrid=0; igrid<Mpi->ngrids; igrid++){
      pgrid = Syn->Grids+igrid;
      for(ilos=0;ilos<Input->nlos;ilos++){
        plos = pgrid->los+ilos;

        lg10T = log10(plos->Para.T);
        for(iatom=0; iatom<Input->Natom; iatom++){
          if((lg10T<Atom[iatom].Ion.frac[0][0]) || \
              (lg10T>Atom[iatom].Ion.frac[0][Atom[iatom].Ion.NT-1])){

            plos->Para.Ion[iatom] = 0;

          }else{     
            for(indx=0; indx<Atom[iatom].Ion.NT-1; indx++){
              if((lg10T>=Atom[iatom].Ion.frac[0][indx]) \
                  && (lg10T<=Atom[iatom].Ion.frac[0][indx+1])){
                DT = (lg10T-Atom[iatom].Ion.frac[0][indx]) \
                    /(Atom[iatom].Ion.frac[0][indx+1] \
                    -Atom[iatom].Ion.frac[0][indx]);
                plos->Para.Ion[iatom] = ((1.0-DT) \
                  *Atom[iatom].Ion.frac[1][indx] \
                  +DT*Atom[iatom].Ion.frac[1][indx+1])*Atom[iatom].Abund;
                break;
              }
            }
          }
        }
      }
    }

    return 0; 
}

/*----------------------------------------------------------------------------*/