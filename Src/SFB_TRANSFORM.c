
#include "SFB_TRANSFORM.h"

/*----------------------------------------------------------------------------*/

    /*######################################################################
     
     revision log:
        8 Sept. 2021.
     
    ######################################################################*/

/*----------------------------------------------------------------------------*/

static void INTEG_REAL(int il, double **Zeros, double R_Max, double R_Min, \
    int Kmax, int Int_Flag, double a, double b, \
    enum boundary_condition BC_Flag, double *Integ);

static void Legendre_Array(int L_Max, double *Theta, int nTheta, \
    double ***har_coffi);

static void Cal_Delta(double *Radius, int nR, double *deltaR);

static void INT_COMP_BC(int il, double **Zeros, double R_Max, \
    double R_Min, int Kmax, int Int_Flag, complex double a, \
    complex double b, enum boundary_condition BC_Flag, \
    complex double *Integ);

/*----------------------------------------------------------------------------*/

static void INTEG_REAL(int il, double **Zeros, double R_Max, double R_Min, \
    int Kmax, int Int_Flag, double a, double b, \
    enum boundary_condition BC_Flag, double *Integ){
    
    /*######################################################################
      Purpose:
        integration of spherical Bessel function (ordr l) from 0 to R_MIN
            (zero-derivative boundary condition).
      Record of revisions:
        20 Otc. 2019.
      Input parameters:
        il, order L.
        Zeros, the zeros of the derivatives of spherical Bessel functions.
        R_MAX, the max value of the radius.
        R_MIN, the min value of the radius
        Kmax, max order of K.
        Int_Flag, a flag. Int_Flag = 1: the field is described by a*R_MIN+b
            (constant); Int_Flag = 2: a*r+b (linear); Int_Flag = 3: a*r^2+b
            (parabolic); Int_Flag = 4: a*r*2-a*r+b(parabolic).
        a, b, coefficients used to describe the field between 0 and R_MIN.
        BC_Flag, boundary condition flag. BC_Flag = enum_zeros for 
            zero-value boundary condition, and = enum_deri_zeros for 
            zero-derivative boundary condition
      Input parameters:
        Integ, the integration.
    ######################################################################*/

    int ik, ir;
    double JJ[il+1], JJp[il+1], norm, K_ln, Rtmp, dR = R_Min/100;
    double field = a*R_Min+b;
    
    for(ik=1; ik<=Kmax; ik++){
      Integ[ik] = 0;
      Spherical_Bessel(Zeros[il][ik], il+1, JJ, JJp);
        
      if(BC_Flag == enum_zeros){
        norm = 2./R_Max/R_Max/R_Max/JJ[il+1]/JJ[il+1];
      }else{
        if(il == 0){
          if(ik == 1){
            norm = 3./R_Max/R_Max/R_Max;
          }else{
            norm = 2./R_Max/R_Max/R_Max/JJ[il]/JJ[il];
          }
        }else{
          norm = 2./R_Max/R_Max/R_Max/JJ[il]/JJ[il]/(1-il*(il+1) \
              /Zeros[il][ik]/Zeros[il][ik]);
        }
      }
        
      K_ln = Zeros[il][ik]/R_Max;
      for(ir=1; ir<=100; ir++){
        Rtmp = (ir-0.5)*dR;
        Spherical_Bessel(K_ln*Rtmp, il, JJ, JJp);
        if(Int_Flag == 1){
          Integ[ik] += Rtmp*Rtmp*dR*JJ[il]*norm*field;
        }else if(Int_Flag == 2){
          Integ[ik] += Rtmp*Rtmp*dR*JJ[il]*norm*(a*Rtmp+b);
        }else if(Int_Flag == 3){
          Integ[ik] += Rtmp*Rtmp*dR*JJ[il]*norm*(a*Rtmp*Rtmp+b);
        }else if(Int_Flag == 4){
          Integ[ik] += Rtmp*Rtmp*dR*JJ[il] *norm*(a*Rtmp*Rtmp-a*Rtmp+b);
        }
      }
    }
    return;
}

/*----------------------------------------------------------------------------*/

extern void SB_DECOMP_SINGLE(double *Field_Array, double *Radius, \
    int Num_R, double **Zeros, int Kmax, int il, int Int_Flag, \
    enum boundary_condition BC_Flag, double *Coef_SB){
    
    /*######################################################################
      Purpose:
        spherical Bessel decomposition (ordr l) of a data array.
      Record of revisions:
        22 Otc. 2019.
      Input parameters:
        Data_Array, a data array.
        Radius, a radius array saved the position of data.
        Num_R, the length of the data and radius arrays.
        Zeros, the zeros of the derivatives of spherical Bessel functions.
        Kmax, max order of K.
        il, order L.
        Int_Flag, a flag. Int_Flag = : the field is described by a*R_MIN+b
            (constant); Int_Flag = 2: a*r+b (linear); Int_Flag = 3: a*r^2+b
            (parabolic); Int_Flag = 4: a*r*2-a*r+b(parabolic).
        BC_Flag, boundary condition flag. BC_Flag = enum_zeros for 
            zero-value boundary condition, and = enum_deri_zeros for 
            zero-derivative boundary condition
      Output parameters:
        Coef_SB, output the spherical Bessel coefficients.
    ######################################################################*/

    int ik, ir;
    double a=0., b=0.;
    double R_Max = Radius[Num_R-1], R_Min = Radius[0];
    double JJ[il+1], JJp[il+1], norm, K_ln;
    double *Integ = (double *)VECTOR(1, Kmax, enum_dbl, true);
    double *dr = (double *)malloc(Num_R*sizeof(double));
    Cal_Delta(Radius, Num_R, dr);

    if(R_Min > 0 && Int_Flag>=1  && Int_Flag<=4){
      if(Int_Flag==1 || Int_Flag==2){
        a = (Field_Array[1]-Field_Array[0])/(Radius[1]-Radius[0]);
        b = Field_Array[0]-a*Radius[0];
      }else if(Int_Flag==3){
        a = (Field_Array[1]-Field_Array[0])/(Radius[1]-Radius[0])/2;
        b = Field_Array[0]-a*Radius[0]*Radius[0];
      }else if(Int_Flag==4){
        a = (Field_Array[1]-Field_Array[0])/(Radius[1]-Radius[0]);
        b = Field_Array[0];
      }
      INTEG_REAL(il, Zeros, R_Max, R_Min, Kmax, Int_Flag, a, b, BC_Flag, \
          Integ);
      for(ik=1; ik<=Kmax; ik++){
        Coef_SB[ik] = Integ[ik];
      }
    }else{
      for(ik=1; ik<=Kmax; ik++){
        Coef_SB[ik] = Integ[ik];
      }
    }
    
    for(ik=1; ik<=Kmax; ik++){
      Spherical_Bessel(Zeros[il][ik], il+1, JJ, JJp);
      K_ln = Zeros[il][ik]/R_Max;

      if(BC_Flag == enum_zeros){
        norm = 2./R_Max/R_Max/R_Max/JJ[il+1]/JJ[il+1];
      }else{
        if(il == 0){
          if(ik == 1){
            norm = 3./R_Max/R_Max/R_Max;
          }else{
            norm = 2./R_Max/R_Max/R_Max/JJ[il]/JJ[il];
          }
        }else{
          norm = 2./R_Max/R_Max/R_Max/JJ[il]/JJ[il]/(1-il*(il+1) \
              /Zeros[il][ik]/Zeros[il][ik]);
        }
      }
        
      for(ir=0; ir<Num_R; ir++){
        Spherical_Bessel(K_ln*Radius[ir], il, JJ, JJp);
        Coef_SB[ik] += Field_Array[ir]*Radius[ir]*Radius[ir] \
            *dr[ir]*JJ[il]*norm;
      }
    }
    
    FREE_VECTOR(Integ, 1, enum_dbl);
    free(dr);
    
    return;
}

/*----------------------------------------------------------------------------*/

extern void SB_RECONSTRUCT_SINGLE(double *Coef_SB, int Kmax, int il, \
    double **Zeros, double R_Max, double *Radius, int Num_R, \
    double *field){
    
    /*######################################################################
      Purpose:
        reconstruct a data array from a spherical Bessel (ordr l)
            coefficient array
      Record of revisions:
        22 Otc. 2019.
      Input parameters:
        Coef_SB, a spherical Bessel coefficient array.
        Kmax, max order of K (e.g. the length of the coefficient array).
        il, order L.
        Zeros, the zeros of the spherical Bessel functions.
        R_Max, max of the radius used in the decompostion.
        Radius, a radius array saved the position of data.
        Num_R, the length of the radius array.
      Output parameters:
        field, the field array.
    ######################################################################*/

    int ik, ir;
    double JJ[il+1], JJp[il+1], K_ln;
    
    for(ir=0; ir<Num_R; ir++){
      field[ir] = 0;
    }
    
    for(ik=1; ik<=Kmax; ik++){
      K_ln = Zeros[il][ik]/R_Max;
      for(ir=0; ir<Num_R; ir++){
        Spherical_Bessel(K_ln*Radius[ir], il, JJ, JJp);
        field[ir] += Coef_SB[ik]*JJ[il];
      }
    }
    
    return;
}

/*----------------------------------------------------------------------------*/

extern void SHB_DECOMP(STRUCT_MODEL *Model, STRUCT_ORDER *Order, \
    double **Zeros, enum boundary_condition BC_Flag, \
    complex double ***Coeff_klm){
    
    /*######################################################################
      Purpose:
        spherical Fourier-Bessel decomposition.
      Record of revisions:
        22 Otc. 2019.
      Input parameters:
        Data, a structure saved the data and the corresponding spherical
            coordinate.
        Order, a structure saved the orders of K, L and R, and total number.
        Zeros, the zeros of the spherical Bessel functions or the derivatives
            thereof.
        BC_Flag, boundary condition flag. BC_Flag = enum_zeros for 
            zero-value boundary condition, and =  enum_deri_zeros for 
            zero-derivative boundary condition.
      Output parameters:
        Coeff_klm[k][l][m], the spherical Fourier-Bessel decomposition
            coefficients.
    ######################################################################*/

    int ik, il, im, itheta, iphi, ir;
    double R_Max=Model->R[Order->R_Max-1];
    double JJ[Order->L_Max+2], JJp[Order->L_Max+2], norm=1, K_ln;
    double tmp1, tmp2;

    double ***Coeff = TENSOR_RHO_DBL(Model->nTheta-1, \
        Order->L_Max, false);
    Legendre_Array(Order->L_Max, Model->Theta, Model->nTheta, Coeff);
    
    double *dr = (double *)malloc(Order->R_Max*sizeof(double));
    double *dtheta = (double *)malloc(Model->nTheta*sizeof(double));
    double *dphi = (double *)malloc(Model->nPhi*sizeof(double));

    Cal_Delta(Model->R, Order->R_Max, dr);
    Cal_Delta(Model->Theta, Model->nTheta, dtheta);
    Cal_Delta(Model->Phi, Model->nPhi, dphi);
    
    complex double **MPHI = (complex double **)MATRIX(0, Order->L_Max, 0, \
        Model->nPhi-1, enum_cplx, false);
    for(im=0; im<=Order->L_Max; im++){
      for(iphi=0; iphi<Model->nPhi; iphi++){
        MPHI[im][iphi] = (cos(im*Model->Phi[iphi]) \
            -sin(im*Model->Phi[iphi])*I);
      }
    }
    
    for(ik=1; ik<=Order->K_Max; ik++){
      for(il=0; il<=Order->L_Max; il++){
        for(im=0; im<=il; im++){
          Coeff_klm[ik][il][im]=0;
        }
      }
    }
    
    for(ik=1; ik<=Order->K_Max; ik++){

      for(il=0; il<=Order->L_Max; il++){
        Spherical_Bessel(Zeros[il][ik], il+1, JJ, JJp);
        K_ln = Zeros[il][ik]/R_Max;
        if(BC_Flag == enum_zeros){
          norm = 2/R_Max/R_Max/R_Max/JJ[il+1]/JJ[il+1];
                
        }else{
          if(il == 0){
            if(ik == 1){
              norm = 3./R_Max/R_Max/R_Max;
            }else{
              norm = 2./R_Max/R_Max/R_Max/JJ[il]/JJ[il];
            }
          }else{
            norm = 2./R_Max/R_Max/R_Max/JJ[il]/JJ[il] \
                /(1-il*(il+1)/Zeros[il][ik] \
                /Zeros[il][ik]);
          }
        }
            
        for(ir=0; ir<Order->R_Max; ir++){
          Spherical_Bessel(K_ln*Model->R[ir], il, JJ, JJp);
          tmp1 = Model->R[ir]*Model->R[ir] \
              *dr[ir]*JJ[il]*norm;
          for(itheta=0; itheta<Model->nTheta; itheta++){
            tmp2 = sin(Model->Theta[itheta]);
            for(iphi=0; iphi<Model->nPhi; iphi++){
              tmp2 *=dtheta[itheta]*dphi[iphi];
              for(im=0; im<=il; im++){
                Coeff_klm[ik][il][im] += \
                    Model->Data[ir][itheta][iphi] \
                    *Coeff[itheta][il][im] \
                    *MPHI[im][iphi]*tmp1*tmp2;
              }
            }
          }
        }
      }
    }
    
    for(ik=0; ik<Order->K_Max; ik++){
      for(il=0; il<=Order->L_Max; il++){
        for(im=1; im<=il; im++){
          if(im%2==0){
            Coeff_klm[ik][il][-im] = \
                creal(Coeff_klm[ik][il][im]) \
                -cimag(Coeff_klm[ik][il][im])*I;
          }else{
            Coeff_klm[ik][il][-im] = \
                -creal(Coeff_klm[ik][il][im]) \
                +cimag(Coeff_klm[ik][il][im])*I;
          }
        }
      }
    }
    
    FREE_TENSOR_RHO_DBL(Coeff);
    free(dr);
    free(dtheta);
    free(dphi);
    FREE_MATRIX(MPHI, 0, 0, enum_cplx);

    return;
}















/*----------------------------------------------------------------------------*/

extern void Sparse_Set(STRUCT_ORDER *Order, double Percent, \
    complex double ***Coeff_klm){
    
    /*######################################################################
      Purpose:
         set spherical Fourier-Bessel coefficients to 0 according to the
            sparcity.
      Record of revisions:
        22 Otc. 2019.
      Input parameters:
        Order, a structure saved the orders of K, L and R, and total number.
        Coeff_klm[k][l][m], the spherical Fourier-Bessel coefficients.
        Percent, the percentage of none-zero values.
      Output parameters:
        Coeff_klm[k][l][m], the spherical Fourier-Bessel coefficients.
    ######################################################################*/

    int ik, il, im, ii=1;
    int num = Order->K_Max*(Order->L_Max+1)*(Order->L_Max+1);
    double *array = (double *)VECTOR(1, num, enum_dbl, false);
    double threshold = array[(int)(num*(1-Percent))];

    for(ik=1; ik<=Order->K_Max; ik++){
      for(il=0; il<=Order->L_Max; il++){
        array[ii] = fabs(creal(Coeff_klm[ik][il][0]));
        ii++;
        for(im=1; im<=il; im++){
          array[ii] = fabs(creal(Coeff_klm[ik][il][im]));
          ii++;
          array[ii] = fabs(cimag(Coeff_klm[ik][il][im]));
          ii++;
        }
      }
    }

    HPSORT(num, array);
    
    /*
     for(ii=1; ii<=num; ii++){
     fprintf(stderr, "%e \n",array[ii]);
     }
     */
    
    fprintf(stderr, "\n Sparse Test. Percent = %f  \n",Percent);
    fprintf(stderr, "Number of Parameters = %d \n",num);
    fprintf(stderr, "min value = %e max value= %e \n",array[1],array[num]);
    fprintf(stderr, "threshold = %e \n",threshold);
    
    for(ik=1; ik<=Order->K_Max; ik++){
      for(il=0; il<=Order->L_Max; il++){
        if(fabs(creal(Coeff_klm[ik][il][0]))<=threshold){
          Coeff_klm[ik][il][0] = 0;
        }
        for(im=1; im<=il; im++){
          if(fabs(creal(Coeff_klm[ik][il][im]))<=threshold \
              && fabs(cimag(Coeff_klm[ik][il][im]))<=threshold ){
            Coeff_klm[ik][il][im] = 0;
            Coeff_klm[ik][il][-im] = 0;
          }else if(fabs(creal(Coeff_klm[ik][il][im]))<=threshold){          
            Coeff_klm[ik][il][im] = cimag(Coeff_klm[ik][il][im])*I;
            if(im%2 == 0){
              Coeff_klm[ik][il][-im] = -Coeff_klm[ik][il][im];
            }else{
              Coeff_klm[ik][il][-im] = Coeff_klm[ik][il][im];
            }
                    
          }else if(fabs(cimag(Coeff_klm[ik][il][im]))<=threshold){        
            Coeff_klm[ik][il][im] = creal(Coeff_klm[ik][il][im]);
            if(im%2 == 0){
              Coeff_klm[ik][il][-im] = Coeff_klm[ik][il][im];
            }else{
              Coeff_klm[ik][il][-im] = -Coeff_klm[ik][il][im];
            }
          }
        }
      }
    }
    
    FREE_VECTOR(array, 1, enum_dbl);
    return;
}

/*----------------------------------------------------------------------------*/

static void Legendre_Array(int L_Max, double *Theta, int nTheta, \
    double ***Coef){
    
    /*######################################################################
      Purpose:
        compute associated Legendre functions (with facors) of a Theta 
            array.
      Record of revisions:
        20 Otc. 2019.
      Input parameters:
        L_Max, max order of L.
        Theta, the theta array.
        nTheta, length of the theta array.
      Output parameters:
        Coef[theta][l][m], the coefficient array.
    ######################################################################*/

    int i, itheta, il, im;
    double tmp;
    double *Theta_Cos = (double *)malloc(nTheta*sizeof(double));
    
    for(i=0; i<nTheta; i++){
      Theta_Cos[i] = cos(Theta[i]);
    }
    
    for(il=0; il<=L_Max; il++){
      for(im=0; im<=il; im++){
        tmp = Harmonic_Coefficient(il, im);
        for(itheta=0; itheta<nTheta; itheta++){
          Coef[itheta][il][im]=tmp*Associated_Legendre \
              (il, im, Theta_Cos[itheta]);
        }
      }
    }
    
    for(itheta=0; itheta<nTheta; itheta++){
      for(il=0; il<=L_Max; il++){
        for(im=1; im<=il; im++){
          if(im%2==0){
            Coef[itheta][il][-im]= Coef[itheta][il][im];
          }else{
            Coef[itheta][il][-im]= -Coef[itheta][il][im];
          }
        }
      }
    }
    
    free(Theta_Cos);
    return;
}

/*----------------------------------------------------------------------------*/

static void Cal_Delta(double *Radius, int nR, double *deltaR){
    
    /*######################################################################
      Purpose:
        compute the delta R from a radius array.
      Record of revisions:
        20 Otc. 2019.
      Input parameters:
        Radius[], the radius array.
        nR, the number of the radius array.
      Output parameters:
        deltaR, the delata r array.
    ######################################################################*/

    int ir;
    deltaR[0] = 0.5*(Radius[1]-Radius[0]);
    deltaR[nR-1] = 0.5*(Radius[nR-1]-Radius[nR-2]);
    for(ir=1; ir<nR-1; ir++){
      deltaR[ir] = 0.5*(Radius[ir+1]-Radius[ir-1]);
    }
    return;
}

/*----------------------------------------------------------------------------*/

static void INT_COMP_BC(int il, double **Zeros, double R_Max, \
    double R_Min, int Kmax, int Int_Flag, complex double a, \
    complex double b, enum boundary_condition BC_Flag, \
    complex double *Integ){
    
    /*######################################################################
      Purpose:
        integration of spherical Bessel function (ordr l) from 0 to R_MIN
            (complex values, zero-derivative boundary condition).
      Record of revisions:
        20 Otc. 2019.
      Input parameters:
        il, order L.
        Zeros, the zeros of the derivatives of spherical Bessel functions.
        Radius_MAX, the max value of the radius.
        Radius_MIN, the min value of the radius
        Kmax, max order of K.
        Int_Flag, a flag. Int_Flag = 1: the field is described by a*R_MIN+b
            (constant); Int_Flag = 2: a*r+b (linear); Int_Flag = 3: a*r^2+b
            (parabolic); Int_Flag = 4: a*r*2-a*r+b(parabolic).
        a, b, coefficients used to describe the field between 0 and R_MIN.
      Input parameters:
        Integ, the integration.
    ######################################################################*/

    int ik, ir;
    double JJ[il+1], JJp[il+1], norm, K_ln, Rtmp, dR = R_Min/100;
    double field = a*R_Min+b;
    
    for(ik=1; ik<=Kmax; ik++){
      Integ[ik] = 0;
      Spherical_Bessel(Zeros[il][ik], il+1, JJ, JJp);
      K_ln = Zeros[il][ik]/R_Max;
  
      if(BC_Flag == enum_zeros){
        norm = 2./R_Max/R_Max/R_Max/JJ[il+1]/JJ[il+1];
      }else{
        if(il == 0){
          if(ik == 1){
            norm = 3./R_Max/R_Max/R_Max;
          }else{
            norm = 2./R_Max/R_Max/R_Max/JJ[il]/JJ[il];
          }
        }else{
          norm = 2./R_Max/R_Max/R_Max/JJ[il]/JJ[il]/(1-il*(il+1) \
              /Zeros[il][ik]/Zeros[il][ik]);
        }
      }

      for(ir=1; ir<=100; ir++){
        Rtmp = (ir-0.5)*dR;
        Spherical_Bessel(K_ln*Rtmp, il, JJ, JJp);
        if(Int_Flag == 1){
          Integ[ik] += Rtmp*Rtmp*dR*JJ[il]*norm*field;
        }else if(Int_Flag == 2){
          Integ[ik] += Rtmp*Rtmp*dR*JJ[il]*norm*(a*Rtmp+b);
        }else if(Int_Flag == 3){
          Integ[ik] += Rtmp*Rtmp*dR*JJ[il]*norm*(a*Rtmp*Rtmp+b);
        }else if(Int_Flag == 4){
          Integ[ik] += Rtmp*Rtmp*dR*JJ[il]*norm*(a*Rtmp*Rtmp-a*Rtmp+b);
        }
      }
    }
    return;
}

/*----------------------------------------------------------------------------*/

extern void SH_DECOMP(STRUCT_MODEL *Model, STRUCT_ORDER *Order, \
                      complex double ***Coeff_rlm){
    
    /*######################################################################
      Purpose:
        spherical Harmonical decomposition of each layer of a 3D spherical
            Field.
      Record of revisions:
        22 Otc. 2019.
      Input parameters:
        Data, a structure saved the data and the corresponding spherical
            coordinate.
        Order, a structure saved the orders of K, L and R, and total number.
        Coef_SB[], a spherical Bessel coefficient array.
      Output parameters:
        Coeff_rlm[r][l][m], output the decomposed Spherical Harmonical
            coefficients of each radius r.
    ######################################################################*/

    int ir, il, im, itheta, iphi;
    double tmp;

    double ***Coeff = (double ***)TENSOR_RHO_DBL(Model->nTheta-1, \
        Order->L_Max, false);
    Legendre_Array(Order->L_Max, Model->Theta, Model->nTheta, Coeff);

    double *dtheta = (double *)malloc(Model->nTheta \
        *sizeof(double));

    dtheta[0] = 0.5*Model->Theta[1];
    dtheta[Model->nTheta-1] = 0.5*(Model->Theta[Model->nTheta-1] \
        -Model->Theta[Model->nTheta-2]);
    for(itheta=1; itheta<Model->nTheta-1; itheta++){
      dtheta[itheta] = 0.5*(Model->Theta[itheta+1] \
          -Model->Theta[itheta-1]);
    }
    
    double *dphi = (double *)malloc(Model->nPhi*sizeof(double));

    dphi[0] = 0.5*Model->Phi[1];
    dphi[Model->nPhi-1] = 0.5*(Model->Phi[Model->nPhi-1] \
        -Model->Phi[Model->nPhi-2]);
    for(iphi=1; iphi<Model->nPhi-1; iphi++){
      dphi[iphi] = 0.5*(Model->Phi[iphi+1]-Model->Phi[iphi-1]);
    }
    
    complex double **MPHI = (complex double **)MATRIX(0, Order->L_Max, 0, \
        Model->nPhi-1, enum_cplx, false);
    for(im=0; im<=Order->L_Max; im++){
      for(iphi=0; iphi<Model->nPhi; iphi++){
        MPHI[im][iphi] = (cos(im*Model->Phi[iphi]) \
            -sin(im*Model->Phi[iphi])*I);
      }
    }
    
    for(ir=0; ir<Order->R_Max; ir++){
      for(il=0; il<=Order->L_Max; il++){
        for(im=0; im<=il; im++){
          Coeff_rlm[ir][il][im] = 0;
        }
      }
    }
    
    for(itheta=0; itheta<Model->nTheta; itheta++){
      tmp = sin(Model->Theta[itheta]);
      for(ir=0; ir<Order->R_Max; ir++){
        for(iphi=0; iphi<Model->nPhi; iphi++){
          for(il=0; il<=Order->L_Max; il++){
            for(im=0; im<=il; im++){
              Coeff_rlm[ir][il][im] += Model->Data[ir][itheta][iphi] \
                  *Coeff[itheta][il][im]*dphi[iphi] \
                  *MPHI[im][iphi]*tmp*dtheta[itheta];
            }
          }
        }
      }
    }
    
    for(ir=0; ir<Order->R_Max; ir++){
      for(il=0; il<=Order->L_Max; il++){
        for(im=1; im<=il; im++){
          if(im%2==0){
            Coeff_rlm[ir][il][-im] = creal(Coeff_rlm[ir][il][im]) \
                -cimag(Coeff_rlm[ir][il][im])*I;
          }else{
            Coeff_rlm[ir][il][-im] = -creal(Coeff_rlm[ir][il][im]) \
                +cimag(Coeff_rlm[ir][il][im])*I;
          }
        }
      }
    }
    
    FREE_TENSOR_RHO_DBL(Coeff);
    free(dtheta);
    free(dphi);
    FREE_MATRIX(MPHI, 0, 0, enum_cplx);
    
    return;
}

/*----------------------------------------------------------------------------*/

extern void SB_DECOMP(complex double ***Coeff_rlm, STRUCT_ORDER *Order, \
    double **Zeros, double *Radius, int Int_Flag, \
    enum boundary_condition BC_Flag, complex double ***Coeff_klm){
    
    /*######################################################################
      Purpose:
        spherical Bessel decomposition of Spherical Harmonical coefficients.
      Record of revisions:
        22 Otc. 2019.
      Input parameters:
        Coeff_rlm[r][l][m], the decomposed Spherical Harmonical coefficients
            of each radius r.
        Order, a structure saved the orders of K, L and R, and total number.
        Zeros, the zeros of the spherical Bessel functions or the derivatives
            thereof.
        Radius, a radius array saved the position of data.
        Int_Flag, a flag. Flag = 0: the field between 0 and Radius[0] is 0;
            (zero); Int_Flag = 1: the field is described by a*R_MIN+b
            (constant); Int_Flag = 2: a*r+b (linear); Int_Flag = 3: a*r^2+b
            (parabolic); Int_Flag = 4: a*r*2-a*r+b(parabolic).
        BC_Flag, boundary condition flag. BC_Flag = enum_zeros for 
            zero-value boundary condition, and = enum_deri_zeros for 
            zero-derivative boundary condition.
      Output parameters:
        Coeff_klm[k][l][m], the spherical Fourier-Bessel decomposition
            coefficients.
    ######################################################################*/

    double R_Max = Radius[Order->R_Max-1], R_Min = Radius[0];
    double JJ[Order->L_Max+2], JJp[Order->L_Max+2], norm = 1, K_ln;
    complex double a = 0, b = 0;
    int ik, il, im, ir;
    double tmp;

    complex double *Integ = \
        (complex double *)VECTOR(1, Order->K_Max, enum_cplx, true);

    double *dr = (double *)malloc(Order->R_Max*sizeof(double));

    dr[0] = 0.5*(Radius[1]-Radius[0]);
    dr[Order->R_Max-1] = 0.5*(Radius[Order->R_Max-1] \
        -Radius[Order->R_Max-2]);    
    for(ir=1; ir<Order->R_Max-1; ir++){
      dr[ir] = 0.5*(Radius[ir+1]-Radius[ir-1]);
    }

    if(R_Min > 0 && Int_Flag>=1 && Int_Flag<=4){
      for(il=0; il<=Order->L_Max; il++){
        for(im=0; im<=il; im++){
          if(Int_Flag==1 || Int_Flag==2){
            a = (Coeff_rlm[1][il][im]-Coeff_rlm[0][il][im]) \
                /(Radius[1]-Radius[0]);
            b = Coeff_rlm[0][il][im]-a*Radius[0];
          }else if(Int_Flag==3){
            a = (Coeff_rlm[1][il][im]-Coeff_rlm[0][il][im]) \
                /(Radius[1]-Radius[0])/2;
            b = Coeff_rlm[0][il][im]-a*Radius[0]*Radius[0];
          }else if(Int_Flag==4){
            a = (Coeff_rlm[1][il][im]-Coeff_rlm[0][il][im]) \
                /(Radius[1]-Radius[0]);
            b = Coeff_rlm[0][il][im];
          }

          INT_COMP_BC(il, Zeros, R_Max, R_Min, Order->K_Max, \
              Int_Flag, a, b, BC_Flag, Integ);

          for(ik=1; ik<=Order->K_Max; ik++){
            Coeff_klm[ik][il][im] = Integ[ik];
          }
        }
      }
    }else{
      for(ik=1; ik<=Order->K_Max; ik++){
        for(il=0; il<=Order->L_Max; il++){
          for(im=0; im<=il; im++){
            Coeff_klm[ik][il][im]=0;
          }
        }
      }
    }

    for(ik=1; ik<=Order->K_Max; ik++){
      for(il=0; il<=Order->L_Max; il++){
        Spherical_Bessel(Zeros[il][ik], il+1, JJ, JJp);
        if(BC_Flag == enum_zeros){
          norm = 2/R_Max/R_Max/R_Max/JJ[il+1]/JJ[il+1];
        }else{
          if(il == 0){
            if(ik == 1){
              norm = 3./R_Max/R_Max/R_Max;
            }else{
              norm = 2./R_Max/R_Max/R_Max/JJ[il]/JJ[il];
            }
          }else{
            norm = 2./R_Max/R_Max/R_Max/JJ[il]/JJ[il] \
                /(1-il*(il+1)/Zeros[il][ik]/Zeros[il][ik]);
          }
        }
        K_ln = Zeros[il][ik]/R_Max;
        for(ir=0; ir<Order->R_Max; ir++){
          Spherical_Bessel(K_ln*Radius[ir], il, JJ, JJp);
          tmp = Radius[ir]*Radius[ir]*dr[ir]*JJ[il]*norm;
          for(im=0; im<=il; im++){
            Coeff_klm[ik][il][im] += Coeff_rlm[ir][il][im]*tmp;
          }
        }
      }
    }
    
    for(ik=1; ik<=Order->K_Max; ik++){
      for(il=0; il<=Order->L_Max; il++){
        for(im=1; im<=il; im++){
          if(im%2==0){
            Coeff_klm[ik][il][-im] = creal(Coeff_klm[ik][il][im]) \
                -cimag(Coeff_klm[ik][il][im])*I;
          }else{
            Coeff_klm[ik][il][-im] = -creal(Coeff_klm[ik][il][im]) \
                +cimag(Coeff_klm[ik][il][im])*I;
          }
        }
      }
    }
    
    free(dr);
    FREE_VECTOR(Integ, 1, enum_cplx);
    return;
}

/*----------------------------------------------------------------------------*/

extern void SFBde(STRUCT_MODEL *Model, STRUCT_ORDER *Order, \
    double **Zeros, int Int_Flag, enum boundary_condition BC_Flag, \
    complex double ***Coeff_klm){
    
    /*######################################################################
      Purpose:
        spherical Harmonical decomposition of each layer of a 3D spherical
            Field.
      Record of revisions:
        22 Otc. 2019.
      Input parameters:
        Data, a structure saved the data and the corresponding spherical
            coordinate.
        Order, a structure saved the orders of K, L and R, and total number.
        Coef_SB[], a spherical Bessel coefficient array.
      Output parameters:
        Coeff_rlm[r][l][m], output the decomposed Spherical Harmonical
            coefficients of each radius r.
    ######################################################################*/

    int ik, il, im, ir, itheta, iphi;
    double tmp, norm4p = 0.;
    double R_Max = Model->R[Order->R_Max-1];
    double JJ[Order->L_Max+2], JJp[Order->L_Max+2], K_ln;

    double **norm = (double **)MATRIX(0, Order->L_Max, 1, \
        Order->K_Max, enum_dbl, false);
    double *dtheta = (double *)malloc(Model->nTheta*sizeof(double));
    double *dphi = (double *)malloc(Model->nPhi*sizeof(double));
    double ***Coeff = (double ***)TENSOR_RHO_DBL(Model->nTheta-1, \
        Order->L_Max, false);
    complex double **MPHI = (complex double **)MATRIX(0, Order->L_Max, 0, \
        Model->nPhi-1, enum_cplx, false);
    double *Rintg = (double *)malloc((Order->R_Max+1)*sizeof(double));

    if(Model->Theta[0]>0){
      dtheta[0] = 0.5*(Model->Theta[0]+Model->Theta[1]) \
          *sin(Model->Theta[0]);
    }else{
      dtheta[0] = 0.5*Model->Theta[1]*sin(Model->Theta[0]);
    }

    if(Model->Theta[Model->nTheta-1]<C_Pi){
      dtheta[Model->nTheta-1] = C_Pi-0.5*(Model->Theta[Model->nTheta-1] \
          +Model->Theta[Model->nTheta-2]) \
          *sin(Model->Theta[Model->nTheta-1]);
    }else{
      dtheta[Model->nTheta-1] = 0.5*(C_Pi-Model->Theta[Model->nTheta-2]) \
          *sin(Model->Theta[Model->nTheta-1]);
    }

    for(itheta=1; itheta<Model->nTheta-1; itheta++){
      dtheta[itheta] = 0.5*(Model->Theta[itheta+1] \
          -Model->Theta[itheta-1])*sin(Model->Theta[itheta]);
    }
    
    dphi[0] = C_Pi+0.5*(Model->Phi[1]-Model->Phi[Model->nPhi-1]);
    dphi[Model->nPhi-1] = C_Pi+0.5*(Model->Phi[0] \
        -Model->Phi[Model->nPhi-2]);
    for(iphi=1; iphi<Model->nPhi-1; iphi++){
      dphi[iphi] = 0.5*(Model->Phi[iphi+1]-Model->Phi[iphi-1]);
    }

    for(itheta=0; itheta<Model->nTheta; itheta++){
      for(iphi=0; iphi<Model->nPhi; iphi++){
        norm4p += dtheta[itheta]*dphi[iphi];
      }
    }
    norm4p = C_Pi*4/norm4p;
    
    Legendre_Array(Order->L_Max, Model->Theta, Model->nTheta, Coeff);

    for(im=0; im<=Order->L_Max; im++){
      for(iphi=0; iphi<Model->nPhi; iphi++){
        MPHI[im][iphi] = cos(im*Model->Phi[iphi]) \
            -sin(im*Model->Phi[iphi])*I;
      }
    }

    for(ik=1; ik<=Order->K_Max; ik++){
      for(il=0; il<=Order->L_Max; il++){
        Spherical_Bessel(Zeros[il][ik], il+1, JJ, JJp);
        if(BC_Flag == enum_zeros){
          norm[il][ik] = 2/R_Max/R_Max/R_Max/JJ[il+1]/JJ[il+1];
        }else{
          if(il == 0){
            if(ik == 1){
              norm[il][ik] = 3./R_Max/R_Max/R_Max;
            }else{
              norm[il][ik] = 2./R_Max/R_Max/R_Max/JJ[il]/JJ[il];
            }
          }else{
            norm[il][ik] = 2./R_Max/R_Max/R_Max/JJ[il]/JJ[il] \
                /(1.-il*(il+1.)/Zeros[il][ik]/Zeros[il][ik]);
          }
        }
        norm[il][ik] = sqrt(norm[il][ik]);
      }
    }


    if(Model->R[0] > 0 && Int_Flag>=1 && Int_Flag<=4){
      for(ik=1; ik<=Order->K_Max; ik++){
        for(il=0; il<=Order->L_Max; il++){
          K_ln = Zeros[il][ik]/R_Max;
          tmp = IntgBessel(0, Model->R[0], il, K_ln)*norm[il][ik];
          for(im=0; im<=il; im++){
            for(itheta=0; itheta<Model->nTheta; itheta++){
              for(iphi=0; iphi<Model->nPhi; iphi++){
                Coeff_klm[ik][il][im] +=  dphi[iphi]*dtheta[itheta] \
                    *tmp*Model->Data[0][itheta][iphi] \
                    *Coeff[itheta][il][im]*MPHI[im][iphi];
              }
            }
          }
        }
      }
    }else{
      for(ik=1; ik<=Order->K_Max; ik++){
        for(il=0; il<=Order->L_Max; il++){
          for(im=0; im<=il; im++){
            Coeff_klm[ik][il][im] = 0;
          }
        }
      }
    }

    Rintg[0] = Model->R[0];
    Rintg[Order->R_Max] = Model->R[Order->R_Max-1];

    for(ir=1; ir<Order->R_Max; ir++){
      Rintg[ir] = 0.5*(Model->R[ir-1]+Model->R[ir]);

 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////     
 //     fprintf(stderr,);
    }

    for(ik=1; ik<=Order->K_Max; ik++){
      fprintf(stderr,"%d \n",ik);
      for(il=0; il<=Order->L_Max; il++){
        K_ln = Zeros[il][ik]/R_Max;
        for(ir=0; ir<Order->R_Max; ir++){
          tmp = IntgBessel(Rintg[ir],Rintg[ir+1], il, K_ln)*norm[il][ik];
          for(im=0; im<=il; im++){
            for(itheta=0; itheta<Model->nTheta; itheta++){
              for(iphi=0; iphi<Model->nPhi; iphi++){
                Coeff_klm[ik][il][im] +=  dphi[iphi]*dtheta[itheta] \
                    *Model->Data[ir][itheta][iphi]*MPHI[im][iphi] \
                    *Coeff[itheta][il][im]*tmp;
              }
            }
          }
        }
      }
    }

    FREE_MATRIX(norm, 0, 1, enum_dbl);
    free(dtheta);
    free(dphi);
    free(Rintg);
    FREE_MATRIX(MPHI, 0, 0, enum_cplx);
    FREE_TENSOR_RHO_DBL(Coeff);
    return;              
}



/*----------------------------------------------------------------------------*/

extern void SFBre(STRUCT_ORDER *Order, complex double ***Coeff_rlm, 
    STRUCT_MODEL *Model, double R_Max, double **Zeros, \
    enum boundary_condition BC_Flag){
    
    /*######################################################################
      Purpose:
        reconstruct 3D data from spherical harmonical coefficients.
      Record of revisions:
        22 Otc. 2019
      Input parameters:
        Coord, a structure saved the coordinates.
        Order, a structure saved the orders of K, L and R, and total number.
        Coeff_rlm[r][l][m], the decomposed Spherical Harmonical coefficients
            of each radius r.
      Output parameters:
        Data, output the reconstructed 3D data.
    ######################################################################*/

    int ik, ir, il, im, itheta, iphi;
    double JJ[Order->L_Max+2], JJp[Order->L_Max+2], K_ln;

    double ***Coeff = (double ***)TENSOR_RHO_DBL(Model->nTheta-1, \
        Order->L_Max, false);
    complex double **MPHI = (complex double **)MATRIX(0, Order->L_Max, 0, \
        Model->nPhi-1, enum_cplx, false);
    double **norm = (double **)MATRIX(0, Order->L_Max, 1, \
        Order->K_Max, enum_dbl, false);

    Legendre_Array(Order->L_Max, Model->Theta, Model->nTheta, Coeff);

    for(im=0; im<=Order->L_Max; im++){
      for(iphi=0; iphi<Model->nPhi; iphi++){
        MPHI[im][iphi] = cos(im*Model->Phi[iphi]) \
            +sin(im*Model->Phi[iphi])*I;
      }
    }

    for(ik=1; ik<=Order->K_Max; ik++){
      for(il=0; il<=Order->L_Max; il++){
        Spherical_Bessel(Zeros[il][ik], il+1, JJ, JJp);
        if(BC_Flag == enum_zeros){
          norm[il][ik] = 2/R_Max/R_Max/R_Max/JJ[il+1]/JJ[il+1];
        }else{
          if(il == 0){
            if(ik == 1){
              norm[il][ik] = 3./R_Max/R_Max/R_Max;
            }else{
              norm[il][ik] = 2./R_Max/R_Max/R_Max/JJ[il]/JJ[il];
            }
          }else{
            norm[il][ik] = 2./R_Max/R_Max/R_Max/JJ[il]/JJ[il] \
                /(1.-il*(il+1.)/Zeros[il][ik]/Zeros[il][ik]);
          }
        }
        norm[il][ik] = sqrt(norm[il][ik]);
      }
    }

    for(ir=0; ir<Order->R_Max; ir++){
      for(itheta=0; itheta<Model->nTheta; itheta++){
        for(iphi=0; iphi<Model->nPhi; iphi++){
          Model->Data[ir][itheta][iphi] = 0;
        }
      }
    }

    for(il=0; il<=Order->L_Max; il++){
      for(ik=1; ik<=Order->K_Max; ik++){
        K_ln = Zeros[il][ik]/R_Max;
        for(ir=0; ir<Order->R_Max; ir++){
          Spherical_Bessel(K_ln*Model->R[ir], il, JJ, JJp);
          for(itheta=0; itheta<Model->nTheta; itheta++){
            for(iphi=0; iphi<Model->nPhi; iphi++){
              for(im=0; im<=il; im++){
                if(im>=0){
                  Model->Data[ir][itheta][iphi] += norm[il][ik]*JJ[il]
                      *creal(Coeff_rlm[ik][il][im]*MPHI[im][iphi]) \
                      *Coeff[itheta][il][im];
                }else{
                  Model->Data[ir][itheta][iphi] += norm[il][ik]*JJ[il]
                      *creal(Coeff_rlm[ik][il][im]*MPHI[im][iphi]) \
                      *Coeff[itheta][il][im]*2;
               
                }
              }
            }
          }
        }
      }
    }

    
    FREE_TENSOR_RHO_DBL(Coeff);

    FREE_MATRIX(MPHI, 0, 0, enum_cplx);
    FREE_MATRIX(norm, 0, 1, enum_dbl);

    return;
}

/*----------------------------------------------------------------------------*/



extern void SHB(STRUCT_MODEL *Model, STRUCT_ORDER *Order, \
                      complex double ***Coeff_rlm){
    
    /*######################################################################
      Purpose:
        spherical Harmonical decomposition of each layer of a 3D spherical
            Field.
      Record of revisions:
        22 Otc. 2019.
      Input parameters:
        Data, a structure saved the data and the corresponding spherical
            coordinate.
        Order, a structure saved the orders of K, L and R, and total number.
        Coef_SB[], a spherical Bessel coefficient array.
      Output parameters:
        Coeff_rlm[r][l][m], output the decomposed Spherical Harmonical
            coefficients of each radius r.
    ######################################################################*/

    int ir, il, im, itheta, iphi;
    double tmp;
    double norm4p = 0.0;

    double *dtheta = (double *)malloc(Model->nTheta*sizeof(double));
    double *dphi = (double *)malloc(Model->nPhi*sizeof(double));
    double ***Coeff = (double ***)TENSOR_RHO_DBL(Model->nTheta-1, \
        Order->L_Max, false);
    complex double **MPHI = (complex double **)MATRIX(0, Order->L_Max, 0, \
        Model->nPhi-1, enum_cplx, false);

    if(Model->Theta[0]>0){
      dtheta[0] = 0.5*(Model->Theta[0]+Model->Theta[1]) \
          *sin(Model->Theta[0]);
    }else{
      dtheta[0] = 0.5*Model->Theta[1]*sin(Model->Theta[0]);
    }

    if(Model->Theta[Model->nTheta-1]<C_Pi){
      dtheta[Model->nTheta-1] = C_Pi-0.5*(Model->Theta[Model->nTheta-1] \
          +Model->Theta[Model->nTheta-2]) \
          *sin(Model->Theta[Model->nTheta-1]);
    }else{
      dtheta[Model->nTheta-1] = 0.5*(C_Pi-Model->Theta[Model->nTheta-2]) \
          *sin(Model->Theta[Model->nTheta-1]);
    }

    for(itheta=1; itheta<Model->nTheta-1; itheta++){
      dtheta[itheta] = 0.5*(Model->Theta[itheta+1] \
          -Model->Theta[itheta-1])*sin(Model->Theta[itheta]);
    }
    
    dphi[0] = C_Pi+0.5*(Model->Phi[1]-Model->Phi[Model->nPhi-1]);
    dphi[Model->nPhi-1] = C_Pi+0.5*(Model->Phi[0] \
        -Model->Phi[Model->nPhi-2]);
    for(iphi=1; iphi<Model->nPhi-1; iphi++){
      dphi[iphi] = 0.5*(Model->Phi[iphi+1]-Model->Phi[iphi-1]);
    }

    for(itheta=0; itheta<Model->nTheta; itheta++){
      for(iphi=0; iphi<Model->nPhi; iphi++){
        norm4p += dtheta[itheta]*dphi[iphi];
      }
    }
    norm4p = C_Pi*4/norm4p;
    
    Legendre_Array(Order->L_Max, Model->Theta, Model->nTheta, Coeff);

    for(im=0; im<=Order->L_Max; im++){
      for(iphi=0; iphi<Model->nPhi; iphi++){
        MPHI[im][iphi] = (cos(im*Model->Phi[iphi]) \
            -sin(im*Model->Phi[iphi])*I);
      }
    }
    


    for(ir=0; ir<Order->R_Max; ir++){
      for(il=0; il<=Order->L_Max; il++){
        for(im=0; im<=il; im++){
          Coeff_rlm[ir][il][im] = 0;
        }
      }
    }
    
    for(itheta=0; itheta<Model->nTheta; itheta++){
      tmp = sin(Model->Theta[itheta]);
      for(ir=0; ir<Order->R_Max; ir++){
        for(iphi=0; iphi<Model->nPhi; iphi++){
          for(il=0; il<=Order->L_Max; il++){
            for(im=0; im<=il; im++){
              Coeff_rlm[ir][il][im] += Model->Data[ir][itheta][iphi] \
                  *Coeff[itheta][il][im]*dphi[iphi] \
                  *MPHI[im][iphi]*tmp*dtheta[itheta];
            }
          }
        }
      }
    }
    
    for(ir=0; ir<Order->R_Max; ir++){
      for(il=0; il<=Order->L_Max; il++){
        for(im=1; im<=il; im++){
          if(im%2==0){
            Coeff_rlm[ir][il][-im] = creal(Coeff_rlm[ir][il][im]) \
                -cimag(Coeff_rlm[ir][il][im])*I;
          }else{
            Coeff_rlm[ir][il][-im] = -creal(Coeff_rlm[ir][il][im]) \
                +cimag(Coeff_rlm[ir][il][im])*I;
          }
        }
      }
    }
    
    FREE_TENSOR_RHO_DBL(Coeff);
    free(dtheta);
    free(dphi);
    FREE_MATRIX(MPHI, 0, 0, enum_cplx);
    
    return;
}

/*----------------------------------------------------------------------------*/

extern void SFB_DECOMP(STRUCT_MODEL *Model, STRUCT_ORDER *Order, \
    double **Zeros, int Int_Flag, enum boundary_condition BC_Flag, \
    complex double ***Coeff_klm){
    
    /*######################################################################
      Purpose:
        spherical Fourier-Bessel decomposition.
      Record of revisions:
        22 Otc. 2019.
      Input parameters:
        Data, a structure saved the data and the corresponding spherical
            coordinate.
        Order, a structure saved the orders of K, L and R, and total number.
        Zeros, the zeros of the spherical Bessel functions or the derivatives
            thereof.
        Int_Flag, a flag. Flag = 0: the field between 0 and Radius[0] is 0;
            (zero); Int_Flag = 1: the field is described by a*R_MIN+b
            (constant); Int_Flag = 2: a*r+b (linear); Int_Flag = 3: a*r^2+b
            (parabolic); Int_Flag = 4: a*r*2-a*r+b(parabolic).
        BC_Flag, boundary condition flag. BC_Flag = 1 for zero-value boundary
            condition, and = 2 for zero-derivative boundary condition.
      Output parameters:
        Coeff_klm[k][l][m], the spherical Fourier-Bessel decomposition
            coefficients.
    ######################################################################*/

    complex double ***Coeff_rlm = TENSOR_RHO_CPLX(Order->R_Max-1, \
        Order->L_Max, true);
    
    SH_DECOMP(Model, Order, Coeff_rlm);
    
    SB_DECOMP(Coeff_rlm, Order, Zeros, Model->R, Int_Flag, BC_Flag, \
        Coeff_klm);
    
    FREE_TENSOR_RHO_CPLX(Coeff_rlm);
    
    return;
    
}

/*----------------------------------------------------------------------------*/

extern void SB_RECONSTRUCT(STRUCT_ORDER *Order, double **Zeros, \
    double *Radius, complex double ***Coeff_klm, \
    complex double ***Coeff_rlm){
    
    /*######################################################################
      Purpose:
        reconstruct Spherical Harmonical coefficients from Spherical Bessel
            coefficients.
      Record of revisions:
        22 Otc. 2019.
      Input parameters:
        Order, a structure saved the orders of K, L and R, and total number.
        Zeros, the zeros of the spherical Bessel functions or the 
            derivatives thereof.
        Radius, a radius array saved the position of data.
        Coeff_klm[k][l][m], the spherical Fourier-Bessel decomposition
            coefficients.
      Output parameters:
        Coeff_rlm[r][l][m], the Spherical Harmonical coefficients
    ######################################################################*/

    int ik, ir, il, im;
    double JJ[Order->L_Max+1], JJp[Order->L_Max+1], K_ln;
    double Radius_MAX = Radius[Order->R_Max-1];
    
    for(il=0; il<=Order->L_Max; il++){
      for(ik=1; ik<=Order->K_Max; ik++){
        K_ln = Zeros[il][ik]/Radius_MAX;
        for(ir=0; ir<Order->R_Max; ir++){
          Spherical_Bessel(K_ln*Radius[ir], il, JJ, JJp);
          for(im=-il; im<=il; im++){
            Coeff_rlm[ir][il][im] += Coeff_klm[ik][il][im]*JJ[il];
          }
        }
      }
    }
    
    return;
}

/*----------------------------------------------------------------------------*/

extern void SH_RECONSTRUCT(STRUCT_ORDER *Order, \
    complex double ***Coeff_rlm, STRUCT_MODEL *Model){
    
    /*######################################################################
      Purpose:
        reconstruct 3D data from spherical harmonical coefficients.
      Record of revisions:
        22 Otc. 2019
      Input parameters:
        Coord, a structure saved the coordinates.
        Order, a structure saved the orders of K, L and R, and total number.
        Coeff_rlm[r][l][m], the decomposed Spherical Harmonical coefficients
            of each radius r.
      Output parameters:
        Data, output the reconstructed 3D data.
    ######################################################################*/

    int ir, il, im, itheta, iphi;
    
    double ***Coeff = (double ***)TENSOR_RHO_DBL(Model->nTheta-1, \
        Order->L_Max, false);
    Legendre_Array(Order->L_Max, Model->Theta, Model->nTheta, Coeff);

    double *dtheta = (double *)malloc(Model->nTheta*sizeof(double));

    dtheta[0] = 0.5*Model->Theta[1];
    dtheta[Model->nTheta-1] = 0.5*(Model->Theta[Model->nTheta-1] \
        -Model->Theta[Model->nTheta-2]);
    for(itheta=1; itheta<Model->nTheta-1; itheta++){
      dtheta[itheta] = 0.5*(Model->Theta[itheta+1] \
          -Model->Theta[itheta-1]);
    }
    
    double *dphi = (double *)malloc(Model->nPhi*sizeof(double));
    dphi[0] = 0.5*Model->Phi[1];
    dphi[Model->nPhi-1] = \
        0.5*(Model->Phi[Model->nPhi-1]-Model->Phi[Model->nPhi-2]);
    for(iphi=1; iphi<Model->nPhi-1; iphi++){
      dphi[iphi] = 0.5*(Model->Phi[iphi+1]-Model->Phi[iphi-1]);
    }
    
    double **MPHI_R = (double **)MATRIX(0, Order->L_Max, 0, Model->nPhi-1, \
        enum_dbl, false);
    double **MPHI_I = (double **)MATRIX(0, Order->L_Max, 0, Model->nPhi-1, \
        enum_dbl, false);
    
    for(im=0; im<=Order->L_Max; im++){
      for(iphi=0; iphi<Model->nPhi; iphi++){
        MPHI_R[im][iphi]=cos(im*Model->Phi[iphi]);
        MPHI_I[im][iphi]=sin(im*Model->Phi[iphi]);
      }
    }
        
    for(ir=0; ir<Order->R_Max; ir++){
      for(itheta=0; itheta<Model->nTheta; itheta++){
        for(iphi=0; iphi<Model->nPhi; iphi++){
          Model->Data[ir][itheta][iphi] = 0;
                
          for(il=0; il<=Order->L_Max; il++){
            for(im=-il; im<=il; im++){
              if(im>=0){
                Model->Data[ir][itheta][iphi] += creal(Coeff_rlm[ir][il][im] \
                    *Coeff[itheta][il][im] \
                    *(MPHI_R[im][iphi]+MPHI_I[im][iphi]*I));
              }else if((-im)%2==0){
                Model->Data[ir][itheta][iphi] += creal(Coeff_rlm[ir][il][im] \
                    *Coeff[itheta][il][-im] \
                    *(MPHI_R[-im][iphi]-MPHI_I[-im][iphi]*I));
              }else{
                Model->Data[ir][itheta][iphi] += creal(Coeff_rlm[ir][il][im] \
                    *Coeff[itheta][il][-im] \
                    *(-MPHI_R[-im][iphi]+MPHI_I[-im][iphi]*I));
              }
            }
          }
        }
      }
    }
    
    free(dtheta);
    free(dphi);
    FREE_MATRIX(MPHI_R, 0, 0, enum_dbl);
    FREE_MATRIX(MPHI_I, 0, 0, enum_dbl);
    FREE_TENSOR_RHO_DBL(Coeff);

    return;
}

/*----------------------------------------------------------------------------*/

extern void SFB_RECONSTRUCT(complex double ***Coeff_klm, \
    STRUCT_ORDER *Order, double **Zeros, STRUCT_MODEL *Model){
    
    /*######################################################################
      Purpose:
        spherical Fourier-Bessel reconstruction.
      Record of revisions:
        22 Otc. 2019.
      Input parameters:
        Coeff_klm[k][l][m], the spherical Fourier-Bessel coefficients.
        Coord, a structure saved the coordinates.
        Order, a structure saved the orders of K, L and R, and total number.
        Zeros, the zeros of the spherical Bessel functions or the derivatives
            thereof.
      Output parameters:
        Reconstruction, the reconstructed data.
    ######################################################################*/

    complex double ***Coeff_rlm = TENSOR_RHO_CPLX(Order->R_Max-1, \
        Order->L_Max, true);
    
    SB_RECONSTRUCT(Order, Zeros, Model->R, Coeff_klm, Coeff_rlm);
    
    SH_RECONSTRUCT(Order, Coeff_rlm, Model);
    
    FREE_TENSOR_RHO_CPLX(Coeff_rlm);
    
    return;
    
}


extern double Xvalue(double x, int N){

    /*######################################################################
      Purpose:
        Compute Xvalue in Bloomfield 2017 'Indefinite Integrals of 
            Spherical Bessel Functions'.
      Record of revisions:
        22 Otc. 2019.
      Input parameters:

      Output parameters:

    ######################################################################*/  

    double X0, Y0, tmpx, tmpy, *xx;
    int in;
    double Sr = sin(x);
    double Cs = cos(x);

    if(N>=0){
      X0 = -Cs;
//###########################      
      if(in==0) return X0;
//###########################      
      Y0 = Sr;
      xx = (double *)malloc((N+1)*sizeof(double));
      xx[1] = x;
      for(in=2;in<=N;in++){
        xx[in] = xx[in-1]*x;
      }

      for(in=1;in<=N;in++){
        tmpx = in*Y0-xx[in]*Cs;
        tmpy = xx[in]*Sr-in*X0;
        X0 = tmpx;
        Y0 = tmpy; 
      }

      free(xx);
      return X0;

    }else{

      cisi(x, &Y0, &X0);
      if(in==-1) return X0;

      xx = (double *)VECTOR(N,-1,enum_dbl,false);

      xx[-1] = 1./x;
      for(in=-2;in>=N;in--){
        xx[in] = xx[in+1]/x;
      }
      for(in=-2;in>=N;in--){
        tmpx = (xx[in+1]*Sr-Y0)/(in+1.);
        tmpy = (xx[in+1]*Cs+X0)/(in+1.);
        X0 = tmpx;
        Y0 = tmpy;
      }
      FREE_VECTOR(xx,N,enum_dbl);
      return X0;
    } 

    return X0; 

}

extern double Inl(double x, int n, int l, double *J){


    double I0nl = Xvalue(x, n-l-1);
    double tmp;
    int il, np;

    double *xx = (double *)VECTOR(n-l-1,n,enum_dbl,false);

    xx[0] = 1.;
    for(il=1;il<=n;il++){
      xx[il] = xx[il-1]*x;
    }
    for(il=-1;il>=n-l-1;il--){
      xx[il] = xx[il+1]/x;
    }

    for(il=1,np=n-l+1;il<=l;il++,np++){
      tmp = (il+np-1)*I0nl-xx[np]*J[il-1];
 //     fprintf(stderr,"%d %d %e %e %e %e\n",il,np,I0nl,tmp,xx[np],pow(x,np));
      I0nl = tmp;
    }

    FREE_VECTOR(xx,n-l-1,enum_dbl);


    return I0nl;

}

extern double IntgBessel(double x0, double x1, int l, double knl){



    double delta0 = 0.01;
    int num = (int)((x1-x0)/delta0);
    double delta = (x1-x0)/num;
    int i;
    double r;
    double sum =0;
    double JJ[l], JJp[l];
    for(i=0;i<num;i++){
      r = x0+(i+0.5)*delta;
      Spherical_Bessel(r*knl, l, JJ, JJp);
      sum += JJ[l]*r*r*delta;

    }



    return sum;
}

/*----------------------------------------------------------------------------*/
