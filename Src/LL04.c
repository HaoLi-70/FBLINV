
#include "LL04.h"

/*----------------------------------------------------------------------------*/

    /*######################################################################
     
      revision log:
      
        15 Aug. 2023.
     
    ######################################################################*/

/*----------------------------------------------------------------------------*/

extern double Geffect(double Gu, double Gl, double Ju, double Jl){
    
    /*######################################################################
      Purpose:
        computes the effective Lande factor.
      Record of revisions:
        25 Apr. 2018.
      Input parameters:
        Gu, The lande factor of upper level.
        Gl, The lande factor of lower level.
        Ju, The total angular momentum of upper level.
        Jl, The total angular momentum of lower level.
      Return:
        return the effect Lande factor.
      References:
        LL04 Chapter 3, Equation 3.44.
    ######################################################################*/
    
    return 0.5*(Gu+Gl)+0.25*(Gu-Gl)*(Ju*(Ju+1)-Jl*(Jl+1));

}

/*----------------------------------------------------------------------------*/

extern double Gfactor(double J, double L, double S){
    
    /*######################################################################
      Purpose:
        computes the Lande factor.
      Record of revisions:
        25 Apr. 2018.
      Input parameters:
        J, The total angular momentum of the electronic cloud.
        L, The total orbital angular momentum of the electronic cloud.
        S, The total spin of the electronic cloud.
      Return:
        return the Lande factor.
      Method:
        L-S coupling asumption, and if the number is 0, return 0 for
            convernience.
    ######################################################################*/

    if(J==0){ 
      return 0;
    }else{  
      return 1.0+(J*(J+1.0)-L*(L+1.0)+S*(S+1.0))/(2.0*J*(J+1.0));
    }

}

/*----------------------------------------------------------------------------*/

extern complex double Djmn(double Alpha, double Beta, double Gamma, \
    double J, double M, double N, STR_FCTSG *fctsg){
    
    /*######################################################################
      Purpose:
        compute the element of the reduced rotation matrices 
            (D^J_MN(Alpha,Beta,Gamma)).
      Record of revisions:
        13 Aug. 2023.
      Input parameters:
        Alpha, Beta, Gamma, the angles.
        J, M, N, quantum numbers.
        fctsg, a structure with factorial of integer.
      Return:
        return the element of the reduced rotation matrices .
      References:
        LL04 Chapter 2, page 54 equation 2.68 and page 56 equation 2.69,
            rotation through an angle Alpha, Beta, and Gamma about the
            Z-axis, Y-axis, and Z-axis.
            page 52 Fig. 2.4
    ######################################################################*/

    double a, b = 0;
    double complex Djmn;
    int t, t1, t2;
    double Sb2 = sin(Beta/2), Cb2 = cos(Beta/2);

    a = sqrt(fctsg->fct[(int)(J+M)]*fctsg->fct[(int)(J-M)] \
        *fctsg->fct[(int)(J+N)]*fctsg->fct[(int)(J-N)]);
    
    t1 = (int)(M-N)<=0?0:(int)(M-N);
    
    t2 = (int)(J+M)<=(int)(J-N)?(int)(J+M):(int)(J-N);
    
    for(t=t1; t<=t2; t++){
      b = b+fctsg->sg[t]*pow(Cb2,(2*J+M-N-2*t)) \
          *pow(Sb2,(2*t-M+N))/fctsg->fct[(int)(J+M-t)] \
          /fctsg->fct[(int)(J-N-t)]/fctsg->fct[t]/fctsg->fct[(int)(t+N-M)];
    }
    
    Djmn = a*b*(cos(M*Alpha+N*Gamma)-sin(M*Alpha+N*Gamma)*I);
    
    return Djmn;
}

/*----------------------------------------------------------------------------*/

extern int Rotmat(double Alpha, double Beta, double Gamma, \
    complex double ***Djmn, int Rank){

    double Sr = sin(Beta);
    double Cr = cos(Beta);
    int j, m, n;

    switch(Rank){

      case 0:
        Djmn[0][0][0] = 1;
        break;

      case 1:

        Djmn[1][-1][-1] = 0.5*(1.0+Cr);
        Djmn[1][-1][0] = Sr/C_sqrt2;
        Djmn[1][-1][1] = 0.5*(1-Cr);

        Djmn[1][0][-1] = -Sr/C_sqrt2;
        Djmn[1][0][0] = Cr;
        Djmn[1][0][1] = -Djmn[1][0][-1];

        Djmn[1][1][-1] = Djmn[1][-1][1];
        Djmn[1][1][0] = -Djmn[1][-1][0];
        Djmn[1][1][1] = Djmn[1][-1][-1];

        for(m=-1;m<=1;m++){
          for(n=-1;n<=1;n++){
            Djmn[1][m][n] = Djmn[1][m][n]\
                *(cos(Alpha*m+Gamma*n)+sin(Alpha*m+Gamma*n)*I);
          }
        }
        break;

      case 2:

        Djmn[2][-2][-2] = 0.25*(1+Cr)*(1+Cr);
        Djmn[2][-2][-1] = 0.5*Sr*(1+Cr);
        Djmn[2][-2][0] = 0.5*C_sqrt3/C_sqrt2*Sr*Sr;
        Djmn[2][-2][1] = 0.5*Sr*(1-Cr);
        Djmn[2][-2][2] = 0.25*(1-Cr)*(1-Cr);

        Djmn[2][-1][-2] = -Djmn[2][-2][-1];
        Djmn[2][-1][-1] = (Cr-0.5)*(Cr+1.0);
        Djmn[2][-1][0] = C_sqrt3/C_sqrt2*Sr*Cr;
        Djmn[2][-1][1] = (Cr+0.5)*(1-Cr);
        Djmn[2][-1][2] = Djmn[2][-2][1];

        Djmn[2][0][-2] = Djmn[2][-2][0];
        Djmn[2][0][-1] = -Djmn[2][-1][0];
        Djmn[2][0][0] = 0.5*(3.0*Cr*Cr-1);
        Djmn[2][0][1] = Djmn[2][-1][0];
        Djmn[2][0][2] = Djmn[2][-2][0];

        Djmn[2][1][-2] = -Djmn[2][-2][1];
        Djmn[2][1][-1] = Djmn[2][-1][1];
        Djmn[2][1][0] = Djmn[2][0][-1];
        Djmn[2][1][1] = Djmn[2][-1][-1];
        Djmn[2][1][2] = Djmn[2][-2][-1];

        Djmn[2][2][-2] = Djmn[2][-2][2];
        Djmn[2][2][-1] =  Djmn[2][1][-2];
        Djmn[2][2][0] = Djmn[2][-2][0];
        Djmn[2][2][1] = Djmn[2][-1][-2];
        Djmn[2][2][2] = Djmn[2][-2][-2];

        for(m=-2;m<=2;m++){
          for(n=-2;n<=2;n++){
            Djmn[2][m][n] = Djmn[2][m][n]\
                *(cos(Alpha*m+Gamma*n)-sin(Alpha*m+Gamma*n)*I);
          }
        }
        break;

      default:

        Djmn[0][0][0] = 1;

        Djmn[1][-1][-1] = 0.5*(1.0+Cr);
        Djmn[1][-1][0] = Sr/C_sqrt2;
        Djmn[1][-1][1] = 0.5*(1-Cr);

        Djmn[1][0][-1] = -Djmn[1][-1][0];
        Djmn[1][0][0] = Cr;
        Djmn[1][0][1] = Djmn[1][-1][0];

        Djmn[1][1][-1] = Djmn[1][-1][1];
        Djmn[1][1][0] = Djmn[1][0][-1];
        Djmn[1][1][1] = Djmn[1][-1][-1];


        Djmn[2][-2][-2] = 0.25*(1+Cr)*(1+Cr);
        Djmn[2][-2][-1] = 0.5*Sr*(1+Cr);
        Djmn[2][-2][0] = 0.5*C_sqrt3/C_sqrt2*Sr*Sr;
        Djmn[2][-2][1] = 0.5*Sr*(1-Cr);
        Djmn[2][-2][2] = 0.25*(1-Cr)*(1-Cr);

        Djmn[2][-1][-2] = -Djmn[2][-2][-1];
        Djmn[2][-1][-1] = (Cr-0.5)*(Cr+1.0);
        Djmn[2][-1][0] = C_sqrt3/C_sqrt2*Sr*Cr;
        Djmn[2][-1][1] = (Cr+0.5)*(1-Cr);
        Djmn[2][-1][2] = Djmn[2][-2][1];

        Djmn[2][0][-2] = Djmn[2][-2][0];
        Djmn[2][0][-1] = -Djmn[2][-1][0];
        Djmn[2][0][0] = 0.5*(3.0*Cr*Cr-1);
        Djmn[2][0][1] = Djmn[2][-1][0];
        Djmn[2][0][2] = Djmn[2][-2][0];

        Djmn[2][1][-2] = -Djmn[2][-2][1];
        Djmn[2][1][-1] = Djmn[2][-1][1];
        Djmn[2][1][0] = Djmn[2][0][-1];
        Djmn[2][1][1] = Djmn[2][-1][-1];
        Djmn[2][1][2] = Djmn[2][-2][-1];

        Djmn[2][2][-2] = Djmn[2][-2][2];
        Djmn[2][2][-1] =  Djmn[2][1][-2];
        Djmn[2][2][0] = Djmn[2][-2][0];
        Djmn[2][2][1] = Djmn[2][-1][-2];
        Djmn[2][2][2] = Djmn[2][-2][-2];

        for(j=1;j<=2;j++){
          for(m=-j;m<=j;m++){
            for(n=-j;n<=j;n++){
              Djmn[j][m][n] *= (cos(Alpha*m+Gamma*n)-sin(Alpha*m+Gamma*n)*I);
            }
          }
        }

        break;
    }

    return 0;
}

/*----------------------------------------------------------------------------*/

extern int TKP (complex double ***T){
    
    /*######################################################################
      Purpose:
        geometry tensor TKP.
      Record of revisions:
        14 Aug. 2023.
      Output parameters:
        T[i][k][p], the tensors. 0<=i<=3,0<=k<=2,-k<=p<=k
      References:
        LL04 Chapter 5, page 210 table 5.5.
    ######################################################################*/

    long Byte = 36*sizeof(complex double);
    memset(T[0][0],0,Byte);

    T[0][0][0] = 1;
    T[0][2][0] = 1/C_sqrt2;

    T[1][2][-2] = -C_sqrt3/2;
    T[1][2][2] = -C_sqrt3/2;

    T[2][2][-2] = C_sqrt3/2*I;
    T[2][2][2] = -C_sqrt3/2*I;

    T[3][1][0] = C_sqrt3/C_sqrt2;

    return 0;

}

/*----------------------------------------------------------------------------*/

extern void TKQ(complex double ***T, double Chi, double Theta, \
    double Gamma){
    
    /*######################################################################
      Purpose:
        geometry tensor TKQ.
      Record of revisions:
        1 Dec. 2019.
      Input parameters:
        Theta, Chi, Gamma, the angles in radians.
        geometry correspond to fig 5.14 at page 196.
        it is actually a rotation of TKP by Djmn(-gamma, -theta, -chi)
      Output parameters:
        T[i][k][p], the tensors.
      References:
        LL04 Chapter 5, page 211 table 5.6.
    ######################################################################*/

    // i=0
    double St = sin(Theta);
    double Ct = cos(Theta);
    double Sc = sin(Chi);
    double Cc = cos(Chi);
    double S2c = sin(2.*Chi);
    double C2c = cos(2.*Chi);
    double S2g = sin(2.*Gamma);
    double C2g = cos(2.*Gamma);

    // i=0
    T[0][0][0] = 1.;
    T[0][1][-1] = 0;
    T[0][1][0] = 0;
    T[0][1][1] = 0;
    T[0][2][-2] = C_sqrt3/4.*St*St*(C2c-S2c*I);
    T[0][2][-1] = C_sqrt3/2.*St*Ct*(Cc-Sc*I);
    T[0][2][0] = 0.5/C_sqrt2*(3*Ct*Ct-1);
    T[0][2][1] = -C_sqrt3/2.*St*Ct*(Cc+Sc*I);
    T[0][2][2] = C_sqrt3/4.*St*St*(C2c+S2c*I);

    // i=1
    T[1][0][0] = 0;
    T[1][1][-1] = 0;
    T[1][1][0] = 0;
    T[1][1][1] = 0;
    T[1][2][-2] = -C_sqrt3/4.*(C2g*(1+Ct*Ct)-2.*S2g*Ct*I)*(C2c-S2c*I);
    T[1][2][-1] = C_sqrt3/2.*St*(C2g*Ct-S2g*I)*(Cc-Sc*I);
    T[1][2][0] = -1.5/C_sqrt2*C2g*St*St;
    T[1][2][1] = -C_sqrt3/2.*St*(C2g*Ct+S2g*I)*(Cc+Sc*I);
    T[1][2][2] = -C_sqrt3/4.*(C2g*(1+Ct*Ct)+2.*S2g*Ct*I)*(C2c+S2c*I);

    // i=2
    T[2][0][0] = 0;
    T[2][1][-1] = 0;
    T[2][1][0] = 0;
    T[2][1][1] = 0;
    T[2][2][-2] = C_sqrt3/4.*(S2g*(1+Ct*Ct)+2.*C2g*Ct*I)*(C2c-S2c*I);
    T[2][2][-1] = -C_sqrt3/2.*St*(S2g*Ct+C2g*I)*(Cc-Sc*I);
    T[2][2][0] = 1.5/C_sqrt2*S2g*St*St;
    T[2][2][1] = C_sqrt3/2.*St*(S2g*Ct-C2g*I)*(Cc+Sc*I);
    T[2][2][2] = C_sqrt3/4.*(S2g*(1+Ct*Ct)-2.*C2g*Ct*I)*(C2c+S2c*I);

    // i=3
    T[3][0][0] = 0;
    T[3][1][-1] = C_sqrt3/2.*St*(Cc-Sc*I);
    T[3][1][0] = C_sqrt3/C_sqrt2*Ct;
    T[3][1][1] = -C_sqrt3/2.*St*(Cc+Sc*I);
    T[3][2][-2] = 0;
    T[3][2][-1] = 0;
    T[3][2][0] = 0;
    T[3][2][1] = 0;
    T[3][2][2] = 0;
    
    return;
}

/*----------------------------------------------------------------------------*/

extern void TKQ90(complex double ***T){
    
    /*######################################################################
      Purpose:
        geometry tensor TKQ for forbidden lines 
            with Chi=0, double Theta=Pi/2, Gamma=0
      Record of revisions:
        14 Aug. 2023.
      Input parameters:
        .
      Output parameters:
        T[i][k][p], the tensors.
      References:
        LL04 Chapter 5, page 211 table 5.6.
    ######################################################################*/

    long Byte = 27*sizeof(complex double);
    memset(T[0][0],0,Byte);

    // i=0
    T[0][0][0] = 1.;
    T[0][2][-2] = C_sqrt3/4.;
    T[0][2][0] = -0.5/C_sqrt2;
    T[0][2][2] = C_sqrt3/4.;

    // i=1
    T[1][2][-2] = -C_sqrt3/4.;
    T[1][2][0] = -1.5/C_sqrt2;
    T[1][2][2] = -C_sqrt3/4.;

    // i=2
    T[2][2][-1] = -C_sqrt3/2.*I;
    T[2][2][1] = -C_sqrt3/2.*I;
    
    return;
}
/*----------------------------------------------------------------------------*/

extern double Planck(double nu, double T){
    
    /*######################################################################
      Purpose:
        compute the Planck function.
      Record of revisions:
        30 Nov. 2019.
      Input parameters:
        nu, the frequency.
        T, the temperature.
      Return:
        return the intensity at frequency nu for temperature T.
    ######################################################################*/

    return (2*C_h*nu*nu*nu/C_c/C_c)/(exp(C_h*nu/C_Kb/T)-1);

}

/*----------------------------------------------------------------------------*/

extern int Jkq_off(double u1, double u2, double r, double *Jkq){
    
    /*######################################################################
      Purpose:
        compute incident radiation tensor J00.
      Record of revisions:
        1 Dec. 2019.
      Input parameters:
        u1, u2, the limb darkening coefficients
        r, distance from solar center, in the unit of the solar radius.
      Output parameters:
        Jkq, J^0_0, and J^2_0;
      References:
        LL04 Chapter 12, page 675 Eqs. 12.34-12.37.
    ######################################################################*/

    double Sr = 1/r;
    double Cr = sqrt(1-Sr*Sr);

    if(u1==0 && u2==0){
      Jkq[0] = 1/2.*(1-Cr);
      Jkq[1] = 1/4./C_sqrt2*Cr*Sr*Sr;
    }else{

      double a[3], c[3];
      
      a[0] = 1-Cr;
      a[1] = Cr-0.5-0.5*Cr*Cr/Sr*log((1+Sr)/Cr);
      a[2] = (Cr+2)*(Cr-1)/3/(Cr+1);
      
      c[0] = Cr*Sr*Sr;
      c[1] = 0.125*(8*Cr*Cr*Cr-3*Cr*Cr-8*Cr+2 \
          +(4-3*Cr*Cr)*Cr*Cr/Sr*log((1+Sr)/Cr));
      c[2] = (Cr-1)/15./(Cr+1)*(9*Cr*Cr*Cr+18*Cr*Cr+7*Cr-4);
      
      Jkq[0] = 0.5*(a[0]+a[1]*u1+a[2]*u2);
      Jkq[1] = 0.25*(c[0]+c[1]*u1+c[2]*u2)/C_sqrt2;//(3Kv-Jv)/2/sqrt(2.);
    }

    return 0;
}

/*----------------------------------------------------------------------------*/

extern double J00_off(double u1, double u2, double r){
    
    /*######################################################################
      Purpose:
        compute incident radiation tensor J00.
      Record of revisions:
        1 Dec. 2019.
      Input parameters:
        u1, u2, the limb darkening coefficients
        r, the radius.
      Return:
        return the radiation tensor J^0_0.
      References:
        LL04 Chapter 12, page 675 Eqs. 12.34-12.37.
    ######################################################################*/

    double Sr, Cr, a[3];
    Sr = 1/r;
    Cr = sqrt(1-Sr*Sr);
    if(u1==0 && u2==0) return 1/2.*(1-Cr);
    
    a[0] = 1-Cr;
    a[1] = Cr-0.5-0.5*Cr*Cr/Sr*log((1+Sr)/Cr);
    a[2] = (Cr+2)*(Cr-1)/3/(Cr+1);
    
    double Jv = 0.5*(a[0]+a[1]*u1+a[2]*u2);
    return Jv;
}

/*----------------------------------------------------------------------------*/

extern double J20_off(double u1, double u2, double r){
    
    /*######################################################################
      Purpose:
        compute incident radiation tensor J20.
      Record of revisions:
        1 Dec. 2019.
      Input parameters:
        u1, u2, the limb darkening coefficients
        r, the radius.
      Return:
        return the radiation tensor J^2_0.
      References:
        LL04 Chapter 12, page 675 Eqs. 12.34-12.37.
    ######################################################################*/

    double Sr, Cr, c[3];
    Sr = 1/r;
    Cr = sqrt(1-Sr*Sr);
    if(u1==0 && u2==0) return 1/4./C_sqrt2*Cr*Sr*Sr;
    
    c[0] = Cr*Sr*Sr;
    c[1] = 1/8.*(8*Cr*Cr*Cr-3*Cr*Cr-8*Cr+2+(4-3*Cr*Cr)*Cr*Cr/Sr*log((1+Sr)/Cr));
    c[2] = (Cr-1)/15./(Cr+1)*(9*Cr*Cr*Cr+18*Cr*Cr+7*Cr-4);
    
    double Kn = 0.25*(c[0]+c[1]*u1+c[2]*u2)/C_sqrt2;//3Kv-Jv;
    return Kn;
}

/*----------------------------------------------------------------------------*/

extern void Limb_Darkening(double Lambda, double *u1, double *u2){
    
    /*######################################################################
      Purpose:
        interpolate the limb darkening coefficients at a wavlegth of Lambda.
      Record of revisions:
        1 Dec. 2019.
      Input parameters:
        Lambda, the wavelength (m^-1).
      Output parameters:
        *u1, *u2, the limb darkening coefficients.
      References:
        Astrophysical Quantities 3rd ed. - C. Allen (Athlone Press, 1973)
    ######################################################################*/

    double Lambda_C[22] = { \
        0.20, 0.22, 0.245, 0.265, 0.28, 0.30, 0.32, 0.35, 0.37, 0.38, 0.40,\
        0.45, 0.50, 0.55, 0.60, 0.80, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0};
    double u1_C[22] = { \
        0.12, -1.3, -0.1, -0.1, 0.38, 0.74, 0.88, 0.98, 1.03, 0.92, 0.91, \
        0.99, 0.97, 0.92, 0.88, 0.73, 0.64, 0.57, 0.48, 0.35, 0.22, 0.15};
    double u2_C[22] = { \
        0.33, 1.6, 0.85, 0.90, 0.57, 0.20, 0.03, -0.1, -0.16, -0.05, -0.05,\
        -0.17, -0.22, -0.23, -0.23, -0.22, -0.20, -0.21, -0.18, -0.12, \
        -0.07, -0.07};
    
    int i, indx = -1;
    double lambda_tmp = Lambda*1e6;
    
    if(lambda_tmp<Lambda_C[0]||lambda_tmp>Lambda_C[21]){
      *u1 = 0;
      *u2 = 0;
      return;
    }
    
    for(i=0; i<21; i++){
      if(lambda_tmp>=Lambda_C[i]&&lambda_tmp<Lambda_C[i+1]){
        indx = i;
        break;
      }
    }
    double l1 = (lambda_tmp-Lambda_C[indx])/(Lambda_C[indx+1]-Lambda_C[indx]);
    double l2 = (Lambda_C[indx+1]-lambda_tmp)/(Lambda_C[indx+1]-Lambda_C[indx]);
    
    *u1 = l1*u1_C[indx+1]+l2*u1_C[indx];
    *u2 = l1*u2_C[indx+1]+l2*u2_C[indx];
    return;
}

/*----------------------------------------------------------------------------*/

extern void Thom_Scat_van(double R, double *PAB, double q){
    
    /*######################################################################
      Purpose:
        compute the Thomson scattering (a and b coefficients) in van de
            Hulst 1950BAN.
      Record of revisions:
        1 Dec. 2019.
      Input parameters:
        R, Radius.
        q, limb darkening coefficient.
      Output parameters:
        PAB[2], A and B in Eqs 11 and 12 (van de Hulst 1950BAN).
      References:
        van de Hulst 1950BAN.
    ######################################################################*/

    double gamma = asin(1/R);
    double temp1, temp2;
    
    temp1 = (1-q)/(1-1/3.*q)*(2*(1-cos(gamma)))+q/(1-1/3.*q) \
        *(1-cos(gamma)*cos(gamma)/sin(gamma)*log((1+sin(gamma))/cos(gamma)));
    temp2 = (1-q)/(1-1/3.*q)*(2/3.*(1-cos(gamma)*cos(gamma)*cos(gamma))) \
        +q/(1-1/3.*q)*(1/4.+sin(gamma)*sin(gamma)/4.-cos(gamma)*cos(gamma) \
        *cos(gamma)*cos(gamma)/4./sin(gamma)*log((1+sin(gamma))/cos(gamma)));
    
    PAB[0] = (temp1+temp2)/4.;
    PAB[1] = (temp1-temp2)/2.;
    
    return;
}

/*----------------------------------------------------------------------------*/

extern void Thom_Scat(double r, double *PAB, double u1, double u2){
    
    /*######################################################################
      Purpose:
        compute the Thomson scattering LL04.
      Record of revisions:
        1 Dec. 2019.
      Input parameters:
        R, Radius.
        u1, u2, limb darkening coefficient.
      Output parameters:
        PB[2], The A and B defined in Eqs 11 and 12 (van de Hulst 1950BAN).
      References:
        LL04, 12.34
      Warning:
        Jv = 2A+B, Kv = 2A-B.
        a and b are difined in van de Hulst 1950 a factor of pi is not
            included.
    ######################################################################*/

    double a[3], b[3], Jv, Kv;
    double Sr = 1/r;
    double Cr = sqrt(1-Sr*Sr);

    a[0] = 1-Cr;
    a[1] = Cr-0.5-0.5*Cr*Cr/Sr*log((1+Sr)/Cr);
    a[2] = (Cr+2)*(Cr-1)/3/(Cr+1);
    
    b[0] = 1./3.*(1-Cr*Cr*Cr);
    b[1] = 1./24.*(8*Cr*Cr*Cr-3*Cr*Cr-2)-1./8.*Cr*Cr*Cr*Cr/Sr*log((1+Sr)/Cr);
    b[2] = (Cr-1)*(3*Cr*Cr*Cr+6*Cr*Cr+4*Cr+2)/15./(Cr+1);
    
    Jv = 0.5*(a[0]+a[1]*u1+a[2]*u2);
    Kv = 0.5*(b[0]+b[1]*u1+b[2]*u2);
        
    PAB[0] = (Jv+Kv);
    PAB[1] = (Jv-Kv)*2.;
    
    return;
}

/*----------------------------------------------------------------------------*/
