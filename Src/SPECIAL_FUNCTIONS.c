
#include "SPECIAL_FUNCTIONS.h"

/*----------------------------------------------------------------------------*/

    /*######################################################################
     
      revision log:

        19 Aug. 2023.
     
    ######################################################################*/

/*----------------------------------------------------------------------------*/

static int sign(double X);

/*----------------------------------------------------------------------------*/

static int sign(double X){
    
    /*######################################################################
      Purpose:
        return the sign of input number X.
      Record of revisions:
        17 Oct. 2019.
      Input parameters:
        X, the input number.
      Return:
        return 1 (X positive), -1 (X negative) 0 (X=0).
    ######################################################################*/

    if(X>0){
      return 1;
    }else if(X<0){
      return -1;
    }else{
      return 0;
    }
}

/*----------------------------------------------------------------------------*/

extern double Legendre_Hoeksema(int L, int M, double Theta, \
    STR_FCTSG *fctsg){
    
    /*######################################################################
      Purpose:
        computes the associated Legendre polynomial Plm(x) with the
            nomalization according to Hoeksema 1984, where Condon–Shortley
            phase factor do not be concerned.
      Record of revisions:
        19 Aug. 2023.
      Input parameters:
        L, M, are integers with L>0 and -L<M<L.
        Theta, the inclination angle.
        fctsg, structure with factorial and double factorial.
      Return:
        return Associated Legendre polynomial Plm(x).
      Reference:
        numerical recipes in C 2ed.
    ######################################################################*/

    if(L<M) return 0.;
    
    int i, Z;
    double Legendre = 0, Coefficent, q;
    double Ct = cos(Theta);
    double St = sin(Theta);
    
    for(Z=0; Z<=(L-M)/2; Z++){
      Coefficent = pow(Ct,(L-M-2.*Z));
      for(i=1; i<=Z; i++){
        Coefficent *= (L-M-2.*i+1.)*(2.*i-2-L+M)/2./i/(2.*L-2.*i+1);
      }
      Legendre += Coefficent;
    }
    
    q = (M==0?1.:2.);
    
    Legendre *= fctsg->fct2[2*L-1]*sqrt(q/fctsg->fct[L-M]/fctsg->fct[L+M]) \
        *pow(St, M);

    return Legendre;
}

/*----------------------------------------------------------------------------*/

extern double Legendre_Hoeksema_P(int L, int M, double Theta, \
    STR_FCTSG *fctsg){
    
    /*######################################################################
      Purpose:
        computes the derivative associated Legendre polynomial Plm(x) with
            the nomalization according to Hoeksema 1984, where 
            Condon–Shortley phase factor do not be concerned
      Record of revisions:
        5 Jun. 2019.
      Input parameters:
        L, M, are integers with L>0 and -L<M<L.
        Theta, the inclination angle.
        fctsg, structure with factorial and double factorial.
      Return:
        return the derivative of the associated Legendre polynomial Plm(x).
      Reference:
        numerical recipes in C 2ed.
    ######################################################################*/

    if(L<M) return 0.;
    
    int i, Z;
    double Legendre1 = 0, Legendre2 = 0, Wlm, Coefficent1 = 0;
    double Coefficent2 = 0, q;
    double Ct = cos(Theta);
    double St = sin(Theta);

    for(Z=0; Z<=(L-M-1)/2; Z++){
      Coefficent1 = -(L-M-2.*Z)*pow(Ct,(L-M-2.*Z-1.))*St;
      for(i=1; i<=Z; i++){
        Coefficent1 *= (L-M-2.*i+1.)*(2.*i-2-L+M)/2./i/(2.*L-2.*i+1);
      }
      Legendre1 += Coefficent1;
    }
        
    if(M==0){
      q = 1.;
      Legendre2 = 0;
    }else{
      Legendre1 *= pow(St,M);
      q=2.;
      for(Z=0; Z<=(L-M)/2; Z++){
        Coefficent2 = pow(Ct,(L-M-2.*Z));
        for(i=1; i<=Z; i++){
          Coefficent2 *= (L-M-2.*i+1.)*(2.*i-2-L+M)/2./i/(2.*L-2.*i+1);
        }
        Legendre2 += Coefficent2;
      }
      Legendre2 *= M*pow(St,(M-1.))*Ct;
    }
    
    Wlm = fctsg->fct2[2*L-1]*sqrt(q/fctsg->fct[L+M]/fctsg->fct[L-M]);
    
    return (Legendre1+Legendre2)*Wlm;
}

/*----------------------------------------------------------------------------*/

extern double Associated_Legendre(int L, int M, double X){
    
    /*######################################################################
      Purpose:
        computes the associated Legendre polynomial Plm(x).
      Record of revisions:
        5 Jun. 2019.
      Input parameters:
        L, M, are integers with L>0 and -L<M<L.
        X, X=cos(theta), theta is the inclination angle -1<X<1.
      Return:
        return the associated Legendre polynomial Plm(x).
      Reference:
        numerical recipes in C 2ed.
    ######################################################################*/

    if(M>L ||-M>L) return 0;
    if(M<0){
      int i;
      double tmp=1.0;
      if( M == 0 ){
        tmp=1.;
      }else if(M > 0){
        for(i=L-M+1; i<=L+M; i++){
          tmp*=i;
        }
      }else{
       for(i=L+M+1; i<=L-M; i++){
          tmp/=i;
        }
      }
      int sgn = M%2==0?1:-1;
      return Associated_Legendre(L,-M,X)*sgn*tmp;
    }
    
    if(fabs(X) > 1.0){
      fprintf(stderr,"gg %d %d %e\n",L,M,X);
      fprintf(stderr,"Bad arguments in routine Associated_Legendre \n");
        
      fprintf(stderr,"...now exiting to system...\n");
        
      exit(1);
    }
    
    double fact,pll=0,pmm=1.0,pmmp1,somx2;
    int i,ll;
    
    somx2=sqrt((1.-X)*(1.+X));
    fact=1.0;
    for(i=1;i<=M;i++){
      pmm *= -fact*somx2;
      fact += 2.0;
    }
    
    if(L == M){
      return pmm;
    }else{
      pmmp1=X*(2*M+1)*pmm;
      if(L == (M+1)){
        return pmmp1;
      }else{
        for(ll=M+2;ll<=L;ll++){
          pll = (X*(2*ll-1)*pmmp1-(ll+M-1)*pmm)/(ll-M);
          pmm = pmmp1;
          pmmp1 = pll;
        }
        return pll;
      }
    }
}

/*----------------------------------------------------------------------------*/

extern double Harmonic_Coefficient(int L, int M){
    
    /*######################################################################
      Purpose:
        calculate the Spherical Harmonics.
      Record of revisions:
        19 Aug. 2023.
      Input parameters:
        L, M, are integers with L>0 and -L<M<L.
      Return:
        return the spherical harmonic.
      Reference:
        numerical recipes in C 2ed.
    ######################################################################*/
    
    int i;
    double tmp = 1.0;
    
    if(M==0){
      tmp = 1.;
    }else if(M > 0){
      for(i=L-M+1; i<=L+M; i++){
        tmp /= i;
      }
    }else{
      for(i=L+M+1; i<=L-M; i++){
        tmp *= i;
      }
    }
    
    tmp = sqrt(tmp*(2.*L+1.)/4./C_Pi);
    
    return tmp;
}

/*----------------------------------------------------------------------------*/

extern complex double Spherical_Harmonic(int L, int M, double X, \
    double Phi){
    
    /*######################################################################
      Purpose:
        calculate the Spherical Harmonics.
      Record of revisions:
        19 Aug. 2023.
      Input parameters:
        L, M, are integers with L>0 and -L<M<L.
        X, X=cos(theta), theta is the inclination angle -1<X<1.
        Phi, the azimuth 0<Phi<2*Pi.
      Return:
        return the spherical harmonic.
      Reference:
        numerical recipes in C 2ed.
    ######################################################################*/

    double tmp = Harmonic_Coefficient(L, M)*Associated_Legendre(L, M, X);
        
    return tmp*cos(M*Phi)+tmp*sin(M*Phi)*I;
}

/*----------------------------------------------------------------------------*/

extern complex double Spherical_Harmonic_conjugate(int L, int M, \
    double X, double Phi){
    
    /*######################################################################
      Purpose:
        calculate the conjugate of the Spherical Harmonics.
      Record of revisions:
        5 Jun. 2019.
      Input parameters:
        L, M, are integers with L>0 and -L<M<L.
        X, X=cos(theta), theta is the inclination angle -1<X<1.
        Phi, the azimuth 0<Phi<2*Pi.
      Return:
        return the spherical harmonic.
      Reference:
        numerical recipes in C 2ed.
    ######################################################################*/

    double tmp = Harmonic_Coefficient(L, M)*Associated_Legendre(L, M, X);
        
    return tmp*cos(M*Phi)-tmp*sin(M*Phi)*I;
}

/*----------------------------------------------------------------------------*/

extern void Spherical_Bessel(double x, int Lmax, double *J, double *Jp){
    
    /*######################################################################
      Purpose:
        steed algorithm to calculate the Spherical Bessel functions and
            their derivatives.
      Record of revisions:
        26 Sept. 2019.
      Input parameters:
        x, the x value of the Spherical Bessel functions.
        Lmax, order of L.
      Output parameters:
        J[], the Spherical Bessel functions.
        Jp[], the derivatives of the Spherical Bessel functions.
      Reference:
        Chapter 9 of 'Computational Atomic Physics' ed. K Bartschat, which
            is called 'The Calculation of Spherical Bessel and Coulomb
            Functions' by A R Barnett.
    ######################################################################*/

    if(Lmax < 0){
      fprintf(stderr, "Errors in calculating the spherical " \
          "Bessel function. \n");
      fprintf(stderr, "Max of L is negative. \n");
      exit(0);
    }
    
    const double ACCUR = 1e-15;
    
    if(Lmax == 0){   
      if(x < sqrt(ACCUR)){
        J[0] = 1;
        Jp[0] = 0;
      }else{
        J[0] = sin(x)/x;
        Jp[0] = (cos(x)-sin(x)/x)/x;
      }
      return;
        
    }else if(Lmax ==1 ){ 
      if(x < sqrt(ACCUR)){
        J[0] = 1;
        Jp[0] = 0;
        J[1] = x/3.;
        Jp[1] = 1./3;
      }else{
        J[0] = sin(x)/x;
        Jp[0] = (cos(x)-sin(x)/x)/x;
        J[1] = -Jp[0];
        Jp[1] = J[0]-2./x*J[1];
      }
      return;
        
    }else{    
      int i = 0;
      if(x == 0){      
        J[0] = 1;
        Jp[0] = 0;
        J[1] = 0;
        Jp[0] = 1./3;
        for(i = 2; i <= Lmax; i++){
          J[i] = 0;
        }
        return;

      }else if(x > sqrt(ACCUR)){
            
        const double Tiny = 1e-15;
        double xinv = 1./x;
        double Sn = Lmax*xinv, twox = 2.*xinv;
        double Tk = Sn*2+xinv;
        double F = Sn, C = Sn, D = 0, Delta, Den = 1.0;
        int limits = 2000;

        if(fabs(F) < Tiny){
          F = Tiny;
        }
            
        do{
          Tk += twox;
          C = Tk-1./C;
          if(fabs(C) < Tiny){
            C = Tiny;
          }
          D=(Tk-D);
          if(fabs(D) < Tiny){
            D = Tiny;
          }
          D = 1./D;
          Delta = C*D;
          F *= Delta;
          if(D <= 0){
            Den = -Den;
          }
          i++;
                
        }while(fabs(Delta-1) > ACCUR && i <= limits);
            
        J[Lmax] = Den;
        Jp[Lmax] = F*Den;
            
        for(i = Lmax; i > 0; i--){
          J[i-1] = (Sn+xinv)*J[i]+Jp[i];
          Sn -= xinv;
          Jp[i-1] = Sn*J[i-1]-J[i];
        }
            
        Den = J[0];
        double Den2=J[1];
        J[0] = xinv*sin(x);
        Jp[0] = (cos(x)-sin(x)*xinv)*xinv;
            
        double OMEGA;
            
        if(fabs(J[0])>1e-5){
          OMEGA = J[0]/Den;
          for(i = 1; i <= Lmax; i++){
            J[i] *= OMEGA;
            Jp[i] *= OMEGA;
          }
          return;
        }else{
          J[1] = -Jp[0];
          Jp[1] = J[0]-twox*J[1];
          OMEGA = J[1]/Den2;
          for(i = 2; i <= Lmax; i++){
            J[i] *= OMEGA;
            Jp[i] *= OMEGA;
          }
          return;
        }
        return;
            
      }else{
            
        J[0] = 1;
        Jp[0] = 0;
        for(i=1; i<=Lmax; i++){
          J[i] = J[i-1]/(2.*i+1)*x;
          Jp[i] = J[i-1]*i/(2.*i+1);
        }
            
      return;
      }
    }
}

/*----------------------------------------------------------------------------*/

extern void Sph_Bessel_Output(char *filename, double x0, double x1, \
    int num, int indx_l){
    
    /*######################################################################
      Purpose:
        write the spherical bessel function and its derivative (order l)
            into a file.
      Record of revisions:
        25 Oct. 2019.
      Input parameters:
        filename, the file nema.
        x0, x1, the range of x.
        num, number of the points,
        indx_l, the order of L.
    ######################################################################*/

    FILE *Fa;
    Fa=fopen(filename, "w");
    double x, delta_x = (x1-x0)/num;
    double JJ[indx_l],JJp[indx_l];
    for(x=x0; x<=x1; x+=delta_x){
      Spherical_Bessel(x, indx_l, JJ, JJp);
      fprintf(Fa, "%e %e %e\n",x,JJ[indx_l],JJp[indx_l]);
    }
    fclose(Fa);
    return;
}

/*----------------------------------------------------------------------------*/

extern int Compute_SB_Zeros(double **Zeros, int Zero_L, int Zero_N, \
    enum boundary_condition type){
    
    /*######################################################################
      Purpose:
        compute the positive zeros of spherical Bessel funtion or of its
            derivative.
      Record of revisions:
        19 Aug. 2023.
      Input parameters:
        Zero_L, The max order L of the Bessel L.
        Zero_N, The max number of the zeros.
        type, If the type = enum_deri_zeros, output the zeros of the 
            derivative, else output the zeros of the function.
      Output parameters:
        Zeros[0...Zero_L][1...Zero_N], The positive zeros of spherical
            Bessel funtion or of its derivative.
            X01 is set to 0 if type = enum_deri_zeros.
    ######################################################################*/

    double JJ[Zero_L+1], JJp[Zero_L+1];
    double delta_x;
    double x0, y0;
    double x1, x2, y1, y2;
    double x_tmp, y_tmp, X_width ;
    int y1_sign, y2_sign, y_tmp_sign;
    int indx_l, indx_n;
    double X_limits = 1e-15, Y_limits = 1e-15;
    for(indx_l=0; indx_l<=Zero_L; indx_l++){
      indx_n = 0;
      if(indx_l==0 && type==enum_deri_zeros){
        indx_n = 1;
        Zeros[0][1] = 0;
      }
      delta_x = 1e-3;
      x0 = delta_x;
        
      x1 = x0;
      Spherical_Bessel(x1, indx_l, JJ, JJp);
      if(type==enum_deri_zeros){
        y1 = JJp[indx_l];
      }else{
        y1 = JJ[indx_l];  
      }
      y1_sign = sign(y1);
        
      do{
        x2 = x1+delta_x;
        Spherical_Bessel(x2, indx_l, JJ, JJp);

        if(type==enum_deri_zeros){
          y2 = JJp[indx_l];
        }else{
          y2 = JJ[indx_l];  
        }
        y2_sign = sign(y2);
          
        if(y1_sign*y2_sign<=0){
          X_width = delta_x;
          do{
            x_tmp = 0.5*(x1+x2);
            Spherical_Bessel(x_tmp, indx_l, JJ, JJp);
            if(type==enum_deri_zeros){
              y_tmp = JJp[indx_l];
            }else{
              y_tmp = JJ[indx_l];  
            }

            y_tmp_sign = sign(y_tmp);
                    
            if(y1_sign*y_tmp_sign<=0){
              x2 = x_tmp;
              y2 = y_tmp;
              y2_sign = y_tmp_sign;
            }else{
              x1 = x_tmp;
              y1 = y_tmp;
              y1_sign = y_tmp_sign;
            }
                    
            if(fabs(y1)<=fabs(y2)){
              x0 = x1;
              y0 = y1;
            }else{
              x0 = x2;
              y0 = y2;
            }
            X_width /= 2.;
          }while(X_width>X_limits &&fabs(y0)>Y_limits);
                
          indx_n++;
          Zeros[indx_l][indx_n] = x0;
          x1 = x0+delta_x;      
          Spherical_Bessel(x1, indx_l, JJ, JJp);
          if(type==enum_deri_zeros){
            y1 = JJp[indx_l];
          }else{
            y1 = JJ[indx_l];  
          }
          y1_sign = sign(y1);
                
        }else{
          x1 = x2;
          y1 = y2;
          y1_sign = y2_sign;
        }
            
      }while(indx_n<Zero_N);
    }
    return 0;
}

/*----------------------------------------------------------------------------*/

extern int cisi(double x, double *ci, double *si){

    /*######################################################################
      Purpose:
        Computes the cosine and sine integrals Ci(x) and Si(x). Ci(0) is 
        returned as a large negative number and no error message is 
        generated. For x < 0 the routine returns Ci(−x) and you must 
        supply the −iπ yourself.
      Record of revisions:
        19 Aug. 2023.
      Input parameters:
        x, The x value.
      Output parameters:
        ci, The Ci(x)
        si, The Si(x)
      Reference:
        numerical recipes in C 2ed P258.
    ######################################################################*/

    double EPS = 6.0e-8;
    double EULER = 0.57721566; 
    int MAXIT = 100; 
    double PIBY2 = 1.5707963; 
    double FPMIN = 1.0e-30; 
    double TMIN = 2.0;

    int i,k,odd;
    double a,err,fact,sign,sum,sumc,sums,t,term; 
    complex double h,b,c,d,del;

    t = fabs(x);
    if(t == 0.0){
      *si = 0.0;
      *ci = -1.0/FPMIN; 
      return 0;
    }

    if(t > TMIN){
      b = 1.0+t*I; 
      c = 1.0/FPMIN; 
      d = 1.0/b;
      h = 1.0/b;
      for(i=2;i<=MAXIT;i++){
        a = -(i-1)*(i-1); 
        b += 2.0; 
        d = 1.0/(a*d+b);
        c = b+a/c;
        del = c*d;
        h *= del;

        if(fabs(creal(del)-1.0)+fabs(cimag(del)) < EPS) break; 
      }
      if(i > MAXIT) nrerror("cf failed in cisi"); 
      h *= cos(t)-sin(t)*I;
      *ci = -creal(h);
      *si = PIBY2+cimag(h);
    }else{
      if(t < sqrt(FPMIN)){
        sumc = 0.0;
        sums = t; 
      }else{
        sum = 0.0;
        sums = 0.0;
        sumc = 0.0; 
        sign = 1.0;
        fact = 1.0;
        odd = true;
        for(k=1;k<=MAXIT;k++){
          fact *= t/k; 
          term = fact/k;
          sum += sign*term;
          err = term/fabs(sum); 
          if(odd){
            sign = -sign; 
            sums = sum; 
            sum = sumc;
          }else{ 
            sumc = sum;
            sum = sums;
          }
          if(err < EPS) break; 
          odd =!odd;
        }
        if(k > MAXIT) nrerror("maxits exceeded in cisi"); 
      }
      *si = sums;
      *ci = sumc+log(t)+EULER; 
    }
    if(x < 0.0) *si = -(*si); 

    return 0;

}

/*----------------------------------------------------------------------------*/
