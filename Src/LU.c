
#include "LU.h"

/*----------------------------------------------------------------------------*/

    /*######################################################################
     
     revision log:
        24 Apr. 2024.
     
    ######################################################################*/

/*----------------------------------------------------------------------------*/

extern int ludcmp_dbl(double **a, int n, int *indx){
    
    /*######################################################################
      Purpose:
        LU decomposition (double matrix).
      Record of revisions:
        30 Nov. 2019.
      Input parameters:
        a[1..n][1..n], the input matrix.
        n, the size of the matrix.
      Output parameters:
        a[1..n][1..n], the output matrix.
        indx[1..n], records the row permutation.
      Return:
        return ±1 depending on whether the number of row interchanges
            was even or odd, respectively.
      Reference:
        numerical recipes in C 2ed.
            Given a matrix a[1..n][1..n], this routine replaces it by
                the LU decomposition of a rowwise permutation of itself.
                a and n are input. a is output, arranged as in equation
                (2.3.14) above; indx[1..n] is an output vector that records
                the row permutation effected by the partial pivoting; d
                is output as ±1 depending on whether the number of row
                interchanges was even or odd, respectively. This routine
                is used in combination with lubksb to solve linear
                equations or invert a matrix.
    ######################################################################*/

    int i, j, k, d, imax = 0;
    double big,dum,sum,temp;
    double *vv = (double *)VECTOR(1, n, enum_dbl, true);
    d = 1;
    
    for(i=1;i<=n;i++){
      big = 0.0;
      for(j=1;j<=n;j++){
        if((temp = fabs(a[i][j])) > big) big=temp;
      }
      if(big == 0.0) nrerror("Singular matrix in routine ludcmp");
      vv[i] = 1.0/big;
    }
    
    for(j=1;j<=n;j++){
      for(i=1;i<j;i++){
        sum = a[i][j];
        for(k=1;k<i;k++) sum -= a[i][k]*a[k][j];
        a[i][j] = sum;
      }
      big = 0.0;
      for(i=j;i<=n;i++){
        sum = a[i][j];
        for(k=1;k<j;k++) sum -= a[i][k]*a[k][j];
        a[i][j] = sum;
        if( (dum=vv[i]*fabs(sum)) >= big){
          big = dum;
          imax = i;
        }
      }
      if(j != imax){
        for(k=1;k<=n;k++){
          dum = a[imax][k];
          a[imax][k] = a[j][k];
          a[j][k] = dum;
        }
        d = -d;
        vv[imax] = vv[j];
      }
      indx[j] = imax;
      if(a[j][j] == 0.0) a[j][j] = TINY_IN_LU;
      if(j != n){
        dum = 1.0/(a[j][j]);
        for(i=j+1;i<=n;i++) a[i][j] *= dum;
      }
    }
    
    FREE_VECTOR(vv, 1, enum_dbl);
    return d;
}

/*----------------------------------------------------------------------------*/

void lubksb_dbl(double **a, int n, int *indx, double *b){
    
    /*######################################################################
      Purpose:
        solve linear equations (double matrix).
      Record of revisions:
        30 Nov. 2019.
      Input parameters:
        a[1..n][1..n], the input matrix (LU decomposition).
        n, the size of the matrix.
        indx[1..n], records the row permutation.
        b[1..n], the right-hand side vector.
      Output parameters:
        b[1..n], the solutions.
      Reference:
        numerical recipes in C 2ed.
            Solves the set of n linear equations A·X = B. Here a[1..n][1..n]
            is input, not as the matrix A but rather as its LU decomposition,
            determined by the routine ludcmp. indx[1..n] is input as the
            permutation vector returned by ludcmp. b[1..n] is input as the
            right-hand side vector B, and returns with the solution vector
            X. a, n, and indx are not modified by this routine and can be
            left in place for successive calls with different right-hand
            sides b. This routine takes into account the possibility that
            b will begin with many zero elements, so it is efficient for
            use in matrix inversion.
    ######################################################################*/

    int i, ip, j, ii = 0;
    double sum;
    
    for(i=1;i<=n;i++){
      ip = indx[i];
      sum = b[ip];
      b[ip] = b[i];
      if(ii){
        for(j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
      }else if(sum){
        ii = i;
      }
      b[i] = sum;
    }
    
    for(i=n;i>=1;i--){
      sum = b[i];
      for(j=i+1;j<=n;j++){
        sum -= a[i][j]*b[j];
      }
      b[i] = sum/a[i][i];
    }
    return;
}

/*----------------------------------------------------------------------------*/

extern void improve_dbl(double **a, double **alud, int n, int *indx, \
    double *b, double *x){
    
    /*######################################################################
      Purpose:
        improves a solution vector x[1..n] of the linear set of equations
            A · X = B (double matrix).
      Record of revisions:
        30 Nov. 2019.
      Input parameters:
        a[1..n][1..n], the A matrix.
        alud[1..n][1..n], the LU decomposition matrix.
        n, the size of the matrix.
        indx[1..n], records the row permutation.
        b[1..n], the right-hand side vector.
        x[1..n], the solutions.
      Output parameters:
        x[1..n], the imporved solutions.
      Reference:
        numerical recipes in C 2ed.
        Improves a solution vector x[1..n] of the linear set of equations
            A · X = B. The matrix a[1..n][1..n], and the vectors b[1..n]
            and x[1..n] are input, as is the dimension n. Also input is
            alud[1..n][1..n], the LU decomposition of a as returned by
            ludcmp, and the vector indx[1..n] also returned by that
            routine. On output, only x[1..n] is modified, to an improved
            set of values.
    ######################################################################*/

    int i, j;
    double sdp;
    double *r = (double *)VECTOR(1, n, enum_dbl, true);
    
    for(i=1;i<=n;i++){
      sdp = -b[i];
      for(j=1;j<=n;j++){
        sdp += a[i][j]*x[j];
      }
      r[i] = sdp;
    }
    
    lubksb_dbl(alud,n,indx,r);
    
    for(i=1;i<=n;i++) x[i] -= r[i];
    
    FREE_VECTOR(r, 1, enum_dbl);
    return;
}

/*----------------------------------------------------------------------------*/

extern int ludcmp_flt(float **a, int n, int *indx){
    
    /*######################################################################
      Purpose:
        LU decomposition (float matrix).
      Record of revisions:
        30 Nov. 2019.
      Input parameters:
        a[1..n][1..n], the input matrix.
        n, the size of the matrix.
      Output parameters:
        a[1..n][1..n], the output matrix.
        indx[1..n], records the row permutation.
      Return:
        return ±1 depending on whether the number of row interchanges
            was even or odd, respectively.
      Reference:
        numerical recipes in C 2ed.
    ######################################################################*/

    int i, j, k, d, imax = 0;
    float big, dum, sum, temp;
    float *vv = (float *)VECTOR(1, n, enum_flt,  true);
    d = 1;
    
    for(i=1;i<=n;i++){
      big = 0.0;
      for(j=1;j<=n;j++){
        if((temp=fabs(a[i][j])) > big){
          big = temp;
        }
      }
      if(big == 0.0) nrerror("Singular matrix in routine ludcmp");
      vv[i] = 1.0/big;
    }
    
    for(j=1;j<=n;j++){
      for(i=1;i<j;i++){
        sum = a[i][j];
        for(k=1;k<i;k++){
          sum -= a[i][k]*a[k][j];
        }
        a[i][j] = sum;
      }
      big = 0.0;
      for(i=j;i<=n;i++){
        sum = a[i][j];
        for(k=1;k<j;k++) sum -= a[i][k]*a[k][j];
        a[i][j] = sum;
        if( (dum=vv[i]*fabs(sum)) >= big){
          big = dum;
          imax = i;
        }
      }
      if(j != imax){
        for(k=1;k<=n;k++){
          dum = a[imax][k];
          a[imax][k] = a[j][k];
          a[j][k] = dum;
        }
        d = -d;
        vv[imax] = vv[j];
      }
      indx[j] = imax;
      if(a[j][j] == 0.0) a[j][j] = TINY_IN_LU;
      if(j != n){
        dum = 1.0/(a[j][j]);
        for(i=j+1;i<=n;i++) a[i][j] *= dum;
      }
    }
    
    FREE_VECTOR(vv, 1, enum_flt);
    return d;
}

/*----------------------------------------------------------------------------*/

extern void lubksb_flt(float **a, int n, int *indx, float *b){
    
    /*######################################################################
      Purpose:
        solve linear equations (float matrix).
      Record of revisions:
        30 Nov. 2019.
      Input parameters:
        a[1..n][1..n], the input matrix (LU decomposition).
        n, the size of the matrix.
        indx[1..n], records the row permutation.
        b[1..n], the right-hand side vector.
      Output parameters:
        b[1..n], the solutions.
      Reference:
        numerical recipes in C 2ed.
    ######################################################################*/

    int i, ip, j, ii = 0;
    float sum;
    
    for(i=1;i<=n;i++){
      ip = indx[i];
      sum = b[ip];
      b[ip] = b[i];
      if(ii){
        for(j=ii;j<=i-1;j++){
          sum -= a[i][j]*b[j];
        }
      }else if(sum){
        ii = i;
      }
      b[i] = sum;
    }
    
    for(i=n;i>=1;i--){
      sum = b[i];
      for(j=i+1;j<=n;j++){
        sum -= a[i][j]*b[j];
      }
      b[i] = sum/a[i][i];
    }
    
    return;
}

/*----------------------------------------------------------------------------*/

extern void improve_flt(float **a, float **alud, int n, int *indx, \
    float *b, float *x){
    
    /*######################################################################
      Purpose:
        improves a solution vector x[1..n] of the linear set of equations
            A · X = B (float matrix).
      Record of revisions:
        30 Nov. 2019.
      Input parameters:
        a[1..n][1..n], the A matrix.
        alud[1..n][1..n], the LU decomposition matrix.
        n, the size of the matrix.
        indx[1..n], records the row permutation.
        b[1..n], the right-hand side vector.
        x[1..n], the solutions.
      Output parameters:
        x[1..n], the imporved solutions.
      Reference:
        numerical recipes in C 2ed.
    ######################################################################*/

    int i, j;
    float sdp;
    float *r = (float *)VECTOR(1, n, enum_flt, true);
    
    for(i=1;i<=n;i++){
      sdp = -b[i];
      for(j=1;j<=n;j++){
        sdp += a[i][j]*x[j];
      }
      r[i] = sdp;
    }
    
    lubksb_flt(alud,n,indx,r);
    
    for(i=1;i<=n;i++) x[i] -= r[i];
    
    FREE_VECTOR(r, 1, enum_flt);
    return;
}

/*----------------------------------------------------------------------------*/

extern int ludcmp_cplx(complex double **a, int n, int *indx){
    
    /*######################################################################
      Purpose:
        LU decomposition (complex double matrix).
      Record of revisions:
        30 Nov. 2019.
      Input parameters:
        a[1..n][1..n], the input matrix.
        n, the size of the matrix.
      Output parameters:
        a[1..n][1..n], the output matrix.
        indx[1..n], records the row permutation.
      Return:
        return ±1 depending on whether the number of row interchanges
            was even or odd, respectively.
      Reference:
        numerical recipes in C 2ed.
    ######################################################################*/

    int i, j, k, d, imax = 0;
    double big, temp;
    complex double sum,dum;
    double *vv = (double *)VECTOR(1, n, enum_dbl, true);
    d = 1;
    
    for(i=1;i<=n;i++){
      big = 0.0;
      for(j=1;j<=n;j++){
        if((temp=creal(a[i][j])*creal(a[i][j])+cimag(a[i][j]) \
            *cimag(a[i][j])) > big){
          big = temp;
        }
      }
      if(big == 0.0) nrerror("Singular matrix in routine ludcmp");
      vv[i] = 1.0/big;
    }
    
    for(j=1;j<=n;j++){
      for(i=1;i<j;i++){
        sum = a[i][j];
        for(k=1;k<i;k++){
          sum -= a[i][k]*a[k][j];
        }
        a[i][j] = sum;
      }
      big = 0.0;
      for(i=j;i<=n;i++){
        sum = a[i][j];
        for(k=1;k<j;k++){
          sum -= a[i][k]*a[k][j];
        }
        a[i][j] = sum;
        if( (temp=vv[i]*(creal(sum)*creal(sum)+cimag(sum)*cimag(sum)))\
            >= big){
          big = temp;
          imax = i;
        }
      }
      if(j != imax){
        for(k=1;k<=n;k++){
          dum = a[imax][k];
          a[imax][k] = a[j][k];
          a[j][k] = dum;
        }
        d = -d;
        vv[imax] = vv[j];
      }
      indx[j] = imax;
      if(a[j][j] == 0.0) a[j][j] = TINY_IN_LU;
      if(j != n){
        dum = 1.0/(a[j][j]);
        for(i=j+1;i<=n;i++) a[i][j] *= dum;
      }
    }
    
    FREE_VECTOR(vv, 1, enum_dbl);
    return d;
}

/*----------------------------------------------------------------------------*/

extern void lubksb_cplx(complex double **a, int n, int *indx, \
    complex double *b){
    
    /*######################################################################
      Purpose:
        solve linear equations (complex double matrix).
      Record of revisions:
        30 Nov. 2019.
      Input parameters:
        a[1..n][1..n], the input matrix (LU decomposition).
        n, the size of the matrix.
        indx[1..n], records the row permutation.
        b[1..n], the right-hand side vector.
      Output parameters:
        b[1..n], the solutions.
      Reference:
        numerical recipes in C 2ed.
    ######################################################################*/

    int i, ip, j, ii = 0;
    complex double sum;
    
    for(i=1;i<=n;i++){
      ip = indx[i];
      sum = b[ip];
      b[ip] = b[i];
      if(ii){
        for(j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
      }else if(sum){
        ii = i;
      }
      b[i] = sum;
    }
    for(i=n;i>=1;i--){
      sum = b[i];
      for(j=i+1;j<=n;j++){
        sum -= a[i][j]*b[j];
      }
      b[i] = sum/a[i][i];
    }
    
    return;
}

/*----------------------------------------------------------------------------*/

extern void improve_cplx(complex double **a, complex double **alud, \
    int n, int *indx, complex double *b, complex double *x){
    
    /*######################################################################
      Purpose:
        improves a solution vector x[1..n] of the linear set of equations
            A · X = B (complex matrix).
      Record of revisions:
        30 Nov. 2019.
      Input parameters:
        a[1..n][1..n], the A matrix.
        alud[1..n][1..n], the LU decomposition matrix.
        n, the size of the matrix.
        indx[1..n], records the row permutation.
        b[1..n], the right-hand side vector.
        x[1..n], the solutions.
      Output parameters:
        x[1..n], the imporved solutions.
      Reference:
        numerical recipes in C 2ed.
    ######################################################################*/

    int i, j;
    complex double sdp, *r = (complex double *)VECTOR(1, n, enum_cplx, true);
    
    for(i=1;i<=n;i++){
      sdp = -b[i];
      for(j=1;j<=n;j++){
        sdp += a[i][j]*x[j];
      }
      r[i] = sdp;
    }
    
    lubksb_cplx(alud,n,indx,r);
    
    for(i=1;i<=n;i++) x[i] -= r[i];
    
    FREE_VECTOR(r, 1, enum_cplx);
    return;
}

/*----------------------------------------------------------------------------*/

