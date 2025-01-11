
#include "WIGNER.h"

/*----------------------------------------------------------------------------*/

    /*######################################################################
     
      revision log:

        17 Aug. 2023.
           --- bugfix: fix a typo in WIGNER_3J

     
    ######################################################################*/

/*----------------------------------------------------------------------------*/
/*
static int ISINTEGER(double val) {

    if (floor(val) == ceil(val)) {
        return 1; 
    } else {
        return 0;
    }
}
*/
/*----------------------------------------------------------------------------*/

extern double WIGNER_3J(double J1, double J2, double J3, double M1, \
    double M2, double M3, STR_FCTSG *fctsg){
    
    /*######################################################################
      Purpose:
        computes the 3j symbol.
      Record of revisions:
        25 Apr. 2018.
      Input parameters:
        J1, J2, J3, The total angular momentum.
        M1, M2, M3, The magnetic quantum number.
      Return:
        return the 3j symbol.
      References:
        LL04 Chapter 2, Page 36, Equation 2.19 and Page 38, Equation 2.22.
    ######################################################################*/

    int indx1 = lround(J1*2);
    int indx2 = lround(J2*2);
    int indx3 = lround(J3*2);
    int indx4 = lround(M1*2);
    int indx5 = lround(M2*2);
    int indx6 = lround(M3*2);

    if((indx4+indx5+indx6)!=0) return 0.0;
    if((indx1+indx2)<indx3||abs(indx1-indx2)>indx3) return 0.0;

    if(abs(indx4)>indx1) return 0.0;
    if(abs(indx5)>indx2) return 0.0;
    if(abs(indx6)>indx3) return 0.0;

    if((indx1+indx2+indx3)%2!=0) return 0.0;
    if((indx1-indx4)%2!=0) return 0.0;
    if((indx2-indx5)%2!=0) return 0.0;

    if(fctsg->memo){
      double *ptr;
      ptr = ele6d(fctsg->J3, indx1, indx2, indx3, indx4, indx5, indx6);

      if(ptr!=NULL){
        return *ptr;
      }
    }

    int IC = lround(J1-M1);
    int IE = lround(J3-J2+M1);
    int IG = lround(J1+J2-J3);
    int IF = lround(J3-J1-M2);
    int IH = lround(J2+M2);
 
    int t, t1 = 0, t2 = IG;

    if(t1<-IE) t1 = -IE;
    if(t1<-IF) t1 = -IF;
        
    if(t2>IC) t2 = IC;
    if(t2>IH) t2 = IH;


    double c, d, e = 0, val;
    
    c = fctsg->fct[IG]*fctsg->fct[lround(J1-J2+J3)] \
        *fctsg->fct[lround(-J1+J2+J3)]/fctsg->fct[lround(J1+J2+J3+1)];
    d = fctsg->fct[lround(J1+M1)]*fctsg->fct[IC] \
        *fctsg->fct[IH]*fctsg->fct[lround(J2-M2)] \
        *fctsg->fct[lround(J3+M3)]*fctsg->fct[lround(J3-M3)];

    for(t=t1; t<=t2; t++){
      e += fctsg->sg[t]/(fctsg->fct[t]*fctsg->fct[IG-t] \
        *fctsg->fct[IC-t]*fctsg->fct[IH-t]*fctsg->fct[IE+t] \
        *fctsg->fct[IF+t]);
    }
    
    val = fctsg->sg[lround(J1-J2-M3)]*sqrt(c*d)*e;

    if(fctsg->memo){
      inster_6d(fctsg->J3, val, indx1, indx2, indx3, indx4, indx5, indx6);
    }

    return val;
}

/*----------------------------------------------------------------------------*/

extern double WIGNER_6J(double J1, double J2, double J3, double J4, \
    double J5, double J6, STR_FCTSG *fctsg){
    
    /*######################################################################
      Purpose:
        Computes the 6j symbol.
      Record of revisions:
        25 Apr. 2018.
      Input parameters:
        J1, J2, J3, J4, J5, J6, the total angular momentum.
      Return:
        return the 6j symbol.
      References:
        LL04 Chapter 2, Page 42, Equation 2.35.
    ######################################################################*/

    int indx1 = lround(J1*2);
    int indx2 = lround(J2*2);
    int indx3 = lround(J3*2);
    int indx4 = lround(J4*2);
    int indx5 = lround(J5*2);
    int indx6 = lround(J6*2);
    
    if((indx1+indx2)<indx3||abs(indx1-indx2)>indx3) return 0.0;
    if((indx1+indx5)<indx6||abs(indx1-indx5)>indx6) return 0.0;
    if((indx4+indx2)<indx6||abs(indx4-indx2)>indx6) return 0.0;
    if((indx4+indx5)<indx3||abs(indx4-indx5)>indx3) return 0.0;

    if((indx1+indx2+indx3)%2!=0) return 0.0;
    if((indx1+indx5+indx6)%2!=0) return 0.0;
    if((indx2+indx4+indx6)%2!=0) return 0.0;
    if((indx3+indx4+indx5)%2!=0) return 0.0;

    if(fctsg->memo){

      double *ptr;
      ptr = ele6d(fctsg->J6, indx1, indx2, indx3, indx4, indx5, indx6);

      if(ptr!=NULL){
        return *ptr;
      }
    }

    int sum1 = lround(J1+J2+J3);
    int sum2 = lround(J1+J5+J6);
    int sum3 = lround(J4+J2+J6);
    int sum4 = lround(J4+J5+J3);
    int II = lround(J1+J2+J4+J5);
    int IJ = lround(J2+J3+J5+J6);
    int IK = lround(J1+J3+J4+J6);

    int t, t1 = sum1, t2 = II;
    if(t1<sum2) t1 = sum2;
    if(t1<sum3) t1 = sum3;
    if(t1<sum4) t1 = sum4;

    if(t2>IJ) t2 = IJ;
    if(t2>IK) t2 = IK;

    double a, b = 0, val;

    a = fctsg->fct[lround(J1+J2-J3)]*fctsg->fct[lround(J1-J2+J3)] \
        *fctsg->fct[lround(-J1+J2+J3)]/fctsg->fct[sum1+1] \
        *fctsg->fct[lround(J1+J5-J6)]*fctsg->fct[lround(J1-J5+J6)] \
        *fctsg->fct[lround(-J1+J5+J6)]/fctsg->fct[sum2+1] \
        *fctsg->fct[lround(J4+J2-J6)]*fctsg->fct[lround(J4-J2+J6)] \
        *fctsg->fct[lround(-J4+J2+J6)]/fctsg->fct[sum3+1] \
        *fctsg->fct[lround(J4+J5-J3)]*fctsg->fct[lround(J4-J5+J3)] \
        *fctsg->fct[lround(-J4+J5+J3)]/fctsg->fct[sum4+1];

    for(t=t1; t<=t2; t++){
      b += fctsg->sg[t]*fctsg->fct[t+1]/(fctsg->fct[t-sum1] \
          *fctsg->fct[t-sum2]*fctsg->fct[t-sum3] \
          *fctsg->fct[t-sum4]*fctsg->fct[II-t] \
          *fctsg->fct[IJ-t]*fctsg->fct[IK-t]);
    }
 
    val = sqrt(a)*b;

    if(fctsg->memo){
      inster_6d(fctsg->J6, val, indx1, indx2, indx3, indx4, indx5, indx6);
    }

    return val;
}

/*----------------------------------------------------------------------------*/

extern double WIGNER_9J(double J1, double J2, double J3, double J4, \
    double J5, double J6, double J7, double J8, double J9, \
    STR_FCTSG *fctsg){
    
    /*######################################################################
      Purpose:
        computes the 9j symbol.
      Record of revisions:
        6 Dec. 2022.
      Input parameters:
        J1, J2, J3, J4, J5, J6, J7, J8, J9, The total angular momentum.
      Return:
        return the 9j symbol.
      References:
        LL04 Chapter 2, Page 47, Equation 2.48.
    ######################################################################*/

    int indx1 = lround(J1*2);
    int indx2 = lround(J2*2);
    int indx3 = lround(J3*2);
    int indx4 = lround(J4*2);
    int indx5 = lround(J5*2);
    int indx6 = lround(J6*2);
    int indx7 = lround(J7*2);
    int indx8 = lround(J8*2);
    int indx9 = lround(J9*2);


    if((indx1+indx2)<indx3||abs(indx1-indx2)>indx3) return 0.0;
    if((indx4+indx5)<indx6||abs(indx4-indx5)>indx6) return 0.0;
    if((indx7+indx8)<indx9||abs(indx7-indx8)>indx9) return 0.0;
    if((indx1+indx4)<indx7||abs(indx1-indx4)>indx7) return 0.0;
    if((indx2+indx5)<indx8||abs(indx2-indx5)>indx8) return 0.0;
    if((indx3+indx6)<indx9||abs(indx3-indx6)>indx9) return 0.0;

    if((indx1+indx2+indx3)%2!=0) return 0.0;
    if((indx4+indx5+indx6)%2!=0) return 0.0;
    if((indx7+indx8+indx9)%2!=0) return 0.0;
    if((indx1+indx4+indx7)%2!=0) return 0.0;
    if((indx2+indx5+indx8)%2!=0) return 0.0;
    if((indx3+indx6+indx9)%2!=0) return 0.0;

    if(fctsg->memo){
      double *ptr;
      ptr = ele9d(fctsg->J9, indx1, indx2, indx3, indx4, indx5, indx6, \
          indx7, indx8, indx9);

      if(ptr!=NULL){
        return *ptr;
      }
    }

    int IA = abs(indx4-indx8);
    int IB = abs(indx2-indx6);

    int IC = indx4+indx8;
    int ID = indx2+indx6;

    int t, t1 = abs(indx1-indx9), t2 = indx1+indx9;

    if(t1<IA) t1 = IA;
    if(t1<IB) t1 = IB;
    if(t2>IC) t2 = IC;
    if(t2>ID) t2 = ID;

    double val = 0.0, ht;

    for(t=t1; t<=t2; t+=2){
      ht = 0.5*t;
      val += (t+1)*WIGNER_6J(J1, J9, ht, J8, J4, J7, fctsg) \
          *WIGNER_6J(J2, J6, ht, J4, J8, J5, fctsg) \
          *WIGNER_6J(J1, J9, ht, J6, J2, J3, fctsg);
    }

    val *= fctsg->sg[t1];

    if(fctsg->memo){
      inster_9d(fctsg->J9, val, indx1, indx2, indx3, indx4, indx5, \
          indx6, indx7, indx8, indx9);
    }
    
    return val;
}

/*----------------------------------------------------------------------------*/
