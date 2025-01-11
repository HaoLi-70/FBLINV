
#include "COEFF.h"

/*----------------------------------------------------------------------------*/

    /*######################################################################
    
      revision log:

        30 Otc. 2024
          --- initial comment 
    
    ######################################################################*/

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
  
    int indx;
    *icoeff = 0;
    for(indx=0; indx<ipara; indx++){
      *icoeff += Input->Para[indx].NLsquare;
    }
    
    *icoeff += (in-1)*Input->Para[ipara].Lsquare+il*il+im*2;
    
    if(real > 0 && im > 0){
      *icoeff = *icoeff-1;
    }
    
    return;
}

/*----------------------------------------------------------------------------*/

extern void array2coeff(double *qk,  STRUCT_INPUT *Input){
  
    /*######################################################################
      Purpose:
        convert qk array to SFB coefficients structure.
      Record of revisions:
        10 Sept. 2021.
      Input parameters:
        qk, a structure saved input configuration.
        Input, a structure saved the inversion information.
      Output parameters:
        Input, a structure saved the inversion information.
    ######################################################################*/
  
    int in, il, im, ipara;
    int ii = 0;

    for(ipara=0; ipara<4; ipara++){
      if(!Input->Para[ipara].invt) continue;
      for(in=1; in<=Input->Para[ipara].N; in++){
        for(il=0; il<=Input->Para[ipara].L; il++){
          for(im=0; im<=il; im++){
            if(im==0){
              Input->Para[ipara].Coeff[in][il][im] = qk[ii];
              ii++;
            }else{
              Input->Para[ipara].Coeff[in][il][im] = qk[ii]+qk[ii+1]*I;
              ii+=2;
            }
          }
        }
      }
    }
    
    return;
}

/*----------------------------------------------------------------------------*/

extern void coeff2array(double *qk,  STRUCT_INPUT *Input){
  
    /*######################################################################
      Purpose:
        convert SFB coefficients structure to qk array.
      Record of revisions:
        10 Sept. 2021.
      Input parameters:
        Input, a structure saved the inversion information.
      Output parameters:
        qk, a structure saved input configuration.
    ######################################################################*/
  
    int in, il, im, ipara;
    int ii = 0;
    
    for(ipara=0; ipara<4; ipara++){
      if(!Input->Para[ipara].invt) continue;
      for(in=1; in<=Input->Para[ipara].N; in++){
        for(il=0; il<=Input->Para[ipara].L; il++){
          for(im=0; im<=il; im++){
            if(im==0){
              qk[ii] = creal(Input->Para[ipara].Coeff[in][il][im]);
              ii++;
            }else{
              qk[ii] = creal(Input->Para[ipara].Coeff[in][il][im]);
              qk[ii+1] = cimag(Input->Para[ipara].Coeff[in][il][im]);
              ii+=2;
            }
          }
        }
      }
    }
    
    return;
}

/*----------------------------------------------------------------------------*/
