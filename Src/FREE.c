
#include "FREE.h"

/*----------------------------------------------------------------------------*/

    /*######################################################################
    
      revision log:

        19 Nov. 2024
          --- update: free the nx in the input structure.

        30 Otc. 2024
          --- initial comment 
    
    ######################################################################*/

/*----------------------------------------------------------------------------*/

extern int FREE_INPUT(STRUCT_INPUT *Input){

    /*######################################################################
      Purpose:
        free the input structure.
      Record of revisions:
        19 Nov. 2023.
      Input parameters:
        Input, a structure with the input information.
    ######################################################################*/

    int ipara;

    if(Input->Mode < 4){

      FREE_MATRIX(Input->Dkmn[2], -2, -2, enum_cplx);

      FREE_TENSOR_RHO_CPLX(Input->Tkq);

      FREE_MATRIX(Input->Path_Atom, 0, 0, enum_char);

      if(Input->NThom > 0){
        free(Input->Thom);
      }
    }

    if(Input->Mode == 0){

      FREE_MATRIX(Input->Zeros, 0, 1, enum_dbl);

      for(ipara=0;ipara<4;ipara++){
        if(Input->Para[ipara].invt){ 
          FREE_TENSOR_RHO_CPLX(Input->Para[ipara].Coeff);
        }
      }
    }

    free(Input->los);
    free(Input->weight);
    if(Input->Nspec>0 ){
      free(Input->Spec);
    }
    if(Input->Mode<=3){
      FREE_MATRIX(Input->nx, 0, 0, enum_int);
    }
    free(Input);
    
    return 0;

}

/*----------------------------------------------------------------------------*/

extern int FREE_COEFF(STRUCT_SFB_COEFF *Para){

    /*######################################################################
      Purpose:
        free the coefficients.
      Record of revisions:
        20 Aut. 2023.
      Input parameters:
        Para, structure with the coefficients.
    ######################################################################*/

    int ipara;

    for(ipara=0;ipara<4;ipara++){
      if(!Para[ipara].invt) continue;
      FREE_TENSOR_RHO_CPLX(Para[ipara].Coeff);
    }

    free(Para);

    return 0;
}

/*----------------------------------------------------------------------------*/

extern void FREE_MODEL(STRUCT_ATMO *Atmo){
  
    /*######################################################################
      Purpose:
        free the memory of the coronal model.
      Record of revisions:
        16 Jan. 2024.
      Input parameters:
        Atmo, a structure with the coronal model.
    ######################################################################*/
  
    free(Atmo->R);
    free(Atmo->Theta);
    free(Atmo->Phi);
    FREE_TENSOR_DBL(Atmo->B1, 0, 0, 0);
    FREE_TENSOR_DBL(Atmo->B2, 0, 0, 0);
    FREE_TENSOR_DBL(Atmo->B3, 0, 0, 0);
    FREE_TENSOR_DBL(Atmo->rho, 0, 0, 0);
    FREE_TENSOR_DBL(Atmo->T, 0, 0, 0);
    if(Atmo->vtype>=0){
      FREE_TENSOR_DBL(Atmo->V1, 0, 0, 0);
      FREE_TENSOR_DBL(Atmo->V2, 0, 0, 0);
      FREE_TENSOR_DBL(Atmo->V3, 0, 0, 0);
    }
    free(Atmo);
        
    return;
}

/*----------------------------------------------------------------------------*/

extern void FREE_ATOM(STRUCT_ATOM *Atom, int Natom){
    
    /*######################################################################
      Purpose:
        free the memory of atomic structure and collisional rates.
      Record of revisions:
        30 Nov. 2019.
      Input parameters:
        Atom, a structure saved the atomic information.
        Natom, number of atoms;
    ######################################################################*/

    int iatom, itrans, ilevel;

    for(iatom=0;iatom<Natom;iatom++){
      // free transitions
      for (itrans=0; itrans<Atom[iatom].Ntrans; itrans++) {
        FREE_MATRIX_RHO_CPLX(Atom[iatom].TR[itrans].JKQ);
        if(Atom[iatom].TwoLv) continue;
        FREE_TENSOR_DBL(Atom[iatom].TR[itrans].TA, 0, 0, 0);
        FREE_TENSOR_DBL(Atom[iatom].TR[itrans].TS, 0, 0, 0);
        FREE_TENSOR_DBL(Atom[iatom].TR[itrans].RA, 0, 0, 0);
        FREE_TENSOR_DBL(Atom[iatom].TR[itrans].RS, 0, 0, 0);
        FREE_VECTOR(Atom[iatom].TR[itrans].TE, 0, enum_dbl);
      }
      free(Atom[iatom].TR);
   
      // free levels
      for (ilevel=0; ilevel<Atom[iatom].Nlevel; ilevel++) {
        FREE_MATRIX_RHO_CPLX(Atom[iatom].LV[ilevel].Rho);
      }
      free(Atom[iatom].LV);    

      // free collision
      if(Atom[iatom].collision){
        FREE_VECTOR(Atom[iatom].col->T, 0, enum_dbl);
        FREE_TENSOR_DBL(Atom[iatom].col->Strength, 0, 0, 0);
        FREE_MATRIX(Atom[iatom].col->Rates, 0, 0, enum_dbl);
        free(Atom[iatom].col);
      }

      free(Atom[iatom].iout);
      free(Atom[iatom].eqindx);
      // free ionization
      FREE_MATRIX(Atom[iatom].Ion.frac, 0, 0, enum_dbl);
    }

    free(Atom);
    
    return;
}

/*----------------------------------------------------------------------------*/

extern int FREE_GRIDS(STRUCT_INPUT *Input, STRUCT_SYN *Syn, \
    STRUCT_MPI *Mpi){

    /*######################################################################
      Purpose:
        free the memory of grids.
      Record of revisions:
        12 Aug. 2023
      Input parameters:
        Input, a structure with the input information.
        Syn, a structure with forward synthesis.
        Mpi, a structure with the Mpi information.
    ######################################################################*/

    int igrid, iatom, ilos;
    STRUCT_GRID *pgrid;
    STRUCT_LOS *plos;

    for(igrid=0; igrid<Mpi->ngrids; igrid++){
      pgrid = Syn->Grids+igrid;

      if(Input->NThom > 0){
        free(pgrid->Kr);
        free(pgrid->Kt);
      }

      if(Input->Nline > 0){
        
        for(iatom=0; iatom<Input->Natom; iatom++){
          free(pgrid->Jkq[iatom].J00);     
          FREE_MATRIX(pgrid->Jkq[iatom].J2q, 0, -Input->JQmax, enum_dbl);
        }
       
        free(pgrid->Jkq);
      }

      for(ilos=0;ilos<Input->nlos;ilos++){
        plos = pgrid->los+ilos;
        if(Input->NThom > 0){
          free(plos->Thom);
        }
        if(Input->Nline>0){
          free(plos->Line);
          free(plos->Para.Ion);
        }
        if(Input->Mode==0){
          free(plos->Line_unperturb);
          free(plos->Thom_unperturb);
          free(plos->Para_unperturb->Ion);
          free(plos->Para_unperturb);
        }
      }

      free(pgrid->los);

      FREE_MATRIX(pgrid->T2Q, 0, -2, enum_cplx);

      if(Input->Mode == 0){
        FREE_VECTOR(pgrid->Phiarray, 0, enum_cplx);
        FREE_MATRIX(pgrid->sBessel, 1, 0, enum_dbl);
        FREE_MATRIX_TRI(pgrid->Legendre, enum_dbl);
        if(Input->Bpotential){
          FREE_MATRIX(pgrid->sBesselD, 1, 0, enum_dbl);
          FREE_MATRIX_TRI(pgrid->LegendreD, enum_dbl);
        }
      }
    }

    free(Syn->Grids);
    free(Syn);

    return 0;
}

/*----------------------------------------------------------------------------*/

extern void FREE_OUTPUT(STRUCT_OUT *Output, STRUCT_INPUT *Input){
    
    /*######################################################################
      Purpose:
        free the output structure.
      Record of revisions:
        30 Otc. 2024.
      Input parameters:
        Output, the output structure.
    ######################################################################*/

    int ilos;

    if(Input->Mode==0){
      FREE_VECTOR(Output->norm,0,enum_dbl);
      FREE_VECTOR(Output->R,0,enum_dbl);  
    }

    if(Input->Mode<=3){

      for(ilos=0;ilos<Input->nlos;ilos++){
        FREE_MATRIX(Output->los[ilos].syn, 0, 0, enum_dbl);
      }
      
      FREE_MATRIX(Output->synloc, 0, 0, enum_dbl);

      if(Input->Nspec>0){

        for(ilos=0;ilos<Input->nlos;ilos++){
          FREE_MATRIX(Output->los[ilos].spec, 0, 0, enum_dbl);
        }  
        FREE_MATRIX(Output->specloc, 0, 0, enum_dbl);
      }

      free(Output->los);     
    }
      
    free(Output);

    return;
}

/*----------------------------------------------------------------------------*/
