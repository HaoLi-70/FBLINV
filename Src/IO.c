
#include "IO.h"

/*----------------------------------------------------------------------------*/

    /*######################################################################
 
      revision log:
        8 Sept. 2021.
 
    ######################################################################*/

/*----------------------------------------------------------------------------*/

extern int READ_ZEROS(char *Filename, STRUCT_INPUT *Input, STRUCT_MPI *Mpi){

    /*######################################################################
      Purpose:
        Read the zeros of Spherical Bessel function of the zeros of the 
            derivative of Spherical Bessel function.
      Record of revisions:
        10 aug. 2023.
      Input parameters:
        Filename[], the input file.
      Output parameters:
        Input, a structure saved the input information.
        Mpi, a structure saved the Mpi information.
    ######################################################################*/

    FILE *fa = fopen(Filename, "rb");

    char tmp[5];
    fread(tmp,sizeof(char),4,fa);
    tmp[4] = '\0';

    if(strcmp(tmp, "zero") == 0){
      Input->BC = enum_zeros;
    }else if(strcmp(tmp, "zdri") == 0){
      Input->BC = enum_deri_zeros;
    }else{
      const char *routine_name = "READ_ZEROS";
      Error(enum_error, routine_name, "not a zeros file");
    }

    fread(&(Input->Lzero),sizeof(int),1,fa);
    fread(&(Input->Nzero),sizeof(int),1,fa);

    if(Input->Verbose && Mpi->rank == 0){
      if(Input->BC == enum_zeros){
        fprintf(stderr, " read the zeros of Spherical Bessel function\n");
      }else{
        fprintf(stderr, " read the zeros of the derivative of Spherical"\
            " Bessel function\n");
      }
      fprintf(stderr, " L = %d, N = %d\n",Input->Lzero,Input->Nzero);
    }

    Input->Zeros = (double **)MATRIX(0, Input->Lzero, 1, \
        Input->Nzero, enum_dbl, false);

    int il;

    for(il=0;il<=Input->Lzero;il++){
      fread(&(Input->Zeros[il][1]),sizeof(double),Input->Nzero,fa);
    }

    fclose(fa);

    return 0;

}

/*----------------------------------------------------------------------------*/

extern int WRITE_ZEROS(char *Filename, int Zero_L, int Zero_N, \
    enum boundary_condition type){
    
    /*######################################################################
      Purpose:
        write positive zeros of spherical Bessel function or of its 
            derivative into a file.
      Record of revisions:
        19 Aug. 2023.
      Input parameters:
        filename, the file nema.
        Zero_L, The max order L of the Bessel L.
        Zero_N, The max number of the zeros.
        type, If the type = enum_deri_zeros, output the zeros of the 
            derivative, else output the zeros of the function.
      Notes:
        X01 is set to 0 if type = enum_deri_zeros.
    ######################################################################*/

    FILE *Fa = fopen(Filename, "wb");

    double **Zeros = (double **)MATRIX(0, Zero_L, 1, Zero_N, enum_dbl, \
        false);
    
    Compute_SB_Zeros(Zeros, Zero_L, Zero_N, type);
    int i;

    if(type==enum_zeros){
      // write zeros of the spherical Bessel functions
      char io[5] = "zero";
      fwrite(io, sizeof(char), 4, Fa);

    }else{
      // write zeros of the derivatives of the spherical Bessel functions
      char io[5] = "zdri";
      fwrite(io, sizeof(char), 4, Fa);
    }

    // order L
    fwrite(&Zero_L, sizeof(int), 1, Fa);

    // number
    fwrite(&Zero_N, sizeof(int), 1, Fa);

    for (i=0; i<=Zero_L; i++) {
      fwrite(&Zeros[i][1], sizeof(double), Zero_N, Fa);
    }

    fclose(Fa);
    FREE_MATRIX(Zeros, 0, 1, enum_dbl);

    return 0;
}

/*----------------------------------------------------------------------------*/

extern int READ_OBSERVATION(STRUCT_INPUT *Input, STRUCT_OBSERVATION *Obs){
  
    /*######################################################################
      Purpose:
        Read the observation.
      Record of revisions:
        12 Aug. 2023.
      Input parameters:
        Input, a structure with the input configuration.
        Obs, structure with the observation data.
      Return:
        .
    ######################################################################*/
  
    if(Input->Mode!=0) return -1;

    const char *routine_name = "READ_OBSERVATION";

    if(!FILE_EXIST(Input->Path_Observation)){
      fprintf(stderr, "%s \n",Input->Path_Observation);
      Error(enum_error, "routine_name", "observation doesn't exist!");
    }
    
    FILE *Fa = fopen(Input->Path_Observation, "rb");

    char header[5];
    fread(header,sizeof(char),4,Fa);
    header[4] = '\0';

    if(strcmp(header, "syth") != 0){
      Error(enum_error, routine_name, "not a observation file");
    }

    int ny, nz, nline, nThom;
    double FOV[4];

    fread(&ny,sizeof(int),1,Fa);
    fread(&nz,sizeof(int),1,Fa);

    if(Input->ny!=ny||Input->nz!=nz){
      Error(enum_warning, routine_name, "update the ny and nz");
      fprintf(stderr,"input ny = %d nz = %d \n",Input->ny, Input->nz);
      fprintf(stderr,"observation ny = %d nz = %d \n",ny, nz);
      Input->ny = ny;
      Input->nz = nz;
    }

    fread(&nline,sizeof(int),1,Fa);
    
    if(Input->Nline!=nline){
      fprintf(stderr,"atom nline = %d observation nline = %d \n", \
          Input->Nline, nline);
      Error(enum_error, routine_name, "nline number mismatch");
    }

    fread(&nThom,sizeof(int),1,Fa);

    if(Input->NThom!=nThom){
      fprintf(stderr,"input nThom = %d observation nThom = %d \n", \
          Input->NThom, nThom);
      Error(enum_error, routine_name, "nThom number mismatch");
    }

    fread(FOV,sizeof(double),4,Fa);

    if(Input->FOV[0][0]!=FOV[0]||Input->FOV[0][1]!=FOV[1] \
        ||Input->FOV[1][0]!=FOV[2]||Input->FOV[1][1]!=FOV[3]){
      Error(enum_warning, routine_name, "update the FOV");
      fprintf(stderr,"input yrange = %e %e , zrange = %e %e  \n", \
          Input->FOV[0][0], Input->FOV[0][1], Input->FOV[1][0], \
          Input->FOV[1][1]);
      fprintf(stderr,"observation yrange = %e %e , zrange = %e %e  \n", \
          FOV[0], FOV[1], FOV[2], FOV[3]);
      fprintf(stderr,"observation ny = %d nz = %d \n",ny, nz);
      Input->FOV[0][0] = FOV[0];
      Input->FOV[0][1] = FOV[1];
      Input->FOV[1][0] = FOV[2];
      Input->FOV[1][1] = FOV[3];
    }

    Input->dy = (Input->FOV[0][1]-Input->FOV[0][0])/(Input->ny-1);
    Input->dz = (Input->FOV[1][1]-Input->FOV[1][0])/(Input->nz-1);

    int iy, iz, ipixel = 0;
    double Y, Z, rsq0, tmp;

    Obs->npixels = 0;

    for(iy=0; iy<Input->ny; iy++){
      Y = Input->FOV[0][0]+iy*Input->dy;
      for(iz=0; iz<Input->nz; iz++){
        Z = Input->FOV[1][0]+iz*Input->dz;
        rsq0 = Y*Y+Z*Z;
        if(rsq0<Input->rsqmin || rsq0>Input->rsqmax) continue;
        tmp = sqrt(Input->rsqint-rsq0);
        if(tmp<Input->dx) continue;
        Obs->npixels++;      
      }
    }

    Obs->num = 3*Input->Nline+2*Input->NThom;

    Obs->Data = (double **)MATRIX(0, Obs->npixels-1, 0, Obs->num, \
        enum_dbl, true);

    for(iy=0; iy<Input->ny; iy++){
      Y = Input->FOV[0][0]+iy*Input->dy;
      for(iz=0; iz<Input->nz; iz++){
        Z = Input->FOV[1][0]+iz*Input->dz;
        rsq0 = Y*Y+Z*Z;
        if(rsq0>=Input->rsqmin && rsq0<=Input->rsqmax){
          tmp = sqrt(Input->rsqint-rsq0);
          if(tmp>=Input->dx){
            fread(Obs->Data[ipixel], sizeof(double), Obs->num, Fa);
            ipixel++;
          }else{
            fseek(Fa, sizeof(double)*Obs->num, SEEK_CUR);
          }
        }else{
          fseek(Fa, sizeof(double)*Obs->num, SEEK_CUR);
        }
      }
    }

    fclose(Fa);
    
    return 0;
}

/*----------------------------------------------------------------------------*/

extern int WRITE_SYNTHESIS(char *Filename, STRUCT_INPUT *Input, \
    STRUCT_ATOM *Atom, STRUCT_OUT *Output, STRUCT_SYN *Syn){
  
    /*######################################################################
      Purpose:
        Output the synthetic data.
      Record of revisions:
        12 Aug. 2023.
      Input parameters:
        Input, a structure with the input configuration.
        Output, output structure.
      Return:
        .
    ######################################################################*/
  
    if(Input->Mode<0 || Input->Mode>3){
      return -1;
    }

    int iy, iz, ipixel = 0, ispec, iatom, iline, itrans, ithom;
    double Y, Z, rsq0, dtmp;

    int num = Input->Nstk*Input->Nline+2*Input->NThom;
    FILE *Fa = fopen(Filename, "wb");

    if(Input->Mode==0 || Input->Mode==1){

      double *zeros = (double *)VECTOR(0, num-1, enum_dbl, true);
      // forward
      char chartmp[5] = "syth";
      fwrite(chartmp,sizeof(char),4,Fa);

      // write the size
      fwrite(&(Input->ny),sizeof(int),1,Fa);
      fwrite(&(Input->nz),sizeof(int),1,Fa);
      fwrite(&(Input->Nstk),sizeof(int),1,Fa);

      // write number of lines and the corresponding wavelength
      fwrite(&(Input->Nline),sizeof(int),1,Fa);
      if(Input->Nline>0){
        for(iatom = 0; iatom <Input->Natom; iatom++){
          for(iline=0; iline<Atom[iatom].Nline; iline++){
            itrans = Atom[iatom].iout[iline];
            dtmp = Atom[iatom].TR[itrans].lambda*1e10;

            fwrite(&dtmp,sizeof(double),1,Fa);
          }
        }
      }

      // write Thomson scattering grid and the corresponding wavelength
      fwrite(&(Input->NThom),sizeof(int),1,Fa);
      if(Input->NThom>0){
        for(ithom=0;ithom<Input->NThom;ithom++){
          dtmp = Input->Thom[ithom].lambda*1e10;
          fwrite(&dtmp,sizeof(double),1,Fa);
        }
      }

      // write the field of view
      fwrite(Input->FOV[0],sizeof(double),4,Fa);

      for(iy=0; iy<Input->ny; iy++){
        Y = Input->FOV[0][0]+iy*Input->dy;
        for(iz=0; iz<Input->nz; iz++){
          Z = Input->FOV[1][0]+iz*Input->dz;
          rsq0 = Y*Y+Z*Z;
          if(rsq0>=Input->rsqmin && rsq0<=Input->rsqmax){
            dtmp = sqrt(Input->rsqint-rsq0);
            if(dtmp>=Input->dx){
              fwrite(Output->syntot[ipixel],sizeof(double),num,Fa);
              ipixel++;
            }else{
              fwrite(zeros,sizeof(double),num,Fa);
            }
          }else{
            fwrite(zeros,sizeof(double),num,Fa);
          }
        }
      }

      // write the profiles if there are
      if(Input->Nspec>0){
      //////////////////////spectrum  
        
        FREE_VECTOR(zeros, 0, enum_dbl);
        zeros = (double *)VECTOR(0, Input->Nl*Input->Nstk-1, enum_dbl, true);

        // write number of spectral lines
        fwrite(&(Input->Nspec),sizeof(int),1,Fa);

        // write the size
        fwrite(&(Input->nys),sizeof(int),1,Fa);
        fwrite(&(Input->nzs),sizeof(int),1,Fa);
        // write the field of view
        fwrite(Input->FOVSPEC[0],sizeof(double),4,Fa);

        // write wavelength
        for(ispec=0;ispec<Input->Nspec;ispec++){
          fwrite(&(Input->Spec[ispec].Nl),sizeof(int),1,Fa);
          fwrite(Input->Spec[ispec].Lambda,sizeof(double), \
              Input->Spec[ispec].Nl,Fa);
        }

        ipixel = 0;
        for(iy=0; iy<Input->ny; iy++){
          Y = Input->FOV[0][0]+iy*Input->dy;
          for(iz=0; iz<Input->nz; iz++){
            Z = Input->FOV[1][0]+iz*Input->dz;
            rsq0 = Y*Y+Z*Z;
            if(rsq0>=Input->rsqmin && rsq0<=Input->rsqmax){
              dtmp = sqrt(Input->rsqint-rsq0);
              if(dtmp>=Input->dx && Y>Input->FOVSPEC[0][0] \
                   && Y<Input->FOVSPEC[0][1] && Z>Input->FOVSPEC[1][0] \
                   && Z<Input->FOVSPEC[1][1]){
                fwrite(Output->spectot[ipixel],sizeof(double), \
                    Input->Nl*Input->Nstk,Fa);
                ipixel++;
              }else{
                fwrite(zeros,sizeof(double),Input->Nl*Input->Nstk,Fa);
              }
            }else{
              fwrite(zeros,sizeof(double),Input->Nl*Input->Nstk,Fa);
            }
          }
        }
      }

      FREE_VECTOR(zeros, 0, enum_dbl);
    
    }else if(Input->Mode==2){
      // single grid
      char chartmp[5] = "sigg";
      fwrite(chartmp,sizeof(char),4,Fa);

      // write coordinate 
      fwrite(&(Input->Ypos),sizeof(double),1,Fa);
      fwrite(&(Input->Zpos),sizeof(double),1,Fa);
      fwrite(&(Input->Xpos),sizeof(double),1,Fa);
      fwrite(&(Input->Nstk),sizeof(int),1,Fa);

      // write number of lines and the corresponding wavelength
      fwrite(&(Input->Nline),sizeof(int),1,Fa);
      if(Input->Nline>0){
        for(iatom = 0; iatom <Input->Natom; iatom++){
          for(iline=0; iline<Atom[iatom].Nline; iline++){
            itrans = Atom[iatom].iout[iline];
            dtmp = Atom[iatom].TR[itrans].lambda*1e10;
            fwrite(&dtmp,sizeof(double),1,Fa);
          }
        }
      }

      // write Thomson scattering grid and the corresponding wavelength
      fwrite(&(Input->NThom),sizeof(int),1,Fa);
      if(Input->NThom>0){
        for(ithom=0;ithom<Input->NThom;ithom++){
          dtmp = Input->Thom[ithom].lambda*1e10;
          fwrite(&dtmp,sizeof(double),1,Fa);
        }
      }

      // write integrated polarization 
      fwrite(Output->syntot[0],sizeof(double),num,Fa);

      // write number of spectral lines
      fwrite(&(Input->Nspec),sizeof(int),1,Fa);

      if(Syn->npspec>0){

        // write wavelength
        for(ispec=0;ispec<Input->Nspec;ispec++){
          fwrite(&(Input->Spec[ispec].Nl),sizeof(int),1,Fa);
          fwrite(Input->Spec[ispec].Lambda,sizeof(double), \
              Input->Spec[ispec].Nl,Fa);
        }
  
        // write spectra
        fwrite(Output->spectot[0],sizeof(double), \
            Input->Nl*Input->Nstk*Syn->npspec,Fa);
      }

    }else if(Input->Mode==3){
      // single pixel
      char chartmp[5] = "sigp";
      fwrite(chartmp,sizeof(char),4,Fa);
      
            // write coordinate 
      fwrite(&(Input->Ypos),sizeof(double),1,Fa);
      fwrite(&(Input->Zpos),sizeof(double),1,Fa);
      fwrite(&(Input->Nstk),sizeof(int),1,Fa);

      // write number of lines and the corresponding wavelength
      fwrite(&(Input->Nline),sizeof(int),1,Fa);
      if(Input->Nline>0){
        for(iatom = 0; iatom <Input->Natom; iatom++){
          for(iline=0; iline<Atom[iatom].Nline; iline++){
            itrans = Atom[iatom].iout[iline];
            dtmp = Atom[iatom].TR[itrans].lambda*1e10;
            fwrite(&dtmp,sizeof(double),1,Fa);
          }
        }
      }

      // write Thomson scattering grid and the corresponding wavelength
      fwrite(&(Input->NThom),sizeof(int),1,Fa);
      if(Input->NThom>0){
        for(ithom=0;ithom<Input->NThom;ithom++){
          dtmp = Input->Thom[ithom].lambda*1e10;
          fwrite(&dtmp,sizeof(double),1,Fa);
        }
      }

      // write integrated polarization 
      fwrite(Output->synloc[0],sizeof(double),num,Fa);

      // write number of spectral lines
      fwrite(&(Input->Nspec),sizeof(int),1,Fa);

      if(Syn->npspec>0){
    
        // write wavelength
        for(ispec=0;ispec<Input->Nspec;ispec++){
          fwrite(&(Input->Spec[ispec].Nl),sizeof(int),1,Fa);
          fwrite(Input->Spec[ispec].Lambda,sizeof(double), \
              Input->Spec[ispec].Nl,Fa);
        }
  
        // write spectra
        fwrite(Output->spectot[0],sizeof(double), \
              Input->Nl*Input->Nstk*Syn->npspec,Fa);
      }

    } 
    
    fclose(Fa);
    
    return 0;
}

/*----------------------------------------------------------------------------*/

extern void WRITE_COEFF(char *Filename, STRUCT_SFB_COEFF *Para, int btype, \
    int rhotype, double Rlim){
    
    /*######################################################################
      Purpose:
        print the spherical Fourier Bessl (or only spherical harmonical)
            coefficients.
      Record of revisions:
        20 Aut. 2023.
      Input parameters:
        Filename, name of the output file.
        Para, structure with the coefficients.
    ######################################################################*/

    FILE *Fa = fopen(Filename, "wb");
    int ipara, in, il, im, indx;
    double coeffr, coeffi;

    char io[5] = "sfbc";
    fwrite(io, sizeof(char), 4, Fa);

    for(ipara=0;ipara<4;ipara++){
      if(Para[ipara].invt){
        indx = 1;
      }else{
        indx = 0;
      }
//      fprintf(stderr,"aa %d %d %d \n ",ipara,indx,Para[ipara].invt);
      fwrite(&indx, sizeof(int), 1, Fa);
    }
//      fprintf(stderr,"bb \n ");

    fwrite(&btype, sizeof(int), 1, Fa);
    fwrite(&rhotype, sizeof(int), 1, Fa);
    fwrite(&Rlim, sizeof(double), 1, Fa);

    for(ipara=0;ipara<4;ipara++){
      if(!Para[ipara].invt) continue;
      fwrite(&Para[ipara].N, sizeof(int), 1, Fa);
      fwrite(&Para[ipara].L, sizeof(int), 1, Fa);

      for(in=1; in<=Para[ipara].N; in++){
        for(il=0; il<=Para[ipara].L; il++){
          for(im=0; im<=il; im++){
            if(im==0){
              coeffr = creal(Para[ipara].Coeff[in][il][im]);
              fwrite(&coeffr, sizeof(double), 1, Fa);
            }else{
              coeffr = creal(Para[ipara].Coeff[in][il][im]);
              coeffi = cimag(Para[ipara].Coeff[in][il][im]);
              fwrite(&coeffr, sizeof(double), 1, Fa);
              fwrite(&coeffi, sizeof(double), 1, Fa);
            }
          }
        }
      }
    }
    
    fclose(Fa);
    return;
}

/*----------------------------------------------------------------------------*/

extern void READ_COEFF(char *Filename, STRUCT_SFB_COEFF *Para, double *Rlim){
    
    /*######################################################################
      Purpose:
        read the spherical Fourier Bessl (or only spherical harmonical)
            coefficients.
      Record of revisions:
        20 Aut. 2023.
      Input parameters:
        Filename, name of the output file.
      Output parameters:
        Para, structure with the coefficients.
    ######################################################################*/

    int btype, rhotype;

    if(!FILE_EXIST(Filename)){
      fprintf(stderr, "%s \n", Filename);
      Error(enum_error, "READ_COEFF", "coefficients file doesn't exist!");
    }

    FILE *Fa = fopen(Filename, "rb");

    char tmp[5];
    fread(tmp,sizeof(char),4,Fa);
    tmp[4] = '\0';

    if(strcmp(tmp, "sfbc") != 0){
      const char *routine_name = "READ_COEFF";
      Error(enum_error, routine_name, "not a coefficient file");
    }

    int ipara, in, il, im, indx;
    double coeffr, coeffi;

    for(ipara=0;ipara<4;ipara++){
      fread(&indx, sizeof(int), 1, Fa);
      Para[ipara].in = (indx==1);
    }

    fread(&btype, sizeof(int), 1, Fa);
    fread(&rhotype, sizeof(int), 1, Fa);
    fread(Rlim, sizeof(double), 1, Fa);

    for(ipara=0;ipara<4;ipara++){
      if(!Para[ipara].in) continue;

      fread(&Para[ipara].N, sizeof(int), 1, Fa);
      fread(&Para[ipara].L, sizeof(int), 1, Fa);
      Para[ipara].Coeff = (complex double ***) \
          TENSOR_RHO_CPLX(Para[ipara].N, Para[ipara].L, false);

      for(in=1; in<=Para[ipara].N; in++){
        for(il=0; il<=Para[ipara].L; il++){
          for(im=0; im<=il; im++){
            if(im==0){
              fread(&coeffr, sizeof(double), 1, Fa);
              Para[ipara].Coeff[in][il][im] = coeffr;
            }else{
              fread(&coeffr, sizeof(double), 1, Fa);
              fread(&coeffi, sizeof(double), 1, Fa);
              Para[ipara].Coeff[in][il][im] = coeffr+coeffi*I;
            }
          }
        }
      }
    }
    
    fclose(Fa);
    return;
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
