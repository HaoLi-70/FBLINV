
#include "IO.h"

/*----------------------------------------------------------------------------*/

    /*######################################################################
 
      revision log:

        19 Nov. 2024
          --- update: output the nx for each pixel.

        30 Otc. 2024
          --- update: moved subroutine FREE_COEFF to FREE.c

        19 Sep. 2024
          --- update: synthesis for multi LOS .
 
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

    size_t nsize;

    char tmp[5];
    nsize = fread(tmp,sizeof(char),4,fa);
    tmp[4] = '\0';

    if(strcmp(tmp, "zero") == 0){
      Input->BC = enum_zeros;
    }else if(strcmp(tmp, "zdri") == 0){
      Input->BC = enum_deri_zeros;
    }else{
      const char *routine_name = "READ_ZEROS";
      Error(enum_error, routine_name, "not a zeros file");
    }

    nsize = fread(&(Input->Lzero),sizeof(int),1,fa);
    nsize = fread(&(Input->Nzero),sizeof(int),1,fa);

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
      nsize = fread(&(Input->Zeros[il][1]),sizeof(double),Input->Nzero,fa);
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
    size_t nsize;

    if(!FILE_EXIST(Input->Path_Observation)){
      fprintf(stderr, "%s \n",Input->Path_Observation);
      Error(enum_error, "routine_name", "observation doesn't exist!");
    }
    
    FILE *Fa = fopen(Input->Path_Observation, "rb");

    char header[5];
    nsize = fread(header,sizeof(char),4,Fa);
    header[4] = '\0';

    if(strcmp(header, "syth") != 0){
      Error(enum_error, routine_name, "not a observation file");
    }

    int ny, nz, nline, nThom;
    double FOV[4];

    nsize = fread(&ny,sizeof(int),1,Fa);
    nsize = fread(&nz,sizeof(int),1,Fa);

    if(Input->ny!=ny||Input->nz!=nz){
      Error(enum_warning, routine_name, "update the ny and nz");
      fprintf(stderr,"input ny = %d nz = %d \n",Input->ny, Input->nz);
      fprintf(stderr,"observation ny = %d nz = %d \n",ny, nz);
      Input->ny = ny;
      Input->nz = nz;
    }

    nsize = fread(&nline,sizeof(int),1,Fa);
    
    if(Input->Nline!=nline){
      fprintf(stderr,"atom nline = %d observation nline = %d \n", \
          Input->Nline, nline);
      Error(enum_error, routine_name, "nline number mismatch");
    }

    nsize = fread(&nThom,sizeof(int),1,Fa);

    if(Input->NThom!=nThom){
      fprintf(stderr,"input nThom = %d observation nThom = %d \n", \
          Input->NThom, nThom);
      Error(enum_error, routine_name, "nThom number mismatch");
    }

    nsize = fread(FOV,sizeof(double),4,Fa);

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
            nsize = fread(Obs->Data[ipixel], sizeof(double), Obs->num, Fa);
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
        19 Nov. 2024.
      Input parameters:
        Input, a structure with the input configuration.
        Output, output structure.
      Return:
        .
    ######################################################################*/
  
    if(Input->Mode<0 || Input->Mode>3){
      return -1;
    }

    int ipixel, ispec, iatom, iline, itrans, ithom, ilos;
    double dtmp;

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

      fwrite(&(Syn->ngrids),sizeof(int),1,Fa);
      fwrite(&(Syn->npixels),sizeof(int),1,Fa);

      // write the number of viewing directions
      fwrite(&(Input->nlos),sizeof(int),1,Fa);

      // write the viewing directions
      for(ilos=0; ilos<Input->nlos; ilos++){
        fwrite(&(Input->los[ilos]),sizeof(double),1,Fa);
      }

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

      fwrite(&Input->rsqmin, sizeof(double), 1, Fa);
      fwrite(&Input->rsqmax, sizeof(double), 1, Fa);
      fwrite(&Input->rsqint, sizeof(double), 1, Fa);
      fwrite(&Input->dx, sizeof(double), 1, Fa);

      // write the field of view
      fwrite(Input->FOV[0], sizeof(double), 4, Fa);

      fwrite(Input->nx[0], sizeof(int), Input->ny*Input->nz, Fa);

      for(ilos=0; ilos<Input->nlos; ilos++){
        ipixel = 0;
        for(ipixel=0; ipixel<Syn->npixels; ipixel++){
          fwrite(Output->los[ilos].syn[ipixel],sizeof(double),num,Fa);
        }
      }
      
      FREE_VECTOR(zeros, 0, enum_dbl);

      // write number of spectral lines
      fwrite(&(Input->Nspec),sizeof(int),1,Fa);

      // write the profiles if there are
      if(Input->Nspec>0){
      //////////////////////spectrum  
        

        // write the size
        fwrite(&(Input->nys),sizeof(int),1,Fa);
        fwrite(&(Input->nzs),sizeof(int),1,Fa);
        fwrite(&(Syn->npspec),sizeof(int),1,Fa);

        // write the field of view
        fwrite(Input->FOVSPECINDX[0],sizeof(int),4,Fa);

        // write wavelength
        for(ispec=0;ispec<Input->Nspec;ispec++){
          fwrite(&(Input->Spec[ispec].Nl),sizeof(int),1,Fa);
          fwrite(Input->Spec[ispec].Lambda,sizeof(double), \
              Input->Spec[ispec].Nl,Fa);
        }

        for(ilos=0; ilos<Input->nlos; ilos++){
          for(ipixel=0; ipixel<Syn->npspec; ipixel++){
            fwrite(Output->los[ilos].spec[ipixel],sizeof(double), \
                Input->Nl*Input->Nstk,Fa);
          }
        }
      }

    
    }else if(Input->Mode==2){
      // single grid
      char chartmp[5] = "sigg";
      fwrite(chartmp,sizeof(char),4,Fa);

      // write coordinate 
      fwrite(&(Input->Ypos),sizeof(double),1,Fa);
      fwrite(&(Input->Zpos),sizeof(double),1,Fa);
      fwrite(&(Input->Xpos),sizeof(double),1,Fa);
      fwrite(&(Input->Nstk),sizeof(int),1,Fa);

      // write the number of viewing directions
      fwrite(&(Input->nlos),sizeof(int),1,Fa);

      // write the viewing directions
      for(ilos=0; ilos<Input->nlos; ilos++){
        fwrite(&(Input->los[ilos]),sizeof(double),1,Fa);
      }

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
      for(ilos=0; ilos<Input->nlos; ilos++){
        fwrite(Output->los[ilos].syn[0],sizeof(double),num,Fa);
      }

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
        for(ilos=0; ilos<Input->nlos; ilos++){
          fwrite(Output->los[ilos].spec[0],sizeof(double), \
              Input->Nl*Input->Nstk*Syn->npspec,Fa);
        }
      }

    }else if(Input->Mode==3){
      // single pixel
      char chartmp[5] = "sigp";
      fwrite(chartmp,sizeof(char),4,Fa);
      
            // write coordinate 
      fwrite(&(Input->Ypos),sizeof(double),1,Fa);
      fwrite(&(Input->Zpos),sizeof(double),1,Fa);
      fwrite(&(Input->Nstk),sizeof(int),1,Fa);

      // write the number of viewing directions
      fwrite(&(Input->nlos),sizeof(int),1,Fa);

      // write the viewing directions
      for(ilos=0; ilos<Input->nlos; ilos++){
        fwrite(&(Input->los[ilos]),sizeof(double),1,Fa);
      }


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

      for(ilos=0; ilos<Input->nlos; ilos++){
        // write integrated polarization 
        fwrite(Output->los[ilos].syn[0],sizeof(double),num,Fa);
      }

      // write number of spectral lines
      fwrite(&(Input->Nspec),sizeof(int),1,Fa);

      if(Syn->npspec>0){
    
        // write wavelength
        for(ispec=0;ispec<Input->Nspec;ispec++){
          fwrite(&(Input->Spec[ispec].Nl),sizeof(int),1,Fa);
          fwrite(Input->Spec[ispec].Lambda,sizeof(double), \
              Input->Spec[ispec].Nl,Fa);
        }
  
        for(ilos=0; ilos<Input->nlos; ilos++){
          // write spectra
          fwrite(Output->los[ilos].spec[0],sizeof(double), \
                Input->Nl*Input->Nstk*Syn->npspec,Fa);
        }
      }

    } 
    
    fclose(Fa);
    
    return 0;
}

/*----------------------------------------------------------------------------*/

extern int WRITE_PARA(STRUCT_INPUT *Input, STRUCT_SYN *Syn, \
    STRUCT_ATOM *Atom, STRUCT_MPI *Mpi){
  
    /*######################################################################
      Purpose:
        Output the parameters.
      Record of revisions:
        19 Nov. 2024.
      Input parameters:
        Input, a structure with the input configuration.
        Syn, a structure with forward synthesis.
        Atmo, a structure with the coronal model.
        Mpi, a structure with the Mpi information.
      Return:
        .
    ######################################################################*/
  

    MPI_File fh;
    MPI_File_open(MPI_COMM_WORLD, Input->Path_Para, \
        MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

    STRUCT_GRID *pgrid;
    STRUCT_LOS *plos;

    int igrid, ilos, iatom;
    int npar = 5+Input->Natom;
    if(Input->Nspec>0){
      npar += 2;
    }
    double h2e = C_H2E;

    size_t header = sizeof(char)*4+sizeof(int)*7+sizeof(double) \
        *(9+Input->Natom)+sizeof(int)*Input->ny*Input->nz;
    size_t offset = Mpi->gindx[0]*(npar*sizeof(double)+sizeof(int));
    size_t offsetlos = Syn->ngrids*(npar*sizeof(double)+sizeof(int));

    if(Mpi->rank==0){

      char chartmp[5] = "para";
      MPI_File_write(fh, chartmp, 4, MPI_CHAR, MPI_STATUS_IGNORE);
      MPI_File_write(fh, &Input->ny, 1, MPI_INT, MPI_STATUS_IGNORE);
      MPI_File_write(fh, &Input->nz, 1, MPI_INT, MPI_STATUS_IGNORE);
      MPI_File_write(fh, &Input->nlos, 1, MPI_INT, MPI_STATUS_IGNORE);
      MPI_File_write(fh, &Syn->ngrids, 1, MPI_INT, MPI_STATUS_IGNORE);
      MPI_File_write(fh, &Syn->npixels, 1, MPI_INT, MPI_STATUS_IGNORE);
      MPI_File_write(fh, &Input->Natom, 1, MPI_INT, MPI_STATUS_IGNORE);
      MPI_File_write(fh, &Input->Nspec, 1, MPI_INT, MPI_STATUS_IGNORE);

      MPI_File_write(fh, &Input->rsqmin, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
      MPI_File_write(fh, &Input->rsqmax, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
      MPI_File_write(fh, &Input->rsqint, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);

      MPI_File_write(fh, &Input->dx, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);

      MPI_File_write(fh, Input->FOV[0], 4, MPI_DOUBLE, MPI_STATUS_IGNORE);

      MPI_File_write(fh, &h2e, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);

      for(iatom=0; iatom<Input->Natom; iatom++){
        MPI_File_write(fh, &(Atom[iatom].Abund), 1, MPI_DOUBLE, \
            MPI_STATUS_IGNORE);
      }

      MPI_File_write(fh, Input->nx[0], Input->ny*Input->nz, MPI_INT, \
          MPI_STATUS_IGNORE);

    }
 
    for(ilos=0; ilos<Input->nlos; ilos++){

      MPI_File_seek(fh, header+offset+offsetlos*ilos, MPI_SEEK_SET);

      for(igrid=0; igrid<Mpi->ngrids; igrid++){
          
        pgrid = Syn->Grids+igrid;
        plos = pgrid->los+ilos;
      
        MPI_File_write(fh, &pgrid->ipixel, 1, \
            MPI_INT, MPI_STATUS_IGNORE);

        MPI_File_write(fh, &plos->Para.B, 1, MPI_DOUBLE, \
            MPI_STATUS_IGNORE);
        MPI_File_write(fh, &plos->Para.ThetaB, 1, MPI_DOUBLE, \
            MPI_STATUS_IGNORE);
        MPI_File_write(fh, &plos->Para.PhiB, 1, MPI_DOUBLE, \
            MPI_STATUS_IGNORE);
        MPI_File_write(fh, &plos->Para.T, 1, MPI_DOUBLE, \
            MPI_STATUS_IGNORE);
        MPI_File_write(fh, &plos->Para.ne, 1, MPI_DOUBLE, \
            MPI_STATUS_IGNORE);

        MPI_File_write(fh, plos->Para.Ion, Input->Natom, \
            MPI_DOUBLE, MPI_STATUS_IGNORE);

        if(Input->Nspec>0){
          MPI_File_write(fh, &plos->Blos, 1, MPI_DOUBLE, \
              MPI_STATUS_IGNORE);
          MPI_File_write(fh, &plos->Vlos, 1, MPI_DOUBLE, \
              MPI_STATUS_IGNORE);
        }

      }
    }

    MPI_File_close(&fh);

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
      fwrite(&indx, sizeof(int), 1, Fa);
    }

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

    size_t nsize;

    FILE *Fa = fopen(Filename, "rb");

    char tmp[5];
    nsize = fread(tmp,sizeof(char),4,Fa);
    tmp[4] = '\0';

    if(strcmp(tmp, "sfbc") != 0){
      const char *routine_name = "READ_COEFF";
      Error(enum_error, routine_name, "not a coefficient file");
    }

    int ipara, in, il, im, indx;
    double coeffr, coeffi;

    for(ipara=0;ipara<4;ipara++){
      nsize = fread(&indx, sizeof(int), 1, Fa);
      Para[ipara].in = (indx==1);
    }

    nsize = fread(&btype, sizeof(int), 1, Fa);
    nsize = fread(&rhotype, sizeof(int), 1, Fa);
    nsize = fread(Rlim, sizeof(double), 1, Fa);

    for(ipara=0;ipara<4;ipara++){
      if(!Para[ipara].in) continue;

      nsize = fread(&Para[ipara].N, sizeof(int), 1, Fa);
      nsize = fread(&Para[ipara].L, sizeof(int), 1, Fa);
      Para[ipara].Coeff = (complex double ***) \
          TENSOR_RHO_CPLX(Para[ipara].N, Para[ipara].L, false);

      for(in=1; in<=Para[ipara].N; in++){
        for(il=0; il<=Para[ipara].L; il++){
          for(im=0; im<=il; im++){
            if(im==0){
              nsize = fread(&coeffr, sizeof(double), 1, Fa);
              Para[ipara].Coeff[in][il][im] = coeffr;
            }else{
              nsize = fread(&coeffr, sizeof(double), 1, Fa);
              nsize = fread(&coeffi, sizeof(double), 1, Fa);
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

extern void WRITE_ATMO(char *Filename, STRUCT_ATMO *Atmo){
  
    /*######################################################################
      Purpose:
        write the coronal model.
      Record of revisions:
        16 Jan. 2024.
      Input parameters:
        Filename, path to the mode atmosphere.
        Atmo, a structure with the coronal model.
    ######################################################################*/

    FILE *fa = fopen(Filename, "wb");

    char tmp[5] = "fbmd";
    fwrite(tmp, sizeof(char), 4, fa);

    fwrite(&(Atmo->nR),sizeof(int),1,fa);
    fwrite(&(Atmo->nTheta),sizeof(int),1,fa);
    fwrite(&(Atmo->nPhi),sizeof(int),1,fa);

    fwrite(Atmo->R,sizeof(double),Atmo->nR,fa);
    fwrite(Atmo->Theta,sizeof(double),Atmo->nTheta,fa);
    fwrite(Atmo->Phi,sizeof(double),Atmo->nPhi,fa);

    int num = Atmo->nR*Atmo->nTheta*Atmo->nPhi;

    fwrite(Atmo->T[0][0],sizeof(double),num,fa);

    fwrite(&(Atmo->rhotype),sizeof(int),1,fa);
    fwrite(Atmo->rho[0][0],sizeof(double),num,fa);

    fwrite(&(Atmo->btype),sizeof(int),1,fa);
    fwrite(Atmo->B1[0][0],sizeof(double),num,fa);
    fwrite(Atmo->B2[0][0],sizeof(double),num,fa);
    fwrite(Atmo->B3[0][0],sizeof(double),num,fa);

    if(Atmo->vtype>=0){
      fwrite(&(Atmo->vtype),sizeof(int),1,fa);
      fwrite(Atmo->V1[0][0],sizeof(double),num,fa);
      fwrite(Atmo->V2[0][0],sizeof(double),num,fa);
      fwrite(Atmo->V3[0][0],sizeof(double),num,fa);
    }

    fclose(fa);
    return;
}

/*----------------------------------------------------------------------------*/