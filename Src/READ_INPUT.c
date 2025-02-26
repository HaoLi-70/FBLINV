
#include "READ_INPUT.h"

/*----------------------------------------------------------------------------*/

    /*######################################################################
    
      revision log:

        19 Nov. 2024
          --- update: add **nx to the input structure.
                      remove the parameter file if exist.

        30 Otc. 2024
          --- update: A keywords to output the model parameters at each 
                      grids.
                      moved subroutine FREE_INPUT to FREE.c

        27 Dec. 2024
          --- update: Capable of inputting the magnetic field for the 
                      grid mode.  
    
    ######################################################################*/

/*----------------------------------------------------------------------------*/

extern int Keywords_Conversion(STRUCT_KEYS Keywords[], \
    STRUCT_INPUT *Input, STRUCT_MPI *Mpi){
  
    /*######################################################################
      Purpose:
        convert the keywords to inversion configuration.
      Record of revisions:
        27 Dec. 2024.
      Input parameters:
        Keywords, a structure saved input configuration.
      Output parameters:
        Input, a structure with the input information.
        Mpi, a structure with the Mpi information.
    ######################################################################*/
    
    const char *routine_name = "Keywords_Conversion";
    
    int indx, nread, itmp;
    char inv[4][Key_Length], parameter[Key_Length], *token, **Path_tmp;
    double tmp[20];
    STRUCT_SPEC *spe_tmp;   

    indx = 0;
    Input->Verbose = atoi(Keywords[indx].line);
    if(Input->Verbose>=1 && Mpi->rank == 0){
      fprintf(stderr, "\n verbose level : %d ", Input->Verbose);
    }
    
    indx = 1;
    String_to_Upper(Keywords[indx].line);
   
    if(strcmp(Keywords[indx].line, "INVERSION") == 0){
      Input->Mode = 0;
    }else if(strcmp(Keywords[indx].line, "FORWARD") == 0){
      Input->Mode = 1;
    }else if(strcmp(Keywords[indx].line, "SINGLE_GRID") == 0){
      Input->Mode = 2;
    }else if(strcmp(Keywords[indx].line, "SINGLE_PIXEL") == 0){
      Input->Mode = 3;
    }else if(strcmp(Keywords[indx].line, "DECOMPOSITION") == 0){
      Input->Mode = 4;
    }else if(strcmp(Keywords[indx].line, "RECONSTRUCTION") == 0){
      Input->Mode = 5;
    }else{
      Error(enum_error, routine_name, "mode error!");
    }

    if(Input->Verbose >= 1 && Mpi->rank == 0){
      if(Input->Mode == 0){
        fprintf(stderr, "\n mode : inversion ");
      }else if(Input->Mode == 1){
        fprintf(stderr, "\n mode : forward synthesis ");
      }else if(Input->Mode == 2){
        fprintf(stderr, "\n mode : single grid calculation ");
      }else if(Input->Mode == 3){
        fprintf(stderr, "\n mode : single pixel calculation ");
      }else if(Input->Mode == 4){
        fprintf(stderr, "\n mode : decomposition ");
      }else if(Input->Mode == 5){
        fprintf(stderr, "\n mode : reconstruction ");
      }
    }

    Path_tmp = Input->Path_Atom;
    Input->Path_Atom = (char **)MATRIX(0, Input->Natom-1, 0, \
        Max_Line_Length-1, enum_char, true);
    
    for(itmp=0; itmp<Input->Natom; itmp++){
      String_Copy(Input->Path_Atom[itmp], Path_tmp[itmp], \
          Max_Line_Length-1, true);
      if(Input->Verbose >= 1 && Mpi->rank == 0){
        fprintf(stderr, "\n Path to the Atom : %d %s", itmp, \
            Input->Path_Atom[itmp]);
      }
      if(!FILE_EXIST(Input->Path_Atom[itmp])){
        if(Mpi->rank == 0){
          Error(enum_error, "routine_name", 
            "atomic file doesn't exist!");
        }
        ABORT();
      }
    }
    FREE_MATRIX(Path_tmp,0,0,enum_char);

    indx = 3;
    if(Input->Mode == 0||Input->Mode == 4||Input->Mode == 5){
      String_Copy(Input->Path_Zero, Keywords[indx].line, \
          strlen(Keywords[indx].line), true);
      if(!FILE_EXIST(Input->Path_Atom[itmp])){
        Error(enum_error, "routine_name", "zeros file doesn't exist!");
      }
      if(Input->Verbose >= 1 && Mpi->rank == 0){
        fprintf(stderr, "\n path_zero : %s ", Input->Path_Zero);
      }
    }

    indx = 4;
    nread = sscanf(Keywords[indx].line,"%d %d %d %d ", \
        &(Input->Para[0].N), &(Input->Para[1].N), \
        &(Input->Para[2].N), &(Input->Para[3].N));
    if(Input->Verbose >= 1 && Mpi->rank == 0){
      fprintf(stderr, "\n Bessel_function_N : %d %d %d %d", \
          Input->Para[0].N, Input->Para[1].N, \
          Input->Para[2].N, Input->Para[3].N);
    }
    
    indx = 5;
    nread = sscanf(Keywords[indx].line,"%d %d %d %d ", \
        &(Input->Para[0].L), &(Input->Para[1].L), \
        &(Input->Para[2].L), &(Input->Para[3].L));
    if(Input->Verbose >= 1 && Mpi->rank == 0){
      fprintf(stderr, "\n harmonic_L : %d %d %d %d", \
          Input->Para[0].L, Input->Para[1].L, \
          Input->Para[2].L, Input->Para[3].L);
    }
    
    indx = 6;
    nread = sscanf(Keywords[indx].line,"%lf %lf %lf %lf ", \
        &(Input->Para[0].sparsity), &(Input->Para[1].sparsity), \
        &(Input->Para[2].sparsity), &(Input->Para[3].sparsity));
    if(Input->Verbose >= 1 && Mpi->rank == 0){
      fprintf(stderr, "\n Sparsity : %f %f %f %f", \
          Input->Para[0].sparsity, Input->Para[1].sparsity, \
          Input->Para[2].sparsity, Input->Para[3].sparsity);
    }

    indx = 7;
    String_to_Upper(Keywords[indx].line);
    if(strcmp(Keywords[indx].line, "YES") == 0){
      Input->Saturated= true;
    }else{
      Input->Saturated= false;
    }
    if(Input->Verbose >= 1 && Mpi->rank == 0){
      if(Input->Saturated){
        fprintf(stderr, "\n Redue the SEEs : YES ");
      }else{
        fprintf(stderr, "\n Redue the SEEs : No ");
      }
    }
    
    indx = 8;
    Input->NThom = 0;
    token = strtok(Keywords[indx].line, ",");
    while(token != NULL){
      nread = sscanf(token,"%lf ", &tmp[Input->NThom]);
      if(nread != 1) break;
      Input->NThom++;
      token = strtok(NULL, " ,");
    }
    if(Input->NThom > 0){
      double COnST = 3./8.*C_sigmaT*1e-10*C_Solarradus;
      Input->Thom = (STRUCT_THOMSON *) \
          malloc(Input->NThom*sizeof(STRUCT_THOMSON));
      for(itmp = 0; itmp < Input->NThom; itmp++){
        Input->Thom[itmp].lambda = tmp[itmp]*1e-10;
        Input->Thom[itmp].Intensity = Planck(C_c/Input->Thom[itmp].lambda, \
            5800.0);
        Limb_Darkening(Input->Thom[itmp].lambda, &(Input->Thom[itmp].u1), \
            &(Input->Thom[itmp].u2));
        Input->Thom[itmp].Const = C_c*COnST*Input->Thom[itmp].Intensity \
            /Input->Thom[itmp].lambda/Input->Thom[itmp].lambda; 
      }  
    }
      
    indx = 9;
    Input->rmin = atof(Keywords[indx].line);
    Input->rsqmin = Input->rmin*Input->rmin;
    if(Input->Verbose >= 1 && Mpi->rank == 0){
      fprintf(stderr, "\n Min of the radius : %e Rs ", Input->rmin);
    }
    
    indx = 10;
    Input->rmax = atof(Keywords[indx].line);
    Input->rsqmax = Input->rmax*Input->rmax;
    if(Input->Verbose >= 1 && Mpi->rank == 0){
      fprintf(stderr, "\n Max of the radius : %e Rs ", Input->rmax);
    }
        
    indx = 11;
    nread = sscanf(Keywords[indx].line,"%lf, %lf, %lf, %lf ", \
        &Input->FOV[0][0], \
        &Input->FOV[0][1], &Input->FOV[1][0], &Input->FOV[1][1]);

    if(Input->Verbose >= 1 && Mpi->rank == 0){
      fprintf(stderr, "\n field of view ");
      fprintf(stderr, "\n X : %e Rs %e Rs ", Input->FOV[0][0], \
          Input->FOV[0][1]);
      fprintf(stderr, "\n Z : %e Rs %e Rs ", Input->FOV[1][0], \
          Input->FOV[1][1]);
    }
    
    indx = 12;
    Input->ny = atoi(Keywords[indx].line);
    Input->dy = (Input->FOV[0][1]-Input->FOV[0][0])/(Input->ny-1);
    if(Input->Verbose >= 1 && Mpi->rank == 0){
      fprintf(stderr, "\n Y size: %d ", Input->ny);
    }

    indx = 13;
    Input->nz = atoi(Keywords[indx].line);
    Input->dz = (Input->FOV[1][1]-Input->FOV[1][0])/(Input->nz-1);
    if(Input->Verbose >= 1 && Mpi->rank == 0){
      fprintf(stderr, "\n Z size: %d ", Input->nz);
    }
    
    indx = 14;
    Input->dx = atof(Keywords[indx].line);
    if(Input->Verbose >= 1 && Mpi->rank == 0){
      fprintf(stderr, "\n delta X : %f ", Input->dx);
    }

    indx = 15;
    Input->Kmax = atoi(Keywords[indx].line);
    if(Input->Verbose >= 1 && Mpi->rank == 0){
      fprintf(stderr, "\n Kmax : %d  ", Input->Kmax);
    }

    indx = 16;
    Input->Kdelta = atoi(Keywords[indx].line);
    if(Input->Verbose >= 1 && Mpi->rank == 0){
      fprintf(stderr, "\n Kdelta : %d  ", Input->Kdelta);
    }

    /////////////////////////////
    indx = 17;
    if(Input->Mode == 0){
      String_Copy(Input->Path_Observation, Keywords[indx].line, \
          strlen(Keywords[indx].line), true);
      if(!FILE_EXIST(Input->Path_Observation)){
        Error(enum_error, "routine_name", "zeros file doesn't exist!");
      }
      if(Input->Verbose >= 1 && Mpi->rank == 0){
        fprintf(stderr, "\n Observation file : %s ", \
            Input->Path_Observation);
      }
    }
    
    indx = 18;
    String_Copy(Input->Path_Output, Keywords[indx].line, \
        strlen(Keywords[indx].line), true);
    if(Input->Verbose >= 1 && Mpi->rank == 0){
      fprintf(stderr, "\n Output file : %s ",Input->Path_Output);
    }
    
    indx = 19;
    String_Copy(Input->Path_coeff, Keywords[indx].line, \
        strlen(Keywords[indx].line), true);
    if(Input->Verbose >= 1 && Mpi->rank == 0){
      fprintf(stderr, "\n SFB Coefficients : %s ", Input->Path_coeff);
    }
    
    indx = 20;
    Input->step_size = atof(Keywords[indx].line);
    if(Input->Verbose >= 1 && Mpi->rank == 0){
      fprintf(stderr, "\n step size: %e ", Input->step_size);
    }
    
    
    indx = 21;
    Input->nweight = String_elements(Keywords[indx].line);
    Input->weight = (double *)malloc(sizeof(double)*Input->nweight);
    for(itmp=0;itmp<Input->nweight;itmp++){
      String_Split(parameter, Keywords[indx].line);
      Input->weight[itmp] = atof(parameter);
    }
    if(Input->Verbose >= 1 && Mpi->rank == 0){
      fprintf(stderr, "\n Weights array : ");
      for(itmp=0;itmp<Input->nweight;itmp++){
        fprintf(stderr, " %f ", Input->weight[itmp]); 
      }
    }
    
    indx = 22;
    String_Copy(Input->Path_Atmo, Keywords[indx].line, \
        strlen(Keywords[indx].line), true);
    if(Input->Verbose >= 1 && Mpi->rank == 0){
      fprintf(stderr, "\n atmosphere model : %s ", Input->Path_Atmo);
    }

    indx = 23;
    Input->Per = atof(Keywords[indx].line);
    if(Input->Verbose >= 1 && Mpi->rank == 0){
      fprintf(stderr, "\n dri - per : %f  ", Input->Per);
    }

    indx = 24;
    nread = sscanf(Keywords[indx].line,"%lf, %lf, %lf, %lf", \
        &(Input->perturb[0]), &(Input->perturb[1]), \
        &(Input->perturb[2]), &(Input->perturb[3]));
    if(Input->Verbose >= 1 && Mpi->rank == 0){
      fprintf(stderr, "\n perturbation : %e %e %e %e ", \
          (Input->perturb[0]), (Input->perturb[1]), \
          (Input->perturb[2]), (Input->perturb[3]));
    }

    indx = 25;
    String_to_Upper(Keywords[indx].line);
    if(strcmp(Keywords[indx].line, "YES") == 0){
      Input->symmetry = true;
    }else{
      Input->symmetry = false;
    }
    if(Input->Verbose >= 1 && Mpi->rank == 0){
      if(Input->symmetry){
        fprintf(stderr, "\n axial symmetry : YES ");
      }else{
        fprintf(stderr, "\n axial symmetry : No ");
      }
    }

    if(Input->Mode == 2 || Input->Mode == 3){
      indx = 26;
      nread = sscanf(Keywords[indx].line,"%lf, %lf, %lf", \
          &(Input->Ypos), &(Input->Zpos), \
          &(Input->Xpos));
      if(Input->Verbose >= 1 && Mpi->rank == 0){
        if(nread == 3){
          fprintf(stderr, "\n Y = %f Z = %f X = %f ", \
              Input->Ypos, Input->Zpos, Input->Xpos);
        }else{
          fprintf(stderr, "\n Y = %f Z = %f ", \
              Input->Ypos, Input->Zpos);
        }
      }
    }

    indx = 27;
    String_to_Upper(Keywords[indx].line);
    nread = sscanf(Keywords[indx].line,"%s %s %s %s", \
          inv[0], inv[1], inv[2], inv[3]);

    for(itmp=0;itmp<4;itmp++){
      Trim(inv[itmp],3);
      if(strcmp(inv[itmp],"NO")==0||strcmp(inv[itmp],"N")==0){
        Input->Para[itmp].invt = false;
      }else{
        Input->Para[itmp].invt = true;
      }
      if(Input->Verbose >= 1 && Mpi->rank == 0){
        if(Input->Para[itmp].invt){
          fprintf(stderr,"\n inv parameter %d: Yes ",itmp);
        }else{
          fprintf(stderr,"\n inv parameter %d: No ",itmp);
        }
      }
    }

    indx = 28;
    Input->inner_flag = atoi(Keywords[indx].line);
    if(Input->Verbose >= 1 && Mpi->rank == 0){
      fprintf(stderr, "\n Inner : %d  ", Input->inner_flag);
    }

    indx = 29;
    Input->Ncut = atoi(Keywords[indx].line);
    if(Input->Verbose >= 1 && Mpi->rank == 0){
      fprintf(stderr, "\n N cut : %d  ", Input->Ncut);
    }

    indx = 30;
    Input->Lcut = atoi(Keywords[indx].line);
    if(Input->Verbose >= 1 && Mpi->rank == 0){
      fprintf(stderr, "\n K cut : %d  ", Input->Lcut);
    }

    indx = 31;
    String_to_Upper(Keywords[indx].line);
    if(strcmp(Keywords[indx].line, "YES") == 0){
      Input->Bpotential = true;
    }else{
      Input->Bpotential = false;
    }
    if(Input->Verbose >= 1 && Mpi->rank == 0){
      if(Input->Bpotential){
        fprintf(stderr, "\n Vector Potential : YES ");
      }else{
        fprintf(stderr, "\n Vector Potential : No ");
      }
    }

    indx = 32;
    Input->Niter = atoi(Keywords[indx].line);
    if(Input->Verbose >= 1 && Mpi->rank == 0){
      fprintf(stderr, "\n Niter : %d  ", Input->Niter);
    }

    indx = 33;
    if(Keywords[33].Set){
      nread = sscanf(Keywords[indx].line,"%lf, %lf, %lf, %lf ", \
          &Input->FOVSPEC[0][0], &Input->FOVSPEC[0][1], \
          &Input->FOVSPEC[1][0], &Input->FOVSPEC[1][1]);
      if(Input->FOVSPEC[0][0]<Input->FOV[0][0]) \
          Input->FOVSPEC[0][0]=Input->FOV[0][0];
      if(Input->FOVSPEC[0][1]>Input->FOV[0][1]) \
          Input->FOVSPEC[0][1]=Input->FOV[0][1];
      if(Input->FOVSPEC[1][0]<Input->FOV[1][0]) \
          Input->FOVSPEC[1][0]=Input->FOV[1][0];
      if(Input->FOVSPEC[1][1]>Input->FOV[1][1]) \
          Input->FOVSPEC[1][1]=Input->FOV[1][1];
    }else{
      Input->FOVSPEC[0][0]=Input->FOV[0][0];
      Input->FOVSPEC[0][1]=Input->FOV[0][1];
      Input->FOVSPEC[1][0]=Input->FOV[1][0];
      Input->FOVSPEC[1][1]=Input->FOV[1][1];
    }
    if(Input->Verbose >= 1 && Mpi->rank == 0){
      fprintf(stderr, "\n Spectrum region ");
      fprintf(stderr, "\n X : %e Rs %e Rs ", Input->FOVSPEC[0][0], \
          Input->FOVSPEC[0][1]);
      fprintf(stderr, "\n Z : %e Rs %e Rs ", Input->FOVSPEC[1][0], \
          Input->FOVSPEC[1][1]);
    }

    indx = 34;
    spe_tmp = Input->Spec;
    if(Input->Nspec>0 && Input->FOVSPEC[0][1]>=Input->FOVSPEC[0][0] \
        && Input->FOVSPEC[1][1]>=Input->FOVSPEC[1][0]){
      Input->Spec = (STRUCT_SPEC *)malloc(Input->Nspec*sizeof(STRUCT_SPEC));

      for(itmp=0; itmp<Input->Nspec; itmp++){
        Input->Spec[itmp].Nl = spe_tmp[itmp].Nl;
        if(spe_tmp[itmp].range[1]>spe_tmp[itmp].range[0]){
          Input->Spec[itmp].range[0] = \
              spe_tmp[itmp].range[0];
          Input->Spec[itmp].range[1] = \
              spe_tmp[itmp].range[1];
        }else{
          Input->Spec[itmp].range[0] = \
              spe_tmp[itmp].range[1];
          Input->Spec[itmp].range[1] = \
              spe_tmp[itmp].range[0];
        }
        if(Input->Verbose >= 1 && Mpi->rank == 0){
          fprintf(stderr, "\n spectrum : %d range: %e %e, grids %d", \
              itmp+1, Input->Spec[itmp].range[0], \
              Input->Spec[itmp].range[1], \
              Input->Spec[itmp].Nl);
        }
      }
    }else{
      Input->Nspec = 0;
      Input->Spec = NULL;
    }
    free(spe_tmp);

    indx = 35;
    String_to_Upper(Keywords[indx].line);
    if(strcmp(Keywords[indx].line, "YES") == 0){
      Input->OutputV = true;
    }else{
      Input->OutputV = false;
    }
    if(Input->Verbose >= 1 && Mpi->rank == 0){
      if(Input->OutputV){
        fprintf(stderr, "\n Output Stokes V image : Yes ");
      }else{
        fprintf(stderr, "\n Output Stokes V image : No ");
      }
    }

    indx = 36;
    Input->rint = atof(Keywords[indx].line);
    Input->rsqint = Input->rint*Input->rint;
    if(Input->Verbose >= 1 && Mpi->rank == 0){
      if(Input->rint<Input->rmax){
        Input->rint = Input->rmax;
        fprintf(stderr, "\n rIntegration is set to rmax ");
      }
     
      fprintf(stderr, "\n Max radius for integration : %e Rs ", \
          Input->rint); 
    }

    indx = 37;
    Input->nlos = String_elements(Keywords[indx].line);
    Input->los = (double *)malloc(sizeof(double)*Input->nlos);
    for(itmp=0;itmp<Input->nlos;itmp++){
      String_Split(parameter, Keywords[indx].line);
      Input->los[itmp] = atof(parameter);
      Input->los[itmp] *= (C_Pi/180.);
    }
    if(Input->Verbose >= 1 && Mpi->rank == 0){
      fprintf(stderr, "\n viewing directions : ");
      for(itmp=0;itmp<Input->nlos;itmp++){
        fprintf(stderr, " %f ", Input->los[itmp]/C_Pi*180.); 
        
      }
    }

    indx = 38;
    String_to_Upper(Keywords[indx].line);
    if(strcmp(Keywords[indx].line, "YES") == 0){
      Input->OutputPara = true;
    }else{
      Input->OutputPara = false;
    }
    if(Input->Verbose >= 1 && Mpi->rank == 0){
      if(Input->OutputPara){
        fprintf(stderr, "\n Output model parameters : YES ");
      }else{
        fprintf(stderr, "\n Output model parameters : No ");
      }
    }

    indx = 39;
    String_Copy(Input->Path_Para, Keywords[indx].line, \
        strlen(Keywords[indx].line), true);
    if(Input->Verbose >= 1 && Mpi->rank == 0){
      fprintf(stderr, "\n path to model parameter : %s ", 
          Input->Path_Para);
    }

    if(Input->Mode == 2){
      indx = 40;
      nread = sscanf(Keywords[indx].line,"%lf, %lf, %lf", \
          &(Input->Bvec[0]), &(Input->Bvec[1]), \
          &(Input->Bvec[2]));
      if(Input->Verbose >= 1 && Mpi->rank == 0){
        if(nread == 3){
          Input->Binput = true;
          Input->Bvec[1] *= 180./C_Pi; 
          Input->Bvec[2] *= 180./C_Pi; 
          fprintf(stderr, "The input magnetic field: %f, %f, %f ", \
              Input->Bvec[0], Input->Bvec[1], \
              Input->Bvec[2]);
        }else{
          Input->Binput = false;
        }
      }
    }

    if(Mpi->rank == 0){
      fprintf(stderr, "\n\n ------------------" \
          "---------------------- \n\n");
    }

    return 0;
}

/*----------------------------------------------------------------------------*/

extern int RDINPUT(char Filename[], STRUCT_INPUT *Input, STRUCT_MPI *Mpi){
  
    /*######################################################################
      Purpose:
        Read the input file for the forbidden line calculation.
      Record of revisions:
        27 Dec. 2024.
      Input parameters:
        Filename[], the input file.
      Output parameters:
        Input, a structure saved the input information.
        Input, a structure saved the inversion information.
        Output, a structure saved the output information.
        Mpi, a structure saved the Mpi information.
    ######################################################################*/
    
    const char *routine_name = "RDINPUT";
    
    FILE *fa = fopen(Filename, "r");
        
    int Num_Keywords = 41;

    STRUCT_KEYS Keywords[] ={
      {"verbose", "1", false, false},  //0
      {"Mode", "Forward", false, false},  //1
      {"atom_path", "", false, true}, //2
      {"path_zero", "", false, true}, //3
      {"Bessel_function_N", "", false, true}, //4
      {"harmonic_L", "", false, true}, //5
      {"sparsity", "", false, true}, //6
      {"saturated", "YES", false, false}, //7
      {"Thomson_scattering", "", false, false}, //8
      {"rmin", "", false, true}, //9     
      {"rmax", "", false, true}, //10
      {"FOV", "", false, true}, //11
      {"Ysize", "", false, true}, //12
      {"Zsize", "", false, true}, //13
      {"deltax", "", false, true}, //14
      {"Kmax", "2", false, false}, //15
      {"Kdelta", "2", false, false}, //16
      {"observation", "", false, true}, //17
      {"output", "", false, true}, //18
      {"path_coeff", "", false, true}, //19
      {"Step_size", "", false, true}, //20
      {"Weights", "1.0, 2.0, 2.0, 5.0, 8.0", false, true}, //21
      {"atmosphere_model", "", false, false}, //22
      {"per", "0.5", false, false}, //23
      {"perturbation", "1e-3, 1e-3, 1e-3, 1e-3", false, false}, //24
      {"symmetry", "YES", false, false},  //25
      {"coordinate", "", false, false},  //26
      {"inv_flag", "YES YES YES YES", false, false},  //27
      {"innerflag", "1", false, false},  //28
      {"Ncut", "100", false, false},  //29
      {"Lcut", "100", false, false},  //30
      {"Bpotential", "No", false, false},  //31
      {"Niter", "10000", false, false},  //32
      {"FOVSPEC", "1.0, -1.0, 1.0, -1.0", false, false}, //33
      {"SPEC", "", false, false}, //34
      {"OutputV", "No", false, false}, //35
      {"rIntegration", "3.0", false, false}, //36
      {"LOS", "0.0", false, false}, //37
      {"OutputPara", "No", false, false}, //38
      {"Parapath", "./para.bin", false, true}, //39
      {"Bvec", "", false, false} //40
    };
    Num_Keywords = sizeof(Keywords)/sizeof(STRUCT_KEYS);
    
    char lines[Max_Line_Length], key[Key_Length], *p;
    int len, itmp;
    long len_tot;
    bool neglect;
    Input->Verbose = 1;

    Input->Path_Atom = (char **)MATRIX(0, 60, 0, Max_Line_Length-1, \
        enum_char, true);
    Input->Spec = (STRUCT_SPEC *)malloc(60*sizeof(STRUCT_SPEC));

    Input->Natom = 0;
    Input->Nspec = 0;

    while(Read_line(lines, fa) > 0){
      len_tot = strlen(lines);
      len = Indx_Char(lines, '=', 1);
      if(len > 1){
        len_tot -= len;
        p = lines+len;
        String_Copy(key, lines, len-1, true);
        Trim(key, 3);
        
        neglect = true;
        if(strcmp(key,Keywords[2].keyword)==0){
          String_Copy(Input->Path_Atom[Input->Natom], p, len_tot, false);
          Input->Natom++;
          neglect = false;
          Keywords[2].Set = true;
        }else if(strcmp(key,Keywords[34].keyword)==0){
          sscanf(p,"%lf %lf %d", &(Input->Spec[Input->Nspec].range[0]), \
              &(Input->Spec[Input->Nspec].range[1]), \
              &(Input->Spec[Input->Nspec].Nl));
          if(Input->Spec[Input->Nspec].Nl>1 \
              &&Input->Spec[Input->Nspec].range[0] \
              !=Input->Spec[Input->Nspec].range[1]) Input->Nspec++;
          Keywords[34].Set = true;
          neglect = false;
        }else{
          for (itmp=0; itmp<Num_Keywords; itmp++){
            if(strcmp(key,Keywords[itmp].keyword)==0){
              if(Keywords[itmp].Set == true && (Mpi->rank) == 0){
                Error(enum_warning, routine_name, "read a keyword twice");
                fprintf(stderr,"keyword = %s \n", key);
              }
              String_Copy(Keywords[itmp].line, p, len_tot, true);
              Keywords[itmp].Set = true;
              if(Input->Verbose >= 2  && (Mpi->rank) == 0){
                fprintf(stderr," %d %s %s \n",itmp, Keywords[itmp].keyword, \
                    Keywords[itmp].line);
              }
              neglect = false;
              break;
            }
          }
        }
        if(Input->Verbose >= 1 && Mpi->rank == 0 && neglect){
          fprintf(stderr, "Warning : Neglect key words %s \n",key);
          fprintf(stderr, "The whole line is %s \n",lines);
        }
      }
    }
    fclose(fa);

    Keywords_Conversion(Keywords, Input, Mpi);

    Input->rJmax = 0;
    Input->rLmax = 0;
    Input->rSmax = 0;

    Input->Nline = 0;
    Input->Ntrans = 0;
    Input->Maxntran = 0;

    if(Input->Kdelta==2) Input->Kmax = (Input->Kmax/2)*2;

    if(Input->Kmax==0) Input->Saturated = true;

    if(Input->Saturated){
      Input->nJKQ = 2;
      Input->Kdelta = 2;
      Input->TQmax = 0;
    }else{
      if(Input->Kdelta==2){
        Input->nJKQ = 6;
      }else{
        Input->nJKQ = 9;
      }
      Input->TQmax = 2;
    }

    if(Input->symmetry){
      Input->JQmax = 0;
    }else{
      Input->JQmax = 2;
    }

    if(Input->Mode==0) Input->Nspec = 0;
    if(Input->Nspec>0) Input->OutputV = true;

    if(Input->OutputV){
      Input->Nstk = 4;
    }else{
      Input->Nstk = 3;
    }

    if(Input->Mode==2||Input->Mode==3){
      Input->ny = 1;
      Input->nz = 1;
    }

    // if mode 0 (inversion), or mode 1 (forward)
    if(Input->Mode<=3){
      Input->nx = (int **)MATRIX(0, Input->ny-1, 0, \
          Input->nz-1, enum_int, true);
    }

    if(Mpi->rank==0){
      if(FILE_EXIST(Input->Path_Para)){
         remove(Input->Path_Para);
      }
    }

    return 0;
}

/*----------------------------------------------------------------------------*/

