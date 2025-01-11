
#include "INIT.h"

/*----------------------------------------------------------------------------*/

    /*######################################################################
    
      revision log:

        19 Nov. 2024
          --- update: output the nx for each pixel.

        30 Otc. 2024
          --- bugfix: fix a typo in multi-direction synthesis .
          --- update: moved subroutine FREE_GRIDS and FREE_OUTPUT to 
                      FREE.c
                      moved subroutine Randomseeds to MPI_CTRL.c

        19 Sep. 2024
          --- update: synthesis for multi LOS.
    
    ######################################################################*/

/*----------------------------------------------------------------------------*/

extern int Init(STRUCT_INPUT *Input, STRUCT_ATOM *Atom, STRUCT_SYN *Syn, \
    STR_FCTSG *fctsg, STRUCT_OUT *Output, STRUCT_MPI *Mpi){
  
    /*######################################################################
      Purpose:
        initialize the inversion or synthesis.
      Record of revisions:
        19 Nov 2024.
      Input parameters:
        Input, a structure with the inputs.
        Atom, a structure with the atoms.
        fctsg, a structure with factoral, signs.
      Output parameters:
        Input, a structure with the inversion information.
        Syn, a structure with the grids.
        Output, a structure with the output matrix
        Mpi, a structure with the Mpi information.
    ######################################################################*/

    int K, au, al, ix, iy, iz, nx = 0, itmp, ilos, iThom;
    int iatom, itrans, igrid, ipixel, ipara, indx, ipspec, ispec, ilambda;
    int Si, q1, q, in, il, im, K_ln, costheta;
    double dtmp1, dtmp2, X, Y, Z, rsq0, rsq1, PAB[2], Jkq[2], \
        Hcoeff, Ju, Jl;


    complex double **Dkmn[3];

    Dkmn[2] = (complex double **)MATRIX(-2, 2, -2, 2, enum_cplx, false);

    complex double ***T_KQ = (complex double ***)TENSOR_RHO_CPLX(2, 2, \
        false);

    STRUCT_GRID *pgrid;
    STRUCT_LOS *plos;

    fctsg->nmax = lround(8*Input->rJmax)> \
        4*lround(Input->rJmax+Input->rLmax+Input->rSmax)? \
        lround(8*Input->rJmax): \
        4*lround(Input->rJmax+Input->rLmax+Input->rSmax);

    Input->Tkq = (complex double ***)TENSOR_RHO_CPLX(2, 2, true);

    Input->Dkmn[2] = (complex double **)MATRIX(-2, 2, -2, 2, \
        enum_cplx, false);

    Input->Nmax=0;
    Input->Lmax=0;

    // compute the coefficient number if mode 0 (inversion)
    if(Input->Mode == 0){

      Randomseeds(Mpi);
      // read the zeros 
      READ_ZEROS(Input->Path_Zero, Input, Mpi);

      Input->ncoeff = 0;
      
      for (ipara=0; ipara<4; ipara++){
        if(!Input->Para[ipara].invt) continue;
        Input->Para[ipara].Lsquare = (Input->Para[ipara].L+1) \
            *(Input->Para[ipara].L+1);
        Input->Para[ipara].NLsquare = Input->Para[ipara].Lsquare \
            *Input->Para[ipara].N;

        Input->Para[ipara].threshold_n = (int)(Input->Para[ipara].NLsquare \
            *Input->Para[ipara].sparsity);

        Input->ncoeff += Input->Para[ipara].NLsquare;
        if(Input->Para[ipara].N>Input->Nmax) Input->Nmax = \
            Input->Para[ipara].N;
        if(Input->Para[ipara].L>Input->Lmax) Input->Lmax = \
            Input->Para[ipara].L;

        if(Input->Verbose>=1 && Mpi->rank == 0){
          fprintf(stderr,"\n parameter %d : \n", ipara);
          fprintf(stderr,"  Lsquare %d NLsquare %d: \n", \
              Input->Para[ipara].Lsquare, Input->Para[ipara].NLsquare);
        }
        Input->Para[ipara].Coeff = (complex double ***) \
            TENSOR_TRI_CPLX(Input->Para[ipara].N, Input->Para[ipara].L, \
            true);
      }

      if(Input->Verbose>=1 && Mpi->rank == 0){
        fprintf(stderr,"\n total coefficient number : %d \n", \
            Input->ncoeff);
      }

      if(fctsg->nmax<Input->Lmax*2) fctsg->nmax=Input->Lmax*2;
    }

    // precompute the sign the factorials
    INIT_FCTSG(fctsg);
    fctsg->memo = false;

    // get the wavelength if output the spectrum
    if(Input->Mode<4&&Input->Mode>0){
      if(Input->Nspec>0){
        for(ispec=0; ispec<Input->Nspec; ispec++){
          dtmp1 = (Input->Spec[ispec].range[1]-Input->Spec[ispec].range[0]) \
              /(Input->Spec[ispec].Nl-1);
          Input->Spec[ispec].Lambda = (double *)malloc( \
              Input->Spec[ispec].Nl*sizeof(double));  
          Input->Spec[ispec].Lambda[0] = Input->Spec[ispec].range[0];
          Input->Spec[ispec].Lambda[Input->Spec[ispec].Nl-1] \
              = Input->Spec[ispec].range[1];
          for(ilambda=1;ilambda<Input->Spec[ispec].Nl-1;ilambda++){
            Input->Spec[ispec].Lambda[ilambda] \
                = Input->Spec[ispec].Lambda[0]+dtmp1*ilambda;
          }
        }
      }
    }

    // compute the w and delta factors
    for(iatom = 0; iatom <Input->Natom; iatom++){
      for(indx=0; indx<Atom[iatom].Nline; indx++){
        itrans = Atom[iatom].iout[indx];

        au = Atom[iatom].TR[itrans].au;
        al = Atom[iatom].TR[itrans].al;
        Ju = Atom[iatom].LV[au].J;
        Jl = Atom[iatom].LV[al].J;

        Atom[iatom].TR[itrans].w[0] = 1;
        dtmp1 = fctsg->sg[(int)(1+Ju+Jl)]*sqrt(3.*Atom[iatom].LV[au].deg);
        for(K=1; K<3; K++){
          Atom[iatom].TR[itrans].w[K] = dtmp1*WIGNER_6J(1.0, 1.0, K, \
              Ju, Ju, Jl,fctsg);
        }

        dtmp1 = sqrt(Ju*(Ju+1.)*(Ju*2.+1.));
        dtmp2 = sqrt(Jl*(Jl+1.)*(Jl*2.+1.));

        Atom[iatom].TR[itrans].delta = -3.0*Atom[iatom].LV[au].sqrt_deg \
            *(fctsg->sg[lround(1+Jl-Ju)]*Atom[iatom].LV[au].g*dtmp1 \
            *WIGNER_6J(2.,1.,1.,Ju,Ju,Ju,fctsg) \
            *WIGNER_6J(1.,1.,1.,Jl,Ju,Ju,fctsg) \
            +Atom[iatom].LV[al].g*dtmp2 \
            *WIGNER_9J(Jl,Ju,1.,Jl,Ju,1.,1.,2.,1.,fctsg));

      }
    }

    // if mode 0 (inversion), or mode 1 (forward)
    if(Input->Mode<=1){
      
      // get the total numbers of pixels and grids
      Syn->ngrids = 0;
      Syn->npixels = 0;

      for(iy=0; iy<Input->ny; iy++){
        Y = Input->FOV[0][0]+iy*Input->dy;
        for(iz=0; iz<Input->nz; iz++){
          Z = Input->FOV[1][0]+iz*Input->dz;
          rsq0 = Y*Y+Z*Z;
          
          if(rsq0<Input->rsqmin || rsq0>Input->rsqmax){ 
            Input->nx[iy][iz] = -1;
            continue;
          }
          dtmp1 = sqrt(Input->rsqint-rsq0);
          if(dtmp1<Input->dx){ 
            Input->nx[iy][iz] = -1;
            continue;
          }
          nx = (int)(dtmp1/Input->dx);
          Input->nx[iy][iz] = nx;
          Syn->ngrids+=2*nx+1;
          Syn->npixels++;      
        }
      }

      // print the numbers
      if(Mpi->rank == 0)fprintf(stderr, " Total grids = %d Total "\
          "pixels = %d\n", Syn->ngrids,  Syn->npixels);

      // if mode 0 (inversion), allocate the nomalization factor 
      if(Input->Mode==0){
        Output->R = (double *)VECTOR(0, Syn->npixels-1, enum_dbl, false);
        Output->norm = (double *)VECTOR(0, Syn->npixels-1, \
            enum_dbl, false);
      }

      // get the numbers of pixels and grids for each processor 
      // and the corresponding grid index
      Mpi->ngrids = Syn->ngrids/Mpi->nprocs;
      Mpi->gindx[0] = Mpi->ngrids*Mpi->rank;

      // get and assign the resdue 
      itmp = Syn->ngrids%Mpi->nprocs;
      if(itmp>0){
        if(Mpi->rank <= itmp && Mpi->rank>0){
          Mpi->ngrids++;
          Mpi->gindx[0] = Mpi->gindx[0]+Mpi->rank;
        }else if(Mpi->rank > itmp){
          Mpi->gindx[0] = Mpi->gindx[0]+itmp;
        }
      }
      Mpi->gindx[1] = Mpi->gindx[0]+Mpi->ngrids-1;
      
      // allocate the RAM for the grids
      Syn->Grids = (STRUCT_GRID *) \
          malloc(sizeof(STRUCT_GRID)*Mpi->ngrids);
      
      // print the grid number and index for each processor 
      fprintf(stderr, "  rank = %d Total grids = %d. "\
          "start from %d to %d\n", Mpi->rank, Mpi->ngrids, \
          Mpi->gindx[0], Mpi->gindx[1]);

      // initialize the indexes of grid, pixel, and spectrum
      igrid = -1;
      ipixel = 0;
      ipspec = 0;

      
      for(iy=0; iy<Input->ny; iy++){
        Y = Input->FOV[0][0]+iy*Input->dy;
        if(Y>=Input->FOVSPEC[0][0]){
          Input->FOVSPECINDX[0][0] = iy;
          break;
        } 
      }
      Input->FOVSPECINDX[0][1] = Input->ny-1;
      for(iy=Input->FOVSPECINDX[0][0]; iy<Input->ny; iy++){
        Y = Input->FOV[0][0]+iy*Input->dy;
        if(Y>Input->FOVSPEC[0][1]){
          Input->FOVSPECINDX[0][1] = iy;
          break;
        } 
      }
      Input->nys = Input->FOVSPECINDX[0][1]-Input->FOVSPECINDX[0][0]+1;


      for(iz=0; iz<Input->nz; iz++){
        Z = Input->FOV[1][0]+iz*Input->dz;
        if(Z>=Input->FOVSPEC[1][0]){
          Input->FOVSPECINDX[1][0] = iz;
          break;
        } 
      }
      Input->FOVSPECINDX[1][1] = Input->nz-1;
      for(iz=Input->FOVSPECINDX[1][0]; iz<Input->nz; iz++){
        Z = Input->FOV[1][0]+iz*Input->dz;
        if(Z>Input->FOVSPEC[1][1]){
          Input->FOVSPECINDX[1][1] = iz;
          break;
        } 
      }
      Input->nzs = Input->FOVSPECINDX[1][1]-Input->FOVSPECINDX[1][0]+1;

      // get the pixel, and spectrum indexes for each grid
      for(iy=0; iy<Input->ny; iy++){
        Y = Input->FOV[0][0]+iy*Input->dy;
        for(iz=0; iz<Input->nz; iz++){
          Z = Input->FOV[1][0]+iz*Input->dz;
          rsq0 = Y*Y+Z*Z;

          // if the grid between rmin and rmax
          nx = Input->nx[iy][iz];
          if(nx<0) continue;

          for(ix=-nx;ix<=nx;ix++){
            igrid++;
            // check the grid index for the processor
            if(igrid<Mpi->gindx[0]||igrid>Mpi->gindx[1]) continue;
                
            X = ix*Input->dx;
            rsq1 = rsq0+X*X;
            itmp = igrid-Mpi->gindx[0];

            // get the coordinate of the grid
            Syn->Grids[itmp].Rsq0 = rsq0;
            Syn->Grids[itmp].R = sqrt(rsq1);
            Syn->Grids[itmp].Theta = acos(Z/Syn->Grids[itmp].R);
            Syn->Grids[itmp].Phi = atan2(Y, X);
            if(Syn->Grids[itmp].Phi < 0){
              Syn->Grids[itmp].Phi += C_Pi*2;
            }

            // get the pixel index for the grid
            Syn->Grids[itmp].ipixel = ipixel;

            // get the spectrum index if output the spectrum 
            if(Input->Nspec>0){
              if(iy>=Input->FOVSPECINDX[0][0] \
                  && iy<=Input->FOVSPECINDX[0][1] \
                  && iz>=Input->FOVSPECINDX[1][0] \
                  && iz<=Input->FOVSPECINDX[1][1]){
                Syn->Grids[itmp].ipspec = ipspec;

              }else{
                Syn->Grids[itmp].ipspec = -1;
              }
            }
          }

          // compute the normalization factor if mode 0 (inversion)
          if(Input->Mode==0){
            Output->R[ipixel] = sqrt(rsq0);
            Output->norm[ipixel] = 1.0/Baumbach(Output->R[ipixel]);
            //Output->norm[ipixel] = 1.0;
          }

          if(Input->Nspec>0){
            if(iy>=Input->FOVSPECINDX[0][0] \
                && iy<=Input->FOVSPECINDX[0][1] \
                && iz>=Input->FOVSPECINDX[1][0] \
                && iz<=Input->FOVSPECINDX[1][1]){
              ipspec++;  
            }
          }
          ipixel++;
        }
      }

      Syn->npspec = ipspec;

    // if mode 2 (single grid)
    }else if(Input->Mode==2){

      // set the grid, pixel, and spectrum numbers
      Syn->ngrids = 1;
      Syn->npixels = 1;
      Input->nx[0][0] = 0;

      if(Input->Nspec>0){
        Syn->npspec = 1;
      }else{
        Syn->npspec = 0;
      }


      // only the first processor does the calculation since 1 grid
      if(Mpi->rank == 0){

        Mpi->ngrids = 1;
        Mpi->gindx[0] = 0;
        Mpi->gindx[1] = 0;

        // allocate the RAM for the grids
        Syn->Grids = (STRUCT_GRID *)malloc(sizeof(STRUCT_GRID));

        rsq1 = Input->Xpos*Input->Xpos+Input->Ypos*Input->Ypos+ \
            Input->Zpos*Input->Zpos;

        // get the coordinate
        Syn->Grids[0].Rsq0 = Input->Ypos*Input->Ypos+ \
            Input->Zpos*Input->Zpos;
        Syn->Grids[0].R = sqrt(rsq1);
        Syn->Grids[0].Theta = acos(Input->Zpos/Syn->Grids[0].R);
        Syn->Grids[0].Phi = atan2(Input->Ypos, Input->Xpos);

        if(Syn->Grids[0].Phi < 0){
          Syn->Grids[0].Phi += C_Pi*2;
        }

        // set the pixel and spectrum index
        Syn->Grids[0].ipixel = 0;
        Syn->Grids[0].ipspec = 0;

      // set the other processors
      }else{

        Mpi->ngrids = 0;
        Mpi->gindx[0] = 0;
        Mpi->gindx[1] = -1;
      }

    // if mode 3 (single pixel)
    }else if(Input->Mode==3){



      rsq0 = Input->Ypos*Input->Ypos+Input->Zpos*Input->Zpos;

      dtmp1 = sqrt(Input->rsqint-rsq0);

      // get the grid number
      if(dtmp1>=Input->dx){
        nx = (int)(dtmp1/Input->dx);
        Input->nx[0][0] = nx;

        Syn->ngrids = 2*nx+1;
        Syn->npixels = 1;

        if(Input->Nspec>0){
          Syn->npspec = 1;
        }else{
          Syn->npspec = 0;
        }

      }else{
        const char *routine_name = "Init";
        Error(enum_error, routine_name, "Mode 3 error!");
      }

      // get the numbers of pixels and grids for each processor 
      // and the corresponding grid index
      Mpi->ngrids = Syn->ngrids/Mpi->nprocs;
      Mpi->gindx[0] = Mpi->ngrids*Mpi->rank;

      // get and assign the resdue 
      itmp = Syn->ngrids%Mpi->nprocs;

      if(itmp>0){
        if(Mpi->rank <= itmp && Mpi->rank>0){
          Mpi->ngrids++;
          Mpi->gindx[0] = Mpi->gindx[0]+Mpi->rank;
        }else if(Mpi->rank > itmp){
          Mpi->gindx[0] = Mpi->gindx[0]+itmp;
        }
      }
      Mpi->gindx[1] = Mpi->gindx[0]+Mpi->ngrids-1;
      
      // allocate the RAM for the grids
      Syn->Grids = (STRUCT_GRID *)malloc(sizeof(STRUCT_GRID) \
          *Mpi->ngrids);
      
      // print the numbers
      fprintf(stderr, "  rank = %d Total grids = %d. "\
          "start from %d to %d\n", Mpi->rank, Mpi->ngrids, \
          Mpi->gindx[0], Mpi->gindx[1]);
        
      // initialize the indexes of grid, pixel, and spectrum
      igrid = -1;

      // get the pixel, and spectrum indexes for each grid
      for(ix=-nx;ix<=nx;ix++){
        igrid++;

        // check the grid index for the processor
        if(igrid<Mpi->gindx[0]||igrid>Mpi->gindx[1]) continue;

        X = ix*Input->dx;
        rsq1 = rsq0+X*X;
        itmp = igrid-Mpi->gindx[0];

        // compute coordinate
        Syn->Grids[itmp].Rsq0 = rsq0;
        Syn->Grids[itmp].R = sqrt(rsq1);

        Syn->Grids[itmp].Theta = acos(Input->Zpos \
            /Syn->Grids[itmp].R);
        Syn->Grids[itmp].Phi = atan2(Input->Ypos, X);

        if(Syn->Grids[itmp].Phi < 0){
          Syn->Grids[itmp].Phi += C_Pi*2;
        }

        Syn->Grids[itmp].ipixel = 0;
        
        if(Input->Nspec>0){
          Syn->Grids[itmp].ipspec = 0;     
        }
      }        
    }

    // get the Tkq for 90 degree
    TKQ90(T_KQ);

    double *JJ = (double *)malloc((Input->Lmax+1)*sizeof(double));
    double *JJp = (double *)malloc((Input->Lmax+1)*sizeof(double));

    // precomputing for each grid
    for(igrid=0; igrid<Mpi->ngrids; igrid++){
      pgrid = Syn->Grids+igrid;
      
      pgrid->los = (STRUCT_LOS *)malloc(sizeof(STRUCT_LOS)*Input->nlos);

      for(ilos=0;ilos<Input->nlos;ilos++){
        plos = pgrid->los+ilos;

        // if output the Thomson scattering
        if(Input->NThom > 0){
          plos->Thom = (STRUCT_THOMSON_GRID *) \
              malloc(sizeof(STRUCT_THOMSON_GRID)*Input->NThom);
        }

        // if output spectral lines
        if(Input->Nline > 0){
        
          // allocate lines and Jkq
          plos->Line = (STRUCT_STOKES_GRID *) \
              malloc(sizeof(STRUCT_STOKES_GRID)*Input->Nline);
                  // allocate the ion fraction
          plos->Para.Ion = (double *)malloc(sizeof(double)*Input->Natom);

        }
      }

      // if output the Thomson scattering
      if(Input->NThom > 0){
        pgrid->Kr = (double *)malloc(sizeof(double)*Input->NThom);
        pgrid->Kt = (double *)malloc(sizeof(double)*Input->NThom);
        for(iThom=0; iThom< Input->NThom; iThom++){
          Thom_Scat(pgrid->R, PAB, Input->Thom[iThom].u1, \
              Input->Thom[iThom].u2);
          rsq0 = pgrid->Rsq0/pgrid->R/pgrid->R;
          pgrid->Kr[iThom] = (1-rsq0)*PAB[0]+rsq0*PAB[1];
          pgrid->Kt[iThom] = PAB[0];
        }
      }

      // if output spectral lines
      if(Input->Nline > 0){

        pgrid->Jkq = (STRUCT_INCIDENT_GRID *)\
              malloc(sizeof(STRUCT_INCIDENT_GRID)*Input->Natom);

        for(iatom=0; iatom<Input->Natom; iatom++){

          // allocate J00 and J2q
          pgrid->Jkq[iatom].J00 = (double *) \
              malloc(sizeof(double)*Atom[iatom].Ntrans);
          pgrid->Jkq[iatom].J2q = (double **)MATRIX(0, \
              Atom[iatom].Ntrans-1, -Input->JQmax, Input->JQmax, \
              enum_dbl, true);

          // comput the J00 and J20
          for(itrans=0; itrans<Atom[iatom].Ntrans; itrans++){
            Jkq_off(Atom[iatom].TR[itrans].u1, Atom[iatom].TR[itrans].u2, \
                pgrid->R, Jkq);
            pgrid->Jkq[iatom].J00[itrans] = Jkq[0] \
                *Atom[iatom].TR[itrans].Intens;
            pgrid->Jkq[iatom].J2q[itrans][0] = Jkq[1] \
                *Atom[iatom].TR[itrans].Intens;

            /*
              to do
              symmetry breaking!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!









            */
          }
        }
      }
      
      // allocate T2q
      pgrid->T2Q = (complex double **)MATRIX(0, 2, -2, 2, enum_cplx, \
          true);

      // get the rotation matrix 
      Rotmat(pgrid->Phi, pgrid->Theta, 0, Dkmn, 2);

      // T^K_P rotated to local vertical coordinate
      for(Si=0;Si<3;Si++){
        for(q=-2;q<=2;q++){
          for(q1=-2;q1<=2;q1++){
            pgrid->T2Q[Si][q] += T_KQ[Si][2][q1]*Dkmn[2][q1][q];
          }
        }
      }

      //if mode 0 (inversion)
      if(Input->Mode == 0){

        costheta = cos(pgrid->Theta);

        // allocate the ion fraction for perterbation
        for(ilos=0;ilos<Input->nlos;ilos++){

          plos = pgrid->los+ilos;

          plos->Para_unperturb = (STRUCT_PARA *) \
              malloc(sizeof(STRUCT_PARA));
          plos->Para_unperturb->Ion = (double *) \
              malloc(sizeof(double)*Input->Natom);

          // if output spectral lines
          if(Input->Nline > 0){
          
            plos->Line_unperturb = (STRUCT_STOKES_GRID *) \
                malloc(sizeof(STRUCT_STOKES_GRID)*Input->Nline);
          }

          // if output the Thomson scattering
          if(Input->NThom > 0){
            plos->Thom = (STRUCT_THOMSON_GRID *) \
                malloc(sizeof(STRUCT_THOMSON_GRID)*Input->NThom);
          }          

        }

        // allocate Bessel, and Legendre functions
        pgrid->sBessel = (double **)MATRIX(1, Input->Nmax, 0, \
            Input->Lmax, enum_dbl, false);
        pgrid->Legendre = (double **)MATRIX_TRI(Input->Lmax, enum_dbl, \
            false);

        // allocate the array for phi
        pgrid->Phiarray = (complex double *)VECTOR(0, Input->Lmax, \
            enum_cplx, false);

        // allocate derivatives of Bessel, and Legendre functions
        // if vector potential
        if(Input->Bpotential){
          pgrid->sBesselD = (double **)MATRIX(1, Input->Nmax, 0, \
              Input->Lmax, enum_dbl, false);
          pgrid->LegendreD = (double **)MATRIX_TRI(Input->Lmax, enum_dbl, \
              false);
        }

        // precomputing the functions
        for(il=0; il<=Input->Lmax; il++){   

          // compute the phi array
          pgrid->Phiarray[il] = cos(il*pgrid->Phi) \
              +sin(il*pgrid->Phi)*I;

          for(in=1; in<=Input->Nmax; in++){
            K_ln = Input->Zeros[il][in]/Input->rint;
            // compute the Bessel function
            Spherical_Bessel(K_ln*pgrid->R, il, JJ, JJp);
            pgrid->sBessel[in][il] = JJ[il];

            // compute derivative of the Bessel function 
            // if vector potential
            if(Input->Bpotential){
              pgrid->sBesselD[in][il] = JJp[il]*K_ln;
            }
          }
           
          for(im=0; im<=il; im++){
            Hcoeff = Harmonic_Coefficient(il, im);

            // compute the associated Legendre function
            pgrid->Legendre[il][im] = Hcoeff \
                *Associated_Legendre(il, im, costheta);

            // compute the derivative of the associated Legendre function
            // if vector potential (forgot the referrence)
            if(Input->Bpotential) \
                pgrid->LegendreD[il][im] = Hcoeff \
                  *0.5*((il+im)*(il-im+1)* \
                  Associated_Legendre(il, im-1, costheta) \
                  +Associated_Legendre(il, im+1, costheta));
          }
        }  
      }
    }

    // free RAM
    FREE_MATRIX(Dkmn[2], -2, -2, enum_cplx);
    FREE_TENSOR_RHO_CPLX(T_KQ);

    Output->los = (STRUCT_OUTLOS *)malloc(sizeof(STRUCT_OUTLOS) \
        *Input->nlos);
    for(ilos=0;ilos<Input->nlos;ilos++){
      Output->los[ilos].syn = (double **)MATRIX(0, Syn->npixels-1, 0, \
          Input->Nstk*Input->Nline+2*Input->NThom-1, enum_dbl, false);
    }
    
    // allocate the images
    Output->synloc = (double **)MATRIX(0, Syn->npixels-1, 0, \
        Input->Nstk*Input->Nline+2*Input->NThom-1, enum_dbl, false);

    // allocate the spectra
    if(Input->Nspec>0 && Syn->npspec>0){
      indx = 0;
      for(ispec=0; ispec<Input->Nspec; ispec++){
        indx += Input->Spec[ispec].Nl;
      }
      Input->Nl = indx;
      Output->specloc = (double **)MATRIX(0, Syn->npspec-1, 0, \
          Input->Nl*Input->Nstk-1, enum_dbl, false);

      for(ilos=0;ilos<Input->nlos;ilos++){
        Output->los[ilos].spec = (double **)MATRIX(0, Syn->npspec-1, 0, \
            Input->Nl*Input->Nstk-1, enum_dbl, false);;
      }  

    }else{
      Input->Nspec = 0;
      for(ispec=0; ispec<Input->Nspec; ispec++){
        free(Input->Spec[ispec].Lambda);
      }
      free(Input->Spec);
    }

    free(JJ);
    free(JJp);

    return 0;

}

/*----------------------------------------------------------------------------*/

extern int INIT_COEFF(STRUCT_SFB_COEFF *Parain, STRUCT_INPUT *Input){
    
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

    int ipara, in, il, im, Nmin, Lmin;

    for(ipara=0;ipara<4;ipara++){
      if(!Input->Para[ipara].invt) continue;
      Input->Para[ipara].in = Parain[ipara].in;
      if(!Parain[ipara].in) continue;
      Nmin = Parain[ipara].N<Input->Para[ipara].N? \
          Parain[ipara].N:Input->Para[ipara].N;
      Lmin = Parain[ipara].L<Input->Para[ipara].L? \
          Parain[ipara].L:Input->Para[ipara].L;

      if(Nmin>Input->Ncut) Nmin = Input->Ncut;
      if(Lmin>Input->Lcut) Lmin = Input->Lcut;
    
      for(in=1; in<=Nmin; in++){
        for(il=0; il<=Lmin; il++){
          for(im=0; im<=il; im++){
            Input->Para[ipara].Coeff[in][il][im] = \
                Parain[ipara].Coeff[in][il][im];
          }
        }
      } 
      FREE_TENSOR_RHO_CPLX(Parain[ipara].Coeff);     
    }

    return 0;
}

/*----------------------------------------------------------------------------*/

extern double Baumbach(double radius){
  
    /*######################################################################
      Purpose:
        Baumbach expressions from Astrophysical Quantities 3rd.
      Record of revisions:
        20 Sept. 2021.
      Input parameters:
        radius, projetcted radial distance.
      Return:
        The coronal intensity.
      Reference:
        Allen, 1973 Astrophysical Quantities 3rd.
    ######################################################################*/
    
    return 0.0532*pow(radius, -0.25)+1.425*pow(radius, -7)\
        +2.565*pow(radius, -17);

}

/*----------------------------------------------------------------------------*/

