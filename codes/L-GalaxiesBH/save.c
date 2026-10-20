/*  Copyright (C) <2016>  <L-Galaxies>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/> */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "allvars.h"
#include "proto.h"


/**@file save.c
 * @brief Copies the relevant properties in Galaxy structure into
 *        Galaxy_Output structure and saves them into the output
 *        files (SA_z**_**) - redshift/filenr.
 *
 *        There are two distinct procedures to write the output depending
 *        on whether GALAXY_TREE option is turned on or off. If it is on
 *        the full galaxy tree is written using save_galaxy_tree. If it
 *        is off, the output is only written for the chosen output snap
 *        numbers using save_galaxies.
 *
 *        If USE_MEMORY_TO_MINIMIZE_IO ON, these routines copy the galaxy
 *        data from the working structures into pointers until that has
 *        been done for all the tree in a given file.
 *
 *        After all the galaxy trees are written finalize_galaxy_file is
 *        called in main.c to include an header in the output files. If
 *        GALAXY_TREE=1 three numbers are written: 1 (int);
 *        size_of_struct(Galaxy_Output) (int);UnitTime_in_years and TotGalCount(int). If
 *        GALAXY_TREE=0 then the header is: Ntrees (int); total number of
 *        galaxies for the snapshot corresponding to the current file ->
 *        TotGalaxies[n] (int); and the number of galaxies on each tree
 *        on the current snapshot -> TreeNgals[n] (int*Ntrees).
 *
 *        If USE_MEMORY_TO_MINIMIZE_IO ON, finalize_galaxy_file also copies
 *        the galaxy data stored in pointers into the output files, so that
 *        it is all done in one go for a given file. This is done using either
 *        write_galaxy_data_snap (for SNAP output), write_all_galaxy_data
 *        (for GALAXYTREE option) or write_galaxy_for_momaf (for MOMAF option).
 *
 *        If UPDATETYPE2 is defined, the positions of type 2 galaxies
 *        (satellites without a dark matter halo) will be updated before
 *        output to contain the subsequent positions of the most bound dark
 *        matter particle at disruption time (using get_coordinates).
 *        */


void create_galaxy_files(int filenr)
{
  // create output files - snapshot option
  int n, i;
  char buf[1000];

  for(n = 0; n < NOUT; n++)
    {
      for(i = 0; i < Ntrees; i++)
    	TreeNgals[n][i] = 0;

      #ifndef NO_SAVE_GALAXY_OUT_SNP
      sprintf(buf, "%s/%s_z%1.2f_%d", OutputDir, FileNameGalaxies, ZZ[ListOutputSnaps[n]], filenr);
      if(!(FdGalDumps[n] = fopen(buf, "w+")))
        {
    	  char sbuf[1000];
    	  sprintf(sbuf, "can't open file `%s'\n", buf);
    	  terminate(sbuf);
        }


      fseek(FdGalDumps[n], (2 + Ntrees) * sizeof(int), SEEK_SET);	/* skip the space for the header */
      #endif
      TotGalaxies[n] = 0;
    }
}

void close_galaxy_files(void)
{
  int n;

  #ifndef NO_SAVE_GALAXY_OUT_SNP
  for(n = 0; n < NOUT; n++)
    {
      fseek(FdGalDumps[n], 0, SEEK_SET);
      myfwrite(&Ntrees, sizeof(int), 1, FdGalDumps[n]);	//Number of trees
      myfwrite(&TotGalaxies[n], sizeof(int), 1, FdGalDumps[n]);	// total number of galaxies
      myfwrite(TreeNgals[n], sizeof(int), Ntrees, FdGalDumps[n]);	// Number of galaxies in each tree
      fclose(FdGalDumps[n]);
    }
  #endif
}


/**@brief Saves the Galaxy_Output structure for all the galaxies in
 *        the current tree into the current output file (one for each
 *        input dark matter file) for the chosen snapshots.
 *
 *        If UPDATETYPETWO=1 then the positions and velocities of type 2
 *        galaxies are updated from the most bound dark matter particle.
 *        After that the GALAXY_OUTPUT structure is created and written.
 *        input: int file number (current file where the output is
 *        being written), int tree number (tree being currently treated).
 *
 *        If USE_MEMORY_TO_MINIMIZE_IO ON, this write statements in this
 *        routine copy the galaxy data from the working structures into
 *        pointers until that has been done for all the tree in a given file.
 */
void save_galaxy_append(int tree, int i, int n)
{
  struct GALAXY_OUTPUT galaxy_output;

  prepare_galaxy_for_output(n, &HaloGal[i], &galaxy_output);
  #ifndef NO_SAVE_GALAXY_OUT_SNP
  myfwrite(&galaxy_output, sizeof(struct GALAXY_OUTPUT), 1, FdGalDumps[n]);
  #endif

  TotGalaxies[n]++;		//this will be written later
  TreeNgals[n][tree]++;		//this will be written later (Number of galaxies in each tree)
}


 /*@brief Copies all the relevant properties from the Galaxy structure
        into the Galaxy output structure, some units are corrected.*/
void prepare_galaxy_for_output(int n, struct GALAXY *g, struct GALAXY_OUTPUT *o)
{
  int j,ibin;
  float ma, macca, maccb , mm , mdisec, mdimmi;
  float  rm , rdisec, rdimmi;
  float mbhhot, mbhgw;
  float mbhtde ,mrestde;

  o->Type = g->Type;
  o->SnapNum = g->SnapNum;
  o->CentralMvir = get_virial_mass(Halo[g->HaloNr].FirstHaloInFOFgroup);
  o->CentralRvir = get_virial_radius(Halo[g->HaloNr].FirstHaloInFOFgroup);
  o->Mvir = g->Mvir;
  o->Rvir = g->Rvir;
  o->Vvir = g->Vvir;

  for(j = 0; j < 3; j++)
    {
      o->Pos[j] = g->Pos[j];
      o->DistanceToCentralGal[j] = wrap(Halo[Halo[g->HaloNr].FirstHaloInFOFgroup].Pos[j] - g->Pos[j],BoxSize);
    }

  o->ColdGas = g->ColdGas;
  o->BulgeMass = g->BulgeMass;
  o->DiskMass = g->DiskMass;
  o->HotGas = g->HotGas;
  o->BlackHoleMass = g->BlHo[0].BlackHoleMass[STEPS-1];   //o->BlackHoleMass = g->BlackHoleMass;
  
  #ifdef EXTENDED_BH_OUTPUT
  o->HeapIndex = g->HeapIndex; //tracer
  o->BlackHoleMass_second = g->BlHo[1].BlackHoleMass[STEPS-1];   //o->BlackHoleMassB = g->BlackHoleMassB;
  /* Black Hole Properties */
  o->QuasarLum = (g->BlHo[0].QuasarLum[STEPS-1] + g->BlHo[1].QuasarLum[STEPS-1]);

  #ifdef RADIOLUMBHS  
  o->RadioLum_TD = (g->BlHo[0].RadioLum_TD[STEPS-1] + g->BlHo[1].RadioLum_TD[STEPS-1])/Hubble_h;
  o->RadioLum_ADAF = (g->BlHo[0].RadioLum_ADAF[STEPS-1] + g->BlHo[1].RadioLum_ADAF[STEPS-1])/Hubble_h;
  #endif
  o->fEDD = g->BlHo[0].fEDD[STEPS-1];
  o->eta_BH = g->BlHo[0].eta_BH[STEPS-1];
  o->spin = g->BlHo[0].spin[STEPS-1];
  ma=0.0; macca =0.0; maccb =0.0; mm =0.0; mdisec=0.0; mdimmi=0.0; 
  //rdimmi =0.0; rdisec=0.0; rm=0.0;
  mbhtde=0.0; mrestde=0.0; 
  mbhhot=0.0; mbhgw=0.0;
  for(j=0; j<STEPS; j++)     // sum all steps together for mass in snap
   {ma += g->BlHo[0].Mass_to_acc[j];
    macca +=  g->BlHo[0].MBH_acc_first_phase[j]; 
    maccb +=  g->BlHo[0].MBH_acc_second_phase[j];  
    mm +=  g->BlHo[0].MBH_acc_Mergers[j];  
    //rm +=  g->BlHo[0].M_res_Mergers[j];  
    mdisec +=  g->BlHo[0].MBH_acc_DI_sec[j];  
    //rdisec +=  g->BlHo[0].M_res_DI_sec[j];  
    mdimmi +=  g->BlHo[0].MBH_acc_DI_mmi[j];  
    //rdimmi +=  g->BlHo[0].M_res_DI_mmi[j];    
    mbhhot += g->BlHo[0].MBH_acc_radiomode[j]; //hot gas accretion
    mbhgw  += g->BlHo[0].MBH_loss_GW[j];             // mass lost in gravitational waves during BH-BH merger
    #ifdef FEED_ON_STARS
    mbhtde += g->BlHo[0].MBH_acc_TDE[j];
    //mrestde+= g->BlHo[0].M_res_TDE[j];
    #endif
   }
  o->M_res_0 = g->BlHo[0].M_res_0;             // reservoir at each phase
  o->M_dot_acc = g->BlHo[0].M_dot_acc;           // accretion rate of the black hole
  o->t_acc_quiescent_phase = g->BlHo[0].t_acc_quiescent_phase;   // time spent in the quiescent phase 
  /*
  o->MBH_0 = g->BlHo[0].MBH_0 ;                  //temporar parameter
  o->MBH_peak = g->BlHo[0].MBH_peak;             // for the second phase
  o->MBH_second_phase = g->BlHo[0].MBH_second_phase;    // for the second phase
  */
  o->Mass_to_acc = ma;  // Reservoir mass around the BH                 // accretion efficiency of the BH
  o->MBH_acc_first_phase = macca ;    // for the first phase
  o->MBH_acc_second_phase = maccb;    // for the second phase
  o->MBH_acc_radiomode = g->BlHo[0].MBH_acc_radiomode[STEPS-1]; //hot gas accretion
  o->MBH_loss_GW = g->BlHo[0].MBH_loss_GW[STEPS-1];             // mass lost in gravitational waves during BH-BH merger
  o->MBH_acc_Mergers = mm ;    // Mass accreted by the BH during the mergers (major and minor)
  o->MBH_acc_DI_sec = mdisec;        // Mass accreted by the BH during the DI secular process
  o->MBH_acc_DI_mmi = mdimmi;  // Mass accreted by the BH during the DI induced by minor merger
  // o->M_res_Mergers = rm;    // Mass accreted during the mergers (major and minor)
  // o->M_res_DI_sec = rdisec ;        // Mass accreted during the DI secular process
  // o->M_res_DI_mmi = rdimmi ; 
  //   o->Mass_to_acc = g->BlHo[0].Mass_to_acc;  // Reservoir mass around the BH                 // accretion efficiency of the BH
  //  o->MBH_acc_first_phase = g->BlHo[0].MBH_acc_first_phase ;    // for the first phase
  //  o->MBH_acc_second_phase = g->BlHo[0].MBH_acc_second_phase;    // for the second phase
  //  o->MBH_acc_radiomode = g->BlHo[0].MBH_acc_radiomode[STEPS-1]; //hot gas accretion
  //  o->MBH_loss_GW = g->BlHo[0].MBH_loss_GW[STEPS-1];             // mass lost in gravitational waves during BH-BH merger
  //  o->MBH_acc_Mergers =  g->BlHo[0].MBH_acc_Mergers ;    // Mass accreted by the BH during the mergers (major and minor)
  //  o->MBH_acc_DI_sec =  g->BlHo[0].MBH_acc_DI_sec;        // Mass accreted by the BH during the DI secular process
  //  o->MBH_acc_DI_mmi =  g->BlHo[0].MBH_acc_DI_mmi;  // Mass accreted by the BH during the DI induced by minor merger
  o->M_res_Mergers = g->BlHo[0].M_res_Mergers[STEPS-1];    // Mass accreted during the mergers (major and minor)
  o->M_res_DI_sec = g->BlHo[0].M_res_DI_sec[STEPS-1] ;        // Mass accreted during the DI secular process
  o->M_res_DI_mmi = g->BlHo[0].M_res_DI_mmi[STEPS-1]  ; 


  o->StellarMass = g->BulgeMass+g->DiskMass;
  o->Bulge_di = g->Bulge_di;
  o->Sersic_index = g->Sersic_index;
  //o-> tBH_merger_clock_dt = g-> BlHo[0].tBH_merger_clock_dt  // clock for the BH to merger with its companion
  o->tBH_merger = g->BlHo[0].tBH_merger;               // time needed for the BH to merger with its companion
  // #ifdef MODELWANDERINGBHS
  // o->t_spent_WPHs = g->BlHo[0].t_spent_WPHs ;     // time spent in the wandering phase: -1 means never was in a WPHs
  // o->t_ejected_WPHs = g->BlHo[0].t_spent_WPHs;   // redshfit in which the BH was ejected: -1 means never was in WPHs
  // o->MBH_WPH = g->BlHo[0].MBH_WPH;       // Mass of the black hole when it was ejected form the galaxy: -1 means no WPH
  // #endif
  //Galaxy properties to store along the BH characteristics 
  o->t_last_DI = g->t_last_DI;
  o->t_last_merger = g->t_last_merger;
  o->t_last_MINOR = g->t_last_MINOR;
  o->t_last_MAJOR = g->t_last_MAJOR;
  o->t_last_insta_minor_induced = g->t_last_insta_minor_induced;
  o->Type_merger= g->Type_merger;           //type of galaxy at merger        
  o->n_Min_Mer = g->n_Min_Mer;                 // Number of the minor mergers in a given SNP (cumulative along the merger tree of that SNP)
  o->n_Maj_Mer = g->n_Maj_Mer;                 // Number of the major mergers in a given SNP (cumulative along the merger tree of that SNP)
  o->n_di = g->n_di;    // Number of disk instabilities per SNP
  o->n_di_MM = g->n_di_MM;
  o->n_di_SP = g->n_di_SP;

  #ifdef FEED_ON_STARS
  o->MBH_acc_TDE = mbhtde ; 
  o->M_res_TDE = g->BlHo[0].M_res_TDE[STEPS-1] ; 
  o->r_inf = g->r_inf[STEPS-1];
  o->sigma_inf = g->sigma_inf[STEPS-1];
  o->rho_inf = g->rho_inf[STEPS-1];
  o->ICM  = g->ICM;
  #endif
    
  #endif //EXTENDED_BH_OUTPUT ended

  #ifdef FEED_ON_STARS
//  o->T_relax = g->BlHo[0].T_relax[STEPS-1]; //simplified 2 variables here
// o->M_core = g->BlHo[0].M_core[STEPS-1];//+g->BlHo[1].M_core[STEPS-1]; //the second one shoud be zero
  o->M_dot_TDE = g->BlHo[0].M_dot_TDE[STEPS-1];
  o->t_TD_phase = g->BlHo[0].t_TD_phase;
  #endif

  /* initialization variables */
  o->SnapOfInit = g->SnapOfInit;
  /* IGM metallicity variables */
  #ifdef SAVE_METALSHELLRADIUS
  for(j=0; j<MAXSNAPS; j++)
    o->MetalShellRadius[j] = g->MetalShellRadius[j];
  #endif
  o->LyWerLuminosity = g->LyWerLuminosity;
  o->LyWerStars = g->LyWerStars;
  o->LyWerSfr = g->LyWerSfr;

  //o->Pristine = g->Pristine;

#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
#ifdef OUTPUT_REST_MAGS 
  /* Luminosities are converted into Mags in various bands */
  for(j = 0; j < NMAG; j++)
 	o->Mag[j] = lum_to_mag(g->Lum[j][n]);
  o->MagDust[j] = lum_to_mag(g->LumDust[j][n]);
#endif
#endif //ndef POST_PROCESS_MAGS
#endif //COMPUTE_SPECPHOT_PROPERTIES

#ifndef LIGHT_OUTPUT
  
#ifdef GALAXYTREE
  o->HaloID = HaloIDs[g->HaloNr].HaloID;
  o->Redshift = ZZ[g->SnapNum];

  int ii = (int) floor(o->Pos[0] * ScaleFactor);
  int jj = (int) floor(o->Pos[1] * ScaleFactor);
  int kk = (int) floor(o->Pos[2] * ScaleFactor);

  o->PeanoKey = peano_hilbert_key(ii, jj, kk, Hashbits);

  o->SubID = calc_big_db_subid_index(g->SnapNum, Halo[g->HaloNr].FileNr, Halo[g->HaloNr].SubhaloIndex);

  int tmpfirst = Halo[g->HaloNr].FirstHaloInFOFgroup;
  int lenmax = 0;
  int next = tmpfirst;
  while(next != -1)
    {
      if(Halo[next].Len > lenmax)
	{
	  lenmax = Halo[next].Len;
	  tmpfirst = next;
	}
      next = Halo[next].NextHaloInFOFgroup;
    }

  o->MMSubID = calc_big_db_subid_index(g->SnapNum, Halo[tmpfirst].FileNr, Halo[tmpfirst].SubhaloIndex);
#endif

  o->LookBackTimeToSnap = NumToTime(g->SnapNum)*UnitTime_in_years/Hubble_h;
  o->InfallVmax = g->InfallVmax;
  o->InfallVmaxPeak = g->InfallVmaxPeak;
  o->InfallSnap = g->InfallSnap;
  o-> InfallHotGas = g-> InfallHotGas;
  o->HotRadius =  g->HotRadius;
#ifdef HALOPROPERTIES
  o->HaloM_Mean200 = g->HaloM_Mean200;
  o->HaloM_Crit200 = g->HaloM_Crit200;
  o->HaloM_TopHat = g->HaloM_TopHat;
  o->HaloVelDisp = g->HaloVelDisp;
  o->HaloVmax = g->HaloVmax;
#endif

  o->Len = g->Len;
  o->Vmax = g->Vmax;

  o->BulgeSize = g->BulgeSize;
  o->EjectedMass = CORRECTDBFLOAT(g->EjectedMass);
  o->ICM = g->ICM;

  for(j = 0; j < 3; j++)
    {
      o->Vel[j] = g->Vel[j];
#ifdef HALOSPIN
      o->HaloSpin[j] = g->HaloSpin[j];
#endif
      o->GasSpin[j] = g->GasSpin[j];
      o->StellarSpin[j] = g->StellarSpin[j];
#ifdef HALOPROPERTIES
      o->HaloPos[j] = g->HaloPos[j];
      o->HaloVel[j] = g->HaloVel[j];
      o->HaloSpin[j] = g->HaloSpin[j];
#endif      
    }

  o->XrayLum = g->XrayLum;
  o->GasDiskRadius = g->GasDiskRadius;
  o->StellarDiskRadius = g->StellarDiskRadius;
  o->CoolingRadius = g->CoolingRadius;

  o->QuasarAccretionRate = g->QuasarAccretionRate * UNITMASS_IN_G / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS;
  o->RadioAccretionRate = g->RadioAccretionRate * UNITMASS_IN_G / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS;
  o->CosInclination = g->CosInclination;

  if(g->Type == 2 || (g->Type == 1 && g->MergeOn == 1)) {
    o->OriMergTime=g->OriMergTime*UnitTime_in_years/Hubble_h;
    o->MergTime = g->MergTime*UnitTime_in_years/Hubble_h;
  }
  else {
    o->OriMergTime=0.0;
    o->MergTime = 0.0;
  }

#ifndef GALAXYTREE
  o->HaloIndex = g->HaloNr;
#endif
#ifdef MBPID
  o->MostBoundID = g->MostBoundID;
#endif


#ifdef GALAXYTREE
  o->DisruptOn = g->DisruptOn;
#endif
  o->MergeOn = g->MergeOn;

//METALS
#ifndef   DETAILED_METALS_AND_MASS_RETURN
  o->MetalsColdGas = CORRECTDBFLOAT(g->MetalsColdGas);
  o->MetalsStellarMass = CORRECTDBFLOAT(g->MetalsDiskMass)+ CORRECTDBFLOAT(g->MetalsBulgeMass);
  o->MetalsDiskMass = CORRECTDBFLOAT(g->MetalsDiskMass);
  o->MetalsBulgeMass = CORRECTDBFLOAT(g->MetalsBulgeMass);
  o->MetalsHotGas = CORRECTDBFLOAT(g->MetalsHotGas);
  o->MetalsEjectedMass = CORRECTDBFLOAT(g->MetalsEjectedMass);
  o->MetalsICM = CORRECTDBFLOAT(g->MetalsICM);
#ifdef METALS_SELF
  o->MetalsHotGasSelf = CORRECTDBFLOAT(g->MetalsHotGasSelf);
#endif
#else
  o->MetalsColdGas = g->MetalsColdGas;
  o->MetalsDiskMass = g->MetalsDiskMass;
  o->MetalsBulgeMass = g->MetalsBulgeMass;
  o->MetalsHotGas = g->MetalsHotGas;
  o->MetalsEjectedMass = g->MetalsEjectedMass;
  o->MetalsICM = g->MetalsICM;
#ifdef METALS_SELF
  o->MetalsHotGasSelf = g->MetalsHotGasSelf;
#endif
#endif

#ifdef TRACK_BURST
  o->BurstMass=g->BurstMass;
#endif

//STAR FORMATION HISTORIES / RATES

#ifdef STAR_FORMATION_HISTORY
  o->sfh_ibin=g->sfh_ibin;
  ibin=0;
  for (j=0;j<=o->sfh_ibin;j++) {
 	  o->sfh_DiskMass[j]=g->sfh_DiskMass[j];
 	  o->sfh_BulgeMass[j]=g->sfh_BulgeMass[j];
 	  o->sfh_ICM[j]=g->sfh_ICM[j];
 	  o->sfh_MetalsDiskMass[j]=g->sfh_MetalsDiskMass[j];
 	  o->sfh_MetalsBulgeMass[j]=g->sfh_MetalsBulgeMass[j];
 	  o->sfh_MetalsICM[j]=g->sfh_MetalsICM[j];
#ifdef INDIVIDUAL_ELEMENTS
	  o->sfh_ElementsDiskMass[j]=g->sfh_ElementsDiskMass[j];
	  o->sfh_ElementsBulgeMass[j]=g->sfh_ElementsBulgeMass[j];
	  o->sfh_ElementsICM[j]=g->sfh_ElementsICM[j];
#endif
#ifdef TRACK_BURST
	  o->sfh_BurstMass[j]=g->sfh_BurstMass[j];
#endif
   }

  //Set all non-used array elements to zero:
  // important if we want to read files in database that all values are valid SQLServer floats
  for (j=o->sfh_ibin+1;j<SFH_NBIN;j++) {
	  o->sfh_DiskMass[j]=0.;
	  o->sfh_BulgeMass[j]=0.;
	  o->sfh_ICM[j]=0.;
	  o->sfh_MetalsDiskMass[j]=metals_init();
	  o->sfh_MetalsBulgeMass[j]=metals_init();
	  o->sfh_MetalsICM[j]=metals_init();
#ifdef INDIVIDUAL_ELEMENTS
	  o->sfh_ElementsDiskMass[j]=elements_init();
	  o->sfh_ElementsBulgeMass[j]=elements_init();
	  o->sfh_ElementsICM[j]=elements_init();
#endif
#ifdef TRACK_BURST
	  o->sfh_BurstMass[j]=0.;
#endif
  }
#endif //STAR_FORMATION_HISTORY

#ifdef INDIVIDUAL_ELEMENTS
  o->DiskMass_elements = g->DiskMass_elements;
  o->BulgeMass_elements = g->BulgeMass_elements;
  o->ColdGas_elements = g->ColdGas_elements;
  o->HotGas_elements = g->HotGas_elements;
  o->EjectedMass_elements = g->EjectedMass_elements;
  o->ICM_elements = g->ICM_elements;
#endif

  o->PrimordialAccretionRate = CORRECTDBFLOAT(g->PrimordialAccretionRate * UNITMASS_IN_G / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS);
  o->CoolingRate = CORRECTDBFLOAT(g->CoolingRate * UNITMASS_IN_G / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS);
  o->CoolingRate_beforeAGN = CORRECTDBFLOAT(g->CoolingRate_beforeAGN * UNITMASS_IN_G / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS);

 //NOTE: in Msun/yr
  o->Sfr = CORRECTDBFLOAT(g->Sfr * UNITMASS_IN_G / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS);
  o->SfrBulge = CORRECTDBFLOAT(g->SfrBulge * UNITMASS_IN_G / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS);

//MAGNITUDES
#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifdef POST_PROCESS_MAGS
    //Convert recorded star formation histories into mags
    post_process_spec_mags(o);
#else //ndef POST_PROCESS_MAGS

#ifdef OUTPUT_REST_MAGS
  // Luminosities are converted into Mags in various bands
  for(j = 0; j < NMAG; j++)
    {
	  //o->Mag[j] = lum_to_mag(g->Lum[j][n]); -> DONE ON TOP FOR LIGHT_OUTPUT AS WELL
	  o->MagBulge[j] = lum_to_mag(g->LumBulge[j][n]);
	  o->MagDust[j] = lum_to_mag(g->LumDust[j][n]);
#ifdef ICL
	  o->MagICL[j] = lum_to_mag(g->ICLLum[j][n]);
#endif
    }
 

#endif //OUTPUT_REST_MAGS
#ifdef OUTPUT_OBS_MAGS
#ifdef COMPUTE_OBS_MAGS
  // Luminosities in various bands
  for(j = 0; j < NMAG; j++)
    {
      o->ObsMag[j] = lum_to_mag(g->ObsLum[j][n]);
      o->ObsMagBulge[j] = lum_to_mag(g->ObsLumBulge[j][n]);
      o->ObsMagDust[j] = lum_to_mag(g->ObsLumDust[j][n]);
#ifdef ICL
      o->ObsMagICL[j] = lum_to_mag(g->ObsICL[j][n]);
#endif

#ifdef OUTPUT_MOMAF_INPUTS
      o->dObsMag[j] = lum_to_mag(g->dObsLum[j][n]);
      o->dObsMagBulge[j] = lum_to_mag(g->dObsLumBulge[j][n]);
      o->dObsMagDust[j] = lum_to_mag(g->dObsLumDust[j][n]);
#ifdef ICL
      o->dObsMagICL[j] = lum_to_mag(g->dObsICL[j][n]);
#endif
#endif
    }
#endif //COMPUTE_OBS_MAGS
#endif //OUTPUT_OBS_MAGS
#endif //ndef POST_PROCESS_MAGS
#endif //COMPUTE_SPECPHOT_PROPERTIES

  if((g->DiskMass+g->BulgeMass)> 0.0)
    {
      o->MassWeightAge = g->MassWeightAge[n] / (g->DiskMass+g->BulgeMass);
      o->MassWeightAge = o->MassWeightAge / 1000. * UnitTime_in_Megayears / Hubble_h;	//Age in Gyr
    }
  else
    o->MassWeightAge = 0.;

#ifdef FIX_OUTPUT_UNITS
  fix_units_for_ouput(o);
#endif
#endif //ndef LIGHT_OUTPUT
}


#ifdef FIX_OUTPUT_UNITS

#ifdef LIGHT_OUTPUT
/**@brief Removes h from units of galaxy properties. If desired (makefile option
 * FIX_OUTPUT_UNITS is set), the output properties of the galaxies can be scaled
 * to physical units excluding any factors of h == Hubble/100 km/s/Mpc. */
void fix_units_for_ouput(struct GALAXY_OUTPUT *o)
{
  o->Pos[0] /= Hubble_h;
  o->Pos[1] /= Hubble_h;
  o->Pos[2] /= Hubble_h;
  o->Mvir /= Hubble_h;
  o->Rvir /= Hubble_h;
  o->ColdGas /= Hubble_h;
  o->DiskMass /= Hubble_h;
  o->BulgeMass /= Hubble_h;
  o->HotGas /= Hubble_h;
  o->BlackHoleMass /= Hubble_h;

}
#else

/**@brief Removes h from units of galaxy properties. If desired (makefile option
 * FIX_OUTPUT_UNITS is set), the output properties of the galaxies can be scaled
 * to physical units excluding any factors of h == Hubble/100 km/s/Mpc. */
void fix_units_for_ouput(struct GALAXY_OUTPUT *o)
{
  int j;

  o->Pos[0] /= Hubble_h;
  o->Pos[1] /= Hubble_h;
  o->Pos[2] /= Hubble_h;
  o->CentralMvir /= Hubble_h;
  o->Mvir /= Hubble_h;
  o->Rvir /= Hubble_h;
#ifdef HALOPROPERTIES
  o->HaloM_Mean200 /= Hubble_h;
  o->HaloM_Crit200 /= Hubble_h;
  o->HaloM_TopHat /= Hubble_h;
  o->HaloPos[0] /= Hubble_h;
  o->HaloPos[1] /= Hubble_h;
  o->HaloPos[2] /= Hubble_h;
  o->HaloSpin[0] /= Hubble_h;
  o->HaloSpin[1] /= Hubble_h;
  o->HaloSpin[2] /= Hubble_h;
#endif
  o->GasSpin[0] /= Hubble_h;
  o->GasSpin[1] /= Hubble_h;
  o->GasSpin[2] /= Hubble_h;
  o->StellarSpin[0] /= Hubble_h;
  o->StellarSpin[1] /= Hubble_h;
  o->StellarSpin[2] /= Hubble_h;
  o->HotRadius /= Hubble_h;
  o->ColdGas /= Hubble_h;
  o->DiskMass /= Hubble_h;
  o->BulgeMass /= Hubble_h;
  o->HotGas /= Hubble_h;
  o->EjectedMass /= Hubble_h;
  o->BlackHoleMass /= Hubble_h;
  #ifdef RADIOLUMBHS 
	o->RadioLum_TD /= pow(Hubble_h,2);
  o->RadioLum_ADAF /= pow(Hubble_h,2);
	#endif
  o->ICM /= Hubble_h;
  o->BulgeSize /= Hubble_h;
  o->StellarDiskRadius /= Hubble_h;
  o->GasDiskRadius /= Hubble_h;
  o->CoolingRadius /= sqrt(Hubble_h);
#ifdef TRACK_BURST
  o->BurstMass /= Hubble_h;
#endif

  o->MetalsColdGas=metals_add(metals_init(),o->MetalsColdGas,1./Hubble_h);
  o->MetalsDiskMass=metals_add(metals_init(),o->MetalsDiskMass,1./Hubble_h);
  o->MetalsBulgeMass=metals_add(metals_init(),o->MetalsBulgeMass,1./Hubble_h);
  o->MetalsHotGas=metals_add(metals_init(),o->MetalsHotGas,1./Hubble_h);
  o->MetalsEjectedMass=metals_add(metals_init(),o->MetalsEjectedMass,1./Hubble_h);
  o->MetalsICM=metals_add(metals_init(),o->MetalsICM,1./Hubble_h);

#ifdef STAR_FORMATION_HISTORY
  for(j=0;j<=o->sfh_ibin;j++) {
    o->sfh_DiskMass[j] /= Hubble_h;
    o->sfh_BulgeMass[j] /= Hubble_h;
    o->sfh_ICM[j] /= Hubble_h;
#ifdef TRACK_BURST
    o->sfh_BurstMass[j] /= Hubble_h;
#endif
    o->sfh_MetalsDiskMass[j]=metals_add(metals_init(),
					   o->sfh_MetalsDiskMass[j],1./Hubble_h);
    o->sfh_MetalsBulgeMass[j]=metals_add(metals_init(),
					 o->sfh_MetalsBulgeMass[j],1./Hubble_h);
    o->sfh_MetalsICM[j]=metals_add(metals_init(),
					 o->sfh_MetalsICM[j],1./Hubble_h);
  }
#endif
}
#endif // NO Light_Output

#endif // FIX_OUTPUT_UNITS

long long calc_big_db_offset(int filenr, int treenr)
{
  long long big;
#ifdef MRII
  big = (((filenr * (long long) 1000000) + treenr) * (long long) 1000000000);
#else
  big = (((filenr * (long long) 1000000) + treenr) * (long long) 1000000);
#endif
  return big;
}

long long calc_big_db_subid_index(int snapnum, int filenr, int subhaloindex)
{
  long long big;
#ifdef MRII
  big = snapnum * (long long) 10000000000000 + filenr * (long long) 1000000000 + subhaloindex;
#else
  big = snapnum * (long long) 1000000000000 + filenr * (long long) 100000000 + subhaloindex;
#endif
  return big;
}


#ifdef STAR_FORMATION_HISTORY
void write_sfh_bins()
{
  FILE* SFH_Bins_File;
  struct SFH_Time *sfh_times;
  char buf[1000];

  int snap,j, nbins,ibin;
  nbins = 0;
  for(snap = 0; snap < MAXSNAPS; snap++)
    for(j=0;j < SFH_ibin[snap][0];j++)
      nbins++;


  sfh_times = (struct SFH_Time *) mymalloc("sfh_times", sizeof(struct SFH_Time) * nbins);
  ibin = 0;
  for(snap = 0; snap < MAXSNAPS; snap++)
    {
      for(j=0;j < SFH_ibin[snap][0];j++)
	{
	  sfh_times[ibin].snapnum = snap;
	  sfh_times[ibin].bin = j;
	  sfh_times[ibin].lookbacktime = (SFH_t[snap][0][j]+SFH_dt[snap][0][j]/2.-NumToTime(snap))*UnitTime_in_years/Hubble_h;
	  sfh_times[ibin].dt=SFH_dt[snap][0][j]*UnitTime_in_years/Hubble_h;
	  sfh_times[ibin].nbins = SFH_Nbins[snap][0][j];
	  ibin++;
	}
    }

  sprintf(buf, "%s/SFH_Bins", FinalOutputDir);
  if(!(SFH_Bins_File = fopen(buf, "w")))
    {
      char sbuf[1000];
      sprintf(sbuf, "can't open file `%s'\n", buf);
      terminate(sbuf);
    }
  else
      printf("writing sfh bins to %s\n", buf);

  // write # bins
  myfwrite(&nbins, sizeof(int), 1, SFH_Bins_File);	// write 1
  myfwrite(sfh_times, sizeof(struct SFH_Time), nbins, SFH_Bins_File);	// size of an output structure (Galaxy_Output)
  fflush(SFH_Bins_File);
  fclose(SFH_Bins_File);


}
#endif

