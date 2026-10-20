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

/*  Created in: 2021
 *      Author: Daniele Spinoso */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"


/** @file model_igm_metallicity.c
 *  @brief This file includes all the routines to deal with the metallicity of the IGM,
 *         which basicaly means:
 *	   1 - associating a metal content to the primordial infall on newly-resolved
 *	       halos (performed in "compute_infall_metalMass" routine). This was
 *             introduced in order to avoid having halos forming at zero-metallicity
 *             as soon as they are resolved (even at very high masses and low redshift,
 *             when they should be already polluted by internal star-formation or local
 *             pollution) and hence to obtain reasonable initial conditions for their
 *             metallic content to which the BH-seeding prescriptions can be anchored.
 *	   2 - computing the local metallic pollution at the position of newly-resolved
 *	       halos, due to SNe feedback of neighboring star-forming halos (this is
 *	       done in "check_metal_pollution" routine)
 *	   All the additional routines in this file are necessary to these two main tasks.
 *
 *
 *
 * @brief Computation of the average IGM metallicity (Z_IGM -> ZIGM variable in the code).
 *         The (average) IGM metallicity as a function of redshift is obtained simply as:
 *         
 *         ZIGM(z) = Total_Metals_Mass_inThe_IGM(z) / Total_Gas_Mass_inThe_IGM(z)
 *
 *         where:
 *         Both the numerator and denominator of ZIGM(z) are obtained by using a run of
 *         the Henriques+2015 LGalaxies version, as follows:
 *         1) Total_Metals_Mass_inThe_IGM(z) is the sum of the mass of metals ejected by
 *            all the galaxies in the whole volume of the simulated box at a given z.
 *            The mass of ejected metals is already computed by LGalaxies and stored in
 *            "EjectedMetalsMass" for each GALAXY or GALAXY_OUTPUT structure. NOTE: the
 *            variable "EjectedMetalsMass" already accounts for the balance between
 *            ejected and reincorporated metals mass, so its value really provides the
 *            mass of metals which is outside the Rvir of each galaxy at each redshift.
 *         2) Total_Gas_Mass_inThe_IGM(z) is the sum of all the baryonic mass which is
 *            not collapsed into haloes, in the whole box, at each redshift. Its value
 *            is obtained as the difference between the total baryonic mass inside the
 *            simulated box, i.e.: fb * particle_DarkMatter_mass * N_particles (with
 *            fb=0.15 being the baryon fraction used by Henriques+2015), and the total
 *            baryonic mass collapsed into halos at each z, i.e. fb * total_DM_inHalos
 *            (with total_DM_inHalos being the sum of the "Mvir" variable of all galaxy
 *            structures in the whole simulated box at each redshift).
 *         To summarize:
 *                              sum[ EjectedMetalsMass(z) ]_box
 *         ZIGM(z) = ___________________________________________________________
 *	               fb * {particle_mass * N_particles - sum[ Mvir(z) ]_box}
 *
 *         Finally, once ZIGM(z) is obtained both for Millennium and Millennium-II at each z,
 *         its values are stored in dedicated tables (one per simulation) inside the folder:
 *         "./ZigmTables/". In particular, the stored values are: Log10( ZIGM(z)/ZSUN ),
 *         where ZSUN is the solar metallicity: ZSUN=0.0177.
 *
 *
 *
 * @brief Computation of the local variations of metallic pollution in the IGM.
 *        The local variations of any background (e.g. radiative or chemical) depend on the
 *        contributions coming from potentially ANY point of the simulated box (especially
 *        in the case of radiative backgrounds with large horizon length). Unfortunately,
 *        since its input tree-files are divided into subvolumes, LGalaxies cannot account
 *        for all the possible contributions to a local field ON THE RUN. In other words,
 *        while the code is running (i.e. it is building and evolving a given merger tree),
 *        it has no info about ALL the halos/galaxies living in the WHOLE simulated box,
 *        outside the FoF group which is currently being analyzed.
 *        To overcome this, we excluded the use of pre-computed grids describing the spatial
 *        variations of a given background in order to avoid: 1) being limited by the grids
 *        resolution, and 2) intense memory usage imposed by the need of having grid per each
 *        SnapNum/redshift during the code run.
 *        Instead, we relied on the following method to check the local metallic pollution at
 *        the position of any newly-initialized galaxy in the box:
 *        1) track in time the expansion of "metallic shells" launched by each actively
 *           star-forming galaxy in the box,
 *        2) check the self-pollution, i.e. the chemical enrichment due to (eventual) internal
 *           star formation of the newly-initialized galaxy.
 *        3) check the FoF-pollution, i.e. the contribution to the local chemical enrichment
 *           provided by the metallic shells of all the galaxies belonging to the FoF group
 *           of the newly-initialized galaxy.
 *        4) check for eventual "external-FoF" contributions to the local pollution by using
 *           a pre-saved list of potential "metal providers" or "metal sources", defined over
 *           the whole box.
 *        Steps 1, 2 and 3 can be easily performed on the run. Instead, to define the external
 *        list of metal-sources we use a preliminary run of the code (with STORE_SOURCES option
 *        active in the "My_Makefile_options") and select as potential "metal sources" all the
 *        galaxies with very extended metal shells (as compared to the galaxies' virial radius,
 *        see note 2, below, and the "store_metal_source" routine). We then store some info
 *        about the metallic shells of these "metal sources" in an external binary file.
 *        NOTES:
 *        1) tracking the expansion of metallic shells is fairly easy, since the variables
 *           "EjectedMass" and "MetalsEjectedMass" already take into account the ejection
 *           of both gas and metals beyond the virial radius of each halo, due to intense
 *           SF episodes. The metallic shell evolution is tracked in "update_metal_shell"
 *           and detailed in the comment above that routine.
 *        2) The "self-pollution" is is likely null during the first substeps after each
 *           galaxy initialization since the primordial infall didn't have time to cool
 *           and form stars yet. Nevertheless, we do not check for metal pollution only
 *           at the substep of initialization, hence it's actually possible for the self
 *           pollution to actively play a role in the local metal enrichment.
 *        2) the potential "metal sources" are galaxies which were able to launch a fairly
 *           big metallic shell and significantly influence their surroundings. The exact
 *           definition of "failry big" metallic shell is controlled by an input parameter
 *           and it's performed in the routine "store_metal_source", below.
 *        3) there are 3 different types of "metal contributors" to the local level of metal
 *           pollution: self-contributors (the newly-initialized galaxy who eventually formed
 *           stars very fast and very efficiently), the FoF-contributors (the galaxies in the
 *           current FoF group, whose metallic shells are able to reach the newly-initialized
 *           galaxy) and the out-of-FOF-contributors (those galaxies in the external list of
 *           metal sources whose metallic shell is able to reach the newly-initialized galaxy).
 *           In case the same galaxy would belong to more than one class of contributors, we
 *           make sure to account for it only once, namely at its first appearence in one of
 *           the classes of metal contributors.
 *           */


/*------  ROUTINES FOR METALLICITY OF PRIMORDIAL INFALL (from average IGM metallicity)  ------*/

/** @brief The read_ZIGM_tables() routine reads the table containing the ZIGM(z) values,
 *         i.e. the average IGM metallicity as a function of redshift. These are computed
 *         in post-processing as described in the documentation here above.
 *         The values read from the tables are stored in the global variable ZIGM[]. These
 *         values are then useed to associate a non-zero "primordial metallicity" to the
 *         primordial infall of gas into halos. */
void read_ZIGM_tables(void)
{
  FILE *frz;
  int ii, dummy;
  double zz, igm;

  if( !(frz = fopen(FileWithZIGMvsZ, "r")) )
    {
      printf("file `%s' not found.\n", FileWithZIGMvsZ);
      exit(0);
    }
  for(ii=0; ii<LastDarkMatterSnapShot+1; ii++)
    {
      dummy = fscanf(frz, "%lf", &zz);
      dummy = fscanf(frz, "%lf", &igm);
      ZIGM[ii] = igm;
    }
  fclose(frz);  
}


/** @brief The compute_infall_metalMass() routine computes the mass of metals which needs
 *         to be associated to the primordial infall onto halos, at the redshift of their
 *         initialization.  This is done by just:
 *         1) extracting ZIGM(z) for the current z (stored in Log10[ZIGM(z)/ZSUN] units)
 *         2) using it to compute the actual mass of metals inside the primordial infall
 *         NOTE: this routine makes use of the GALAXY.SnapNum clock (i.e. Halo.SnapNum - 1),
 *               since that's when the infall actually happens, physically, onto the galaxy.
 * */
double compute_infall_metalMass(double infallingGas, int halonr, int centralgal)
{
  double primordialMetals = 0;
  double metalRatio = 0;

  metalRatio = ZSUN * pow(10, ZIGM[Halo[halonr].SnapNum-1]);
  primordialMetals = infallingGas * metalRatio;   // [10^10 Msun/h]. Mass of metals in the infalling gas.

  if(primordialMetals*infallingGas < 0)
    printf("Error in 'compute_infall_metal': opposite signs! infall=%f, metals=%f\n",infallingGas,primordialMetals);

  return primordialMetals;
}




/*------  ROUTINES FOR LOCAL VARIATIONS OF IGM POLLUTION (especially needed for BH seeding)  ------*/

/** @brief The following routine opens the binary file where the list of metal sources is stored
 *         (it should have been saved there in a preliminary run of the model), counts the number
 *         of items in the list, allocates the necessary memory according to this number and then
 *         stores the file content in a global array of metal sources (used for BH seeding). The
 *         latter, will then be used in a following run to compute the local pollution of eventual
 *         Black-Hole seeding candidates (if the "SEEDING" option of the Makefile is ON). The list
 *         of Metal-Sources is defined in the documentation at the beginning of this module and in
 *         the comments of "store_metal_source" routine.
 * */
void read_metal_sources(void)
{
  METALSOURCESFILE = fopen(FileWithMetalSources, "rb");
    
  if( METALSOURCESFILE != NULL )
    {
      /*---  first we count the elements in the external list of MetalSources ---*/
      /*---  in order to properly allocate the memory to store and use them.  ---*/
      METALSOURCESNR = 0;
      struct SHELLINFO source;   // dummy structure to read the external list
      
      while( fread(&source, sizeof(struct SHELLINFO), 1, METALSOURCESFILE)!=0 )
	METALSOURCESNR++;

      fclose(METALSOURCESFILE);
      printf("Metal Sources counted in model_igm_metallicity.c/read_metal_sources. Number of Metal Sources found: N=%d\n",METALSOURCESNR);
      
      /*---  now we can allocate the memory we need ---*/
      METALSOURCE = mymalloc("MetalSource", sizeof(struct SHELLINFO)*METALSOURCESNR);

      /*-- now we can read the external list again and store its content in the allocated array ---*/
      METALSOURCESFILE = fopen(FileWithMetalSources, "rb");
      int ii=0;
      while( fread(&METALSOURCE[ii], sizeof(struct SHELLINFO), 1, METALSOURCESFILE)!=0 )
	ii++;
      fclose(METALSOURCESFILE);

      /*---  now we check if the reading was fine ---*/
      /* int num; */
      /* num = 0; */
      /* print_metal_source(num); */
      /* num = METALSOURCESNR-3; */
      /* print_metal_source(num); */
      /* num = METALSOURCESNR-2; */
      /* print_metal_source(num); */
      /* num = METALSOURCESNR-1; */
      /* print_metal_source(num); */
      /* terminate("\nHope everything is fine with the numbers\n"); */
    }
  else
    terminate("\nMetal shell sources not read in model_igm_metallicity.c --> read_metal_sources\n");
}



//-----------------------------------------------------------------------------------------
//  This routine just prints a single metals source
//  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  ONCE MADE SURE THAT THE WRITING/READING OF THE METAL SOURCES IS OK, REMOVE THIS ROUTINE !!!!!!!!!!
//-----------------------------------------------------------------------------------------
void print_metal_source(int num)
{
  printf("\nHalo %d: mvir=%e [Msun]", num , METALSOURCE[num].mvir*1.e10/Hubble_h);
  printf("\nHalo %d: rvir=%e [Kpc]", num , METALSOURCE[num].rvir*1.e3/Hubble_h);
  printf("\nHalo %d: X=%e [Mpc]", num , METALSOURCE[num].xpos/Hubble_h);
  printf("\nHalo %d: Y=%e [Mpc]", num , METALSOURCE[num].ypos/Hubble_h);
  printf("\nHalo %d: Z=%e [Mpc]", num , METALSOURCE[num].zpos/Hubble_h);
  printf("\nHalo %d: Rshell=%f [Kpc]", num , METALSOURCE[num].rshell*1.e3/Hubble_h);
  printf("\nHalo %d: Metals=%e [Msun]", num , METALSOURCE[num].shellmetals*1.e10/Hubble_h);
  printf("\nHalo %d: ShellGas=%e [Msun]", num , METALSOURCE[num].shellgas*1.e10/Hubble_h);
  printf("\nHalo %d: l.b.-time/redshift=%f [int.units]", num , METALSOURCE[num].time);
  printf("\nHalo %d: HaloID=%d", num , METALSOURCE[num].haloID);
  printf("\nHalo %d: pclesNr=%d", num , METALSOURCE[num].pcles);
  printf("\nHalo %d: boundID=%llu", num , METALSOURCE[num].boundID);
  printf("\n\n");
}



/** @brief "store_metal_source" is used to save to external binary files some info about galaxies
 *         whose associated metallic shell is actively polluting their surroundings with metals,
 *         These galaxies can be of two kinds:
 *         1) Metal Sources: defined as galaxies whose metallic shell has a radius of "N" times
 *            their virial radius. The expansion of the metallic shell is powered by the ejecta
 *            from star formation (see "update_metal_shell_radius" routine below). The factor "N"
 *            used to control the definition of "metal sources" is a free-parameter given as input
 *            in the input-parameter file. These "metal sources" should be thought as "potential
 *            providers" of metals, or "potential contributors" to the local pollution at any
 *            point in the box.
 *         2) Metal Contributors: these are those galaxies which are actually contributing to the
 *            local level of chemical pollution at a specific point of the box at a given redshift.
 *            As an example, they can be galaxies from the list of "metal sources" whose metallic
 *            shell is actually able to reach the specific point at which the local amount of metal
 *            pollution is being evaluated. Another example of "metal contributors" is provided by
 *            galaxies belonging to the FoF-group of a newly-initialized galaxy, which are actively
 *            providing metals to the point in space where the new galaxy is being initialized.
 *         The following "store_metal_source" routine is in charge of saving-out both "metal sources"
 *         or "metal contributors" to external binary files. In the first case, it writes out the
 *         external list of "metal sources" defined over the whole box, in the second case it writes
 *         out the smaller list of "metal contributors" actually providing metals to a single, very
 *         specific location inside the box.
 *         NOTE:
 *         The input pointer is actually deciding whether we are writing to the external list of
 *         metal sources, or to the smaller list of metal-contributors. The logic to write out is
 *         as follows: by default, we would write out the input galaxy to the input file. Indeed,
 *         we always want to save-out the metal contributors which were specifically-selected by
 *         the routine "check_metal_pollution". For the case of metal sources, instead, we first
 *         need to check whether their shell radius satisfies a given threshold (defined in terms
 *         of the "RshellMin" input parameter). If this condition is not verified, we deactivate
 *         the writing-out of the current galaxy. This check is only performed if the option
 *         "STORE_SOURCES" of the makefile is active for the current run.
 * */
void store_metal_source(int p, int nstep, FILE * OutFilePointer)
{
  double old_ShellRadius, new_ShellRadius, ShellRadius;
  float writeOut;

  /*--- The shell radii of each galaxy are saved only at each SnapNum, so here the values are ---*/
  /*--- linearly interpolated at the current substep:    oldRad + dRad/deltat * nstep*dt      ---*/
  old_ShellRadius = Gal[p].MetalShellRadius[Gal[p].SnapNum] * Hubble_h;         // radius at this Gal[p].SnapNum, [Mpc/h]
  new_ShellRadius = Gal[p].MetalShellRadius[Gal[p].SnapNum+1] * Hubble_h;       // radius at next Gal[p].SnapNum, [Mpc/h]
  ShellRadius = (old_ShellRadius + (new_ShellRadius-old_ShellRadius) * nstep / STEPS);  // radius at this substep [Mpc/h]
  
  /*--- here we decide whether to write-out the current galaxy ---*/
  writeOut = 1.;
#ifdef STORE_SOURCES
  if (ShellRadius < RshellMin * Gal[p].Rvir)
    writeOut = -1.;  // ShellRadius did not satisfy the threshold, so we do not save out
#endif

  if (writeOut > 0.)
    {
      struct SHELLINFO source;  // dummy structure to save out.
      double mytime, dt_int, deltat, prevtime, swept, rhocrit, hubble_z, metalRatio;

      /*--- computing the time at current substep using Galaxy clock ---*/
      prevtime = NumToTime(Gal[p].SnapNum);
      deltat = (prevtime - NumToTime(Gal[p].SnapNum+1));
      dt_int = deltat / STEPS;                      // substep dt [int.units]
      mytime = prevtime - (nstep + 0.5) * dt_int;   // lookback time at the middle of the substep
      
      /*--- computing the IGM mass swept by the shell propagation, assuming:  ---*/
      /*--- loading factor==1 and flat,uniform universe with density=critical ---*/
      hubble_z = hubble_of_z(Gal[p].HaloNr);
      rhocrit = 3 * hubble_z * hubble_z / (8 * M_PI * G);     // internal units --> [1.e10Msun * h^2 / Mpc^3]
      swept = 4 * M_PI * (pow(ShellRadius,3.0) - pow(Gal[p].Rvir,3.0)) * rhocrit;  // shell swept-mass. [1.e10 Msun/h]
      metalRatio = ZSUN * pow(10, ZIGM[Gal[p].SnapNum]);      // IGM metal fraction (given the average IGM metallicity)

      /*--- All the following quantities are in internal units. ---*/
      /*--- Mass:                [1.e10 M_sun/h]     ---*/
      /*--- Position or Radius:  [Mpc/h]             ---*/
      /*--- Time:                [(Mpc/h) / (Km/s)]  ---*/
      source.mvir = Gal[p].Mvir;
      source.rvir = Gal[p].Rvir;
      source.xpos = Gal[p].Pos[0];
      source.ypos = Gal[p].Pos[1];
      source.zpos = Gal[p].Pos[2];
      source.rshell = ShellRadius;
      source.shellmetals = Gal[p].MetalsEjectedMass + swept*metalRatio;
      source.shellgas = Gal[p].EjectedMass + swept;
      source.time = mytime;
      source.haloID = Gal[p].HaloNr; // Integer: halo ID.
      source.pcles = Gal[p].Len;     // Integer: number of particles in the associated subhalo
      source.boundID = Gal[p].MostBoundID;  // longlong integer: ID of the most bound particle
      
      fwrite(&source, sizeof(struct SHELLINFO), 1, OutFilePointer);
    } 
}



/** @brief "update_metal_shell_radius" tracks the evolution in time of the metallic bubble/shell
 *         launched by the SNe feedback associated to strong star-formation events. More in detail,
 *         it associates a "shell-radius" to any galaxy who is able to produce "EjectedMass" in a
 *         given SF episode. From the moment of the SF episode, it evolves the radius of the shell
 *         in time, computing one value per SnapNum until z=0. The radius evolution is modelled
 *         exactly as in Dijkstra et al. 2014 (eq. 4) and it does not differentiate between kinetic
 *         or thermal modes to support the bubble/shell expansion (see e.g. Dijkstra et al. 2014).
 *         NOTES:
 *         1) the shell radius is launched starting from "HotRadius". The computation of this value
 *            is already implemented in the code, and it represents the radius of the hot atmosphere
 *            surrounding galaxies. In particular: HotRadius = Rvir if Gal.Type==0, else its value
 *            is reduced due to stripping; see model_stripping.c
 *         2) every time a galaxy produces ejecta, the update of the shell radius is triggered and
 *            a new radius evolution is computed. If, during this new evolution, the new radius of
 *            the shell overcomes the previously-computed values, then the shell radius is updated
 *            with the newly-computed values. In this way, we account for subsequent SF episodes
 *            over the whole "shell life".
 *         3) The Dijkstra model depends on the following physical parameters and quantities:
 *            e0       : energy of 1 supernova explosion in cgs, i.e:  1 Bethe = 10^51 erg.
 *            h0       : Hubble constant at z=0, needed to compute critical density at z=0. Units: [(km/s)/Mpc]
 *            Hubble   : Cosmology-independent hubble constant. Units: [(km/s)/Mpc]. Global variable in LGalaxies. Should be 100.00 but is 100.000066361520
 *            rhocrit  : critical density at z=0. Units: [1e10 Msun]/[Mpc^3]
 *            G        : gravitational constant. Global variable in Lgalaxies. Int.units: [Mpc^3]/([10^10Msun]*[Mpc/(km/s)]^2) (see init.c -> set_units())
 *            Ob       : omega_baryon at z=0 -> omega matter * average cosmic baryon fraction (f_b=0.15 in Planck15 cosmology)
 *            n        : MASS density of the IGM inside which the metallic shell expands. See Eq.4 in Dijkstra_Ferrara_Mesinger-2014
 *            dtsep    : time step at the current SnapNum. Units: [Mpc/(km/s)]
 *            SNfreq   : number of supernovae produced by a unit of stellar mass. Units: 1 / ([1e10 Msun]/[h])
 *            IMPORTANT NOTE on SNfreq:
 *            Eq.4 of Dijkstra et al. 2014 needs the "supernovae frequency per units of stellar mass"
 *            i.e. the number of SNe produced per unit of the galaxy's total stellar mass. To guarantee
 *            self-consistency in the code, we use the "RecycleFraction" parameter of LGalaxies, since
 *            this can be thought as the number of SNe exploded per unit of newly-formed stellar mass).
 * */
//-----------------------------------------------------------------------------------------
void update_metal_shell_radius(int p, int nstep, double time, double stars)
{
  float r_dummy, SNfreq, zz;
  float e0, h0, rhocrit, Ob, n, cdijkstra;         //variables to compute factors of eq. 4 in Dijkstra et al. 2014
  int snap;                                        //loop variables for shell radius computation
  float tt, dtstep, time_now, time_pre, time_nex;  //loop variables for shell radius computation
 
  e0 = EnergySN / UnitEnergy_in_cgs;                         // [10^10Msun]*([Mpc]^2) / [Mpc/(Km/s)]^2
  h0 = Hubble * Hubble_h;                                    // [(Km/s)/Mpc]
  rhocrit = 3 * h0 * h0 / (8 * M_PI * G);                    // [1e10 Msun]/([Mpc]^3)
  Ob = Omega * BaryonFrac;                                   // [adim]
  n = 60. * Ob * rhocrit * pow((1.+ZZ[Gal[p].SnapNum]),3.);  // [1e10 Msun]/([Mpc]^3)

  /*--- lookback (l.b.) time quantities without the "h" given by "NumToTime" routine   ---*/
  /*--- (and inherited from "time_to_present" routine). All quantities in [(Km/s)/Mpc] ---*/
  time_pre = NumToTime(Gal[p].SnapNum)/Hubble_h;    // at previous SnapNum.
  time_nex = NumToTime(Gal[p].SnapNum+1)/Hubble_h;  // at next SnapNum.
  dtstep = (time_pre - time_nex) / STEPS;           // delta-t associated to current substep.
  time_now = time_pre - nstep*dtstep;               // time at current substep.

  /*--- Eq.4 of Dijkstra et al. 2014; we use the SNfreq given by "RecycleFraction" ---*/
  SNfreq = RecycleFraction * 1e10/Hubble_h;  // units: 1 / ([1e10 Msun]/[h])
  cdijkstra = ( e0 * SNfreq * stars ) / n;   // units: ([Mpc]^5 / [Mpc/(Km/s)]^2)
    
  if (cdijkstra <= 0.)
    terminate("Problem in update_metal_shell_radius: Null or negative constant C!!\n\n");

  /*--- Once the constant "cdijkstra" is obtained (sort of an "expansion speed" for ---*/
  /*--- the shell), we can compute the radius evolution and eventually update it.   ---*/
  for(snap=Gal[p].SnapNum; snap<MAXSNAPS; snap++)
    {
      if (ZZ[snap] >= 0.)  // avoids snapshots with z<0 in the rescaled cosmology
        {
          tt = time_now - NumToTime(snap+1)/Hubble_h;  //[Mpc/(km/s)]
          if (tt < 0.)
	    terminate("Problem in update_metal_shell_radius: negative time?!?!?!\n\n");

	  // Shell radius computation. If the new values are bigger
	  // than the old ones, we will update the old values
          r_dummy = pow(cdijkstra, 1./5.) * pow(tt, 2./5.) + Gal[p].HotRadius;  // [Mpc].
          if (r_dummy > Gal[p].MetalShellRadius[snap])
	    Gal[p].MetalShellRadius[snap] = r_dummy;
        }
    }
}
