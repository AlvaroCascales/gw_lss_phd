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


/** @file model_lyman_werner.c
 *  @brief This module collects all the routines devoted to compute the local variations
 *         of the Lyman-Werner background, i.e. a specific UV band necessary to model the
 *         formation of Black-Holes at very high-redshift. Every galaxy in the code which
 *         is actively forming stars will have an associated Lyman-Werner luminosity, in
 *         relation to the amount of young, massive and hot stars produced in the last SF
 *         event. Indeed, LyWer photons are very efficiently produced by massive OB stars
 *         (see e.g. Schaerer D. 2002 or Greif_Bromm 2006), which are instantly recycled
 *         in the LGalaxies code. This means that each SF event will update the LyWer
 *         luminosity (hence considered as "instantaneous luminosity") of each SF galaxy.
 *
 *
 *  @brief Computation of the spatially varying Lyman-Werner background.
 *         The Lyman-Werner (in brief: LyWer) band is a specific component of the UV
 *         spectrum, which is able to photo-ionize H2 molecules. The ionising potential
 *         of H2 molecules implies that LyWer photons have energies in the interval:
 *         [11.2 - 13.6] eV. These photons are particularly important to suppress gas
 *         cooling and fragmentation (hence to hinder star-formation) in young collapsing
 *         halos with a very pristine chemical composition. Indeed, in these conditions
 *         (especially frequent at high redhisft), the gas cooling via atomic transitions
 *         lines is not efficient, while the roto-vibrational transitions of H2 molecules
 *         are the main cooling channel. The LyWer photons act on these "mini-halos" by
 *         specifically suppressing their H2 content, and hence reducing their ability to
 *         cool, fragment and form stars. This, in turn, paves the way for the existence
 *         of "monolithic" clouds of pristine gas, which can collapse into a single object
 *         (via e.g. gravitational instabilities) and form a so-called "direct-collapse"
 *         black-hole (DCBH) seed, with a mass of ~ 10^4 - 10^5 M_sun. In order to model
 *         the formation of these heavy seeds, it is thus crucial to track the spatial
 *         variations of the LyWer background at the potential sites of high-z Black-Hole
 *         formation (in our case: where new halos are initialized during the code runs).
 *         As for the case of the local variation of the metallic pollution (see the
 *         "model_igm_metallicity.c" module) we do not rely on pre-computed grids of
 *         spatially varying LyWer background. Instead, we pre-select galaxies with very
 *         intense Star-Formation rates in a preliminary run of the code, and store them
 *         into an external list of potential "LyWer providers" (or just "LyWer sources")
 *         defined over the whole simulated box. During the following runs of the code,
 *         for each newly-initialized galaxy whose chemical composition is pristine
 *         (i.e. lower than a given input threshold), we check the TOTAL LyWer illumination
 *         received by the newly-initialized galaxy from:
 *         1) its eventual neighbors belonging to its FoF-group,
 *         2) ALL the LyWer sources in the external list.
 *         NOTES:
 *         1) at very high redshift, the horizon-length of LyWer photons is of the order
 *            of 100 Mpc (see e.g. Ahn et al. 2009), hence it is FUNDAMENTAL that the
 *            list of LyWer sources is defined over the WHOLE simulated box, since a
 *            non-negligible contribute to the local LyWer illumination can be provided
 *            by relatively far-away actively star-forming galaxies.
 *         2) Due to their very big horizon length, LyWer photons emitted by a source can
 *            be affected by cosmological redshift and hence "drop-out" from the LyWer band.
 *            At the same time, more energetic photons emitted by the same source are also
 *            affected by cosmological redshift and can "enter" into the LyWer band. As a
 *            consequence, to compute the exact local LyWer illumination at a given point
 *            of a simulated box, it should be necessary to implement a radiative-transfer
 *            treatment (which should also account for the opacity towards UV photons of
 *            the non-reionized IGM at z>8). Fortunately, Ahn et al. 2009 demonstrated that
 *            this exact approach can be well approximated by a fully-analytic approach.
 *            We then base our computation on the Ahn et al. 2009 theoretical background.
 *         3) The total local LyWer illumination for a newly initialized galaxy is the sum
 *            of all the contributions inside its past light-cone. Each of these ones will
 *            depend on the distance "D" between the source and the newly-initialized halo
 *            as: D^-2. In order to avoid time-consuming loops-back in the past, we stop to
 *            account for contributions to the local LyWer FLUX when the distance between
 *            an hypothetic very-strong LyWer source and the newly-initialized galaxy is
 *            such that the hypothetic very-strong source would only minimally contribute
 *            to the local LyWer illumination. For the exact definitions of these conditions
 *            see the comments in the "check_LyWer_illumination" routine.
 * */

/*------  ROUTINES FOR DEFINING AND STORING LyWer SOURCES AND CONTRIBUTORS (especially needed for BH seeding)  ------*/

/** @brief The following routine opens the binary file where the list of LyWer sources is stored
 *         (it should have been saved there in a preliminary run of the model), counts the number
 *         of items in the list, allocates the necessary memory according to this number and then
 *         stores the file content in a global array of LyWer sources (used for BH seeding). The
 *         latter, will then be used in a following run to compute the local illumination at eventual
 *         Black-Hole seeding candidates (if the "SEEDING" option of the Makefile is ON). The list
 *         of LyWer-Sources is defined in the comments of the "store_LyWer_source" routine.
 * */
void read_LyWer_sources(void)
{
  LYWERSOURCESFILE = fopen(FileWithLyWerSources, "rb");

  if (LYWERSOURCESFILE != NULL)
    {
      /*---  first we count the elements in the external list of LyWerSources ---*/
      /*---  in order to properly allocate the memory to store and use them.  ---*/
      LYWERSOURCESNR = 0;
      struct LYWERINFO source;   // dummy structure to read the external list
      
      while( fread(&source, sizeof(struct LYWERINFO), 1, LYWERSOURCESFILE)!=0 )
	LYWERSOURCESNR++;

      fclose(LYWERSOURCESFILE);
      printf("LyWer Sources counted in model_lyman_werner.c/read_LyWer_sources. Number of LyWer Sources found: N=%d\n",LYWERSOURCESNR);
      
      /*---  now we can allocate the memory we need ---*/
      LYWERSOURCE = mymalloc("LyWerSource", sizeof(struct LYWERINFO)*LYWERSOURCESNR);

      /*-- now we can read the external list again and store its content in the allocated array ---*/
      LYWERSOURCESFILE = fopen(FileWithLyWerSources, "rb");
      int ii=0;
      while( fread(&LYWERSOURCE[ii], sizeof(struct LYWERINFO), 1, LYWERSOURCESFILE)!=0 )
	ii++;
      fclose(LYWERSOURCESFILE);
      
      //CHECK FILE READING
      /* int num; */
      /* num = LYWERSOURCESNR-4; */
      /* print_LyWer_sources(num); */
      /* num = LYWERSOURCESNR-3; */
      /* print_LyWer_sources(num); */
      /* num = LYWERSOURCESNR-2; */
      /* print_LyWer_sources(num); */
      /* num = LYWERSOURCESNR-1; */
      /* print_LyWer_sources(num); */
      /* terminate("\nHope everything is fine with these numbers\n"); */
    }
  else
    terminate("\nLyWer sources not read in model_lyman_werner.c --> read_LyWer_sources\n");  
}


//-----------------------------------------------------------------------------------------
//  This routine just prints a single LyWers source
//  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  ONCE MADE SURE THAT THE WRITING/READING OF THE LYWER SOURCES IS OK, REMOVE THIS ROUTINE !!!!!!!!!!
//-----------------------------------------------------------------------------------------
void print_LyWer_sources(int num)
{
  printf("\nHalo %d: mvir=%e [i.u.]  (%e  [Msun])", num, LYWERSOURCE[num].mvir,LYWERSOURCE[num].mvir*1.e10/Hubble_h);
  printf("\nHalo %d: rvir=%f [Mpc/h]  (%f  [Mpc])", num, LYWERSOURCE[num].rvir,LYWERSOURCE[num].rvir/Hubble_h);
  printf("\nHalo %d: ColdGas=%e [i.u.]  (%e  [Msun])", num, LYWERSOURCE[num].coldgas, LYWERSOURCE[num].coldgas*1.e10/Hubble_h);
  printf("\nHalo %d: Stars=%e [i.u.]  (%e  [Msun])", num, LYWERSOURCE[num].stars, LYWERSOURCE[num].stars*1.e10/Hubble_h);
  printf("\nHalo %d: InstSfr=%e [Msun/yr]", num, LYWERSOURCE[num].InstSfr);
  printf("\nHalo %d: RecycledStars=%e [i.u.]  (%e  [Msun])", num, LYWERSOURCE[num].RecycledStars, LYWERSOURCE[num].RecycledStars*1.e6);
  printf("\nHalo %d: x=%f [Mpc/h]  (%2.2f  [Mpc])", num, LYWERSOURCE[num].xpos, LYWERSOURCE[num].xpos/Hubble_h);
  printf("\nHalo %d: y=%f [Mpc/h]  (%2.2f  [Mpc])", num, LYWERSOURCE[num].ypos, LYWERSOURCE[num].ypos/Hubble_h);
  printf("\nHalo %d: z=%f [Mpc/h]  (%2.2f  [Mpc])", num, LYWERSOURCE[num].zpos, LYWERSOURCE[num].zpos/Hubble_h);
  printf("\nHalo %d: time/redshift=%f", num , LYWERSOURCE[num].time);
  printf("\nHalo %d: Rshell=%f [Mpc/h]  (%f  [Mpc])", num, LYWERSOURCE[num].rshell, LYWERSOURCE[num].rshell/Hubble_h);
  printf("\nHalo %d: Metals=%e [i.u.]  (%e  [Msun])", num, LYWERSOURCE[num].metals, LYWERSOURCE[num].metals*1.e10/Hubble_h);
  printf("\nHalo %d: ShellGas=%e [i.u.]  (%e  [Msun])", num, LYWERSOURCE[num].shellgas, LYWERSOURCE[num].shellgas*1.e10/Hubble_h);
  printf("\n\n");
}



/** @brief The "store_LyWer_routine" is used to save a list of galaxies who either act as potential
 *         sources of Lyman-Werner photons (due to their high instantaneous SF rates) or galaxies
 *         who are actually able to provide a non-negligible flux of Ly-Wer photons to the location
 *         of a newly-initialized galaxy in the box. For these two cases, we speak respectively of
 *         Lyman-Werner "sources" and Lyman-Werner "contributors".
 *         1) LyWer "sources" are defined at every SUBSTEP during a preliminary run of the code,
 *            if the option "STORE_SOURCES" is active in the Makefile. Basically, LyWer sources
 *            are galaxies with very strong SF rates, which are collected into an "external list"
 *            of LyWer "sources", or "potential LyWer contributors". More details about how the
 *            SFR of galaxies are converted to UV (and Lyman-Werner) fluxes are provided in the
 *            documentation of the routines devoted to compute the LyWer flux.
 *         2) LyWer "contributors" are defined when checking for the local LyWer illumination at
 *            a specific location in the box, where a new galaxy is initialized. The "contributors"
 *            of this newly-initialized galaxy are either those members of its FoF-group or those
 *            LyWer sources (from the external list of LyWer sources) which are actively providing
 *            a non-negligible flux of LyWer photons to its location.
 *         The following routine is in charge of saving-out the properties of either "conributors"
 *         or "sources". The input pointer of the routine and whether the option "STORE_SOURCES" is
 *         active in the Makefile are the factors who decide where the input galaxy is going to be
 *         saved. If the input galaxy is a potential "source", then the routine checks the value
 *         of its instantaneous SFR, in order to effectively include it among the LyWer sources list.
 * */
void store_LyWer_source(int p, int nstep, FILE * OutFilePointer)
{
  struct LYWERINFO source;  // we'll need this in any case
  
#ifdef STORE_SOURCES
  /*
    in this case we have to decide whether to write-out the galaxy as a LyWer source
    hence we first check if some conditions are verified, then we write the data out
  */
  double mytime, dt_int, dt_yrs, prevtime;
  double rvir, lwsf;

  /*--- time quantities according to galaxy "clock" ---*/
  prevtime = NumToTime(Gal[p].SnapNum);
  dt_int = (prevtime - NumToTime(Gal[p].SnapNum+1)) / STEPS; // substep dt [int.units] -> [(Mpc/h) / (Km/s)]
  dt_yrs = dt_int * (UnitTime_in_years / Hubble_h);          // substep dt [years]
  
  /*--- galaxy properties to check the instantaneous SFR of massive, recycled stars (those producing LW flux) ---*/
  rvir = Gal[p].Rvir/Hubble_h;                   // rvir in [Mpc]
  lwsf = Gal[p].RecycledStars * 1.e6 / dt_yrs;   // iSFR in [Msun/yr]

  /*--- Here we select galaxies with J_21(Rvir) >~ 10 before z~6.5 (DEPENDS ON REDSHIFT!!! CHANGE THIS!!!) ---*/
  if ((lwsf > 0.0) &&  (lwsf/(rvir * rvir) > 2.e3) && (ZIGM[Gal[p].SnapNum]<=ZIGMcrit))
    {
      /*--- Computing the current properties of metallic shells of this LyWer source, just for post-processing analysis ---*/
      double old_ShellRadius, new_ShellRadius, ShellRadius;
      old_ShellRadius = Gal[p].MetalShellRadius[Gal[p].SnapNum] * Hubble_h;         // radius at this Gal[p].SnapNum, [Mpc/h]
      new_ShellRadius = Gal[p].MetalShellRadius[Gal[p].SnapNum+1] * Hubble_h;       // radius at next Gal[p].SnapNum, [Mpc/h]
      ShellRadius = (old_ShellRadius + (new_ShellRadius-old_ShellRadius) * nstep / STEPS);  // radius at this substep [Mpc/h]
      double hubble_z, rhocrit, swept, metalRatio;
      hubble_z = hubble_of_z(Gal[p].HaloNr);
      rhocrit = 3 * hubble_z * hubble_z / (8 * M_PI * G);     // internal units --> [1.e10Msun * h^2 / Mpc^3]
      swept = 4 * M_PI * (pow(ShellRadius,3.0) - pow(Gal[p].Rvir,3.0)) * rhocrit;  // shell swept-mass. [1.e10 Msun/h]
      metalRatio = ZSUN * pow(10, ZIGM[Gal[p].SnapNum]);      // IGM metal fraction (given the average IGM metallicity)
      
      source.mvir    = Gal[p].Mvir;
      source.rvir    = Gal[p].Rvir;
      source.coldgas = Gal[p].ColdGas;
      source.stars   = Gal[p].DiskMass+Gal[p].BulgeMass;
      source.InstSfr = lwsf;
      source.RecycledStars = Gal[p].RecycledStars;  // [1e6 Msun]
      source.xpos = Gal[p].Pos[0];
      source.ypos = Gal[p].Pos[1];
      source.zpos = Gal[p].Pos[2];
      source.time = prevtime - (nstep + 0.5) * dt_int;  // lookback time at the middle of the substep
      source.jlw = 0.0;  // this is non-zero only for LyWer contributors. Here we are saving a LyWer source.
      source.rshell = ShellRadius;
      source.metals = Gal[p].MetalsEjectedMass + swept*metalRatio;
      source.shellgas = Gal[p].EjectedMass + swept;
      
      fwrite(&source, sizeof(struct LYWERINFO), 1, OutFilePointer);
    }
#else

  terminate("\n\nwe are trying to save-out a LyWer contributor... is this ok???");   ////// remove after seeding is implemented!!!!!!!!!!!!!!!!!!!
    
  /*
    in this case we write out the LyWer contributor without further checks.
    We only need to change a couple of fields in the LyWer structure.
  */
  source.mvir    = LYWERSOURCE[LyWer_contrib_IDs[p]].mvir;
  source.rvir    = LYWERSOURCE[LyWer_contrib_IDs[p]].rvir;
  source.coldgas = LYWERSOURCE[LyWer_contrib_IDs[p]].coldgas;
  source.stars   = LYWERSOURCE[LyWer_contrib_IDs[p]].stars;
  source.InstSfr = LYWERSOURCE[LyWer_contrib_IDs[p]].InstSfr;
  source.RecycledStars = LYWERSOURCE[LyWer_contrib_IDs[p]].RecycledStars;
  source.xpos = LYWERSOURCE[LyWer_contrib_IDs[p]].xpos;
  source.ypos = LYWERSOURCE[LyWer_contrib_IDs[p]].ypos;
  source.zpos = LYWERSOURCE[LyWer_contrib_IDs[p]].zpos;
  source.time = LyWer_contrib_Zed[p];  // redshift of the Lywer contributor (stored in check_LyWer_illumination)
  source.jlw  = LyWer_contrib_Flu[p];  // actual flux received by the DCBH candidate from this LyWer contributor (computed in check_LyWer_illumination)
  source.rshell = LYWERSOURCE[LyWer_contrib_IDs[p]].rshell;
  source.metals = LYWERSOURCE[LyWer_contrib_IDs[p]].metals;
  source.shellgas = LYWERSOURCE[LyWer_contrib_IDs[p]].shellgas;  // [10^10 Msun / h]
  
  fwrite(&source, sizeof(struct LYWERINFO), 1, OutFilePointer);
  
#endif
  
}


/** @brief The following routine keeps track of the amount of massive, short-lived stars produced
 *         during the last SF event. Indeed, LyWer photons are most efficiently produced by these
 *         young, massive and short-lived stars (see e.g. Leitherer et al. 1999). After ~10^6 yrs
 *         it is reasonable to assume that the LyWer flux produced by a single starburst event is
 *         zero (see e.g. Leitherer et al. 1999). Since the LGalaxies substeps are much longer than
 *         10^6 years (usually >10^7 yrs), we assume that the mass of stars actually producing the
 *         LyWer photons is the one who will be instantly recycled (because these stars will explode
 *         as SNe). As a consequence, RecycledStars is the mass of stars producing LyWer photons.
 *         NOTE:
 *         1) Recycled stars are produced in any SF event, i.e. both by teh routines: starformation
 *            and collisional_starburst_recipe. The first routine is called if the conditions for
 *            star formation are met (over-critical mass of ColdGas), the second one is called if
 *            there is a merger which sends the ColdGas mass to over-critical values. IF both are
 *            called, then they are SURELY called during the SAME substep. In this case, then, we
 *            need to account for BOTH episodes of star formation to compute the TOTAL LyWer flux
 *            produced in the current substep. Therefore, if we are calling the routine for the
 *            second time (i.e. during a SF burst induced by a merger) we will INTEGRATE the new
 *            stars and the old ones, and recompute all the LyWer-related quantities.
 *         2) we do NOT directly compute and store the flux produced by each source in a peculiar
 *            frequency band, since this quantity depends on the distance between the photon-source
 *            and the observer. Indeed this distance also determines whether the UV photons will
 *            exit from the LyWer band before reaching the observer (due to cosmological redshift).
 *            Hence it is more convenient to simply keep track of the iSFR of massive, short lived
 *            stars and later use this quantity to compute the flux received by a specific observer
 *            at a well-defined distance and redshift.
 *         3) The following routine is called by "update_from_starformation" routine which itself
 *            is called only TWICE in the code:
 *            2.1) during an ordinary SF event (model_starformation_and_feedback.c --> starformation)
 *            2.2) during a merger-induced burst (model_mergers.c --> collisional_starburst_recipe)
 *            These two calls are done in the same substep, on after the other, so we keep track of
 *            BOTH events by summing the amount of stars produced by both.
 * */
void update_LyWer_luminosity(int p, bool flag_burst, double dt, double stars)
{
  double dist;
  double dnu, fesc;
  double Erate, mstars;

  /*--- Physical parameters(eventually set from input_parameter_file) ---*/
  Erate = 1.e44;                     //energy-emission rate in the LW band [erg/(s*10^6Msun)]. IMPORTANT: THIS IS NOW FIXED BUT IT MUST DEPEND ON AGE AND METALLICITY OF THE STELLAR POPULATION!!!!!!!!
  fesc = 1.;                         //escape fraction of Lyman-Werner photons  !!!! SET AS FREE PARAMETER? INCLUDE IGM ABSORPTION (e.g. Agarwal 2012)?????????
  dnu = (3.288 - 2.708) * 1.e15;     //Lyman-Werner bandwidth in [Hz] (Hertz). This corresponds to the energy range: [13.6 - 11.2] eV
    
  if (flag_burst)
    Gal[p].RecycledStars += (stars * RecycleFraction * 1.e4/Hubble_h);  // [1e6 Msun]
  else
    Gal[p].RecycledStars = (stars * RecycleFraction * 1.e4/Hubble_h);   // [1e6 Msun]

  /*--- updating LyWer-related variables to keep track of LyWer quantities in the outputs ---*/
  Gal[p].LyWerStars = Gal[p].RecycledStars;
  Gal[p].LyWerLuminosity = (fesc/dnu) * (Erate * Gal[p].LyWerStars);            //now in [erg/(s*Hz)]
  Gal[p].LyWerSfr = Gal[p].LyWerStars / (dt/(Hubble_h/UnitTime_in_Megayears));  //actual Sfr at the moment of LW luminosity production in 10^6 Msun / Myr (i.e. Msun/yr)
  
}




/*------  ROUTINES FOR ACTUALLY COMPUTING LOCAL VARIATIONS OF LyWer ILLUMINATION (through a past-light-cone computation)  ------*/

/** @brief This routine follows the analytic approach of Anh et al.2009 (eqs. 22-23) to compute
 *         the Lyman-Werner horizon at a given redshift. This is the maximum distance at which
 *         ionizing UV photons are expected to travel before being either received by an observer
 *         as LyWer photons or absorbed by a resonant line of neutral hydrogen in the IGM.
 * */
float get_LyWer_horizon(float zz)
{
  float alpha, lw_hor;

  alpha = (0.7/Hubble_h) * pow(0.27/Omega, 0.5) * pow(21.0/(1+zz), 0.5);
  lw_hor = 97.39 * alpha * Hubble_h;
  
  return lw_hor;
}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// continue from here:
//  - check all the new variables implementation!!!!  The RecycledStars were barely implemented...............
//  - final checks and then push to git



/** @brief This routine follows the analytic approaches of Anh et al.2009, Agarwal et al. 2012
 *         and Dijkstra et al. 2014 to compute the Lyman-Werner flux received by an observer
 *         at a distance "d" from a given source of LyWer photons. This computation is based
 *         on several approximations assumed to avoid a full radiative-transfer approach (see
 *         Ahn et al. 2009 for details). Indeed, the latter should describe the energy shift
 *         of the whole spectrum of UV-ionising photons emitted by the source, as caused by
 *         cosmological redshift along their path from the source until the observer.
 *         The main feature of the Ahn et al. 2009 approximation is the use of a "modulation
 *         factor" (fmod) which condensates all the action of the IGM along the photons path
 *         and allows to compute the received LyWer-flux (from the observer point of view)
 *         just as a function of the photons emitted by the source IN THE LY-WER BAND.
 *         The fmod factor depends directly on redshift (due to the evolving properties of
 *         the IGM with time) and, obviously, on the distance between source and observer.
 *         Finally, due to the IGM absorption, UV-ionising photons can only travel up to a
 *         finite distance before being absorbed. The exact distance they can travel depend
 *         primarily on their energy (at the emission moment) and on the IGM density. To
 *         describe this complex phenomenon, the Ahn+2009 approximation introduces a LyWer
 *         "horizon" which represents the maximum distance at which UV-ionising photons can
 *         travel and being received as LyWer photons (independently of their initial energy).
 *         IMPORTANT NOTES:
 *         1) due to the "instantaeous recycling" approximation, LGalaxies assumes that a
 *            fixed fraction of newly-formed stars (namely, the "RecycleFraction" variable)
 *            is composed of hot, massive stars which immediately explode (i.e. in the same
 *            substep in which they are created). At the same time, these stars are the most
 *            efficient producers of UV photons. Therefore, we directly connect the SFR of
 *            these massive stars to the production of UV photons. A full treatment should
 *            connect the production of UV flux (and the "energy content" of this flux) to
 *            the actual properties of the current stellar population (age and metallicity).
 *            For the moment, this is not implemented (but it could be implemented).
 *         2) Another important aspect regarding the emission of UV-ionizing radiation from
 *            galaxies is the modelling of the escape fraction. This should account for the
 *            absorption of UV-ionising photons in the inter-stellar medium (e.g. by dust).
 *            We greatly simplify this model by simply using and calibrating a (scalar)
 *            escape-fraction parameter.
 *        PHYSICAL PARAMETERS:
 *        Erate: energy-emission rate in the LW band. Units: [erg/(s*10^6Msun)], as in Agarwal+2012.
 *        fesc:  escape fraction of LyWer photons from the source. It's a free-parameter
 *        dnu:   Lyman-Werner bandwidth in Hertz. Corresponds to the energy range: [13.6 - 11.2] eV
 * */
double compute_LyWer_flux(double starmass, double zz, double d, float LyWerHor)
{
  double LyWerFlux, Erate, fesc, dnu, dd;
  double alpha, fmod;
  /*--- Physical parameters ---*/
  Erate = 1.e44;   // IMPORTANT: THIS IS NOW FIXED BUT IT MUST DEPEND ON AGE AND METALLICITY OF THE STELLAR POPULATION!!!!!!!!
  fesc = 1.;       // IMPORTANT: SET AS FREE PARAMETER? INCLUDE IGM ABSORPTION (e.g. Agarwal 2012)?????????
  dnu = (3.288 - 2.708) * 1.e15;
  
  /*--- Here we compute the Ahn+2009 attenuation factor (see comment above for details) ---*/
  alpha = (0.7/Hubble_h) * pow(0.27/Omega, 0.5) * pow(21.0/(1+zz), 0.5);
  printf("Alpha = %2.10f\n", alpha);
  fmod = 0.0;
  if (d/alpha <= 97.39)
    fmod = (1.7 * exp(-1.*pow(d/(116.29*alpha),0.68))) - 0.7;
  printf("fmod = %2.10f\n",fmod);

  /*--- Here we actually compute the LyWer flux received by the current observer from the current LyWer-source ---*/
  LyWerFlux = 0.0;  // if the source is outside the LyWer horizon of current observer, the LyWer flux is null
  dd = d * 3.0857;  // in units of [10^24 cm]. 3.0857 = conversion factor from Mpc to 10^24 cm (weird units otherwise numbers too big)
  if (d/alpha < LyWerHor)
    LyWerFlux = fmod * (fesc / dnu) * (Erate * starmass) / ( pow(4. * M_PI * dd * (1.+zz),2.0) );
  LyWerFlux = LyWerFlux * 1.e-27;  //now in units of 10^-21 [erg]/[s][Hz][sr][cm2]. Obtanied from: jlw * [10^-48 --> cm^2] * 10^21
  printf("LyWerFlux = %2.10f\n\n",LyWerFlux);
  
  if ((isnan(LyWerFlux)==1) || (isinf(LyWerFlux)==1))
    {
      printf("LW flux=%E, zz=%f, d=%E, LyWerHor=%f\n", LyWerFlux, zz, d, LyWerHor);
      terminate("Error: Lyman-Werner flux either INF or NAN in function compute_LyWer_flux\n");
    }

  return LyWerFlux;
}




//-----------------------------------------------------------------------------------------
// COMPUTE ILLUMINATION BY LYMAN-WERNER PHOTONS
// Here we compute the UV illumination from Lyman-Werner photons at the location
// of a given seeding-candidate. LW photons travel as far as ~100 cMpc (see Ahn et al. 2009),
// hence the contribution to the local UV flux can likely come from outisde the
// FoF of the seeding halo. To compute the local LW flux at each redshift, one
// could build grids of LW illumination and read from there the local LW whenever
// needed. This proved to be time consuming hence a "on-the-run" computation was
// chosen. Hence, we look for sources of LW-photons (LW sources) within a list of
// STRONG LW-emitters (see model_lyman_werner.c  -->  store_LyWer_sources) which
// is saved during a "preliminary" run of LGalaxies (without seeding). For every
// DCBH candidate, then, we integrate the local LW flux received by all the LW
// sources in the past light-cone of the DCBH candidate and then check if the
// cumulative LW flux is over a given threshold "Jcrit". If so, we put a massive
// DCBH (M_bh ~ 10^5 Msun) in the halo.
// NOTE: the extent of the light-cone (towards the past) is defined by the
//       Lyman-Werner horizon, which is the maximum distance (i.e. maximum
//       lookback time) from which we can receive LW photons (Ahn et al. 2009).
//       Usually it's not necessary to loop-back until the lightcone limit,
//       since, way before that, LW sources are too distant (i.e. faint) to
//       actively contribute to the local LW flux at the DCBH site. So we
//       stop the LW integration when the contribution of a conservatively
//       strong LW source is 1/100 of the Jcrit threshold (i.e. even the
//       strongest LW sources are too faint to contribute).
// NOTE: sometimes the link between Halo.SnapNum and Gal.SnapNum is lost
//       (i.e. they differ by more than 1). Hence, mixing time quantities
//       measured by Halo.SnapNum and Gal.SanpNum results in errors. Here
//       we perform all the computation in the Gal.SnapNum "clock" reference,
//       hence we re-define all the time-related quantities.
//-----------------------------------------------------------------------------------------
void check_LyWer_illumination(int p, int nstep)
{
  // variables for loop activation
  double masslim = 5.e5;  //minimum cold gas mass for the galaxy to be potential DCBH formation site  !!!!!!!!!!!!!!!!!!!!!!!!!!! <--- USE AN EXTERNAL VARIABLE READ FROM INPUT_PARAMETERS_FILE?
  double gasfrac, Munit;

  Munit = 1.e-10*Hubble_h;    // mass units converter (multiply by it to go from Msun to internal units)
  gasfrac = Gal[p].ColdGas / (Gal[p].ColdGas + Gal[p].DiskMass + Gal[p].BulgeMass); // current galaxy's gas fraction

  //--- LIGHT-CONE COMPUTATION: only if conditions for DCBH seeding are met ---//
  if ((gasfrac > 0.99) && (Gal[p].ColdGas > (masslim * Munit)) && ((Gal[p].BulgeMass + Gal[p].DiskMass) < 1.e2 * Munit))
    {
      // lightcone loop variables
      double prevtime, nexttime, dtstep, diti, tdum, tdum2, t0, t1;  // time loop computation
      int snp, stp, nstep_dum, ns;               // loop on SnapNums, Steps and LyWer sources
      double rmin, rmax, d, z_sour, z_rel;          // distance and redshift in the time loop
      double LyWerTot, LyWerDummy, m_dum, d_dum, f_dum, sfr_thresh;      // LyWer integration
      double t_cand, z_cand, LyWerHor, t_lw;         // seeding candidate reference variables
        
      // time quantities in the "galaxy clock"
      t_cand = 0.;
      prevtime = NumToTime(Gal[p].SnapNum);      // previous snapshot
      nexttime = NumToTime(Gal[p].SnapNum+1);    // next snapshot
      dtstep = (prevtime - nexttime) / STEPS;    // current substep
      t_cand = prevtime - nstep * dtstep;        // current l.b. time (at SnapNum&nstep where the candidate DCBH is)
      z_cand = redshift_at_t(ZZ[Gal[p].SnapNum], nstep, dtstep);  // redshift at the current lookback time

      // lyman-werner horizon: from current lookback-time we go back of deltaT = LyWerHor_light_travel_time
      LyWerHor = get_LyWer_horizon(z_cand);      // [Mpc/h]
      t_lw = LyWerHor/(C_LIGHT * 1.e-5) + t_cand;      // [Mpc/h]/[km/s] - here we go back of deltaT from t_cand
      
      // Initialisations of loop variables
      nstep_dum = nstep;            // nstep_dum loops-back on steps over the past light-cone of the DCBH candidate
      rmin = 0.;                    // [Mpc/h]
      if( nstep==0 )
	dtstep = (NumToTime(Gal[p].SnapNum-1) - NumToTime(Gal[p].SnapNum)) / STEPS;
      rmax = dtstep * (C_LIGHT * 1.e-5);  // [Mpc/h]
      LyWerTot = 0.0;               // integral of Lyman-Werner flux. initialized to zero for each potential DCBH
      tdum2 = 0.0;                  // time-threshold; stops the loop-back on light-cone when LW sources are faint
      sfr_thresh = 200.;            // conservatively-high SFR threshold [Msun/yr]  !!!!!!!!!!!!!!!!!!!!!!!!!!! <--- USE AN EXTERNAL VARIABLE READ FROM INPUT_PARAMETERS_FILE?
	
      //--- lookback loop on past-lightcone ---//
      for(snp=Gal[p].SnapNum; snp>=0; snp--)
	{
	  t0 = NumToTime(snp);
	  t1 = NumToTime(snp+1);
	  diti = (t0 - t1) / STEPS; // substep for the current "snp" (substep changes while looping on SnapNums)
	  
	  // here we loop "back" on steps to look for Lyman-Werner sources in the current Galaxy's light-cone
	  for(stp=nstep_dum; stp>=0; stp--)
	    {
	      if( stp==0 )
		diti = (NumToTime(snp-1) - NumToTime(snp)) / STEPS;
	      tdum = t0 - (diti * stp);   // lookback time from DCBH candidate (while looping back)
	      
	      z_sour = redshift_at_t(ZZ[snp], stp, diti);   // redshift of (eventual) LW sources
	      z_rel = (1+z_sour) / (1+z_cand) - 1.;         // relative redshift of (eventual) DCBH and LyWer source
			      
	      if (fabs(z_rel) < 1.e-5)
		z_rel = 0.;
	      
	      if (z_rel < 0.)
		terminate("Negative redshift in model_seeding.c/DCBH_seeding!!!\n");
	      
	      //--- Computations to stop the lookback-loop ---//
	      // The logic is: if we are too far in the past, we either go out of the LW horizon or we are at a
	      // distance too big for even very strong sources to contribute significantly. If so, we can stop
	      // the LW integration and save time. For this, we check that both: 1) the incremental lookback-time
	      // at each step-back is smaller than the time-extent of the LW horizon and 2) a very star-forming
	      // galaxy at the current look-back distance would produce a significant LW flux. For the latter check
	      // we use a high Star-Formation threshold (see sfr_thresh above) and compute the mass of stars formed
	      // in the current "diti", with that sfr_thresh (first line below here). Then, we compute the LW flux
	      // produced produced by those stars (i.e. in a source which keeps a SFR=sfr_thresh for a dt=diti),
	      // as recieved at a distance: d_lum = (1+zrel)*d_dum  (second line below here)
	      tdum2 += diti;                                                    // [int.Units] - incremental lookback time from the DCBH
	      d_dum = tdum2 * (C_LIGHT * 1.e-5) / Hubble_h;                           // [Mpc] - distance associated to tdum2	      
	      m_dum = sfr_thresh * diti * (UnitTime_in_years/Hubble_h) * 1.e-6; // [1e6 Msun]; Stars formed in diti
	      f_dum = compute_LyWer_flux(m_dum, z_rel, d_dum, LyWerHor);        // [J21]; LW flux at distance d_dum
	      
	      // if the following is true, we are out of the light-cone (either earlier
	      // than the LW horizon lookback-time or too far to get significant LW flux)
	      if ((tdum2 >= t_lw) || (f_dum < Jcrit/100.))
		{
		  stp = -1;      // breaks the loop on steps
		  snp = -1;      // breaks the loop on snaps
		}
	      // if the above is false, we are inside the LW lightcone, hence we integrate LW contributions
	      else
		{
		  for(ns=0; ns<LYWERSOURCESNR; ns++) 
		    {
		      // we first look for sources at the correct lookback time from the observer
		      if ((LYWERSOURCE[ns].time >= tdum) && (LYWERSOURCE[ns].time < tdum+diti))
			{
			  d = pow(Gal[p].Pos[0]-LYWERSOURCE[ns].xpos,2.0) + pow(Gal[p].Pos[1]-LYWERSOURCE[ns].ypos,2.0) + pow(Gal[p].Pos[2]-LYWERSOURCE[ns].zpos,2.0);
			  d = pow(d,0.5);  // [Mpc/h]
			  
			  // now we also look for sources at the correct "lookback distance". d>0 excludes the DCBH candidate
			  if ((d >= rmin) && (d < rmax) && (d > 0.))
			    {			      
			      // contribution of current source to local LF flux received
			      LyWerDummy = compute_LyWer_flux(LYWERSOURCE[ns].RecycledStars, z_rel, d, LyWerHor);
			      LyWerTot += LyWerDummy;
			      
			      // storing some data of the current source (now a contributor), to write it out later
			      LyWer_contrib_IDs[counter_LyWer_contrib] = ns;          // the ID of the contributor
			      LyWer_contrib_Flu[counter_LyWer_contrib] = LyWerDummy;  // actual LyWer contribution of the current contributor
			      LyWer_contrib_Zed[counter_LyWer_contrib] = z_sour;      // redshift of the contributor
			      
			      counter_LyWer_contrib++;
			    }//end space-selection of external LyWer sources
			}//end time-selection of external LyWer sources
		    }//end loop on LYWERSOURCEs
		}// end else (looking for LYWERSOURCEs within lw_horizon)
	      
	      rmin += (C_LIGHT * 1.e-5) * diti;  // in [Mpc/h]
	      
	      if( stp==1 )  // in this case rmax must be computed with the new dt from the "next" snp in the loop
		rmax += (C_LIGHT * 1.e-5) * ((NumToTime(snp-1) - NumToTime(snp)) / STEPS);  // in [Mpc/h]
	      else
		rmax += (C_LIGHT * 1.e-5) * diti;  // in [Mpc/h]
	      
	    }//end loop in "lookback" substeps
	  nstep_dum = STEPS-1;
	}// end loop on "lookback" SnapNums
      
      //--- Evaluating the local LyWer flux ---//
      if (LyWerTot >= Jcrit)
	{
	  Gal[p].Pristine = 5;  // above the critical threshold for local LW flux. We'll form a DCBH
	  
	  // storing some properties of the DCBH-candidate to be saved later
	  if (counter_LyWer_contrib <= 498)
	    {
	      LyWer_contrib_Flu[counter_LyWer_contrib+1] = LyWerTot; // we append the total LyWer received by the seeding candidate after those of the LyWer contributors
	      LyWer_contrib_Zed[counter_LyWer_contrib+1] = z_cand;   // we append the seeding-candidate redshift after those of the LyWer contributors
	    }
	  else
	    terminate("\n\nProblem in check_LyWer_illumination() routine: found more than 499 Lyman-Werner contributors! There is no space to store the DCBH properties\n");
	  
	}
      else
	Gal[p].Pristine = 6;  // below the critical threshold for local LW flux. We'll form a PopIII remnant     
      
    }//end if, i.e. end loop looking for DCBH candidates for the current "p" galaxy    
}
