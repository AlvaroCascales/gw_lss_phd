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
 *      Author: David Izquierdo-Villalba */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

/** This module deals with the Pseudobulge, Calssical bulge and elliptical treatment
    It also includes the Smooth accretion
    It also includes the disipation during merger to properly fit the observed bulge sizes
*/

void init_Bulge_variables(int p, int Option)
{
 // Here I init all the properties that I will use for the bulge tracking
 if(Option < 0) 
  {
 	int j;
 	for(j=0;j<3;j++) Gal[p].StellarBulgeSpin[j] = 0.0;
 	Gal[p].Smooth_minor_merger = -2.;  // -2 means no smooth no normal minor merger, i.e just initialize the variable
 	Gal[p].induced = -1;
 	Gal[p].v_sigma_now = -1.;
 	Gal[p].t_dyn_last_disk = 0;
 	Gal[p].t_dyn_disk_merger = 0; 
 	Gal[p].counter_Delta_M = 0.;
 	Gal[p].Delta_M = 0.;
 	Gal[p].step_started_DI = -1;
 	Gal[p].snp_started_DI = -1;
 	Gal[p].M_stellar_started_DI = 0;
 	Gal[p].is_the_first_after_minor_Merger = -1.;
 	Gal[p].Mcrit_merger = -1.0;
 	Gal[p].Mburst_merger = 0.;
 	Gal[p].increment_due_to_burst = -99;
 	Gal[p].happened_minor_merger = -1.;
 	Gal[p].Mstellar_sec = 0.;
 	Gal[p].sfr_Burst_merger = -1.;
 	Gal[p].sfr_merger = -1.;
 	Gal[p].t_last_disk = 0.;
 	Gal[p].t_last_merger = 0.;
 	Gal[p].t_last_MINOR = 0.;
 	Gal[p].t_last_MAJOR = 0.;
 	Gal[p].t_last_DI = 0.;
 	Gal[p].Snap_last_merger = -1;
 	Gal[p].step_last_merger = -1;
 	Gal[p].Snap_last_DI = -1;
 	Gal[p].step_last_DI = -1;
 	Gal[p].z_last_merger = -1.;
 	Gal[p].z_last_Minor_merger = -1.;
 	Gal[p].z_last_Major_merger = -1.;
 	Gal[p].z_last_Smooth_accretion = -1.;
 	Gal[p].mass_ratio_mm = -1.;
 	Gal[p].mass_ratio_Mm = -1.;
 	Gal[p].n_Min_Mer = 0;
 	Gal[p].n_Maj_Mer = 0;
 	Gal[p].n_di_SP = 0;
 	Gal[p].n_di_MM = 0;
 	Gal[p].n_di = 0;
 	Gal[p].merger_during_DI = -1;
 	Gal[p].Type_merger = -1;
 	Gal[p].Bulge_Minor_merger = 1e-9;
 	Gal[p].Bulge_Major_merger = 1e-9;
 	Gal[p].Bulge_di = 1e-9;
 	Gal[p].t_last_insta_minor_induced = -1000.;
 	Gal[p].type_v_sigma = -1.; // Type of bulge formation 1: Major merger, 2: Minor merger, 3: Secular
 	Gal[p].v_sigma_now = -1.0;
 	Gal[p].v_sigma = -1.0;
    }
 
  else{ Gal[p].induced = -1; Gal[p].t_dyn_disk_merger = 0.;
        Gal[p].t_dyn_last_disk = 0.; Gal[p].t_last_merger = 0.;
        Gal[p].t_last_disk = 0.;}
}	


void join_N_merger_DI(int t, int p)
{
  // t --> Central
  // p --> Satellite
  
  // Here I have the compute of the number of XXXX trought the whole hostory of the galaxy
  // Note that the remnent galaxy inherit all the event of the satellite and central galaxy

  Gal[t].n_di_SP += Gal[p].n_di_SP;
  Gal[p].n_di_SP = 0;

  Gal[t].n_di += Gal[p].n_di;
  Gal[p].n_di = 0;

  Gal[t].n_di_MM +=Gal[p].n_di_MM;
  Gal[p].n_di_MM = 0;

  Gal[t].n_Min_Mer += Gal[p].n_Min_Mer;
  Gal[p].n_Min_Mer = 0;

  Gal[t].n_Maj_Mer += Gal[p].n_Maj_Mer;
  Gal[p].n_Maj_Mer = 0;

}

void deal_with_SatelliteMass_during_Smooth_accretion_Minor_merger(int t, int p, float tdisk_mass, float tbulge_mass, float pdisk_mass, float pbulge_mass, int Smooth_accretion)
{

  // Here I decide where to put the total mass of the satellite galaxy
  // In case Minor merger --> Bulge + Disk of satellite go to the bulge of the central galaxy
  // In case Smooth accretion --> Bulge + Disk of satellite go to the disk of the central galaxy

  int j;
  if(Smooth_accretion>0)
   {
      transfer_stars(t,"Disk",p,"Disk",1.);
      transfer_stars(t,"Disk",p,"Bulge",1.);
      // Here I update the bulge spin
      if((tdisk_mass + pdisk_mass + pbulge_mass) > 1e-9)
         {
           for(j=0; j<3; j++) Gal[t].StellarSpin[j]=(Gal[t].StellarSpin[j]*tdisk_mass + pdisk_mass*Gal[p].StellarSpin[j] + pbulge_mass*Gal[p].StellarBulgeSpin[j] )/(tdisk_mass + pdisk_mass + pbulge_mass);
           if(isnan(Gal[t].StellarDiskRadius) == 1) {printf("The stellar disk radii is wrong in model_mergers.c : Gal[t].StellarSpin[0] = %E, Gal[t].StellarSpin[1] = %E, Gal[t].StellarSpin[2] = %E, t_disk_mass = %f, pdisk_mass = %f\n",Gal[t].StellarSpin[0],Gal[t].StellarSpin[1],Gal[t].StellarSpin[2],tdisk_mass,pdisk_mass); exit(0);}
         }
    }

  else
    {
      // Here is the fiducual bulge treatement during minor mergers of Henriques et al. 2015
      compute_Sersic_index_of_galaxy(t,(pdisk_mass+pbulge_mass),1);
      if(BulgeFormationInMinorMergersOn)
        {Gal[t].Bulge_Minor_merger += Gal[p].DiskMass; transfer_stars(t,"Bulge",p,"Disk",1.);}
      else
        transfer_stars(t,"Disk",p,"Disk",1.);

      Gal[t].Bulge_Minor_merger += Gal[p].BulgeMass;
      transfer_stars(t,"Bulge",p,"Bulge",1.);
      // Here I update the bulge spin
      if((tbulge_mass + pbulge_mass)> 1e-9)
       {
        for(j=0;j<3;j++) Gal[t].StellarBulgeSpin[j] = (Gal[t].StellarBulgeSpin[j]*tbulge_mass + pbulge_mass*Gal[p].StellarBulgeSpin[j])/(tbulge_mass + pbulge_mass);
        if(isnan(Gal[t].StellarBulgeSpin[0]) == 1 || isnan(Gal[t].StellarBulgeSpin[1]) == 1 || isnan(Gal[t].StellarBulgeSpin[2]) == 1) {printf("The bulge spin is wrong in model_mergers.c : Gal[t].StellarBulgeSpin[0] = %E, Gal[t].StellarBulgeSpin[1] = %E, Gal[t].StellarBulgeSpin[2] = %E, tbulge_mass = %f, pbulge_mass = %f\n",Gal[t].StellarBulgeSpin[0],Gal[t].StellarBulgeSpin[1],Gal[t].StellarBulgeSpin[2],tbulge_mass,pbulge_mass); exit(0);}
       }
    }



}

int determine_if_Smooth_accretion(int t, int p, float mass_ratio, float Mcstar, float Mcbulge, float Mcgas, float Mpstar, float Mpbulge, float Mpgas, float *Mass_added_disk)
{
  // t --> central galaxy
  // p --> satellite galaxy
  // McXXX --> central
  // MpXXX --> satellite
  // This code determines according to the merger ratio and the binding energy if the minor merger is a smooth accretion

  int Smooth_accretion;

  // First: I compute the binding energy of the two interacting systems

  float f_binding, R_sat, R_cent;
  R_cent=(Gal[t].StellarDiskRadius/3.*1.68*(Mcstar-Mcbulge)+Gal[t].GasDiskRadius*Mcgas/3.*1.68+Gal[t].GasDiskRadius*Mpgas/3.*1.68)/(Mpgas + Mcgas+(Mcstar-Mcbulge));
  R_sat = (Gal[p].StellarDiskRadius/3.*1.68*(Mpstar-Mpbulge)+Gal[p].BulgeSize*Mpbulge)/(Mpstar);
  if( (Mcgas+(Mcstar-Mcbulge)) < 1e-9 || Mpstar < 1e-9) f_binding = 9999.;
  else f_binding = (R_cent/R_sat) * pow((Mpstar*Mpstar)/(Mpgas + Mcgas + (Mcstar-Mcbulge)),2.0);

  // Second: Now I check if the conditions for a smooth accretion are satisfied

  if((mass_ratio<ThreshMajorMerger) &&  (f_binding<ThreshSmoothaccretion))
     {Smooth_accretion = +1; *Mass_added_disk = Mpstar; Gal[p].Smooth_minor_merger = 1.0;}
  else
     {Smooth_accretion = -1; *Mass_added_disk = 1e-10;
      if(mass_ratio<ThreshMajorMerger) Gal[p].Smooth_minor_merger = -1.0;
      else Gal[p].Smooth_minor_merger = -2.0;
     }

  return Smooth_accretion;
}

void set_galaxy_merger_properties(int p, int nstep, int Smooth_accretion, float mass_ratio)
{
    // Here I assigne values to the properties that I will use for the bulge (pseudobulge) tracking, i.e merger or secular DI     

  float dt;
  dt = ( (NumToTime(Gal[p].SnapNum) - NumToTime(Gal[p].SnapNum+1) ) /STEPS );
  Gal[p].t_last_merger = (NumToTime(Gal[p].SnapNum)) - ( 1.0 * nstep * dt );
  Gal[p].Mstellar_sec = 0.0;
  Gal[p].Snap_last_merger = Gal[p].SnapNum;
  Gal[p].step_last_merger = nstep;
  Gal[p].z_last_merger = redshift_at_t(ZZ[Gal[p].SnapNum],nstep,dt); //ZZ[Gal[p].SnapNum];


  if(mass_ratio < ThreshMajorMerger){
       Gal[p].t_last_MINOR = Gal[p].t_last_merger/ (Hubble_h/UnitTime_in_Megayears); // Myr
       Gal[p].mass_ratio_mm = mass_ratio; // mass ratio of the las minor merger
       Gal[p].n_Min_Mer++;
       if(Smooth_accretion<0.) Gal[p].z_last_Minor_merger = Gal[p].z_last_merger;//ZZ[Gal[p].SnapNum];
       else Gal[p].z_last_Smooth_accretion = Gal[p].z_last_merger;
       if(Gal[p].Type == 0) Gal[p].t_dyn_disk_merger =  (Gal[p].StellarDiskRadius/3.0)/Gal[p].Vmax;
       else Gal[p].t_dyn_disk_merger = (Gal[p].StellarDiskRadius/3.0)/Gal[p].InfallVmax;
       Gal[p].happened_minor_merger = +1;
       Gal[p].merger_during_DI = check_if_previous_step_DI(p,nstep,Gal[p].t_last_merger);
       Gal[p].sfr_merger = Gal[p].Sfr;
       Gal[p].is_the_first_after_minor_Merger = +1;
       Gal[p].Type_merger = Gal[p].Type;
       //checks
       if(Gal[p].t_dyn_disk_merger<=0.)
               {printf("Error in set_galaxy_merger_properties: The t_dyn_disk_merger is 0 or negative\n"); exit(0);}
       }

  else{
        Gal[p].t_last_MAJOR = Gal[p].t_last_merger / (Hubble_h/UnitTime_in_Megayears); // Myr
        Gal[p].z_last_Major_merger = Gal[p].z_last_merger;
        Gal[p].n_Maj_Mer++;
        Gal[p].mass_ratio_Mm = mass_ratio; // mass ratio of the last major merger
        Gal[p].increment_due_to_burst = -99;
        Gal[p].t_last_merger = 0.;
        Gal[p].t_last_disk = 0.;
        Gal[p].t_dyn_disk_merger = 0.;
        Gal[p].t_dyn_last_disk = 0.;
        Gal[p].merger_during_DI = -1;
        Gal[p].is_the_first_after_minor_Merger = -1;
      }

}

void  determine_quantities_needed_to_distinguish_DI_induced_secular(int p, float mass_ratio, float Mdisk_before, float deltaT, float New_stars, float Extra_Mass_disk)
{

  // Here I set all the properties needed to determine during the code run, if the DI is secular induced or 
  // minor merger induced
  // NOTE: New_stars is the amount of new stars that are generated after the collisional starburst
  float Mcrit;


  if(mass_ratio<ThreshMajorMerger)
        {
            if (Gal[p].Type != 0)
                Mcrit = Gal[p].InfallVmax * Gal[p].InfallVmax * (Gal[p].StellarDiskRadius/1.0) / ( G_GRAVITY * pow(Efstathiou_crit, 2.) );
            else
                Mcrit = Gal[p].Vmax * Gal[p].Vmax * (Gal[p].StellarDiskRadius/1.0) / (G_GRAVITY * pow(Efstathiou_crit,2.) );
            Gal[p].Mcrit_merger = Mcrit;
            Gal[p].increment_due_to_burst = ( Extra_Mass_disk + ( (1.0 - RecycleFraction) * New_stars) )/ ( Gal[p].Mcrit_merger - Mdisk_before); // TODO add the mass of the satellite galaxie in case
            Gal[p].sfr_Burst_merger = (New_stars / deltaT);
            if(New_stars <= 0.0) Gal[p].Mburst_merger = -2.;
            else  Gal[p].Mburst_merger = (1.0 - RecycleFraction) *New_stars;
            if(fabs(Gal[p].Mcrit_merger - Mdisk_before )<=1e-4)
               {
                   if(New_stars == 0) Gal[p].increment_due_to_burst = 1e-11;
                   else Gal[p].increment_due_to_burst = -3.;
               }
        }
   else {Gal[p].increment_due_to_burst = -99.0; Gal[p].Mburst_merger = -1; Gal[p].Mcrit_merger = -1; Gal[p].sfr_Burst_merger  = -1;}

}



int check_if_previous_step_DI(int p, int step, float time_now)
{
  int dummy, dummy_inside_tdyn, dummy_consecutive;

  // Given that in LGalaxies is divided in substeps, one disk instability can expand several consecutive substeps. Therefore
  // we performe here acode to determine if two DI in different substeps are the same disk instability or are two different
  

  if(Gal[p].Snap_last_DI < 0) dummy =  -1;
  else
       {
         // 1) First I check if the DI is a consecutive one
         if(step == 0)
                {
                   if( ( (Gal[p].SnapNum-1) == Gal[p].Snap_last_DI)  && ( Gal[p].step_last_DI == (STEPS-1) ) ) dummy_consecutive =  1;
                   else if( (Gal[p].SnapNum == Gal[p].Snap_last_DI)  && ( Gal[p].step_last_DI == 0 ) ) dummy_consecutive =  1;
                   else dummy_consecutive =  -1;
                }
        else
                {
                   if( (Gal[p].SnapNum == Gal[p].Snap_last_DI) && ( Gal[p].step_last_DI == (step-1) ) ) dummy_consecutive =  1;
                   if( (Gal[p].SnapNum == Gal[p].Snap_last_DI) && ( Gal[p].step_last_DI == step ) ) dummy_consecutive =  1;
                   else dummy_consecutive =  -1;
                }

         // 2) Second I check if the DI is inside of the DI dynamical time
        if(Gal[p].t_last_disk <=0.) dummy_inside_tdyn = -1;
        else{
                if( time_now < (Gal[p].t_last_disk - (Gal[p].t_dyn_last_disk * DelayBetweenDI ) )  ) dummy_inside_tdyn = -1;
                else dummy_inside_tdyn = +1;
            }

        // 3) I join the two conditions
        if( (dummy_inside_tdyn < 0 ) && ( dummy_consecutive < 0 ) ) dummy = -1;
        else dummy = +1;
        //if(dummy <0) printf("Dummy negative!\n");

       }
   return dummy;
}



void check_if_DI_is_minor_merger_induced(int p, double dt_merger, double time_DI, float Bulge_after, float Bulge_before, float t_merger, int n_step)
{
    // Here we take into account if the bulge is formed by
    // a minor merger consecuence or a secular process

    float to_be_consider_SAME_DI;
    float v_sigma_before;
    to_be_consider_SAME_DI = DelayBetweenDI;
    int induced_before;
    int dummy;
    induced_before = Gal[p].induced;
    v_sigma_before = Gal[p].v_sigma_now ;

    if(Gal[p].counter_Delta_M <=1 || Gal[p].is_the_first_after_minor_Merger > 0)
        {
           if(Gal[p].is_the_first_after_minor_Merger > 0 ) {
                if(Gal[p].induced > 0)  Gal[p].induced = +1;
                else {
                          if(dt_merger <= 1e-9) {
                                   if( Gal[p].sfr_Burst_merger > Gal[p].sfr_merger ) {Gal[p].induced = +1;}// printf("Induced entered from sfr_inst > sfr_burst. Increment = %f\n",Gal[p].increment_due_to_burst);}
                                   else Gal[p].induced = -1;
                                }
                          else {
                                   if(Gal[p].increment_due_to_burst < 0. && Gal[p].increment_due_to_burst!=-3.0){printf("Weird!, f_m = %f Snp_merger = %d, step_merger = %d, Snp_DI = %d, step_DI = %d\n",
                                                        Gal[p].increment_due_to_burst, Gal[p].Snap_last_merger, Gal[p].step_last_merger,  Gal[p].Snap_last_DI, Gal[p].step_last_DI );}
                                   if( (Gal[p].increment_due_to_burst > 0.5) && ((dt_merger/Gal[p].t_dyn_disk_merger)<= DelayMergerDI) && (dt_merger>=0.)) {Gal[p].induced = +1;}// printf("Induced entered from normal DI\n"); }
                                   else Gal[p].induced = -1;
                               }
                        }
                }
          else {
                 if( (Gal[p].t_dyn_last_disk > 0.) && (((Gal[p].t_last_disk - time_DI)/Gal[p].t_dyn_last_disk) <= to_be_consider_SAME_DI) && ( (Gal[p].t_last_disk - time_DI) >= 0.) ) Gal[p].induced = Gal[p].induced;
                 else Gal[p].induced = -1;
               }
          Gal[p].is_the_first_after_minor_Merger = -1;
          Gal[p].increment_due_to_burst= -99;
        }
    // Finish the asignation of minor merger induced or not minor merger induced
    if( Gal[p].induced > 0 ) Gal[p].Bulge_Minor_merger += (Bulge_after - Bulge_before);
    else Gal[p].Bulge_di += (Bulge_after - Bulge_before);
   // This is to save when was the last minor merger which induced a disk instability
   if(induced_before != Gal[p].induced && Gal[p].induced>0) Gal[p].t_last_insta_minor_induced = Gal[p].t_last_MINOR;



  // Lets compute the bulge Sersic index
  float Delta_Mbulge;
  Delta_Mbulge = Bulge_after-Bulge_before;
  if(isnan(Delta_Mbulge) == 1 || isinf(Delta_Mbulge) == 1) Delta_Mbulge = 1e-9;
  if(Gal[p].induced > 0)  compute_Sersic_index_of_galaxy(p,Delta_Mbulge,1);
  else compute_Sersic_index_of_galaxy(p,Delta_Mbulge,2);


  dummy = -1;
  if(Gal[p].type_v_sigma<0)
    {
     dummy = -1;
     if(Gal[p].induced > 0)
        {Gal[p].type_v_sigma = 1; Gal[p].v_sigma_now = altura_number(Gal[p].type_v_sigma,p,n_step);}
     else
        {Gal[p].type_v_sigma = 2; Gal[p].v_sigma_now = altura_number(Gal[p].type_v_sigma,p,n_step);}
    }
   else
    {
        if(Gal[p].induced > 0)
           {
               dummy = 0;
               Gal[p].n_di_MM++;
               Gal[p].type_v_sigma = 1;
               Gal[p].v_sigma_now = altura_number(Gal[p].type_v_sigma,p,n_step);
               Gal[p].v_sigma_now = ((Bulge_after - Bulge_before)*Gal[p].v_sigma_now + Gal[p].Bulge_di*v_sigma_before)/((Bulge_after - Bulge_before) + Gal[p].Bulge_di);
               if(isnan(Gal[p].v_sigma_now) == 1) Gal[p].v_sigma_now = altura_number(Gal[p].type_v_sigma,p,n_step); // Usually happens when 1) Gal[p].Bulge_di = 0 and 2) the ammount of mass trasnfered during the DI is almost ~ 0 so: Bulge_after ~=Bulge_after
           }
         else
           {
              if( ( (Gal[p].t_last_disk - time_DI) < 0 ) )
                {
                   Gal[p].n_di_SP++;
                   dummy = 1;
                   Gal[p].type_v_sigma = 2;
                   Gal[p].v_sigma_now = ((Bulge_after - Bulge_before)*Gal[p].v_sigma_now + Gal[p].Bulge_di*v_sigma_before)/((Bulge_after - Bulge_before) + Gal[p].Bulge_di);
                   if(isnan(Gal[p].v_sigma_now) == 1) Gal[p].v_sigma_now = altura_number(Gal[p].type_v_sigma,p,n_step); // Usually happens when 1) Gal[p].Bulge_di = 0 and 2) the ammount of mass trasnfered during the DI is almost ~ 0 so: Bulge_after ~=Bulge_after
                   //Gal[p].v_sigma_now = altura_number(Gal[p].type_v_sigma);
                }
               else
                {
                    if( (Gal[p].counter_Delta_M <= 1) &&  (((Gal[p].t_last_disk - time_DI)/Gal[p].t_dyn_last_disk) > to_be_consider_SAME_DI) )
                     {
                           dummy = 2;
                           Gal[p].n_di_SP++;
                           Gal[p].type_v_sigma = 2;
                           Gal[p].v_sigma_now = ((Bulge_after - Bulge_before)*Gal[p].v_sigma_now + Gal[p].Bulge_di*v_sigma_before)/((Bulge_after - Bulge_before) + Gal[p].Bulge_di);
                           if(isnan(Gal[p].v_sigma_now) == 1) Gal[p].v_sigma_now = altura_number(Gal[p].type_v_sigma,p,n_step); // Usually happens when 1) Gal[p].Bulge_di = 0 and 2) the ammount of mass trasnfered during the DI is almost ~ 0 so: Bulge_after ~=Bulge_after
                           //Gal[p].v_sigma_now = altura_number(Gal[p].type_v_sigma);
                      }
                    else {Gal[p].v_sigma_now = Gal[p].v_sigma_now;  Gal[p].type_v_sigma = Gal[p].type_v_sigma; dummy = 3;}
                }
            }
     }
 if(isnan(Gal[p].v_sigma_now) == 1 || Gal[p].v_sigma_now < 0) {printf("Error in check_if_DI_is_minor_merger_induced(). The value of Galaxy v_sigma_now is nan or negative! dummy = %d, Bulge_after = %E, Bulge_before = %E, Gal[p].Bulge_di = %E, ((Bulge_after - Bulge_before) + Gal[p].Bulge_di) = %E\n",dummy, Bulge_after, Bulge_before, Gal[p].Bulge_di,(Bulge_after - Bulge_before) + Gal[p].Bulge_di); exit(0);}


}




double altura_number(int type_now, int p, int n_step)
{
double p_x,x,y,x_Max,y_Max,v_sigma, counter, xarr[100], yarr[100];
int acepted, i;
float aa, bb, cc, xx, xx0, zz, sigma, alpha, med, beta, dt;
aa = 0.6949;
bb = 0.4548;
cc = 8.5255;
beta =  3.2126;
alpha = -0.8680;
acepted = -1;
x_Max = 0;
y_Max = 0;

if(type_now==0) { x_Max = 0; y_Max = type_of_galaxy_hybrid_model(type_now, x_Max,p,n_step);}
else if(type_now==1) { x_Max = 0.5196; y_Max = type_of_galaxy_hybrid_model(type_now, x_Max,p,n_step);}
//else if(type_now==1) { x_Max = 0.49777175;  y_Max = type_of_galaxy_hybrid_model(type_now, x_Max,p,n_step);}
else if(type_now==2)
        {
            dt = (NumToTime(Gal[p].SnapNum) - NumToTime(Gal[p].SnapNum + 1))/STEPS;
            zz = redshift_at_t(ZZ[Gal[p].SnapNum],n_step,dt);
            if(zz<0.2) zz = 0.2;
            sigma = 0.34;
            xx = log10(1.0+zz);
            xx0 = bb*(log10((Gal[p].DiskMass + Gal[p].BulgeMass)*1e10/Hubble_h) - cc);
            med = aa * pow(xx/xx0,alpha) * pow(1.0+pow(xx/xx0,beta),-alpha/beta);
            if(log10((Gal[p].DiskMass + Gal[p].BulgeMass)*1e10/Hubble_h) < cc ) med = 1.0;
            if(isnan(med) == 1) {printf("Error in altura_number() : type_of_galaxy_hybrid_model, med is nan\n"); exit(0);}
            counter = 0.01;
            for(i=0;i<100;i++)
              {
                    xarr[i] = counter;
                    counter = counter + 0.1;
                    yarr[i] = exp(-pow(log(xarr[i])-med,2.0)/(2.0*pow(sigma,2.0)))/(xarr[i]*sigma*pow(2.0*M_PI,0.5));
              }
            i = 0;
            max_and_index(yarr,100,&y_Max,&i);
            x_Max = xarr[i];
      }

while(acepted < 0)
{
        x = 10.5*((double)rand()/((double)RAND_MAX+1));
        //p_x = type_of_galaxy(type_now, x, p);
        p_x = type_of_galaxy_hybrid_model(type_now, x, p,n_step);
        y = y_Max*((double)rand()/((double)RAND_MAX+1));
        if(y<=p_x){ v_sigma = x; acepted = +1; }
}
return v_sigma;
}


void max_and_index(double arr[200], int dim, double *xmax, int *indx_xmax)
{
  double maximum;
  int c, location;
  maximum = arr[0];
  location = 0;
  for (c = 1; c < dim; c++)
  {
    if (arr[c] > maximum)
    {
       maximum  = arr[c];
       location = c;
    }
  }
 *xmax = maximum;
 *indx_xmax = location;
}


double type_of_galaxy(int tt, double x, int p)
{
double val;
  if(tt == 0) val =  14.59*exp(-3.98*pow(x,0.96)); // elliptical
  else if(tt == 1) val = (2.10/0.24)*pow(x,0.24)* exp(-0.5*pow((x-0.47)/0.24,2.)); // bulge MERGER spiral
  else if(tt == 2) val = 0.083*(x+0.53)/(pow(x,2.) -1.38*x + 0.5); // bulge SECULAR (DI) spirals
  else {printf("Error in type_of_galaxy(): Not type recognized in type_of_galaxy\n"); exit(0);}
return val;
}


double type_of_galaxy_hybrid_model(int tt, double x, int p, int n_step)
{
  // This hybrid model assumes that the pseudobulge v/sigma value is similar to the one of the disk
  // for furhter details check Izquierdo-Villalba et al. 2020 and Sesana et al. 2016
  double val;
  float aa, bb, cc, xx, xx0, zz, sigma, alpha, med, beta, dt;
  aa = 0.6949;
  bb = 0.4548;
  cc = 8.5255;
  beta =  3.2126;
  alpha = -0.8680;
  dt = (NumToTime(Gal[p].SnapNum) - NumToTime(Gal[p].SnapNum + 1))/STEPS;
  zz = redshift_at_t(ZZ[Gal[p].SnapNum],n_step,dt);
  if(zz<0.2) zz = 0.2;
  //printf("x == %f\n",x);
  if(tt == 0) val =  14.59*exp(-3.98*pow(x,0.96)); // elliptical
  else if(tt == 1) val = 42.24*pow(x,1.29)* exp(-2.95*pow(x,2.81)); // bulge MERGER spiral. The first factor is the normalization
  else if(tt == 2)
        {
           sigma = 0.34;
           //val = 0.083*(x+0.53)/(pow(x,2.) -1.38*x + 0.5); // bulge SECULAR (DI) spirals
           xx = log10(1.0+zz);
           xx0 = bb*(log10((Gal[p].DiskMass + Gal[p].BulgeMass)*1e10/Hubble_h) - cc);
           med = aa * pow(xx/xx0,alpha) * pow(1.0+pow(xx/xx0,beta),-alpha/beta);
           if(log10((Gal[p].DiskMass + Gal[p].BulgeMass)*1e10/Hubble_h) < cc ) med = 1.0;
           if(isnan(med) == 1) {printf("Error in type_of_galaxy_hybrid_model(): med is nan, tt = %d, x = %f, p = %d, xx = %f, xx0 = %f\n",tt,x,p,xx,xx0); exit(0);}
           val = exp(-pow(log(x)-med,2.0)/(2.0*pow(sigma,2.0)))/(x*sigma*pow(2.0*M_PI,0.5));
        }
  else {printf("Error in type_of_galaxy(): Not type recognized in type_of_galaxy\n"); exit(0);}

return val;
}



float altura_number_Sersic(int Bulge_type)
{
  float p_x,x,y,x_Max,y_Max,n_sersic;
  int acepted;
  acepted = -1;
  x_Max = 0;
  // Bulge_type == 0  Elliptical
  // Bulge type == 1  Clasical bulges
  // Bulge type == 2  Pseudobulges
  if(Bulge_type==0) x_Max = 3.885;//0.116146972585;
  else if(Bulge_type ==1) x_Max = 3.385;//0.0746627141931;
  else if(Bulge_type==2) x_Max = 1.635;//0.0945435956907;
  //if(Bulge_type==0)  return 1.0;
  //srand(time(NULL));
  y_Max = type_of_galaxy_Sersic(Bulge_type, x_Max);
  //printf("y_Max = %f\n",y_Max);
  while(acepted < 0)
  {
          x = 8.0*((float)rand()/((float)RAND_MAX+1));
          p_x = type_of_galaxy_Sersic(Bulge_type, x);
          y = y_Max*((float)rand()/((float)RAND_MAX+1));
          if(y<=p_x){ n_sersic = x; acepted = +1; }
  }
  return n_sersic;
}

float type_of_galaxy_Sersic(int tt, float z)
{
  double val, A, z0, beta;
  //printf("x == %f\n",x);
  if(tt == 0) {A = 1.15498593; z0 = 4.24177435; beta = 5.75390158;}
  else if(tt == 1) {A = 0.60460942; z0 = 2.46891546; beta = 1.87983007;}
  else if(tt == 2) {A = 0.02104979; z0 = 0.16567257; beta = 0.71373334;}
  else {printf("Error in type_of_galaxy(): Not type recognized in type_of_galaxy\n"); exit(0);}
  val = A*pow(z/z0,3.5) * exp(-1.0*pow(z/z0,beta)); // Extracted from 0 . Le Fevre et al. 2013
  return val;
}


void compute_Sersic_index_of_galaxy(int p, float Mbulge_new, int Bulge_type)
{
  float n_random;
  //printf("Sersic Index Before = %f ", Gal[p].Sersic_index);
  n_random = altura_number_Sersic(Bulge_type);
  if(Gal[p].Sersic_index < 0. || Gal[p].BulgeMass < 1e-9 ||  Bulge_type < 1) Gal[p].Sersic_index = n_random;
  else Gal[p].Sersic_index = ( ( Gal[p].Sersic_index * Gal[p].BulgeMass ) + ( n_random * Mbulge_new ) ) / (Gal[p].BulgeMass + Mbulge_new);
  //if(Bulge_type == 0) Gal[p].Sersic_index = altura_number_Sersic(2);
  //printf("BulgeType = %d, n_random = %f, Gal[p].Sersic_index = %f\n",Bulge_type, n_random, Gal[p].Sersic_index);
}


void  bulgesize_from_Disipative_merger(double mass_ratio,int merger_centralgal,int p, double Mcstar,double Mcbulge,double Mcgas, double Mpstar,double Mpbulge,double Mpgas,double frac, int n_step, int filenr)
{
 
 // This function computes the bulge size when we include the disipation during merger

  double Mc,Rc;
  double Mp,Rp;
  double fint,c;
  double fgas;
  double Crad;
  double disipation_factor, B_T_1, B_T_2;
  float M_DM;

  c = 0.5; // just in case
  Crad = 0.; //just in case
  if((Mcstar+Mcgas) > (Mpstar + Mpgas)){B_T_1 = Mcbulge/Mcstar; B_T_2=Mpbulge/Mpstar;}
  else {B_T_1 = Mpbulge/Mpstar; B_T_2=Mcbulge/Mcstar;}

  if(isnan(B_T_1) == 1) B_T_1 = 0.;
  if(isnan(B_T_2) == 1) B_T_2 = 0.;

  Crad = 2.75;
  c=0.5;
  fint=0.;
  fgas = (Mcgas + Mpgas) / ( Mcgas + Mpgas + Mpstar + Mcstar );
  if(fgas > 1.) {printf("fgas in bulge_size() model mergers is > 1 (%f)\n",fgas); exit(0);}
  disipation_factor = (Crad*fgas);

  /* calculate radius for the object that will form the new bulge - Rc and Rp */
  /* Minor Merger */
   if(mass_ratio < ThreshMajorMerger) {
    /* In a minor merger only the stars of the satellite galaxy are moved
     * to the bulge of the central galaxy, therefore only stellar
     * components are used to compute radius and masses. */
    frac=0.0;
    disipation_factor = 0.;
    /* in a minor merger only consider the bulge mass of the central galaxy */
    Mc=Mcbulge;
    Rc=Gal[merger_centralgal].BulgeSize;
    //M_DM = DM_mass_inside_radi(merger_centralgal,Rc);
    /* and stellarmass of satellite*/
    Mp=Mpstar;
    if (Mp >0.0)
      Rp=(Gal[p].StellarDiskRadius/3.*1.68*(Mpstar-Mpbulge)+Gal[p].BulgeSize*Mpbulge)/Mpstar;
    else
      Rp=0.0;
  }

  /* Major Merger */
  else {
    /* on a major merger both stellar and gas (after a burst) components
     * from the final bulge and need to be considered */
    /* Mc = bulge mass + burst of central*/

    Mc=Mcstar+frac*Mcgas;
    if (Mc > 0.0)
      { Rc=(Gal[merger_centralgal].StellarDiskRadius/3.*1.68*(Mcstar-Mcbulge)+Gal[merger_centralgal].BulgeSize*Mcbulge+Gal[merger_centralgal].GasDiskRadius*Mcgas/3.*1.68)/(Mcgas+Mcstar); M_DM = DM_mass_inside_radi(merger_centralgal,Rc,n_step,filenr); }
    else
      {Rc=0.0; M_DM = 0.0;}
    Mc = Mc + Mcgas + 2.0*M_DM;
    /* and satellite Mp */
    Mp=Mpstar+frac*Mpgas;
    if (Mp > 0.0){
      Rp=(Gal[p].StellarDiskRadius/3.*1.68*(Mpstar-Mpbulge)+Gal[p].BulgeSize*Mpbulge+Gal[p].GasDiskRadius*Mpgas/3.*1.68)/(Mpgas+Mpstar); M_DM = DM_mass_inside_radi(p,Rp,n_step,filenr);}
    else{
      Rp=0.0; M_DM = 0.0;}
    Mp = Mp + Mpgas;
  }

  if(Rp>0. && Rp<1.e-8)
      Rp=1.e-8;
  if(Rc>0. && Rc<1.e-8)
      Rc=1.e-8;
  /* If both original radius are bigger then 0 then this is Eq. 33 in Guo 2010
   * solved for R_new,bulge with all terms divided by G and C. */
  if(Rc >= 1.e-8 && Rp >= 1.e-8)
      Gal[merger_centralgal].BulgeSize=(Mp+Mc)*(Mp+Mc)/( ( (1.+disipation_factor)*(Mp*Mp/Rp+Mc*Mc/Rc) ) + (fint/c*Mp*Mc/(Rp+Rc)) );
  if(Rc >= 1.e-8 && Rp <= 1.e-8)
      Gal[merger_centralgal].BulgeSize=(Mp+Mc)*(Mp+Mc)/( ( (1.+disipation_factor)*(Mc*Mc/Rc) ) + (fint/c*Mp*Mc/(Rp+Rc)) );
  if(Rc <= 1.e-8 && Rp <=1.e-8)
      Gal[merger_centralgal].BulgeSize=0.0;
  if(Rc <= 1.e-8 && Rp >= 1.e-8)
      Gal[merger_centralgal].BulgeSize=(Mp+Mc)*(Mp+Mc)/( ( (1.+disipation_factor)*(Mp*Mp/Rp) ) + (fint/c*Mp*Mc/(Rp+Rc)) );

  if(isnan(Gal[merger_centralgal].BulgeSize)==1 && Gal[merger_centralgal].SnapNum<58) {
    printf("In bulgesize_from_dissipative_merger, the bulge size is nan! check it!\n"); exit(0);}
  if( Gal[merger_centralgal].BulgeSize == 0.) printf("Error\n");

  if ((Mp+Mc > 0.0 && Gal[merger_centralgal].BulgeSize == 0.0 && mass_ratio > 0.) || (Mp+Mc == 0.0 && Gal[merger_centralgal].BulgeSize> 0.0))
    {
      char sbuf[1000];
      printf("halonr=%d, merger_centralgal %d\n\n", Gal[merger_centralgal].HaloNr, merger_centralgal);

      printf("New Bulge Mass from Central (MC)=%e\n New Bulge Mass from Satellite (Mp)=%e\n NewBulge size=%e\n\n",
             Mp,Mc,Gal[merger_centralgal].BulgeSize);

      if(mass_ratio < ThreshMajorMerger)
        {
          printf("minor merger, new mass=original mass\n");
          printf("New Bulge Mass From Central (Mc)   = Mcbulge = %f\n",Mcbulge);
          printf("New Bulge Mass From Satellite (Mp) = Mpstar  = %f\n",Mpstar);
        }
      else
        {
          printf("New Bulge From Central (Mc)   = Mcstar+frac*Mcgas = %f+%f*%f\n", Mcstar, frac, Mcgas);
          printf("New Bulge From Satellite (Mp) = Mpstar+frac*Mpgas = %f+%f*%f\n", Mpstar, frac, Mpgas);
        }

      printf("BulgeSize from Central (Rc)=%e\nBulgeSize from Satellite (Rp)=%e\nmass ratio=%f\n\n",Rc,Rp, mass_ratio);

      printf("the following masses don't tell a lot because things have been merged already!!!\n");
      printf("    sat: BulgeMass=%0.7f, BulgeSize=%0.7f, GasMass=%0.7f, GasSize=%0.7f, DiskMass=%0.7f StellarSize=%0.7f \n",
             Gal[p].BulgeMass, Gal[p].BulgeSize, Gal[p].ColdGas, Gal[p].GasDiskRadius,
             Gal[p].DiskMass, Gal[p].StellarDiskRadius);
      printf(   "central: BulgeMass=%0.7f, BulgeSize=%0.7f, GasMass=%0.7f, GasSize=%0.7f, DiskMass=%0.7f StellarSize=%0.7f \n",
                Gal[merger_centralgal].BulgeMass,  Gal[merger_centralgal].BulgeSize, Gal[merger_centralgal].ColdGas,
                Gal[merger_centralgal].GasDiskRadius, Gal[merger_centralgal].DiskMass, Gal[merger_centralgal].StellarDiskRadius);

      sprintf(sbuf,"\n bulgesize wrong in merger");
      terminate(sbuf);
      exit(0);
  }

}


void Generate_file_to_link_DI(int p, float tim_DI)
{
  FILE *fileout;
    // Generate the file for the delay between ALL (secular + induced) DI. If you want only the Secular add in the if: Gal[p].induced < 0
      if( (Gal[p].counter_Delta_M == 1 ) && ( (Gal[p].t_last_disk - tim_DI) >0. ) )
        {
              //fileout = fopen("/global/users/dizquierdo/thesis/LGalaxies/DI/MR/epsilon_1_5/All_DI","ab");
              #ifndef MRII
                   fileout = fopen("/global/users/dizquierdo/thesis/LGalaxies/DI/MR/epsilon_1_5/All_DI","ab");
              #else
                   fileout = fopen("/global/users/dizquierdo/thesis/LGalaxies/DI/MRII/epsilon_1_5/All_DI","ab");
              #endif
              float arr1[2] = {(Gal[p].t_last_disk - tim_DI)/(Hubble_h/UnitTime_in_Megayears),Gal[p].t_dyn_last_disk/(Hubble_h/UnitTime_in_Megayears)};
              fwrite(arr1,sizeof(arr1),1,fileout);
              fclose(fileout);
        }
}
void Generate_file_number_density_DI(int p, float tim_DI)
{
  FILE *fileout;
  int dummy;
        // Generate the file for the DI number density
        //if( ((Gal[p].t_last_disk - tim_DI)/Gal[p].t_dyn_last_disk) <= DelayBetweenDI )
        if( ( ((Gal[p].t_last_disk - tim_DI)/Gal[p].t_dyn_last_disk) <= DelayBetweenDI) && (Gal[p].t_last_disk - tim_DI)>0.)  dummy = 0.;
        else
        {
              #ifndef MRII
                   fileout = fopen("/global/users/dizquierdo/thesis/LGalaxies/DI/MR/epsilon_1_5/Number_density_DI","ab");
              #else
                   fileout = fopen("/global/users/dizquierdo/thesis/LGalaxies/DI/MRII/epsilon_1_5/Number_density_DI","ab");
              #endif
              float arr1[4] = {Gal[p].SnapNum,ZZ[Gal[p].SnapNum],log10((Gal[p].BulgeMass + Gal[p].DiskMass)/(Hubble_h*1.e-10)),Gal[p].induced};
              fwrite(arr1,sizeof(arr1),1,fileout);
              fclose(fileout);
        }
}
void Generate_file_DI_Dleta_M(int p, float tim_DI, int induced_b, int induced_a, float Delta_M_transfere, int step)
{

  // This function generates the files to determine the amount of mass that has gone form the disk to the bulge during a DI (i.e can be several DI events given that I link them) 
  FILE *fileout;
  //int dummy;
        if( ( ((Gal[p].t_last_disk - tim_DI)/Gal[p].t_dyn_last_disk) <= DelayBetweenDI) && (Gal[p].t_last_disk - tim_DI)>0. && (induced_b==induced_a) ) Gal[p].Delta_M += Delta_M_transfere;
        else
        {
              if((Gal[p].t_last_disk - tim_DI)<0.)
                 {
                        Gal[p].snp_started_DI = Gal[p].SnapNum;
                        Gal[p].M_stellar_started_DI = (Gal[p].BulgeMass + Gal[p].DiskMass);
                        Gal[p].Delta_M += Delta_M_transfere;
                        Gal[p].step_started_DI = step;
                 }
                else
                {
                        #ifndef MRII
                           fileout = fopen("/global/users/dizquierdo/thesis/LGalaxies/DI/MR/epsilon_1_5/Delta_M_DI","ab");
                        #else
                           fileout = fopen("/global/users/dizquierdo/thesis/LGalaxies/DI/MRII/epsilon_1_5/Delta_M_DI","ab");
                        #endif
                        float time;
                        time = ( (NumToTime(Gal[p].snp_started_DI) - ( 1.0 * Gal[p].step_started_DI * ( (NumToTime(Gal[p].snp_started_DI) - NumToTime(Gal[p].snp_started_DI+1))/STEPS) ))   -
                               (NumToTime(Gal[p].SnapNum) -        ( 1.0 * step * ( (NumToTime(Gal[p].SnapNum) - NumToTime(Gal[p].SnapNum+1))/STEPS) ) ) ) /(Hubble_h/UnitTime_in_Megayears);;
                        if(time<0.) {printf("Generate_file_DI_Dleta_M() error: time is negative\n"); exit(0);}
                        float arr1[7] = {Gal[p].snp_started_DI,
                                        log10((Gal[p].M_stellar_started_DI)/(Hubble_h*1.e-10)),
                                        log10((Gal[p].Delta_M)/(Hubble_h*1.e-10)),
                                        induced_b,
                                        Gal[p].SnapNum,
                                        time,
                                        log10((Gal[p].BulgeMass + Gal[p].DiskMass)/(Hubble_h*1.e-10))};
                        fwrite(arr1,sizeof(arr1),1,fileout);
                        fclose(fileout);
                        Gal[p].snp_started_DI = Gal[p].SnapNum;
                        Gal[p].step_started_DI = step;
                        Gal[p].M_stellar_started_DI = (Gal[p].BulgeMass + Gal[p].DiskMass);
                        Gal[p].Delta_M = 0.;
                        Gal[p].Delta_M += Delta_M_transfere;
                }
        }
}