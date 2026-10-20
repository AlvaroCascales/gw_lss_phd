/*  Copyright (C) <2016>  <L-Galaxies>
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

/** This module deals with the wandering black holes (i.e ejected and orphan) 
 *  All the physics included can be found summarized in Izquierdo-Villalba et al. 2020 (and references therein)
 *  IMPORTANT NOTE -- The orbits of ejected black holes are going to be assumed radial
 *                     while the orbits of the orphan black hole are going to be tracked in the 3D cartesian coordinates
 *                 -- When a galaxy becomes Type 2 we transfer all its wandering black holes to the new central galaxy (Type 0)
*/

#ifdef MODELWANDERINGBHS
void restart_the_wandering_BH(int p, int indx_wBH)
{
  // This function re-init a wandering black hole spot
  // if indx_wBH is negative I re-init all the slots of wBHs in the galaxy
  // in case indx_wBH is positive, I only re-init the slot whose position coincides with indx_wBH

  int i,j;
  if(indx_wBH<0)
   {
       for(i=0; i<WBHS; i++)
        {  Gal[p].WBlHo[i].e_BlackHoleMass = 1e-15; Gal[p].WBlHo[i].e_Mass_to_acc = 1e-15; Gal[p].WBlHo[i].e_QuasarLum = 0.0;
           Gal[p].WBlHo[i].e_spin = 0.0;  Gal[p].WBlHo[i].e_eta_BH = 0.0; Gal[p].WBlHo[i].e_fEDD = 0.0;
           Gal[p].WBlHo[i].e_M_Res_0 = 0.0; Gal[p].WBlHo[i].e_MBH_0 = 0.0; Gal[p].WBlHo[i].e_MBH_peak = 0.0;
           Gal[p].WBlHo[i].e_MBH_second_phase = 0.0; Gal[p].WBlHo[i].e_t_acc_quiescent_phase = 0.0;
           for(j=0; j<3; j++) {Gal[p].WBlHo[i].e_BH_pos[j] = 3.e-4; Gal[p].WBlHo[i].e_BH_vel[j] = 1e-3;}
           Gal[p].WBlHo[i].e_v_sigma_now = -1.;
           //Gal[p].WBlHo[i].to_be_deleted = -1.0;  
	   Gal[p].WBlHo[i].time_since_recoil = 0.0;
           Gal[p].WBlHo[i].v_kick = -1.0;   Gal[p].WBlHo[i].z_ejected_WPHs = -1.0;
           Gal[p].WBlHo[i].disrupted_galaxy = -1.0;  Gal[p].WBlHo[i].v_escape = -1.0;
           Gal[p].WBlHo[i].inherit = -1.0;
        }
    }
  else
    {
        Gal[p].WBlHo[indx_wBH].e_BlackHoleMass = 1e-15; Gal[p].WBlHo[indx_wBH].e_Mass_to_acc = 1e-15; Gal[p].WBlHo[indx_wBH].e_QuasarLum = 0.0;
        Gal[p].WBlHo[indx_wBH].e_spin = 0.0;  Gal[p].WBlHo[indx_wBH].e_eta_BH = 0.0; Gal[p].WBlHo[indx_wBH].e_fEDD = 0.0;
        Gal[p].WBlHo[indx_wBH].e_M_Res_0 = 0.0; Gal[p].WBlHo[indx_wBH].e_MBH_0 = 0.0; Gal[p].WBlHo[indx_wBH].e_MBH_peak = 0.0;
        Gal[p].WBlHo[indx_wBH].e_MBH_second_phase = 0.0; Gal[p].WBlHo[indx_wBH].e_t_acc_quiescent_phase = 0.0;
        for(j=0; j<3; j++) {Gal[p].WBlHo[indx_wBH].e_BH_pos[j] = 3.e-4; Gal[p].WBlHo[indx_wBH].e_BH_vel[j] = 1e-3;}
        Gal[p].WBlHo[indx_wBH].e_v_sigma_now = -1.;
        //Gal[p].WBlHo[indx_wBH].to_be_deleted = -1.0;  
	Gal[p].WBlHo[indx_wBH].time_since_recoil = 0.0;
        Gal[p].WBlHo[indx_wBH].v_kick = -1.0; Gal[p].WBlHo[indx_wBH].z_ejected_WPHs = -1.0;
        Gal[p].WBlHo[indx_wBH].disrupted_galaxy = -1.0;  Gal[p].WBlHo[indx_wBH].v_escape = -1.0;
        Gal[p].WBlHo[indx_wBH].inherit = -1.0;
    }
}




void exchange_Wandering_BHs(int gal_targ, int wbh_targ, int gal_obj, int wbh_obj) 
{
 
  // This functions echange the wandering black holes of the galaxy: gat_targ to the galaxy: gal_obj targ 

  int i;
  float New_cartesianPos[3], New_cartesianVel[3];

  // First case is basically the wBH comes form a galaxy that had an ejected wBH
  if( fabs(Gal[gal_obj].WBlHo[wbh_obj].e_BH_vel[0]) > 1e-2 && fabs(Gal[gal_obj].WBlHo[wbh_obj].e_BH_vel[1])>1e-4  
   && fabs(Gal[gal_obj].WBlHo[wbh_obj].e_BH_vel[1])<1e-2 && fabs(Gal[gal_obj].WBlHo[wbh_obj].e_BH_vel[2])>1e-4 && fabs(Gal[gal_obj].WBlHo[wbh_obj].e_BH_vel[2])<1e-2)
              {compute_offset_EJECTED_wBH_respect_central_gal(gal_targ,gal_obj,wbh_obj,New_cartesianPos,New_cartesianVel);}
  // Second case is basically the WBH comes form the deposit of a central BH in the DM after a disruption
  // Third case is when the wBH is already in the DM halo and its central galaxy became a type 2
  else
    for(i=0;i<3;i++) { New_cartesianPos[i] = wrap((Gal[gal_obj].Pos[i] + Gal[gal_obj].WBlHo[wbh_obj].e_BH_pos[i]) - Gal[gal_targ].Pos[i],BoxSize);
                       New_cartesianVel[i] = (Gal[gal_obj].Vel[i] + Gal[gal_obj].WBlHo[wbh_obj].e_BH_vel[i]) - Gal[gal_targ].Vel[i]; }

  Gal[gal_targ].WBlHo[wbh_targ].e_BlackHoleMass = Gal[gal_obj].WBlHo[wbh_obj].e_BlackHoleMass;
  Gal[gal_targ].WBlHo[wbh_targ].e_Mass_to_acc = Gal[gal_obj].WBlHo[wbh_obj].e_Mass_to_acc;
  Gal[gal_targ].WBlHo[wbh_targ].e_spin = Gal[gal_obj].WBlHo[wbh_obj].e_spin;
  Gal[gal_targ].WBlHo[wbh_targ].e_eta_BH = Gal[gal_obj].WBlHo[wbh_obj].e_eta_BH;
  Gal[gal_targ].WBlHo[wbh_targ].e_fEDD = Gal[gal_obj].WBlHo[wbh_obj].e_fEDD;
  Gal[gal_targ].WBlHo[wbh_targ].e_M_Res_0 = Gal[gal_obj].WBlHo[wbh_obj].e_M_Res_0;
  Gal[gal_targ].WBlHo[wbh_targ].e_MBH_0 = Gal[gal_obj].WBlHo[wbh_obj].e_MBH_0;
  Gal[gal_targ].WBlHo[wbh_targ].e_MBH_peak = Gal[gal_obj].WBlHo[wbh_obj].e_MBH_peak;
  Gal[gal_targ].WBlHo[wbh_targ].e_MBH_second_phase = Gal[gal_obj].WBlHo[wbh_obj].e_MBH_second_phase;
  Gal[gal_targ].WBlHo[wbh_targ].e_t_acc_quiescent_phase = Gal[gal_obj].WBlHo[wbh_obj].e_t_acc_quiescent_phase;
  for(i=0;i<3;i++) { Gal[gal_targ].WBlHo[wbh_targ].e_BH_pos[i] = New_cartesianPos[i]; Gal[gal_targ].WBlHo[wbh_targ].e_BH_vel[i] = New_cartesianVel[i];}
  //Gal[gal_targ].WBlHo[wbh_targ].to_be_deleted = Gal[gal_obj].WBlHo[wbh_obj].to_be_deleted;
  Gal[gal_targ].WBlHo[wbh_targ].time_since_recoil = Gal[gal_obj].WBlHo[wbh_obj].time_since_recoil;
  Gal[gal_targ].WBlHo[wbh_targ].e_v_sigma_now = Gal[gal_obj].WBlHo[wbh_obj].e_v_sigma_now;
  Gal[gal_targ].WBlHo[wbh_targ].v_kick = Gal[gal_obj].WBlHo[wbh_obj].v_kick;
  Gal[gal_targ].WBlHo[wbh_targ].v_escape = Gal[gal_obj].WBlHo[wbh_obj].v_escape;
  Gal[gal_targ].WBlHo[wbh_targ].z_ejected_WPHs = Gal[gal_obj].WBlHo[wbh_obj].z_ejected_WPHs;
  Gal[gal_targ].WBlHo[wbh_targ].disrupted_galaxy = Gal[gal_obj].WBlHo[wbh_obj].disrupted_galaxy;
  Gal[gal_targ].WBlHo[wbh_targ].inherit = +1;// I put this by hand because by seting equal to Gal[gal_obj].WBlHo[wbh_obj].inherit it might have problems
}


void obtain_index_of_wBHs_ordered_by_mass_and_pos(int p)
{
   // This function basically obtains the order of the wBHs of a galaxy ordered by mass
   // In case two BHs have the same mass, we always put before the one which is closer to the central galaxy which is orbiting
   
   int i, x, temp_indx, counter, index[WBHS];
   float temp, dummy_WBHs[WBHS], dummy_pos_WBHs[WBHS];
   float dummy_posX_WBHs[WBHS], dummy_posY_WBHs[WBHS], dummy_posZ_WBHs[WBHS];
   float r1, r2;
   counter = 0;
   for(i=0;i<WBHS;i++) index[i] = i;
   //#ifdef CATALOGUESEEDING
   //for(i=0;i<WBHS;i++) Gal[p].indx_Order_by_most_massive_WBHs[i] = -1;
   //#endif
   for(i=0;i<WBHS;i++) {dummy_WBHs[i] = Gal[p].WBlHo[i].e_BlackHoleMass;}//; if(Gal[p].WBlHo[i].e_BlackHoleMass>1e-9) counter++;}
   for(i=0;i<WBHS;i++) {dummy_posX_WBHs[i] = Gal[p].WBlHo[i].e_BH_pos[0];
                        dummy_posY_WBHs[i] = Gal[p].WBlHo[i].e_BH_pos[1];
                        dummy_posZ_WBHs[i] = Gal[p].WBlHo[i].e_BH_pos[2];}

   for(x=0;x<WBHS;x++)
    {
      for(i=0;i<WBHS-x-1;i++)
       {
         if(dummy_WBHs[i] < dummy_WBHs[i+1])
          {temp = dummy_WBHs[i+1]; dummy_WBHs[i+1] = dummy_WBHs[i]; dummy_WBHs[i] = temp;
           temp = dummy_posX_WBHs[i+1]; dummy_posX_WBHs[i+1] = dummy_posX_WBHs[i]; dummy_posX_WBHs[i] = temp;
           temp = dummy_posY_WBHs[i+1]; dummy_posY_WBHs[i+1] = dummy_posY_WBHs[i]; dummy_posY_WBHs[i] = temp;
           temp = dummy_posZ_WBHs[i+1]; dummy_posZ_WBHs[i+1] = dummy_posZ_WBHs[i]; dummy_posZ_WBHs[i] = temp;
           temp_indx = index[i+1]; index[i+1] = index[i]; index[i] = temp_indx;}
         else
          {
           if(dummy_WBHs[i] == dummy_WBHs[i+1])
                {
                   r1 = pow(pow(dummy_posX_WBHs[i],2.0) + pow(dummy_posY_WBHs[i],2.0) + pow(dummy_posZ_WBHs[i],2.0),0.5);
                   r2 = pow(pow(dummy_posX_WBHs[i+1],2.0) + pow(dummy_posY_WBHs[i+1],2.0) + pow(dummy_posZ_WBHs[i+1],2.0),0.5);
                   if(r1 < r2)
                     {temp = dummy_WBHs[i+1]; dummy_WBHs[i+1] = dummy_WBHs[i]; dummy_WBHs[i] = temp;
                      temp = dummy_posX_WBHs[i+1]; dummy_posX_WBHs[i+1] = dummy_posX_WBHs[i]; dummy_posX_WBHs[i] = temp;
                      temp = dummy_posY_WBHs[i+1]; dummy_posY_WBHs[i+1] = dummy_posY_WBHs[i]; dummy_posY_WBHs[i] = temp;
                      temp = dummy_posZ_WBHs[i+1]; dummy_posZ_WBHs[i+1] = dummy_posZ_WBHs[i]; dummy_posZ_WBHs[i] = temp;
                      temp_indx = index[i+1]; index[i+1] = index[i]; index[i] = temp_indx;}
                }// end if
          }//end else
       }
    }


  // TODO UNCOMMENT THIS WHEN I INCLUDE THE SEEDING_MR
  //#ifdef CATALOGUESEEDING
  // for(i=0;i<WBHS;i++) Gal[p].indx_Order_by_most_massive_WBHs[i] = index[i];
  //#endif

}


void wandering_BH_Growth_phase(int gal_num, int filenr)
{
 // This function attempts to deal with the growth of wandering black holes
 // It is never used in the code given that we are not sure yet how is the growth of black holes in a wandering phase
 // Therefore, this code has been never tested
 
 int i;
 float dt, fEdd, a, eta, epsilon, MBH_i, MBH_f, MBH_teo, M_dot;
 dt = (NumToTime(Gal[gal_num].SnapNum) - NumToTime(Gal[gal_num].SnapNum + 1))/STEPS;

 for(i=0; i<WBHS; i++)
   {
      if(Gal[gal_num].WBlHo[i].e_BlackHoleMass>1e-9 && Gal[gal_num].WBlHo[i].e_Mass_to_acc>1e-13)
        {

          // Setting the initial conditions
          fEdd = 1.0;
          a = Gal[gal_num].WBlHo[i].e_spin;
          eta = Gal[gal_num].WBlHo[i].e_eta_BH;
          epsilon =  Gal[gal_num].WBlHo[i].e_eta_BH;
          MBH_i = Gal[gal_num].WBlHo[i].e_BlackHoleMass;
          MBH_f = MBH_i + (1.-eta)*Gal[gal_num].WBlHo[i].e_Mass_to_acc; // Max BH mass
          // end initial conditions

          MBH_teo = MBH_i + ( (1.-eta) * MBH_i * fEdd * dt / (epsilon * T_EDD) );
          M_dot = (MBH_teo - MBH_i)/(dt*(1.-eta)); // Mass[internal] / seg[int]
          Gal[gal_num].WBlHo[i].e_BlackHoleMass = MBH_teo;
          Gal[gal_num].WBlHo[i].e_Mass_to_acc -= M_dot*dt;

          if(MBH_teo>=MBH_f)
                {
                      Gal[gal_num].WBlHo[i].e_BlackHoleMass = MBH_f;
                      Gal[gal_num].WBlHo[i].e_Mass_to_acc = 0.;
                      M_dot = 0.0;
                      fEdd = 0.0;
                }

          if(Gal[gal_num].WBlHo[i].e_Mass_to_acc<0) {printf("Error : Mass_to_acc_total in model_merger(): black hole growth is negative!"); exit(0);}
          Gal[gal_num].WBlHo[i].e_fEDD = fEdd;
          Gal[gal_num].WBlHo[i].e_QuasarLum = fEdd * L_EDDINGTON_ERG * (MBH_i / Hubble_h);  //1e40;
	  #ifdef SPINCOMPUTATION
         // The spin only changes in the accretion phase
         if(fEdd > 0.) Self_Gravity_Accretion_Sesana_Model(gal_num,MBH_i,M_dot,fEdd,&a,&eta,epsilon,Gal[gal_num].WBlHo[i].e_v_sigma_now,dt,filenr);
         if(fEdd > 0.) if(a>=1. || a<-1. || isnan(a) ==1) {printf("Error in wandering_BH_Growth_phase(): During the gas accretion, the spin parameter is bigger than one or nan. Please check in model_mergers.c (David functions) , spin = %f, Gal[merger_centralgal].e_BlackHoleMass[i] = %E, MBH_i = %E,m_dot = %E snp = %d, BH_ejected index = %d\n",a,Gal[gal_num].WBlHo[i].e_BlackHoleMass,MBH_i,M_dot,Gal[gal_num].SnapNum,i); exit(0);}
          #else
          epsilon = 0.1; eta = 0.1; a = 0.0;
          #endif
          Gal[gal_num].WBlHo[i].e_spin = a;
          Gal[gal_num].WBlHo[i].e_eta_BH = eta;
          MBH_i = Gal[gal_num].WBlHo[i].e_BlackHoleMass;
        }// if there is a BH there

    else {Gal[gal_num].WBlHo[i].e_Mass_to_acc = 0.0; Gal[gal_num].WBlHo[i].e_QuasarLum = 0.0; Gal[gal_num].WBlHo[i].e_fEDD = 0.0; Gal[gal_num].WBlHo[i].e_M_Res_0 = 1e-15; Gal[gal_num].WBlHo[i].e_MBH_0=1.e-5;}

  }// loop in the ejected BH slots 
}



float compute_escape_velocity_halo(int p, int n_step)
{
  // This function computes the escape velocity of a given halo of Mass, Mvir and redshift zz
  float Vvir, concentration, v_escape_halo;
  concentration = compute_concentration_DM_halo(p, n_step);
  Vvir = pow(G_GRAVITY * Gal[p].Mvir / Gal[p].Rvir,1./2.); // Units [km/s]
  v_escape_halo = pow(2.0*Vvir*Vvir*concentration/NFW_g_function(concentration),1./2.); // km/s
  if(v_escape_halo<0. || isnan(v_escape_halo)==1) {printf("Error in compute_escape_velocity_halo( ): v_escape_halo<0 or isnan(v_escape_halo)==1 (v_escape_halo = %f) in Snapnum = %d, concentration = %f\n",v_escape_halo,Gal[p].SnapNum,concentration); exit(0);}
  return v_escape_halo;

}


float compute_escape_velocity(int p, float position)
{
   // This function computes the escape velocity taking into account the whole galaxy. 
   // Notice that we do not divide between pseudo-bulge or classical one. We  assume the Hernquist profile.
   // Vescape = pow(-2*Potential,0.5)


   float V_escape,V_escape_sdisk,V_escape_gdisk,V_escape_bulge,a;
   a = 0.;

   if( (Gal[p].BulgeMass<1e-9 && Gal[p].BulgeMass>1e-9) || (Gal[p].BulgeMass>1e-9 && Gal[p].BulgeMass<=1e-9)  )
        {printf("Something weird in the bulge when I compute the escape velocity\n"); exit(0);}

   // Bulge 
   if (Gal[p].BulgeMass<1.e-9 || Gal[p].BulgeSize<1.e-9) V_escape_bulge = 0.0;
   else{ a = Gal[p].BulgeSize / (1.0+pow(2.0,0.5)); V_escape_bulge = (2.0 * G_GRAVITY * Gal[p].BulgeMass) / (fabs(position) + a);}
  
   // Stellar Disk
   if(Gal[p].DiskMass<1e-9) V_escape_sdisk = 0.0;
   else V_escape_sdisk = (3.36 * G_GRAVITY * Gal[p].DiskMass )/(1.68*Gal[p].StellarDiskRadius/3.);

   // Gas Disk
   if(Gal[p].ColdGas<1e-9) V_escape_gdisk = 0.0;
   else V_escape_gdisk = (3.36* G_GRAVITY * Gal[p].ColdGas )/(1.68*Gal[p].GasDiskRadius/3.);

   // Total escape velocity
   V_escape = pow(V_escape_bulge + V_escape_sdisk + V_escape_gdisk,0.5);

   // Some checks
   if(isinf(V_escape) == 1 || isnan(V_escape) == 1) {printf("V_escape is infinity/nan: %f, V_escape_bulge = %f, V_escape_sdisk = %f. V_escape_gdisk = %f, abulge = %f, BulgeSize = %E, Mbulge = %E\n",V_escape,V_escape_bulge,V_escape_sdisk,V_escape_gdisk,a,Gal[p].BulgeSize,Gal[p].BulgeMass); exit(0);}
   return V_escape;
}


float kick_velocity(float a1, float a2, float beta, float gamma, float qq, float mu, int snapnum, int filenr)
{
  // This function computes the kick velocity of a black hole according to Lousto et al. (2012)

  double zeta, vm, v_perp, v_par, H, Am, Bm, V11, VA, VB, VC, cos_phi_delta_minus_phi_1;
  double vx, vy, vz;
  double S_par, a1_par, a2_par, a1_perp, a2_perp;
  float V_rec;
  zeta = 145.0 * M_PI / 180.0; // rads
  Am = 1.2*1e4;
  Bm = -0.93;
  H = 6.9*1e3;
  V11 = 3677.76; // km/s
  VA = 2481.21; //  km/h
  VB = 1792.45; //  km/h
  VC = 1506.52; //  km/h
  //beta should be in radians
  //gamma should be in radians
  //a1(beta) according to the BH_BH_spin_mass function
  //a2(gamma) according to the BH_BH_spin_mass funtion
  a1_par = fabs(a1)*cos(beta);
  a2_par = fabs(a2)*cos(gamma);
  a1_perp = fabs(a1)*sin(beta);
  a2_perp = fabs(a2)*sin(gamma);
  S_par = 2*(a2_par+(qq*qq*a1_par))/pow(1.0+qq,2.0);
  cos_phi_delta_minus_phi_1 =  (2.*M_PI) * ((double)rand()/((double)RAND_MAX+1)); // random between 0 and 2pi
  cos_phi_delta_minus_phi_1 = cos(cos_phi_delta_minus_phi_1);
  vm = Am*mu*mu*(1.0-qq)*(1.0+(Bm*mu))/(1.0+qq);
  v_perp = H*mu*mu*(a2_par - (qq*a1_par))/(1.0+qq);
  v_par = (16.0*mu*mu/(1.0+qq)) * ( V11 + (VA*S_par) + (VB*S_par*S_par) + (VC*S_par*S_par*S_par)) * (a2_perp - (qq*a1_perp)) * cos_phi_delta_minus_phi_1;

  vx = vm + (v_perp*cos(zeta));
  vy = v_perp*sin(zeta);
  vz = v_par;
  V_rec = pow(pow(vx,2.0) + pow(vy,2.0) + pow(vz,2.0),0.5); // km/s

  if(isnan(V_rec) == 1 || isinf(V_rec) == 1 || V_rec<=0.){printf("Error in kick_velocity() the value of V_rec = %f has a weird value!: spin1 = %f, spin2 = %f, qq = %f, snapnum = %d, filenr = %d\n",V_rec,a1,a2,qq, snapnum,filenr);exit(0);}
  return V_rec;

}





float mass_of_gas_retained_by_ejected_BH(int p, int nstep, float Vkick)
{
    // This function computes the amount of mass reteined by the black hole after an ejection from its host galaxy
    // The equations have been extracted from Loeb (2007)
    float fedd, Vmin, alpha, Mdisk;
    fedd = 1.0; // We assume that we start by Eddington limited phase
    alpha = 0.1;
    Vmin = 550 * pow(alpha/0.1,-0.28) * pow(Gal[p].BlHo[0].eta_BH[nstep],-0.21) * pow(fedd,0.21) * pow(Gal[p].BlHo[0].BlackHoleMass[nstep]*1.e3,0.43); // Units [km/s]
    if(Vkick>=Vmin) Mdisk = (1.9e-4) * pow(alpha/0.1,-0.8) * pow(Gal[p].BlHo[0].eta_BH[nstep],-0.6) * pow(fedd,0.6) * pow(Gal[p].BlHo[0].BlackHoleMass[nstep]*1.e3,2.2) * pow(Vkick/1e3,-2.8); // Units = [10^10 M_sun]
    else Mdisk = Gal[p].BlHo[0].BlackHoleMass[nstep]; // Units = [10^10/h M_sun]
    return Mdisk;
}




int check_recoil_velocity_BH(int p, int nstep, float V_recoil, float *Mass_can_retain, int *index_BH_ejected)
{
  // This function determines if the recoiled black hole is ejected from the galaxy (or even halo) according to 
  // the escape velocity of the galaxy and the modulus of the recoil velocity
  // Notation +1 is ejected; -1 is not ejected ; +2 is ejected form the halo


  float V_escape,V_escape_halo;
  V_escape = compute_escape_velocity(p,0.0);
  V_escape_halo = compute_escape_velocity_halo(p,nstep);

  if(((1.0*V_recoil)/(1.0*V_escape)) > 1.) 
      {
        Gal[p].BlackHoleMass_ejected = Gal[p].BlHo[0].BlackHoleMass[nstep];
        Gal[p].Snp_ejected = Gal[p].SnapNum;

        // FIRST CASE: The BH is ejected from the halo 
        if(((1.0*V_recoil)/(1.0*V_escape_halo))>1.)
         {
           Gal[p].ColdGas += Gal[p].BlHo[0].Mass_to_acc[nstep];
           restart_the_BH(p,0,nstep);
           return +2;
         }

        #ifdef FOLLOWWANDERINGORBITS
        float Velocity[3] = {V_recoil,1e-3,1e-3};
        *Mass_can_retain = mass_of_gas_retained_by_ejected_BH(p,nstep,V_recoil); // The disk that a ejected BH can retain accroding to Loeb et al 2007 and Guedes et al. 2009 
	// I only save the BH in a WBH slot if and only if the BHMass ejected is above the threshold WBH_MASS_TH
	if(Gal[p].BlackHoleMass_ejected > (WBH_MASS_TH*1e-10/Hubble_h) )  *index_BH_ejected = add_BH_in_halo(p,nstep,Velocity,-1); // Here I transfere the properties of the central BH to one of the N-allowed (index_BH_ejected) wandering BHs
        // Notice that the code how it is written avoids the growth of the wBHs because we do not know yet how it grows in the wandering phase
        *Mass_can_retain = 1e-15;
        Gal[p].ColdGas += Gal[p].WBlHo[*index_BH_ejected].e_Mass_to_acc;
        Gal[p].WBlHo[*index_BH_ejected].e_Mass_to_acc = *Mass_can_retain;
        Gal[p].ColdGas -= *Mass_can_retain;
        #else
        Gal[p].ColdGas += Gal[p].BlHo[0].Mass_to_acc[nstep];
        #endif
        restart_the_BH(p,0,nstep);
        return +1; // This mean that is ejected
       }//end case of ejection

  else { *index_BH_ejected = -1; *Mass_can_retain = -1.; return -1;} // This means that is not ejected

}


float add_BH_in_halo(int p, int nstep, float Vel[3], float form_disruption)
{
 
  // This function adds a black hole in a given halo (of the galaxy with index p)
  // The value of Vel is the initial black hole velocity when it is added in the halo
  // form_disruption means if the black hole comes form a galaxy disruption (orphan wandering black hole, +1) 
  // of if it comes from an ejection (-1)  --> This difference is important because we track in a different way the 
  // orbits of orphan black holes (3D orbit) and ejected black holes (radial orbit)
  // The return is the position of the black hole in the wBH array

  int i,find, pos_BH_ejected_min;
  float BH_ejected_min;
  BH_ejected_min = 1e8;
  find = -1;
  pos_BH_ejected_min = -1.0;

  for(i=0; i<WBHS; i++) // loop for all the wandering black holes
    {
	// This loop is done to find a free slot in the array of wandering black hole
        if(Gal[p].WBlHo[i].e_BlackHoleMass<1e-9) // Find the free slot!
           {
                find = +1;
                pos_BH_ejected_min = i;
                Gal[p].WBlHo[i].e_BlackHoleMass = Gal[p].BlHo[0].BlackHoleMass[nstep];
                Gal[p].WBlHo[i].e_Mass_to_acc = Gal[p].BlHo[0].Mass_to_acc[nstep];
                Gal[p].WBlHo[i].e_spin = Gal[p].BlHo[0].spin[nstep];
                Gal[p].WBlHo[i].e_eta_BH = Gal[p].BlHo[0].eta_BH[nstep];
                Gal[p].WBlHo[i].e_BH_vel[0] = Vel[0];
                Gal[p].WBlHo[i].e_BH_vel[1] = Vel[1];
                Gal[p].WBlHo[i].e_BH_vel[2] = Vel[2];
                // This three following lines are a bit redundant but just in case...
                Gal[p].WBlHo[pos_BH_ejected_min].e_BH_pos[0] = 3.e-4;
                Gal[p].WBlHo[pos_BH_ejected_min].e_BH_pos[1] = 3.e-4;
                Gal[p].WBlHo[pos_BH_ejected_min].e_BH_pos[2] = 3.e-4;
                Gal[p].WBlHo[i].e_MBH_0 = Gal[p].BlHo[0].MBH_0;
                Gal[p].WBlHo[i].e_v_sigma_now = Gal[p].v_sigma_now;
                Gal[p].WBlHo[i].time_since_recoil = 0.0;
                Gal[p].WBlHo[i].v_kick = pow( pow(Vel[0],2) + pow(Vel[1],2) + pow(Vel[2],2), 0.5);
                Gal[p].WBlHo[i].v_escape = compute_escape_velocity(p,0.0);
                Gal[p].WBlHo[i].z_ejected_WPHs =  redshift_at_t(ZZ[Gal[p].SnapNum],nstep,(NumToTime(Gal[p].SnapNum)-NumToTime(Gal[p].SnapNum+1) ) / STEPS);
                break;
           }
    }
  if(find<0) 
   {
	// Is all the slots of wBHs in the galaxy are  booked: If the new ejected one is bigger than other wandering BH I keep this new one a I remove the smallest!
        for(i=0; i<WBHS; i++) // loop for all the wandering black holes
          if(Gal[p].WBlHo[i].e_BlackHoleMass<BH_ejected_min)
           {BH_ejected_min = Gal[p].WBlHo[i].e_BlackHoleMass; pos_BH_ejected_min = i;}// This is the smaller BH

        if(BH_ejected_min<Gal[p].BlHo[0].BlackHoleMass[nstep])
          {
             restart_the_wandering_BH(p,pos_BH_ejected_min);
             Gal[p].WBlHo[pos_BH_ejected_min].e_BlackHoleMass = Gal[p].BlHo[0].BlackHoleMass[nstep];
             Gal[p].WBlHo[pos_BH_ejected_min].e_Mass_to_acc = Gal[p].BlHo[0].Mass_to_acc[nstep];
             Gal[p].WBlHo[pos_BH_ejected_min].e_spin = Gal[p].BlHo[0].spin[nstep];
             Gal[p].WBlHo[pos_BH_ejected_min].e_eta_BH = Gal[p].BlHo[0].eta_BH[nstep];
             Gal[p].WBlHo[pos_BH_ejected_min].e_BH_vel[0] = Vel[0];
             Gal[p].WBlHo[pos_BH_ejected_min].e_BH_vel[1] = Vel[1];
             Gal[p].WBlHo[pos_BH_ejected_min].e_BH_vel[2] = Vel[2];
             Gal[p].WBlHo[pos_BH_ejected_min].e_BH_pos[0] = 3.e-4;
             Gal[p].WBlHo[pos_BH_ejected_min].e_BH_pos[1] = 3.e-4;
             Gal[p].WBlHo[pos_BH_ejected_min].e_BH_pos[2] = 3.e-4;
             Gal[p].WBlHo[pos_BH_ejected_min].e_MBH_0 = Gal[p].BlHo[0].BlackHoleMass[nstep];
             Gal[p].WBlHo[pos_BH_ejected_min].e_v_sigma_now = Gal[p].v_sigma_now;
             Gal[p].WBlHo[pos_BH_ejected_min].time_since_recoil = 0.0;
             Gal[p].WBlHo[pos_BH_ejected_min].v_kick = pow( pow(Vel[0],2) + pow(Vel[1],2) + pow(Vel[2],2) ,0.5);
             Gal[p].WBlHo[pos_BH_ejected_min].z_ejected_WPHs =  redshift_at_t(ZZ[Gal[p].SnapNum],nstep,(NumToTime(Gal[p].SnapNum)-NumToTime(Gal[p].SnapNum+1) ) / STEPS);
          }
   }
   //TODO --> NOTICE THIS: 
   Gal[p].WBlHo[pos_BH_ejected_min].inherit = -1.0;
   Gal[p].ColdGas += Gal[p].WBlHo[pos_BH_ejected_min].e_Mass_to_acc; // all the mass of the black hole (in case it has) is deposited again in the galaxy
   Gal[p].WBlHo[pos_BH_ejected_min].e_Mass_to_acc = 1e-15;
   Gal[p].WBlHo[pos_BH_ejected_min].disrupted_galaxy = form_disruption;
 return pos_BH_ejected_min;
}




void add_in_gal_the_wondering_BH(int p, int indx_wondering, int nstep, int filenr, float dt)
{

  int indx_BH, j;
  indx_BH = -1;

  if( ( Gal[p].BlHo[1].BlackHoleMass[nstep] >= Gal[p].WBlHo[indx_wondering].e_BlackHoleMass ) && ( Gal[p].BlHo[0].BlackHoleMass[nstep] >= Gal[p].WBlHo[indx_wondering].e_BlackHoleMass ) )
     {restart_the_wandering_BH(p,indx_wondering);}
  else
   {
     if( (Gal[p].BlHo[0].BlackHoleMass[nstep] < 1e-14) && (Gal[p].BlHo[1].BlackHoleMass[nstep] <1e-14) )
      {
        // Basically there is no BH in the galaxy. So the new BH is the one that is going to be the main one
        for(j=nstep; j<STEPS; j++)
         {
           Gal[p].BlHo[0].BlackHoleMass[j] = Gal[p].WBlHo[indx_wondering].e_BlackHoleMass;
           Gal[p].BlHo[0].Mass_to_acc[j] += Gal[p].WBlHo[indx_wondering].e_Mass_to_acc;
           Gal[p].BlHo[0].QuasarLum[j] = 0.0;// Gal[p].e_QuasarLum[i][j];       
           Gal[p].BlHo[0].spin[j] = Gal[p].WBlHo[indx_wondering].e_spin;
           Gal[p].BlHo[0].eta_BH[j] = Gal[p].WBlHo[indx_wondering].e_eta_BH;
           Gal[p].BlHo[0].fEDD[j] = Gal[p].WBlHo[indx_wondering].e_fEDD;
          }
           Gal[p].BlHo[0].M_res_0 = 0.; Gal[p].BlHo[0].MBH_0 = 0.; Gal[p].BlHo[0].MBH_peak = 0.; Gal[p].BlHo[0].MBH_second_phase = 0.;
           Gal[p].BlHo[0].M_dot_acc = 0.0; Gal[p].BlHo[0].t_acc_quiescent_phase = -1.;
           //Gal[p].BlHo[0].M_acc_Mergers = 0.; Gal[p].BlHo[0].M_acc_DI_sec = 0.; Gal[p].BlHo[0].M_acc_DI_mmi = 0.;
           Gal[p].BlHo[0].tBH_merger_clock_dt = -1.0; Gal[p].BlHo[0].tBH_merger = -1.0; Gal[p].fgas = -1.0; Gal[p].z_binary_form = -1.0;
           Gal[p].t_binary_form = -1.0;
           Gal[p].BlHo[0].t_spent_WPHs = Gal[p].WBlHo[indx_wondering].time_since_recoil;  Gal[p].BlHo[0].z_ejected_WPHs = Gal[p].WBlHo[indx_wondering].z_ejected_WPHs;
           Gal[p].BlHo[0].MBH_WPH = Gal[p].WBlHo[indx_wondering].e_BlackHoleMass;

           // Re-start the slot of ejected BH   
           restart_the_wandering_BH(p,indx_wondering);
       }//end if 
     else
      {
           if(Gal[p].WBlHo[indx_wondering].e_BlackHoleMass > Gal[p].BlHo[1].BlackHoleMass[nstep])
                {indx_BH = 1;}
           if(Gal[p].WBlHo[indx_wondering].e_BlackHoleMass > Gal[p].BlHo[0].BlackHoleMass[nstep])
                {indx_BH  = 0; Gal[p].BlHo[1] = Gal[p].BlHo[0]; Gal[p].BlHo[1].tBH_merger_clock_dt = -1.0; Gal[p].BlHo[1].tBH_merger = -1.0;}
           if(indx_BH<0) {printf("Error in add_in_gal_the_wondering_BH(). Filenr = %d, WNB_index = %d, The index_BH is negative (%d), BH0 = %E, BH1 = %E, WBH = %E, nstep = %d\n",filenr,indx_wondering,indx_BH,Gal[p].BlHo[0].BlackHoleMass[nstep],Gal[p].BlHo[1].BlackHoleMass[nstep],Gal[p].WBlHo[indx_wondering].e_BlackHoleMass,nstep); exit(0);}



           #ifdef TRIPLETS_EJECTED
                  if(Gal[p].BlHo[indx_BH].BlackHoleMass[nstep]>1e-9)
                  {
                    FILE *fileout;
                    char file_BH[5000];
                    float zzz,dt_dummy;
                    dt_dummy = (NumToTime(Gal[p].SnapNum) - NumToTime(Gal[p].SnapNum + 1))/STEPS;
                    zzz = redshift_at_t(ZZ[Gal[p].SnapNum],nstep,dt_dummy);
                    sprintf(file_BH, "/global/users/dizquierdo/thesis/LGalaxies/Recoil_BH/Hybrid_model/Triplets_ejected/Snap_%d", filenr);
                    //printf("zzz = %f ZZ[Gal[p].SnapNum] = %f, Gal[p].SnapNum=%d, dt = %E, nstep = %d, dt_dummy = %E\n",zzz,ZZ[Gal[p].SnapNum],Gal[p].SnapNum,dt,nstep,dt_dummy);
                    fileout = fopen(file_BH,"ab");
                    float arr[8] = {Gal[p].BulgeMass+Gal[p].DiskMass,
                                Gal[p].BulgeMass,
                                Gal[p].DiskMass,
                                Gal[p].Mvir,
                                Gal[p].Bulge_di/Gal[p].BulgeMass,
                                Gal[p].BlHo[indx_BH].BlackHoleMass[nstep],
                                zzz,-1.0};
                    fwrite(arr,sizeof(arr),1,fileout);
                    fclose(fileout);
                  }
           #endif

           for(j=nstep; j<STEPS; j++)
              {
                Gal[p].BlHo[indx_BH].BlackHoleMass[j] = Gal[p].WBlHo[indx_wondering].e_BlackHoleMass;
                Gal[p].BlHo[indx_BH].Mass_to_acc[j] += Gal[p].WBlHo[indx_wondering].e_Mass_to_acc;
                Gal[p].BlHo[indx_BH].QuasarLum[j] = 0.0;// Gal[p].e_QuasarLum[i][j];       
                Gal[p].BlHo[indx_BH].spin[j] = Gal[p].WBlHo[indx_wondering].e_spin;
                Gal[p].BlHo[indx_BH].eta_BH[j] = Gal[p].WBlHo[indx_wondering].e_eta_BH;
                Gal[p].BlHo[indx_BH].fEDD[j] = Gal[p].WBlHo[indx_wondering].e_fEDD;
              }
                Gal[p].BlHo[indx_BH].M_res_0 = 0.; Gal[p].BlHo[indx_BH].MBH_0 = 0.; Gal[p].BlHo[indx_BH].MBH_peak = 0.; Gal[p].BlHo[indx_BH].MBH_second_phase = 0.;
                Gal[p].BlHo[indx_BH].M_dot_acc = 0.0; Gal[p].BlHo[indx_BH].t_acc_quiescent_phase = -1.;
                //Gal[p].BlHo[indx_BH].M_acc_Mergers = 0.; Gal[p].BlHo[indx_BH].M_acc_DI_sec = 0.; Gal[p].BlHo[indx_BH].M_acc_DI_mmi = 0.;
                Gal[p].BlHo[indx_BH].tBH_merger_clock_dt = -1.0; Gal[p].BlHo[indx_BH].tBH_merger = -1.0; Gal[p].fgas = -1.0; Gal[p].z_binary_form = -1.0;
                Gal[p].t_binary_form = -1.0;
                Gal[p].BlHo[indx_BH].t_spent_WPHs = Gal[p].WBlHo[indx_wondering].time_since_recoil;  Gal[p].BlHo[indx_BH].z_ejected_WPHs = Gal[p].WBlHo[indx_wondering].z_ejected_WPHs;
                Gal[p].BlHo[indx_BH].MBH_WPH = Gal[p].WBlHo[indx_wondering].e_BlackHoleMass;


           if(Gal[p].BlHo[0].BlackHoleMass[nstep]<Gal[p].BlHo[1].BlackHoleMass[nstep]) {printf("In reincirporation of the BH function the BH0<BH1 (file = %d)\n",filenr); exit(0);}
          // Re-start the slot of ejected BH
           restart_the_wandering_BH(p,indx_wondering);

          if(Gal[p].BlHo[0].BlackHoleMass[nstep] > 1e-14 && Gal[p].BlHo[1].BlackHoleMass[nstep] > 1e-14)
            {
                // I compute the delay in the BH - BH merger
                float dt;
                dt = (NumToTime(Gal[p].SnapNum) - NumToTime(Gal[p].SnapNum + 1))/STEPS;
                compute_BH_BH_merger_delay(p,nstep,filenr);
                //printf("I'm in the reincirporation of the wBH\n");
                if(Gal[p].BlHo[0].tBH_merger_clock_dt < dt)
                      merger_BHs(p,nstep,filenr);
            }
        } // end internal else 
   } //end external else
}




void compute_offset_EJECTED_wBH_respect_central_gal(int gal_targ, int gal_obj, int wbh_obj, float New_cartesianPos[3], float New_cartesianVel[3])
{
  // Computes the new off between wandering black hole and its new galaxy (the previous gaalxy was converted into type 2 galaxy)
  // The wandering black holes has to be type EJECTED!
  // Variables description: 
  // gal_targ is the id of the galaxy which is going to recieve the wBH (normally central)
  // gal_obj is the id of the galaxy who is going to give the wBH (normally satellite)
  // wbh_obj is the indx of the given wBH from galaxy gal_obj

  float u, v, theta, phi;
  float pos_wBH[3], vel_wBH[3];
  int j;

  //The offset that I'm going to use is not that one!

  u = ((1.0)*(((double)rand()/(double)RAND_MAX)));
  v = ((1.0)*(((double)rand()/(double)RAND_MAX)));
  phi = 2.0*M_PI*u; // rad
  theta = acos((2.0*v)-1.0); // rad
  //x = fabs(Gal[gal_obj].WBlHo[wbh_obj].e_BH_pos) * sin(theta) * cos(phi);
  //y = fabs(Gal[gal_obj].WBlHo[wbh_obj].e_BH_pos) * sin(theta) * sin(phi);
  //z = fabs(Gal[gal_obj].WBlHo[wbh_obj].e_BH_pos) * cos(theta);

  pos_wBH[0] = Gal[gal_obj].Pos[0] + (fabs(Gal[gal_obj].WBlHo[wbh_obj].e_BH_pos[0]) * sin(theta) * cos(phi));//x;//BoxSize
  pos_wBH[1] = Gal[gal_obj].Pos[1] + (fabs(Gal[gal_obj].WBlHo[wbh_obj].e_BH_pos[1]) * sin(theta) * sin(phi));//BoxSize
  pos_wBH[2] = Gal[gal_obj].Pos[2] + (fabs(Gal[gal_obj].WBlHo[wbh_obj].e_BH_pos[2]) * cos(theta));//BoxSize

  vel_wBH[0] = Gal[gal_obj].Vel[0] + (fabs(Gal[gal_obj].WBlHo[wbh_obj].e_BH_vel[0]) * sin(theta) * cos(phi));//x;//BoxSize
  vel_wBH[1] = Gal[gal_obj].Vel[1] + (fabs(Gal[gal_obj].WBlHo[wbh_obj].e_BH_vel[1]) * sin(theta) * sin(phi));//BoxSize
  vel_wBH[2] = Gal[gal_obj].Vel[2] + (fabs(Gal[gal_obj].WBlHo[wbh_obj].e_BH_vel[2]) * cos(theta));//BoxSize

  for(j=0;j<3;j++)
        {
           New_cartesianPos[j] = pos_wBH[j]- Gal[gal_targ].Pos[j];
           New_cartesianVel[j] = vel_wBH[j] - Gal[gal_targ].Vel[j];
        }
}



void join_wandering_BHs(int central, int satellite, int filenr)
{

   //  Here I move all the wandering black holes from the satellite galaxy (with id satellite)
   //  to the central one (central)

   int i, k;
   int n_satellite, n_central;

   // First I count the wandering BHs that are in the central and satellite galaxy
   n_satellite = 0;
   n_central = 0;
   for(i=0; i<WBHS; i++)
      {
        if(Gal[central].WBlHo[i].e_BlackHoleMass>1e-9) n_central++;
        if(Gal[satellite].WBlHo[i].e_BlackHoleMass>1e-9) n_satellite++;
      }

   if((n_central + n_satellite) > WBHS)
    {
      reorganize_wanderingBHs(central,satellite,filenr);
      restart_the_wandering_BH(satellite,-1);
    }
    else
     {
        if(n_satellite>0)
         {
           for(k=0; k<WBHS; k++) // loop in all the WBHs in the satellite galaxy
             if(Gal[satellite].WBlHo[k].e_BlackHoleMass>1e-9) // if in this position of the array there is a WBHs
              for(i=0; i<WBHS; i++) // loop in all the WBHs in the central galaxy
                if(Gal[central].WBlHo[i].e_BlackHoleMass<1e-9) // if in this position of the array there is a free slot WBHs
                 {
                  exchange_Wandering_BHs(central,i,satellite,k);
                  restart_the_wandering_BH(satellite,k);
                  break;
                 } // if there is a free slot in the central galaxy 
             }
        }
}




void reorganize_wanderingBHs(int central, int satellite, int filenr)
{

  // This function for a given set of wandering black holes, coming form the satellite and central galaxy
  // it re-order all of them in the [WBHS]-length array, saving only the ones that are the most massive
  // In case you have two equal mass black hole, I save always the one that is closer to the galaxy

  int i;
  int indx_smallest_central;
  int indx_bigest_satellite;
  int to_be_ended;
  float max, min;
  float r1, r2;
  to_be_ended = -1;

  while(to_be_ended<0)
   {
     max = -1.0; min = 1e3; indx_smallest_central = -1; indx_bigest_satellite = -1;
     for(i=0; i<WBHS; i++)
       {

        if(max<Gal[satellite].WBlHo[i].e_BlackHoleMass)
            {max = Gal[satellite].WBlHo[i].e_BlackHoleMass; indx_bigest_satellite = i;}
        else
          if(max == Gal[satellite].WBlHo[i].e_BlackHoleMass && indx_bigest_satellite>0)
            {
	       r1 = pow( pow(wrap(Gal[satellite].WBlHo[i].e_BH_pos[0]-Gal[central].Pos[0],BoxSize),2.0) + pow(wrap(Gal[satellite].WBlHo[i].e_BH_pos[1]-Gal[central].Pos[1],BoxSize),2.0) + pow(wrap(Gal[satellite].WBlHo[i].e_BH_pos[2]-Gal[central].Pos[2],BoxSize),2.0), 0.5);
               r2 = pow( pow(wrap(Gal[satellite].WBlHo[indx_bigest_satellite].e_BH_pos[0]-Gal[central].Pos[0],BoxSize),2.0) + pow(wrap(Gal[satellite].WBlHo[indx_bigest_satellite].e_BH_pos[1]-Gal[central].Pos[1],BoxSize),2.0) + pow(wrap(Gal[satellite].WBlHo[indx_bigest_satellite].e_BH_pos[2]-Gal[central].Pos[2],BoxSize),2.0), 0.5);

               if(r1<r2)
                        {max = Gal[satellite].WBlHo[i].e_BlackHoleMass; indx_bigest_satellite = i;}
            }
        if(min>Gal[central].WBlHo[i].e_BlackHoleMass){min = Gal[central].WBlHo[i].e_BlackHoleMass; indx_smallest_central = i;}
        else
          if(min == Gal[central].WBlHo[i].e_BlackHoleMass && indx_smallest_central>0)
            {
                r1 = pow( pow(Gal[central].WBlHo[i].e_BH_pos[0],2.0) + pow(Gal[central].WBlHo[i].e_BH_pos[1],2.0) + pow(Gal[central].WBlHo[i].e_BH_pos[2],2.0), 0.5);
                r2 = pow( pow(Gal[central].WBlHo[indx_smallest_central].e_BH_pos[0],2.0) + pow(Gal[central].WBlHo[indx_smallest_central].e_BH_pos[1],2.0) + pow(Gal[central].WBlHo[indx_smallest_central].e_BH_pos[2],2.0), 0.5);
                if(r1>r2)
                   {min = Gal[central].WBlHo[i].e_BlackHoleMass; indx_smallest_central = i;}
            }
       }

     if(indx_smallest_central < 0 || indx_bigest_satellite < 0)
      {printf("Error in  reorganize_wonderingBHs() (File = %d). indx_smallest_central<0 or indx_bigest_satellite < 0",filenr); exit(0);}
     if(max>min)
      {
        exchange_Wandering_BHs(central,indx_smallest_central,satellite,indx_bigest_satellite);
        to_be_ended = -1;
        restart_the_wandering_BH(satellite,indx_bigest_satellite);

      }
     else to_be_ended = +1;

   }

  restart_the_wandering_BH(satellite,-1);

}




void deal_with_wondering_BHs_in_satellite(int ngal, int filenr)
{
  // This function moves ALL the wandering black holes of Type 2 galaxies in the 
  // wandering balck hole array of its central galaxy (Gal[i].CentralGal)
  int i;
  float offset;
  for(i=0;i<ngal;i++)
      if(Gal[i].Type == 2)
        join_wandering_BHs(Gal[i].CentralGal,i,filenr);
}

void checks_wandering_BHs_in_type2_gal(int p, int filenr)
{ 

  // This function checks if all the wandering black holes of the Type 2 galaxies
  // have been moved in their respective central galaxy (Gal[p].CentralGal)
  int i,n_wBH_type2;
  float offset;
  n_wBH_type2 = 0;

  for(i=0;i<WBHS;i++)
      if(Gal[p].WBlHo[i].e_BlackHoleMass>1e-9)
             n_wBH_type2++;
   if(n_wBH_type2>0)
        join_wandering_BHs(Gal[p].CentralGal,p,filenr);
}



/*  This part is the part of the code in which I implement all the equations used 
 *  for tracking the orbit of a wandering black hole (see the details in Izquierdo-Villalba et al. 2020
 */


void update_position_wandering_BHs(int p, float dt, int nstep, int filenr)
{
  // This function deals with the integration of the wandering black hole orbit
  // This functions calls the Runge-Kutta numerical integration
  // In this case the time-step used to integrate the orbits are dt_Rauge_kutta set to 1 Myr of resolution

  float v_n,v_ini,p_ini,p_n,zz,dt_dummy,time,dt_Rauge_kutta, Vescape_halo;
  float r,v_tot;
  int i, k, Ejected;
  float Rgal;
  Vescape_halo = compute_escape_velocity_halo(p,nstep);
  dt_Rauge_kutta = 1.0; // [Myr]
  // Effective radius of the galaxy --> Mass weighted
  Rgal =((Gal[p].BulgeMass*Gal[p].BulgeSize)+(Gal[p].DiskMass*Gal[p].StellarDiskRadius/3.*1.68)+(Gal[p].ColdGas*Gal[p].GasDiskRadius/3.*1.68))/(Gal[p].ColdGas + Gal[p].DiskMass + Gal[p].BulgeMass);
  if(Rgal<3e-4) Rgal = 3e-4;
  dt_dummy = ( NumToTime(Gal[p].SnapNum) -  NumToTime(Gal[p].SnapNum+1) ) / STEPS;
  zz = redshift_at_t(ZZ[Gal[p].SnapNum],nstep,dt_dummy);


  // FIRST we loop in all the ejected BHs
  for(i=0;i<WBHS;i++) 
   {
     if(Gal[p].WBlHo[i].e_BlackHoleMass>1e-9)
      {
         // Now, according to the velocity of the wandering black hole, we can deduce if it is Ejected (the orbit will be radial) 
	 // or if it is orphan (the owrbit will be computed in the 3D)
         if( fabs(Gal[p].WBlHo[i].e_BH_vel[0]) > 1e-2 && fabs(Gal[p].WBlHo[i].e_BH_vel[1])>1e-4 && fabs(Gal[p].WBlHo[i].e_BH_vel[1])<1e-2 && fabs(Gal[p].WBlHo[i].e_BH_vel[2])>1e-4 && fabs(Gal[p].WBlHo[i].e_BH_vel[2])<1e-2) Ejected = +1;
         else  Ejected = -1;

          for(time = dt_Rauge_kutta; time<dt;time=time+dt_Rauge_kutta) // time should be dt/division_time [Myr]
           {
	      // We set the initial conditions of the differential equation of the orbit
              if(Ejected >0)
                 p_ini = Gal[p].WBlHo[i].e_BH_pos[0];
              else
                 p_ini = pow( pow(Gal[p].WBlHo[i].e_BH_pos[0],2) + pow(Gal[p].WBlHo[i].e_BH_pos[1],2) + pow(Gal[p].WBlHo[i].e_BH_pos[2],2),0.5);

              v_ini = pow( pow(Gal[p].WBlHo[i].e_BH_vel[0],2) + pow(Gal[p].WBlHo[i].e_BH_vel[1],2) + pow(Gal[p].WBlHo[i].e_BH_vel[2],2),0.5);
              Gal[p].WBlHo[i].time_since_recoil+= dt_Rauge_kutta;
              integrate_Runge_kutta(p,dt_Rauge_kutta,Gal[p].WBlHo[i].e_BH_pos,Gal[p].WBlHo[i].e_BH_vel,Gal[p].WBlHo[i].e_BlackHoleMass,Gal[p].Mvir,zz,nstep,Ejected,filenr);
              v_n = pow( pow(Gal[p].WBlHo[i].e_BH_vel[0],2) + pow(Gal[p].WBlHo[i].e_BH_vel[1],2) + pow(Gal[p].WBlHo[i].e_BH_vel[2],2),0.5);
              if(Ejected >0) p_n = Gal[p].WBlHo[i].e_BH_pos[0]; // The coordinate is just the radial one
              else p_n = pow( pow(Gal[p].WBlHo[i].e_BH_pos[0],2) + pow(Gal[p].WBlHo[i].e_BH_pos[1],2) + pow(Gal[p].WBlHo[i].e_BH_pos[2],2),0.5);

              //Here we do some tests, i.e we check if the integration was good or there are weird things
              if(isnan(Gal[p].WBlHo[i].e_BH_vel[0])==1 || isnan(Gal[p].WBlHo[i].e_BH_vel[1])==1 || isnan(Gal[p].WBlHo[i].e_BH_vel[1])==1 || isnan(Gal[p].WBlHo[i].e_BH_pos[0])==1 || isnan(Gal[p].WBlHo[i].e_BH_pos[1])==1 || isnan(Gal[p].WBlHo[i].e_BH_pos[2])==1) {printf("Error in update_position_wandering_BHs(): the value of v_x(%f) v_y(%f) v_z(%f) or p_x(%f) p_y(%f) p_z(%f)is nan!\n",Gal[p].WBlHo[i].e_BH_vel[0],Gal[p].WBlHo[i].e_BH_vel[1],Gal[p].WBlHo[i].e_BH_vel[2],Gal[p].WBlHo[i].e_BH_pos[0],Gal[p].WBlHo[i].e_BH_pos[1],Gal[p].WBlHo[i].e_BH_pos[2]); exit(0);}

              r = pow(pow(Gal[p].WBlHo[i].e_BH_pos[0],2.0) + pow(Gal[p].WBlHo[i].e_BH_pos[1],2.0) + pow(Gal[p].WBlHo[i].e_BH_pos[2],2.0) , 0.5);
              v_tot = pow(pow(Gal[p].WBlHo[i].e_BH_vel[0],2.0) + pow(Gal[p].WBlHo[i].e_BH_vel[1],2.0) + pow(Gal[p].WBlHo[i].e_BH_vel[2],2.0) , 0.5);

	      // If the BH is at 3*Rvir and moving away ( p_ini < p_n (end) ) from the center of the galaxy, we assume that the BH will never return to the center of the galaxy
              if( ( r  > 3.0*Gal[p].Rvir) && (fabs(p_n)>fabs(p_ini)) ) {restart_the_wandering_BH(p,i); break;}
	      // Now I check if after the the n+1 [Delta_dt = dt_Rauge_kutta] we had a re-incorporation of the BH
              if( ( r < Rgal ) && ( v_tot  <= compute_escape_velocity(p,r) ) )
               { add_in_gal_the_wondering_BH(p, i, nstep, filenr, dt); break; }
           } // end integration range kutta
      } // end if it exits a wBH
   } // end in the loop of number of wandering BHs

} 




void integrate_Runge_kutta(int p, float deltaT, float pos[3], float vel[3], float MBH, float Mvir, float redshift, int n_step, int Ejected, int filenr)
{
 
  // Here I used the expression of the Runge-Kutta expression to integrate the orbit ( i.e the function called aceleration_ejected_BH(...) ) 
 //deltaT [Myr]
 double k1n,k2n,k3n,k4n;
 double q1n,q2n,q3n,q4n;
 double c1n,c2n,c3n,c4n;
 double l1n,l2n,l3n,l4n;
 double m1n,m2n,m3n,m4n;
 double d1n,d2n,d3n,d4n;
 float r, v_tot;


  if(Ejected < 0)
   {
        r = pow(  pow(pos[0],2) + pow(pos[1],2) + pow(pos[2],2) , 0.5);
        v_tot = pow(  pow(vel[0],2) + pow(vel[1],2) + pow(vel[2],2) , 0.5);
     // X-axis
        k1n = (deltaT*SECONDS_TO_MEGAYRS)*(vel[0])*KM_TO_MPC*Hubble_h; // Mpc/h 
        l1n = deltaT*aceleration_ejected_BH(p,pos[0],r,v_tot,vel[0],Mvir,redshift,MBH,n_step,filenr)/(Hubble_h*KM_TO_MPC*SECONDS_TO_MEGAYRS); // km/s
     // Y-axis
        q1n = (deltaT*SECONDS_TO_MEGAYRS)*(vel[1])*KM_TO_MPC*Hubble_h; // Mpc/h 
        m1n = deltaT*aceleration_ejected_BH(p,pos[1],r,v_tot,vel[1],Mvir,redshift,MBH,n_step,filenr)/(Hubble_h*KM_TO_MPC*SECONDS_TO_MEGAYRS); // km/s
     // Z-axis
        c1n = (deltaT*SECONDS_TO_MEGAYRS)*(vel[2])*KM_TO_MPC*Hubble_h; // Mpc/h 
        d1n = deltaT*aceleration_ejected_BH(p,pos[1],r,v_tot,vel[2],Mvir,redshift,MBH,n_step,filenr)/(Hubble_h*KM_TO_MPC*SECONDS_TO_MEGAYRS); // km/s


        r = pow(  pow(pos[0]+0.5*k1n,2) + pow(pos[1]+0.5*q1n,2) + pow(pos[2]+0.5*c1n,2) , 0.5);
        v_tot = pow(  pow(vel[0]+(0.5*l1n),2) + pow(vel[1]+(0.5*m1n),2) + pow(vel[2]+(0.5*d1n),2) , 0.5);
     // X-axis
        k2n = (deltaT*SECONDS_TO_MEGAYRS)*(vel[0]+(0.5*l1n))*KM_TO_MPC*Hubble_h; // Mpc/h
        l2n = deltaT*aceleration_ejected_BH(p,pos[0]+0.5*k1n,r,v_tot,vel[0]+(0.5*l1n),Mvir,redshift,MBH,n_step,filenr)/(Hubble_h*KM_TO_MPC*SECONDS_TO_MEGAYRS); // km/s
     // Y-axis
        q2n = (deltaT*SECONDS_TO_MEGAYRS)*(vel[1]+(0.5*m1n))*KM_TO_MPC*Hubble_h; // Mpc/h
        m2n = deltaT*aceleration_ejected_BH(p,pos[1]+0.5*q1n,r,v_tot,vel[1]+(0.5*m1n),Mvir,redshift,MBH,n_step,filenr)/(Hubble_h*KM_TO_MPC*SECONDS_TO_MEGAYRS); // km/s
     // Z-axis
        c2n = (deltaT*SECONDS_TO_MEGAYRS)*(vel[2]+(0.5*d1n))*KM_TO_MPC*Hubble_h; // Mpc/h
        d2n = deltaT*aceleration_ejected_BH(p,pos[2]+0.5*c1n,r,v_tot,vel[2]+(0.5*d1n),Mvir,redshift,MBH,n_step,filenr)/(Hubble_h*KM_TO_MPC*SECONDS_TO_MEGAYRS); // km/s


        r = pow(  pow(pos[0]+0.5*k2n,2) + pow(pos[1]+0.5*q2n,2) + pow(pos[2]+0.5*c2n,2) , 0.5);
        v_tot = pow(  pow(vel[0]+(0.5*l2n),2) + pow(vel[1]+(0.5*m2n),2) + pow(vel[2]+(0.5*d2n),2) , 0.5);
     // X-axis
        k3n = (deltaT*SECONDS_TO_MEGAYRS)*(vel[0]+(0.5*l2n))*KM_TO_MPC*Hubble_h; // Mpc/h
        l3n = deltaT*aceleration_ejected_BH(p,pos[0]+0.5*k2n,r,v_tot,vel[0]+(0.5*l2n),Mvir,redshift,MBH,n_step,filenr)/(Hubble_h*KM_TO_MPC*SECONDS_TO_MEGAYRS); // km/s
     // Y-axis
        q3n = (deltaT*SECONDS_TO_MEGAYRS)*(vel[1]+(0.5*m2n))*KM_TO_MPC*Hubble_h; // Mpc/h
        m3n = deltaT*aceleration_ejected_BH(p,pos[1]+0.5*q2n,r,v_tot,vel[1]+(0.5*m2n),Mvir,redshift,MBH,n_step,filenr)/(Hubble_h*KM_TO_MPC*SECONDS_TO_MEGAYRS); // km/s
     // Z-axis
        c3n = (deltaT*SECONDS_TO_MEGAYRS)*(vel[2]+(0.5*d2n))*KM_TO_MPC*Hubble_h; // Mpc/h
        d3n = deltaT*aceleration_ejected_BH(p,pos[2]+0.5*c2n,r,v_tot,vel[2]+(0.5*d2n),Mvir,redshift,MBH,n_step,filenr)/(Hubble_h*KM_TO_MPC*SECONDS_TO_MEGAYRS); // km/s



        r = pow(  pow(pos[0]+k3n,2) + pow(pos[1]+q3n,2) + pow(pos[2]+c3n,2) , 0.5);
        v_tot = pow(  pow(vel[0]+l3n,2) + pow(vel[1]+m3n,2) + pow(vel[2]+d3n,2) , 0.5);
     // X-axis
        k4n = (deltaT*SECONDS_TO_MEGAYRS)*(vel[0]+l3n)*KM_TO_MPC*Hubble_h; // Mpc/h
        l4n = deltaT*aceleration_ejected_BH(p,pos[0]+k3n,r,v_tot,vel[0]+l3n,Mvir,redshift,MBH,n_step,filenr)/(Hubble_h*KM_TO_MPC*SECONDS_TO_MEGAYRS); // km/s
     // Y-axis
        q4n = (deltaT*SECONDS_TO_MEGAYRS)*(vel[1]+m3n)*KM_TO_MPC*Hubble_h; // Mpc/h
        m4n = deltaT*aceleration_ejected_BH(p,pos[1]+q3n,r,v_tot,vel[1]+m3n,Mvir,redshift,MBH,n_step,filenr)/(Hubble_h*KM_TO_MPC*SECONDS_TO_MEGAYRS); // km/s
     // Z-axis
        c4n = (deltaT*SECONDS_TO_MEGAYRS)*(vel[2]+d3n)*KM_TO_MPC*Hubble_h; // Mpc/h
        d4n = deltaT*aceleration_ejected_BH(p,pos[2]+c3n,r,v_tot,vel[2]+d3n,Mvir,redshift,MBH,n_step,filenr)/(Hubble_h*KM_TO_MPC*SECONDS_TO_MEGAYRS); // km/s


        pos[0] += (1./6.)*(k1n+(2.*k2n)+(2.*k3n)+k4n);
        pos[1] += (1./6.)*(q1n+(2.*q2n)+(2.*q3n)+q4n);
        pos[2] += (1./6.)*(c1n+(2.*c2n)+(2.*c3n)+c4n);

        vel[0] += (1./6.)*(l1n+(2.*l2n)+(2.*l3n)+l4n);
        vel[1] += (1./6.)*(m1n+(2.*m2n)+(2.*m3n)+m4n);
        vel[2] += (1./6.)*(d1n+(2.*d2n)+(2.*d3n)+d4n);
    }

   else
    {
         // For the Ejected wBHS we assume that they are moving in a radial orbit (i.e just one direction, the r-axis)
         k1n = (deltaT*SECONDS_TO_MEGAYRS)*(vel[0])*KM_TO_MPC*Hubble_h; // Mpc/h
         l1n = deltaT*aceleration_ejected_BH(p,pos[0],fabs(pos[0]),fabs(vel[0]),vel[0],Mvir,redshift,MBH,n_step,filenr)/(Hubble_h*KM_TO_MPC*SECONDS_TO_MEGAYRS); // km/s
         k2n = (deltaT*SECONDS_TO_MEGAYRS)*(vel[0]+(0.5*l1n))*KM_TO_MPC*Hubble_h; // Mpc/h
         l2n = deltaT*aceleration_ejected_BH(p,pos[0]+(0.5*k1n),fabs(pos[0]+(0.5*k1n)),fabs(vel[0]+(0.5*l1n)),vel[0]+(0.5*l1n),Mvir,redshift,MBH,n_step,filenr)/(Hubble_h*KM_TO_MPC*SECONDS_TO_MEGAYRS); // km/s
         k3n = (deltaT*SECONDS_TO_MEGAYRS)*(vel[0]+(0.5*l2n))*KM_TO_MPC*Hubble_h; // Mpc/h
         l3n = deltaT*aceleration_ejected_BH(p,pos[0]+(0.5*k2n),fabs(pos[0]+(0.5*k2n)),fabs(vel[0]+(0.5*l2n)),vel[0]+(0.5*l2n),Mvir,redshift,MBH,n_step,filenr)/(Hubble_h*KM_TO_MPC*SECONDS_TO_MEGAYRS); // km/s
         k4n = (deltaT*SECONDS_TO_MEGAYRS)*(vel[0]+l3n)*KM_TO_MPC*Hubble_h; // Mpc/h
         l4n = deltaT*aceleration_ejected_BH(p,pos[0]+k3n,fabs(pos[0]+k3n),fabs(vel[0]+(0.5*l2n)),vel[0]+l3n,Mvir,redshift,MBH,n_step,filenr)/(Hubble_h*KM_TO_MPC*SECONDS_TO_MEGAYRS); // km/s
         pos[0] += (1./6.)*(k1n+(2.*k2n)+(2.*k3n)+k4n);
         vel[0] += (1./6.)*(l1n+(2.*l2n)+(2.*l3n)+l4n);
     }
}


float aceleration_ejected_BH(int p, float coor, float pos, float vel, float vi, float Mvir, float redsfhit, float MBH, int n_step, int filenr)
{
    // pos --> Is the radial distance
    // coor --> Is the (x or , y or , z) position
    // vel --> Is the total velocity of the BH
    // vi --> Is the velocity of the BH in the coor x, or y, or z
    float aceleration_Gravity, aceleration_DF_STARS_DM, aceleration_GAS, aceleration_Total;
    aceleration_Gravity = compute_gravity_force(p,Mvir,coor,pos,redsfhit, n_step, filenr);
    aceleration_DF_STARS_DM = comute_a_DF_DM(p,pos,Mvir,redsfhit,MBH,vel,vi,n_step,filenr);
    aceleration_GAS = comute_a_DF_gas(p,pos,Mvir,redsfhit,MBH,vel,vi,filenr);
    aceleration_Total = aceleration_Gravity + aceleration_DF_STARS_DM + aceleration_GAS;
    if(isnan(aceleration_Total)==1) {printf("Error in the function aceleration_ejected_BH(): aceleration is nan (%E)! Term 1: compute_gravity_force = %E, Term 2: comute_a_DF_DM = %E, Term 3: comute_a_DF_gas = %E, pos = %E, vel = %f, Type of gal = %d\n",aceleration_Total, aceleration_Gravity, aceleration_DF_STARS_DM, aceleration_GAS, pos, vel, Gal[p].Type); exit(0);}
    return aceleration_Total;

}


float compute_gravity_force(int p, float Mvir, float coor, float pos, float redsfhit, int n_step, int filenr)
{

   // Here I compute the term of aceleration caused by the gravity
   // To understand the difference between coor and pos, check the comment in the function aceleration_ejected_BH
   float a_gravity, SF, Mt_r;
   SF = 1e-6; // Mpc "Softening length"

   // Total mass inside the radius pos
   Mt_r = DM_mass_inside_radi(p,fabs(pos),n_step,filenr)/Hubble_h;
   Mt_r += Stellar_mass_inside_radi_Bulge(p,fabs(pos),filenr);
   Mt_r += mass_inside_radii_disk(p,fabs(pos),filenr,1);
   Mt_r += mass_inside_radii_disk(p,fabs(pos),filenr,-1);

   a_gravity = ( -1.0 * (fabs(coor)/Hubble_h) * G_GRAVITY * Mt_r / pow(pow(pos/Hubble_h,2.0) + pow(SF/Hubble_h,2.0),3.0/2.0) ) * pow(KM_TO_MPC,2.0) * pow(SECONDS_TO_MEGAYRS,2.0) * Hubble_h; // (Mpc/h)/Myr^2

   if(isnan(a_gravity)==1) {printf("Error in the function compute_gravity_force variable: a_gravity is nan (%E)! DM_mass_inside_radi = %E, pos = %E, n_step = %d, filenr = %d, Mvir = %E\n",a_gravity,DM_mass_inside_radi(p,fabs(pos),n_step,filenr),pos,n_step,filenr,Mvir); exit(0);}
   if(isinf(a_gravity)==1) {printf("Error in the function compute_gravity_force variable: a_gravity is inf (%E)! DM_mass_inside_radi = %E, pos = %E, n_step = %d, filenr = %d, Mvir = %E\n",a_gravity,DM_mass_inside_radi(p,fabs(pos),n_step,filenr),pos,n_step,filenr,Mvir); exit(0);}

   return a_gravity;


}


float comute_a_DF_DM(int p, float pos, float Mvir, float redsfhit, float MBH, float vel, float vi, int n_step, int filenr)
{
   // Here I compute the aceleration caused by the dynamical friction caused by collisionless particles (stars and DM)
   // See Eq. 33 of Izquierdo-Villalba et al. 2020
   float a_DF;
   float sigma_DM, rho_DM, rho_ICM, rho_Tot, L_Lambda, X;

   sigma_DM = compute_sigma_DM_halo(p); // Units = [km/s]
   rho_DM = NFW_density_inside_radi(p,pos,n_step,filenr);// Units = [10^10 M_sun]/[Mpc]**3
   rho_ICM =density_of_ISOTHERMAL(p,Gal[p].ICM,fabs(pos),filenr);
   rho_Tot = rho_DM + rho_ICM;
   L_Lambda =3.1;
   X = fabs(vel)/(pow(2.0,0.5)*sigma_DM); // Units = [Adimensional]

   a_DF = -4. * M_PI *pow(G_GRAVITY,2.0) * (MBH/Hubble_h) * rho_Tot * vi * L_Lambda * ( erf(X) - (2.*X*exp(-X*X)/pow(M_PI,1./2.))) * pow(KM_TO_MPC,2.0) * pow(SECONDS_TO_MEGAYRS,2.0) * Hubble_h / pow(fabs(vel),3.0); // Units: (Mpc/h)/Myr^2

   // I do some checks to test if the dynamical friction has been computed correctly
   if(isnan(a_DF)==1) {printf("Error in file = %d step = %d: function comute_a_DF_DM variable: a_DF is nan (%E)! rho_DM = %E, L_Lambda = %E, X = %Evel = %f, rho_DM = %E, rho_ICM = %E, (2.*X*exp(-X*X)/pow(M_PI,1./2.)) = %f, erf(X) = %f \n",filenr,n_step,a_DF,rho_DM,rho_ICM,L_Lambda,X, vel, rho_DM,(2.*X*exp(-X*X)/pow(M_PI,1./2.)), erf(X)); exit(0);}
   if(isinf(a_DF)==1) {printf("Error in file = %d step = %d: function comute_a_DF_DM variable: a_DF is INF (%E)! rho_DM = %E, L_Lambda = %E, X = %E vel = %f, rho_DM = %E, rho_ICM = %E, (2.*X*exp(-X*X)/pow(M_PI,1./2.)) = %f, erf(X) = %f \n",filenr,n_step,a_DF,rho_DM,rho_ICM,L_Lambda,X, vel, rho_DM,(2.*X*exp(-X*X)/pow(M_PI,1./2.)), erf(X)); exit(0);}
   return a_DF;
}



float comute_a_DF_gas(int p, float pos, float Mvir, float redsfhit, float MBH, float vel, float vi, int filenr)
{
   float a_DF, cs, M, M_eq, f_M;
   float rho_gas,L_Lambda, sigma_DM;
   M_eq = 1.5;

   cs = 1.8*pow(1.+redsfhit,1./2.)*pow(1.e3*Mvir/Hubble_h,1./3.) * (Omega*Hubble_h*Hubble_h/0.14); // km/s
   L_Lambda =3.1;
   rho_gas = density_of_ISOTHERMAL(p,Gal[p].HotGas,fabs(pos),filenr);
   sigma_DM = compute_sigma_DM_halo(p); // Units = [km/s]
   M = fabs(vel)/cs;
   if(Gal[p].HotGas < 1e-9) return 0.;

   if(M<=0.8) f_M = 0.5*L_Lambda * ( erf(M/pow(2.,0.5)) - ( (pow(2./M_PI,0.5))*M*exp(-M*M/2.) ) );
   else if(0.8<M && M<=M_eq) f_M = 1.5*L_Lambda * ( erf(M/pow(2.,0.5)) - ( (pow(2./M_PI,0.5))*M*exp(-M*M/2.) ) );
   else if(M>M_eq) f_M = (0.5*(1.-pow(M,-2.)) ) + L_Lambda;
   else { printf("Error in comute_a_DF_gas(): M (%f) is not understtod in comute_a_DF_gas, M_eq = %f\n",M,M_eq); exit(0); }

    a_DF = -4. * M_PI *pow(G_GRAVITY,2.0) * (MBH/Hubble_h) * rho_gas * vi * f_M * pow(KM_TO_MPC,2.0) * pow(SECONDS_TO_MEGAYRS,2.0) * Hubble_h / pow(fabs(vel),3.0); // Units: (Mpc/h)/Myr^2
    if(isnan(a_DF)==1) {printf("Error in file = %d : function comute_a_DF_gas variable: comute_a_DF_gas is nan (%E)! \n",filenr,a_DF); exit(0);}
  return a_DF;
}




void save_file_with_kick_BH(int p, int nstep, float V_recoil, float a1, float a2, float BH1, float BH2, float qq, float BHfinal, float fgas_BHs, int filenr)
{
  float V_escape,zz,dt;
  V_escape = compute_escape_velocity(p,0.0);
  if((V_recoil/V_escape) > 1. && BHfinal > 0.)
      {
        FILE *fileout;
        dt = ( NumToTime(Gal[p].SnapNum) -  NumToTime(Gal[p].SnapNum+1) ) / STEPS;
        zz = redshift_at_t(ZZ[Gal[p].SnapNum],nstep,dt);

        float arr1[16] = {zz,
                          qq,
                          a1,
                          a2,
                          BH1,
                          BH2,
                          BHfinal,
                          fgas_BHs,
                          (Gal[p].BulgeMass + Gal[p].DiskMass),
                          Gal[p].BulgeMass,
                          Gal[p].DiskMass,
                          Gal[p].ColdGas,
                          Gal[p].Bulge_di,
                          V_recoil,
                          V_escape,
                          Gal[p].Mvir};
        char file_ejected_BH[1000];
        sprintf(file_ejected_BH, "/Recoil_BH/Hybrid_model/Ejected/Snap_%d", filenr);
        fileout = fopen(file_ejected_BH,"ab");
        fwrite(arr1,sizeof(arr1),1,fileout);
        fclose(fileout);
      }

}

#endif
