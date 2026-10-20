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



/** This module deals with the nuclear black holes
 *  All the physics included can be found summarized in Izquierdo-Villalba et al. 2020 (and references therein)
*/

//###########################################################################
//################# These are the Growth functions ##########################
//###########################################################################


void BH_growth_Hot_Gas(int p, int BH_indx, float AGNaccreted, float dt, int nstep, int filenr)
{
  // This function computes the mass growth of the black hole during the hot gas accretion
  // Given that there is not yet a model to deal with the spin evolution during the ADAF accretion, we assume the same spin evolution with the Thin disk accretion
  
  float m_dot, a, eta, epsilon, M_Edd_dot, MBH_i, M_dot_EDDrate, a_new;
  float fEdd, M_dot;
  int j;

  MBH_i = Gal[p].BlHo[BH_indx].BlackHoleMass[nstep];  // initial BH mass
  #ifdef SPINCOMPUTATION
    // Here I asumme that the first accretion is prograd! (See papers of Volonteri)
    a = fabs(Gal[p].BlHo[BH_indx].spin[nstep]);
    if(Gal[p].BlHo[BH_indx].fEDD[nstep] <= 1e-15) eta = compute_eta(a);
    else eta = Gal[p].BlHo[BH_indx].eta_BH[nstep];
    epsilon = eta;
  #else
     eta = 0.1; epsilon = eta;
  #endif

  // All this is based on the paper of Merloni and Heinz 2005
  
  M_Edd_dot = compute_M_dot_Edd(MBH_i,eta); // Notice that it should give the same as  M_dot_EDDrate = ((1.-eta)/eta)*MBH_i/T_EDD;
  M_dot = AGNaccreted/dt;
  a = Gal[p].BlHo[BH_indx].spin[nstep];
  //MBH_i = Gal[p].BlHo[BH_indx].BlackHoleMass[nstep];
  m_dot = M_dot / M_Edd_dot;
  if(m_dot<=m_DOT_CRIT) epsilon = eta * (m_dot/m_DOT_CRIT); // ADAF accretion
  else epsilon =  eta; // Thin disk acretion
  fEdd = m_dot * epsilon / eta; // Volonteri Silik Dubus 2015 // / eta;

  #ifdef RADIOLUMBHS
  Compute_Radio_Luminosities(p,nstep,BH_indx,MBH_i,M_dot,a,fEdd,eta);
  #endif
  Gal[p].BlHo[BH_indx].fEDD[nstep] = fEdd; //m_dot * epsilon / eta;
  Gal[p].BlHo[BH_indx].QuasarLum[nstep] += (Gal[p].BlHo[BH_indx].fEDD[nstep] * L_EDDINGTON_ERG * (MBH_i / Hubble_h));
  //a_new =  Gal[p].BlHo[BH_indx].spin[nstep] - ((2.0*Gal[p].BlHo[BH_indx].spin[nstep]*(AGNaccreted/dt)/MBH_i)*dt);
  #ifdef SPINCOMPUTATION
  // The spin only changes in the accretion phase
  // Notice that in this case, for the spin computation, we put v_sigma = 0
  // this is because the hot gas accretion is not linked with the bulge and we assume that is an isotropic accretion
  // i.e it is a chaotic accretion (as happens in Elliptical galaxies)
  Self_Gravity_Accretion_Sesana_Model(p,MBH_i,M_dot,fEdd,&a,&eta,epsilon,0.0,dt,filenr);
  #endif
  if(a<= 0.) a_new = 0.0;
  if(a>=1.) a_new = 0.998;
#ifdef EXTENTED_BH_OUTPUT
save_extended_BHproperties(p, nstep, BH_indx, filenr, 3, MBH_i ,MBH_i+((1.-eta)*AGNaccreted), MBH_i+((1.-eta)*AGNaccreted), AGNaccreted, eta);
#endif
  for(j=nstep; j<STEPS; j++)
   {  
   Gal[p].BlHo[BH_indx].BlackHoleMass[j] += ((1.-eta)*AGNaccreted);  
   Gal[p].BlHo[BH_indx].spin[j]  = a; 
   Gal[p].BlHo[BH_indx].eta_BH[j] = eta;
   #ifdef EXTENDED_BH_OUTPUT
   Gal[p].BlHo[BH_indx].MBH_acc_radiomode[j] += ((1.-eta)*AGNaccreted);
   #endif
   }

}


void BH_growth_Cold_Gas(int gal_num, int BH_indx, int n_step, float BHaccrete, int continue_growth, int filenr)
{

  // This function deals with the accretion of cold gas coming from the mergers or disk instabilities
  // We have three ways in which the black hole can consume the gas
  // The first model of growth is the established by Henriques2015, basically the BH accretes the matter instantaneously 
  // The second one is Eddington limited model in whcih the black hole grows fast (fedd = 1.0)
  // The third one is two phases model. The first one is the Eddington limited one but when the BH consumes a percentage of reservoir mass, the BH enters in a low-accretion regimen emulating the AGN accretion
  //
  if( (Gal[gal_num].BlHo[BH_indx].BlackHoleMass[n_step] > 1e-9) && ( (Gal[gal_num].BlHo[BH_indx].Mass_to_acc[n_step] + BHaccrete)>1e-11 ) )
  {
  // Note that we put a threshold of 1e-11 (0.1M_sun/h) to be sure of only BHs with enought mass in the reservoir grow

  float MBH_i, MBH_f, MaxBHMass_without_epsilon_losses, Mass_to_acc_total;
  float dt;
  float a, eta;

  dt = (NumToTime(Gal[gal_num].SnapNum) - NumToTime(Gal[gal_num].SnapNum + 1))/STEPS;
  #ifdef SPINCOMPUTATION
  // Here I asumme that the first accretion is prograd! (See papers of Volonteri)
  a = Gal[gal_num].BlHo[BH_indx].spin[n_step];
  eta = Gal[gal_num].BlHo[BH_indx].eta_BH[n_step];
  #else
   a = 0.0; eta = 0.1;
  #endif

  if(continue_growth<0 || BHaccrete>1e-12) {
          Gal[gal_num].BlHo[BH_indx].M_res_0  = Gal[gal_num].BlHo[BH_indx].Mass_to_acc[n_step] + BHaccrete;
          Gal[gal_num].BlHo[BH_indx].MBH_0 = Gal[gal_num].BlHo[BH_indx].BlackHoleMass[n_step];
          
          #ifdef SPINCOMPUTATION
          a = Gal[gal_num].BlHo[BH_indx].spin[n_step]; eta = compute_eta(a);  // Here I asumme that the first accretion is prograd!
          #endif
          Gal[gal_num].BlHo[BH_indx].t_acc_quiescent_phase = -1;
          //Gal[gal_num].BlHo[BH_indx].MBH_peak = 0.;
          Gal[gal_num].BlHo[BH_indx].MBH_second_phase = 0.;
        }
   // # Some global parameters #
  //MBH_i = Gal[gal_num].BlHo[BH_indx].BlackHoleMass[n_step]; // Mass of the BH after a merger 
   //MBH_f = MBH_i + (1.-eta)*Mass_to_acc_total; // Max BH mass
   Mass_to_acc_total = Gal[gal_num].BlHo[BH_indx].Mass_to_acc[n_step] + BHaccrete;
  //MaxBHMass_without_epsilon_losses = MBH_i + Mass_to_acc_total; // This will be used only for doing checkings!]

#ifdef EXTENDED_BH_OUTPUT
  if(continue_growth==-1)
  {
     if(Gal[gal_num].induced > 0) {
      Gal[gal_num].BlHo[BH_indx].M_acc_DI_mmi[n_step] += BHaccrete; 
     Gal[gal_num].BlHo[BH_indx].M_res_DI_mmi[n_step]  += BHaccrete;} //this will be constant
     else { 
      Gal[gal_num].BlHo[BH_indx].M_acc_DI_sec[n_step] += BHaccrete; 
      Gal[gal_num].BlHo[BH_indx].M_res_DI_sec[n_step] += BHaccrete; } 
     }
   if(continue_growth==-2){
    Gal[gal_num].BlHo[BH_indx].M_acc_Mergers[n_step] += BHaccrete;
    Gal[gal_num].BlHo[BH_indx].M_res_Mergers[n_step] +=BHaccrete;}
   #ifdef FEED_ON_STARS 
   //if accreting with the two phase model, otherwise directly in BH_growth_TDE function
   if(continue_growth==-3)  {
    Gal[gal_num].BlHo[BH_indx].M_acc_TDE[n_step] += BHaccrete; 
   Gal[gal_num].BlHo[BH_indx].M_res_TDE[n_step] +=BHaccrete;}
     if(continue_growth<-3 || continue_growth>5) printf("Wrong continue growth parameter; check extended_bh_output option");

  #else
  if(continue_growth<-2 || continue_growth>4) printf("Wrong continue growth parameter; check extended_bh_output option");
  #endif
#endif

  // CHECKS FOR THE SPIN
  if(a >=1) {printf("Error in BH_Growth_phase() initial condition a >=1: %f, MBH_i = %E, Gal[p].Snp = %d, n_step = %d Snp_last_merg = %d Step_last_merg = %d \n",a,MBH_i,Gal[gal_num].SnapNum,n_step,Gal[gal_num].Snap_last_merger,Gal[gal_num].step_last_merger); exit(0);}
  if( Gal[gal_num].BlHo[BH_indx].spin[n_step]>=1){printf("Error in BH_Growth_phase() initial condition a >=1: %f, MBH_i = %f, Gal[p].Snp = %d, n_step = %d Snp_last_merg = %d Step_last_merg = %d \n",Gal[gal_num].BlHo[BH_indx].spin[n_step],MBH_i,Gal[gal_num].SnapNum,n_step,Gal[gal_num].Snap_last_merger,Gal[gal_num].step_last_merger);}

   if(BlackHoleAccretionModel==0) Growth_BH_as_Henriques(gal_num, n_step, BHaccrete);
   //else if(BlackHoleAccretionModel==1) Growth_BH_Eddington_Limited(gal_num, n_step, BH_indx, eta, a, MBH_i, MBH_f, Mass_to_acc_total, MaxBHMass_without_epsilon_losses, dt, filenr);
   //else if(BlackHoleAccretionModel==2) Growth_BH_Two_Phases_Model( gal_num, n_step, BH_indx, eta, a, MBH_i, MBH_f, Mass_to_acc_total, MaxBHMass_without_epsilon_losses, dt, continue_growth, filenr);
   else if(BlackHoleAccretionModel==1) Growth_BH_Eddington_Limited(gal_num, n_step, BH_indx, eta, a,Mass_to_acc_total, dt, filenr);
   else if(BlackHoleAccretionModel==2) Growth_BH_Two_Phases_Model( gal_num, n_step, BH_indx, eta, a,Mass_to_acc_total, dt, continue_growth, filenr);
   else {printf("Error in BH_growth_Cold_Gas: The accretion model is not recognized!\n"); exit(0);}
 } // end the condition of the mass of BH and mass of the reservoir
}


void Merger_BH_as_Henriques(int t, int p, int n_step)
{
   // This function basically does as the version of L-Galaxies of Henriques 2015 does: An instantaneous merger!
   int j;
   for(j = n_step; j<STEPS; j++)
    {
        Gal[t].BlHo[0].BlackHoleMass[j] += Gal[p].BlHo[0].BlackHoleMass[j];
        Gal[t].BlHo[0].Mass_to_acc[j] = 0.0;
        Gal[t].BlHo[0].spin[j] = 0.0;
        Gal[p].BlHo[0].BlackHoleMass[j] = 0.0;
        Gal[p].BlHo[0].Mass_to_acc[j]   = 0.0;
        Gal[p].BlHo[0].spin[j] = 0.0;
     }
}

void Growth_BH_as_Henriques(int p, int n_step, float BHaccrete)
{
  // This function computes the instantaneous growth (Henriques et al. 2015 model)

  int j;
  for(j=n_step; j<STEPS; j++)
   {
       Gal[p].BlHo[0].BlackHoleMass[j] += BHaccrete;
       Gal[p].BlHo[0].Mass_to_acc[j] = 0.0;
       //if(j==n_step) Gal[p].QuasarLum[j] = log10(L_EDDINGTON_ERG * Gal[p].BlackHoleMass[j] / Hubble_h);
       //else Gal[p].QuasarLum[j] = 0.;
   }
}


void Growth_BH_Eddington_Limited(int p, int n_step, int BH_indx, float eta, float a, float Mass_to_acc_total, float dt, int filenr)
{
   // This model of growth assumes that the BH comsumes all the gas in the reservoir in an Eddington limited phase
   // For more details see Izquierdo-Villalba et al. 2020 or Volonteri et al. 2005
   int j, ii;
   float fEdd = 1.0; // The Eddington limited model assumes that all the time that the BH is accreting mass, it does it at Eddington rate, i.e fedd = 1.0
   float epsilon, M_dot, MBH_teo, MBH_i, dMBH_i, MBH_f , MaxBHMass_without_epsilon_losses;
   dMBH_i = 0.0;
   for(j=n_step; j<STEPS; j++)
    {  
       MBH_i  = Gal[p].BlHo[BH_indx].BlackHoleMass[j] + dMBH_i;
       MaxBHMass_without_epsilon_losses = MBH_i + Mass_to_acc_total;
       MBH_f = MBH_i + (1.-eta)*Mass_to_acc_total; // Notice that I re-compute the value of MBH_f given that the value of eta has changed after the gas accretion
       fEdd = 1.0;
       epsilon =  eta;
       MBH_teo = MBH_i + ( (1.-eta) * MBH_i * fEdd * dt / (epsilon * T_EDD) ); // Should be equivalent to MBH_teo = MBH_i*exp( dt * fEdd * (1.-eta) / (eta * T_EDD) );
       M_dot = (MBH_teo - MBH_i)/(dt*(1.-eta)); // Mass[internal] / seg[internal]  Accretion rate of the BH
       Gal[p].BlHo[BH_indx].BlackHoleMass[j] += MBH_teo - MBH_i + dMBH_i;
      #ifdef EXTENDED_BH_OUTPUT
      Gal[p].BlHo[BH_indx].MBH_acc_first_phase[j] = MBH_teo - MBH_i;
      #endif
       
       Mass_to_acc_total -= M_dot*dt; // This is the mass that the reservoir has lost given the accretion 
       if(MBH_teo>=MBH_f)
        {
           Gal[p].BlHo[BH_indx].BlackHoleMass[j] += MBH_f - MBH_teo;
           dMBH_i+= Gal[p].BlHo[BH_indx].BlackHoleMass[j]  - MBH_i;
           MaxBHMass_without_epsilon_losses = MBH_i+(MBH_f-MBH_i)/(1-eta);
           Mass_to_acc_total = 0.;
           M_dot = 0.0;
           fEdd = 0.0;
           for(ii=j+1; ii<STEPS; ii++)
           {Gal[p].BlHo[BH_indx].BlackHoleMass[ii] += dMBH_i;
            Gal[p].BlHo[BH_indx].fEDD[ii] = 0.0;
            Gal[p].BlHo[BH_indx].QuasarLum[ii] = 0.0;
            Gal[p].BlHo[BH_indx].Mass_to_acc[ii] = 0.0;
            Gal[p].BlHo[BH_indx].spin[j] = a;
            Gal[p].BlHo[BH_indx].eta_BH[j] = eta;
           }
           break;
        }
      dMBH_i+= Gal[p].BlHo[BH_indx].BlackHoleMass[j]  - MBH_i;
       if(Mass_to_acc_total<0) {printf("Error : Mass_to_acc_total in model_merger(): black hole growth is negative!"); exit(0);}
       if(Gal[p].BlHo[BH_indx].BlackHoleMass[j] > MaxBHMass_without_epsilon_losses) {
      printf("In accretion model 1: Snap = %d, Error in model mergers BH. Exceded the max mass,\n MBH = %E,MaxMass_without_epsilon_losses = %E, eta = %f, spin = %f, Mass_to_acc_total = %e ", Gal[p].SnapNum, log10(Gal[p].BlHo[BH_indx].BlackHoleMass[j]*1e10/Hubble_h),log10(MaxBHMass_without_epsilon_losses*1e10/Hubble_h), eta,a,Mass_to_acc_total);
      printf(" \nMBH_f = %f, MBH_teo = %f, \n", log10(MBH_f*1e10/Hubble_h), log10(MBH_teo*1e10/Hubble_h)); 
      exit(0);}

       Gal[p].BlHo[BH_indx].fEDD[j] = fEdd;
       Gal[p].BlHo[BH_indx].QuasarLum[j] = fEdd * L_EDDINGTON_ERG * (MBH_i / Hubble_h); // [1e40 erg/s];
       Gal[p].BlHo[BH_indx].Mass_to_acc[j] = Mass_to_acc_total;
       #ifdef SPINCOMPUTATION
         // The spin only changes in the accretion phase
         if(fEdd > 0.) Self_Gravity_Accretion_Sesana_Model(p,MBH_i,M_dot,fEdd,&a,&eta,epsilon,Gal[p].v_sigma_now,dt,filenr);
         if(fEdd > 0.) if(a>=1. || a<-1. || isnan(a) ==1) {printf("Error: During the gas accretion, the spin parameter is bigger than one or nan. Please check in model_mergers.c (David functions) , spin = %f, BlackHoleMass[j] = %E, MBH_i = %E,m_dot = %E snp = %d, j = %d, BH_idx = %d\n",a,Gal[p].BlHo[BH_indx].BlackHoleMass[j],MBH_i,M_dot,Gal[p].SnapNum,j,BH_indx); exit(0);}
       #else
   	 epsilon = 0.1; eta = 0.1; a = 0.0;
       #endif
       Gal[p].BlHo[BH_indx].spin[j] = a;
       Gal[p].BlHo[BH_indx].eta_BH[j] = eta;
       //MBH_i = Gal[p].BlHo[BH_indx].BlackHoleMass[j];
       #ifdef RADIOLUMBHS
       Compute_Radio_Luminosities(p,n_step,BH_indx,MBH_i,M_dot,a,fEdd,eta);
       #endif
      if (Gal[p].BlHo[BH_indx].BlackHoleMass[j]- MBH_i<0.0){
            printf("\n Start at nstep = %d",n_step);
            printf("\n Negative mass added in Two_Phase_Growth  prev MBH/ MBH_teo - MBH_i / dMBH =  %e, %e\n", Gal[p].BlHo[BH_indx].BlackHoleMass[j], Gal[p].BlHo[BH_indx].BlackHoleMass[j]- MBH_i);
            printf("\n At step %d log prev MBH_teo, MBH_f , eta, a = %e , %e, %e, %e\n",j, MBH_teo, MBH_f, eta,a);
            exit(0);}
       
      }
}

void Growth_BH_Two_Phases_Model(int p, int n_step, int BH_indx, float eta, float a, float Mass_to_acc_total,  float dt, int continue_growth, int filenr)
{
   // This model is based on  Volonteri et al. 2013  Bonoli et al. 2009  Marulli et al. 2008 Hopkins et al. 2006 model
   // The black hole first undergoes a Eddignton limited phase, until it consumes the Ff% of the mass in the reservoir. when this happens, the BH undesgoes a AGN phase (low-accretion rates)
   // We set the parameters of the model
   float Ff, eta_Hop, beta, to, t_quasar, fEdd, epsilon;
   float time_now, time_acc, M_res_min;
   float MBH_teo,  M_dot, M_acc_step, dMBH_i, MBH_i, MBH_f, MaxBHMass_without_epsilon_losses;
   int j, times_substep, kk ,ii;
   Ff = 0.7;                                                              // Fraction of the mass in the reservoir that the blak hole has to consume to reach the low-accretion regimen
   eta_Hop = 0.3;                                                         // Parameter eta for the self-regulated phase
   beta = 0.4;                                                            // Parameter beta for the self-regulated phase
   to =1.26e8*(1.-(1.8*(beta-0.6)))* YRS_TO_SEC * Hubble_h/UnitTime_in_s; // [Internal units]
   t_quasar = to*pow(eta_Hop,beta)/(beta*log(10));                        // tQuasasr for the self-regulated phase. [Internal units]
   times_substep = 1;                                                     // This helps you to do a thinner bin discretitazion in STEPS, by defalult it is 1 , i.e no thinner discretizations.
       // a value of 2 would mean reducing by a factor of 2 the size of time-STEP.

   //MaxBHMass_without_epsilon_losses = MBH_i + Mass_to_acc_total; // This will be used only for doing checkings!]
   if(to<0) {printf("Error in BH growth. to(%E) < 0!\n",to); exit(0);}
   //if(Gal[p].BlHo[BH_indx].t_acc_quiescent_phase<=0) Gal[p].BlHo[BH_indx].MBH_peak = Gal[p].BlHo[BH_indx].MBH_0  + Ff*(Gal[p].BlHo[BH_indx].M_res_0*(1-eta)); // This means that the black hole entered in the growth for the first time. So we set the time properties of the growth
   if(Gal[p].BlHo[BH_indx].t_acc_quiescent_phase<=0) M_res_min = (1-Ff)*Gal[p].BlHo[BH_indx].M_res_0; // This means that the black hole entered in the growth for the first time. So we set the time properties of the growth

   dMBH_i = 0.0;
    for(j=n_step; j<STEPS; j++)
    {
     kk=0;
    //for(kk=0; kk<times_substep;kk++) // This loop is to do more integrations inside the nstep. Just to be sure that the integration is super well done. I checked that with kk == 1 we obtain the same as with kk == times_substep
    //{
          MBH_i = Gal[p].BlHo[BH_indx].BlackHoleMass[j] + dMBH_i;
          MaxBHMass_without_epsilon_losses = MBH_i + Mass_to_acc_total;
          MBH_f = MBH_i + (1.-eta)*Mass_to_acc_total; // Notice that I re-compute the value of MBH_f given that the value of eta has changed after the gas accretion
          
          //if(MBH_i<Gal[p].BlHo[BH_indx].MBH_peak)
          if( Gal[p].BlHo[BH_indx].Mass_to_acc[j]>M_res_min)
           {
              // Phase of Eddington growth
              fEdd = 1.0;
              epsilon =  eta;
              MBH_teo = MBH_i + ((1.-eta) * MBH_i * fEdd * (dt/(1.0*times_substep)) / (epsilon * T_EDD));
              #ifdef EXTENDED_BH_OUTPUT
              Gal[p].BlHo[BH_indx].MBH_acc_first_phase[j] = MBH_teo - MBH_i;
              #endif
           }
          else // Phase of self-regulation 
           {  //printf("\n I am in the second phase");
              if(Gal[p].BlHo[BH_indx].t_acc_quiescent_phase < 0.)
                 { 
                    Gal[p].BlHo[BH_indx].t_acc_quiescent_phase = (NumToTime(Gal[p].SnapNum)) -  ( 1.0 * j * ( (NumToTime(Gal[p].SnapNum) - NumToTime(Gal[p].SnapNum+1))/STEPS) ) - (1.0 * kk * dt/(1.0*times_substep));
                    Gal[p].BlHo[BH_indx].MBH_second_phase = MBH_i;
                 }
                 time_now =  ((NumToTime(Gal[p].SnapNum)) -  ( 1.0 * j * ( (NumToTime(Gal[p].SnapNum) - NumToTime(Gal[p].SnapNum+1))/STEPS) ))  - (1.0 * kk * dt/(1.0*times_substep));
                 time_acc = (Gal[p].BlHo[BH_indx].t_acc_quiescent_phase - time_now) + (dt/(2.0*times_substep));
                 
                 #ifdef SEEDINGWITHMR2
                 if(time_acc<0.0) {Gal[p].BlHo[BH_indx].t_acc_quiescent_phase = time_now; time_acc = 0.0;} // This takes into account the slightly difference in time resolution of MS and MSII 
                 #endif
                 
                 if(time_acc < 0) {printf("Error in BH_Growth_phase() the accretion time in the second phase is negative! time_acc = %E, time_now = %E, t_acc_quiescent_phase = %E, MBH_ini = %E, BH_indx = %d\n",time_acc,time_now,Gal[p].BlHo[BH_indx].t_acc_quiescent_phase,log10(Gal[p].BlHo[BH_indx].BlackHoleMass[n_step]/(Hubble_h*1.e-10)),BH_indx); exit(0);}
                 fEdd = pow(1.+ pow(time_acc/t_quasar,1./2.),-2./beta);
                 if(fEdd<=FEDD_DOT_CRIT) epsilon = eta * (fEdd/FEDD_DOT_CRIT); // ADAF acretion
                 else epsilon =  eta; // Thin disk acretion
                 MBH_teo = MBH_i + ((1.-eta) * MBH_i * fEdd * (dt/(1.0*times_substep)) / (epsilon * T_EDD));
                 //printf("\n MBH_teo, MBH_i, added = %e , %e, %e",MBH_teo, MBH_i, ((1.-eta) * MBH_i * fEdd * (dt/(1.0*times_substep)) / (epsilon * T_EDD)));
                 #ifdef EXTENDED_BH_OUTPUT
                 Gal[p].BlHo[BH_indx].MBH_acc_second_phase[j] = MBH_teo - MBH_i;
                 #endif
           }
          M_dot = (MBH_teo - MBH_i)/((dt/(1.0*times_substep))*(1.-eta)); // Perego et al. 2009 [internal units], i.e (Mass[internal]/Time[internal])
          M_acc_step = M_dot*(dt/(1.0*times_substep));
          Mass_to_acc_total -= M_acc_step;
          Gal[p].BlHo[BH_indx].BlackHoleMass[j] += MBH_teo - MBH_i+dMBH_i; //prev #MP  = MBH_teo;
         //Correction of MBH_teo -> MBH_f in the save accretion event function
         #ifdef EXTENDED_BH_OUTPUT
         save_extended_BHproperties(p, j, BH_indx, filenr, continue_growth, MBH_i, MBH_f, MBH_teo, M_acc_step, eta);
         #endif

          if(MBH_teo>=MBH_f)
          {
           Gal[p].BlHo[BH_indx].BlackHoleMass[j] += MBH_f - MBH_teo;
           dMBH_i+= Gal[p].BlHo[BH_indx].BlackHoleMass[j]  - MBH_i;
           MaxBHMass_without_epsilon_losses = MBH_i+(MBH_f-MBH_i)/(1-eta);
           Mass_to_acc_total = 0.;
           M_dot = 0.0;
           fEdd = 0.0;
           //Gal[p].BlHo[BH_indx].MBH_peak = MBH_f;
           for(ii=j+1; ii<STEPS; ii++)
           {Gal[p].BlHo[BH_indx].BlackHoleMass[ii] += dMBH_i;
            Gal[p].BlHo[BH_indx].fEDD[ii] = 0.0;
            Gal[p].BlHo[BH_indx].QuasarLum[ii] = 0.0;
            Gal[p].BlHo[BH_indx].Mass_to_acc[ii] = 0.0;
            Gal[p].BlHo[BH_indx].spin[j] = a;
            Gal[p].BlHo[BH_indx].eta_BH[j] = eta;
           }
           break;
           }

          #ifdef RADIOLUMBHS
           Compute_Radio_Luminosities(p,n_step,BH_indx,MBH_i,M_dot,a,fEdd,eta);
          #endif

          if(Mass_to_acc_total<0) {printf("Error : Mass_to_acc_total in model_merger(): black hole growth is negative!, Mass_to_acc_total = %E\n",Mass_to_acc_total); exit(0);}
          if(Gal[p].BlHo[BH_indx].BlackHoleMass[j] > MaxBHMass_without_epsilon_losses) {printf("Error in Black Hole Accretion model 2: Error in model mergers BH. Exceeded the max mass, MBH = %E, MaxBHMass_without_epsilon_losses = %E\n", Gal[p].BlHo[BH_indx].BlackHoleMass[j], MaxBHMass_without_epsilon_losses); exit(0);}
          Gal[p].BlHo[BH_indx].fEDD[j] = fEdd;
          Gal[p].BlHo[BH_indx].QuasarLum[j] = fEdd * L_EDDINGTON_ERG * (MBH_i / Hubble_h);///1e40;
          Gal[p].BlHo[BH_indx].Mass_to_acc[j] = Mass_to_acc_total;
          #ifdef SPINCOMPUTATION
          // The spin only changes in the accretion phase
          if(fEdd > 0.) Self_Gravity_Accretion_Sesana_Model(p,MBH_i,M_dot,fEdd,&a,&eta,epsilon,Gal[p].v_sigma_now,dt/(1.0*times_substep),filenr);
          if(fEdd > 0.) if(a>=1. || a<-1. || isnan(a) ==1) {printf("Error: During the gas accretion, the spin parameter is bigger than one or nan. Please check in model_mergers.c (David functions) , spin = %f, BlackHoleMass[j] = %E, MBH_i = %E,m_dot = %E snp = %d, j = %d, BH_idx = %d\n",a,Gal[p].BlHo[BH_indx].BlackHoleMass[j],MBH_i,M_dot,Gal[p].SnapNum,j,BH_indx); exit(0);}
          #else
          eta = 0.1; a = 0.0;
          #endif

          Gal[p].BlHo[BH_indx].spin[j] = a;
          Gal[p].BlHo[BH_indx].eta_BH[j] = eta;
          dMBH_i += Gal[p].BlHo[BH_indx].BlackHoleMass[j]- MBH_i;
         if (Gal[p].BlHo[BH_indx].BlackHoleMass[j]- MBH_i<0.0){
             printf("\n Start at nstep = %d",n_step);
             printf("\n Negative mass added in Two_Phase_Growth  prev MBH/ MBH_teo - MBH_i / dMBH = %e , %e \n", Gal[p].BlHo[BH_indx].BlackHoleMass[j], Gal[p].BlHo[BH_indx].BlackHoleMass[j]- MBH_i);
             printf("\n At step %d log prev MBH_teo, MBH_f , eta, a = %e , %e, %e, %e\n",j, MBH_teo, MBH_f, eta,a);
             exit(0);}
         //if (fEdd>0.0) MBH_f = MBH_i+ (1.-eta)*Mass_to_acc_total;
          //if(Gal[p].BlHo[BH_indx].t_acc_quiescent_phase < 0) Gal[p].BlHo[BH_indx].MBH_peak = Gal[p].BlHo[BH_indx].MBH_0  + Ff*(Gal[p].BlHo[BH_indx].M_res_0*(1.0-eta));
          if(a>=1. || isnan(a) ==1) {printf("Error : During the BH gas accretion, the spin parameter is bigger than one or nan. Please check in model_mergers.c (David functions) , spin = %f, Gal[p].BlackHoleMass[j] = %E, MBH_i = %E,m_dot = %E snp = %d, j = %d, BH_indx=%d\n",a,Gal[p].BlHo[BH_indx].BlackHoleMass[j],MBH_i,M_dot,Gal[p].SnapNum,j,BH_indx); exit(0);}
    //} end the kk
  }// End substeps
}



#ifdef FEED_ON_STARS
void BH_growth_TDE(int p, int j, int BH_indx, float event_rate, float m_monochromatic, int filenr)
{
  // According to Antonini et al. 2018
  float zz, dt ,MBH;
  float dt_acc, t_TD, T_relax, time_now , unbound_debris_frac_escape , Mass_unbound_debris;
  float Mbh_teo, MBH_i,MBH_f,MBH_teo,M_dot, Mass_to_acc, eta, epsilon, fEdd, MaxBHMass_without_epsilon_losses;


  dt = ((NumToTime(Gal[p].SnapNum) - NumToTime(Gal[p].SnapNum + 1))/STEPS);
  zz = redshift_at_t(ZZ[Gal[p].SnapNum],j,dt);
  //MBHtot = (Gal[p].BlHo[0].BlackHoleMass[j]+ Gal[p].BlHo[1].BlackHoleMass[j])/Hubble_h; //[10^10 Msun]
  
  //T_relax = Gal[p].BlHo[0].T_relax[j];  // out of this function now
  T_relax=Gal[p].BlHo[0].t_TD_phase;
  //if (T_relax==0.0)
  //{
   //T_relax = compute_T_relax(p, j, m_monochromatic); // out if this function
  //}
  // for TDEs here only, later for EMRIS and plunges
  //Gal[p].BlHo[BH_indx].BlackHoleMass += event_rate*m_monochromatic/2.0*dt_acc; no Eddington limi, no loss other than initial disruption maximal scenario, no reservoir added
   time_now =  ((NumToTime(Gal[p].SnapNum)) -  ( 1.0 * j * ( (NumToTime(Gal[p].SnapNum) - NumToTime(Gal[p].SnapNum+1))/STEPS) ));
   if(Gal[p].BlHo[BH_indx].t_TD_phase < 0.)
      {
         Gal[p].BlHo[BH_indx].t_TD_phase = time_now;
         #ifdef EXTENDED_BH_OUTPUT
         Gal[p].BlHo[BH_indx].MBH_acc_TDE[j] = 0.0;
         #endif
      }
   // #ifdef SEEDINGWITHMR2
   // if(time_tde<0.0) {Gal[p].BlHo[BH_indx].t_TD_phase = time_now; time_tde = 0.0;} // This takes into account the slightly difference in time resolution of MS and MSII 
   // #endif
   
   // #ifdef SPINCOMPUTATION
   //    // Here I asumme that the first accretion is prograd! (See papers of Volonteri)
   //    a = fabs(Gal[p].BlHo[BH_indx].spin[j]);
   //    if(Gal[p].BlHo[BH_indx].fEDD[j] <= 1e-15) eta = compute_eta(a);
   //    else eta = Gal[p].BlHo[BH_indx].eta_BH[j];
   //    //epsilon = eta;
   //    epsilon = eta;//Dai 2019
   // #else
   //    eta = 0.1; epsilon = eta;//Dai 2019
   // #endif
   //Eddington accretion caping here possible if alternative feeding scenario is preffered. For now super-eddington accretion is allowed in TDEs
   //This is to not form permanent disk from regular TDE rates
   fEdd = 1.0; eta = 0.027; epsilon = 0.027; //Dai 2019, if eta=epsilon=0.1 adapted only 5% of the stars are consumed (according to observations)

   if ((Gal[p].BlHo[BH_indx].BlackHoleMass[j]> 1e-2/Hubble_h)) 
   {eta=0.0; epsilon =  dt_acc/T_EDD; fEdd=2.0;} //Reached hill mass limit the stars to be swallowed whole

   t_TD = Gal[p].BlHo[BH_indx].t_TD_phase; 
   dt_acc = dt;
   if(dt > T_relax) dt_acc = T_relax;
   if( t_TD-time_now > T_relax)
      {
      dt_acc = dt -(t_TD-time_now - T_relax);
      Gal[p].BlHo[BH_indx].t_TD_phase = -1.0;
      }
   MBH_i = Gal[p].BlHo[BH_indx].BlackHoleMass[j];  // initial BH mass
   Mass_to_acc  =  event_rate * m_monochromatic/2.0 * dt_acc;
   Mass_unbound_debris = event_rate * m_monochromatic/2.0 * dt_acc;
   MBH_f = MBH_i + (1.-eta)*Mass_to_acc; // Max BH mass
   MaxBHMass_without_epsilon_losses = MBH_f;
   MBH_teo = MBH_i + ((1.-eta) * MBH_i * fEdd * dt_acc / (epsilon * T_EDD));

  Gal[p].BlHo[BH_indx].BlackHoleMass[j] += MBH_teo-MBH_i;

  M_dot = (MBH_teo - MBH_i)/((dt_acc)*(1.-eta));
  Mass_to_acc -= M_dot*dt_acc;

  if((MBH_teo>=MBH_f) || (Gal[p].BlHo[BH_indx].BlackHoleMass[j]> 1e-15*Hubble_h))
   {
      Gal[p].BlHo[BH_indx].BlackHoleMass[j] += MBH_f-MBH_teo;
      Gal[p].BlHo[BH_indx].Mass_to_acc[j] += MBH_f-MBH_teo;
      #ifdef EXTENDED_BH_OUTPUT
      Gal[p].BlHo[BH_indx].M_res_TDE[j] += (MBH_f-MBH_teo)/(1-eta);
      #endif
      MaxBHMass_without_epsilon_losses = MBH_f;
      Mass_to_acc = 0.;
      Mass_unbound_debris = 0.;
      M_dot = 0.0; fEdd = 0.0; 
   }

  #ifdef EXTENDED_BH_OUTPUT
  save_extended_BHproperties(p, j, BH_indx, filenr, -3 ,MBH_i, MBH_f, MBH_teo, Mass_to_acc, Mass_unbound_debris);
  #endif

   // #ifdef RADIOLUMBHS
   // Compute_Radio_Luminosities(p,j,BH_indx,MBH_i,M_dot, Gal[p].BlHo[BH_indx].spin[j],fEdd,eta);
   // #endif
   if(Mass_to_acc<0) {printf("Error : Mass_to_acc in BH_growth_TDE(): black hole growth is negative!, Mass_to_acc = %f\n",Mass_to_acc); exit(0);}
   if(Gal[p].BlHo[BH_indx].BlackHoleMass[j] > MaxBHMass_without_epsilon_losses) {printf("Error in Black Hole Accretion model 2: Error in model mergers BH. Exceeded the max mass, MBH = %f, MaxBHMass_without_epsilon_losses = %f\n", Gal[p].BlHo[BH_indx].BlackHoleMass[j], MaxBHMass_without_epsilon_losses); exit(0);}
   
   // do feedback or preserve mass
   #ifdef TDE_FEEDBACK
   if (Mass_unbound_debris>0.0)  TDE_feedback(p, BH_indx, Mass_unbound_debris);
   #else
   Gal[p].EjectedMass += Mass_unbound_debris;
   #endif
   //Gal[p].BlHo[BH_indx].fEDD[j] = fEdd;
   //Gal[p].BlHo[BH_indx].QuasarLum[j] = fEdd * L_EDDINGTON_ERG * (MBH_i / Hubble_h);///1e40;
   //Gal[p].BlHo[BH_indx].Mass_to_acc[j] = Mass_to_acc_total;
   // End substeps
//printf("Finish\n");
}

  // char file_starscore[3000];
  // FILE *fileout;
  // sprintf(file_starscore,"%s/Starscores_z%1.2f_%d",OutputDir, ZZ[Gal[p].SnapNum],filenr);
  // fileout = mymalloc(fopen(file_starscore,"ab"), sizeof(int) * (8*ngal));
  //  fileout =fopen(file_starscore,"ab");

// void BH_growth_TDE_direct(int p, int nstep, int BH_indx, float  Mincluded, float events_fraction, float time, float m_monochromatic, int filenr)
// {
//       int jj;
//       float MBHtot ,Mdot_TDE, dt_step ,dM_over,Mdot_TDE_prev ,event_rate;
//       // global now  Stars2TDEeff = 1.0; //how many of the nuclear stars (Mincluded) are able to be consumed in total at the maximal rate
//       float eta_TDE = 0.3; //swallows have of the star, epsilon_rad=0.0
      
//       MBHtot = Gal[p].BlHo[0].BlackHoleMass[nstep]+Gal[p].BlHo[1].BlackHoleMass[nstep];
//       #ifndef TIME_DEPENDENT_TDE
//       if (Gal[p].BlHo[BH_indx].BlackHoleMass[nstep]>1.e-2*Hubble_h) {eta_TDE =1.0 ;} //Hills mass criterion at 10^8 Modot
//       event_rate = compute_TDE_rate(p ,jj, BH_indx, m_monochromatic, time);
//       //printf("\n Rate loaded %f [yr^-1] ", event_rate_TDE/1e12);
//       #endif


//       dt_step =(NumToTime(Gal[p].SnapNum) - NumToTime(Gal[p].SnapNum+1))/STEPS;
//       if(Gal[p].BlHo[BH_indx].t_TD_phase < 0.)
//          {//if((Gal[p].BlHo[BH_indx].M_core[nstep] > 0.) || (Gal[p].BlHo[BH_indx].T_relax[nstep] > 0.)) 
//          //flag initializing of TDE in the snapshot
//          Gal[p].BlHo[BH_indx].t_TD_phase = time;}
//       else  
//       {printf("Error in BH_growth_TDE_direct> check time-switch for TDE \n should be off on this mode of FEED_ON_STARS, t");}

//       // outside of function for now
//       //Central region is allowed 1) to be consumed 2) to not be disrupted > mass checks later so that it never exceeds rest of stellar mass components
//       if (Gal[p].BulgeMass< MBHtot)
//          {Gal[p].DiskMass += -1.0*(Mincluded-Gal[p].BlHo[0].M_core[nstep])*Hubble_h; //give back the core that is overwritten to the stars
//          if (Gal[p].DiskMass<0.0) {Gal[p].DiskMass=0.0; printf("\n DiskMass exhausted by TDEs!, yet continuing");}
//          }
//       else
//          {Gal[p].BulgeMass += -1.0*(Mincluded-Gal[p].BlHo[0].M_core[nstep])*Hubble_h; //give back the core that is overwritten to the stars
//          if (Gal[p].BulgeMass<0.0) {Gal[p].BulgeMass=0.0;printf("\n BulgeMass exhausted by TDEs!, yet continuing");}
//          }
//          //establish feeding from the core in the future
//       for (jj=nstep;jj<STEPS;jj++)
//       { //Relaxation time used at creation to propagate
//          Mdot_TDE_prev = Gal[p].BlHo[BH_indx].M_dot_TDE[jj]; 
//          #ifdef TIME_DEPENDENT_TDE
//          event_rate = compute_TDE_rate(p ,jj, BH_indx, m_monochromatic, time);
//          if (Gal[p].BlHo[BH_indx].BlackHoleMass[jj]>1.e-2*Hubble_h) {eta_TDE =1.0 ;} //Hills mass criterion at 10^8 Modot
//          Mdot_TDE = events_fraction *event_rate * m_monochromatic *eta_TDE;
//          //printf("\n Rate loaded %f [yr^-1] ", event_rate_TDE/1e12);
//          #endif
//         if((Mincluded>0.0)&&(Mdot_TDE>=Mdot_TDE_prev))
//          {  
//             Gal[p].BlHo[BH_indx].T_relax[jj] = Gal[p].BlHo[0].T_relax[nstep];   
//             Gal[p].BlHo[BH_indx].M_dot_TDE[jj] = Mdot_TDE;
//             Gal[p].BlHo[BH_indx].BlackHoleMass[jj] += (jj-nstep+1)*dt_step*(Mdot_TDE-Mdot_TDE_prev);//subtract also what was added previously, TDE_rates do not add up
//             #ifdef EXTENDED_BH_OUTPUT
//             Gal[p].BlHo[BH_indx].MBH_acc_TDE[jj] += dt_step*(Mdot_TDE-Mdot_TDE_prev);//subtract also what was added previously, TDE_rates do not add up
//             #endif
//             Gal[p].BlHo[0].M_core[jj] = Mincluded - (jj-nstep+1)*dt_step*Mdot_TDE/Stars2TDEefficiency;
//             if (Gal[p].BlHo[0].M_core[jj] <=0.0)
//             {
//                dM_over = -1.0* Gal[p].BlHo[0].M_core[jj]*Stars2TDEefficiency;  //excess should be given back to BH/accreted mass   
//                if (dM_over<dt_step*Mdot_TDE) 
//                {
//                Gal[p].BlHo[BH_indx].BlackHoleMass[jj] += -1.0*dM_over*(Mdot_TDE-Mdot_TDE_prev)/Mdot_TDE; //remove the extra mass falsfully added
//                #ifdef EXTENDED_BH_OUTPUT
//                Gal[p].BlHo[BH_indx].MBH_acc_TDE[jj] +=  -1.0*dM_over*(Mdot_TDE-Mdot_TDE_prev)/Mdot_TDE;
//                #endif
//                }
//                else
//                {
//                Gal[p].BlHo[BH_indx].BlackHoleMass[jj] += -1.0*(jj-nstep+1)*dt_step*(Mdot_TDE-Mdot_TDE_prev);//subtract also what was added previously, TDE_rates do not add up
//                #ifdef EXTENDED_BH_OUTPUT
//                Gal[p].BlHo[BH_indx].MBH_acc_TDE[jj] += -1.0*dt_step*(Mdot_TDE-Mdot_TDE_prev);//subtract also what was added previously, TDE_rates do not add up
//                #endif
//                //printf("\n Bad treatment in BH_growth_TDE_direct.\n log Mleftover = %f , logMBH = %f!",log10(dM_over),log10(Gal[p].BlHo[BH_indx].BlackHoleMass[jj]));
//                }
//                Gal[p].BlHo[BH_indx].M_dot_TDE[jj] = 0.0;
//                Gal[p].BlHo[0].M_core[jj] = 0.0; 
//                Gal[p].BlHo[0].T_relax[jj]=0.0;
//                Mincluded = 0.0;
               
//                //Gal[p].BlHo[BH_indx].t_TD_phase  = -1.0; uncomment to set flag to zero if the star-feeding ends within this snapshot
//                //break; //break here to allow for remaining M_core, T_relax and M_dot_TDE from past activity to remain in the following times
//             }
//             // do feedback or preserve mass
//             #ifdef TDE_FEEDBACK
//             if ((1-eta_TDE/eta_TDE)*Gal[p].BlHo[BH_indx].M_dot_TDE[jj]>0.0)  TDE_feedback(p, BH_indx, (1-eta_TDE/eta_TDE)*Gal[p].BlHo[BH_indx].M_dot_TDE[jj]);
//             #else
//             Gal[p].EjectedMass += (1-eta_TDE/eta_TDE)*Gal[p].BlHo[BH_indx].M_dot_TDE[jj];
//             #endif
//          }
//             // else{printf("\nThis should not be happening: new TDE rate less than half the previous one (if binary created, rate splits); \n Check BH_growth_TDE_direct");
//             //      printf("\nold (log) new (log): %f, %f ", log10(Mdot_TDE_prev) , log10(Mdot_TDE));}
//       }
//      if(Gal[p].BlHo[BH_indx].BlackHoleMass[nstep] <0.0) printf("\n MBH<0.0 %f after subtracting dM_leftover=%f,%f",log10(Gal[p].BlHo[BH_indx].BlackHoleMass[nstep]),log10(dM_over));
//      if(isnan(Gal[p].BlHo[BH_indx].BlackHoleMass[nstep])==1) printf("\n MBH nan after subtracting dM_leftover=%f, %f  ",log10(Gal[p].BlHo[BH_indx].BlackHoleMass[nstep]), log10(dM_over));
//   } //end BH_growth_TDE_direct


void BH_growth_TDE_instant(int p, int nstep, int BH_indx, float  Mincluded, float events_fraction, float time, float m_monochromatic, int filenr)
{
      float MBHtot ,Mdot_TDE, dt_step, event_rate, unbound_debris , MBH_f, MBH_teo;
      // global now  Stars2TDEeff = 1.0; //how many of the nuclear stars (Mincluded) are able to be consumed in total at the maximal rate
      float eff_TDE = 0.3; //swallows have of the star, epsilon_rad=0.0
      float frac_trapped = 0.0;

      event_rate = compute_TDE_rate(p , nstep, BH_indx, m_monochromatic, time) ;
      if (Gal[p].BlHo[BH_indx].BlackHoleMass[nstep]>1.e-2*Hubble_h) {eff_TDE =1.0;} //Hills mass criterion at 10^8 Modot

      Mdot_TDE = events_fraction *event_rate * m_monochromatic *eff_TDE;
      MBHtot = Gal[p].BlHo[0].BlackHoleMass[nstep]+Gal[p].BlHo[1].BlackHoleMass[nstep];


      dt_step =(NumToTime(Gal[p].SnapNum) - NumToTime(Gal[p].SnapNum+1))/STEPS;
      if(Gal[p].BlHo[BH_indx].t_TD_phase < 0.) Gal[p].BlHo[BH_indx].t_TD_phase = time;
      MBH_teo = Gal[p].BlHo[BH_indx].BlackHoleMass[nstep]+Mdot_TDE*dt_step;
      MBH_f =  Gal[p].BlHo[BH_indx].BlackHoleMass[nstep]+Mincluded*eff_TDE;
      if (Gal[p].BulgeMass<= MBHtot)
         {
         Gal[p].DiskMass += -1.0*Mdot_TDE/eff_TDE*dt_step; //give back the core that is overwritten to the stars
         if (Gal[p].DiskMass<0.0) {
         MBH_teo = MBH_f;
         printf("\n DiskMass exhausted by TDEs!, yet continuing,\n time= %f : MBH=%e , event_rate s=%e, dt_step=%f \n",time, MBHtot, event_rate/1e6 * MYRs_TO_INTERNAL,dt_step);
          //exit(0);
         Mdot_TDE += Gal[p].DiskMass/dt_step*eff_TDE;  //regulate Mdot of TDE for that instance
         Gal[p].DiskMass = 0.0;
         Gal[p].StellarDiskRadius = 0.0;
         }
         else{
         Gal[p].StellarDiskRadius = Gal[p].StellarDiskRadius*pow( 1-Mdot_TDE/eff_TDE*dt_step/(Gal[p].DiskMass+Mdot_TDE/eff_TDE*dt_step),0.5); //reduce the radius, maintain density, Hiils 75
         }
         }
      else
      {  
         Gal[p].BulgeMass += -1.0*Mdot_TDE/eff_TDE*dt_step; //give back the core that is overwritten to the stars
         if (Gal[p].BulgeMass<0.0) {
            MBH_teo = MBH_f;
            //printf("\n BulgeMass exhausted by TDEs!, yet continuing,\n time=%f : MBH=%e , event_rate =%e, dt_step=%f \n",time, MBHtot, event_rate/1e6 * MYRs_TO_INTERNAL,dt_step);
            //exit(0);
            Mdot_TDE += Gal[p].BulgeMass/dt_step*eff_TDE; //regulate Mdot of TDE for that instance
            Gal[p].BulgeMass = 0.0;
            Gal[p].BulgeSize = 0.0;
         }
         else{
         Gal[p].BulgeSize = Gal[p].BulgeSize*pow(1 - Mdot_TDE/eff_TDE*dt_step/(Gal[p].BulgeMass + Mdot_TDE/eff_TDE*dt_step),0.5); //reduce the radius, maintain density, Hiils 75, it could be from 0.5 to 0.33, but use 0.5 if -1.5 used for scaling relation co-relation
         }  
      }
      if (Mdot_TDE<0.0) printf("\n Negative Mdot in TDE growth: %f",Mdot_TDE);
      #ifdef EXTENDED_BH_OUTPUT
      //Gal[p].BlHo[BH_indx].MBH_acc_TDE[nstep] += dt_step*Mdot_TDE; //subtract also what was added previously, TDE_rates do not add up
      //printf("\n MdotTDE*dt = %e , %e",dt_step*Mdot_TDE, Mdot_TDE);
      // void save_extended_BHproperties(int p, int n_step, int BH_indx, int filenr, int index_growth , float MBH_i, float MBH_f, float MBH_teo, float M_acc_step, float eta)
      save_extended_BHproperties(p, nstep, BH_indx,filenr,-3, Gal[p].BlHo[BH_indx].BlackHoleMass[nstep], MBH_f, MBH_teo, (Mdot_TDE/eff_TDE)*dt_step, eff_TDE);
      #endif
      Gal[p].BlHo[BH_indx].M_dot_TDE[nstep] = Mdot_TDE; //subtract also what was added previously, TDE_rates do not add up
      Gal[p].BlHo[BH_indx].BlackHoleMass[nstep] += dt_step*Mdot_TDE; 
     // else{printf("\nThis should not be happening: new TDE rate less than half the previous one (if binary created, rate splits); \n Check BH_growth_TDE_direct");
     // printf("\nold (log) new (log): %f, %f ", log10(Mdot_TDE_prev) , log10(Mdot_TDE));}

   // unbound_debris = (1-eta_TDE/eta_TDE)*Gal[p].BlHo[BH_indx].M_dot_TDE[nstep]*dt_step;
   // #ifdef TDE_FEEDBACK
   // if (frac_trapped>0.0) BH_growth_Cold_Gas(p,0,nstep,unbound_debris*frac_trapped,-3,filenr); //different channel, gas accretion
   // if ((unbound_debris >0.0)&&(frac_trapped<1.0)) TDE_feedback(p, BH_indx,unbound_debris*(1-frac_trapped));
   // #else
   // if (frac_trapped>0.0) BH_growth_Cold_Gas(p,0,nstep,unbound_debris*frac_trapped,-3,filenr); 
   // Gal[p].EjectedMass += (1-frac_trapped)*unbound_debris;
   // #endif

   if(Gal[p].BlHo[BH_indx].BlackHoleMass[nstep] <0.0) printf("\n MBH<0.0  after subtracting , MBH=%e , Mdot_TDE=%e",MBHtot, Mdot_TDE);
   if(isnan(Gal[p].BlHo[BH_indx].BlackHoleMass[nstep])==1) printf("\n MBH nan after subtracting , MBH=%e , Mdot_TDE=%e",MBHtot, Mdot_TDE);
  } //end BH_growth_TDE_instant


/* The unbound material escapes back to the intergalactic medium ,  usually with high temperature T=10^6 K 
   The amount of mater and energy from this process should be conserved


/* The unbound material escapes back to the intergalactic medium ,  usually with high temperature T=10^6 K 
   The amount of mater and energy from this process should be conserved
*/
#ifdef TDE_FEEDBACK
void TDE_feedback(int p, int j, int BH_indx, double unbound_debris)
{
  double CentralVvir, MergeCentralVvir=0., EjectVmax, EjectVvir, TDE_Energy, Reheat_Energy, fac;
  double reheated_mass=0., ejected_mass=0.;
  double TDE_EjectionEfficiency,  TDE_ReheatingEpsilon , TDEPreVelocity ,Vdebris, Vvir;
  /* TDE FEEDBACK MODEL  based on SN*/ 
  /* with feedback location cold gas only and no central galaxy considered*/

	Vvir = Gal[p].Vvir; // main halo Vvir
   Vdebris = TDEPreVelocity * (Gal[p].BlHo[BH_indx].BlackHoleMass[j]*1.0e-10/Hubble_h)/1e6; //dependence on stellar mass is too weak (m )^-0.07 Jonker 2021 Review book page 399

   if (Gal[Gal[p].CentralGal].Type == 0)
   reheated_mass = TDE_ReheatingEpsilon * unbound_debris *
   (.5+1./pow(Gal[Gal[p].CentralGal].Vmax/TDEPreVelocity,ReheatSlope));
   else
   reheated_mass = TDE_ReheatingEpsilon * unbound_debris *
   (.5+1./pow(Gal[Gal[p].CentralGal].InfallVmax/TDEPreVelocity,ReheatSlope));

   if (reheated_mass * Gal[Gal[p].CentralGal].Vvir * Gal[Gal[p].CentralGal].Vvir > unbound_debris * Vdebris * Vdebris )
   reheated_mass =  unbound_debris * Vdebris * Vdebris / (Gal[Gal[p].CentralGal].Vvir * Gal[Gal[p].CentralGal].Vvir);
	 
	if(reheated_mass > Gal[p].ColdGas)
	  reheated_mass = Gal[p].ColdGas;

   if(reheated_mass > Gal[p].BulgeMass)
	  reheated_mass = Gal[p].BulgeMass;

   if(reheated_mass > Gal[p].DiskMass)
	  reheated_mass = Gal[p].DiskMass;

  /* Determine ejection (for FeedbackEjectionModel 2 we have the dependence on Vmax)
   * Guo2010 - eq 22
   * Note that satellites can now retain gas and have their own gas cycle*/

  if (Gal[Gal[p].CentralGal].Type == 0)
    {
      EjectVmax=Gal[p].Vmax;
      EjectVvir=Gal[p].Vvir;// main halo Vvir
    }
  else
    {
      EjectVmax=Gal[Gal[p].CentralGal].InfallVmax;
      EjectVvir=Gal[Gal[p].CentralGal].Vvir; //central subhalo Vvir
    }

  if(TDE_FeedbackEjectionModel == 0)
    {
      ejected_mass = (TDE_EjectionEfficiency * unbound_debris *
	      min(1./TDE_EjectionEfficiency, .5+1/pow(EjectVmax/Vdebris,EjectSlope)) -
	      reheated_mass*EjectVvir*EjectVvir) /(EjectVvir*EjectVvir);
    }
  else if(TDE_FeedbackEjectionModel == 1)//the ejected material is assumed to have V_tde
    {
      TDE_Energy = .5 * unbound_debris * Vdebris * Vdebris;
      Reheat_Energy = .5 * reheated_mass * EjectVvir * EjectVvir;

      ejected_mass = (TDE_Energy - Reheat_Energy)/(0.5 * TDE_EjectionEfficiency*TDE_Energy);

      //if VSN^2<Vvir^2 nothing is ejected
      if(TDE_EjectionEfficiency * TDE_Energy< TDE_Energy - Reheat_Energy)
	   ejected_mass =0.0;
    }

  // Finished calculating mass exchanges, so just check that none are negative
  if (reheated_mass < 0.0) reheated_mass = 0.0;
  if (ejected_mass < 0.0) ejected_mass = 0.0;

  /* Update For Feedback */
  /* update cold, hot, ejected gas fractions and respective metallicities
   * there are a number of changes introduced by Guo2010 concerning where
   * the gas ends up */

  //ejected_mass = 0.01*Gal[centralgal].HotGas;
  if (reheated_mass + ejected_mass > 0.)
    update_from_feedback(p, Gal[p].CentralGal, reheated_mass, ejected_mass);
}
#endif // END TDE_FEEDBACK



// float compute_r_ini_BlackHole(int p, int Smooth_accretion)
// {  
//    float Reff;
//    if(Smooth_accretion > 0) 
//     {
//       if(Gal[p].DiskMass > 1e-9) Reff =  1.68*Gal[p].StellarDiskRadius/3.; //[Mpc/h]
//       else Reff =  1.68*Gal[p].ColdGasRadius/3.; // [Mpc/h];
//     }
//    else 
//     { if(Gal[p].BulgeMass > 1e-9) Reff = Gal[p].BulgeSize; // [Mpc/h]
//       else if(Gal[p].DiskMass > 1e-9) Reff =  1.68*Gal[p].StellarDiskRadius/3.; //[Mpc/h]
//       else Reff =  1.68*Gal[p].ColdGasRadius/3.; // [Mpc/h];
//     }
//    return Reff;
//}

// float extract_tidal_radii_from_log_normal(float r_t)
// {
//  float x, p_x, sigma, mu, y, y_Max, binsize, r_t_random;
//  int acepted;
//  sigma = 0.3;//#pow(0.4,0.5);
//  y_Max = -1;
//  mu = r_t;
//  acepted = -1;
//  binsize = 10*r_t/100.;

//  for(x=0;x<(10*r_t);x = x + binsize)
//   { 
//     //p_x = pow(sigma * x * pow(2.0*M_PI,0.5),-1.0) * exp(-pow(log(x)-mu,2.0)/(2.0*pow(sigma,2.0)));
//     p_x = pow(sigma * pow(2.0*M_PI,0.5),-1.0) * exp(-pow(x-mu,2.0)/(2.0*pow(sigma,2.0)));
//     if(p_x > y_Max) y_Max = p_x;
//   }

//  while(acepted < 0)
//   {
//         x = (10*r_t)*((double)rand()/((double)RAND_MAX+1));
//         //p_x = pow(sigma * x * pow(2.0*M_PI,0.5),-1.0) * exp(-pow(log(x)-mu,2.0)/(2.0*pow(sigma,2.0)));
//         p_x = pow(sigma * pow(2.0*M_PI,0.5),-1.0) * exp(-pow(x-mu,2.0)/(2.0*pow(sigma,2.0)));
//         y = y_Max*((double)rand()/((double)RAND_MAX+1));
//         if(y<=p_x){ r_t_random = x; acepted = +1; }
//   }

//  return r_t_random;
// }


#endif //FEED_ON_STARS


#ifdef EXTENDED_BH_OUTPUT 
void save_extended_BHproperties(int p, int n_step, int BH_indx, int filenr, int index_growth , float MBH_i, float MBH_f, float MBH_teo, float M_acc_step, float eta)
{
   int IND;
   double MBH_acc;
   float  Mres , Mbh;
   float frac_mergers = 0.0;
   float frac_di = 0.0;

   char file_extBH[3000];
   FILE *fileout;
   //sprintf(file_extBH,"%s/Extended_BH",OutputDir);
   //fileout = fopen(file_extBH,"ab");
   sprintf(file_extBH,"%s/extBH_%d_z%1.2f_%d",OutputDir, n_step, ZZ[Gal[p].SnapNum],filenr);
   

   if(index_growth==3){  //radiomode accretion
      IND=4; 
      Mres = 0.0;
      MBH_acc= MBH_teo - MBH_i;}

   if (index_growth==4){  //GW loss
      IND=5; 
      Mres = 0.0;
      MBH_acc= MBH_teo - MBH_i;
   }

   #ifdef FEED_ON_STARS
   float frac_tde = 0.0;
   if (index_growth==5){  //Direct star accretion, no gas reservoir TDEs
      IND=6;
      Mres = 0.0;
      MBH_acc= MBH_teo - MBH_i;
   }
   #endif

   MBH_acc= MBH_teo - MBH_i;
   if(MBH_teo>=MBH_f) MBH_acc= MBH_f - MBH_i;
   if((MBH_acc<0.0)&&(index_growth!=4)){printf("\n\n !!Error in save_extended_BHproperties / negative MBH_acc: %f, MBH_i = %E, MBH_f = %E, MBH_teo = %E,  M_acc_step = %E, Mbh_acc = %E,  type growth = %d, n_step = %d \n", MBH_acc, MBH_i,MBH_f, MBH_teo,   M_acc_step, MBH_acc   , index_growth,n_step); exit(0);}
   // #ifndef BASICSEEDING
   // if((Mbh<BlackHoleSeedMass*0.5*Hubble_h*1e-10)&&(index_growth!=4)){printf("\n\n !!Error in save_extended_BHproperties / Mbh less than half the seed mass: %E, MBH_i = %E, MBH_f = %E, MBH_teo = %E, M_acc_step = %f, Mbh_acc = %f,  type growth = %d, n_step = %d \n", BlackHoleSeedMass*0.99*Hubble_h*1e-10, MBH_i,MBH_f, MBH_teo, M_acc_step, MBH_acc , index_growth,n_step); exit(0);}
   // #endif

   //for saving
   Mbh = MBH_teo;
   if(MBH_teo>MBH_f) Mbh = MBH_f;

  //track reservoirs. in installed versions of this function check also update_BH_properties so that M_res_TYPE[0] = M_res_TYPE[STEPS-1] in snapshot transition
  if (n_step>0)
  { Gal[p].BlHo[BH_indx].M_res_Mergers[n_step] = Gal[p].BlHo[BH_indx].M_res_Mergers[n_step-1];
    Gal[p].BlHo[BH_indx].M_res_DI_mmi[n_step] = Gal[p].BlHo[BH_indx].M_res_DI_mmi[n_step-1];
    Gal[p].BlHo[BH_indx].M_res_DI_sec[n_step] = Gal[p].BlHo[BH_indx].M_res_DI_sec[n_step-1];
    #ifdef FEED_ON_STARS
    Gal[p].BlHo[BH_indx].M_res_TDE[n_step] = Gal[p].BlHo[BH_indx].M_res_TDE[n_step-1];
    #endif
  }
   if(index_growth==1){ //continue growth
      IND=0;
      if (Gal[p].BlHo[BH_indx].M_res_Mergers[n_step]>0.0) frac_mergers = Gal[p].BlHo[BH_indx].M_res_Mergers[n_step] / Gal[p].BlHo[BH_indx].Mass_to_acc[n_step];
      Gal[p].BlHo[BH_indx].M_res_Mergers[n_step] -=  M_acc_step * frac_mergers;
      Gal[p].BlHo[BH_indx].MBH_acc_Mergers[n_step] += MBH_acc * frac_mergers;
      ////we allow for negative values on the TOTAL reservoir since we can borrow from other sub-reservoirs
      if(MBH_teo>=MBH_f) {
         Gal[p].BlHo[BH_indx].M_res_Mergers[n_step] = 0.0;
         Gal[p].BlHo[BH_indx].MBH_acc_Mergers[n_step] = (MBH_f - MBH_i)* frac_mergers; // to guarantee mass conservation
      }
      Mres=Gal[p].BlHo[BH_indx].M_res_Mergers[n_step];

      float arr_acca[8] = {Gal[p].HaloNr,Gal[p].HeapIndex,IND, MBH_acc*frac_mergers, Mres,Gal[p].BlHo[BH_indx].fEDD[n_step],eta};

      if(Gal[p].induced > 0){  //induced DI from minor mergers
         IND=1;
      if (Gal[p].BlHo[BH_indx].M_res_DI_mmi[n_step]>0.0) 
      {frac_di = Gal[p].BlHo[BH_indx].M_res_DI_mmi[n_step] / Gal[p].BlHo[BH_indx].Mass_to_acc[n_step];
      Gal[p].BlHo[BH_indx].M_res_DI_mmi[n_step] -= M_acc_step *frac_di;
      Gal[p].BlHo[BH_indx].MBH_acc_DI_mmi[n_step] += MBH_acc* frac_di ;
      Mres=Gal[p].BlHo[BH_indx].M_res_DI_mmi[n_step];
      }
      }
      else{ IND=2;   // secular DI 
      if (Gal[p].BlHo[BH_indx].M_res_DI_sec[n_step]>0.0) 
            {frac_di = Gal[p].BlHo[BH_indx].M_res_DI_sec[n_step] / Gal[p].BlHo[BH_indx].Mass_to_acc[n_step];
            Gal[p].BlHo[BH_indx].M_res_DI_sec[n_step] -= M_acc_step* frac_di;
            Gal[p].BlHo[BH_indx].MBH_acc_DI_sec[n_step] += MBH_acc*frac_di;
            Mres=Gal[p].BlHo[BH_indx].M_res_DI_sec[n_step];
            }
      }

      if(MBH_teo>=MBH_f) 
      {  Mres= 0.0;
      if(Gal[p].induced > 0){
         Gal[p].BlHo[BH_indx].M_res_DI_mmi[n_step]  = 0.0;
         Gal[p].BlHo[BH_indx].MBH_acc_DI_mmi[n_step]  = (MBH_f - MBH_i)*frac_di;//to guarantee mass conservation
         }
      else {Gal[p].BlHo[BH_indx].M_res_DI_sec[n_step] = 0.0;
         Gal[p].BlHo[BH_indx].MBH_acc_DI_sec[n_step]  = (MBH_f - MBH_i)*frac_di;// to guarantee mass conservation
      }
      }
     float arr_accb[8] = {Gal[p].HaloNr,Gal[p].HeapIndex, IND, MBH_acc*frac_di, Mres, Mbh, Gal[p].BlHo[BH_indx].fEDD[n_step], eta};

      #ifdef FEED_ON_STARS
      IND=3;   // secular DI 
     if (Gal[p].BlHo[BH_indx].M_res_TDE[n_step] >0.0) frac_tde = Gal[p].BlHo[BH_indx].M_res_TDE[n_step] / Gal[p].BlHo[BH_indx].Mass_to_acc[n_step];
          Gal[p].BlHo[BH_indx].M_res_TDE[n_step] -= M_acc_step* frac_tde;
          Gal[p].BlHo[BH_indx].MBH_acc_TDE[n_step] += MBH_acc*frac_tde;

      if(MBH_teo>=MBH_f) {Gal[p].BlHo[BH_indx].M_res_TDE[n_step]  = 0.0;
            Gal[p].BlHo[BH_indx].MBH_acc_TDE[n_step] = (MBH_f - MBH_i)*frac_tde;//to guarantee mass conservation
            }
      Mres=Gal[p].BlHo[BH_indx].M_res_TDE[n_step];
      float arr_accc[8] = {Gal[p].HaloNr,Gal[p].HeapIndex, IND, MBH_acc*frac_tde, Mres, Mbh, Gal[p].BlHo[BH_indx].fEDD[n_step], eta};
      #endif

     fileout = fopen(file_extBH,"ab");
     fwrite(arr_acca,sizeof(arr_acca),1,fileout);
     fwrite(arr_accb,sizeof(arr_accb),1,fileout);
     #ifdef FEED_ON_STARS 
     fwrite(arr_accc,sizeof(arr_accc),1,fileout);
     #endif 
     fclose(fileout);
      }
else{
      if(index_growth==-2) //accretion via mergers
      { 
      IND=0;
      Gal[p].BlHo[BH_indx].M_res_Mergers[n_step] -= M_acc_step;
      Gal[p].BlHo[BH_indx].MBH_acc_Mergers[n_step] += MBH_acc ;
      if(MBH_teo>=MBH_f) {
         Gal[p].BlHo[BH_indx].M_res_Mergers[n_step] = 0.0;
         Gal[p].BlHo[BH_indx].MBH_acc_Mergers[n_step] = (MBH_f - MBH_i); //to guarantee mass conservation
      }
      Mres=Gal[p].BlHo[BH_indx].M_res_Mergers[n_step];
      }
    
      if(index_growth==-1){
      if(Gal[p].induced > 0){  //induced DI from minor mergers
         IND=1;
      Gal[p].BlHo[BH_indx].M_res_DI_mmi[n_step] -= M_acc_step;
      Gal[p].BlHo[BH_indx].MBH_acc_DI_mmi[n_step] += MBH_acc ;
      Mres=Gal[p].BlHo[BH_indx].M_res_DI_mmi[n_step];
      }
      else{ IND=2;   // secular DI 
            Gal[p].BlHo[BH_indx].M_res_DI_sec[n_step] -= M_acc_step;
            Gal[p].BlHo[BH_indx].MBH_acc_DI_sec[n_step] += MBH_acc ;
            Mres=Gal[p].BlHo[BH_indx].M_res_DI_sec[n_step];
      }
      if(MBH_teo>=MBH_f) 
      {  Mres= 0.0;
      if(Gal[p].induced > 0){Gal[p].BlHo[BH_indx].M_res_DI_mmi[n_step]  = 0.0;
         Gal[p].BlHo[BH_indx].MBH_acc_DI_mmi[n_step]  = (MBH_f - MBH_i);//to guarantee mass conservation
         }
      else {Gal[p].BlHo[BH_indx].M_res_DI_sec[n_step] = 0.0;
         Gal[p].BlHo[BH_indx].MBH_acc_DI_sec[n_step]  = (MBH_f - MBH_i); //to guarantee mass conservation
       }
      }
      }

      #ifdef FEED_ON_STARS
      if(index_growth==-3) //accretion via mergers
      { 
      IND=3;
      Gal[p].BlHo[BH_indx].M_res_TDE[n_step] -=  M_acc_step;
      Gal[p].BlHo[BH_indx].MBH_acc_TDE[n_step] = MBH_acc;
      if(MBH_teo>=MBH_f){
      Gal[p].BlHo[BH_indx].M_res_TDE[n_step] = 0.0;
      Gal[p].BlHo[BH_indx].MBH_acc_TDE[n_step] = MBH_f - MBH_i;// to guarantee mass conservation
      }
      Mres=Gal[p].BlHo[BH_indx].M_res_TDE[n_step];
      }
      #endif

     float arr_acc[8] = {Gal[p].HaloNr,Gal[p].HeapIndex,IND,MBH_acc,Mres,Mbh,Gal[p].BlHo[BH_indx].fEDD[n_step],eta};
     fileout = fopen(file_extBH,"ab");
     fwrite(arr_acc,sizeof(arr_acc),1,fileout);
     fclose(fileout);
} //end else
}
#endif

/*
###########################################################################
###################### End growth functions ##############################
###########################################################################
*/


void restart_the_BH(int p, int BH_indx, int nstep)
{
   // This function inits all the variables of a BH
   // p --> Galaxy id
   // index_BH --> if is the main BH of the binary system (i = 0) or is the satellite (i = 1)
   // nstep --> is the starting step from which we are going to (re)init everything
   int j;
   for(j=nstep; j<STEPS; j++)
      {
        Gal[p].BlHo[BH_indx].BlackHoleMass[j] = 1e-15; Gal[p].BlHo[BH_indx].Mass_to_acc[j] = 1e-15;
        Gal[p].BlHo[BH_indx].QuasarLum[j] = 0.0;       Gal[p].BlHo[BH_indx].spin[j]=0.0;
        Gal[p].BlHo[BH_indx].eta_BH[j] = 0.0;          Gal[p].BlHo[BH_indx].fEDD[j] = 0.0;

        #ifdef FEED_ON_STARS
         Gal[p].BlHo[BH_indx].M_dot_TDE[j] = 0.0;
       //Gal[p].BlHo[BH_indx].T_relax[j] = 0.0;     Gal[p].BlHo[BH_indx].M_core[j] = 0.0;
         #endif

        #ifdef RADIOLUMBHS
        Gal[p].BlHo[BH_indx].JetLum_TD[j] = 0.0;  Gal[p].BlHo[BH_indx].JetLum_ADAF[j] = 0.0;
        Gal[p].BlHo[BH_indx].RadioLum_TD[j] = 0.0; Gal[p].BlHo[BH_indx].RadioLum_ADAF[j] = 0.0;
        #endif

        #ifdef EXTENDED_BH_OUTPUT
        Gal[p].BlHo[BH_indx].MBH_acc_first_phase[j]  = 0.;   Gal[p].BlHo[BH_indx].MBH_acc_second_phase[j]  = 0.;
        Gal[p].BlHo[BH_indx].MBH_loss_GW[j]  = 0.0;     Gal[p].BlHo[BH_indx].MBH_acc_radiomode[j]  = 0.0; 
        Gal[p].BlHo[BH_indx].M_acc_Mergers[j]  = 0.;      Gal[p].BlHo[BH_indx].M_acc_DI_sec[j]  = 0.;
        Gal[p].BlHo[BH_indx].M_acc_DI_mmi[j]  = 0.;      Gal[p].BlHo[BH_indx].MBH_acc_Mergers[j]  = 0.; 
        Gal[p].BlHo[BH_indx].MBH_acc_DI_sec[j]  = 0.;     Gal[p].BlHo[BH_indx].MBH_acc_DI_mmi[j]  = 0.;
        Gal[p].BlHo[BH_indx].M_res_Mergers[j]  = 0.0;     Gal[p].BlHo[BH_indx].M_res_DI_sec[j]  = 0.0; 
        Gal[p].BlHo[BH_indx].M_res_DI_mmi[j]  = 0.0; 
        #ifdef FEED_ON_STARS
        Gal[p].BlHo[BH_indx].MBH_acc_TDE[j] = 0.0;   Gal[p].BlHo[BH_indx].M_res_TDE[j] = 0.0; Gal[p].BlHo[BH_indx].M_acc_TDE[j]  = 0.;   
        #endif 
        #endif
       }
        //Gal[p].BlHo[BH_indx].M_res_0 = 0.;               
        Gal[p].BlHo[BH_indx].MBH_0 = 0.;
        //Gal[p].BlHo[BH_indx].MBH_peak = 0.;              
        Gal[p].BlHo[BH_indx].MBH_second_phase = 0.;
        Gal[p].BlHo[BH_indx].M_dot_acc = 0.0; 
        Gal[p].BlHo[BH_indx].t_acc_quiescent_phase = -1.;
        Gal[p].BlHo[BH_indx].tBH_merger_clock_dt = -1.0; 
        Gal[p].BlHo[BH_indx].tBH_merger = -1.0;
        Gal[p].BlHo[BH_indx].t_spent_WPHs = -1.0;        
        Gal[p].BlHo[BH_indx].z_ejected_WPHs = -1.0;
        Gal[p].fgas = -1.0;                              
        Gal[p].z_binary_form = -1.0;                     
        Gal[p].t_binary_form = -1.0;
        Gal[p].BlHo[BH_indx].MBH_WPH = -1.0;
        #ifdef FEED_ON_STARS
        Gal[p].BlHo[BH_indx].t_TD_phase = -1.0;  
        #endif

}


void init_nuclear_BHs(int p)
{
  // This function inits all the nuclear black holes allowed inside of a galaxy
  // p --> is the id of the galaxy
  int i, j;
  float ao;
  #if defined MRII && defined BASICSEEDING
  float seedMin, seedMax, seedMass;

  if(ZZ[Gal[p].SnapNum]>=12.0)
    {float Ocupation;
    float Probability, ProbabilityTh, GQdGraftingProbability;
        seedMin = BlackHoleSeedMass;
        float seedMax = 1.0e4;
        if (BlackHoleSeedMass>seedMax) seedMax = BlackHoleSeedMass;
        seedMass = seedMin*pow(seedMax/seedMin,((double)rand()/(double)RAND_MAX));
       //Ocupation = gsl_rng_uniform(random_generator);
       //if(Ocupation>0.5) {Gal[p].BlHo[i].BlackHoleMass[j] = BlackHoleSeedMass*(Hubble_h*1.e-10);}
       // NOW WE COMPUTE THE OCUPATION PROBABILITY
       ProbabilityTh = (float)rand()/((float)RAND_MAX);
       GQdGraftingProbability = 0.45;
       Probability = GQdGraftingProbability * (Gal[p].Mvir*1.e10/Hubble_h)/1e8;
       if(Probability>=1.) {Probability = 1.0;}
       if(Probability<ProbabilityTh) {seedMass = 1e-15*1e10/Hubble_h;}
    }
  else {seedMass = 1e-15*1e10/Hubble_h;}
  #endif
  
  for(i=0;i<2;i++) 
     {
         ao = gsl_rng_uniform(random_generator); //the spin value is random
        for(j=0;j<STEPS;j++) // loop in the substeps
          {
             if(i==0) // This properties should only be added in the BH that is the main black hole. The galaxy does not have yet a Binary
                {
                 #if defined MRII && defined BASICSEEDING
                Gal[p].BlHo[i].BlackHoleMass[j] = seedMass*(Hubble_h*1.e-10); // internal units
                 #else
                 Gal[p].BlHo[i].BlackHoleMass[j] = BlackHoleSeedMass*(Hubble_h*1.e-10); // internal units
                  #endif
                  Gal[p].BlHo[i].spin[j] = ao;
                  if(Gal[p].BlHo[i].spin[j]>=0.998) Gal[p].BlHo[i].spin[j] = 0.998;
                  Gal[p].BlHo[i].eta_BH[j] = compute_eta(Gal[p].BlHo[i].spin[j]);
                  #ifndef SPINCOMPUTATION
                  Gal[p].BlHo[i].spin[j] = 0.0;
                  Gal[p].BlHo[i].eta_BH[j] = 0.1;
                  #endif
                }
             else
                 { Gal[p].BlHo[i].BlackHoleMass[j] = 1e-15;
                   Gal[p].BlHo[i].spin[j] = 0.0;
                   Gal[p].BlHo[i].eta_BH[j] = 0.1;
                 }
             Gal[p].BlHo[i].Mass_to_acc[j] = 1e-15;
             Gal[p].BlHo[i].QuasarLum[j] = 0.0;
             Gal[p].BlHo[i].fEDD[j] = 0.0;
             #ifdef RADIOLUMBHs
             Gal[p].BlHo[i].JetLum_TD[j] = 0.0;
             Gal[p].BlHo[i].JetLum_ADAF[j] =0.0;
             Gal[p].BlHo[i].RadioLum_TD[j] = 0.0;
             Gal[p].BlHo[i].RadioLum_ADAF[j] =  0.0;
             #endif
            #ifdef FEED_ON_STARS
            //Gal[p].BlHo[i].M_core[j] = 0.0;
            //Gal[p].BlHo[i].T_relax[j] = 0.0;
            Gal[p].BlHo[i].M_dot_TDE[j] = 0.0;
            #endif
            #ifdef EXTENDED_BH_OUTPUT
            Gal[p].BlHo[i].MBH_acc_first_phase[j]  = 0.0;
            Gal[p].BlHo[i].MBH_acc_second_phase[j]  = 0.0;
            Gal[p].BlHo[i].MBH_acc_radiomode[j]  = 0.0;  
            Gal[p].BlHo[i].MBH_loss_GW[j]  = 0.0;  
            Gal[p].BlHo[i].MBH_acc_Mergers[j]  = 0.0; 
            Gal[p].BlHo[i].MBH_acc_DI_sec[j]  = 0.0;
            Gal[p].BlHo[i].MBH_acc_DI_mmi[j]  = 0.0;
            Gal[p].BlHo[i].M_acc_Mergers[j]  = 0.0; 
            Gal[p].BlHo[i].M_acc_DI_sec[j]  = 0.0;
            Gal[p].BlHo[i].M_acc_DI_mmi[j]  = 0.0;
            Gal[p].BlHo[i].M_res_Mergers[j]  =  0.0; 
            Gal[p].BlHo[i].M_res_DI_sec[j]  =   0.0; 
            Gal[p].BlHo[i].M_res_DI_mmi[j]  =  0.0; 
            #ifdef FEED_ON_STARS
            Gal[p].BlHo[i].MBH_acc_TDE[j] = 0.0;
            Gal[p].BlHo[i].M_res_TDE[j] = 0.0;
            Gal[p].BlHo[i].M_acc_TDE[j] = 0.0;
            #endif 
            #endif
               } // end loop substeps
      //Gal[p].BlHo[i].M_res_0 = 0.; 
      Gal[p].BlHo[i].MBH_0 = 0.; 
      //Gal[p].BlHo[i].MBH_peak = 0.;
      Gal[p].BlHo[i].M_dot_acc = 0.0;
      Gal[p].BlHo[i].t_acc_quiescent_phase = -1.;
      Gal[p].BlHo[i].t_spent_WPHs = -1.0;  
      Gal[p].BlHo[i].z_ejected_WPHs = -1.0;
      Gal[p].BlHo[i].tBH_merger_clock_dt = -1.0; 
      Gal[p].BlHo[i].tBH_merger = -1.0;
      Gal[p].BlHo[i].MBH_WPH = -1.0;
      #ifdef FEED_ON_STARS
      Gal[p].BlHo[i].t_TD_phase = -1.0;
      #endif
     } // end loop in allowed BHs inside a galaxy
}

void update_BH_properties(int p)
{
  // The black holes in this version of L-Galaxies are treated as vectors of lenght STEPS
  // Therefore, when the loop in a STEPS of a given SNP is finished, we initialize the 
  // step = 0 of the new SNP (SNP+1) with the last values of steps of the SNP -1 
  int j,i ;
  //#ifdef FEED_ON_STARS //simplified
  //int mcore_neg=1;
  //if(Gal[p].BlHo[0].M_core[STEPS-1]>0.0) {mcore_neg=-1;}
  //#endif
  for(i=0;i<2;i++)
    {
       // SOME CHECKS IN THE SPIN COMPUTATION! 0 < spin < 1.0 
       if(fabs(Gal[p].BlHo[i].spin[STEPS-1])>=1.0)
         {
            printf("Error in update_BH_properties(int p), THE BH num %d the spin is > 1 (%f)\n",i,Gal[p].BlHo[i].spin[STEPS-1]);
            for(j=0;j<STEPS;j++)
                printf("For the BH num %d in the galaxy a[j] = %f, BlackHoleMass = %f\n",i,Gal[p].BlHo[i].spin[j],log10(Gal[p].BlHo[i].BlackHoleMass[j]/(Hubble_h*1.e-10)));
            exit(0);
          }
       // END CHECKS!
       for(j=0;j<STEPS;j++)
         {
             Gal[p].BlHo[i].BlackHoleMass[j] = Gal[p].BlHo[i].BlackHoleMass[STEPS-1];
             Gal[p].BlHo[i].Mass_to_acc[j] = Gal[p].BlHo[i].Mass_to_acc[STEPS-1];
             #ifdef EXTENDED_BH_OUTPUT
             Gal[p].BlHo[i].M_res_Mergers[j] = Gal[p].BlHo[i].M_res_Mergers[STEPS-1];
             Gal[p].BlHo[i].M_res_DI_sec[j] = Gal[p].BlHo[i].M_res_DI_sec[STEPS-1];
             Gal[p].BlHo[i].M_res_DI_mmi[j] = Gal[p].BlHo[i].M_res_DI_mmi[STEPS-1];
             #endif
             Gal[p].BlHo[i].QuasarLum[j] = Gal[p].BlHo[i].QuasarLum[STEPS-1];
             #ifdef RADIOLUMBHs
             Gal[p].BlHo[i].JetLum_TD[j] = Gal[p].BlHo[i].JetLum_TD[STEPS-1];
             Gal[p].BlHo[i].JetLum_ADAF[j] = Gal[p].BlHo[i].JetLum_ADAF[STEPS-1];
             Gal[p].BlHo[i].RadioLum_TD[j] = Gal[p].BlHo[i].RadioLum_TD[STEPS-1];
             Gal[p].BlHo[i].RadioLum_ADAF[j] = Gal[p].BlHo[i].RadioLum_ADAF[STEPS-1];
             #endif
             Gal[p].BlHo[i].fEDD[j] = Gal[p].BlHo[i].fEDD[STEPS-1];
             Gal[p].BlHo[i].spin[j] = Gal[p].BlHo[i].spin[STEPS-1];
             Gal[p].BlHo[i].eta_BH[j] = Gal[p].BlHo[i].eta_BH[STEPS-1];
            #ifdef FEED_ON_STARS
            // float dt_step;
            // dt_step =(NumToTime(Gal[p].SnapNum) - NumToTime(Gal[p].SnapNum+1))/STEPS;
            //if (mcore_neg<0){
            //Gal[p].BlHo[i].M_dot_TDE[j] =Gal[p].BlHo[i].M_dot_TDE[STEPS-1];
            //Gal[p].BlHo[i].T_relax[j] =Gal[p].BlHo[0].T_relax[STEPS-1];   
            //Gal[p].BlHo[i].BlackHoleMass[j] += (j+1)*dt_step*(Gal[p].BlHo[i].M_dot_TDE[STEPS-1]);
            //#ifdef EXTENDED_BH_OUTPUT
            //Gal[p].BlHo[i].MBH_acc_TDE[j] += dt_step*Gal[p].BlHo[i].M_dot_TDE[STEPS-1];//subtract also what was added previously, TDE_rates do not add up
            //#endif
            //Gal[p].BlHo[i].M_core[j] = Gal[p].BlHo[0].M_core[STEPS-1] -(j+1)*dt_step*(Gal[p].BlHo[i].M_dot_TDE[STEPS-1]/Stars2TDEefficiency);
            //if (Gal[p].BlHo[0].M_core[j]<0.0){
            //Gal[p].BlHo[i].BlackHoleMass[j] += Gal[p].BlHo[0].M_core[j]*Stars2TDEefficiency;
            //#ifdef EXTENDED_BH_OUTPUT
            //Gal[p].BlHo[i].MBH_acc_TDE[j] +=  Gal[p].BlHo[0].M_core[j]*Stars2TDEefficiency;//subtract also what was added previously, TDE_rates do not add up
            //#endif
            //Gal[p].BlHo[i].M_dot_TDE[j] =0.0;
            // Gal[p].BlHo[0].M_core[j] =0.0;
            // Gal[p].BlHo[0].T_relax[j] =0.0;
            // mcore_neg=1;
            //         }
            //}
            //  else{
           // Gal[p].BlHo[i].M_dot_TDE[j] =0.0;
            //  Gal[p].BlHo[0].M_core[j] =0.0;
            //  Gal[p].BlHo[0].T_relax[j] =0.0;
            //  }
            //simplified
            Gal[p].BlHo[i].M_res_TDE[j] = Gal[p].BlHo[i].M_res_TDE[STEPS-1];
            #endif
    }//end substeps
    }// end loop in the two BHs
}


void exchange_central_BH_positions(int p)
{
   // This function has just one porpuse, keeping always the BH[0] (i.e the main BH of the binary) as the most massive one
   // So, if we notice that BlackHoleMass[1] > BlacHoleMass[0] --> We need to exhange the black holes otherwise the code will complain  
   struct BH_struct BHs_aux;
   BHs_aux = Gal[p].BlHo[0];
   Gal[p].BlHo[0] = Gal[p].BlHo[1];
   Gal[p].BlHo[1] = BHs_aux;
   Gal[p].BlHo[0].tBH_merger_clock_dt = Gal[p].BlHo[1].tBH_merger_clock_dt;
   Gal[p].BlHo[0].tBH_merger = Gal[p].BlHo[1].tBH_merger;
   Gal[p].BlHo[1].tBH_merger_clock_dt = -1.0;
   Gal[p].BlHo[1].tBH_merger = -1.0;
}

float BH_Mass_accreted_via_merger(int p, double mass_ratio, int n_step, int filenr)
{
  /** @brief Grows black hole through accretion from cold gas during mergers,
   *         as in Kauffmann & Haehnelt (2000).
   *                              in main.c */
  // p --> index of the galaxy
  // mass_ratio --> baryonic mass ratio between the two interacting galaxies
  // n_step --> step in which the merger takes place
  // filenr --> Number of the merger-tree file that L-Galaxies is computing. 
  double zz;
  float dt, BHaccrete;

  if(Gal[p].ColdGas > 1e-9)
    {
        dt = ( NumToTime(Gal[p].SnapNum) -  NumToTime(Gal[p].SnapNum+1) ) / STEPS;
        zz = redshift_at_t(ZZ[Gal[p].SnapNum],n_step,dt);

        if(GassReservoirRedshiftDepend==1)
                BHaccrete = (1.0+ zz) * BlackHoleGrowthRate * mass_ratio / (1.0 + pow2((BlackHoleCutoffVelocity / Gal[p].Vvir))) * Gal[p].ColdGas;
        else if (GassReservoirRedshiftDepend==2)
                BHaccrete = pow(1.0+ zz,3.0/2.0) * BlackHoleGrowthRate * mass_ratio / (1.0 + pow2((BlackHoleCutoffVelocity / Gal[p].Vvir))) * Gal[p].ColdGas;
        else
                BHaccrete = BlackHoleGrowthRate * mass_ratio / (1.0 + pow2((BlackHoleCutoffVelocity / Gal[p].Vvir))) * Gal[p].ColdGas;

      //counter dMgas ,  Mass ratio, Mhalo, cold gas
   //    #ifdef EXTENDED_BH_OUTPUT

   //    char file_BHaccmerg[3000];
   //    FILE *fileout;
   //    //sprintf(file_extBH,"%s/Extended_BH",OutputDir);
   //    //fileout = fopen(file_extBH,"ab");
   //    //sprintf(file_BHaccmerg,"%s/BHacc_mergers_%d_z%1.2f_%d",OutputDir, n_step, ZZ[Gal[p].SnapNum],filenr);
   //    sprintf(file_BHaccmerg,"%s/BHaccmerg_%d",OutputDir,filenr);

   //   float arr_acc[9] = {zz,BHaccrete,Gal[p].BlHo[0].BlackHoleMass[n_step], Gal[p].BulgeMass+Gal[p].DiskMass, mass_ratio, Gal[p].ColdGas, Gal[p].Mvir, Gal[p].BlHo[0].fEDD[n_step], Gal[p].BlHo[0].eta_BH[n_step]};
   //   fileout = fopen(file_BHaccmerg,"ab");
   //   fwrite(arr_acc,sizeof(arr_acc),1,fileout);
   //   fclose(fileout);

   //    #endif

      /* cannot accrete more gas than is available! */
      if(BHaccrete > Gal[p].ColdGas){ 
      BHaccrete = Gal[p].ColdGas;
      }
      if(BHaccrete<1e-13) BHaccrete = 1e-13;


  }

 else {BHaccrete = 1e-13;}
      //printf("Mergers => BHaccrete = %E, mass_ratio = %E, n_step = %d, BlackHoleGrowthRate = %E, zz = %f\n",BHaccrete,mass_ratio, n_step,BlackHoleGrowthRate,zz);
      //exit(0);

 return  BHaccrete;
}


float BH_Mass_accreted_via_DI(int p, float fraction, float stars, int step, int filenr)
{
  /** @brief Grows black hole through accretion from cold gas during DI,
   *         as in Izquierdo-Villalba et al. (2020) */
  // p --> index of the galaxy
  // fraction --> baryonic mass ratio between the two interacting galaxies
  // stars --> Amount of stars that trigger the disk instability 
  // step --> step in which the merger takes place
  // filenr --> Number of the merger-tree file that L-Galaxies is computing. 

  float dt_step, BHaccrete;
  double zz;


  if(Gal[p].ColdGas > 1e-9)
    {
        dt_step =  ( NumToTime(Gal[p].SnapNum) -  NumToTime(Gal[p].SnapNum+1) ) / STEPS;
        zz = redshift_at_t(ZZ[Gal[p].SnapNum],step,dt_step);
        if( Gal[p].induced > 0 )
           {
             if(GassReservoirRedshiftDepend==1) BHaccrete = (1.0+zz) * BlackHoleGrowthRate / (1.0 + pow2((BlackHoleCutoffVelocity / Gal[p].Vvir))) * stars;
             else if (GassReservoirRedshiftDepend==2) BHaccrete = pow(1.0+zz,3.0/2.0) * BlackHoleGrowthRate / (1.0 + pow2((BlackHoleCutoffVelocity / Gal[p].Vvir))) *stars;
             else BHaccrete = BlackHoleGrowthRate / (1.0 + pow2((BlackHoleCutoffVelocity / Gal[p].Vvir))) * stars;
           }

        else
           {
             if(GassReservoirRedshiftDepend==1) BHaccrete = (1.0+zz) * BlackHoleGrowthRateDI / (1.0 + pow2((BlackHoleCutoffVelocity / Gal[p].Vvir))) * stars;
             else if (GassReservoirRedshiftDepend==2) BHaccrete = pow(1.0+zz,3.0/2.0) * BlackHoleGrowthRateDI / (1.0 + pow2((BlackHoleCutoffVelocity / Gal[p].Vvir))) * stars;
             else BHaccrete = BlackHoleGrowthRateDI / (1.0 + pow2((BlackHoleCutoffVelocity / Gal[p].Vvir)))* stars;
           }

   //    #ifdef EXTENDED_BH_OUTPUT

   //    char file_BHaccdi[3000];
   //    FILE *fileout;
   //    //sprintf(file_extBH,"%s/Extended_BH",OutputDir);
   //    //fileout = fopen(file_extBH,"ab");
   //    //sprintf(file_BHaccmerg,"%s/BHacc_mergers_%d_z%1.2f_%d",OutputDir, n_step, ZZ[Gal[p].SnapNum],filenr);
   //    sprintf(file_BHaccdi,"%s/BHaccdi_%d",OutputDir,filenr);

   //   float arr_acc[9] = {zz,BHaccrete,Gal[p].BlHo[0].BlackHoleMass[step], Gal[p].BulgeMass+Gal[p].DiskMass, stars, Gal[p].ColdGas, Gal[p].Mvir, Gal[p].BlHo[0].fEDD[step],Gal[p].BlHo[0].eta_BH[step]};
   //   fileout = fopen(file_BHaccdi,"ab");
   //   fwrite(arr_acc,sizeof(arr_acc),1,fileout);
   //   fclose(fileout);

   //    #endif

        //printf("DI => BHaccrete = %E\n",BHaccrete);
        if(BlackHoleAccretionModel==0) BHaccrete = Gal[p].ColdGas*fraction;
        if(BHaccrete > Gal[p].ColdGas) {
      char file_BHaccdi[3000];
      FILE *fileout;
      //sprintf(file_extBH,"%s/Extended_BH",OutputDir);
      //fileout = fopen(file_extBH,"ab");
      //sprintf(file_BHaccmerg,"%s/BHacc_mergers_%d_z%1.2f_%d",OutputDir, n_step, ZZ[Gal[p].SnapNum],filenr);
         sprintf(file_BHaccdi,"%s/BHexcessdi_%d",OutputDir,filenr);
         float arr_acc[9] = {zz,BHaccrete,Gal[p].BlHo[0].BlackHoleMass[step], Gal[p].BulgeMass+Gal[p].DiskMass, stars, Gal[p].ColdGas, Gal[p].Mvir,Gal[p].BlHo[0].fEDD[step], Gal[p].BlHo[0].eta_BH[step]};
         fileout = fopen(file_BHaccdi,"ab");
         fwrite(arr_acc,sizeof(arr_acc),1,fileout);
         fclose(fileout);
          BHaccrete = Gal[p].ColdGas;
        }
        if(BHaccrete<1e-13) BHaccrete = 1e-13;
        if(BHaccrete<0) printf("Check in DI the value of BHaccrete = %f\n",BHaccrete);
        //if(BHaccrete<1e-15) BHaccrete = 1e-16;
     }//end Gal[p].ColdGas
  else {BHaccrete = 1e-13;}

  return BHaccrete;
}




float compute_M_dot_Edd(float MBH, float epsilon)
{ 
  // This function computes the M_dot_Edd of an accreting BH in INTERNAL UNITS!!
  float c, M_Edd_dot;
  c = 2.99e8;//speed of light [m/s]
  M_Edd_dot = 6.28*(UnitTime_in_s)*MBH/(epsilon*c*c*Hubble_h); // [internal units of mass and time]
  return M_Edd_dot;

}


#ifdef RADIOLUMBHS
void Compute_Radio_Luminosities(int p, int nstep, int BH_indx, float MBH, float M_dot, float a, float fEdd, float eta)
{
  // This function computes the radio luminosities follwoing a similar approach than Fanidakis et al. 201X
  float M_Edd_dot;
  M_Edd_dot = compute_M_dot_Edd(MBH,eta);
 
   if(fEdd>FEDD_DOT_CRIT) 
      {
         // Thin diks geometry
         Gal[p].BlHo[BH_indx].JetLum_TD[nstep] = pow(10.,2.7)*( pow(0.1/0.01,-0.1)*pow(MBH*10/Hubble_h,0.9)*pow(M_dot/(0.1*M_Edd_dot),6./5.) * (1.0+(1.1*a)+(0.29*a*a))); // CAUTION in units of 10^40 erg/s
         Gal[p].BlHo[BH_indx].RadioLum_TD[nstep] = 1e-4*pow(M_dot/(0.1*M_Edd_dot),-1.2)*pow(MBH*1e10/Hubble_h,0.32)*Gal[p].BlHo[BH_indx].JetLum_TD[nstep];
         Gal[p].BlHo[BH_indx].JetLum_ADAF[nstep] = 0.0; Gal[p].BlHo[BH_indx].RadioLum_ADAF[nstep] = 0.0;
         //printf("HOT Thin disk --> JetLum_TD = %f, RadioLum_TD = %f\n", 40.0+log10(Gal[p].BlHo[BH_indx].JetLum_TD[nstep]), 40.0+log10(Gal[p].BlHo[BH_indx].RadioLum_TD[nstep]));
      }
   else 
      {
          // ADAF geometry
          Gal[p].BlHo[BH_indx].JetLum_ADAF[nstep] = pow(10.,5.1)*( pow(0.1/0.3,-1.)*(MBH*10/Hubble_h)*(M_dot/(0.1*M_Edd_dot)) * pow(2.3,2.0) * (0.55+(1.5*a)+(a*a))); // CAUTION in units of 10^40 erg/s
          Gal[p].BlHo[BH_indx].RadioLum_ADAF[nstep] = 2e-8*pow((M_dot/(0.1*M_Edd_dot))*MBH*1e10/Hubble_h,0.42)*Gal[p].BlHo[BH_indx].JetLum_ADAF[nstep];
          Gal[p].BlHo[BH_indx].JetLum_TD[nstep] = 0.0; Gal[p].BlHo[BH_indx].RadioLum_TD[nstep] = 0.0;
          //printf("HOT ADAF --> P = %f\n",log10(Gal[p].BlHo[BH_indx].RadioLum_ADAF[nstep])+40.-7-log10(1.4*1e9));
          //exit(0);
      }

   if(fEdd <= 1e-11 || fEdd == 0.) {Gal[p].BlHo[BH_indx].JetLum_TD[nstep] = 0.; Gal[p].BlHo[BH_indx].RadioLum_TD[nstep] = 0.; Gal[p].BlHo[BH_indx].JetLum_ADAF[nstep] = 0.; Gal[p].BlHo[BH_indx].RadioLum_ADAF[nstep] = 0.;}
}
#endif

void Determine_V_Sigma_for_BH_after_galaxy_interaction(int t, int p, float m_ratio, float f_central, float f_satellite, float tmass, float pmass, float tbulge_mass, int n_step, int Smooth_accretion)
{
  // In this function I determine what is the v_sigma that will determine the spin evolution of the BH
  //
  // t central, p satellite
  
  if(m_ratio >= ThreshMajorMerger) // Major merger
        {
           // Satellite galaxy
           Gal[p].v_sigma_now = -1;
           Gal[p].type_v_sigma = -1;
           // Central galaxy
           Gal[t].type_v_sigma = 0;
           Gal[t].v_sigma_now = altura_number(Gal[t].type_v_sigma,t,n_step);
           if(isnan(Gal[t].v_sigma_now) == 1) {printf("Error in add_galaxies_together(). The value of central v_sigma_now is nan!\n"); exit(0);}

        }
  else //Minor merger Central galaxy
        {
           int zzz = -1;
           // Both of them do not have bulge
           if(Gal[t].v_sigma_now < 0 && Gal[p].v_sigma_now < 0)
              {
                  if(Smooth_accretion>0){Gal[t].type_v_sigma = 1; Gal[t].v_sigma_now = altura_number(2,t,n_step);}
                  else{Gal[t].type_v_sigma = 1; Gal[t].v_sigma_now = altura_number(Gal[t].type_v_sigma,t,n_step);}
               }
           // The central galaxy has a bulge but the other one no
           else if(Gal[t].v_sigma_now > 0 && Gal[p].v_sigma_now < 0)
               {
                   if(Smooth_accretion>0){Gal[t].type_v_sigma = 1; Gal[t].v_sigma_now =  Gal[t].v_sigma_now;}
                   else {Gal[t].type_v_sigma = 1; Gal[t].v_sigma_now = altura_number(Gal[t].type_v_sigma,t,n_step);}
               }
           // If the satellite galaxy has a bulge but the other one no
           else if(Gal[t].v_sigma_now <0 && Gal[p].v_sigma_now > 0)
               {       
		   if(Smooth_accretion>0){Gal[t].type_v_sigma = 1; Gal[t].v_sigma_now =  Gal[p].v_sigma_now;}
                   else {Gal[t].type_v_sigma = 1; Gal[t].v_sigma_now = altura_number(Gal[p].type_v_sigma,t,n_step);}
               }
           // The central and satellite galaxy have a bulge
           else if(Gal[t].v_sigma_now >0 && Gal[p].v_sigma_now > 0)
               {  Gal[t].type_v_sigma = 1;
                  if(Smooth_accretion>0)
                      {Gal[t].v_sigma_now = ( ((f_central + tmass)* Gal[t].v_sigma_now) + ( (f_satellite+pmass)*Gal[p].v_sigma_now) ) / (pmass + tmass + f_central+f_satellite); zzz = +1;}
                  else
                      {   Gal[t].type_v_sigma = 1;
                          zzz = +1;
                          if((tbulge_mass + pmass)<=2e-9) Gal[t].v_sigma_now = altura_number(Gal[t].type_v_sigma,t,n_step);
                          else Gal[t].v_sigma_now = ( (tbulge_mass* Gal[t].v_sigma_now) + (pmass*Gal[p].v_sigma_now) ) / (tbulge_mass+pmass);
                       }

                }
           // of not the case.... be carefull, smoething is wrong!
           else {printf("Error in add galaxies together v_sigma in minor merger, Gal[p].type_v_sigma = %d, Gal[t].type_v_sigma = %d\n",Gal[p].type_v_sigma,Gal[t].type_v_sigma); exit(0);}

           if(Gal[t].v_sigma_now < 0) printf("Here is the error in Determine_V_Sigma_for_BH_after_galaxy_interaction, zzz = %d\n",zzz);
           if(isnan(Gal[t].v_sigma_now) == 1) {printf("Error in add_galaxies_together(). The value of central v_sigma_now is nan!\n"); exit(0);}
           // Satellite galaxy
           Gal[p].v_sigma_now = -1;
           Gal[p].type_v_sigma = -1;
        }

  //printf("n_step = %d, m_ratio = %f, f_central = %f, tmass = %f, f_satellite = %f, pmass = %f, tbulge_mass = %f, Gal[t].v_sigma_now = %f\n",n_step,m_ratio,f_central,tmass,f_satellite,pmass,tbulge_mass,Gal[t].v_sigma_now);
  //exit(0);
}


void order_and_Exchange_Black_Holes_By_Mass(int t, int p, int n_step, int order[4], int *exchange, int filenr)
{

 // This function puts the two most massive black holes as the ones that will forme a binary black hole
 float first, second;
 int i, idx_first, idx_second, idx_third, idx_forth;
 first = second = -1.0;
 idx_first = idx_second = idx_third = idx_forth = -1;
 struct BH_struct BlHo_aux[4];
 BlHo_aux[0] = Gal[t].BlHo[0];
 BlHo_aux[1] = Gal[t].BlHo[1];
 BlHo_aux[2] = Gal[p].BlHo[0];
 BlHo_aux[3] = Gal[p].BlHo[1];

 // I look for the first and second most massive BHs
 // Notice that the way that I use here, if two BHs of coming from different galaxies have the same maes
 //  the one that will be firstly saved will be always the one of the central galaxy
 for(i = 0; i<4; i++)
   {
      if (BlHo_aux[i].BlackHoleMass[n_step] > first)
	  {second = first; first = BlHo_aux[i].BlackHoleMass[n_step]; idx_second = idx_first; idx_first = i;}
      else if (BlHo_aux[i].BlackHoleMass[n_step] > second && BlHo_aux[i].BlackHoleMass[n_step] <= first)
	  {second = BlHo_aux[i].BlackHoleMass[n_step]; idx_second = i;}
   }

  if(idx_second<0 || idx_first<0) {printf("Error: Something weird in add_and_update_BHs_mergering_galaxies(): idx_second and idx_first should be > 0 and now are idx_first = %d idx_second = %d\n",idx_first,idx_second); exit(0);}


  if(idx_first > 1 ) *exchange = +1;
  else if(*exchange < 0 &&  idx_second > 1) *exchange = +1; // Here we consider that if the the primary galaxy has not a 2nd BH, then there is an exchange 
  //else if(*exchange < 0 &&  idx_second > 1 && BlHo_aux[i].BlackHoleMass[n_step] > 1e-14) *exchange = +1; // Here we consider that if the the primary galaxy has not a 2nd BH, then there is not an exchange 
  else *exchange = -1;


  order[0] = idx_first;
  order[1] = idx_second;
  int order_aux[2], kk = 0;
  for(i=0;i<4;i++) if(i!=order[0] && i!=order[1]) {order_aux[kk] = i; kk++;}
  if (BlHo_aux[order_aux[0]].BlackHoleMass[n_step] >= BlHo_aux[order_aux[1]].BlackHoleMass[n_step]) {idx_third = order_aux[0]; idx_forth = order_aux[1];}
  else {idx_third = order_aux[1]; idx_forth = order_aux[0];}

  order[2] = idx_third;
  order[3] = idx_forth;
  Gal[t].BlHo[0] = BlHo_aux[idx_first];
  Gal[t].BlHo[1] = BlHo_aux[idx_second];
  Gal[p].BlHo[0] = BlHo_aux[idx_third];
  Gal[p].BlHo[1] = BlHo_aux[idx_forth]; 
	  
  // SOME CHECKS
  if(Gal[t].BlHo[0].BlackHoleMass[n_step] < 1e-14 && Gal[t].BlHo[1].BlackHoleMass[n_step]<1e-14) {printf("Error: Something weird in add_and_update_BHs_mergering_galaxies(): In the central galaxy: BH0<1e-14((%E), BH1<1e-14(%E)\n",Gal[t].BlHo[0].BlackHoleMass[n_step],Gal[t].BlHo[1].BlackHoleMass[n_step]); exit(0);}
  if(Gal[t].BlHo[0].BlackHoleMass[n_step]<Gal[t].BlHo[1].BlackHoleMass[n_step]) {printf("In add_and_update_BHs_mergering_galaxies of the BH function the BH0<BH1 (file = %d)\n",filenr); exit(0);}
}

void add_and_update_BHs_mergering_galaxies(int t, int p, int n_step, float m_ratio, int Smooth_accretion, float pmass, float tmass, float f_central, float f_satellite, float pbulge_mass, float tbulge_mass, int filenr)
{
  // t central, p satellite
  int j,i;

  Determine_V_Sigma_for_BH_after_galaxy_interaction(t,p,m_ratio,f_central, f_satellite, tmass, pmass, tbulge_mass, n_step,Smooth_accretion); //  What is the type of v/sigma corresponding to this process
 
  // Here I join the wandering black holes of the two galaxies: I transfer all the wandering BHs of the satellite to the central galaxy 
  // Given that the merger happens between Type2 and Type0 (most of of the cases), the transfer should have happened when the satellite galaxy
  // becamed from Type0/Type1 to Type2
  
  #ifdef MODELWANDERINGBHS
  join_wandering_BHs(t,p,filenr);
  #endif
  
  if(BlackHoleAccretionModel == 0) Merger_BH_as_Henriques(t,p,n_step);
  else
      {
        if( (Gal[t].BlHo[0].BlackHoleMass[n_step] + Gal[t].BlHo[1].BlackHoleMass[n_step] + Gal[p].BlHo[0].BlackHoleMass[n_step] + Gal[p].BlHo[1].BlackHoleMass[n_step]) >= 1e-14)
         {
            int order[4], exchange;
	    order_and_Exchange_Black_Holes_By_Mass(t, p, n_step, order, &exchange, filenr); // I order the possible 4 BHs by mass (notice that it can be just two black holes)
	    for(j=n_step; j<STEPS; j++) {
       Gal[t].BlHo[0].Mass_to_acc[j] += (Gal[p].BlHo[0].Mass_to_acc[j] + Gal[p].BlHo[1].Mass_to_acc[j]);
       #ifdef FEED_ON_STARS
       //Gal[t].BlHo[0].M_core[j] =  0.0; //Gal[p].BlHo[0].M_core[j]; disruption during mergers, if new bulge satisfies conditions will form new 
       #ifdef TIME_DEPENDENT_TDE
         Gal[p].BlHo[0].t_TD_phase  = -1.0; //simplified
       #endif
       #endif
       }
            // I add in the reservoir of the most massive one the reservoir of the rest of BHs except the second one
	    /////////for(j = n_step; j<STEPS; j++) Gal[t].BlHo[0].Mass_to_acc[j] += Gal[t].BlHo[0].Mass_to_acc[j] + Gal[t].BlHo[1].Mass_to_acc[j];
            // I reset the satellite BHS  (galaxy p)
            restart_the_BH(p,0,n_step);
            restart_the_BH(p,1,n_step);
            // I check what is the delay in the BH-BH merger. Notice that I nedd that BH1 > 0 and BH2>0, i.e there are two BHs to merge!
	    // Notice as well that I'm computing the delay as in Izquierdo-Villalba et al. 2020, future versions will deal with the BBH in a more sophisticated way
            if(Gal[t].BlHo[0].BlackHoleMass[n_step] > 1e-14 && Gal[t].BlHo[1].BlackHoleMass[n_step] > 1e-14)
              {
                 float dt, BH_accreted_by_merger;
                 dt = (NumToTime(Gal[t].SnapNum) - NumToTime(Gal[t].SnapNum + 1))/STEPS;
                 compute_BH_BH_merger_delay(t,n_step,filenr);
		 if(Gal[t].BlHo[0].tBH_merger_clock_dt < dt)
		   {
		      BH_accreted_by_merger = BH_Mass_accreted_via_merger(t,m_ratio, n_step, filenr); // Note that I do not count this twice eh!
                      Gal[t].BlHo[0].Mass_to_acc[n_step] += BH_accreted_by_merger; // TODO: NOTE THIS AND THE NEXT TWO LINES
		      merger_BHs(t,n_step,filenr);
                      if(Gal[t].BlHo[0].BlackHoleMass[n_step]<1e-9)  Gal[t].ColdGas -= BH_accreted_by_merger; // This is because in add_BH_in_halo I transfered the BH reservoir to the Cold gas
		      else  Gal[t].BlHo[0].Mass_to_acc[n_step] -= BH_accreted_by_merger; // TODO: NOTE THAT HERE I AM REMOVINF THE GAS THAT LATTER WILL BE ADDED
                   }
              } // end if there are two BHs to merge!

          } // end if MBH1 + MBH2 + MBH3 + MBH4 > 0
        else  for(i=0;i<2;i++)  {restart_the_BH(p,i,n_step); restart_the_BH(t,i,n_step);}
      }// end else Accretion Model
}


void merger_BHs(int p, int n_step, int filenr)
{

  // This function deals with the properties of the remament black hole after the merger	
  int j,i;
  float M_bh_f_GW, a_fin,V_recoil;
  int index_BH_ejected, is_ejected;
  float Mdisk_ejection, qq, BH1, BH2, a1, a2;
  float fgas_BHs, Mres_tot, MBH_tot;
  float BHf, af;
  index_BH_ejected = -1;
  is_ejected = -1;
  Mdisk_ejection = 0.0;
  
//   #ifdef FEED_ON_STARS
//   if(Gal[p].BlHo[1].BlackHoleMass[n_step] >  Gal[p].BlHo[0].BlackHoleMass[n_step])
//   {
//    struct BH_struct BHs_aux;
//    BHs_aux = Gal[p].BlHo[0];
//    Gal[p].BlHo[0] = Gal[p].BlHo[1];
//    Gal[p].BlHo[1] = BHs_aux;
//    Gal[p].BlHo[0].tBH_merger_clock_dt = Gal[p].BlHo[1].tBH_merger_clock_dt;
//    Gal[p].BlHo[0].tBH_merger = Gal[p].BlHo[1].tBH_merger;
//    Gal[p].BlHo[1].tBH_merger_clock_dt = -1.0;
//    Gal[p].BlHo[1].tBH_merger = -1.0;
//   }
//   #endif
  qq = Gal[p].BlHo[1].BlackHoleMass[n_step] / Gal[p].BlHo[0].BlackHoleMass[n_step]; // Remember that the BH 0 is always more massive than the BH 1
  
  //Re-set this two values always of tBH_merger and tBH_merger_clock_dt here!
  Gal[p].BlHo[0].tBH_merger = -1.0;
  Gal[p].BlHo[0].tBH_merger_clock_dt = -1.0;
  Gal[p].BlHo[1].tBH_merger = -1.0;
  Gal[p].BlHo[1].tBH_merger_clock_dt = -1.0;

  //In case the qq value is > one I show the issue and I exit the code
  if(qq>1.) {printf("Error in merger_BHs: File = %d qq>1 (%f). In the step =%d, BH0 = %E, BH1 = %E \n",filenr,qq,n_step,Gal[p].BlHo[0].BlackHoleMass[n_step],Gal[p].BlHo[1].BlackHoleMass[n_step]);
             for(i=0;i<STEPS;i++) printf("step = %d, BH0 = %E, BH1 = %E\n",i,Gal[p].BlHo[0].BlackHoleMass[i],Gal[p].BlHo[1].BlackHoleMass[i]); exit(0);}

  BH1 = Gal[p].BlHo[0].BlackHoleMass[n_step];
  BH2 = Gal[p].BlHo[1].BlackHoleMass[n_step];
  a1 = Gal[p].BlHo[0].spin[n_step];
  a2 = Gal[p].BlHo[1].spin[n_step];
  Mres_tot = Gal[p].BlHo[0].Mass_to_acc[n_step] + Gal[p].BlHo[1].Mass_to_acc[n_step];
  MBH_tot  = Gal[p].BlHo[0].BlackHoleMass[n_step] + Gal[p].BlHo[1].BlackHoleMass[n_step];
  fgas_BHs = (Mres_tot)/(Mres_tot+MBH_tot);

  #ifdef SPINCOMPUTATION
  BH_BH_spin_mass(p, n_step, &a_fin, &M_bh_f_GW,&V_recoil,filenr); // Here I compute the spin of the final black hole
  //printf("BH1 = %E, BH2 = %E, a1 = %f, a2 = %f, M_bh_f_GW = %E\n",BH2,BH1,a1,a2,M_bh_f_GW);
  if(isnan(M_bh_f_GW) ==1) {printf("Error in M_bh_f_GW is nan, check it.\n"); exit(0);}
#ifdef EXTENDED_BH_OUTPUT
save_extended_BHproperties(p, n_step, 0, filenr, 4, M_bh_f_GW, MBH_tot, MBH_tot, MBH_tot - M_bh_f_GW, 0.0);
#endif
// add together accretion masses from each timestep for both black holes before they merge
#ifdef EXTENDED_BH_OUTPUT
Gal[p].BlHo[0].MBH_loss_GW[n_step] += MBH_tot - M_bh_f_GW;  
for(j=0; j<n_step; j++)
{    
     //Gal[p].BlHo[0].M_acc_Mergers[j] += Gal[p].BlHo[1].M_acc_Mergers[j];
     //Gal[p].BlHo[0].M_acc_DI_sec[j] += Gal[p].BlHo[1].M_acc_DI_sec[j];
     //Gal[p].BlHo[0].M_acc_DI_mmi[j] += Gal[p].BlHo[1].M_acc_DI_mmi[j];
     Gal[p].BlHo[0].MBH_acc_Mergers[j] += Gal[p].BlHo[1].MBH_acc_Mergers[j];
     Gal[p].BlHo[0].MBH_acc_DI_sec[j] += Gal[p].BlHo[1].MBH_acc_DI_sec[j];
     Gal[p].BlHo[0].MBH_acc_DI_mmi[j] += Gal[p].BlHo[1].MBH_acc_DI_mmi[j];
     Gal[p].BlHo[0].M_res_Mergers[j] += Gal[p].BlHo[1].M_res_Mergers[j];
     Gal[p].BlHo[0].M_res_DI_sec[j] += Gal[p].BlHo[1].M_res_DI_sec[j];
     Gal[p].BlHo[0].M_res_DI_mmi[j] += Gal[p].BlHo[1].M_res_DI_mmi[j];
     Gal[p].BlHo[0].MBH_acc_first_phase[j] += Gal[p].BlHo[1].MBH_acc_first_phase[j];
     Gal[p].BlHo[0].MBH_acc_second_phase[j] += Gal[p].BlHo[1].MBH_acc_second_phase[j];
      #ifdef FEED_ON_STARS
      Gal[p].BlHo[0].MBH_acc_TDE[j] += Gal[p].BlHo[1].MBH_acc_TDE[j];
      Gal[p].BlHo[0].M_res_TDE[j] += Gal[p].BlHo[1].M_res_TDE[j];
     #endif
}
#endif
  for(j = n_step; j<STEPS; j++)
   {
     Gal[p].BlHo[0].BlackHoleMass[j] = M_bh_f_GW;
     Gal[p].BlHo[0].Mass_to_acc[j] += Gal[p].BlHo[1].Mass_to_acc[j];
   //   #ifdef FEED_ON_STARS
   //    Gal[p].BlHo[0].M_core[j] += Gal[p].BlHo[1].M_core[j]; //add cores together //simplified
   //   #endif
     #ifdef EXTENDED_BH_OUTPUT
   //option here if cumulative accretion mass is stored
      Gal[p].BlHo[0].M_acc_Mergers[j] += Gal[p].BlHo[1].M_acc_Mergers[n_step];
      Gal[p].BlHo[0].M_acc_DI_sec[j] += Gal[p].BlHo[1].M_acc_DI_sec[n_step];
      Gal[p].BlHo[0].M_acc_DI_mmi[j] += Gal[p].BlHo[1].M_acc_DI_mmi[n_step];
      #ifdef FEED_ON_STARS
      Gal[p].BlHo[0].M_acc_TDE[j] += Gal[p].BlHo[1].M_acc_TDE[j];
      #endif
     #endif
     Gal[p].BlHo[0].spin[j] = a_fin;
     Gal[p].BlHo[0].eta_BH[j] = compute_eta(a_fin);}
  #else
  for(j = n_step; j<STEPS; j++)
    {Gal[p].BlHo[0].BlackHoleMass[j] += Gal[p].BlHo[1].BlackHoleMass[j];
     Gal[p].BlHo[0].Mass_to_acc[j] += Gal[p].BlHo[1].Mass_to_acc[j];
   //    #ifdef FEED_ON_STARS
   //    Gal[p].BlHo[0].M_core[j] += Gal[p].BlHo[1].M_core[j];
   //   #endif
     Gal[p].BlHo[0].spin[j] = 0.;}
  V_recoil = 0.0;
  return 
  #endif

  // I re-set the second BH of the central galaxy
  restart_the_BH(p,1,n_step);
  // Some checks
  if(Gal[p].BlHo[0].spin[n_step] >= 1.)
        {printf("Error in add_galaxies_together(): File (%d) Gal[p].BlHo[0].spin[n_step]> 1. (%f with a_fin = %f) , Gal[p].BlHo[0].BlackHoleMass[n_step] = %E\n",filenr,Gal[p].BlHo[0].spin[n_step],a_fin, Gal[p].BlHo[0].BlackHoleMass[n_step]); exit(0);}

  BHf = Gal[p].BlHo[0].BlackHoleMass[n_step];
  af = Gal[p].BlHo[0].spin[n_step];

  #ifdef MODELWANDERINGBHS
   // is_ejected = -1 --> No ejection ; is_ejected == 1 ---> Is ejected
   // If the BH is ejected all the properties of the BH 0 are re-set in the function check_recoil_velocity_BH()
   is_ejected = check_recoil_velocity_BH(p,n_step,V_recoil,&Mdisk_ejection,&index_BH_ejected);
   //if(is_ejected>0) save_file_with_kick_BH(p,n_step,V_recoil,a1,a2,BH1,BH2,qq,BHf,fgas_BHs,filenr);
  #endif

  #ifdef SAVE_ALL_BH_MERGER
  save_merger_events(p,n_step,V_recoil,a1,a2,af,BH1,BH2,qq,BHf,(Mres_tot)/(MBH_tot),filenr,is_ejected);
  #endif
}


void compute_BH_BH_merger_delay(int p, int n_step, int filenr)
{
  // Here I compute the merger delay according to Izquierdo-Villalba et al. 2020
  float fgas, qq;
  float time_binary, t_dealy;

  time_binary = NumToTime(Gal[p].SnapNum) - (((1.0*n_step)/(1.0*STEPS))*(NumToTime(Gal[p].SnapNum) - NumToTime(Gal[p].SnapNum+1))); // internal units
  qq = Gal[p].BlHo[1].BlackHoleMass[n_step] / Gal[p].BlHo[0].BlackHoleMass[n_step]; // Remember that the BH 0 is always more massive than the BH 1
  if(qq>1.) {printf("Error in compute_BH_BH_merger_delay(): qq>1 (%f), Gal[p].t_last_MERGER = %f, tnow = %f \n",qq,Gal[p].t_last_merger,time_binary/(Hubble_h/UnitTime_in_Megayears)); exit(0);}
  Gal[p].fgas = Gal[p].ColdGas/(Gal[p].ColdGas + Gal[p].DiskMass + Gal[p].BulgeMass ); // notice that all the quantities of the satellite galaxy (p) are 0 (already trasnfered in add_galaxies_together 
  fgas = Gal[p].ColdGas/(Gal[p].ColdGas + Gal[p].DiskMass + Gal[p].BulgeMass );

  Gal[p].z_binary_form = redshift_at_t(ZZ[Gal[p].SnapNum],n_step,(NumToTime(Gal[p].SnapNum)-NumToTime(Gal[p].SnapNum+1))/STEPS);
  Gal[p].t_binary_form = time_binary/(3.16e13*Hubble_h/UnitTime_in_s); // [Myr]
  t_dealy =  0.01*(0.1/qq)*(0.3/fgas); // [Gyr]  <--- This is the correct TODO

  Gal[p].BlHo[0].tBH_merger = t_dealy;
  Gal[p].BlHo[0].tBH_merger_clock_dt = Gal[p].BlHo[0].tBH_merger; // [Gyr];
  Gal[p].BlHo[0].tBH_merger_clock_dt = Gal[p].BlHo[0].tBH_merger_clock_dt * (3.16e16*Hubble_h/UnitTime_in_s); // [internal]; Note the 3.16e16
  Gal[p].BlHo[1].tBH_merger = -1.0;
  Gal[p].BlHo[1].tBH_merger_clock_dt = -1.0;
  #ifdef TIMEBHSMERGERALL
   char file_BH_time_merger[1000];
   float t_binary_form,z_binary_form,dt;
   FILE *fileout;
   dt = ( NumToTime(Gal[p].SnapNum) -  NumToTime(Gal[p].SnapNum+1) ) / STEPS;
   t_binary_form = NumToTime(Gal[p].SnapNum) - (((1.0*n_step)/(1.0*STEPS))*(NumToTime(Gal[p].SnapNum) - NumToTime(Gal[p].SnapNum+1))); // internal units
   t_binary_form = t_binary_form/(3.16e13*Hubble_h/UnitTime_in_s); // [Myr]
   z_binary_form = redshift_at_t(ZZ[Gal[p].SnapNum],n_step,dt);
   sprintf(file_BH_time_merger, "/home/david/Documentos/GitGalaxies/LGalaxies/code/Recoil_BH/Hybrid_model/time_BH_mergers/All/Snap_%d", filenr);
   fileout = fopen(file_BH_time_merger,"ab");
   float arr1[8] = {fgas,t_dealy,qq,z_binary_form,Gal[p].BlHo[0].BlackHoleMass[n_step],Gal[p].BlHo[1].BlackHoleMass[n_step],1.0,t_binary_form};
   fwrite(arr1,sizeof(arr1),1,fileout);
   fclose(fileout);
  #endif
}