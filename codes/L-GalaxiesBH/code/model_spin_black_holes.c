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

/** This module deals with black hole spin evolution
 *  The equations of risco, Lisco, Eisco can be found in Jefremov et al. (2015) or in Fiacconi et al. 2017
*/



float compute_rlso(float sp)
{
     float Z1, Z2, rlso, spin, rotacion;
     if(sp >= 0.) rotacion = +1.;
     else rotacion = -1.;
     //if(fabs(sp) >= 0.99998) sp = 0.99998; // Avoiding the nan in (1-spin)^(1./3). See Thorne (1974)
     spin = fabs(sp);
     Z1 = 1. + (pow(1.-pow(spin,2.),1./3.) * (pow(1.+spin,1./3.) + pow(1.-spin,1./3.)));
     Z2 = pow(3.*pow(spin,2.) + pow(Z1,2),1./2.);
     rlso = (3.+Z2) - (rotacion * pow((3.-Z1)*(3.+Z1 + 2.*Z2),0.5));
     return rlso;
}


float compute_eta(float sp)
{
    float rlso, eta;
    rlso = compute_rlso(sp); // Here we are going to take into account the signature of the spin
    eta = 1. - pow(1 - ((2./3.) * (1/rlso)),1./2.);
    return eta;
}


float Compute_Spec_Angular(float sp)
{
   float rlso, spin, rotacion, L;
   if(sp >= 0.) rotacion = 1.;
   else rotacion = -1.;
   spin = fabs(sp);
   rlso = compute_rlso(sp);
   L = rotacion * ( pow(rlso,2.) - (2. * spin * rotacion * pow(rlso,0.5)) + pow(spin,2.)) / (rlso * pow(rlso - 3. + (2.*rotacion*spin*pow(rlso,-0.5)),0.5));
   if(isnan(L) == 1) {printf("Error in Compute_Spec_Angular() the energy is nan. Check why! (sp = %f, spin = %f, rotation = %f, L = %f\n)",sp,spin,rotacion,L);}
   return L;
}


float Compute_Spec_Energy(float sp)
{
   float rlso,spin, rotacion, E;
   if(sp >= 0.) rotacion = 1.;
   else rotacion = -1.;
   spin = fabs(sp);
   rlso = compute_rlso(sp);
   E = ( pow(rlso,3./2.) - 2.*pow(rlso,1./2.) + (rotacion * spin) ) / (rlso * pow(rlso - 3. + (2.*rotacion*spin*pow(rlso,-0.5)),0.5));
   if(isnan(E) == 1) {printf("Error in Compute_Spec_Energy() the energy is nan. Check why! (sp = %f, spin = %f, rotation = %f, E = %f\n)",sp,spin,rotacion,E);}
   return E;
}


float update_spin(float Mini,float Mfin, float sp)
{
        float rlso, spin_final;
        rlso = compute_rlso(sp);
        if( (Mfin/Mini) > pow(rlso,0.5) ) spin_final = 0.998;
        else spin_final = (1./3.) * pow(rlso,0.5) * (Mini/Mfin) * ( 4. - pow((3.*rlso*pow(Mini/Mfin,2.)) -2.,0.5));
        //if(spin_final<0.) {printf("In update_spin: Spin final negative!\n"); exit(0);}
        if(spin_final>=1.) spin_final = 0.998; //{printf("In update_spin: Spin final >= 1. (%f)!\n",spin_final); exit(0);} //spin_final = 0.998;
        return spin_final;

}


void Self_Gravity_Accretion_Sesana_Model(int p, float BlackHole, float M_dot, float fEdd, float *a, float *eta, float epsilon, float v_sigma, double dt, int filenr)
{

  // This function tracks the spin evolution during gas accretion. It is based on the model of Dotti et al. 2013 and Sesana et al. 2014. 
  // Basically links the spin evolution (number of prograds accretion disks) with the bulge morphology of the galaxy
  // The larger is the bulge coherence (pseudobulges/bars) the larger is the final spin of the black hole after the gas accreiton
  
  float alpha_visco1, Rdisk,Rsch, Av,ratio;
  float Msg, Mcloud, w, F, a_dot;
  int j;
  alpha_visco1 = 0.1;
  Mcloud = 7.e3 * 1e-10 * Hubble_h; // M_sun; Now are internal units --> This parameter controls the mass where the BHs change form maximum spinning to no maximaly spinning (See Dotti et al. 2013 and Izquierdo-Villalba et al. 2020)
  F = 0.;
  if((*a) < 0 || *a >= 1) {printf("Error in self_gravity_accretion_Sesana() the spin is negative or 1 (a* = %f)...but should not. The spin is always positive ...Wierd....BH = %f\n",*a,log10(BlackHole/(Hubble_h*1.e-10))); exit(0);}
  if(isnan(*a) == 1) {printf("Error its nan\n"); printf("Input parameters = BlakHole = %E, M_dot = %E, fEdd = %f, spin = %f, eta = %f, epsilon = %f M_dot = %E, v_sigma = %E, dt = %E, \n",BlackHole,M_dot,fEdd,*a,*eta,epsilon,M_dot,v_sigma,dt);}
  if(M_dot<=0.){ printf("Error in self_gravity_accretion: The M_dot < 0. (%E), MBH = %E\n",M_dot,BlackHole*1e10/Hubble_h); exit(0);}

  Msg = 2.e4 * pow(alpha_visco1 / 0.1, -1./45.) * pow(0.1/(*eta), 4./45.)* pow(fEdd,4./45.) * pow(1.e4 * BlackHole / Hubble_h, 34./45.); // M_sun
  Msg = Msg * 1e-10 * Hubble_h; // Internal units
  Rsch = 2.96e13 * 100. * BlackHole / Hubble_h; // cm
  if(Msg < Mcloud)
           Rdisk = 1.e5 * Rsch * pow(alpha_visco1/0.1,28./45.) * pow(1.e4 * BlackHole / Hubble_h,-52./45.) * pow(fEdd,-22./45.) * pow(0.1/(*eta),-22./45.);
  else
           Rdisk = 4.e4 * Rsch * pow(1.e6*Mcloud/Hubble_h,4./5.) * pow(alpha_visco1/0.1,16./25.) * pow(1.e4 * BlackHole / Hubble_h,-44./25.)* pow(fEdd,-14./25.) * pow(0.1/(*eta),-14./25.);
  Av = 9.14e6*pow(alpha_visco1 / 0.1,4./5.) * pow(1.e4 * BlackHole / Hubble_h,1./20.) * pow(fEdd,3./10.) * pow(0.1/(*eta),3./10.); // cm^5/4 s^-1
  Av = Av * pow(Hubble_h/UnitTime_in_s,-1.); // cm^5/4 s^-1 --> cm^5/4 s_internal^-1
  ratio = 0.5*(8./21.)*pow(2,1./2.)* M_dot * pow(Rdisk,7./4.) / (BlackHole*(*a)*Av*pow(Rsch,1./2.));//Jdisk/JBH;

  if(isnan(ratio) == 1) {printf("Error in function Self_Gravity_Accretion_Sesana_Model() the value of the ratio is nan ! ######\n "); exit(0);}
  if(ratio>1.) w = 1.;
  else
        {
                if(v_sigma>= 1.) F = 1.;
                else if(v_sigma == 0.) F = 0.5;
                else if(v_sigma >0. && v_sigma<1.) F = (1.+v_sigma)/2.;
                else {printf("Error in v_sigma in Self_Gravity_Accretion_Sesana_Model(): The value of v_sigma is not undestood. v_sigma = %f, filenr = %d, spin = %f, eta = %f, M_dot = %E, fedd = %f, BHmass = %f, dt = %E, BulgeMass = %E, ColdGas = %E, DiskMass = %E\n",v_sigma,filenr,*a,*eta,M_dot,fEdd,log10(BlackHole*1e10/Hubble_h),dt,Gal[p].BulgeMass,Gal[p].ColdGas,Gal[p].DiskMass); for(j=0;j<STEPS;j++) printf("BH1 Mass_to_acc[%d]=%E, BH2 Mass_to_acc[%d] = %E\n",j,Gal[p].BlHo[0].Mass_to_acc[j],j,Gal[p].BlHo[1].Mass_to_acc[j]); exit(0);}
                #ifdef CHAOTICBHSCENARIO
                F = 0.5;
                #endif
                w = F + (ratio*(1.-F));
                if(w < 0 || isnan(w) == 1) {printf("Error in Self_Gravity_Accretion_Sesana_Model(): The value of w is 0 or nan (w = %f, ratio = %f, F = %f). Check it.\n",w,ratio,F); exit(0);}
                //w = 0.5; //With this value we get a similar model to King's one!
        }

  if(M_dot <=0) {printf("Error in the Self_Gravity_Accretion_Sesana_Model(). The accretion rate onto the BH is negative: M_dot(%f)<=0",M_dot); exit(0);}
  if(dt<=0) {printf("Error in the Self_Gravity_Accretion_Sesana_Model(). The value of dt is negative: dt = %f <=0",dt); exit(0);}
  a_dot = ( ( w*Compute_Spec_Angular(*a) + (1.-w)*Compute_Spec_Angular(-1.*(*a)) ) - (2.*fabs(*a)*( w*Compute_Spec_Energy(*a) + (1.-w)*Compute_Spec_Energy(-1.*(*a)) ) ) ) * (M_dot/BlackHole);
  *a = fabs(*a) + (a_dot * dt);

  if(*a>=0.998) *a = 0.998; // The BH spin can NOT reach the a == 1 according to Kip Thorn 
  if(*a<0.) *a = 0.01;     // Negative means that the accretion decreased the spin up to remove it
  if(*a>=1. || *a<-1. || isnan(*a) == 1) {printf("Error in Self_Gravity_Accretion_Sesana_Model(): The spin bigger/sammler +-1 or nan (a = %f): SelfGravity mass = %E, BHMass_ini = %E, M_dot, = %E delta_episode = %E, fEdd = %f, v_sigma = %f\n",*a, Msg, BlackHole, M_dot, Msg/BlackHole,fEdd,v_sigma);}// exit(0);}
  *eta = w*compute_eta(*a) + ( (1.-w)*compute_eta(-1.*(*a)) );
}


void BH_BH_spin_mass(int p, int n_step, float *a_fin, float *M_bh_f_GW, float *V_recoil, int filenr)
{
  // Here I compute the final spin of the final black hole after a BH-BH coalescence
  // The equation used here are Lousto et al. (2012) 
  // We do an approximation: If the two BH involced in the merger are surrounded by a lot of gas (Mres1+Mres2>MBH1+MBH2) the gas torque will aling the black holes (i.e the spins point towards the same direction)
  // in case the two black holes are surrounded by low amount of gas, the direction of the two BH spins will be random
  
  float Mres_tot, MBH_tot, qq, mu_spin, a1, a2,beta, gamma, alpha, l_orbital_momentum ,t0,t2,t3,s4,s5,p0,p1,p2,m0,m1;
  t0 = -2.8904;
  t2 = -3.5171;
  t3 = 2.5763;
  s4 = -0.1229;
  s5 = 0.4537;
  p0 = 4.826e-2;
  p1 = 1.559e-2;
  p2 = 0.485e-2;
  m0 = 0.9515;
  m1 = -0.013;

  if(Gal[p].BlHo[1].BlackHoleMass[n_step] < 1e-14 || Gal[p].BlHo[0].BlackHoleMass[n_step] < 1e-14)
        {printf("Error in BH_BH_spin_mass(): One of the BHs has mass == 0. MBH1 = %E, MBH2 = %E\n",Gal[p].BlHo[0].BlackHoleMass[n_step],Gal[p].BlHo[1].BlackHoleMass[n_step]); exit(0);}
  qq = Gal[p].BlHo[1].BlackHoleMass[n_step] / Gal[p].BlHo[0].BlackHoleMass[n_step];
  mu_spin = qq/pow(1.+qq,2.);

  if(qq>1.0) {printf("Error in BH_BH_spin_mass(): The value of qq > 1 (%f). Filenr = %d",qq,filenr); exit(0);}
  if(isnan(qq) == 1){printf("Error in BH_BH_spin_mass: The qq value is nan!  (MBH1 = %E, MBH2 = %E)\n",Gal[p].BlHo[0].BlackHoleMass[n_step],Gal[p].BlHo[1].BlackHoleMass[n_step]); exit(0);}

  a1 = fabs(Gal[p].BlHo[0].spin[n_step]);
  a2 = fabs(Gal[p].BlHo[1].spin[n_step]);

  Mres_tot = Gal[p].BlHo[0].Mass_to_acc[n_step] + Gal[p].BlHo[1].Mass_to_acc[n_step];
  MBH_tot  = Gal[p].BlHo[0].BlackHoleMass[n_step] + Gal[p].BlHo[1].BlackHoleMass[n_step];


  if( Mres_tot > MBH_tot )
     {
	//Dotti et al. 2010 said the aligment is not perfect but it might have an offset of [0-10] deg
        beta  = 0.98481 + ((1.0-0.98481)*(((double)rand()/(double)RAND_MAX)));//cos(0.174533 * (((double)rand()/(double)RAND_MAX))); // [0,10] Dotti 2010
        gamma = 0.98481 + ((1.0-0.98481)*(((double)rand()/(double)RAND_MAX)));//cos(0.174533 * (((double)rand()/(double)RAND_MAX))); // [0,10] Dotti 2010
        alpha = cos(acos(beta) - acos(gamma)) + ( (cos(acos(beta) + acos(gamma)) - cos(acos(beta) - acos(gamma)) ) *  (((double)rand()/(double)RAND_MAX)) );
     }
   else
    {
        beta  = -1. + ( 2.* (((double)rand()/(double)RAND_MAX)) );
        gamma = -1. + ( 2.* (((double)rand()/(double)RAND_MAX)) );
        //alpha = -1. + ( 2.* (((double)rand()/(double)RAND_MAX)) );//(M_PI)*gsl_rng_uniform(random_generator); //(2.*gsl_rng_uniform(random_generator)) -1.; // 0 < cos gamma < -1 //(2.*M_PI)*gsl_rng_uniform(random_generator);
        alpha = cos(acos(beta) - acos(gamma)) + ( (cos(acos(beta) + acos(gamma)) - cos(acos(beta) - acos(gamma)) ) *  (((double)rand()/(double)RAND_MAX)) );
     }

  if( ( Gal[p].BlHo[0].BlackHoleMass[n_step] == Gal[p].BlHo[1].BlackHoleMass[n_step] )  && ( Mres_tot >= MBH_tot ) )
  *M_bh_f_GW = (1. - (p0 + 2.*p1*(a1+a2) + 4*p2*pow(a1+a2,2.)) ) * MBH_tot; // Gravitational wave losses Reisswig et al. (2009)
  else
      *M_bh_f_GW = MBH_tot * ((1. + 4.*mu_spin*(m0-1.) + 16.*m1*pow(mu_spin,2.)*(a1*(beta) + a2*(gamma))));// Gravitational wave losses Tichy & Marronetti (2008)

   l_orbital_momentum = 2.*pow(3,1./2.)+ ( t2*mu_spin ) + ( t3*pow(mu_spin,2.) ) +  (  (s4/pow(1.+pow(qq,2.),2.) )  * ( pow(a1,2.) + ( pow(a2,2.)*pow(qq,4.) )  +  ( 2*a1*a2*pow(qq,2.)*(alpha)) ) ) + ( ( ( (s5*mu_spin) + t0 + 2.)/(1.+pow(qq,2.)) )* ( a1*(beta) + a2*(gamma)*pow(qq,2.) )) ;

   *a_fin = (1./pow(1.+qq,2.))* pow( pow(a1,2.) + ( pow(a2,2.)*pow(qq,4.) ) + ( 2*a1*a2*pow(qq,2) * (alpha)) + 2*( (a1*(beta) + a2*pow(qq,2.)*(gamma)))*l_orbital_momentum*qq  + ( pow(l_orbital_momentum,2.)*pow(qq,2.) ),1./2. );

  if(*a_fin >=1.) *a_fin = 0.998; // See Bert & Volonteri (2008)
  if(*a_fin >=1.) {printf("Error in the mergers of two BHs, the spin resulting is >1.. Check why?. Some information: Gal[p].BlHo[0].BlackHoleMass[n_step] = %E, Gal[p].BlHo[1].BlackHoleMass[n_step] = %E, a1 = %f, a2 = %f, qq = %E, mu = %E, l_orbital = %E, a_fin = %E, beta = %f, gamma = %E, alpha = %f, t0 = %E, t1 = %E, t2 = %E; For a_fin = a0 = %E, a1 = %E, (a0 + a1)^1/2 = %E\n, ",Gal[p].BlHo[0].BlackHoleMass[n_step],Gal[p].BlHo[1].BlackHoleMass[n_step],a1,a2,qq,mu_spin,l_orbital_momentum, *a_fin, beta, gamma, alpha, 2.*pow(3,1./2.)+ ( t2*mu_spin ) + ( t3*pow(mu_spin,2.) )  , (  (s4/pow(1.+pow(qq,2.),2.) )  * ( pow(a1,2.) + ( pow(a2,2.)*pow(qq,4.) )  +  ( 2*a1*a2*pow(qq,2.)*(alpha)) ) ), ( ( ( (s5*mu_spin) + t0 + 2.)/(1+pow(qq,2.)) ) * ( a1*(beta) + a2*(gamma)*pow(qq,2.) )), pow(a1,2.) + ( pow(a2,2.)*pow(qq,4.) ) + ( 2*a1*a2*pow(qq,2) * (alpha)), 2*( (a1*(beta) + a2*pow(qq,2.)*(gamma)))*l_orbital_momentum*qq  + ( pow(l_orbital_momentum,2.)*pow(qq,2.) ), pow( pow(a1,2.) + ( pow(a2,2.)*pow(qq,4.) ) + ( 2*a1*a2*pow(qq,2) * (alpha)) + 2*( (a1*(beta) + a2*pow(qq,2.)*(gamma)))*l_orbital_momentum*qq  + ( pow(l_orbital_momentum,2.)*pow(qq,2.) ),1./2. ) ); exit(0);}

  if(isnan(*a_fin) == 1) {printf("Error in the mergers of two BHs, the spin resulting is NAN. Check why?. Some information: Gal[p].BlHo[0].BlackHoleMass[n_step] = %E, Gal[p].BlHo[1].BlackHoleMass[n_step] = %E, a1 = %f, a2 = %f, qq = %E, mu = %E, l_orbital = %E, a_fin = %E, beta = %f, gamma = %E, alpha = %f, t0 = %E, t1 = %E, t2 = %E; For a_fin = a0 = %E, a1 = %E, (a0 + a1)^1/2 = %E\n, ",Gal[p].BlHo[0].BlackHoleMass[n_step],Gal[p].BlHo[1].BlackHoleMass[n_step],a1,a2,qq,mu_spin,l_orbital_momentum, *a_fin, beta, gamma, alpha, 2.*pow(3,1./2.)+ ( t2*mu_spin ) + ( t3*pow(mu_spin,2.) )  , (  (s4/pow(1.+pow(qq,2.),2.) )  * ( pow(a1,2.) + ( pow(a2,2.)*pow(qq,4.) )  +  ( 2*a1*a2*pow(qq,2.)*(alpha)) ) ), ( ( ( (s5*mu_spin) + t0 + 2.)/(1+pow(qq,2.)) )* ( a1*(beta) + a2*(gamma)*pow(qq,2.) )), pow(a1,2.) + ( pow(a2,2.)*pow(qq,4.) ) + ( 2*a1*a2*pow(qq,2) * (alpha)), 2*( (a1*(beta) + a2*pow(qq,2.)*(gamma)))*l_orbital_momentum*qq  + ( pow(l_orbital_momentum,2.)*pow(qq,2.) ), pow( pow(a1,2.) + ( pow(a2,2.)*pow(qq,4.) ) + ( 2*a1*a2*pow(qq,2) * (alpha)) + 2*( (a1*(beta) + a2*pow(qq,2.)*(gamma)))*l_orbital_momentum*qq  + ( pow(l_orbital_momentum,2.)*pow(qq,2.) ),1./2. ) ); exit(0);}
  //if(isnan(*a_fin) == 1) {printf("in nan, beta = %f, gamma = %f, alpha = %f, l_orbital_momentum  %f, *a_fin = %E\n",beta,gamma,alpha,l_orbital_momentum,*a_fin); exit(0);}

  if(*M_bh_f_GW > (MBH_tot) || *M_bh_f_GW < 0. || isnan(*M_bh_f_GW) == 1) {printf("Error in the BH_BH_spin_mass(): The final BH mass is larger after the GW losses than without (or is negative): MBHf = %E\n",*M_bh_f_GW); exit(0);}
 #ifdef SPINCOMPUTATION
 #ifdef MODELWANDERINGBHS
  *V_recoil = kick_velocity(a1,a2,beta,gamma,qq,mu_spin,Gal[p].SnapNum,filenr);
 #else
  *V_recoil = 1e-9;
 #endif
 #else
  *V_recoil = 1e-9;
 #endif
}

