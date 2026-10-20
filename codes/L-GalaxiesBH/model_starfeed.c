#ifdef FEED_ON_STARS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"


  
void Feed_On_Stars(int p,int nstep ,int BH_indx,int filenr,float time)
{
  float Minnerstars,  nb, Re, MBHtot, Mincluded; //#MP
  float event_rate_TDE , max_stellar_mass_removed, rho_0;
  float rho_inf, r_inf; //sigma_inf computed in T_relax, these are used #MP
  float Ms, Rs, corr; // simplified
  #ifndef TIME_DEPENDENT_TDE
  double  T_relax;
  #endif
  //simplified rho_TDEcrit  ->0.0 rcore_min_lim ->1.e-9
  float rcore_min_lim = 1.e-9/Hubble_h; // h^-1 Mpc , a condition does not go into calculations
  float m_monochromatic = 1e-10*Hubble_h; //10^10 Msun
  float rho_TDEcrit = 0.0;// simplified (30.0*m_monochromatic)/pow(rcore_min_lim,3); //mass in central pc^3 must be above that to consider the process
  float qa, qb; //for constant rate of events, if binary split rate to half

  if(((Gal[p].BulgeMass>Gal[p].BlHo[0].BlackHoleMass[nstep])&&((Gal[p].BulgeSize>rcore_min_lim))) || ((Gal[p].DiskMass>Gal[p].BlHo[0].BlackHoleMass[nstep])&&(Gal[p].StellarDiskRadius> rcore_min_lim)))
    {
      MBHtot = Gal[p].BlHo[0].BlackHoleMass[nstep]+Gal[p].BlHo[1].BlackHoleMass[nstep];
      Minnerstars = Gal[p].BulgeMass;
      Re  = Gal[p].BulgeSize;// [Mpc]
      nb = Gal[p].Sersic_index;
      if (nb <0.0) nb=1.0;
      if (nb < 0.01) nb=0.01;
      if(((Minnerstars <= Gal[p].BlHo[0].BlackHoleMass[nstep]) || (Re<=rcore_min_lim)))
      {
          Minnerstars = Gal[p].DiskMass; // [10^10 Msun]
          Re = (1.68/3.)*Gal[p].StellarDiskRadius; // [Mpc] The 1.68 is to transform the scale length into effective radius, I leave the h because in the mass is also the h
          //if(Gal[p].Sersic_index < 0.) 
          nb = 1.0;
          Gal[p].Sersic_index = 1.0;
      }
      //rho_0 = compute_rho_0_Sersic_Bulge(Minnerstars,nb,Re);
      //METHOD 2 (METHOD 1 COMMENTED OUT BELOW used BH_Growth_TDE)
      //r_inf =  compute_influence_radius(MBHtot,rho_0,nb,Re); //can be temporal  or directly store to Gal[p]

      //if((rinf<=0.0) || isnan(rinf)==1) printf("\n r_inf= %f pc",rinf*1e6);
      //r_inf =  3.0e-3 * MBHtot;
      //if (MBHtot>1e-5) printf("\n sigma = %f", Gal[p].sigma_inf[nstep]);
      //rho_inf = density_bulge_Sersic_model(rho_0, nb, Re, r_inf); //can be temporal or directly store to Gal[p]
      rho_inf = 1.0;//simplified from above , rhoTDEcrit = 0.0
      if(rho_inf > rho_TDEcrit){
      //   if (ZZ[Gal[p].SnapNum]<1.0) {
      //   printf("\n just out of M_include time= %f (z=%f), rho_inf = %e Msun/pc^3, nb = %e , rinf/Re= %f",previoustime - (nstep + 0.5) * (deltaT / STEPS), ZZ[Gal[p].SnapNum],rho_inf/1.e8,nb,r_inf/Re);  
      //   fflush(stdout);
      //  }
        Mincluded = Minnerstars; //compute_included_mass(20*r_inf,rho_0,nb,Re); //simplified
        max_stellar_mass_removed  = Mincluded;
        //#ifdef EXTENDED_BH_OUTPUT
        //Gal[p].sigma_inf[nstep]  = sigma_Sersic_profile_r(r_inf ,nb, Re, rho_0, rho_inf);  // put in rho>rho_crit loop
        //Gal[p].rho_inf[nstep] = rho_inf;//can be temporal  or directly store to Gal[p]
        //Gal[p].r_inf[nstep]  =  r_inf; 
        //#endif
        #ifndef TIME_DEPENDENT_TDE
          #ifdef EXTENDED_BH_OUTPUT
          Gal[p].BlHo[0].T_relax[nstep] = compute_T_relax(p, nstep, m_monochromatic);
          #else
          Gal[p].BlHo[0].T_relax[nstep] = compute_T_relax(p, r_inf,  rho_inf, rho_0,  nb, Re, m_monochromatic);
          #endif
          //printf("\n T_relax = %f ",log10(Gal[p].BlHo[0].T_relax[nstep]));

          event_rate_TDE = compute_TDE_rate(p,nstep,0, m_monochromatic); //internal units 
          if(isnan(max_stellar_mass_removed)==1 && Gal[p].SnapNum<58) {printf("\n Maximum Stellar_mass_removed nan");}
          if(isnan(event_rate_TDE)==1 && Gal[p].SnapNum<58) {printf("\n event_rate_TDE nan");}
          max_stellar_mass_removed  = event_rate_TDE*m_monochromatic*Gal[p].BlHo[0].T_relax[nstep];
        #endif

          if(Gal[p].BulgeMass>Gal[p].BlHo[0].BlackHoleMass[nstep])
          {
          if ( max_stellar_mass_removed > Minnerstars*MaxFracBulgeStarFeed)
          { max_stellar_mass_removed  =  Minnerstars*MaxFracBulgeStarFeed;}
          }
          else{
          if (max_stellar_mass_removed > Minnerstars*MaxFracDiskStarFeed)
          {  max_stellar_mass_removed  =  Minnerstars*MaxFracDiskStarFeed;}
          }
          if (Mincluded > max_stellar_mass_removed) {Mincluded = max_stellar_mass_removed;} //reduce if relation with stellar mass is set

     // if (Mincluded > 2*Gal[p].BlHo[0].M_core[nstep]) //simplified
      if (Mincluded > 0.0)
       
        { 
          
      // Sersic profile with total mass (Greene+2020):
        Ms = 3.0 *Hubble_h* pow(10,(log10(MBHtot/Hubble_h)+10.0-7.43)/1.61);
      //and the scale radius depends on the above mass as (Shen+2003):
        Rs = 1.0e-3*Hubble_h*pow(10,(0.14*log10(Ms*1e10/Hubble_h)-1.21));
        corr = (Mincluded/Ms) * pow(Re/Rs, -1.5);
        if (corr<1.0) corr=1.0; //printf("\n corr = %f",corr);
        if (corr>1000.0) {printf("\n corr = %f",corr); corr=1000.0;}
          //printf("\n  T_relax = %f, sigma_inf=%f , r_inf = %f , rho_inf = %f",log10(Gal[p].BlHo[0].T_relax[nstep]),log10(Gal[p].sigma_inf[nstep]) ,log10(r_inf), log10(rho_inf));
          Gal[p].BlHo[0].M_res_TDE[nstep] = 0.0;//Mincluded; 
          Gal[p].BlHo[0].M_acc_TDE[nstep] = Mincluded;  
          BH_growth_TDE_instant(p, nstep, BH_indx, Mincluded, corr, time, m_monochromatic,filenr);
         //if (ZZ[Gal[p].SnapNum]<6.0) 
         //{printf("\n z=%f TDE is happening, Mincluded =%f , Minner = %f,  max_stellar_mass_removed  %f, event_rate  %f ",ZZ[Gal[p].SnapNum],log10(Mincluded*1.0e10/Hubble_h), log10(Minnerstars*1.0e10/Hubble_h),  log10(max_stellar_mass_removed*1.0e10/Hubble_h), log10(event_rate_TDE*1e-12*Hubble_h));
        //  fflush(stdout);
        // }
          //printf("\n %f ",  max_stellar_mass_removed );
        
        //if(Mincluded>5e-1*Minnerstars/Hubble_h)  //simplified
        //{ printf("\n Not a problem, just curious: Core is 1/2 of the StellarMass with log10 Mcore %f at z=%f",log10(Mincluded), ZZ[Gal[p].SnapNum]);}


        //if (Gal[p].BlHo[1].BlackHoleMass[nstep] > 1e-9*Hubble_h) { qa=1.0;  qb=1.0;} 
        //else{ qa=1.0; qb=0.0;}
        // qa=1.0; 
        // qb=0.0;
        // Gal[p].BlHo[0].t_TD_phase  = -1.0; //turn off switch (useful if TDE feeding happening elsewhere in the code,
        // BH_growth_TDE_direct(p, nstep,0 , Mincluded, qa, time, m_monochromatic, filenr );
        // if (qb>0.0) {
        //   Gal[p].BlHo[1].t_TD_phase  = -1.0; //especially for later versions where T_relax is recomputed 
        // BH_growth_TDE_direct(p, nstep,1, Mincluded, qb, time, m_monochromatic, filenr);
        // }
        // if((Gal[p].BlHo[1].BlackHoleMass[0]>1e-14) && (Gal[p].BlHo[1].BlackHoleMass[nstep]>Gal[p].BlHo[0].BlackHoleMass[nstep]))
        // exchange_central_BH_positions(p);
        //printf("\n TDE feeding is activated at redshift = %f,new logMBH = %f ", ZZ[Gal[p].SnapNum], log10(Gal[p].BlHo[0].BlackHoleMass[nstep]));
        //MASS CHECK (if we have borrowed mass from the bulge to finish the >M_dot_star accretion)
        //M_core in step nstep should be Mincluded - stars accreted by BH immediately,
        // we do not subtract event_rate_TDE*Mstar because the rate might be lower or higher in BH_growth_TDE
        //if greater than event_rate_TDE rates are allowed (in the future), then we have to make bulge's artificially greater

        //Central region is allowed 1) to be consumed 2) to not be disrupted > mass checks later so that it never exceeds rest of stellar mass components
         //establish feeding from the core in the future

         //simplified
        // if ((Gal[p].BlHo[0].M_dot_TDE[nstep]>event_rate_TDE)|| (Gal[p].BlHo[1].M_dot_TDE[nstep]>event_rate_TDE))
        // {
        // if(Gal[p].BlHo[0].M_core[nstep] > Gal[p].DiskMass/Hubble_h) {  
        //   transfer_stars(p,"Disk",p,"ICM",(Gal[p].BlHo[0].M_core[nstep]*Hubble_h - Gal[p].DiskMass)/(Gal[p].ICM*Hubble_h));
        //   printf("Transfering from ICM to disk !");}
        // if(Gal[p].BlHo[0].M_core[nstep] > Gal[p].BulgeMass/Hubble_h) { 
        //    transfer_stars(p,"Bulge",p,"ICM",(Gal[p].BlHo[0].M_core[nstep]*Hubble_h - Gal[p].BulgeMass)/(Gal[p].ICM*Hubble_h));
        //    printf("Transfering from ICM to bulge !");}
        // }

        //if (ZZ[Gal[p].SnapNum]<0.05) printf("Core log Mcore= %f, z=%f",log10(Gal[p].BlHo[0].M_core[nstep]), ZZ[Gal[p].SnapNum]);
            // printf("log Mb , log Mbh , nb , log Re =  %f , %f , %f , %f \n",log10(Minnerstars),log10(MBHtot), nb,log10(Re));
            //rho_0 = compute_rho_0_Sersic_Bulge(Minnerstars,nb,Re);
            //rinf = compute_influence_radius(MBHtot,rho_0,nb,Re);
            //if((rinf<=0.0) || isnan(rinf)==1) 
            //printf("\n r_inf= %f pc",rinf*1e6);
            //if (rinf<=0.0) rinf =  3.0e-3 * (Gal[p].BlHo[0].BlackHoleMass[STEPS-1]+Gal[p].BlHo[1].BlackHoleMass[STEPS-1]);
            
          //  float arr_acc[8] = {ZZ[Gal[p].SnapNum],Gal[p].BlHo[0].BlackHoleMass[STEPS-1], Gal[p].BulgeMass, Gal[p].DiskMass, rho_0, rinf, sigma,T_relax};
          //  fwrite(arr_acc,sizeof(arr_acc),1,fileout);
        } //end if Mincluded > 2Mcore
      }  //end if rho_crit
      else{ Mincluded=0.0;}
    
      //#ifdef FEED_ON_STARS       //METHOD 1 
      //event_rate_TDE = compute_TDE_rate(p,nstep,0, m_monochromatic); //internal units 
      // #ifdef EXTENDED_BH_OUTPUT
      // Gal[p].sigma_inf[nstep]  = sigma_Sersic_profile_r(Gal[p].r_inf[nstep] ,nb, Re, rho_0);
      // Gal[p].BlHo[0].T_relax[nstep] = compute_T_relax(p, nstep, m_monochromatic);
      // #else
      // Gal[p].BlHo[0].T_relax[nstep] = compute_T_relax(p, MBHtot,  rho_0,  n,  Re, m_monochromatic);
      // #endif     
      //BH_growth_TDE(p, nstep,0,event_rate_TDE*qq, m_monochromatic, filenr);
      //if (qb>0.0) BH_growth_TDE(p, nstep,1, event_rate_TDE*qq, m_monochromatic, filenr);

      //for (Gal[p].BlHo[1].BlackHoleMass[nstep] > BlackHoleSeedMass*1e-10*Hubble_h){
      //qq = Gal[p].BlHo[bh_idx].BlackHoleMass[nstep] / Gal[p].BlHo[0].BlackHoleMass[nstep];
      // BH_growth_TDE(p, nstep,1, event_rate_TDE*qq/(1+qq), m_monochromatic, filenr);
      // BH_growth_TDE(p, nstep,0, event_rate_TDE*1/(1+qq), m_monochromatic, filenr);
      //Grow second black hole as well but with a different rate
      // event_rate_plunges = compute_plunges_rate(p,nstep,0, m_monochromatic);
      // BH_growth_plunges(p, nstep, BH_indx, event_rate_TDE, m_monochromatic)
      // event_rate_EMRIs = compute_EMRIS_rate(p,nstep,0, m_monochromatic);
      // BH_growth_EMRIs(p, nstep, BH_indx, event_rate_TDE, m_monochromatic) }
      //else{
      // BH_growth_TDE(p, nstep,0,, event_rate_TDE, m_monochromatic, filenr);}
      //#endif      //end FEED_ON_STARS
    
    }//end if bulge/disk/bh  >bhseed and bulgesize,disksize >1e-6
  }

#ifndef TIME_DEPENDENT_TDE
float compute_TDE_rate(int p , int j, int BH_indx, float m_monochromatic)   //function of redshift and primary/secondary dependence too 
{
return  1e-4 /(MYRs_TO_INTERNAL/1.0e6) ; // number of events * h/Mpc * Km/s - code units
}
#else
float compute_TDE_rate(int p , int j, int BH_indx, float m_monochromatic, float time_now)   //function of redshift and primary/secondary dependence too 
{
int kk,ll, ns; 
int inarray=-1;
float dt, td1, td2, td3, td4, mbhtde, rate;
dt = Gal[p].BlHo[0].t_TD_phase - time_now;
mbhtde = Gal[p].BlHo[BH_indx].BlackHoleMass[j];
ns = round(Gal[p].Sersic_index);
if (ns>6) ns=6;
if (ns<1) ns=1;

if(Gal[p].BlHo[BH_indx].BlackHoleMass[j]>1e-6*Hubble_h) mbhtde =1.0001*1.0e-6*Hubble_h;
if(Gal[p].BlHo[BH_indx].BlackHoleMass[j]<1e-2*Hubble_h) mbhtde = 0.9999*1.0e-2*Hubble_h;
for (ll=0;ll<NUM_MBH_BINS-1;ll++)
  {
    for (kk=0;kk<NUM_MBH_BINS-1;kk++)
      {
      if ((TDEtimes[kk]< dt < TDEtimes[kk+1])&&(TDEmbh[ll]< mbhtde < TDEmbh[ll+1]))
      inarray=1;
      break;
      }
   if(inarray==1) break;
 }
 if(inarray!=1){printf("\n Apparently, the TDE rate array does not cover the parameter space inputed in the compute_TDE_events. \n dt= %f Mbh = %f ",dt,Gal[p].BlHo[BH_indx].BlackHoleMass[j]); exit(0);}
 td1 = TDErates[ns][kk][ll];
 td2 = TDErates[ns][kk][ll+1];
 td3 = TDErates[ns][kk+1][ll];
 td4 = TDErates[ns][kk+1][ll+1];
 // either interpolate between them linearly, or in logscale or use maximum
 rate = (0.5*(td1+td2)*log(mbhtde/TDEmbh[kk])+0.5*(td3+td4)*log(TDEmbh[kk+1]/mbhtde))/log(TDEmbh[kk+1]/TDEmbh[kk]);
 //printf("\n Interpolated TDE rate (internal units): %e",rate);
 return  rate;
   }
#endif

#ifdef EXTENDED_BH_OUTPUT
float compute_T_relax(int p,int j, float m_monochromatic)
{ float norm_Trelax = 0.081;
  float log_Lambda_Coulomb =10.0;
  //printf("rho_0 = %E rinf = %E\n",rho_0, rinf);
  return norm_Trelax* pow(Gal[p].sigma_inf[j],3)/(pow(G_GRAVITY,2)*Gal[p].rho_inf[j] * m_monochromatic *log_Lambda_Coulomb); 
}
#else
float compute_T_relax(int p, float r_inf, float rho_inf, float rho_0, float n, float Re, float m_monochromatic)
{ 
  float sigma_inf;
  float norm_Trelax = 0.081;
  float log_Lambda_Coulomb =10.0;
  //printf("rho_0 = %E rinf = %E\n",rho_0, rinf);
  sigma_inf =  sigma_Sersic_profile_r(r_inf, n,  Re,  rho_0 ,rho_inf);
  return norm_Trelax* pow(sigma_inf,3)/(pow(G_GRAVITY,2)* rho_inf* m_monochromatic *log_Lambda_Coulomb);  
}
#endif
#endif