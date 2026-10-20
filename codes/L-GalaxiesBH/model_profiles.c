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

/** This module deals with galaxy and halo profiles: NFW, Sersic, Exponential
*/

float NFW_g_function(float value)
{
   float g;
   g = log(1.0+value) - ((1.0*value)/(1.0+value));
   return g;
}


float compute_concentration_DM_halo(int p, int n_step)
{
   // Extracted form the paper Dutton & Maccio 2014: Cold dark matter halos in the Planck era
   float concentration, zz, a, b, dt;
   dt = ( NumToTime(Gal[p].SnapNum) -  NumToTime(Gal[p].SnapNum+1) ) / STEPS;
   zz = redshift_at_t(ZZ[Gal[p].SnapNum],n_step,dt);
   zz = ZZ[Gal[p].SnapNum];
   if(zz<=0) zz = 0.0;
   if(zz>5.0) zz = 5.0;
   a = 0.537 + (1.025-0.537)*exp(-0.718*pow(zz,1.08));
   b = -0.097+(0.024*zz);
   concentration = a + (b * log10(Gal[p].Mvir/100.)); // log_10 c = a + b log_10 (M/[10^12 M_sun/h])
   concentration = pow(10.,concentration);
   if(isnan(concentration) == 1) {printf("Error in compute_concentration_DM_halo(). The value of concenrtration isn nan (%f). p = %d, n_step = %d, zz = %f,dt = %E,a =%f, b = %f, Mvir = %E\n",concentration,p,n_step,zz,dt,a,b,Gal[p].Mvir); exit(0);}
   return concentration;
}

float DM_mass_inside_radi(int p, float radius, int n_step, int filenr) // radius in [Mpc/h]
{
   //Concentration: log10(cvir) = 1.025 - 0.097 log10(Mvir/1e12 h^-1 M_sun)
   float x, Mhalo_radius,concentration, A_NFW;
   concentration = compute_concentration_DM_halo(p, n_step);
   x = fabs(radius)/Gal[p].Rvir; // Units = [Adimensional] [Mpc/h]/[Mpc/h]
   A_NFW = (Gal[p].Mvir) / NFW_g_function(concentration);
   Mhalo_radius = A_NFW * ( log(1+(x*concentration)) + (1./(1+(x*concentration))) -1);
   if(Mhalo_radius < 0.0) Mhalo_radius = 1e-21;
   if(isnan(Mhalo_radius)==1 || Mhalo_radius<0.){
        if(ZZ[Gal[p].SnapNum]>0.) {printf("Error in file = %d: function: DM_mass_inside_radi(): Mhalo_radius is nan! or negative, redsfhit = %f, Mhalo_radius = %E, x (radius/Rvir) = %E, radius = %E, concentration = %f\n",filenr,ZZ[Gal[p].SnapNum],Mhalo_radius,x,radius,concentration); exit(0);}}
   return Mhalo_radius;
}




float Stellar_mass_inside_radi_Bulge(int p, float radius, int filenr) // radius in [Mpc/h]
{
   // NOTE: Here we assume a Hernquist profile
   float a, Mb_r;
   int i;
   a = Gal[p].BulgeSize/(1.+pow(2.0,0.5));
   if(Gal[p].BulgeMass<1e-9 || Gal[p].BulgeSize <=0.)
        return 0.0;
   else{
          Mb_r = (Gal[p].BulgeMass/Hubble_h) * (pow(fabs(radius),2.0) / pow(fabs(radius)+a,2.0));
          if(isnan(Mb_r)==1) {printf("Error in file = %d in Stellar_mass_inside_radi_Bulge(): Mb_r == NAN (%f), BulgeMass = %f, radius = %f, a = %f, BulgeSize = %f\n",filenr,Mb_r,log10(Gal[p].BulgeMass*1e10/Hubble_h),radius,a,Gal[p].BulgeSize); exit(0);}
          if(isinf(Mb_r)==1) {printf("Error in file = %d in Stellar_mass_inside_radi_Bulge(): Mb_r == INF (%f), BulgeMass = %f, radius = %f, a = %f, BulgeSize = %f\n",filenr,Mb_r,log10(Gal[p].BulgeMass*1e10/Hubble_h),radius,a,Gal[p].BulgeSize); exit(0);}
          if(Mb_r<0.) {printf("Error in file = %d in Stellar_mass_inside_radi_Bulge(): Mb_r <0 (%f), BulgeMass = %f, radius = %f, a = %f, BulgeSize = %f\n",filenr,Mb_r,log10(Gal[p].BulgeMass*1e10/Hubble_h),radius,a,Gal[p].BulgeSize); exit(0);}
          return Mb_r;
        }
}


float mass_inside_radii_disk(int p, float radius, int filenr, int gas)
{
    //This function computes the total mass inside a radius for a expoential disk
    // For a gas disk: +1
    // For a stellar disk : -1 
    float Md_r, Md_tot, Rsl;
    if(gas > 0) {Md_tot = Gal[p].ColdGas; Rsl = Gal[p].GasDiskRadius/3.;}
    else        {Md_tot = Gal[p].DiskMass; Rsl = Gal[p].StellarDiskRadius/3.;}
    if(Md_tot<1e-9 || Rsl <=0.)
        return 0.;
    else{
         Md_r = (Md_tot/Hubble_h) * (1.- ( ( 1. + (fabs(radius)/Rsl) )*exp(-fabs(radius)/Rsl) ) );
	 // Some checks
         if(isnan(Md_r)==1) {printf("Error in file = %d in mass_inside_radii_disk(): Md_r == NAN (%f), gas = %d, radius = %f, DiskSize = %f\n",filenr,Md_r,gas,radius,Rsl); exit(0);}
         if(isinf(Md_r)==1) {printf("Error in file = %d in mass_inside_radii_disk(): Md_r == INF (%f), gas = %d, radius = %f, DiskSize = %f\n",filenr,Md_r,gas,radius,Rsl); exit(0);}
         if(Md_r<0.){printf("Error in file = %d in mass_inside_radii_disk(): Md_r < 0 (%f), gas = %d, radius = %f, DiskSize = %f\n",filenr,Md_r,gas,radius,Rsl); exit(0);}
         return Md_r;
        }
}


float NFW_density_inside_radi(int p, float radius, int n_step, int filenr) 
{
  float x, delta_c, rho_NFW, rho_radius, concentration, dt, zz;
  float zplus1, H_z0, Omega_zz;//, H_z;
  dt = ( NumToTime(Gal[p].SnapNum) -  NumToTime(Gal[p].SnapNum+1) ) / STEPS;
  //zz = redshift_at_t(ZZ[Gal[p].SnapNum],n_step,dt);
  zz = ZZ[Gal[p].SnapNum];
  concentration = compute_concentration_DM_halo(p, n_step);

  float A_NFW, rho_Halo;
  x = fabs(radius/Gal[p].Rvir); // Units = [Adimensional] [Mpc/h]/[Mpc/h]
  concentration = compute_concentration_DM_halo(p,n_step);
  A_NFW = (Gal[p].Mvir) / NFW_g_function(concentration);
  rho_Halo = Gal[p].Mvir/(4.*M_PI*pow(Gal[p].Rvir,3.0)/3.0);
  rho_radius = rho_Halo/(3.*A_NFW*x*pow(x + pow(concentration,-1.),2.0));
  return rho_radius;

}



float compute_sigma_DM_halo(int p)
{
  float sigma_DM;
  sigma_DM = pow(0.5*G*Gal[p].Mvir/Gal[p].Rvir,1./2.); // sigma_DM units: [km/s]
  if(sigma_DM<0. || isnan(sigma_DM)==1) {printf("Error in compute_sigma_DM_halo(int p): sigma_DM<0 or isnan(sigma_DM)==1 (sigma_DM = %f)",sigma_DM); exit(0);}
  return sigma_DM;

}

float density_of_ISOTHERMAL(int p, float mass, float radius, int filenr)
{
   // Density for a given simple isothermal model. It is similiar to the Equation S.4 in Henriques et al 2015
   // radius in [Mpc/h]
   float prho, Rgal;
   Rgal =((Gal[p].BulgeMass*Gal[p].BulgeSize)+(Gal[p].DiskMass*Gal[p].StellarDiskRadius/3.*1.68)+(Gal[p].ColdGas*Gal[p].GasDiskRadius/3.*1.68))/(Gal[p].ColdGas + Gal[p].DiskMass + Gal[p].BulgeMass);
   if(mass <= 1.e-9) return 0.;
   else
     {
        prho = (mass/Hubble_h)/(4 * M_PI * (Gal[p].Rvir/Hubble_h) * pow(fabs(radius)/Hubble_h,2.0));
        if(isinf(prho)){printf("Error in density_of_ISOTHERMAL() file = %d. prho is INF (%f) at radius = %f, Galaxy = %d\n", filenr,prho,radius,p); exit(0);}
        if(isnan(prho)){printf("Error in density_of_ISOTHERMAL() file = %d. prho is NAN (%f) at radius = %f, Galaxy = %d\n", filenr,prho,radius,p); exit(0);}
        if(prho<0.){printf("Error in density_of_ISOTHERMAL() file = %d. prho_ICM is negative (%f) at radius = %f, Galaxy = %d\n",filenr,prho,radius,p); exit(0);}
        return prho;
     } 
}



float incomplete_gamma(float a, float x)
{

   int CHECK = -1;
  //if(CHECK > 0) printf("a = %f, x = %f\n",a,x);

  int n;
  float sum = 1e-12;
  float next_term;
  float tg_a = 1.0;
  if (a!=1.0){ tg_a = tgamma(a);}
  for(n=0;n<50;n=n+1)
        {

          //sum += tg_a * pow(x,n)/ tgamma(a + 1.0 + n);
          next_term = tg_a * pow(x,n)/ tgamma(a + 1.0 + n);
          sum += tg_a * pow(x,n)/ tgamma(a + 1.0 + n);
          //if(isinf(tgamma(a + 1.0 + n)) == 1 || next_term < 1e-2) break;
          if(isinf(tgamma(a + 1.0 + n)) == 1 || (next_term/sum) < 1e-3 || next_term > 1e9) break;
          //if(CHECK > 0) printf("n = %f, next_term = %f, sum = %f\n",n,next_term,sum);
          //printf("n = %d, x = %f, pow(x,n) = %E, T(a+1+n)=%E\n", n, x, pow(x,n),tgamma(a + 1.0 + n));
        }
  sum =  exp(-x) * pow(x,a) * sum; // / tg_a; 
  if(isinf(sum) == 1 || isnan(sum) == 1) sum = 1e-20;

  ////U. Blahak 2010: Efficient approximation of the incomplete gamma function  , good till n<=0.25
//   int n;
//   float sum,  sumprev, gamma;
//   if (a>1.0)
//   { gamma= 1.0;
//     sumprev = 1-exp(-x); 
//     for(n=2;n<12;n=n+1)
//     { sum = (n-1)*sumprev - pow(x,n-1)*exp(-1.0*x);
//      if (n>a) break; 
//     gamma = gamma*n;
//     }
//   sum = gamma -(sum- (n-a)*(sum-sumprev));// (sumprev +(n-a)*(sum-sumprev));
//   if (sum<0.0) sum = 0.0;
//   }
//   else{
//   sum = a;
//   }

// if (CHECK>0)
//   {
//   float sum_approx = 0.0;
//   sum_approx = sum;
//   sum = 1e-12;
//   float next_term;
//   float tg_a = 1.0;
//   if (a!=1.0){ tg_a = tgamma(a);}
//   for(n=0;n<50;n=n+1)
//         {

//           //sum += tg_a * pow(x,n)/ tgamma(a + 1.0 + n);
//           next_term = tg_a * pow(x,n)/ tgamma(a + 1.0 + n);
//           sum += tg_a * pow(x,n)/ tgamma(a + 1.0 + n);
//           //if(isinf(tgamma(a + 1.0 + n)) == 1 || next_term < 1e-2) break;
//           if(isinf(tgamma(a + 1.0 + n)) == 1 || (next_term/sum) < 1e-3 || next_term > 1e9) break;
//           //if(CHECK > 0) printf("n = %f, next_term = %f, sum = %f\n",n,next_term,sum);
//           //printf("n = %d, x = %f, pow(x,n) = %E, T(a+1+n)=%E\n", n, x, pow(x,n),tgamma(a + 1.0 + n));
//         }
//   sum =  exp(-x) * pow(x,a) * sum; // / tg_a; 
//   if(isinf(sum) == 1 || isnan(sum) == 1) sum = 1e-20;
//   printf("\n Approx /gfunction computation  = %f,sum_real=%e, a=%f, x=%e, sumprev= %e ,gamma=%e",sum_approx/sum, sum, a, x, sumprev, gamma);
//   }

  //printf("Sum = %f, a = %f, x = %f\n", sum,a,x);
  //if(isinf(sum)== 1 || isnan(sum)== 1) {sum = 1.0;}
  //if(isinf(sum)== 1) {printf("Error in incomplete_gamma\n"); exit(0);}
  //exit(0);
  //printf("n = %f Integral = %E\n",n, sum);
  //if(CHECK > 0) printf("sum = %f ---------------\n",sum);
  return sum;

}


float compute_bn_Sersic(float n)
{
  float bn;
  if(n>0.36) bn = (2.0*n) - (1./3.) + (4./(405.0*n));
  else bn = 0.01945 - (0.8902*n) + (10.95*n*n) - (19.67*n*n*n) + (13.43*n*n*n*n);
  return bn;
}

float compute_p_Sersic(float n)
{
 return 1.0 - (0.6097/n) + (0.05563/pow(n,2.0));
}

float density_bulge_Sersic_model(float rho_0, float n, float Re, float r)
{
 float bn, p;
 //printf("density_bulge_Sersic_model: rho_0 = %E n = %f, Re = %E, r = %E\n",rho_0, n, Re, r);
 //exit(0);
 if(r < 0){printf("At r = 0 the Sersic density profile is not determined. Check in the code where is the r = 0\n Values = rho_0 = %E, n = %f, Re = %E, r = %E\n",rho_0,n,Re,r); exit(0);}
 //printf("rho_0 = %E, n = %f, Re = %f, r = %f\n",rho_0,n, Re,r);
 bn = compute_bn_Sersic(n);
 p = compute_p_Sersic(n);
 //printf("bn = %f , p = %f return = %f\n", bn, p, rho_0 * pow(r/Re,-p) * exp(-bn*pow(r/Re,1.0/n)));
 //printf("pow(r/Re,-p) = %f, exp(-bn*pow(r/Re,1.0/n)) = %f\n",pow(r/Re,-p) ,exp(-bn*pow(r/Re,1.0/n)));
 return rho_0 * pow((1.0*r)/(1.0*Re),-p) * exp(-1.0*bn*pow((1.0*r)/(1.0*Re),1./n));

}


float compute_rho_0_Sersic_Bulge(float Mbulge, float n, float Re)
{
 float rho_0, sum, h, bn, i, p, dummy1, dummy2;
 int stp;
 //printf("compute_rho_0_Bulge: Mbulge = %f, n = %f, Re = %f\n",Mbulge,n,Re);
 bn = compute_bn_Sersic(n);
 p = compute_p_Sersic(n);
 rho_0 = Mbulge/(4.*M_PI*pow(Re,3.0) * n * pow(bn,n*(p-3.))*tgamma(n*(3.0-p)));//gamma(n*(3.0-p))) 
 //printf("compute_rho_0_Bulge: rho_0 = %f\n",rho_0);
 return rho_0;

}


float sigma_Sersic_profile_r(float r, float n, float Re, float rho_0, float rho_r)
{
 float sum, sigma, h, bn, z, p, dummy1, dummy2, x_mn, x_mx;

  bn = compute_bn_Sersic(n);
  p = compute_p_Sersic(n);
int stp = 10;
 x_mn = bn*pow(r/Re,1.0/n);
 if(r<=Re)  x_mx = bn*pow(1.0,1.0/n);
 else
    {
      float alpha;
      alpha = pow(Re/r,-1);
      x_mx = bn*pow(5.0*alpha,1.0/n);
    }


 h = fabs(x_mx-x_mn)/(1.0*stp);
 sum = 1e-14;

 for(z=x_mn;z<=x_mx;z=z+h)
  {
    dummy1 = pow(z,  (-n*(p+1.))-1.) * exp(-z) * incomplete_gamma(n*(2.0-p),z);
    dummy2 = pow(z+h,(-n*(p+1.))-1.) * exp(-(z+h))*incomplete_gamma(n*(2.0-p),z+h);
    if(isnan(dummy1) == 1 || isnan(dummy2) == 1) {break; printf("EXXXXXXXXXXXXXXXIT\n");}
    sum += 0.5*h*(dummy1+dummy2);
  }
//Approximation , #MP
// z = bn*pow(r/Re,1.0/n);
// sum = pow(z,-1.0*n*(p+1.0)) ;// integral over angles approx unity
 // Eq. A5 paper Terzic and Graham 2005 
// sigma = pow(4.0*M_PI*G*pow(rho_0,2.0)*pow(Re,2.0)*pow(n,2.0)*pow(bn,2.0*n*(p-1.0)) * sum / density_bulge_Sersic_model(rho_0,n,Re,r) , 0.5 ); // km/s
 sigma = pow(4.0*M_PI*G_GRAVITY*pow(rho_0,2.0)*pow(Re,2.0)*pow(n,2.0)*pow(bn,2.0*n*(p-1.0)) * sum / rho_r , 0.5 ); // km/s
 //sigma = rho_0*Re*n*pow(bn,n*(p-1.0)) * pow( 4.0*M_PI*G_GRAVITY *sum/ rho_r, 0.5 ); // km/s

 return sigma;
}


float compute_influence_radius(float MBHTot, float rho_0, float n, float Re)
{
   float bn , p, z;
   float xmax, xmin, xmed, Diff1, Diff2, rinf;
   float Mbulge_xmax, Mbulge_xmed, Mbulge_10, prec;
   int Check, times, find;
   Check = -1;
   if(Check > 0) printf("------------------\n");
   bn = compute_bn_Sersic(n);
   p = compute_p_Sersic(n);

   xmax = Re;
   //xmax = 1e-7 + (3*Re-1e-7)*((double)rand()/((double)RAND_MAX+1));//3*Re;
   xmin = 1e-9;
   find = -1;
   times = 0;
   //float r;
   prec =0.03; //precision default 0.001
   //float dummy;
   //dummy = incomplete_gamma(n*(3.0-p),z);
   //return 1e-9;
   z = bn*pow(10.,1./n);
   Mbulge_10 = 4*M_PI* rho_0 * pow(Re,3.0) * n * pow(bn,n*(p-3.)) *incomplete_gamma(n*(3.0-p),z);
   if(Mbulge_10 < (2.*MBHTot)) {return 0.023*Re;} // See Fig. 2 of Biava et al. 2019
   //if(Mbulge_30 < (2.*MBHTot)) {return Re;} // See Fig. 2 of Biava et al. 2019
   //if(Mbulge_30 < (2.*MBHTot)) {return (0.0023 + (0.23-0.0023)*((double)rand()/((double)RAND_MAX+1)))*Re;} // See Fig. 2 of Biava et al. 2019
   while(find<0.)
    {
       xmed = xmin + ((xmax - xmin)/2.);
       if(Check > 0)
       {
         printf("xmax = %E (%ERe), xmin = %E (%ERe), xmed = %E(%ERe) \t",xmax, xmax/Re, xmin, xmin/Re, xmed, xmed/Re);
         //printf("xmax/Re = %E, xmin/Re = %E, xmed/Re = %E     \t",xmax/Re, xmin/Re, xmed/Re);
       }
       //printf("xmin/2. = %E \t ",xmin/2.);
       z = bn*pow(xmax/Re,1./n);
       Mbulge_xmax = 4*M_PI* rho_0 * pow(Re,3.0) * n * pow(bn,n*(p-3.))*incomplete_gamma(n*(3.0-p),z);
       z = bn*pow(xmed/Re,1./n);
       Mbulge_xmed = 4*M_PI* rho_0 * pow(Re,3.0) * n * pow(bn,n*(p-3.))*incomplete_gamma(n*(3.0-p),z);
       //z = bn*pow(5.864e-05/Re,1./n);
       //Mbulge_xmed = 4*M_PI* rho_0 * pow(Re,3.0) * n * pow(bn,n*(p-3.))*incomplete_gamma(n*(3.0-p),z);

       //printf("\n --> %E, 2.*MBHTot = %E, ratio = %f\n",Mbulge_xmed,2.*MBHTot,Mbulge_xmed/(2.*MBHTot));
        //exit(0);

       //if(Check > 0) printf("Mbulge_xmax = %E, Mbulge_xmed = %E", Mbulge_xmax, Mbulge_xmed);

       Diff1 = (Mbulge_xmax / (2.*MBHTot));
       Diff2 = (Mbulge_xmed / (2.*MBHTot));

       //printf("Mbulge_xmax = %E, Mbulge_xmed = %E", Mbulge_xmax, Mbulge_xmed);

       if(isnan(xmin) == 1 || isinf(xmin) == 1 || isnan(xmax) == 1 || isinf(xmax) == 1 || isnan(xmed) == 1 || isinf(xmed) == 1) exit(0);

       if(Check < 0)
       {
        if(Diff1 > 1 && Diff2 > 1) {xmax = xmed;}
        else if(Diff1 > 1 && Diff2 < 1) {xmin = xmed;}
        else if(Diff1 < 1 && Diff2 < 1) {xmin = xmax; xmax = 3.0*xmax;}
        else if(fabs(Diff1) >  (1.0-prec) && fabs(Diff1) < (1.0+prec)) {find = 1; rinf = xmax; break;}
        else if(fabs(Diff2) >  (1.0-prec) && fabs(Diff2) < (1.0+prec)) {find = 1; rinf = xmed; break;}
        else {
          printf("MBHTot = %f, rho_0 = %E, n = %f, Re = %E\n",MBHTot,rho_0,n,Re);
          printf("Error in compute_influence_radius(): Diff1 = %f Diff2 = %f\n",Diff1,Diff2); exit(0);}
       }
       else
       {
        printf("  Diff1 = %f, Diff2 = %f D2-1 = %E, D2-1 = %E",Diff1,Diff2, fabs(Diff1 - 1.0), fabs(Diff2 - 1.0));
        if(Diff1 > 1 && Diff2 > 1) {printf(" \t Entre 1)\n"); xmax = xmed;}
        else if(Diff1 > 1 && Diff2 < 1) {printf(" \t Entre 2)\n"); xmin = xmed;}
        else if(Diff1 < 1 && Diff2 < 1) {printf(" \t Entre 3)\n"); xmin = xmax; xmax = 2.0*xmax;}
        else {printf("Something weird is happening (MBHt = %E)\n",MBHTot); exit(0);}
       }
       if(fabs(Diff1) > (1.0-prec) && fabs(Diff1) < (1.0+prec)) {find = 1; rinf = xmax;}
       if(fabs(Diff2) >  (1.0-prec)&& fabs(Diff2) < (1.0+prec)) {find = 1; rinf = xmed;}
        times = times + 1;
        if(times > 1000) {printf("\nProbably there is a mistake in the compute_influence_radius (MBHT = %E)\n",MBHTot); 
                        printf("\n MBulge_10 = %E , n=%f, Re=%E,  Mbulge_xmed =%E ,Mbulge_xmax = %E \n",Mbulge_10, n, Re, Mbulge_xmed, Mbulge_xmax );
                        printf("Error in compute_influence_radius(): Diff1 = %f Diff2 = %f\n",Diff1,Diff2); 
                        exit(0);}
       //printf("xmax = %E, x_med = %E, Diff1 = %f, Diff2 = %f\n",xmax, xmed, Diff1,Diff2);
    }

   //if(Mbulge_30 < (2.*MBHTot))
  if(Check > 0) {printf("rinf = %E\n----\n",rinf);}
  //exit(0);
  return rinf;
}


float compute_included_mass(float r, float rho_0, float n, float Re)
{
   float bn , ps, z, M;
   //printf("\n starting M_include rinf = %e, rho_0 = %e, nb = %e , Re= %e",r,rho_0,n,Re);  
   bn = compute_bn_Sersic(n);
   ps = compute_p_Sersic(n);
   z = bn*pow(r/Re,1./n);
   M = 4*M_PI* rho_0 * pow(Re,3.0) * n * pow(bn,n*(ps-3.))* incomplete_gamma(n*(3.0-ps),z);
  //printf("\n M=%e",M);
  return M;
}

