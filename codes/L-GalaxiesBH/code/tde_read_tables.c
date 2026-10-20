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

/*
 *  Created in: 2022
 *  #MP
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

void read_TDErates_tables(void)
{	
	static char *name1 = "mbh_m0.txt";
	static char *name2 = "times.txt";
	//------------------------------------------
	//READ MASSES OF BLACKHOLES  (LOGSCALE)
	//------------------------------------------
	FILE *fd1;
	char buf1[100];
	int i1;
	float m1;

	sprintf(buf1, "./TDErateTables/%s", name1);

	if(!(fd1 = fopen(buf1, "r")))
        {
          printf("file `%s' not found.\n", buf1);
          exit(0);
        }

	for(i1=0; i1< NUM_MBH_BINS ; i1++)
        {
	  fscanf(fd1, "%f", &m1);
	   TDEmbh[i1] = m1*1.0e-10/Hubble_h; //10^10 Msun h^-1 to compare with BH mass
        }
	fclose(fd1);
	//printf("Lifetime masses read.\n");

	//------------------------------------------
	//READ TIME OF SIMULATION (LOGSCALE) 
	//------------------------------------------
	FILE *fd2;
	char buf2[100];
	int i2;
	float m2;

	sprintf(buf2, "./TDErateTables/%s", name2);

	if(!(fd2 = fopen(buf2, "r")))
        {
          printf("file `%s' not found.\n", buf2);
          exit(0);
        }

	for(i2=0; i2<NUM_LOGTIME_BINS; i2++)
        {
	  fscanf(fd2, "%f", &m2);
	  TDEtimes[i2] = (1.4e10- m2)* 1.0e-6 *MYRs_TO_INTERNAL; //iternal time, little h included in units, also reversed to lookback time
        }
	fclose(fd2);
	//printf("Lifetime metallicities read.\n");



	//------------------------------------------
	//READ AGB YIELD TABLES:
	//------------------------------------------
	FILE *fd;
	char buf[100];
	int i,j ,ns;
	float m;

    for (ns=1;ns<7;ns++)
	{
	  sprintf(buf, "./TDErateTables/rates_ind_%d.txt",ns); 

	  	if(!(fd = fopen(buf, "r")))
          {
            printf("file `%s' not found, check for theSersic index TDErates.\n", buf);
			printf("ns = %d",ns);
            exit(0);
          }
	  for(i=0; i<NUM_MBH_BINS; i++) //Number of element species (inc H and He) = 11
          {
            for(j=0; j<NUM_LOGTIME_BINS; j++)
            {
            	fscanf(fd, "%f", &m);
            	TDErates[ns][i][j] = m * 1e6 /MYRs_TO_INTERNAL; //internal time, little h included in units
				//if (m>1.0) printf("\n At (  nsersic, i , j) = (%d, %d, %d) value seems high : %f", ns, i ,j, TDErates[ns][i][j]);
            }
          }
	}

	printf("TDE Rates tables read.\n");
	//------------------------------------------
}
