#include <cstdio>
#include <iostream>
#include "kmc_global.h"
#include "kmc_events.h"

using namespace std;

double class_events::ecal_otf(bool is_itl, int x1, int y1, int z1, int x2, int y2, int z2) const{ // corrected H on the fly 
	int pos[2][3]= {{x1, y1, z1}, {x2, y2, z2}};

	// AB-A, AB-B bonds, types of the 2 atoms (itl only)
	int N_corr1[5]={0}; // []: 0:BB, 1:B, 2:ABorV, 3:A, 4:AA
	int N_corr2[5]={0}; // 2nd nn

	// sums of 1st nn
	int s1um44=0, s1um42=0, s1um22=0;		// class 1
	int s1um43=0, s1um41=0, s1um32=0, s1um21=0;	// class 2
	int s1um33=0, s1um31=0, s1um11=0;		// class 3
	int s1um40=0, s1um30=0, s1um20=0, s1um10=0, s1um00=0;
	// sums of 2nd nn
	int s2um44=0, s2um42=0, s2um22=0;		// class 1
	int s2um43=0, s2um41=0, s2um32=0, s2um21=0;	// class 2
	int s2um33=0, s2um31=0, s2um11=0;		// class 3
	int s2um40=0, s2um30=0, s2um20=0, s2um10=0, s2um00=0;
	for(int i=0; i<2; i ++){
		int xi= pbc(pos[i][0], nx);
		int yi= pbc(pos[i][1], ny);
		int zi= pbc(pos[i][2], nz);

		int state_i= states[xi][yi][zi];

		for(int a=0; a<n1nbr; a ++){ // 1st neighbors
			int xj= pbc(xi+(*(v1nbr+a))[0], nx);
			int yj= pbc(yi+(*(v1nbr+a))[1], ny);
			int zj= pbc(zi+(*(v1nbr+a))[2], nz);

			if((xj==x1) && (yj==y1) && (zj==z1)) continue; // avoid double count the (x1, y1, z1)-(x2, y2, z2) bond

			int state_j1= states[xj][yj][zj];
			if(0==state_i)		N_corr1[state_j1+2] ++;
			else if(0==state_j1)	N_corr1[state_i +2] ++;

			// calculating the energy for ABVI system. 1st nn
			s1um44 += powc(state_i, 4) * powc(state_j1, 4);
			s1um43 += powc(state_i, 4) * powc(state_j1, 3) + powc(state_i, 3) * powc(state_j1, 4);
			s1um42 += powc(state_i, 4) * powc(state_j1, 2) + powc(state_i, 2) * powc(state_j1, 4);
			s1um41 += powc(state_i, 4) * powc(state_j1, 1) + powc(state_i, 1) * powc(state_j1, 4);
			s1um33 += powc(state_i, 3) * powc(state_j1, 3);
			s1um32 += powc(state_i, 3) * powc(state_j1, 2) + powc(state_i, 2) * powc(state_j1, 3);
			s1um31 += powc(state_i, 3) * powc(state_j1, 1) + powc(state_i, 1) * powc(state_j1, 3);
			s1um22 += powc(state_i, 2) * powc(state_j1, 2);
			s1um21 += powc(state_i, 2) * powc(state_j1, 1) + powc(state_i, 1) * powc(state_j1, 2);
			s1um11 += powc(state_i, 1) * powc(state_j1, 1);
			
			s1um40 += powc(state_i, 4) + powc(state_j1, 4);
			s1um30 += powc(state_i, 3) + powc(state_j1, 3);
			s1um20 += powc(state_i, 2) + powc(state_j1, 2);
			s1um10 +=      state_i     +      state_j1    ;
			s1um00 ++;
		}
				
		if(! is_e2nbr) continue;

		for(int b=0; b<n2nbr; b ++){ // 2nd neighbors
			int xj= pbc(xi+(*(v2nbr+b))[0], nx);
			int yj= pbc(yi+(*(v2nbr+b))[1], ny);
			int zj= pbc(zi+(*(v2nbr+b))[2], nz);

			int state_j2= states[xj][yj][zj];
			if(0==state_i)		N_corr2[state_j2+2] ++;
			else if(0==state_j2)	N_corr2[state_i +2] ++;

			// calculating the energy: 2nd nn
			s2um44 += powc(state_i, 4) * powc(state_j2, 4);
			s2um43 += powc(state_i, 4) * powc(state_j2, 3) + powc(state_i, 3) * powc(state_j2, 4);
			s2um42 += powc(state_i, 4) * powc(state_j2, 2) + powc(state_i, 2) * powc(state_j2, 4);
			s2um41 += powc(state_i, 4) * powc(state_j2, 1) + powc(state_i, 1) * powc(state_j2, 4);
			s2um33 += powc(state_i, 3) * powc(state_j2, 3);
			s2um32 += powc(state_i, 3) * powc(state_j2, 2) + powc(state_i, 2) * powc(state_j2, 3);
			s2um31 += powc(state_i, 3) * powc(state_j2, 1) + powc(state_i, 1) * powc(state_j2, 3);
			s2um22 += powc(state_i, 2) * powc(state_j2, 2);
			s2um21 += powc(state_i, 2) * powc(state_j2, 1) + powc(state_i, 1) * powc(state_j2, 2);
			s2um11 += powc(state_i, 1) * powc(state_j2, 1);
			
			s2um40 += powc(state_i, 4) + powc(state_j2, 4);
			s2um30 += powc(state_i, 3) + powc(state_j2, 3);
			s2um20 += powc(state_i, 2) + powc(state_j2, 2);
			s2um10 +=      state_i     +      state_j2    ;
			s2um00 ++;
		}
	}

	double sum_energy= 
		unc1_44 * s1um44 + unc1_43 * s1um43 + unc1_42 * s1um42 + unc1_41 * s1um41 +
		unc1_33 * s1um33 + unc1_32 * s1um32 + unc1_31 * s1um31 +
		unc1_22 * s1um22 + unc1_21 * s1um21 +
		unc1_11 * s1um11 +
		unc1_40 * s1um40 + unc1_30 * s1um30 + unc1_20 * s1um20 + unc1_10 * s1um10 + unc1_00 * s1um00 +

		unc2_44 * s2um44 + unc2_43 * s2um43 + unc2_42 * s2um42 + unc2_41 * s2um41 +
		unc2_33 * s2um33 + unc2_32 * s2um32 + unc2_31 * s2um31 +
		unc2_22 * s2um22 + unc2_21 * s2um21 +
		unc2_11 * s2um11 +
		unc2_40 * s2um40 + unc2_30 * s2um30 + unc2_20 * s2um20 + unc2_10 * s2um10 + unc2_00 * s2um00;

	if(is_itl) sum_energy -= (N_corr1[3]*eA1V  + N_corr1[2]*eV1V   + N_corr1[1]*eV1B)  + (N_corr2[3]*eA2V  + N_corr2[2]*eV2V   + N_corr2[1]*eV2B);
	else	   sum_energy -= (N_corr1[3]*eA1AB + N_corr1[2]*eAB1AB + N_corr1[1]*eAB1B) + (N_corr2[3]*eA2AB + N_corr2[2]*eAB2AB + N_corr2[1]*eAB2B);

	return sum_energy;
}

