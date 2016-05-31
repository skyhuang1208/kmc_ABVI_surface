#include <iostream>
#include <vector>
#include "kmc_par.h"
#include "kmc_global.h"
#include "kmc_events.h"

using namespace std;

double class_events::cal_ratesVsp(vector <int> &etype, vector <double> &rates, vector <int> &ilist, vector <int> &inbr){
	double sum_rate= 0;
	if(nV != list_vcc.size()) error(2, "(cal_ratesV) vcc number inconsistent", 2, nV, list_vcc.size());

	for(int ivcc=0; ivcc < nV; ivcc ++){
		int i= (int) (list_vcc[ivcc].ltcp/nz)/ny;
		int j= (int) (list_vcc[ivcc].ltcp/nz)%ny;
		int k= (int)  list_vcc[ivcc].ltcp%nz;

		if(states[i][j][k] != 0) error(2, "(cal_ratesV) there's an non-vacancy in the vacancy list", 2, states[i][j][k], timestep);
		
		for(int a=0; a<n1nbr; a ++){
			int x= pbc(i+v1nbr[a][0], nx);
			int y= pbc(j+v1nbr[a][1], ny);
			int z= pbc(k+v1nbr[a][2], nz);
		
			if(1==states[x][y][z] || -1==states[x][y][z]){
                double ediff= ecal_sp(states[x][y][z], a, i, j, k) - ecal_bond(i, j, k, x, y, z);

                double e0_nonb= ecal_nonb(i, j, k, x, y, z); // non-broken e0
				states[i][j][k]= states[x][y][z];
                if(srf[x][y][z]) states[x][y][z]= 4;
                else             states[x][y][z]= 0;
                ediff += ecal_nonb(i, j, k, x, y, z) - e0_nonb; // non-broken e1
				states[x][y][z]= states[i][j][k];
				states[i][j][k]= 0;

                if((ediff)<0) error(2, "(cal_rateV) got minus sign ", 1, ediff);
                
                double mu= (1==states[x][y][z]) ? muvA:muvB;
                rates.push_back(mu * exp(-beta*ediff));
							
				etype.push_back(1);
				ilist.push_back(ivcc);
				inbr.push_back(a);

				sum_rate += rates.back();
			}
		}
	}

	return sum_rate;
}

double class_events::ecal_sp(int stateA1, int inbr, int i, int j, int k) const{ // calculate saddle-point e 
    if(stateA1 != 1 && stateA1 != -1) error(2, "(ecal_sp) stateA1 is not an atom (type)", 1, stateA1);
    int Nbonds1[3][3]= {0}; // should be {A or B}{A or B}
    int Nbonds2[3][3]= {0}; // should be {A or B}{A or B}
    
    for(int a=0; a<n1sp; a ++){ // 1st SP nbrs
	    int x= pbc(i+v1sp[inbr][a][0], nx);
		int y= pbc(j+v1sp[inbr][a][1], ny);
		int z= pbc(k+v1sp[inbr][a][2], nz);
		int stateA2= states[x][y][z];
        if(abs(stateA2)>1) error(2, "(ecal_sp) stateA2 is an itl", 1, stateA2);
        
        Nbonds1[stateA1+1][stateA2+1] ++;
    }
    
    for(int a=0; a<n2sp; a ++){ // 1st SP nbrs
	    int x= pbc(i+v2sp[inbr][a][0], nx);
		int y= pbc(j+v2sp[inbr][a][1], ny);
		int z= pbc(k+v2sp[inbr][a][2], nz);
		int stateA2= states[x][y][z];
        if(abs(stateA2)>1) error(2, "(ecal_sp) stateA2 is an itl", 1, stateA2);
		
        Nbonds2[stateA1+1][stateA2+1] ++;
    }

    return Nbonds1[2][2]*eSPA1A + Nbonds1[2][1]*eSPA1V + Nbonds1[2][0]*eSPA1B + Nbonds1[0][2]*eSPB1A + Nbonds1[0][1]*eSPB1V + Nbonds1[0][0]*eSPB1B +
           Nbonds2[2][2]*eSPA2A + Nbonds2[2][1]*eSPA2V + Nbonds2[2][0]*eSPA2B + Nbonds2[0][2]*eSPB2A + Nbonds2[0][1]*eSPB2V + Nbonds2[0][0]*eSPB2B;
}
