#include <iostream>
#include <cmath>
#include <vector>
#include "kmc_global.h"
#include "kmc_events.h"

using namespace std;

double class_events::cal_ratesI(vector <int> &etype, vector <double> &rates, vector <int> &ilist, vector <int> &inbr){
	double sum_rate= 0;
	if(nAA + nAB + nBB != list_itl.size()) error(2, "(cal_ratesI) itl number inconsistent");
	
	for(int ii=0; ii < list_itl.size(); ii ++){ // ii: index of interstitial
		int ltcp= list_itl[ii].ltcp;
		int i= (int) (ltcp/nz)/ny;
		int j= (int) (ltcp/nz)%ny;
		int k= (int)  ltcp%nz;
		int stateI= states[i][j][k]; // state of the itl

        int dir= list_itl[ii].dir; // neighbor index
		    
		int x= pbc(i+v1nbr[dir][0], nx);
		int y= pbc(j+v1nbr[dir][1], ny);
		int z= pbc(k+v1nbr[dir][2], nz);
		int stateA= states[x][y][z]; // state of the itl
        if(stateA != 1) continue;
    	
        rates.push_back(muiAA*exp(-beta*emiAA));
                
        sum_rate += rates.back();
        etype.push_back(0);
   		ilist.push_back(ii);
    	inbr.push_back(dir);
    }
        
	return sum_rate;
}
