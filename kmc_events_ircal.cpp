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
		if(list_itl[ii].trapped) continue; // if trapped; cant move; skip
        
		int i= list_itl[ii].x;
		int j= list_itl[ii].y;
		int k= list_itl[ii].z;
		int stateI= states[i][j][k]; // state of the itl

        if(2==stateI){
            int dir1= list_itl[ii].dir; // neighbor index
            int dir2= (dir1<n1nbr/2) ? (dir1+n1nbr/2):(dir1-n1nbr/2);
		    
            for(int a=0; a<n1nbr; a ++){
		    	int x= pbc(i+v1nbr[a][0], nx);
		    	int y= pbc(j+v1nbr[a][1], ny);
		    	int z= pbc(k+v1nbr[a][2], nz);
		        int stateA= states[x][y][z]; // state of the itl
                if(-1==stateA) error(2, "(cal_rateI) an B atom surrounded by AA itl");
                if(stateA != 1) continue;
    	
                double e0= ecal_bond(i, j, k, x, y, z) + ecal_nonb(i, j, k, x, y, z);
                states[i][j][k]= 1; 
                states[x][y][z]= 2;
		        double ediff= ecal_bond(i, j, k, x, y, z) + ecal_nonb(i, j, k, x, y, z) - e0; 
                states[i][j][k]= 2; 
                states[x][y][z]= 1;
               
                double cf= 1.0; // correlation factor (only 1/8 has not-1.0 value)
                double e;
                if(a==dir2)             cf= corrfac;
                if(a==dir1 || a==dir2)  e= emiAA;
                else                    e= emiAA+erAA;
                if(ediff>0)             e+= ediff;
                rates.push_back(muiAA*exp(-beta*e)/cf); // rate_back= rate_for/cf
                
                sum_rate += rates.back();
                etype.push_back(0);
   		        ilist.push_back(ii);
    	        inbr.push_back(a);
            }
        }
        else{ // AB itl
            for(int a=0; a<n2nbr; a ++){
		    	int x= pbc(i+v2nbr[a][0], nx);
		    	int y= pbc(j+v2nbr[a][1], ny);
		    	int z= pbc(k+v2nbr[a][2], nz);
		        int stateA= states[x][y][z]; // state of the itl
                if(stateA != 1) continue;
    	
                double e0= ecal_bond(i, j, k, x, y, z) + ecal_nonb(i, j, k, x, y, z);
                states[i][j][k]= 1; 
                states[x][y][z]= 3;
		        double ediff= ecal_bond(i, j, k, x, y, z) + ecal_nonb(i, j, k, x, y, z) - e0; 
                states[i][j][k]= 3;
                states[x][y][z]= 1;
                
                double e;
                if(ediff<0) e= emiAB;
                else        e= emiAB+ediff;
                rates.push_back(muiAB*exp(-beta*e));
                
                sum_rate += rates.back();
                etype.push_back(0);
   		        ilist.push_back(ii);
    	        inbr.push_back(a);
            }
        }
    }
        
	return sum_rate;
}
