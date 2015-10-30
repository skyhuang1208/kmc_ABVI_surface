#include <iostream>
#include <vector>
#include "kmc_global.h"
#include "kmc_events.h"

using namespace std;

double class_events::cal_ratesC(vector <int> &etype, vector <double> &rates, vector <int> &ilist, vector <int> &nltcp, vector <int> &jatom){
	double sum_rate= 0;
	
	for(int a=0; a < nx*ny*nz; a ++){
		if(*(&srf[0][0][0]+a)){
            int i= (int) (a/nz)/ny;
		    int j= (int) (a/nz)%ny;
		    int k= (int)  a%nz;

            double em, mu;
			switch(states[i][j][k]){
                case  1: em= emvA; mu= muvA; break; 
                case -1: em= emvB; mu= muvB; break; 
                default: error(2, "(cal_ratesC) an srf atom isnt 1 or -1", 2, states[i][j][k], timestep);
            }

		    for(int b=0; b<n1nbr; b ++){
			    int x= pbc(i+v1nbr[b][0], nx);
			    int y= pbc(j+v1nbr[b][1], ny);
			    int z= pbc(k+v1nbr[b][2], nz);

			    if(4==states[x][y][z]){
				    double e0= cal_energy(false, i, j, k, x, y, z);

				    states[x][y][z]= states[i][j][k];
                    states[i][j][k]= 0;

				    double ediff= cal_energy(false, i, j, k, x, y, z) - e0;
                    if(ediff<0) error(2, "(cal_rateC) minus ediff in creation of vcc", 1, ediff);

				    states[i][j][k]= states[x][y][z];
				    states[x][y][z]= 4;
			
                    rates.push_back(mu * exp(-beta*(em+0.5*ediff)));
							
				    etype.push_back(8);
				    ilist.push_back(a);
				    nltcp.push_back(x*ny*nz+y*nz+z);
				    jatom.push_back(states[i][j][k]);
				
				    sum_rate += rates.back();
                }
            }
		}
	}

	return sum_rate;
}
