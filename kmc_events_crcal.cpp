#include <iostream>
#include <vector>
#include "kmc_global.h"
#include "kmc_events.h"

using namespace std;

double class_events::update_ratesC(int ltcp_in){ // search for the ltcp input and the 1st-nn of it; remove rate info of them and rebuilt it
	double rate_change= 0;

    int i0= (int) (ltcp_in/nz)/ny;
    int j0= (int) (ltcp_in/nz)%ny;
    int k0= (int)  ltcp_in%nz;

	for(int a= -1; a<n1nbr; a ++){
        int i, j, k, ltcp;
        if(-1==a){ i= i0; j= j0; k= k0; ltcp= ltcp_in;}
        else{
	        i= pbc(i0+v1nbr[a][0], nx);
	        j= pbc(j0+v1nbr[a][1], ny);
	        k= pbc(k0+v1nbr[a][2], nz);
    	    ltcp= i*ny*nz + j*nz + k;
        }

        if(1==cvcc.count(ltcp)){
            if(cvcc[ltcp].step == timestep) continue; // if has been modified, skip it

            for(int b=0; b<cvcc[ltcp].rates.size(); b ++) rate_change -= cvcc[ltcp].rates[b];
                
            if(*(&srf[0][0][0]+ltcp)){
                cvcc[ltcp].mltcp.clear();
                cvcc[ltcp].rates.clear();
            }
            else cvcc.erase(ltcp);
        }
        else if(*(&srf[0][0][0]+ltcp)) cvcc[ltcp]= cvcc_info(); 
    
        if(*(&srf[0][0][0]+ltcp)){
            cvcc[ltcp].step= timestep;
            
            double em, mu;
			switch(states[i][j][k]){
                case  1: em= emvA; mu= muvA; break; 
                case -1: em= emvB; mu= muvB; break; 
                default: error(2, "(cal_ratesC) an srf atom isnt 1 or -1", 2, states[i][j][k], timestep);
            }

            int count_AM= 0; // count how many A-M bonds.               if count > n1nbr/2, erase it.
            double rate_temp= 0; // unless meet criterion, rates_change += rate_temp
	        for(int b= 0; b<n1nbr; b ++){
	            int x= pbc(i+v1nbr[b][0], nx);
	            int y= pbc(j+v1nbr[b][1], ny);
	            int z= pbc(k+v1nbr[b][2], nz);
    	        int ltcp2= x*ny*nz + y*nz + z;

                if(4==states[x][y][z]){
                    count_AM ++;

                    double e0= ecal_bond(i, j, k, x, y, z) + ecal_nonb(i, j, k, x, y, z);

			        states[x][y][z]= states[i][j][k];
                    states[i][j][k]= 0;

		            double ediff= ecal_bond(i, j, k, x, y, z) + ecal_nonb(i, j, k, x, y, z) - e0; 

				    states[i][j][k]= states[x][y][z];
				    states[x][y][z]= 4;
			
				    cvcc[ltcp].mltcp.push_back(ltcp2);
                    cvcc[ltcp].rates.push_back(mu * exp(-beta*(em+0.5*ediff)));
                    if((em+0.5*ediff)<0) error(2, "(update_ratesC) minus e in creation of vcc", 1, em+0.5*ediff); // delete it
				
				    rate_temp += cvcc[ltcp].rates.back();
                }
		    }

            if(double( count_AM) > RATIO_NOCVCC * n1nbr) cvcc.erase(ltcp);
            else                                         rate_change += rate_temp;
	    }
    }

	return rate_change;
}

double class_events::init_ratesC(){
	double sum_rate= 0;

	for(int a=0; a < nx*ny*nz; a ++)
		if(*(&srf[0][0][0]+a)) sum_rate += update_ratesC(a);

	return sum_rate;
}
