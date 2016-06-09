#include <iostream>
#include <vector>
#include "kmc_global.h"
#include "kmc_events.h"
#include <unordered_map>

using namespace std;

double class_events::update_ratesC(int ltcp_in, bool is_recb){ // search for the ltcp input and the 1st-nn of it; remove rate info of them and rebuilt it
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

        if(cvcc.find(ltcp) != cvcc.end()){
            if(! is_recb && cvcc[ltcp].step == timestep) continue; // if has been modified, skip it

            for(int b=0; b<cvcc[ltcp].rates.size(); b ++) rate_change -= cvcc[ltcp].rates[b];
            cvcc.erase(ltcp);
        }
    
        if(*(&srf[0][0][0]+ltcp)){
            if(states[i][j][k] != 1 && states[i][j][k] != -1) error(2, "(update_ratesC) a srf atom not A nor B", 1, states[i][j][k]);
            cvcc[ltcp].step= timestep;
            
            double mu= (1==states[i][j][k]) ? muvA:muvB;
            double em= (1==states[i][j][k]) ? emvA:emvB;

            int count_AM= 0; // count how many A-M bonds.               if count > n1nbr/2, erase it.
            double rate_temp= 0; // unless meet criterion, rates_change += rate_temp
	        bool isne= false; // is neg energy
            for(int b= 0; b<n1nbr; b ++){
	            int x= pbc(i+v1nbr[b][0], nx);
	            int y= pbc(j+v1nbr[b][1], ny);
	            int z= pbc(k+v1nbr[b][2], nz);
    	        int ltcp2= x*ny*nz + y*nz + z;

                if(4==states[x][y][z]){
                    count_AM ++;

                    int MA= 0; // count how many M-A bonds if count > n1nbr/2, skip it.
                    for(int c= 0; c<n1nbr; c++){
	                    int x2= pbc(x+v1nbr[c][0], nx);
	                    int y2= pbc(y+v1nbr[c][1], ny);
	                    int z2= pbc(z+v1nbr[c][2], nz);
                        if(1==abs(states[x2][y2][z2])) MA ++;
                    }
                    if(MA > (RATIO_NOCVCC*n1nbr)+0.1) continue; // a M inside Atom matrix. SKIP

                    double e0= ecal_bond(x, y, z, i, j, k) + ecal_nonb(x, y, z, i, j, k);
				    states[x][y][z]= states[i][j][k];
                    states[i][j][k]= 0;
		            double ediff= ecal_bond(x, y, z, i, j, k) + ecal_nonb(x, y, z, i, j, k) - e0; 
				    states[i][j][k]= states[x][y][z];
				    states[x][y][z]= 4;
                    if(ediff<0) isne= true;
                    
				    cvcc[ltcp].mltcp.push_back(ltcp2);
                    cvcc[ltcp].rates.push_back(mu * exp(-beta*(ediff+em)));
				
				    rate_temp += cvcc[ltcp].rates.back();
                }
		    }

            if((count_AM > (RATIO_NOCVCC*n1nbr)+0.1) || 0==cvcc[ltcp].mltcp.size()) cvcc.erase(ltcp);
            else if(isne) error(2, "(update_ratesC) e is neg but not srf atom jump case. nAM", 1, count_AM);
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
