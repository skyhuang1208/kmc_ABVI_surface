#include <iostream>
#include <vector>
#include "kmc_global.h"
#include "kmc_events.h"
#include <unordered_map>

using namespace std;

double class_events::update_ratesC(int ltcp_in, bool is_recb){ // remove rate info of them and rebuilt it
	double rate_change= 0;

    int i0= (int) (ltcp_in/nz)/ny;
    int j0= (int) (ltcp_in/nz)%ny;
    int k0= (int)  ltcp_in%nz;

	for(int a= -1; a<n1nbr; a ++){ // loop tho itself & nbrs
        int i, j, k, ltcp;
        if(-1==a){ i= i0; j= j0; k= k0; ltcp= ltcp_in;}
        else{
	        i= pbc(i0+v1nbr[a][0], nx);
	        j= pbc(j0+v1nbr[a][1], ny);
	        k= pbc(k0+v1nbr[a][2], nz);
    	    ltcp= i*ny*nz + j*nz + k;
        }

        if(cvcc.find(ltcp) != cvcc.end()){
            if(! is_recb && cvcc[ltcp].step == timestep) continue; // if has been modified, dont do again

            for(int b=0; b<cvcc[ltcp].rates.size(); b ++) rate_change -= cvcc[ltcp].rates[b];
            cvcc.erase(ltcp); // re-calculate it
        }
    
        if(srf[i][j][k]){
            if(states[i][j][k] != 1 && states[i][j][k] != -1) error(2, "(update_ratesC) a srf atom not A nor B", 1, states[i][j][k]);
            cvcc[ltcp].step= timestep;
            
            double mu= (1==states[i][j][k]) ? muvA:muvB;
            double em= (1==states[i][j][k]) ? emvA:emvB;

            double rate_temp= 0; // temporary rate
	        bool isne= false;    // is neg energy
            int AM= 0;           // count of A-M bonds
            for(int b= 0; b<n1nbr; b ++){
	            int x= pbc(i+v1nbr[b][0], nx);
	            int y= pbc(j+v1nbr[b][1], ny);
	            int z= pbc(k+v1nbr[b][2], nz);
    	        int ltcp2= x*ny*nz + y*nz + z;

                if(4==states[x][y][z]){
                    AM ++;

                    states[x][y][z]= 0; // considering as V+atom->atom+V, not M+atom->atom+V
                    double e0= ecal_bond(x, y, z, i, j, k) + ecal_nonb(x, y, z, i, j, k);
				    states[x][y][z]= states[i][j][k];
                    states[i][j][k]= 0;
		            double ediff= ecal_bond(x, y, z, i, j, k) + ecal_nonb(x, y, z, i, j, k) - e0; 
				    states[i][j][k]= states[x][y][z];
				    states[x][y][z]= 4;
                    
                    ediff += em; // ediff becomes ediff+em
                    if(ediff<0){
                        isne= true;
                        ediff= 0;
                    }
                    
				    cvcc[ltcp].mltcp.push_back(ltcp2);
                    cvcc[ltcp].rates.push_back(mu * exp(-beta*ediff));
				
				    rate_temp += cvcc[ltcp].rates.back();
                }
		    }

            if(0==AM) error(2, "(cvcc) an srf atom does not bond with vacuum");
            else if(AM >= N_NOCVCC) cvcc.erase(ltcp); // either extrusion or intrusion
            else{
                rate_change += rate_temp;
                if(isne) cout << "*** neg e in cvcc ***";
            }
        }
    }

	return rate_change;
}

double class_events::init_ratesC(){
	double sum_rate= 0;

	for(int i=0; i < nx; i ++){
	    for(int j=0; j < ny; j ++){
	        for(int k=0; k < nz; k ++){
		        if(srf[i][j][k]) sum_rate += update_ratesC(i*ny*nz+j*nz+k);
    }}}

	return sum_rate;
}
