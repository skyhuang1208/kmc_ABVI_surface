#include <cstdio>
#include <iostream>
#include "kmc_par.h"
#include "kmc_global.h"
#include "kmc_events.h"

using namespace std;

double class_events::ecal_bond(int x1, int y1, int z1, int x2, int y2, int z2) const{ // calculate the energy of the input 2 ltcps 
	int pos[2][3]= {{x1, y1, z1}, {x2, y2, z2}};
	int bonds1[7][7]= {0}; // all shifted +2, e.g., AA: 4, BB: 0, AB: 5
	int bonds2[7][7]= {0}; // 2nd nn
    vector <int> AB1; // AB for diff B concentrations for 1st-nn
    vector <int> AB2; // AB for diff B concentrations for 2nd-nn
    vector <int> BV1; // BV (1st-nn)
    vector <int> BV2; // BV (2nd-nn)
	
    for(int i=0; i<2; i ++){
		int xi= pbc(pos[i][0], nx);
		int yi= pbc(pos[i][1], ny);
		int zi= pbc(pos[i][2], nz);
        int nBnbr_i= 0;
		int state_i= states[xi][yi][zi];

        if(4==state_i) continue;

		for(int a=0; a<n1nbr; a ++){ // 1st neighbors
			int xj= pbc(xi+(*(v1nbr+a))[0], nx);
			int yj= pbc(yi+(*(v1nbr+a))[1], ny);
			int zj= pbc(zi+(*(v1nbr+a))[2], nz);
		    int nBnbr_j1= 0;
			int state_j1= states[xj][yj][zj];
		
			if((xj==x1) && (yj==y1) && (zj==z1)) continue; // avoid double count the (x1, y1, z1)-(x2, y2, z2) bond
            if(4==state_j1) continue; // bond with vacuum equals to 0

			bonds1[state_i+2][state_j1+2] ++;
                    
            if((-1==state_i && 1==state_j1) || (1==state_i && -1==state_j1)){
                nBnbr_i = cal_Bnbr(nBnbr_i, xi, yi, zi);
                nBnbr_j1= cal_Bnbr(nBnbr_j1, xj, yj, zj);
                AB1.push_back(nBnbr_i);
                AB1.push_back(nBnbr_j1);
            }
            if((-1==state_i && 0==state_j1) || (0==state_i && -1==state_j1)){
                nBnbr_i = cal_Bnbr(nBnbr_i, xi, yi, zi);
                nBnbr_j1= cal_Bnbr(nBnbr_j1, xj, yj, zj);
                BV1.push_back(nBnbr_i);
                BV1.push_back(nBnbr_j1);
            }
		}
	
		if(! is_e2nbr) continue;

		for(int a=0; a<n2nbr; a ++){ // 2nd neighbors
			int xj= pbc(xi+(*(v2nbr+a))[0], nx);
			int yj= pbc(yi+(*(v2nbr+a))[1], ny);
			int zj= pbc(zi+(*(v2nbr+a))[2], nz);
		    int nBnbr_j2= 0;
			int state_j2= states[xj][yj][zj];
		
			if((xj==x1) && (yj==y1) && (zj==z1)) continue; // should not happen in 2nd neighbors
            if(4==state_j2) continue; // bond with vacuum equals to 0

			bonds2[state_i+2][state_j2+2] ++;
            
            if((-1==state_i && 1==state_j2) || (1==state_i && -1==state_j2)){
                nBnbr_i = cal_Bnbr(nBnbr_i, xi, yi, zi);
                nBnbr_j2= cal_Bnbr(nBnbr_j2, xj, yj, zj);
                AB2.push_back(nBnbr_i);
                AB2.push_back(nBnbr_j2);
            }
            if((-1==state_i && 0==state_j2) || (0==state_i && -1==state_j2)){
                nBnbr_i = cal_Bnbr(nBnbr_i, xi, yi, zi);
                nBnbr_j2= cal_Bnbr(nBnbr_j2, xj, yj, zj);
                BV2.push_back(nBnbr_i);
                BV2.push_back(nBnbr_j2);
            }
		}
	}
  
    const int locals= n1nbr + n2nbr +1; // local sites
    double e= 0;
    for(int i=0; i < AB1.size(); i ++) e += e0A1B + e1A1B * (AB1[i]*1.0/locals);
    for(int i=0; i < AB2.size(); i ++) e += e0A2B + e1A2B * (AB2[i]*1.0/locals);
    for(int i=0; i < BV1.size(); i ++) e += e0B1V + e1B1V / (BV1[i]*1.0/locals); // !! different formula
    for(int i=0; i < BV2.size(); i ++) e += e0B2V + e1B2V * (BV2[i]*1.0/locals);
    e /= 2; // a bond energy is avg of eAB of 2 atoms

    e += 
    eM1M * bonds1[6][6]                                                 + eM1A  * bonds1[6][3] + eM1V * bonds1[6][2] + eM1B  * bonds1[6][1]                         + 
                          eAB1AB * bonds1[5][5] + eAA1AB * bonds1[5][4] + eA1AB * bonds1[5][3]                       + eAB1B * bonds1[5][1] + eAB1BB * bonds1[5][0] +
	                      eAA1AB * bonds1[4][5] + eAA1AA * bonds1[4][4] + eAA1A * bonds1[4][3]                       + eAA1B * bonds1[4][1] + eAA1BB * bonds1[4][0] +
	eM1A * bonds1[3][6] + eA1AB  * bonds1[3][5] + eAA1A  * bonds1[3][4] + eA1A  * bonds1[3][3] + eA1V * bonds1[3][2] + 0                    + eA1BB  * bonds1[3][0] +
	eM1V * bonds1[2][6]                                                 + eA1V  * bonds1[2][3] + eV1V * bonds1[2][2] + 0                                            +
	eM1B * bonds1[1][6] + eAB1B  * bonds1[1][5] + eAA1B  * bonds1[1][4] + 0                    + 0                   + eB1B  * bonds1[1][1] + eB1BB  * bonds1[1][0] +
	                      eAB1BB * bonds1[0][5] + eAA1BB * bonds1[0][4] + eA1BB * bonds1[0][3]                       + eB1BB * bonds1[0][1] + eBB1BB * bonds1[0][0] +
// 2nd-nn
    eM2M * bonds2[6][6]                                                 + eM2A  * bonds2[6][3] + eM2V * bonds2[6][2] + eM2B  * bonds2[6][1]                         + 
	                      eAB2AB * bonds2[5][5] + eAA2AB * bonds2[5][4] + eA2AB * bonds2[5][3]                       + eAB2B * bonds2[5][1] + eAB2BB * bonds2[5][0] +
	                      eAA2AB * bonds2[4][5] + eAA2AA * bonds2[4][4] + eAA2A * bonds2[4][3]                       + eAA2B * bonds2[4][1] + eAA2BB * bonds2[4][0] +
	eM2A * bonds2[3][6] + eA2AB  * bonds2[3][5] + eAA2A  * bonds2[3][4] + eA2A  * bonds2[3][3] + eA2V * bonds2[3][2] + 0                    + eA2BB  * bonds2[3][0] +
	eM2V * bonds2[2][6]                                                 + eA2V  * bonds2[2][3] + eV2V * bonds2[2][2] + 0                                            +
	eM2B * bonds2[1][6] + eAB2B  * bonds2[1][5] + eAA2B  * bonds2[1][4] + 0                    + 0                   + eB2B  * bonds2[1][1] + eB2BB  * bonds2[1][0] +
	                      eAB2BB * bonds2[0][5] + eAA2BB * bonds2[0][4] + eA2BB * bonds2[0][3]                       + eB2BB * bonds2[0][1] + eBB2BB * bonds2[0][0];
	
	return e;
}

double class_events::ecal_nonb(int x1, int y1, int z1, int x2, int y2, int z2) const{ // cal non-broken A-B bonds 
    int pos[2][3]= {{x1, y1, z1}, {x2, y2, z2}};
    vector <int> AB1; // AB for diff B concentrations for 1st-nn
    vector <int> AB2; // AB for diff B concentrations for 2nd-nn
    vector <int> BV1; // BV (1st-nn)
    vector <int> BV2; // BV (2nd-nn)
	
    for(int o=0; o<2; o ++){
		int x0= pbc(pos[o][0], nx);
		int y0= pbc(pos[o][1], ny);
		int z0= pbc(pos[o][2], nz);

		for(int i=0; i<n12nbr; i ++){ // 1st neighbors
			int xi= pbc(x0+v12nbr[i][0], nx);
			int yi= pbc(y0+v12nbr[i][1], ny);
			int zi= pbc(z0+v12nbr[i][2], nz);
            if((xi==x1) && (yi==y1) && (zi==z1)) continue; // connected to r1 or r2 are not non-broken
			if((xi==x2) && (yi==y2) && (zi==z2)) continue;
			
            int state_i= states[xi][yi][zi];
		    int nBnbr= 0;
                  
		    for(int a=0; a<n1nbr; a ++){ // 1st neighbors
			    int xj= pbc(xi+(*(v1nbr+a))[0], nx);
			    int yj= pbc(yi+(*(v1nbr+a))[1], ny);
			    int zj= pbc(zi+(*(v1nbr+a))[2], nz);
                if((xj==x1) && (yj==y1) && (zj==z1)) continue; // connected to r1 or r2 are not non-broken
			    if((xj==x2) && (yj==y2) && (zj==z2)) continue;
			    
                int state_j1= states[xj][yj][zj];
                
                if((-1==state_i && 1==state_j1) || (1==state_i && -1==state_j1)){
                    nBnbr= cal_Bnbr(nBnbr, xi, yi, zi);
                    AB1.push_back(nBnbr);
                }
                if((-1==state_i && 0==state_j1) || (0==state_i && -1==state_j1)){
                    nBnbr= cal_Bnbr(nBnbr, xi, yi, zi);
                    BV1.push_back(nBnbr);
                }
            }

		    for(int a=0; a<n2nbr; a ++){ // 2nd neighbors
		    	int xj= pbc(xi+(*(v2nbr+a))[0], nx);
			    int yj= pbc(yi+(*(v2nbr+a))[1], ny);
			    int zj= pbc(zi+(*(v2nbr+a))[2], nz);
                if((xj==x1) && (yj==y1) && (zj==z1)) continue; // connected to r1 or r2 are not non-broken
			    if((xj==x2) && (yj==y2) && (zj==z2)) continue;
			    
                int state_j2= states[xj][yj][zj];
            
                if((-1==state_i && 1==state_j2) || (1==state_i && -1==state_j2)){
                    nBnbr= cal_Bnbr(nBnbr, xi, yi, zi);
                    AB2.push_back(nBnbr);
                }
                if((-1==state_i && 0==state_j2) || (0==state_i && -1==state_j2)){
                    nBnbr= cal_Bnbr(nBnbr, xi, yi, zi);
                    BV2.push_back(nBnbr);
                }
            }
		}
	}

    const int    locals= n1nbr + n2nbr +1; // local sites
    double e= 0;
    for(int i=0; i < AB1.size(); i ++) e += e0A1B + e1A1B * (AB1[i]*1.0/locals);
    for(int i=0; i < AB2.size(); i ++) e += e0A2B + e1A2B * (AB2[i]*1.0/locals);
    for(int i=0; i < BV1.size(); i ++) e += e0B1V + e1B1V / (BV1[i]*1.0/locals); // !! different formula
    for(int i=0; i < BV2.size(); i ++) e += e0B2V + e1B2V * (BV2[i]*1.0/locals);
    
    return e/2; // a bond energy is avg of eAB of 2 atoms
}

double class_events::ecal_range(int xlo, int xhi, int ylo, int yhi, int zlo, int zhi, bool is_corr){
    if(is_corr){ // do correction on boundary
	    if(xlo>xhi){int temp= xlo; xlo= xhi; xhi= temp;}
	    if(ylo>yhi){int temp= ylo; ylo= yhi; yhi= temp;}
	    if(zlo>zhi){int temp= zlo; zlo= zhi; zhi= temp;}
        xlo -= 2; xhi += 2;
        ylo -= 2; yhi += 2;
        zlo -= 2; zhi += 2;
    }
    
    int bonds1[7][7]= {0}; // all shifted +2, e.g., AA: 4, BB: 0, AB: 5
	int bonds2[7][7]= {0}; // 2nd nn
    int AB1[25]= {0}; // AB for diff B concentrations for 1st-nn
    int AB2[25]= {0}; // AB for diff B concentrations for 2nd-nn
    int BV1[25]= {0}; // BV 1st-nn
    int BV2[25]= {0}; // BV 2nd-nn

	for(int ii= xlo; ii<= xhi; ii ++){ // includes
		for(int jj= ylo; jj<= yhi; jj ++){
			for(int kk= zlo; kk<= zhi; kk ++){
				int i= pbc(ii, nx);
				int j= pbc(jj, ny);
				int k= pbc(kk, nz);
                int state0= states[i][j][k];
                int nBnbr= 0;

                if(4==state0) continue;

				for(int a=0; a<n1nbr; a ++){ // 1st neighbors
					int x= pbc(i+(*(v1nbr+a))[0], nx);
					int y= pbc(j+(*(v1nbr+a))[1], ny);
					int z= pbc(k+(*(v1nbr+a))[2], nz);
					int state1a= states[x][y][z];

                    if(4==state1a) continue;

					bonds1[state0+2][state1a+2] ++;
                    
                    if((-1==state0 && 1==state1a) || (1==state0 && -1==state1a)){
                        nBnbr= cal_Bnbr(nBnbr, i, j, k);
                        AB1[nBnbr] ++;
                    }
                    if((-1==state0 && 0==state1a) || (0==state0 && -1==state1a)){
                        nBnbr= cal_Bnbr(nBnbr, i, j, k);
                        BV1[nBnbr] ++;
                    }
				}
	
				for(int b=0; b<n2nbr; b ++){ // 2nd neighbors
					int x= pbc(i+(*(v2nbr+b))[0], nx);
					int y= pbc(j+(*(v2nbr+b))[1], ny);
					int z= pbc(k+(*(v2nbr+b))[2], nz);
					int state1b= states[x][y][z];
                    
                    if(4==state1b) continue;

					bonds2[state0+2][state1b+2] ++;
                    
                    if((-1==state0 && 1==state1b) || (1==state0 && -1==state1b)){
                        nBnbr= cal_Bnbr(nBnbr, i, j, k);
                        AB2[nBnbr] ++;
                    }
                    if((-1==state0 && 0==state1b) || (0==state0 && -1==state1b)){
                        nBnbr= cal_Bnbr(nBnbr, i, j, k);
                        BV2[nBnbr] ++;
                    }
				}
	}}}
    
    // if ABV	
//	double e= eA1A * bonds1[3][3] + eA1V * bonds1[3][2] + eA1B * 0 +
//	          eA1V * bonds1[2][3] + eV1V * bonds1[2][2] + eV1B * bonds1[2][1] +
//	          eA1B * 0            + eV1B * bonds1[1][2] + eB1B * bonds1[1][1] +
//	          eA2A * bonds2[3][3] + eA2V * bonds2[3][2] + eA2B * 0 +
//	          eA2V * bonds2[2][3] + eV2V * bonds2[2][2] + eV2B * bonds2[2][1] +
//	          eA2B * 0            + eV2B * bonds2[1][2] + eB2B * bonds2[1][1];
	
    double e= 
    eM1M * bonds1[6][6]                                                 + eM1A  * bonds1[6][3] + eM1V * bonds1[6][2] + eM1B  * bonds1[6][1]                         + 
                          eAB1AB * bonds1[5][5] + eAA1AB * bonds1[5][4] + eA1AB * bonds1[5][3]                       + eAB1B * bonds1[5][1] + eAB1BB * bonds1[5][0] +
	                      eAA1AB * bonds1[4][5] + eAA1AA * bonds1[4][4] + eAA1A * bonds1[4][3]                       + eAA1B * bonds1[4][1] + eAA1BB * bonds1[4][0] +
	eM1A * bonds1[3][6] + eA1AB  * bonds1[3][5] + eAA1A  * bonds1[3][4] + eA1A  * bonds1[3][3] + eA1V * bonds1[3][2] + 0                    + eA1BB  * bonds1[3][0] +
	eM1V * bonds1[2][6]                                                 + eA1V  * bonds1[2][3] + eV1V * bonds1[2][2] + 0                                            +
	eM1B * bonds1[1][6] + eAB1B  * bonds1[1][5] + eAA1B  * bonds1[1][4] + 0                    + 0                   + eB1B  * bonds1[1][1] + eB1BB  * bonds1[1][0] +
	                      eAB1BB * bonds1[0][5] + eAA1BB * bonds1[0][4] + eA1BB * bonds1[0][3]                       + eB1BB * bonds1[0][1] + eBB1BB * bonds1[0][0] +
    // 2nd-nn
    eM2M * bonds2[6][6]                                                 + eM2A  * bonds2[6][3] + eM2V * bonds2[6][2] + eM2B  * bonds2[6][1]                         + 
	                      eAB2AB * bonds2[5][5] + eAA2AB * bonds2[5][4] + eA2AB * bonds2[5][3]                       + eAB2B * bonds2[5][1] + eAB2BB * bonds2[5][0] +
	                      eAA2AB * bonds2[4][5] + eAA2AA * bonds2[4][4] + eAA2A * bonds2[4][3]                       + eAA2B * bonds2[4][1] + eAA2BB * bonds2[4][0] +
	eM2A * bonds2[3][6] + eA2AB  * bonds2[3][5] + eAA2A  * bonds2[3][4] + eA2A  * bonds2[3][3] + eA2V * bonds2[3][2] + 0                    + eA2BB  * bonds2[3][0] +
	eM2V * bonds2[2][6]                                                 + eA2V  * bonds2[2][3] + eV2V * bonds2[2][2] + 0                                            +
	eM2B * bonds2[1][6] + eAB2B  * bonds2[1][5] + eAA2B  * bonds2[1][4] + 0                    + 0                   + eB2B  * bonds2[1][1] + eB2BB  * bonds2[1][0] +
	                      eAB2BB * bonds2[0][5] + eAA2BB * bonds2[0][4] + eA2BB * bonds2[0][3]                       + eB2BB * bonds2[0][1] + eBB2BB * bonds2[0][0];

    const int    locals= n1nbr + n2nbr +1; // local sites
    for(int x=1; x<=locals; x ++) e += AB1[x] * (e0A1B + e1A1B * (x*1.0/locals)); // 1st-nn
    for(int x=1; x<=locals; x ++) e += AB2[x] * (e0A2B + e1A2B * (x*1.0/locals)); // 2nd-nn
    for(int x=1; x<=locals; x ++) e += BV1[x] * (e0B1V + e1B1V / (x*1.0/locals));
    for(int x=1; x<=locals; x ++) e += BV2[x] * (e0B2V + e1B2V * (x*1.0/locals));
    
	return e/2;
}

double class_events::ecal_whole() const{ // DELETE: WRONG FUNC; REPLACED BY ECAL_RANGE 
	int bonds1[7][7]= {0}; // all shifted +2, e.g., AA: 4, BB: 0, AB: 5
	int bonds2[7][7]= {0}; // 2nd nn
	
	for(int i=0; i<nx; i ++){
		for(int j=0; j<ny; j ++){
			for(int k=0; k<nz; k ++){
				int state0= states[i][j][k];

				for(int a=0; a<n1nbr; a ++){ // 1st neighbors
					int x= pbc(i+(*(v1nbr+a))[0], nx);
					int y= pbc(j+(*(v1nbr+a))[1], ny);
					int z= pbc(k+(*(v1nbr+a))[2], nz);
					int state1a= states[x][y][z];

					bonds1[state0+2][state1a+2] ++;
				}
	
	if(! is_e2nbr) continue;

				for(int b=0; b<n2nbr; b ++){ // 2nd neighbors
					int x= pbc(i+(*(v2nbr+b))[0], nx);
					int y= pbc(j+(*(v2nbr+b))[1], ny);
					int z= pbc(k+(*(v2nbr+b))[2], nz);
					int state1b= states[x][y][z];

					bonds2[state0+2][state1b+2] ++;
				}
				
	}}}
	
	double e=             eM1M   * bonds1[6][6] + 0      * bonds1[6][5] + 0     * bonds1[6][4] + eM1A  * bonds1[6][3] + eM1V  * bonds1[6][2] + eM1B   * bonds1[6][1] +0 + 
                          eAB1AB * bonds1[5][5] + eAA1AB * bonds1[5][4] + eA1AB * bonds1[5][3] + 0     * bonds1[5][2] + eAB1B * bonds1[5][1] + eAB1BB * bonds1[5][0] +
	                      eAA1AB * bonds1[4][5] + eAA1AA * bonds1[4][4] + eAA1A * bonds1[4][3] + 0     * bonds1[4][2] + eAA1B * bonds1[4][1] + eAA1BB * bonds1[4][0] +
	eM1A * bonds1[3][6] + eA1AB  * bonds1[3][5] + eAA1A  * bonds1[3][4] + eA1A  * bonds1[3][3] + eA1V  * bonds1[3][2] + eA1B  * bonds1[3][1] + eA1BB  * bonds1[3][0] +
	eM1V * bonds1[2][6] + 0      * bonds1[2][5] + 0      * bonds1[2][4] + eA1V  * bonds1[2][3] + eV1V  * bonds1[2][2] + eV1B  * bonds1[2][1] + 0      * bonds1[2][0] +
	eM1B * bonds1[1][6] + eAB1B  * bonds1[1][5] + eAA1B  * bonds1[1][4] + eA1B  * bonds1[1][3] + eV1B  * bonds1[1][2] + eB1B  * bonds1[1][1] + eB1BB  * bonds1[1][0] +
	                      eAB1BB * bonds1[0][5] + eAA1BB * bonds1[0][4] + eA1BB * bonds1[0][3] + 0     * bonds1[0][2] + eB1BB * bonds1[0][1] + eBB1BB * bonds1[0][0];

	if(is_e2nbr)
	      e+=             eM2M   * bonds2[6][6] + 0      * bonds2[6][5] + 0     * bonds2[6][4] + eM2A  * bonds2[6][3] + eM2V  * bonds2[6][2] + eM2B   * bonds2[6][1] +0 + 
	                      eAB2AB * bonds2[5][5] + eAA2AB * bonds2[5][4] + eA2AB * bonds2[5][3] + 0     * bonds2[5][2] + eAB2B * bonds2[5][1] + eAB2BB * bonds2[5][0] +
	                      eAA2AB * bonds2[4][5] + eAA2AA * bonds2[4][4] + eAA2A * bonds2[4][3] + 0     * bonds2[4][2] + eAA2B * bonds2[4][1] + eAA2BB * bonds2[4][0] +
	eM2A * bonds2[3][6] + eA2AB  * bonds2[3][5] + eAA2A  * bonds2[3][4] + eA2A  * bonds2[3][3] + eA2V  * bonds2[3][2] + eA2B  * bonds2[3][1] + eA2BB  * bonds2[3][0] +
	eM2V * bonds2[2][6] + 0      * bonds2[2][5] + 0      * bonds2[2][4] + eA2V  * bonds2[2][3] + eV2V  * bonds2[2][2] + eV2B  * bonds2[2][1] + 0      * bonds2[2][0] +
	eM2B * bonds2[1][6] + eAB2B  * bonds2[1][5] + eAA2B  * bonds2[1][4] + eA2B  * bonds2[1][3] + eV2B  * bonds2[1][2] + eB2B  * bonds2[1][1] + eB2BB  * bonds2[1][0] +
	                      eAB2BB * bonds2[0][5] + eAA2BB * bonds2[0][4] + eA2BB * bonds2[0][3] + 0     * bonds2[0][2] + eB2BB * bonds2[0][1] + eBB2BB * bonds2[0][0];

	return e/2;
}
