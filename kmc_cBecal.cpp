#include "kmc_par.h"
#include "kmc_global.h"
using namespace std;

int cal_Bnbr(int N_Bnbr, int x, int y, int z){
    if(0==N_Bnbr){
        int n= 0;
        if(-1==states[x][y][z]) n ++;
        for(int a=0; a<n1nbr; a ++){ // search 1st-nn for B
		    int x1= pbc(x+(*(v1nbr+a))[0], nx);
		    int y1= pbc(y+(*(v1nbr+a))[1], ny);
		    int z1= pbc(z+(*(v1nbr+a))[2], nz);
		    if(-1==states[x1][y1][z1]) n ++;
	    }
	    for(int b=0; b<n2nbr; b ++){ // search 2nd-nn for B
		    int x2= pbc(x+(*(v2nbr+b))[0], nx);
		    int y2= pbc(y+(*(v2nbr+b))[1], ny);
		    int z2= pbc(z+(*(v2nbr+b))[2], nz);
		    if(-1==states[x2][y2][z2]) n ++;
        }

        return n;
    }
    else return N_Bnbr;
}

double ecal_one(int xi, int yi, int zi){ // use only for ABV metropolis MC 
	int state_i= states[xi][yi][zi];
	int bonds1[6][6]= {0}; // all shifted +2, e.g., AA: 4, BB: 0, AB: 5
	int bonds2[6][6]= {0}; // 2nd nn
    int N_Bnbr[nx][ny][nz]= {0}; // local B concentration
    int N_AB1[25]= {0}; // Nbonds for diff B concentrations for 1st-nn
    int N_AB2[25]= {0}; // Nbonds for diff B concentrations for 2nd-nn
    int N_BV1[25]= {0}; // BV 1st-nn
    int N_BV2[25]= {0}; // BV 2nd-nn

    for(int a=0; a<n1nbr; a ++){ // 1st neighbors
		int xj= pbc(xi+(*(v1nbr+a))[0], nx);
		int yj= pbc(yi+(*(v1nbr+a))[1], ny);
		int zj= pbc(zi+(*(v1nbr+a))[2], nz);
		int state_j1= states[xj][yj][zj];
		
        bonds1[state_i+2][state_j1+2] ++;
		
        if((1==state_i && -1==state_j1) || (-1==state_i && 1==state_j1)){ // AB bond for ij1
            N_Bnbr[xi][yi][zi]= cal_Bnbr(N_Bnbr[xi][yi][zi], xi, yi, zi);
            N_AB1[N_Bnbr[xi][yi][zi]] ++;
        }
        if((0==state_i && -1==state_j1) || (-1==state_i && 0==state_j1)){ 
            N_Bnbr[xi][yi][zi]= cal_Bnbr(N_Bnbr[xi][yi][zi], xi, yi, zi);
            N_BV1[N_Bnbr[xi][yi][zi]] ++;
        }
        for(int b=0; b<n1nbr; b ++){ // search 1st-nn for AB bonds
		    int xk= pbc(xj+(*(v1nbr+b))[0], nx);
		    int yk= pbc(yj+(*(v1nbr+b))[1], ny);
		    int zk= pbc(zj+(*(v1nbr+b))[2], nz);
		    int state_k1= states[xk][yk][zk];
        
            if((1==state_j1 && -1==state_k1) || (-1==state_j1 && 1==state_k1)){ // AB bond for j1k1
                N_Bnbr[xj][yj][zj]= cal_Bnbr(N_Bnbr[xj][yj][zj], xj, yj, zj);
                N_AB1[N_Bnbr[xj][yj][zj]] ++;
            }
            if((0==state_j1 && -1==state_k1) || (-1==state_j1 && 0==state_k1)){
                N_Bnbr[xj][yj][zj]= cal_Bnbr(N_Bnbr[xj][yj][zj], xj, yj, zj);
                N_BV1[N_Bnbr[xj][yj][zj]] ++;
            }
	    }
	    for(int b=0; b<n2nbr; b ++){ // search 2nd-nn for AB bonds
		    int xk= pbc(xj+(*(v2nbr+b))[0], nx);
		    int yk= pbc(yj+(*(v2nbr+b))[1], ny);
		    int zk= pbc(zj+(*(v2nbr+b))[2], nz);
		    int state_k2= states[xk][yk][zk];
            
            if((1==state_j1 && -1==state_k2) || (-1==state_j1 && 1==state_k2)){ // AB bond for j1k2
                N_Bnbr[xj][yj][zj]= cal_Bnbr(N_Bnbr[xj][yj][zj], xj, yj, zj);
                N_AB2[N_Bnbr[xj][yj][zj]] ++;
            }
            if((0==state_j1 && -1==state_k2) || (-1==state_j1 && 0==state_k2)){
                N_Bnbr[xj][yj][zj]= cal_Bnbr(N_Bnbr[xj][yj][zj], xj, yj, zj);
                N_BV2[N_Bnbr[xj][yj][zj]] ++;
            }
        }
    }
	
	for(int a=0; a<n2nbr; a ++){ // 2nd neighbors
		int xj= pbc(xi+(*(v2nbr+a))[0], nx);
		int yj= pbc(yi+(*(v2nbr+a))[1], ny);
		int zj= pbc(zi+(*(v2nbr+a))[2], nz);
		int state_j2= states[xj][yj][zj];
	
		bonds2[state_i+2][state_j2+2] ++;
        
        if((1==state_i && -1==state_j2) || (-1==state_i && 1==state_j2)){ // AB bond for ij2
            N_Bnbr[xi][yi][zi]= cal_Bnbr(N_Bnbr[xi][yi][zi], xi, yi, zi);
            N_AB2[N_Bnbr[xi][yi][zi]] ++;
        }
        if((0==state_i && -1==state_j2) || (-1==state_i && 0==state_j2)){
            N_Bnbr[xi][yi][zi]= cal_Bnbr(N_Bnbr[xi][yi][zi], xi, yi, zi);
            N_BV2[N_Bnbr[xi][yi][zi]] ++;
        }
        for(int b=0; b<n1nbr; b ++){ // search 1st-nn for AB bonds
		    int xk= pbc(xj+(*(v1nbr+b))[0], nx);
		    int yk= pbc(yj+(*(v1nbr+b))[1], ny);
		    int zk= pbc(zj+(*(v1nbr+b))[2], nz);
		    int state_k1= states[xk][yk][zk];
        
            if((1==state_j2 && -1==state_k1) || (-1==state_j2 && 1==state_k1)){ // AB bond for j2k1
                N_Bnbr[xj][yj][zj]= cal_Bnbr(N_Bnbr[xj][yj][zj], xj, yj, zj);
                N_AB1[N_Bnbr[xj][yj][zj]] ++;
            }
            if((0==state_j2 && -1==state_k1) || (-1==state_j2 && 0==state_k1)){
                N_Bnbr[xj][yj][zj]= cal_Bnbr(N_Bnbr[xj][yj][zj], xj, yj, zj);
                N_BV1[N_Bnbr[xj][yj][zj]] ++;
            }
	    }
	    for(int b=0; b<n2nbr; b ++){ // search 2nd-nn for AB bonds
		    int xk= pbc(xj+(*(v2nbr+b))[0], nx);
		    int yk= pbc(yj+(*(v2nbr+b))[1], ny);
		    int zk= pbc(zj+(*(v2nbr+b))[2], nz);
		    int state_k2= states[xk][yk][zk];
            
            if((1==state_j2 && -1==state_k2) || (-1==state_j2 && 1==state_k2)){ // AB bond for j2k2
                N_Bnbr[xj][yj][zj]= cal_Bnbr(N_Bnbr[xj][yj][zj], xj, yj, zj);
                N_AB2[N_Bnbr[xj][yj][zj]] ++;
            }
            if((0==state_j2 && -1==state_k2) || (-1==state_j2 && 0==state_k2)){
                N_Bnbr[xj][yj][zj]= cal_Bnbr(N_Bnbr[xj][yj][zj], xj, yj, zj);
                N_BV2[N_Bnbr[xj][yj][zj]] ++;
            }
        }
	}

    double e= 0; // eAB(X)(1st)+eAB(X)(2nd)= w0+w1*X +0.5(eAA(1)+eBB(1)) +0.5(eAA(2)+eBB(2))
   
    const int    locals= n1nbr + n2nbr +1; // local sites
    for(int x=1; x<=locals; x ++) e += N_AB1[x] * (e0A1B + e1A1B*(x*1.0/locals)); // 1st-nn
    for(int x=1; x<=locals; x ++) e += N_AB2[x] * (e0A2B + e1A2B*(x*1.0/locals)); // 2nd-nn
    for(int x=1; x<=locals; x ++) e += N_BV1[x] * (e0B1V + e1B1V/(x*1.0/locals)); // 1st-nn(!! different formula)
    for(int x=1; x<=locals; x ++) e += N_BV2[x] * (e0B2V + e1B2V*(x*1.0/locals)); // 2nd-nn
	e /= 2.0; // AB bond contribution
	
    e+= eA1A * bonds1[3][3] + eA1V * bonds1[3][2] + 
	    eA1V * bonds1[2][3] + eV1V * bonds1[2][2] + 0                    +
                              0                   + eB1B  * bonds1[1][1] +
        eA2A * bonds2[3][3] + eA2V * bonds2[3][2] +
	    eA2V * bonds2[2][3] + eV2V * bonds2[2][2] + 0                    +
	                          0                   + eB2B  * bonds2[1][1];

	return e;
}

double ecal_range(int xlo, int xhi, int ylo, int yhi, int zlo, int zhi){
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

				for(int a=0; a<n1nbr; a ++){ // 1st neighbors
					int x= pbc(i+(*(v1nbr+a))[0], nx);
					int y= pbc(j+(*(v1nbr+a))[1], ny);
					int z= pbc(k+(*(v1nbr+a))[2], nz);
		
					int state1a= states[x][y][z];
					if(itlAB[x][y][z]) state1a= 3;

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
					if(itlAB[x][y][z]) state1b= 3;

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
	
	double e= eA1A * bonds1[3][3] + eA1V * bonds1[3][2] + 0                   +
	          eA1V * bonds1[2][3] + eV1V * bonds1[2][2] + 0                   +
	          0                   + 0                   + eB1B * bonds1[1][1] +
	          eA2A * bonds2[3][3] + eA2V * bonds2[3][2] + 0                   +
	          eA2V * bonds2[2][3] + eV2V * bonds2[2][2] + 0                   +
	          0                   + 0                   + eB2B * bonds2[1][1];
    
    const int    locals= n1nbr + n2nbr +1; // local sites
    for(int x=1; x<=locals; x ++) e += AB1[x] * (e0A1B + e1A1B*(x*1.0/locals)); // 1st-nn
    for(int x=1; x<=locals; x ++) e += AB2[x] * (e0A2B + e1A2B*(x*1.0/locals)); // 2nd-nn
    for(int x=1; x<=locals; x ++) e += BV1[x] * (e0B1V + e1B1V/(x*1.0/locals)); // 1st-nn(!! different formula)
    for(int x=1; x<=locals; x ++) e += BV2[x] * (e0B2V + e1B2V*(x*1.0/locals)); // 2nd-nn

	return e/2;
}

double ecal_swap(int x1, int y1, int z1, int x2, int y2, int z2){
#define NDIS 2 // define neighbor distance
    int dx= x2 - x1;
    int dy= y2 - y1;
    int dz= z2 - z1;
    int xlo, xhi, ylo, yhi, zlo, zhi;
    bool ISnear= true;

    int type[3]; // 0:OUT; 1:IN(pbc)&2>1; 2:IN&2>1; 3:IN&1>2; 4:IN(pbc)&1>2

    // x
    if(dx >= (nx-NDIS))             { xlo= x2-2; xhi= x1+nx+2;}
    else if(dx <= NDIS && dx >= 0)  { xlo= x1-2; xhi= x2+2;}
    else if(dx >= -NDIS && dx < 0)  { xlo= x2-2; xhi= x1+2;}
    else if(dx <= -(nx-NDIS))       { xlo= x1-2; xhi= x2+nx+2;}
    else                            ISnear= false;
    // y
    if(dy >= (ny-NDIS))             { ylo= y2-2; yhi= y1+ny+2;}
    else if(dy <= NDIS && dy >= 0)  { ylo= y1-2; yhi= y2+2;}
    else if(dy >= -NDIS && dy < 0)  { ylo= y2-2; yhi= y1+2;}
    else if(dy <= -(ny-NDIS))       { ylo= y1-2; yhi= y2+ny+2;}
    else                            ISnear= false;
    // z
    if(dz >= (nz-NDIS))             { zlo= z2-2; zhi= z1+nz+2;}
    else if(dz <= NDIS && dz >= 0)  { zlo= z1-2; zhi= z2+2;}
    else if(dz >= -NDIS && dz < 0)  { zlo= z2-2; zhi= z1+2;}
    else if(dz <= -(nz-NDIS))       { zlo= z1-2; zhi= z2+nz+2;}
    else                            ISnear= false;

    if(ISnear)  return ecal_range(xlo, xhi, ylo, yhi, zlo, zhi); 
    else        return ecal_one(x1, y1, z1)+ecal_one(x2, y2, z2); 
}
