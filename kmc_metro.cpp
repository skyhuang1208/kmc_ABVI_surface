#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include "kmc_par.h"
#include "kmc_global.h"
#include "kmc_initial.h"
using namespace std;

double ecal_metro(int x1, int y1, int z1, int x2, int y2, int z2); // use only for ABV metropolis MC 
double ecal_metro_whole();

int main(int nArg, char *Arg[]){
	int t0cpu= time(0);

	long long int ts_bg;
	double time_bg;

    cout << "********** Metropolis MC simulation **********\n" << endl;
	cout << "########## Initializing System... ##########" << endl;
	class_initial init(ts_bg, time_bg, nArg, Arg[1]);

	cout << "\n########## The Simulation Begins !! ##########" << endl;
	timestep=  ts_bg;
	cout << "TIMESTEP() TIME(s) GENR()	NA() NB()	NV() NAA() NAB() NBB()	AJUMP_V% AJUMP_I%";
	printf("\n%d %e %d     %d %d     %d %d %d %d", 0, 0.0, 0, nA, nB, nV, nAA, nAB, nBB);
    long long int N_attempt= 0;
    int N_stuck= 0;
    while(timestep != ts_bg+par_step){
        int s1= ran_generator()*nx*ny*nz; // site 1, 2
        int s2= ran_generator()*nx*ny*nz;
        int t1= *(&states[0][0][0]+s1);   // type 1, 2
        int t2= *(&states[0][0][0]+s2);
        if(t1==t2) continue;
        
        N_attempt ++;

		int x1= (s1/nz)/ny, y1= (s1/nz)%ny, z1= s1%nz;
		int x2= (s2/nz)/ny, y2= (s2/nz)%ny, z2= s2%nz;
	    double e0= ecal_metro(x1, y1, z1, x2, y2, z2);
	    
        states[x1][y1][z1]= t2;
	    states[x2][y2][z2]= t1;
	    double ediff= ecal_metro(x1, y1, z1, x2, y2, z2) - e0;

        if(ediff>1e-10){
            double prob= exp(-beta*ediff);
            if(ran_generator()>=prob){
                states[x1][y1][z1]= t1;
                states[x2][y2][z2]= t2;
                continue;
            }
        }
        
        timestep ++;
#define DEF_STUCK 100000
        if(N_attempt>DEF_STUCK) N_stuck ++;
        N_attempt= 0;
    
        if(nA+nB+nV+nAA+nBB+nAB+nM != nx*ny*nz) error(2, "(jump) numbers of ltc points arent consistent, diff=", 1, nA+nB+nV+nAA+nBB+nAB+nM-nx*ny*nz); // check
        if(2*nAA+nA-nB-2*nBB       != sum_mag)  error(2, "(jump) magnitization isnt conserved", 2, 2*nAA+nA-nB-2*nBB, sum_mag);

		// OUTPUT DATA
		if(0==timestep%step_log) printf("\n%lld %.10e %d     %d %d     %d %d %d %d", timestep, totaltime, N_genr, nA, nB, nV, nAA, nAB, nBB); 
		if(0==timestep%step_his) write_metrohis();
		if(0==timestep%step_out) fprintf(out_engy, "%lld %e %f\n", timestep, totaltime, ecal_metro_whole());
		if(0==timestep%step_conf){ write_conf(); cout << "   <Output conf files at: " << timestep << ">";}

#define TOL_STUCK 10000
        if(N_stuck>TOL_STUCK){
            cout << "Reach the criterion of breaking Metropolis MC" << endl;
            cout << "Criterion: " << DEF_STUCK << " # of attempts for " << TOL_STUCK << " times." << endl;
            write_metrohis();
            break;
        }
	}

	// finalizing
	if(timestep%step_log != 0) printf("\n%lld %f %d	%d %d	%d %d %d %d ", timestep, totaltime, N_genr, nA, nB, nV, nAA, nAB, nBB);
	if(timestep%step_out != 0) fprintf(out_engy, "%lld %e %f\n", timestep, totaltime, ecal_metro_whole());
	write_conf(); cout << "<Output conf files at: " << timestep << ">";
    
	int tfcpu= time(0);
	cout << "\n**** The simulation is done! Total CPU time: " << tfcpu - t0cpu << " secs ****" << endl;

	return 0;
}

double ecal_metro(int x1, int y1, int z1, int x2, int y2, int z2){ // use only for ABV metropolis MC 
	int pos[2][3]= {{x1, y1, z1}, {x2, y2, z2}};

	int bonds1[6][6]= {0}; // all shifted +2, e.g., AA: 4, BB: 0, AB: 5
	int bonds2[6][6]= {0}; // 2nd nn
	for(int i=0; i<2; i ++){
		int xi= pos[i][0];
		int yi= pos[i][1];
		int zi= pos[i][2];
		int state_i= states[xi][yi][zi];

		for(int a=0; a<n1nbr; a ++){ // 1st neighbors
			int xj= pbc(xi+(*(v1nbr+a))[0], nx);
			int yj= pbc(yi+(*(v1nbr+a))[1], ny);
			int zj= pbc(zi+(*(v1nbr+a))[2], nz);
			int state_j1= states[xj][yj][zj];
		
			if((xj==x1) && (yj==y1) && (zj==z1)) continue; // avoid double count the (x1, y1, z1)-(x2, y2, z2) bond

			bonds1[state_i+2][state_j1+2] ++;
		}
	
		if(! is_e2nbr) continue;

		for(int a=0; a<n2nbr; a ++){ // 2nd neighbors
			int xj= pbc(xi+(*(v2nbr+a))[0], nx);
			int yj= pbc(yi+(*(v2nbr+a))[1], ny);
			int zj= pbc(zi+(*(v2nbr+a))[2], nz);
			int state_j2= states[xj][yj][zj];
		
			if((xj==x1) && (yj==y1) && (zj==z1)) continue; // should not happen in 2nd neighbors

			bonds2[state_i+2][state_j2+2] ++;
		}
				
	}
	
	double e= eAB1AB * bonds1[5][5] + eAA1AB * bonds1[5][4] + eA1AB * bonds1[5][3] + 0     * bonds1[5][2] + eAB1B * bonds1[5][1] + eAB1BB * bonds1[5][0] +
	          eAA1AB * bonds1[4][5] + eAA1AA * bonds1[4][4] + eAA1A * bonds1[4][3] + 0     * bonds1[4][2] + eAA1B * bonds1[4][1] + eAA1BB * bonds1[4][0] +
	          eA1AB  * bonds1[3][5] + eAA1A  * bonds1[3][4] + eA1A  * bonds1[3][3] + eA1V  * bonds1[3][2] + eA1B  * bonds1[3][1] + eA1BB  * bonds1[3][0] +
	          0      * bonds1[2][5] + 0      * bonds1[2][4] + eA1V  * bonds1[2][3] + eV1V  * bonds1[2][2] + eV1B  * bonds1[2][1] + 0      * bonds1[2][0] +
	          eAB1B  * bonds1[1][5] + eAA1B  * bonds1[1][4] + eA1B  * bonds1[1][3] + eV1B  * bonds1[1][2] + eB1B  * bonds1[1][1] + eB1BB  * bonds1[1][0] +
	          eAB1BB * bonds1[0][5] + eAA1BB * bonds1[0][4] + eA1BB * bonds1[0][3] + 0     * bonds1[0][2] + eB1BB * bonds1[0][1] + eBB1BB * bonds1[0][0];

	if(is_e2nbr)
	      e+= eAB2AB * bonds2[5][5] + eAA2AB * bonds2[5][4] + eA2AB * bonds2[5][3] + 0     * bonds2[5][2] + eAB2B * bonds2[5][1] + eAB2BB * bonds2[5][0] +
	          eAA2AB * bonds2[4][5] + eAA2AA * bonds2[4][4] + eAA2A * bonds2[4][3] + 0     * bonds2[4][2] + eAA2B * bonds2[4][1] + eAA2BB * bonds2[4][0] +
	          eA2AB  * bonds2[3][5] + eAA2A  * bonds2[3][4] + eA2A  * bonds2[3][3] + eA2V  * bonds2[3][2] + eA2B  * bonds2[3][1] + eA2BB  * bonds2[3][0] +
	          0      * bonds2[2][5] + 0      * bonds2[2][4] + eA2V  * bonds2[2][3] + eV2V  * bonds2[2][2] + eV2B  * bonds2[2][1] + 0      * bonds2[2][0] +
	          eAB2B  * bonds2[1][5] + eAA2B  * bonds2[1][4] + eA2B  * bonds2[1][3] + eV2B  * bonds2[1][2] + eB2B  * bonds2[1][1] + eB2BB  * bonds2[1][0] +
	          eAB2BB * bonds2[0][5] + eAA2BB * bonds2[0][4] + eA2BB * bonds2[0][3] + 0     * bonds2[0][2] + eB2BB * bonds2[0][1] + eBB2BB * bonds2[0][0];

	return e;
}

double ecal_metro_whole(){
	int bonds1[7][7]= {0}; // all shifted +2, e.g., AA: 4, BB: 0, AB: 5
	int bonds2[7][7]= {0}; // 2nd nn
	
	for(int i=0; i<nx; i ++){
		for(int j=0; j<ny; j ++){
			for(int k=0; k<nz; k ++){
		
				int state0= states[i][j][k];
				if(itlAB[i][j][k]) state0= 3;

				for(int a=0; a<n1nbr; a ++){ // 1st neighbors
					int x= pbc(i+(*(v1nbr+a))[0], nx);
					int y= pbc(j+(*(v1nbr+a))[1], ny);
					int z= pbc(k+(*(v1nbr+a))[2], nz);
		
					int state1a= states[x][y][z];
					if(itlAB[x][y][z]) state1a= 3;

					bonds1[state0+2][state1a+2] ++;
				}
	
	if(! is_e2nbr) continue;

				for(int b=0; b<n2nbr; b ++){ // 2nd neighbors
					int x= pbc(i+(*(v2nbr+b))[0], nx);
					int y= pbc(j+(*(v2nbr+b))[1], ny);
					int z= pbc(k+(*(v2nbr+b))[2], nz);
		
					int state1b= states[x][y][z];
					if(itlAB[x][y][z]) state1b= 3;

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

	return e;
}
