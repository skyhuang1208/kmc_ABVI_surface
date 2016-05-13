#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include "kmc_par.h"
#include "kmc_global.h"
#include "kmc_initial.h"
using namespace std;

double cal_sro();

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
    long long int N_attempt= 0, sum_att= 0;
    int N_stuck= 0;
    while(timestep != par_step){ // !! Warning: when restart it runs until par_step not plus !!
        N_attempt ++;
        
        int s= ran_generator()*nx*ny*nz; // site
        int t= *(&states[0][0][0]+s);    // type
		int x= (s/nz)/ny, y= (s/nz)%ny, z= s%nz;
        int s2, t2, x2, y2, z2;
        double ediff;
	    
        if(0==t){ // VCC SWAP
            do{
                s2= ran_generator()*nx*ny*nz; // site
                t2= *(&states[0][0][0]+s2);   // type
            } while(0==t2);
		    x2= (s2/nz)/ny; y2= (s2/nz)%ny; z2= s2%nz;
        
            double e0= ecal_swap(x, y, z, x2, y2, z2);
            states[x][y][z]= t2;
            states[x2][y2][z2]= t;
            ediff= ecal_swap(x, y, z, x2, y2, z2) - e0;
        }
        else{ // A-B FLIP
            double e0= ecal_one(x, y, z);
            states[x][y][z]= (1==t) ? -1:1;
            int dmag= (1==t) ? -2:2;
            ediff= ecal_one(x, y, z) - e0 + mu*dmag; // effective ediff for semi-grand canonical
        }
            
        if(ediff>1e-10){
            double prob= exp(-beta*ediff);
            if(ran_generator()>=prob){
                states[x][y][z]= t;
                if(0==t) states[x2][y2][z2]= t2;
                continue;
            }
        }

        if(t != 0){
            switch(t){
                case -1: nA ++; nB --; break;
                case  1: nA --; nB ++; break;
                default: error(2, "(flip) an unknown type", 1, t);
            }
        }
        
        timestep ++;

#define DEF_STUCK 100000
        if(N_attempt>DEF_STUCK) N_stuck ++;
        sum_att+= N_attempt; N_attempt= 0;
   
        if(nA+nB+nV+nAA+nBB+nAB+nM != nx*ny*nz) error(2, "(jump) numbers of ltc points arent consistent, diff=", 1, nA+nB+nV+nAA+nBB+nAB+nM-nx*ny*nz); // check
        mag= nA-nB; // magnitization

		// OUTPUT DATA
		if(0==timestep%step_log){
            cout << endl;
            printf("%lld %.10e %d     %d %d     %d %d %d %d  %f", timestep, totaltime, N_genr, nA, nB, nV, nAA, nAB, nBB, sum_att*1.0/step_log);
            sum_att= 0;
        }
		if(0==timestep%step_his) write_metrohis();
		if(0==timestep%step_out){ 
            double etotal= ecal_range();
            fprintf(out_engy, "%lld %e %f\n", timestep, totaltime, etotal); fflush(out_engy);
            fprintf(out_eSGC, "%lld %e %f\n", timestep, totaltime, etotal+mu*(nA-nB)); fflush(out_eSGC);
            fprintf(out_sro, "%lld %e %f\n", timestep, totaltime, cal_sro()); fflush(out_sro);
        }
		if(0==timestep%step_conf){ write_conf(); cout << "   <Output conf files at: " << timestep << ">";}

#define TOL_STUCK 1000
        if(N_stuck>TOL_STUCK){
            cout << "Reach the criterion of breaking Metropolis MC" << endl;
            cout << "Criterion: " << DEF_STUCK << " # of attempts for " << TOL_STUCK << " times." << endl;
            write_metrohis();
            break;
        }
	}

	// finalizing
	if(timestep%step_log != 0) printf("\n%lld %f %d	%d %d	%d %d %d %d ", timestep, totaltime, N_genr, nA, nB, nV, nAA, nAB, nBB);
	if(timestep%step_out != 0) fprintf(out_engy, "%lld %e %f\n", timestep, totaltime, ecal_range());
	if(N_stuck != 0) cout << "<N of stuck: " << N_stuck << ">";
	write_conf(); cout << "<Output conf files at: " << timestep << ">";
    
	int tfcpu= time(0);
	cout << "\n**** The simulation is done! Total CPU time: " << tfcpu - t0cpu << " secs ****" << endl;

	return 0;
}

double cal_sro(){
    double cA= nA*1.0/(nx*ny*nz);
    double sro= 0;

    int ncheck= 0;
	for(int i=0; i<nx; i ++){
		for(int j=0; j<ny; j ++){
			for(int k=0; k<nz; k ++){
				int state0= states[i][j][k];

                if(-1==state0){
                    ncheck ++;
                    int nAnbr= 0;

				    for(int a=0; a<n1nbr; a ++){ // 1st neighbors
					    int x= pbc(i+(*(v1nbr+a))[0], nx);
					    int y= pbc(j+(*(v1nbr+a))[1], ny);
					    int z= pbc(k+(*(v1nbr+a))[2], nz);
                        int state1= states[x][y][z];

                        if(1==state1) nAnbr ++;
                    }
				    
                    for(int a=0; a<n2nbr; a ++){ // 1st neighbors
					    int x= pbc(i+(*(v2nbr+a))[0], nx);
					    int y= pbc(j+(*(v2nbr+a))[1], ny);
					    int z= pbc(k+(*(v2nbr+a))[2], nz);
                        int state1= states[x][y][z];

                        if(1==state1) nAnbr ++;
                    }

                    sro += 1.0 - nAnbr*1.0/(n1nbr+n2nbr)/cA;
                }
    }}}
    
    if(ncheck != nB) error(2, "(cal_sro) number inconsistent", 2, ncheck, nB);
    return sro/nB;
}
