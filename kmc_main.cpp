#include <iostream>
#include <fstream>
#include <ctime>
#include "kmc_par.h"
#include "kmc_global.h"
#include "kmc_initial.h"
#include "kmc_events.h"
using namespace std;

int main(int nArg, char *Arg[]){
	int t0cpu= time(0);

	long long int ts_bg;
	double time_bg;

	cout << "########## Initializing System... ##########" << endl;
	class_initial init(ts_bg, time_bg, nArg, Arg[1]);
	
	cout << "\n########## Initializing Events ... ##########" << endl;
	class_events events; 

	cout << "\n########## The Simulation Begins !! ##########" << endl;
	timestep=  ts_bg; totaltime= time_bg;
	cout << "TIMESTEP() TIME(s) GENR()    NA() nV() NAA()    NVD() NVS()    nIrecb nIsink";
	printf("\n%lld %.10e %d    %d %d %d    %d %d    %d %d", timestep, totaltime, N_genr, nA, nV, nAA, nVD, nVs, nIrecb, nIsink);
    int N_0def= 0;
    int N_conf= 0;
	while((totaltime<= time_bg+par_time) && (timestep != ts_bg+par_step)){
        timestep ++;
		
		// CALCULATIONS
		double dt= events.main(); // Defect jumps
        totaltime += dt;
        if(abs(dt-1.0/rate_genr)<1e-10) N_0def ++;

		// OUTPUT DATA
		if(0==timestep%step_log){
            cout << endl;
	        printf("%lld %.10e %d    %d %d %d    %d %d    %d %d", timestep, totaltime, N_genr, nA, nV, nAA, nVD, nVs, nIrecb, nIsink);
            if(N_0def != 0){            cout << "  *** 0-defect genr: " << N_0def << " ***";            N_0def= 0;}
            if(events.N_nediff != 0){   cout << "  *** neg e: " << events.N_nediff << " ***";           events.N_nediff= 0;}
		}

		if(0==timestep%step_conf || totaltime>=(N_conf+1)*time_conf){
			if(0==timestep%step_conf)   write_conf(2);
            else                        write_conf(1);
            cout << "   <Output conf files at: " << timestep << ">";
            N_conf= (int) (totaltime / time_conf);
        }

		if(0==timestep%step_out){
//			fprintf(out_engy, "%lld %e %f\n", timestep, totaltime, events.ecal_range());
//            fflush(out_engy);
			fprintf(out_sro, "%lld %e %f\n", timestep, totaltime, cal_sro());
	        fprintf(out_log, "%lld %.10e %d    %d %d %d    %d %d    %d %d\n", timestep, totaltime, N_genr, nA, nV, nAA, nVD, nVs, nIrecb, nIsink);
            fflush(out_sro); fflush(out_log);
			nIrecb= 0; nIsink= 0;
            if(par_isOUTrestart) write_conf(3);
        }
		
		if(0==timestep%step_his){
 			write_hissol(); // write sol first to construct list_srf
			write_hisdef(); // and then output his_srf here
		}
	}

	// finalizing
	if(timestep%step_log != 0) printf("\n%lld %f %d	%d %d	%d %d %d   %d %d %d ", timestep, totaltime, N_genr, nA, nB, nV, nVD, nVs, nAA, nAB, nBB);
	write_conf(false); cout << "<Output conf files at: " << timestep << ">";
    
	int tfcpu= time(0);
	cout << "\n**** The simulation is done! Total CPU time: " << tfcpu - t0cpu << " secs ****" << endl;

	return 0;
}
