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

    double init_sro= cal_sro();

    if(0==timestep){ fprintf(out_engy, "%lld %e %f\n", timestep, totaltime, events.ecal_range()); 
                     fprintf(out_sro,  "%lld %e %e %e %e %e %e %e\n", timestep, totaltime, init_sro, 0.0, 0.0, 0.0, 0.0, 0.0); }
    cout << "TIMESTEP() TIME(s) GENR()	NA() NB()	NV() NAA() NAB() NBB()	AJUMP_V% AJUMP_I%";
	printf("\n%lld %.10e %d     %d %d     %d %d %d %d     %f %f", timestep, totaltime, 0, nA, nB, nV, nAA, nAB, nBB, 0.0, 0.0);
    int N_0def= 0; // genr performed when no other event
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
			printf("%lld %.10e %d     %d %d     %d %d %d %d     %f %f", timestep, totaltime, N_genr, nA, nB, nV, nAA, nAB, nBB, 
                    ((double) Vja[0])/(Vja[0]+Vja[1]), ((double) Ija[0])/(Ija[0]+Ija[1]));
			if(N_0def != 0){
				cout << "  *** 0-defect genr: " << N_0def << " ***"; 
                N_0def= 0;
			}
            if(events.N_nediff != 0){
				cout << "  *** neg e: " << events.N_nediff << " ***"; 
                events.N_nediff= 0;
            }
            fflush(stdout);
		}

		if(0==timestep%step_conf || totaltime>=(N_conf+1)*time_conf){
			if(0==timestep%step_conf) write_conf(1);
            else                      write_conf(2);
            cout << "   <Output conf files at: " << timestep << ">";
            N_conf= (int) (totaltime / time_conf);
        }

		if(0==timestep%step_out){
			fprintf(out_engy, "%lld %e %f\n", timestep, totaltime, events.ecal_range()); fflush(out_engy);
            double sro= cal_sro();
            double acc_dsroI= sro - init_sro - acc_dsroV - acc_dsroRi - acc_dsroRv - acc_dsroG;
			fprintf(out_sro,  "%lld %e %e %e %e %e %e %e\n", timestep, totaltime, sro, acc_dsroI, acc_dsroV, acc_dsroRi, acc_dsroRv, acc_dsroG); fflush(out_sro);
            if(1==(nV+nAA+nAB+nBB) && (! is_genr)){ fprintf(out_msd,  "%lld %e %f\n", timestep, totaltime, cal_msd()); fflush(out_msd); }
//            write_vdep();
            if(par_isOUTrestart) write_conf(3);
        }
		
		if(0==timestep%step_his){
 			write_hissol(); // write sol first to construct list_srf
			write_hisdef(); // and then output his_srf here
		}
	}

	// finalizing
	if(timestep%step_log != 0) printf("\n%lld %f %d	%d %d	%d %d %d %d ", timestep, totaltime, N_genr, nA, nB, nV, nAA, nAB, nBB);
	write_conf(1); cout << "<Output conf files at: " << timestep << ">";
    cout << "\nAccumulative sro change from vcc: " << acc_dsroV << endl;

	int tfcpu= time(0);
	cout << "**** The simulation is done! Total CPU time: " << tfcpu - t0cpu << " secs ****" << endl;

	return 0;
}
