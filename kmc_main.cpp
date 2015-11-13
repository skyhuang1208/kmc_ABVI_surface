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
	cout << "TIMESTEP() TIME(s) GENR()	NA() NB()	NV() NAA() NAB() NBB()	AJUMP_V% AJUMP_I%";
	printf("\n%d %e %d     %d %d     %d %d %d %d     %f %f", 0, 0.0, 0, nA, nB, nV, nAA, nAB, nBB, 0.0, 0.0);
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
			printf("\n%lld %.10e %d     %d %d     %d %d %d %d     %f %f    %f", timestep, totaltime, N_genr, nA, nB, nV, nAA, nAB, nBB, 
                                                        ((double) Vja[0])/(Vja[0]+Vja[1]), ((double) Ija[0])/(Ija[0]+Ija[1]), events.cvcc_rates);
			if(N_0def != 0){
				cout << "  *** 0-defect genr: " << N_0def << " ***"; 
                N_0def= 0;
			}
		}

		if(0==timestep%step_conf || totaltime>=(N_conf+1)*time_conf){
			write_conf(); 
            cout << "   <Output conf files at: " << timestep << ">";
            N_conf= (int) (totaltime / time_conf);
        }

		if(0==timestep%step_out){
			fprintf(out_engy, "%lld %e %f\n", timestep, totaltime, events.ecal_whole());
            write_vdep();
        }
		
		if(0==timestep%step_his){
 			write_hissol(); // write sol first to construct list_srf
			write_hisdef(); // and then output his_srf here
		}
	}

	// finalizing
	if(timestep%step_log != 0) printf("\n%lld %f %d	%d %d	%d %d %d %d ", timestep, totaltime, N_genr, nA, nB, nV, nAA, nAB, nBB);
	write_conf(); cout << "<Output conf files at: " << timestep << ">";
	if(timestep%step_out != 0){
        fprintf(out_engy, "%lld %e %f\n", timestep, totaltime, events.ecal_whole());
        write_vdep();
    }
    
    for(int i=8; i>=0; i --) njump[i] += njump[i+1]; // output ratio of effective vcc creation
    cout << "\n## vcc creation from srf ##" << endl;
    cout << "total: " << njump[0] << endl;
    cout << ", >=1: " << njump[1] << ", >=2: " << njump[2] << ", >=3: " << njump[3] << ", >=4: " << njump[4] << ", >=5: " << njump[5] << endl;
    cout << ", >=6: " << njump[6] << ", >=7: " << njump[7] << ", >=8: " << njump[8] << ", >=9: " << njump[9] << endl;
    cout << "recbED at 0: " << njump[10] << endl;

	int tfcpu= time(0);
	cout << "\n**** The simulation is done! Total CPU time: " << tfcpu - t0cpu << " secs ****" << endl;

	return 0;
}
