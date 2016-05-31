#ifndef KMC_SYSTEM_INCLUDED
#define KMC_SYSTEM_INCLUDED
#include <iostream>
#include <cstring>
#include "kmc_global.h"
#include "kmc_par.h"

#define MAX_NSPE 20

using namespace std;

class class_initial{
	public:
		class_initial(long long int &ts_bg, double &time_bg, int nArg, char name_restart[]){
			stpcpy(type_ltc, par_ltc);

			ltc_constructor();	// execute constructor for vbra, n*nbr, and v*nbr
			
			// Printf the parameters
			cout << "Setting: " << endl << "nx= " << nx << ", ny= " << ny << ", nz= " << nz << endl;
			cout << "The crystal structure is: " << type_ltc << endl;
			for(int i=0; i<3; i++)
				cout << "v" << i << ": " << vbra[i][0] << " " << vbra[i][1] << " " << vbra[i][2] << endl;
			cout << "And the number of 1st neighbors: " << n1nbr << endl;
			cout << "    the number of 2nd neighbors: " << n2nbr << endl;
		
			cout << "\n####### Generating configuration... #######" << endl;
            if(par_isrestart){
				cout << "RESTART FROM restart file..." << endl;
				if(nArg != 2) error(0, "when restart flag is true, nArg must be 2");
				
				read_restart(name_restart, ts_bg, time_bg);
				
			    his_sol = fopen(par_name_sol,  "a");
			    his_def = fopen(par_name_def,  "a");
			    his_srf = fopen(par_name_srf,  "a");
			    out_engy= fopen(par_name_engy, "a");
			    out_vdep= fopen(par_name_vdep, "a");
			    out_sro = fopen(par_name_sro,  "a");
				cout << "Open " << par_name_sol << " & " << par_name_def << " with append mode" << endl;
            }
			else{
				cout << "START FROM a random configuration..." << endl;
				ts_bg= 0; time_bg= 0;
				
				init_states_array(par_compV, par_compA, par_nMlayer);
				write_conf();
				cout << "Output t0 conf files" << endl;
				
			    his_sol = fopen(par_name_sol,  "w");
			    his_def = fopen(par_name_def,  "w");
			    his_srf = fopen(par_name_srf,  "w");
			    out_engy= fopen(par_name_engy, "w");
			    out_vdep= fopen(par_name_vdep, "w");
			    out_sro = fopen(par_name_sro,  "w");
                cout << "Open " << par_name_sol << " & " << par_name_def << " with write mode" << endl;
			}
			if(NULL==his_sol)  error(2, "(class_events) the solute  history file was not opened!");
			if(NULL==his_def)  error(2, "(class_events) the vacancy history file was not opened!");
			if(NULL==his_srf)  error(2, "(class_events) the surface history file was not opened!");
			if(NULL==out_engy) error(2, "(class_events) the          energy file was not opened!");
			if(NULL==out_vdep) error(2, "(class_events) the          vdepth file was not opened!");
		
			sum_mag= 2*nAA + 1*nA -1*nB -2*nBB;

			init_par();
			init_uncorrH();
		}
		
	private:
		char type_ltc[4];	// crystal structure type

		// functions
		void ltc_constructor();
		void init_states_array(double compV, double compA, int nMlayer);
		void read_restart(char name_restart[], long long int &ts_initial, double &time_initial);
		void init_par();
		void init_uncorrH();
};
#endif // KMC_SYSTEM_INCLUDED
