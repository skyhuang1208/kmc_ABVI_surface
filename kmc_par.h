#ifndef KMC_PAR_INCLUDED
#define KMC_PAR_INCLUDED

const char   par_ltc[4]=		"BCC";

const int    par_nx=                       64; // 620 + 2*20
const int    par_ny=                       64;
const int    par_nz=                       64;
const int    par_nMlayer=                   0;

const int    par_Nv=                     1000;
const double par_mu=                     0.16; // chemical potential

const double	    par_time=           1e300; // toal time (s)
const double        time_conf=          1e300; // time that output a conf file for restart later

const long long int par_step=             2e9; // toal timestep (give a minus step to ignore this quiterior to end the simulation)
const long long int step_log=	          1e5; 
const long long int step_conf=            5e8; // timestep that output a conf file for restart later
const long long int step_out=	          1e5; 
const long long int step_his=             1e5;
const int           tcpu_his=            7200;

const bool par_isrestart=		     true;
const bool par_isreshis=            false; // restart from history file
const bool par_isnoflckr=		    false; // use a straight jumps for cvcc (no flickering)
const int  par_N_NFjumps=		        1; // number of straight jumps (no flickering)
const double par_NFratio=		        1; // ratio for the cvcc to not flickering (no flickering)

const char   par_name_sol[20]=      "history.sol";
const char   par_name_def[20]=      "history.def";
const char   par_name_srf[20]=      "history.srf";
const char   par_name_engy[20]=      "out.energy";
const char   par_name_vdep[20]=      "out.vdepth";
const char   par_name_eSGC[20]=        "out.eSGC";
const char   par_name_sro[20]=          "out.sro";
const char   par_name_log[20]=          "out.log";

// Parameters for events
const double par_dis_rec=0.866*3; // recombination distance
const double par_dpasm1=     1.0; 

// Ising model energy calculation parameters
const double par_temp=                      1200.0;
const double par_beta= 1.0/par_temp/8.617332478e-5; // 1/kbT, units: eV, K

// bonding energy parameters
// AB bond & ratio 2/1
const double r21=                     0.421875;

const double e0A1B= -1.5090; // eA1B: e0AB + e1AB*XB
const double e1A1B= -0.0219;
const double e0A2B= -0.6366; // eA1B: e0AB + e1AB*XB
const double e1A2B= -0.0092;

const double e0B1V= -0.4898; // eB1V: e0B1V + e1B1V/XB
const double e1B1V= -0.009432;
const double e0B2V= -0.3311; // eB2V: e0B2V + e1B2V*XB
const double e1B2V=  0.036;

// 1st nn
const double eA1A=                     -1.5815;
const double eB1B=                     -1.4067;

const double eA1V=                     -0.4898;
const double eV1V=                      0.5873;

// 2nd nn
const double eA2A=                     -0.6672;
const double eB2B=                     -0.5935;

const double eA2V=                     -0.2067;
const double eV2V=                      0.5566;

#endif // KMC_PAR_INCLUDED
