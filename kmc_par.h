#ifndef KMC_PAR_INCLUDED
#define KMC_PAR_INCLUDED

const char   par_ltc[4]=		"BCC";

const int    par_nx=                       64; // 620 + 2*20
const int    par_ny=                       64;
const int    par_nz=                       64;
const int    par_nMlayer=                   0;

const int    par_Nv=                     1000;
const double par_mu=                      0.6; // chemical potential

const double	    par_time=           1e300; // toal time (s)
const double        time_conf=          1e300; // time that output a conf file for restart later

const long long int par_step=             5e6; // toal timestep (give a minus step to ignore this quiterior to end the simulation)
const long long int step_log=	          1e3; 
const long long int step_conf=            1e6; // timestep that output a conf file for restart later
const long long int step_out=	          1e4; 
const long long int step_his=             1e6;

const bool par_isrestart=		    false;
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

// Parameters for events
const double par_dis_rec=0.866*3; // recombination distance
const double par_dpasm1=     1.0; 

// Ising model energy calculation parameters
const double par_temp=                      1500.0;
const double par_beta= 1.0/par_temp/8.617332478e-5; // 1/kbT, units: eV, K

const double par_nuvA=			   6.0e+12; // units: s^-1     !!!!! A: Au; B: Cu !!!!!
const double par_nuvB=			   6.0e+12; // units: s^-1 
const double par_nuiA=			         0; // units: s^-1 
const double par_nuiB=			         0; // units: s^-1 

const double par_emvA=			      1.54;
const double par_emvB=			      1.57;
const double par_emiA=			         0; 
const double par_emiB=			         0;

const double par_erAA=                   0; // rotation energy barrier
const double par_erAB=                   0;
const double par_erBB=                   0;

const double par_emiAA=			         0; // !!! these three para only for Dubey, CMS 2015
const double par_emiAB=			         0; // !!! 
const double par_emiBB=			         0; // !!!
const double par_eciAAtAB=               0; // !!!
const double par_eciABtAA=               0; // !!!
const double par_eciABtBB=               0; // !!!
const double par_eciBBtAB=               0; // !!!

// bonding energy parameters
// AB bond & ratio 2/1
const double w0AB=                 -0.01963602;
const double w1AB=                 -0.02888548;
const double r21=                          0.5;

// 1st nn
const double eAA1AA=                         0;
const double eAA1A=                          0;
const double eAA1AB=                         0;
const double eAA1B=                          0;
const double eAA1BB=                         0;
// ---
const double eA1A=          -12.5590/(4+3*r21);
const double eA1V=                           0;
const double eA1AB=                          0;
const double eA1B=                           0;
const double eA1BB=                          0;
// ---
const double eV1V=                           0;
const double eV1B=                           0;
// ---
const double eAB1AB=                         0;
const double eAB1B=                          0;
const double eAB1BB=                         0;
// ---
const double eB1B=          -11.0662/(4+3*r21);
const double eB1BB=                          0;
// ---
const double eBB1BB=                         0;
// ...
const double eM1M=                           0; //vacuum
const double eM1A=                           0;
const double eM1V=                           0;
const double eM1B=                           0;

// 2nd nn
const double eAA2AA=                         0;
const double eAA2A=                          0;
const double eAA2AB=                         0;
const double eAA2B=                          0;
const double eAA2BB=                         0;
// ---
const double eA2A=                    eA1A*r21;
const double eA2V=                    eA1V*r21;
const double eA2AB=                          0;
const double eA2B=                    eA1B*r21;
const double eA2BB=                          0;
// ---
const double eV2V=                    eV1V*r21;
const double eV2B=                    eV1B*r21;
// ---
const double eAB2AB=                         0;
const double eAB2B=                          0;
const double eAB2BB=                         0;
// ---
const double eB2B=                    eB1B*r21;
const double eB2BB=                          0;
// ---
const double eBB2BB=                         0;
// ...
const double eM2M=                           0; //vacuum
const double eM2A=                           0;
const double eM2V=                           0;
const double eM2B=                           0;

// trapping number: solute atom trapping and intersitial trapping							 
// 	const int par_trNsol= 3; // trapping number by solute atoms
// 	const int par_trNint= 3; // trapping number by interstitials
#endif // KMC_PAR_INCLUDED
