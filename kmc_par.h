#ifndef KMC_PAR_INCLUDED
#define KMC_PAR_INCLUDED

const char   par_ltc[4]=		"FCC";

const int    par_nx=                      276;
const int    par_ny=                       64;
const int    par_nz=                       64;
const int    par_nMlayer=                  10;

const double par_compA=                   0.9; // composition of A atoms
const int    par_nV=                        0;

const double	    par_time=            1e09; 	// toal time (s)
const long long int par_step=            5e09; 	// toal timestep (give a minus step to ignore this quiterior to end the simulation)

const long long int step_log=	          1e3; 
const long long int step_conf=   par_step/5e0;	// timestep that output a conf file for restart later
const long long int time_conf=   par_time/5e0;	//     time that output a conf file for restart later
const long long int step_out=	 par_step/5e2; 
const long long int step_his=             1e6;

const bool par_isrestart=		    false;

const char   par_name_sol[20]=      "history.sol";
const char   par_name_def[20]=      "history.def";
const char   par_name_srf[20]=      "history.srf";
const char   par_name_engy[20]=      "out.energy";

// Parameters for events
const double par_dis_rec=0.866*3; // recombination distance
const double par_dpasm1=     0.1; 

// Ising model energy calculation parameters
const double par_temp=                       700.0;
const double par_beta= 1.0/par_temp/8.617332478e-5; // 1/kbT, units: eV, K

const double par_muvA=			  1.68e+15; // units: s^-1     !!!!! A: Au; B: Cu !!!!!
const double par_muvB=			  4.20e+14; // units: s^-1 
const double par_muiA=			   5.0e+12; // units: s^-1 
const double par_muiB=			   5.0e+12; // units: s^-1 

const double par_emvA=			      0.88;
const double par_emvB=			      0.76;
const double par_emiA=			         0; 
const double par_emiB=			         0;

const double par_erAA=                   0; // rotation energy barrier
const double par_erAB=                   0;
const double par_erBB=                   0;

const double par_emiAA=			     0.300; // !!! these three para only for Dubey, CMS 2015
const double par_emiAB=			     0.377; // !!! 
const double par_emiBB=			     0.120; // !!!
const double par_eciAAtAB=             0.3; // !!!
const double par_eciABtAA=             0.5; // !!!
const double par_eciABtBB=            0.32; // !!!
const double par_eciBBtAB=            0.12; // !!!

// bonding energy parameters

// 1st nn
const double eAA1AA=                   0.25000;
const double eAA1A=                    0.24625;
const double eAA1AB=                   0.25000;
const double eAA1B=                    0.24625;
const double eAA1BB=                   0.25000;
// ---
const double eA1A=                    -0.14250;
const double eA1V=                    -0.01625;
const double eA1AB=                    0.12875;
const double eA1B=                    -0.14250;
const double eA1BB=                    0.14625;
// ---
const double eV1V=                           0;
const double eV1B=                    -0.01625;
// ---
const double eAB1AB=                   0.25000;
const double eAB1B=                    0.12875;
const double eAB1BB=                   0.25000;
// ---
const double eB1B=                    -0.14250;
const double eB1BB=                    0.14625;
// ---
const double eBB1BB=                   0.25000;
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
const double eA2A=                           0;
const double eA2V=                           0;
const double eA2AB=                          0;
const double eA2B=                           0;
const double eA2BB=                          0;
// ---
const double eV2V=                           0;
const double eV2B=                           0;
// ---
const double eAB2AB=                         0;
const double eAB2B=                          0;
const double eAB2BB=                         0;
// ---
const double eB2B=                           0;
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
