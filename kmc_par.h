#ifndef KMC_PAR_INCLUDED
#define KMC_PAR_INCLUDED

const char   par_ltc[4]=		"BCC";

const int    par_nx=                       64; // 620 + 2*20
const int    par_ny=                       64;
const int    par_nz=                       64;
const int    par_nMlayer=                   0;

const double par_compA=                   0.9; // composition of A atoms
const double par_compV=                     0; // vcc; set >1.0 to get only 1 vcc

const double	    par_time=           1e300; // toal time (s)
const double        time_conf=          1e300; // time that output a conf file for restart later

const long long int par_step=             1e7; // toal timestep (give a minus step to ignore this quiterior to end the simulation)
const long long int step_log=	          1e3; 
const long long int step_conf=            1e6; // timestep that output a conf file for restart later
const long long int step_out=	          1e3; 
const long long int step_his=             1e5;

const bool par_isrestart=		    false;
const bool par_isnoflckr=		    false; // use a straight jumps for cvcc (no flickering)
const int  par_N_NFjumps=		        1; // number of straight jumps (no flickering)
const double par_NFratio=		        1; // ratio for the cvcc to not flickering (no flickering)

const char   par_name_sol[20]=      "history.sol";
const char   par_name_def[20]=      "history.def";
const char   par_name_srf[20]=      "history.srf";
const char   par_name_engy[20]=      "out.energy";
const char   par_name_vdep[20]=      "out.vdepth";
const char   par_name_sro[20]=      "out.sro";

// Parameters for events
const double par_dis_rec= 0.866*3; // recombination distance
const bool   par_isgenr=    false;
const double par_dpasm1=      1.0;
const double par_corrfac=     0.7;

// Ising model energy calculation parameters
const double par_temp=                      1500.0;
const double par_beta= 1.0/par_temp/8.617332478e-5; // 1/kbT, units: eV, K

const double par_muvA=			   6.0e+12; // units: s^-1     !!!!! A: Au; B: Cu !!!!!
const double par_muvB=			   6.0e+12; // units: s^-1 
const double par_muiA=			         0; // units: s^-1 
const double par_muiB=			         0; // units: s^-1 
const double par_muiAA=            6.0e+12;
const double par_muiAB=            6.0e+12;

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

// Saddle-point parameters    
const double par_eSPA1A= 0;
const double par_eSPA1B= 0;
const double par_eSPA1V= 0;
const double par_eSPB1A= 0;
const double par_eSPB1B= 0;
const double par_eSPB1V= 0;
const double par_eSPA2A= 0;
const double par_eSPA2B= 0;
const double par_eSPA2V= 0;
const double par_eSPB2A= 0;
const double par_eSPB2B= 0;
const double par_eSPB2V= 0;

// bonding energy parameters
// AB bond & ratio 2/1
const double w0AB=                 -0.01963602; // change it!
const double w1AB=                 -0.02888548; // chgnge it!
const double r21=                          0.5;

// 1st nn
const double eAA1AA=                         1;
const double eAA1A=                          2;
const double eAA1AB=                         3;
const double eAA1B=                          4;
const double eAA1BB=                         5;
// ---
const double eA1A=      -8.327616192/(4+3*r21);
const double eA1V=                     -0.3320;
const double eA1AB=                          6;
const double eA1B=                           7;
const double eA1BB=                          8;
// ---
const double eV1V=                           9;
const double eV1B=                     -0.4000;
// ---
const double eAB1AB=                        10;
const double eAB1B=                         11;
const double eAB1BB=                        12;
// ---
const double eB1B=           -7.4070/(4+3*r21);
const double eB1BB=                         13;
// ---
const double eBB1BB=                        14;
// ...
const double eM1M=                           0; //vacuum
const double eM1A=                           0;
const double eM1V=                           0;
const double eM1B=                           0;

// 2nd nn
const double eAA2AA=                        15;
const double eAA2A=                         16;
const double eAA2AB=                        17;
const double eAA2B=                         18;
const double eAA2BB=                        19;
// ---
const double eA2A=                    eA1A*r21;
const double eA2V=                     -0.3320;
const double eA2AB=                         20;
const double eA2B=                          21;
const double eA2BB=                         22;
// ---
const double eV2V=                          23;
const double eV2B=                        -0.4;
// ---
const double eAB2AB=                        24;
const double eAB2B=                         25;
const double eAB2BB=                        26;
// ---
const double eB2B=                    eB1B*r21;
const double eB2BB=                         27;
// ---
const double eBB2BB=                        28;
// ...
const double eM2M=                           0; //vacuum
const double eM2A=                           0;
const double eM2V=                           0;
const double eM2B=                           0;

// trapping number: solute atom trapping and intersitial trapping							 
// 	const int par_trNsol= 3; // trapping number by solute atoms
// 	const int par_trNint= 3; // trapping number by interstitials
#endif // KMC_PAR_INCLUDED
