#ifndef KMC_PAR_INCLUDED
#define KMC_PAR_INCLUDED

const char   par_ltc[4]=		"BCC";

const int    par_nx=                      256; // 620 + 2*20
const int    par_ny=                       64;
const int    par_nz=                       64;
const int    par_nMlayer=                   0;

const double par_compA=                  0.99; // composition of A atoms
const double par_compV=                     0; // vcc; set >1.0 to get only 1 vcc

const double	    par_time=           1e300; // toal time (s)
const double        time_conf=          1e300; // time that output a conf file for restart later

const long long int par_step=             1e7; // toal timestep (give a minus step to ignore this quiterior to end the simulation)
const long long int step_log=	          1e4; 
const long long int step_conf=            1e6; // timestep that output a conf file for restart later
const long long int step_out=	          1e4; 
const long long int step_his=             1e4;

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
const bool   par_isgenr=     true;
const double par_dpasm1=     1e-3;

// Ising model energy calculation parameters
const double par_temp=                       900.0;
const double par_beta= 1.0/par_temp/8.617332478e-5; // 1/kbT, units: eV, K
const double par_corrfac=  2.93 - 0.00055*par_temp;

const double par_muvA=			  6.46e+12; // vcc mu and Em
const double par_muvB=			  6.46e+12;  
const double par_emvA=			     1.699;
const double par_emvB=			     1.638;

const double par_muiAA=           6.46e+12; // itl mu and Em and Er
const double par_muiAB=           6.46e+12;
const double par_emiAA=			     0.003;
const double par_emiAB=			      0.12; 
const double par_erAA=                0.43; 

const double par_muiA=			         0; // parameters that is not using 
const double par_muiB=			         0;  
const double par_emiA=			         0; 
const double par_emiB=			         0;
const double par_erAB=                   0;
const double par_erBB=                   0;
const double par_emiBB=			         0; // !!! these three para only for Dubey, CMS 2015
const double par_eciAAtAB=               0; // !!!
const double par_eciABtAA=               0; // !!!
const double par_eciABtBB=               0; // !!!
const double par_eciBBtAB=               0; // !!!

// Saddle-point parameters    
const double par_eSPA1A= -2.4801; // (A-A, A-B, A-V), A-M
const double par_eSPA1B= -2.5923;
const double par_eSPA1V=  0.5539;
const double par_eSPA1M=  0.1640;
const double par_eSPB1A= -2.4155; // (B-A, B-B, B-V), B-M
const double par_eSPB1B= -2.4615;
const double par_eSPB1V=  0.3493;
const double par_eSPB1M=  0.1380;
const double par_eSPA2A= -0.62;   // (A-A, A-B, A-V), A-M
const double par_eSPA2B= -0.6491;
const double par_eSPA2V=  0.1385;
const double par_eSPA2M=  0.041;
const double par_eSPB2A= -0.6039; // (B-A, B-B, B-V), B-M
const double par_eSPB2B= -0.6154;
const double par_eSPB2V=  0.0873;
const double par_eSPB2M=  0.0345;

// bonding energy parameters
// AB bond & ratio 2/1
const double r21=                          0.5;
const double e0AB=                     -1.4447;
const double e1AB=                     -0.0210;
// BV bond
const double e0BV=  -0.6029; // eB1V= e0BV + e1BV*conc
const double e1BV=        0; // change

// 1st nn
const double eA1A=                     -1.5141;
const double eB1B=                     -1.3467;

const double eA1V=                     -0.4690;
const double eV1B=                           0;
const double eV1V=                      0.7614;

const double eAA1A=                     0.1666;
const double eA1AB=                     0.1057;
const double eAA1B=                    -0.2854;
const double eAB1B=                     0.1737;

const double eAA1AA=                   -0.2727;
const double eAA1AB=                   -0.3336;
const double eAB1AB=                   -1.4745;
// ---
// ---
const double eA1B=                           0;
const double eA1BB=                          0;
const double eB1BB=                          0;
const double eAA1BB=                         0;
const double eAB1BB=                         0;
const double eBB1BB=                         0;
// ...
const double eM1M=                           0; //vacuum
const double eM1A=                           0;
const double eM1V=                           0;
const double eM1B=                           0;

// 2nd nn
const double eA2A=                     -0.7571;
const double eB2B=                     -0.6734;

const double eA2V=                     -0.2345;
const double eV2B=                           0;
const double eV2V=                      0.7850;

const double eAA2A=                     0.0833;
const double eA2AB=                     0.0529;
const double eAA2B=                    -0.1427;
const double eAB2B=                     0.0829;

const double eAA2AA=                   -0.1364;
const double eAA2AB=                   -0.1668;
const double eAB2AB=                   -0.7373;
// ---
// ---
const double eA2B=                           0;
const double eA2BB=                          0;
const double eB2BB=                          0;
const double eAA2BB=                         0;
const double eAB2BB=                         0;
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
