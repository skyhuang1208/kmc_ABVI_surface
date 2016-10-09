#ifndef KMC_PAR_INCLUDED
#define KMC_PAR_INCLUDED

const char   par_ltc[4]=		"BCC";

const int    par_nx=                       64; // 620 + 2*20
const int    par_ny=                       64;
const int    par_nz=                       64;
const int    par_nMlayer=                   0;

const double par_compA =                 0.85;
const double par_compV=                     0; // vcc; set >1.0 to get only 1 defect
const double par_typeD=                     0; // if sinlge defect, what is the type of defect

const double	    par_time=           100.0; // toal time (s)
const double        time_conf=           10.0; // time that output a conf file for restart later

const long long int par_step=             5e9; // toal timestep (give a minus step to ignore this quiterior to end the simulation)
const long long int step_log=	          1e5;
const long long int step_conf=            5e8; // timestep that output a conf file for restart later
const long long int step_out=	          1e6;
const long long int step_his=             1e6;

const bool par_isrestart=		    false;
const bool par_isOUTrestart=         true;

const bool par_isnoflckr=		    false; // NO USE: use a straight jumps for cvcc (no flickering)
const int  par_N_NFjumps=		        1; // NO USE: number of straight jumps (no flickering)
const double par_NFratio=		        1; // NO USE: ratio for the cvcc to not flickering (no flickering)

const char   par_name_sol[20]=      "history.sol";
const char   par_name_def[20]=      "history.def";
const char   par_name_srf[20]=      "history.srf";
const char   par_name_engy[20]=      "out.energy";
const char   par_name_vdep[20]=      "out.vdepth";
const char   par_name_sro[20]=          "out.sro";
const char   par_name_msd[20]=          "out.msd";

// Parameters for events
const double par_dis_rec=     0.866*3; // recombination distance
const bool   par_isgenr=         true;
const double par_dpasm1=         1e-3;
const bool   par_trap_included= false;

// Ising model energy calculation parameters
const double par_temp = 2100;
const double par_beta= 1.0/par_temp/8.617332478e-5; // 1/kbT, units: eV, K
const double par_corrfac=  2.93 - 0.00055*par_temp;

const double par_muvA=			  6.46e+12; // vcc mu and Em
const double par_muvB=			  6.46e+12;  
const double par_emvA=			     1.623;
const double par_emvB=			     1.651;

const double par_muiAA=           6.46e+12; // itl mu and Em and Er
const double par_muiAB=           6.46e+12;
const double par_emiAA=			     0.003;
const double par_emiAB=			      0.12; 
const double par_erAA=                0.43; 

// Saddle-point parameters    
const double par_eSPA1A=               -2.5975; // A-A
const double par_eSPA2A=               -0.5041;

const double par_eSPA1B=               -2.6451; // A-B
const double par_eSPA2B=               -0.5532;

const double par_eSPA1V=                0.5465; // A-V
const double par_eSPA2V=                0.1060;

const double par_eSPB1A=               -2.5188; // B-A
const double par_eSPB2A=               -0.4888;

const double par_eSPB1B=               -2.5417; // B-B
const double par_eSPB2B=               -0.4943;

const double par_eSPB1V=                0.2902; // B-V
const double par_eSPB2V=                0.0563;

// bonding energy parameters
// AB bond & ratio 2/1
const double r21=                     0.421875;

const double e0A1B=                    -1.5090; // eA1B: e0AB + e1AB*XB
const double e1A1B=                    -0.0219;
const double e0A2B=                    -0.6366; // eA1B: e0AB + e1AB*XB
const double e1A2B=                    -0.0092;

const double e0B1V=                    -0.4898; // eB1V: e0B1V + e1B1V/XB
const double e1B1V=                  -0.009432;
const double e0B2V=                    -0.3311; // eB2V: e0B2V + e1B2V*XB
const double e1B2V=                      0.036;

// bonds
const double eA1A=                     -1.5815;
const double eA2A=                     -0.6672;
const double eB1B=                     -1.4067;
const double eB2B=                     -0.5935;

const double eA1V=                     -0.4898;
const double eA2V=                     -0.2067;
const double eV1V=                      0.5873;
const double eV2V=                      0.5566;

const double eAA1A=                     0.1740;
const double eAA2A=                     0.0734;
const double eA1AB=                     0.1104;
const double eA2AB=                     0.0466;

const double eAA1B=                    -0.2750;
const double eAA2B=                    -0.1160;
const double eAB1B=                    -0.3486;
const double eAB2B=                    -0.1470;

const double eAA1AA=                   -0.1905;
const double eAA2AA=                   -0.0804;
const double eAA1AB=                   -0.2505;
const double eAA2AB=                   -0.1057;
const double eAB1AB=                   -1.3977;
const double eAB2AB=                   -0.5897;

const double eM1M=                           0; //vacuum
const double eM1A=                           0;
const double eM1V=                           0;
const double eM1B=                           0;
const double eM2M=                           0; //vacuum
const double eM2A=                           0;
const double eM2V=                           0;
const double eM2B=                           0;


// PARAMETERS THAT DON'T USE
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
const double par_eSPA1M=                 0;
const double par_eSPB1M=                 0;
const double par_eSPA2M=                 0;
const double par_eSPB2M=                 0;

const double eV1B=                           0;
const double eV2B=                           0;
const double eA1B=                           0;
const double eA2B=                           0;
const double eA1BB=                          0;
const double eA2BB=                          0;
const double eB1BB=                          0;
const double eB2BB=                          0;
const double eAA1BB=                         0;
const double eAA2BB=                         0;
const double eAB1BB=                         0;
const double eAB2BB=                         0;
const double eBB1BB=                         0;
const double eBB2BB=                         0;
#endif // KMC_PAR_INCLUDED
