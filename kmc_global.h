#ifndef KMC_GLOBAL_INCLUDED
#define KMC_GLOBAL_INCLUDED
#include <cstring>
#include <vector>
#include "kmc_par.h"

using namespace std;

// Time parameters
extern long long int timestep;
extern double totaltime;

// ltc variables
#define MAX_NNBR 20
extern int n1nbr, n2nbr, n3nbr;	// number of neighbors
extern int v1nbr[MAX_NNBR][3];	// indexes vectors of 1st neighbors
extern int v2nbr[MAX_NNBR][3];	// indexes vectors of 2nd neighbors
extern int v3nbr[MAX_NNBR][3];	// indexes vectors of 3rd neighbors
extern double vbra[3][3];	// coordinate vectors of bravice lattice
extern int n1sp, n2sp;
extern int v1sp[MAX_NNBR][MAX_NNBR][3]; // [index of jump nbr][index of sp nbr][xyz]
extern int v2sp[MAX_NNBR][MAX_NNBR][3];
extern int n12nbr, n123nbr;
extern vector < vector<int> > v12nbr;
extern vector < vector<int> > v123nbr;

// System
const int nx=  par_nx;
const int ny=  par_ny;
const int nz=  par_nz;
const int x_sink= (int) (nx/2);
extern int nA, nB, nV, nAA, nBB, nAB, nM;
extern int sum_mag; // sum of magnitization; should be conserved
extern int  states[nx][ny][nz];
extern bool    srf[nx][ny][nz];

// Files
extern FILE * his_sol;		// history file of solute atoms
extern FILE * his_def;		// history file of vacancies and interstitials: record every several steps
extern FILE * his_srf;		// history file of surface atoms 
extern FILE * out_engy;		// out file for energy calculations
extern FILE * out_vdep;		// out file for store how long created vcc go
extern FILE * out_sro;		// out file of sro
extern FILE * out_msd;		// out file of msd

// Parameters for mechanisms
const double dis_rec=       par_dis_rec; // recombination distance
const bool   is_genr=       par_isgenr;
const bool   trap_included= par_trap_included;
const double rate_genr=     par_dpasm1*(nx-2*par_nMlayer)*ny*nz;
const double corrfac=       par_corrfac;
const bool   iscaldsro=     par_iscaldsro;

// Defect lists
struct vcc{ // information of an vacancy
    vcc(): njump(-1) {}	
    int njump;
    int ltcp;
	int ix, iy, iz;
};
struct itl{ // information of an interstitial
    itl(): trapped(false) {}
	int ltcp;
	int dir; // direction
	int head; // the atom of the itl that in the front along the dir; useful for AB itl
	int ix, iy, iz;
    bool trapped; // in some cases(low T) AB-itl trapped around B or AB, make it fixed
};
extern vector <vcc> list_vcc;	// A list contains information of all vacancies
extern vector <itl> list_itl;  	// A list contains information of all interstitials
extern vector <int> list_sink;  // A list contains atoms in the sink

extern int N_genr;
extern int njump[10];
extern long long int Vja[2], Ija[2];
extern double init_sro;
extern double acc_dsroV; // accumulative sro change from vcc
extern double acc_dsroG; // accumulative sro change from genr
extern double acc_dsroRi; // accumulative sro change from recb
extern double acc_dsroRv; // accumulative sro change from recb

// Migration parameters
const double temp= par_temp; 
const double beta= par_beta;

const double muvA=  par_muvA; 
const double muvB=  par_muvB; 
const double muiA=  par_muiA; 
const double muiB=  par_muiB; 
const double muiAA= par_muiAA;
const double muiAB= par_muiAB;

const double emvA= par_emvA;
const double emvB= par_emvB;
const double emiA= par_emiA;
const double emiB= par_emiB;

const double emiAA= par_emiAA;       // !!! only for DubeyCMS2015
const double emiAB= par_emiAB;       // !!!
const double emiBB= par_emiBB;       // !!!
const double eciAAtAB= par_eciAAtAB; // !!!
const double eciABtAA= par_eciABtAA; // !!!
const double eciABtBB= par_eciABtBB; // !!!
const double eciBBtAB= par_eciBBtAB; // !!!

const double erAA= par_erAA; // rotation energy barrier
const double erAB= par_erAB;
const double erBB= par_erBB;

const double eSPA1A= par_eSPA1A;
const double eSPA1B= par_eSPA1B;
const double eSPA1V= par_eSPA1V;
const double eSPA1M= par_eSPA1M;
const double eSPB1A= par_eSPB1A;
const double eSPB1B= par_eSPB1B;
const double eSPB1V= par_eSPB1V;
const double eSPB1M= par_eSPB1M;
const double eSPA2A= par_eSPA2A;
const double eSPA2B= par_eSPA2B;
const double eSPA2V= par_eSPA2V;
const double eSPA2M= par_eSPA2M;
const double eSPB2A= par_eSPB2A;
const double eSPB2B= par_eSPB2B;
const double eSPB2V= par_eSPB2V;
const double eSPB2M= par_eSPB2M;

// Ising model energy constants
extern double h0;
// Corrected H
extern double c1_44, c1_43, c1_42, c1_41, c1_33, c1_32, c1_31, c1_22, c1_21, c1_11, c1_40, c1_30, c1_20, c1_10;
extern double c1_00, c1_0_ABA, c1_0_ABB, c1_0_A, c1_0_B, c1_0_V, c1_0_AA, c1_0_AB, c1_0_BB;
extern double c2_44, c2_43, c2_42, c2_41, c2_33, c2_32, c2_31, c2_22, c2_21, c2_11, c2_40, c2_30, c2_20, c2_10;
extern double c2_00, c2_0_ABA, c2_0_ABB, c2_0_A, c2_0_B, c2_0_V, c2_0_AA, c2_0_AB, c2_0_BB;
extern double c1_MM, c1_MA, c1_MV, c1_MB;
extern double c2_MM, c2_MA, c2_MV, c2_MB;
// Uncorrected H
extern double unc1_44, unc1_43, unc1_42, unc1_41, unc1_33, unc1_32, unc1_31, unc1_22, unc1_21, unc1_11, unc1_40, unc1_30, unc1_20, unc1_10, unc1_00;
extern double unc2_44, unc2_43, unc2_42, unc2_41, unc2_33, unc2_32, unc2_31, unc2_22, unc2_21, unc2_11, unc2_40, unc2_30, unc2_20, unc2_10, unc2_00;

extern bool is_e2nbr;

// Global functions
extern void error(int nexit, string errinfo, int nnum=0, double num1=0, double num2=0);
extern void error(int nexit, string errinfo, char c[]);
extern double ran_generator();
extern int pbc(int x_, int nx_);
void write_conf(int flag);
void write_hissol();
void write_hisdef();
void write_metrohis();
void write_rcvcc();
void write_vdep(); // record how far created vcc go
int cal_Bnbr(int N_Bnbr, int x, int y, int z);
double cal_sro();
double cal_msd();
double cal_dsro(int xv, int yv, int zv, int xa, int ya, int za); // cal sro change during vcc jump

#endif // KMC_GLOBAL_INCLUDED
