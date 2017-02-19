#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <vector>
#include "kmc_global.h"

using namespace std;

////////// GLOBAL VARIABLES //////////
long long int timestep;
double totaltime;

#define MAX_NNBR 20
double vbra[3][3];
int n1nbr, n2nbr, n3nbr;
int v1nbr[MAX_NNBR][3];
int v2nbr[MAX_NNBR][3];
int v3nbr[MAX_NNBR][3];
int n1sp, n2sp;
int v1sp[MAX_NNBR][MAX_NNBR][3];
int v2sp[MAX_NNBR][MAX_NNBR][3];
int n12nbr, n123nbr;
int v12nbr[MAX_NNBR*2][3];
int v123nbr[MAX_NNBR*3][3];

int nA, nB, nV, nAA, nBB, nAB, nM;
int sum_mag;
vector<vector<vector<int>>> states(nx, vector<vector<int>>(ny, vector<int>(nz)));
vector<vector<vector<bool>>> srf(nx, vector<vector<bool>>(ny, vector<bool>(nz)));

FILE * his_sol;
FILE * his_def;
FILE * his_srf;
FILE * out_engy;
FILE * out_vdep;
FILE * out_sro;
FILE * out_msd;

vector <vcc> list_vcc;
vector <itl> list_itl;
vector <int> list_sink;

int N_genr= 0;
int njump[10]= {0};
long long int Vja[2]= {0};
long long int Ija[2]= {0};
double acc_dsroV=  0;
double acc_dsroG=  0;
double acc_dsroRi= 0;
double acc_dsroRv= 0;

double h0;
double c1_44, c1_43, c1_42, c1_41, c1_33, c1_32, c1_31, c1_22, c1_21, c1_11, c1_40, c1_30, c1_20, c1_10, c1_00;
double c1_0_ABA, c1_0_ABB, c1_0_A, c1_0_B, c1_0_V, c1_0_AA, c1_0_AB, c1_0_BB;
double c2_44, c2_43, c2_42, c2_41, c2_33, c2_32, c2_31, c2_22, c2_21, c2_11, c2_40, c2_30, c2_20, c2_10, c2_00;
double c2_0_ABA, c2_0_ABB, c2_0_A, c2_0_B, c2_0_V, c2_0_AA, c2_0_AB, c2_0_BB;
double c1_MM, c1_MA, c1_MV, c1_MB;
double c2_MM, c2_MA, c2_MV, c2_MB;
double unc1_44, unc1_43, unc1_42, unc1_41, unc1_33, unc1_32, unc1_31, unc1_22, unc1_21, unc1_11, unc1_40, unc1_30, unc1_20, unc1_10, unc1_00;
double unc2_44, unc2_43, unc2_42, unc2_41, unc2_33, unc2_32, unc2_31, unc2_22, unc2_21, unc2_11, unc2_40, unc2_30, unc2_20, unc2_10, unc2_00;
bool is_e2nbr;

////////// GLOBAL FUNCTIONS //////////
void error(int nexit, string errinfo, int nnum, double num1, double num2){
	// exit number represents:
	// 0: main; 1: class_system; 2: class_events
	
	cout << "\nError: ";
	switch(nexit){
		case 0:  cout << "In main function "; break;
		case 1:  cout << "In class_initial ";  break;
		case 2:  cout << "In class_events ";  break;
		default: cout << "UNKNOWN NEXIT" << endl;
	}

	cout << errinfo;
	switch(nnum){
		case 0:  cout << endl; break;
		case 1:  cout << ": " << num1 << endl; break;
		case 2:  cout << ": " << num1 << " " << num2 << endl; break;
		default: cout << "!!!ERROR FUNCTION MALFUNCTION!!! WRONG NNUM!!!" << endl;
	}
	cout << endl;
	exit(1); 
}

void error(int nexit, string errinfo, char c[]){
	// exit number represents:
	// 0: main; 1: class_system; 2: class_events
	
	cout << "\nError: ";
	switch(nexit){
		case 0:  cout << "In main function "; break;
		case 1:  cout << "In class_system ";  break;
		case 2:  cout << "In class_events ";  break;
		default: cout << "UNKNOWN NEXIT" << endl;
	}

	cout << errinfo << " " << c << endl;
	exit(1); 
}

double ran_generator(){
	static bool first= true;
	if(first){
		srand(time(NULL));
		first= false;
	}
	
	return rand()/((double) RAND_MAX+1.0);
}

int pbc(int x_, int nx_){ // Periodic Boundary Condition
	if	(x_<-nx_ || x_>=2*nx_)
		error(1, "(pbc) input x is out of bound", 2, x_, nx_);

	if	(x_<0)          return (x_ + nx_);
	else if	(x_<nx_)    return  x_;
	else                return (x_ - nx_);
}

void write_conf(int flag){
// A atoms are omitted
// flag: 0: t0; 1: timestep; 2: time; 3: RESTART
	ofstream of_xyz;
	ofstream of_ltcp;

	// determine the names of conf files
	char name_xyz[40], name_ltcp[40];
	
    if(0==flag){      // t0
		of_xyz.open("t0.xyz");
		of_ltcp.open("t0.ltcp");
	}
	else if(1==flag){ // step
		sprintf(name_xyz, "%lld.xyz", timestep);    of_xyz.open(name_xyz);
		sprintf(name_ltcp, "%lld.ltcp", timestep);  of_ltcp.open(name_ltcp);
	}
    else if(2==flag){ // time
		sprintf(name_xyz, "time%.2f.xyz", totaltime);   of_xyz.open(name_xyz);
		sprintf(name_ltcp, "time%.2f.ltcp", totaltime); of_ltcp.open(name_ltcp);
    }
    else if(3==flag) of_ltcp.open("RESTART");
    else error(0, "wrong flag: ", 1, flag);

	if(flag != 3 && !of_xyz.is_open()) error(1, "(write_conf) xyz file is not opened!");// check
	if(!of_ltcp.is_open()) error(1, "(write_conf) ltcp file is not opened!");           // check
	
	// write out data
	of_ltcp << nx*ny*nz-nA << "\n" << "ltcp " << timestep << " ";
	of_ltcp << setprecision(15) << totaltime << "\n";
    if(flag != 3){
	    of_xyz << nx*ny*nz-nA << "\n" << "xyz " << timestep << " ";
	    of_xyz << setprecision(15) << totaltime << "\n";
    }

	for(int i=0; i<nx; i ++){
		for(int j=0; j<ny; j ++){
			for(int k=0; k<nz; k ++){
				double x= i*vbra[0][0] + j*vbra[1][0] + k*vbra[2][0];
				double y= i*vbra[0][1] + j*vbra[1][1] + k*vbra[2][1];
				double z= i*vbra[0][2] + j*vbra[1][2] + k*vbra[2][2];
		
				if(1==states[i][j][k]) continue;
                else if(-1==states[i][j][k]){   // B atom
					if(flag != 3) of_xyz << states[i][j][k] << " " << x << " " << y << " " << z << endl;
					of_ltcp << states[i][j][k] << " " << i << " " << j << " " << k << " ";
				    if(srf[i][j][k]) of_ltcp << "1" << endl;
                    else             of_ltcp << "0" << endl;
                }
				else if(0==states[i][j][k]){    // vcc
					int id; for(id=0; list_vcc[id].x != i && list_vcc[id].y != j && list_vcc[id].z != k; id ++);
					if(flag != 3) of_xyz << states[i][j][k] << " " << x << " " << y << " " << z << " " << endl;
					of_ltcp << states[i][j][k] << " " << i << " " << j << " " << k << " " 
					        << list_vcc[id].ix << " " << list_vcc[id].iy << " " << list_vcc[id].iz << endl;
				}
                else if(4==states[i][j][k]){    // srf
					if(flag != 3) of_xyz << states[i][j][k] << " " << x << " " << y << " " << z << endl;
					of_ltcp << states[i][j][k] << " " << i << " " << j << " " << k << endl;
				}
				else{                           // itl
					int id; for(id=0; list_itl[id].x != i && list_itl[id].y != j && list_itl[id].z != k; id ++);
					int type= states[i][j][k];
					if(flag != 3) of_xyz  << type << " " << x << " " << y << " " << z << " " << endl; 
					of_ltcp << type << " " << i << " " << j << " " << k << " "
					        << list_itl[id].ix << " " << list_itl[id].iy << " " << list_itl[id].iz << " "
					        << list_itl[id].dir << " " << list_itl[id].head << endl;
				}
	}}}
	
	if(flag != 3) of_xyz.close();
	of_ltcp.close();
}

void write_hissol(){
	int ncheck= 0;
    vector <vector<int>> list_srf; // A list store srf info

    // OUTPUT his_sol
	fprintf(his_sol, "%d\n", nB);
	fprintf(his_sol, "T: %lld %e\n", timestep, totaltime);
	for(int i=0; i<nx; i++){
	    for(int j=0; j<ny; j++){
	        for(int k=0; k<nz; k++){
		        if(-1== states[i][j][k]){
			        ncheck ++;
                    int ltcp= i*ny*nz+j*nz+k;
			        fprintf(his_sol, "%d\n", ltcp);
                }

		        if(srf[i][j][k]) list_srf.push_back({i, j, k}); // know # before write
    }}}
	if(ncheck != nB) error(0, "(write_hissol) nB inconsistent", 2, ncheck, nB);
    
    // OUTPUT his_srf
    if(list_srf.size()==0) return;
	
    fprintf(his_srf, "%lu\n", list_srf.size());
	fprintf(his_srf, "T: %lld %e\n", timestep, totaltime);
    for(int i=0; i<list_srf.size(); i++){
        int x= list_srf[i][0];
        int y= list_srf[i][1];
        int z= list_srf[i][2];
		fprintf(his_srf, "%d %d %d %d\n", states[x][y][z], x, y, z);
	}

    fflush(his_sol);
    fflush(his_srf);
}

void write_hisdef(){
	// OUTPUT his_def
    fprintf(his_def, "%lu\n", list_vcc.size()+list_itl.size());
	fprintf(his_def, "T: %lld %e\n", timestep, totaltime);
    
    for(int i=0; i<list_vcc.size(); i++){
        int ltcp= list_vcc[i].x*ny*nz + list_vcc[i].y*nz + list_vcc[i].z;
		fprintf(his_def, "0 %d %d %d %d\n", ltcp, list_vcc[i].ix, list_vcc[i].iy, list_vcc[i].iz);
	}
	for(int i=0; i<list_itl.size(); i++){
        int ltcp= list_itl[i].x*ny*nz + list_itl[i].y*nz + list_itl[i].z;
		int type= states[list_itl[i].x][list_itl[i].y][list_itl[i].z];
		fprintf(his_def, "%d %d %d %d %d\n", type, ltcp, list_itl[i].ix, list_itl[i].iy, list_itl[i].iz);
	}
    
    fflush(his_def);
}

int cal_Bnbr(int N_Bnbr, int x, int y, int z){
    if(0==N_Bnbr){
        int n= 0;
        if(-1==states[x][y][z]) n ++;
        for(int a=0; a<n1nbr; a ++){ // search 1st-nn for B
		    int x1= pbc(x+v1nbr[a][0], nx);
		    int y1= pbc(y+v1nbr[a][1], ny);
		    int z1= pbc(z+v1nbr[a][2], nz);
		    if(-1==states[x1][y1][z1]) n ++;
	    }
	    for(int b=0; b<n2nbr; b ++){ // search 2nd-nn for B
		    int x2= pbc(x+v2nbr[b][0], nx);
		    int y2= pbc(y+v2nbr[b][1], ny);
		    int z2= pbc(z+v2nbr[b][2], nz);
		    if(-1==states[x2][y2][z2]) n ++;
        }

        return n;
    }
    else return N_Bnbr;
}

double cal_sro(){
    double cA= nA*1.0/(nx*ny*nz);
    double sro= 0;

    int ncheck= 0;
	for(int i=0; i<nx; i ++){
		for(int j=0; j<ny; j ++){
			for(int k=0; k<nz; k ++){
				int state0= states[i][j][k];

                if(-1==state0 || 3==state0){ // AB are considered B too
                    ncheck ++;
                    int nAnbr= 0;

				    for(int a=0; a<n1nbr; a ++){ // 1st neighbors
					    int x= pbc(i+(*(v1nbr+a))[0], nx);
					    int y= pbc(j+(*(v1nbr+a))[1], ny);
					    int z= pbc(k+(*(v1nbr+a))[2], nz);
                        int state1= states[x][y][z];

                        if(1==state1) nAnbr ++;
                    }
				    
                    for(int a=0; a<n2nbr; a ++){ // 1st neighbors
					    int x= pbc(i+(*(v2nbr+a))[0], nx);
					    int y= pbc(j+(*(v2nbr+a))[1], ny);
					    int z= pbc(k+(*(v2nbr+a))[2], nz);
                        int state2= states[x][y][z];

                        if(1==state2) nAnbr ++;
                    }

                    sro += 1.0 - nAnbr*1.0/(n1nbr+n2nbr)/cA;
                }
    }}}
    
    if(ncheck != (nAB+nB)) error(2, "(cal_sro) number inconsistent", 2, ncheck, nAB+nB);
    return sro/(nAB+nB);
}

/*
double cal_msd(){ // OBSELETE: USE DIFFUSION KMC CODE IF MSD IS NEEDED
    int i, j, k;
    if((! is_genr) && nV==1){
        i= list_vcc[0].ix*nx + ((int) (list_vcc[0].ltcp/nz)/ny); // vcc position
	    j= list_vcc[0].iy*ny + ((int) (list_vcc[0].ltcp/nz)%ny);
	    k= list_vcc[0].iz*nz + ((int)  list_vcc[0].ltcp%nz);
    }
    else if((! is_genr) && (nAA+nAB+nBB)==1){
        i= list_itl[0].ix*nx + ((int) (list_itl[0].ltcp/nz)/ny); // itl position
	    j= list_itl[0].iy*ny + ((int) (list_itl[0].ltcp/nz)%ny);
	    k= list_itl[0].iz*nz + ((int)  list_itl[0].ltcp%nz);
    }
    else return -1;

    double x= i*vbra[0][0] + j*vbra[1][0] + k*vbra[2][0];
	double y= i*vbra[0][1] + j*vbra[1][1] + k*vbra[2][1];
	double z= i*vbra[0][2] + j*vbra[1][2] + k*vbra[2][2];

    return x*x + y*y + z*z;
}

void write_vdep(){
    // format: step time total >=1 >=2 >=3 >=4 >=5 >=6 >=7 >=8 >=9
    // total 12 columns

    fprintf(out_vdep, "%lld %e ", timestep, totaltime);
    
    for(int a=0; a<10; a ++){
        int n=0;
        for(int b=a; b<10; b++){
            n += njump[b];
        }
        fprintf(out_vdep, "%d ", n);
    }
    fprintf(out_vdep, "\n");

    fflush(out_vdep);
}
*/

/*
// something's wrong with cal dsro method; check carefully if want to use
double cal_dsro(int xv, int yv, int zv, int xa, int ya, int za){ // cal sro change during vcc jump
    int stateV= states[xv][yv][zv];
    int stateA= states[xa][ya][za];
    if(stateV != 0)                 error(1, "(cal_dsro) supposed to be a vcc", 1, stateV);
    if(stateA != 1 && stateA != -1) error(1, "(cal_dsro) supposed to be a atom", 1, stateA);

    vector < vector<int> > list_cal(3);
    if(-1==stateA){ // A itself
        list_cal[0].push_back(xv); list_cal[0].push_back(xa);
        list_cal[1].push_back(yv); list_cal[1].push_back(ya);
        list_cal[2].push_back(zv); list_cal[2].push_back(za);
    }
    for(int a=0; a<n12nbr; a ++){ // V's neighbors
	    int x= pbc(xv+v12nbr[a][0], nx);
		int y= pbc(yv+v12nbr[a][1], ny);
		int z= pbc(zv+v12nbr[a][2], nz);
        if(-1==states[x][y][z]){
            list_cal[0].push_back(x);
            list_cal[1].push_back(y);
            list_cal[2].push_back(z);
        }
    }
    for(int a=0; a<n12nbr; a ++){ // A's neighbors
	    int x= pbc(xa+v12nbr[a][0], nx);
		int y= pbc(ya+v12nbr[a][1], ny);
		int z= pbc(za+v12nbr[a][2], nz);
        if(-1==states[x][y][z]){
            list_cal[0].push_back(x);
            list_cal[1].push_back(y);
            list_cal[2].push_back(z);
        }
    }

    double cA= nA*1.0/(nx*ny*nz);
    double dsro= 0;
    for(int a=-1; a<=1; a +=2){ // before(-1) & after(+1)
        for(int b= 0; b<list_cal[0].size(); b ++){
            int i= list_cal[0].at(b);
            int j= list_cal[1].at(b);
            int k= list_cal[2].at(b);
    	
            if(states[i][j][k] != -1) continue;

            int nAnbr= 0;
	        for(int c=0; c<n12nbr; c ++){ // 1st, 2nd neighbors
	            int x= pbc(i+v12nbr[c][0], nx);
		        int y= pbc(j+v12nbr[c][1], ny);
		        int z= pbc(k+v12nbr[c][2], nz);
                int state1= states[x][y][z];

                if(1==state1) nAnbr ++;
            }
    
            dsro += a * (1.0 - nAnbr*1.0/(n1nbr+n2nbr)/cA);
        }
        
        states[xv][yv][zv]= stateA; // do jump
        states[xa][ya][za]= stateV;
    } 
    states[xv][yv][zv]= stateV; // restore
    states[xa][ya][za]= stateA;
    
    return dsro / (nAB+nB);
}
*/

