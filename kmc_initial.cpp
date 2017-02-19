#include <cstdlib>
#include <cstring>
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include "kmc_initial.h"
#include "kmc_par.h"

using namespace std;

void class_initial::ltc_constructor(){
	double (*ptr_vbra)[3]; 
	int   (*ptr_v1nbr)[3];
	int   (*ptr_v2nbr)[3];
			
	// coordinate vectors of bravais lattices
		
	// BCC
	double vbra_bcc[3][3]= {{-0.5,  0.5,  0.5}, { 0.5, -0.5,  0.5}, { 0.5,  0.5, -0.5}};
	
    int   v1nbr_bcc[8][3]= {{ 1,  0,  0}, { 0,  1,  0}, { 0,  0,  1}, { 1,  1,  1},
                            {-1,  0,  0}, { 0, -1,  0}, { 0,  0, -1}, {-1, -1, -1}};
	
    int   v2nbr_bcc[6][3]= {{ 0,  1,  1}, { 1,  0,  1}, { 1,  1,  0},
                            { 0, -1, -1}, {-1,  0, -1}, {-1, -1,  0}};
	
    // FCC
	double vbra_fcc[3][3]= {{0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};

	int   v1nbr_fcc[12][3]= {{ 1,  0,  0}, { 0,  1,  0}, { 0,  0,  1}, { 1, -1,  0}, { 1,  0, -1}, { 0,  1, -1},
                             {-1,  0,  0}, { 0, -1,  0}, { 0,  0, -1}, {-1,  1,  0}, {-1,  0,  1}, { 0, -1,  1}};

	int   v2nbr_fcc[6][3]=  {{ 1,  1, -1}, { 1, -1,  1}, {-1,  1,  1},
                             {-1, -1,  1}, {-1,  1, -1}, { 1, -1, -1}};

	// Choose ltc structure
	if     (strcmp(type_ltc, "SC ")==0){}
	else if(strcmp(type_ltc, "BCC")==0){ptr_vbra= vbra_bcc;	n1nbr= 8; ptr_v1nbr= v1nbr_bcc; n2nbr= 6; ptr_v2nbr= v2nbr_bcc;}
	else if(strcmp(type_ltc, "FCC")==0){ptr_vbra= vbra_fcc; n1nbr=12; ptr_v1nbr= v1nbr_fcc; n2nbr= 6; ptr_v2nbr= v2nbr_fcc;}
	else if(strcmp(type_ltc, "HCP")==0){}
	else	error(1, "(ltc_constructor) coldn't match the lattice type", type_ltc);
			
	// assign array values
	for(int i=0; i<3; i ++){     // Bravice vectors
		for(int j=0; j<3; j ++){
			vbra[i][j]= (*(ptr_vbra+i))[j];
	}}
	for(int i=0; i<n1nbr; i ++){ // v1nbr
		for(int j=0; j<3; j ++){
			v1nbr[i][j]= (*(ptr_v1nbr+i))[j];
			
            if(i>=n1nbr/2)
				if(v1nbr[i][j] != -v1nbr[i-n1nbr/2][j]) error(1, "(ltc_constructor) v1nbr isn't symmetry");
		}
	}
	for(int i=0; i<n2nbr; i ++){ // v2nbr
		for(int j=0; j<3; j ++){
			v2nbr[i][j]= (*(ptr_v2nbr+i))[j];
			
            if(i>=n2nbr/2)
				if(v2nbr[i][j] != -v2nbr[i-n2nbr/2][j]) error(1, "(ltc_constructor) v1nbr isn't symmetry");
		}
	}
}

void class_initial::init_states_array(){
	// STATE 0: vacancy, 1: A atom, -1: B atom, 4: Vacuum

    for(int i=0; i<nx; i ++) // fill system with A
	    for(int j=0; j<ny; j ++)	
	        for(int k=0; k<nz; k ++)
                states[i][j][k]= 1;
    
    int count= 0;
    while(count != par_Nv){
        int site= ran_generator()*nx*ny*nz;
        if(*(&states[0][0][0]+site)!=0) count ++;
        *(&states[0][0][0]+site)= 0;
    }

	nV= 0; nA= 0; nB= 0; nAA= 0; nBB= 0; nAB= 0; nM= 0;
	
    ////////// CHECK //////////
	for(int i=0; i<nx; i ++){ 
	    for(int j=0; j<ny; j ++){ 
	        for(int k=0; k<nz; k ++){ 
		        switch(states[i][j][k]){
                    case  0:
                        nV ++; break;
		            case  1:
                        nA ++; break;
                    case -1:
                        nB ++; break;
                    case  4:
                        for(int a=0; a<n1nbr; a ++){ // mark srf atoms
                            int x= pbc(i+v1nbr[a][0], nx);
                            int y= pbc(j+v1nbr[a][1], ny);
                            int z= pbc(k+v1nbr[a][2], nz);

                            if(     states[x][y][z] == 0) error(1, "(init_states_array) vcc adjacent to vacuum");
                            else if(states[x][y][z] != 4) srf[x][y][z]= true;
                        } 
                        nM ++; break;
                    default: error(1, "(init_states_array) a state type is unrecognizable", 1, states[i][j][k]);
	            }
    }}}
	int nAtotal= nA + nB;
	////////// CHECK //////////
	
	cout << "The random solution configuration has been generated!" << endl;
	cout << "Vacancy: " << nV << endl;
    cout << "Vacuum:  " << nM << endl;
	cout << "Atype A: " << nA << ", pct: " << 100* (double) nA/(nAtotal) << "%" << endl;
	cout << "Atype B: " << nB << ", pct: " << 100* (double) nB/(nAtotal) << "%" << endl;
}

void class_initial::read_restart(char name_restart[], long long int &ts_initial, double &time_initial){
	ifstream if_re(name_restart, ios::in);
	if(!if_re.is_open()) error(1, "(read_restart) the file is not opened!");
	
	long long int timestep;
	double time;

	int ntotal;
	if_re >> ntotal;
//	if(ntotal != nx*ny*nz) error(1, "(read_restart) the input total ltc number isnt consistent", 2, ntotal, nx*ny*nz);
	
	char c_ltcp[5];
	if_re >> c_ltcp >> timestep >> time;
	if(strcmp(c_ltcp, "ltcp") !=0) error(1, "(read_restart) please input a ltcp file (ltcp at the second line)"); // check
	ts_initial= timestep;
	time_initial= time;

    for(int i=0; i<nx*ny*nz; i ++) *(&states[0][0][0]+i)= 1; // RESTART not contain A atoms; init with As

	nV= 0; nA= nx*ny*nz; nB= 0; nAA= 0; nBB= 0; nAB= 0;
	for(int index=0; index<ntotal; index ++){	
		int type, i, j, k, is_srf, ix, iy, iz, dir, head;
		if_re >> type >> i >> j >> k;
//		if(index != i*ny*nz+j*nz+k) error(1, "(read_restart) the input index inconsistent");
	
		if( 0==type){
			if_re >> ix >> iy >> iz;
			nV ++; nA --;
		}
        else if(1==type || -1==type){
		    if_re >> is_srf;
            
		    if(-1==type){ nB ++; nA --;}
        }
        else error(0, "(read_restart) a weird type", 1, type);

		int pid= i*ny*nz+j*nz+k; // position
		*(&states[0][0][0] + pid)= type;
	}
	if(nV+nA+nB+nAA+nBB+nAB+nM != nx*ny*nz) error(1, "(read_restart) the number inconsistent", 2, nV+nA+nB+nAA+nBB+nAB+nM, nx*ny*nz);

	cout << "The configuration has been generated from the restart file!" << endl;
	cout << "Vacancy: " << nV << endl;
	cout << "Atype A: " << nA << ", pct: " << 100* (double)nA /(nx*ny*nz) << "%" << endl;
	cout << "Atype B: " << nB << ", pct: " << 100* (double)nB /(nx*ny*nz) << "%" << endl;
	cout << " Itl AA: " << nAA << endl;
	cout << " Itl AB: " << nAB << endl;
	cout << " Itl BB: " << nBB << endl;
	
	if_re.close();
}

void class_initial::read_res_his(long long int &ts_initial, double &time_initial){ // read restart from history.sol
	error(0, "(read_res_his) do not read history.def, need modify");
    ifstream if_his(par_name_sol, ios::in);
	if(!if_his.is_open()) error(1, "(read_restart) the history.sol is not opened!");

	int ns;
	double time;
	long long int timestep;
    vector <int> sltcp;
	while(if_his >> ns){
		if_his.ignore();

	    char c_T[5];
	    if_his >> c_T >> timestep >> time;
        if(strcmp(c_T, "T:") !=0) error(1, "(read_restart) not T: on second line"); // check

        sltcp.clear();
        for(int a=0; a<ns; a ++){
            int data;
            if_his >> data;
            sltcp.push_back(data);
        }
    }
	ts_initial= timestep;
    time_initial= time;

    for(int i=0; i<nx; i ++) 
        for(int j=0; j<ny; j ++) 
            for(int k=0; k<nz; k ++) states[i][j][k]= 1;
    for(int a=0; a<ns; a ++) *(&states[0][0][0]+sltcp[a])= -1;
	
    nA= nx*ny*nz - ns; nB= ns; nV= 0; nAA= 0; nAB= 0; nBB= 0;
    cout << "Restart from history file (!! Warning: AB system only !!)" << endl;
	cout << "Vacancy: " << nV << endl;
	cout << "Atype A: " << nA << ", pct: " << 100* (double)nA/(nx*ny*nz) << "%" << endl;
	cout << "Atype B: " << nB << ", pct: " << 100* (double)nB/(nx*ny*nz) << "%" << endl;
	cout << " Itl AA: " << nAA << endl;
	cout << " Itl AB: " << nAB << endl;
	cout << " Itl BB: " << nBB << endl;

    if_his.close();
}

void class_initial::init_par(){
	is_e2nbr= true;

	// print out the parameters to log file
	cout << "\n##### Energy calculation parameters #####" << endl; 
	
	cout << "temperature= " << temp << ", beta= " << beta << endl;
	
	cout << "\n##### Input epsilons: #####" << endl;
	cout << "(1st neigbor)" << endl;
	printf("A-A: %f, B-B: %f, A-V: %f, V-V: %f\n", eA1A, eB1B, eA1V, eV1V);
    printf("A-B: %f+%f*XB, B-V: %f+%f/XB", e0A1B, e1A1B, e0B1V, e1B1V);
	cout << "(2nd neigbor)" << endl;
	printf("A-A: %f, B-B: %f, A-V: %f, V-V: %f\n", eA2A, eB2B, eA2V, eV2V);
    printf("A-B: %f+%f*XB, B-V: %f+%f*XB", e0A2B, e1A2B, e0B2V, e0B2V);
	
	if(is_e2nbr) cout << "\n2nd nn parameters are non-zero" << endl;
	else         cout << "\n2nd nn are 0, skip 2nd-nn calculations" << endl;
}

