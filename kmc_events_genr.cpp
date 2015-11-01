#include <cstdio>
#include <iostream>
#include <vector>
#include "kmc_global.h"
#include "kmc_events.h"

using namespace std;

void class_events::genr(){
	int ltcp[2]; // the chosen 2 ltcps in generating frenkel pair; [0]->itl; [1]->vcc
    double ran;

	do{ // generating itl position: not on ltc occupied by AA, AB, BB, surface atom
		ran= ran_generator();
		ltcp[0]= (int) (ran*nx*ny*nz);
	}while( (*(&states[0][0][0]+ltcp[0]) !=-1 && *(&states[0][0][0]+ltcp[0]) !=1) || *(&srf[0][0][0]+ltcp[0]) );
	do{ // generating vcc position: not on ltc occupied by AA, AB, BB, surface atom, and ltcp[0]
		ran= ran_generator();
		ltcp[1]= (int) (ran*nx*ny*nz);
	}while( (*(&states[0][0][0]+ltcp[1]) !=-1 && *(&states[0][0][0]+ltcp[1]) !=1) || *(&srf[0][0][0]+ltcp[1]) || ltcp[0]==ltcp[1]);

	int iid= list_itl.size();
	int vid= list_vcc.size();

	// initialize the itl in the list_itl
	list_itl.push_back(itl());
	
	list_itl[iid].ltcp= ltcp[0]; 
	ran= ran_generator();
	list_itl[iid].dir= (int) (ran*n1nbr);
	list_itl[iid].head= *(&states[0][0][0]+ltcp[0]); // choose the ltcp[0] because dir is randomly selected
	list_itl[iid].ix= 0; 
	list_itl[iid].iy= 0; 
	list_itl[iid].iz= 0; 

	// initialize the vcc in the list_vcc
	list_vcc.push_back(vcc());
	
	list_vcc[vid].ltcp= ltcp[1];
	list_vcc[vid].ix= 0;
	list_vcc[vid].iy= 0;
	list_vcc[vid].iz= 0;

	// Update numbers (before)
	switch(*(&states[0][0][0]+ltcp[0])){
		case  1: nA --; break;
		case -1: nB --; break;
		default: error(2, "(genr) an unknown atom type", 1, *(&states[0][0][0]+ltcp[0]));
	}
	switch(*(&states[0][0][0]+ltcp[1])){
		case  1: nA --; break;
		case -1: nB --; break;
		default: error(2, "(genr) an unknown atom type", 1, *(&states[0][0][0]+ltcp[1]));
	}
	
	// Update states
	*(&states[0][0][0]+ltcp[0])= *(&states[0][0][0]+ltcp[0]) + *(&states[0][0][0]+ltcp[1]);
	*(&states[0][0][0]+ltcp[1])= 0;
	if(0==*(&states[0][0][0]+ltcp[0])) *(&itlAB[0][0][0]+ltcp[0])= true;

	// update numbers (after)
	switch(*(&states[0][0][0]+ltcp[0])){
		case  2: nAA ++; break;
		case  0: nAB ++; break;
		case -2: nBB ++; break;
		default: error(2, "(genr) an unknown itl type", 1, *(&states[0][0][0]+ltcp[0]));
	}
	nV ++;

    genr_1strecb(iid, vid);
}

void class_events::genr_1strecb(int iid, int vid){
	vector <int> list_rec1; // recombination candidates for itl (searching for vcc)
	vector <int> list_rec2; //                              vcc (searching for itl)

	int xi= (int) (list_itl[iid].ltcp/nz)/ny; // The itl position
	int yi= (int) (list_itl[iid].ltcp/nz)%ny;
	int zi= (int)  list_itl[iid].ltcp%nz;
	int xv= (int) (list_vcc[vid].ltcp/nz)/ny; // The vcc position
	int yv= (int) (list_vcc[vid].ltcp/nz)%ny;
	int zv= (int)  list_vcc[vid].ltcp%nz;

	for(int a=0; a<n1nbr; a ++){
		int x1= pbc(xi+v1nbr[a][0], nx);
		int y1= pbc(yi+v1nbr[a][1], ny);
		int z1= pbc(zi+v1nbr[a][2], nz);
        if(states[x1][y1][z1]== 0 && ! itlAB[x1][y1][z1])
            list_rec1.push_back(x1*ny*nz+y1*nz+z1);
		
        int x2= pbc(xv+v1nbr[a][0], nx);
		int y2= pbc(yv+v1nbr[a][1], ny);
		int z2= pbc(zv+v1nbr[a][2], nz);
        if((x2*ny*nz+y2*nz+z2)== list_itl[iid].ltcp) continue; // avoid put iltcp into list_rec2(list_rec1 would do)
        if(states[x2][y2][z2]== 2 || states[x2][y2][z2]==-2 || itlAB[x2][y2][z2]) list_rec2.push_back(x2*ny*nz+y2*nz+z2);
    }

    bool has_recb= false;
    if(list_rec1.size() != 0){
		double ran= ran_generator();
		int vltcp_recb= list_rec1[(int) (ran*list_rec1.size())];
        
        for(int i= 0; i<list_vcc.size(); i ++){ // brutal search for the vcc in the list
    		if(vltcp_recb==list_vcc[vid].ltcp) has_recb= true;

            if(vltcp_recb==list_vcc[i].ltcp){
                rules_recb(false, iid, i);
                break;
            }
        }
	}

    if(list_rec2.size() != 0 && ! has_recb){
		double ran= ran_generator();
		int iltcp_recb= list_rec2[(int) (ran*list_rec2.size())];
        
        for(int i= 0; i<list_itl.size(); i ++){ // brutal search
    		if(iltcp_recb==list_itl[i].ltcp){
                rules_recb(false, i, vid);
                break;
            }
        }
	}
}
