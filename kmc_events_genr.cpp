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
	}while( (*(&states[0][0][0]+ltcp[1]) !=-1 && *(&states[0][0][0]+ltcp[1]) !=1) || *(&srf[0][0][0]+ltcp[1]) || 
            ltcp[0]==ltcp[1] || -2==(*(&states[0][0][0]+ltcp[0])+*(&states[0][0][0]+ltcp[1])) );

	int iid= list_itl.size();
	int vid= list_vcc.size();

	// initialize the itl in the list_itl
	list_itl.push_back(itl());
	
	list_itl[iid].ltcp= ltcp[0]; 
	list_itl[iid].dir= (int) (ran_generator()*n1nbr);
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
	if(0==*(&states[0][0][0]+ltcp[0])) *(&states[0][0][0]+ltcp[0])= 3; // AB itl

	// update numbers (after)
	switch(*(&states[0][0][0]+ltcp[0])){
		case  2: nAA ++; break;
		case  3: nAB ++; break;
		default: error(2, "(genr) an unknown itl type", 1, *(&states[0][0][0]+ltcp[0]));
	}
	nV ++;
    
    recb_checki(iid);
    if(list_vcc.size() != 0) recb_checkv(list_vcc.size()-1);
}
