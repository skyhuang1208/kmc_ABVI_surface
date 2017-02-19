#include <cstdio>
#include <iostream>
#include <vector>
#include "kmc_global.h"
#include "kmc_events.h"

using namespace std;

void class_events::genr(){
    int x1, y1, z1; // the chosen 2 ltcps in generating frenkel pair; 1->itl; 2->vcc
    int x2, y2, z2;
    double ran;

	do{ // generating itl position: not on ltc occupied by AA, AB, BB, surface atom
		x1= (int) (ran_generator()*nx);
		y1= (int) (ran_generator()*ny);
		z1= (int) (ran_generator()*nz);
	}while( (states[x1][y1][z1] !=-1 && states[x1][y1][z1] !=1) || srf[x1][y1][z1] );
	do{ // generating vcc position: not on ltc occupied by AA, AB, BB, surface atom, and ltcp[0], become BB
		x2= (int) (ran_generator()*nx);
		y2= (int) (ran_generator()*ny);
		z2= (int) (ran_generator()*nz);
	}while( (states[x2][y2][z2] !=-1 && states[x2][y2][z2] !=1) || srf[x2][y2][z2] || 
            (x1==x2 && y1==y2 && z1==z2) || -2==(states[x1][y1][z1]+states[x2][y2][z2]) );

	int iid= list_itl.size();
	int vid= list_vcc.size();

	// initialize the itl in the list_itl
	list_itl.push_back(itl());
	
	list_itl[iid].x= x1;
	list_itl[iid].y= y1;
	list_itl[iid].z= z1;
	list_itl[iid].dir= (int) (ran_generator()*n1nbr);
	list_itl[iid].head= states[x1][y1][z1]; // choose the ltcp[0] because dir is randomly selected
	list_itl[iid].ix= 0; 
	list_itl[iid].iy= 0; 
	list_itl[iid].iz= 0; 

	// initialize the vcc in the list_vcc
	list_vcc.push_back(vcc());
	
	list_vcc[vid].x= x2;
	list_vcc[vid].y= y2;
	list_vcc[vid].z= z2;
	list_vcc[vid].ix= 0;
	list_vcc[vid].iy= 0;
	list_vcc[vid].iz= 0;

	// Update numbers (before)
	switch(states[x1][y1][z1]){
		case  1: nA --; break;
		case -1: nB --; break;
		default: error(2, "(genr) an unknown atom type", 1, states[x1][y1][z1]);
	}
	switch(states[x2][y2][z2]){
		case  1: nA --; break;
		case -1: nB --; break;
		default: error(2, "(genr) an unknown atom type", 1, states[x2][y2][z2]);
	}
	
	// Update states
	states[x1][y1][z1]= states[x1][y1][z1] + states[x2][y2][z2];
	states[x2][y2][z2]= 0;
	if(0==states[x1][y1][z1]) states[x1][y1][z1]= 3; // AB itl

	// update numbers (after)
	switch(states[x1][y1][z1]){
		case  2: nAA ++; break;
		case  3: nAB ++; break;
		default: error(2, "(genr) an unknown itl type", 1, states[x1][y1][z1]);
	}
	nV ++;
    
    recb_checki(iid);
    if(list_vcc.size() != 0) recb_checkv(list_vcc.size()-1);
}
