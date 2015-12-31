#include <iostream>
#include <cmath>
#include <vector>
#include "kmc_global.h"
#include "kmc_events.h"
#include "kmc_par.h"

using namespace std;

void class_events::rules_recb(bool is_vcm, int ii, int iv, int ja){ // execute the recombination: vcc or vacuum
	int iltcp= list_itl[ii].ltcp;
	int vltcp;
    if(is_vcm)  vltcp= iv; // vacuum don't hv a list
    else        vltcp= list_vcc[iv].ltcp;

	if(*(&states[0][0][0]+vltcp) != 0 && *(&states[0][0][0]+vltcp) != 4) error(2, "(rules_recb) the vacancy isnt 0 or 4", 1, *(&states[0][0][0]+vltcp));
	double ran;
	switch(*(&states[0][0][0]+iltcp)){
		case  2:
			nAA --; nA +=2;
			ja= 1;
			break;
		case  0:
			nAB --; nA ++; nB ++;
			*(&itlAB[0][0][0]+iltcp)= false;
			if(0==ja){ // did not assign value when calling
			    ran= ran_generator();
                if(ran<0.5) ja= 1; 
			    else	    ja=-1;
            }    
            break;
		case -2:
			nBB --; nB +=2;
			ja=-1;
			break;
		default: error(2, "(rules_recb) an unknown itl type", 1, *(&states[0][0][0]+iltcp));
	}
			
	*(&states[0][0][0]+iltcp) -=ja;
	*(&states[0][0][0]+vltcp)  =ja; // should work for both vcc and vcm
	
	list_itl.erase(list_itl.begin()+ii);
    
    if(is_vcm){
        nM --;
        srf_check(vltcp);
    }
    else{
        nV --;
        list_vcc.erase(list_vcc.begin()+iv);
    }
}

void class_events::srf_check(int mltcp){ // when vacuum changed, check if srf array changes
    int i= (int) (mltcp/nz)/ny; // the ltcp where is used to be a vacuum
    int j= (int) (mltcp/nz)%ny;
    int k= (int)  mltcp%nz;

    for(int a= -1; a<n1nbr; a ++){ // here we check if some surface atoms become non-surface atoms (hv no bond with vacuum)
        int x, y, z;
        if(-1==a){x= i; y= j; z= k;} // check itself
        else{
            x= pbc(i+v1nbr[a][0], nx);
            y= pbc(j+v1nbr[a][1], ny);
            z= pbc(k+v1nbr[a][2], nz);
        }
        
        if(states[x][y][z] != 1 && states[x][y][z] != -1){
            srf[x][y][z]= false;
            continue;
        }

        bool is_srf= false;
        for(int b=0; b<n1nbr; b ++){
            int d= pbc(x+v1nbr[b][0], nx);
            int e= pbc(y+v1nbr[b][1], ny);
            int f= pbc(z+v1nbr[b][2], nz);
            if(4==states[d][e][f]){
                is_srf= true;
                break;
            }
        }

        if((is_srf != srf[x][y][z]) && is_noflckr) list_update.push_back(x*ny*nz+y*nz+z);
        
        srf[x][y][z]= is_srf; 
    }
}

// THE functions that have been deleted. To find them please go to ABVI_fixSINK //
// void class_events::recb_dir(int index){
// bool class_events::cal_dis(int d1, int d2, int d3){
// void class_events::recb_randomV(int index){
// bool class_events::recb_randomI(int index){
