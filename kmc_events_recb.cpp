#include <iostream>
#include <cmath>
#include <vector>
#include "kmc_global.h"
#include "kmc_events.h"
#include "kmc_par.h"

using namespace std;

void class_events::rules_recb(int ii, int xi, int yi, int zi, int iv, int xv, int yv, int zv){ // execute the recombination: vcc or vacuum
    if(-1==ii){
        int ltcp= xi*ny*nz + yi*nz + zi;
        for(int a=0; a<list_itl.size(); a ++) 
            if(list_itl[a].ltcp==ltcp) ii= a;
        if(-1==ii) error(2, "(rules_recb) cant find the itl");
    }
    if(-1==iv && 0==states[xv][yv][zv]){
        int ltcp= xv*ny*nz + yv*nz + zv;
        for(int a=0; a<list_vcc.size(); a ++)
            if(list_vcc[a].ltcp==ltcp) iv= a;
        if(-1==iv) error(2, "(rules_recb) cant find the vcc");
    }

    if(-1==states[xv][yv][zv]){
        if(states[xi][yi][zi] != 2) error(2, "(rules_recb) a recb with B atom not AA itl (type)", 1, states[xi][yi][zi]);
        nAA --; nB --;
        states[xi][yi][zi]= 1;
        states[xv][yv][zv]= 3;
        nAB ++; nA ++;
        list_itl[ii].ltcp= xv*ny*nz+yv*nz+zv;
        recb_checki(ii);
    }
    else{
	    if(states[xv][yv][zv] != 0 && states[xv][yv][zv] != 4) error(2, "(rules_recb) the vcc isnt 0 or 4", 1, states[xv][yv][zv]);
        switch(states[xi][yi][zi]){
		    case 2: // AA
			    nAA --; nA +=2;
	            states[xi][yi][zi]= 1;
	            states[xv][yv][zv]= 1;
			    break;
            case 3: // AB
			    nAB --; nA ++; nB ++;
	            states[xi][yi][zi]= 1;
	            states[xv][yv][zv]=-1;
                double e1= ecal_bond(xi, yi, zi, xv, yv, zv) + ecal_nonb(xi, yi, zi, xv, yv, zv);
	            states[xi][yi][zi]=-1;
	            states[xv][yv][zv]= 1;
                double e2= ecal_bond(xi, yi, zi, xv, yv, zv) + ecal_nonb(xi, yi, zi, xv, yv, zv);
                if(e1<e2){
	                states[xi][yi][zi]= 1;
	                states[xv][yv][zv]=-1;
                }
                break;
		    default: error(2, "(rules_recb) an unknown itl type", 1, states[xi][yi][zi]);
	    }
			
	    list_itl.erase(list_itl.begin()+ii);
    
        if(4==states[xv][yv][zv]) nM --;
        else{
            nV --;
            list_vcc.erase(list_vcc.begin()+iv);
        }

        if(nM !=0){
            srf_check(xi, yi, zi);
            srf_check(xv, yv, zi);
        }
    }
}

bool class_events::recb_checki(int id){
    vector<vector<int>> list_recb;
    vector<vector<int>> list_AAtoAB;
    
    int ltcp= list_itl[id].ltcp;
    int i= (int) (ltcp/nz)/ny;
	int j= (int) (ltcp/nz)%ny;
	int k= (int)  ltcp%nz;
    int stateI= states[i][j][k];
    int x, y, z;

	for(int a=0; a<n1nbr; a ++){
        x= pbc(i+v1nbr[a][0], nx);
		y= pbc(j+v1nbr[a][1], ny);
		z= pbc(k+v1nbr[a][2], nz);
        if(0==states[x][y][z] || 4==states[x][y][z]) list_recb.push_back({x, y, z});
        if(-1==states[x][y][z] && 2==stateI)         list_AAtoAB.push_back({x, y, z}); // AA+B->A+AB
    }
    if(list_recb.size() !=0) goto letsRECB;
	
    for(int a=0; a<n2nbr; a ++){
        x= pbc(i+v2nbr[a][0], nx);
		y= pbc(j+v2nbr[a][1], ny);
		z= pbc(k+v2nbr[a][2], nz);
        if(0==states[x][y][z] || 4==states[x][y][z]) list_recb.push_back({x, y, z});
        if(-1==states[x][y][z] && 2==stateI)         list_AAtoAB.push_back({x, y, z}); // AA+B->A+AB
    }
    if(list_recb.size() !=0) goto letsRECB;
	
    for(int a=0; a<n3nbr; a ++){
        x= pbc(i+v3nbr[a][0], nx);
		y= pbc(j+v3nbr[a][1], ny);
		z= pbc(k+v3nbr[a][2], nz);
        if(0==states[x][y][z] || 4==states[x][y][z]) list_recb.push_back({x, y, z});
    }
    if(list_recb.size() !=0) goto letsRECB;

    if(list_AAtoAB.size() !=0){ // AA+B->A+AB
        int ran= (int) ran_generator()*list_AAtoAB.size();
        x= list_AAtoAB[ran].at(0);
        y= list_AAtoAB[ran].at(1);
        z= list_AAtoAB[ran].at(2);
        rules_recb(id, i, j, k, -1, x, y, z); // vccID is unknown, give -1
        return true;
    }
    else return false;

letsRECB: // recb
    int ran= (int) ran_generator()*list_recb.size();
    x= list_recb[ran].at(0);
    y= list_recb[ran].at(1);
    z= list_recb[ran].at(2);
    rules_recb(id, i, j, k, -1, x, y, z); // vccID is unknown, give -1
    if(3==stateI) recb_checkv(x*ny*nz + y*nz + z, true);
    return true;
}

bool class_events::recb_checkv(int id, bool isB){
    vector<vector<int>> list_recb;
    vector<vector<int>> list_AAtoAB;
    
    int ltcp;
    if(isB) ltcp= id;
    else    ltcp= list_vcc[id].ltcp;
    int i= (int) (ltcp/nz)/ny;
	int j= (int) (ltcp/nz)%ny;
	int k= (int)  ltcp%nz;
    if(isB && states[i][j][k] != -1) error(2, "(recb_checkv) isB but not -1", 1, states[i][j][k]);
    int x, y, z;

	for(int a=0; a<n1nbr; a ++){
        x= pbc(i+v1nbr[a][0], nx);
		y= pbc(j+v1nbr[a][1], ny);
		z= pbc(k+v1nbr[a][2], nz);
        if(! isB && (2==states[x][y][z] || 3==states[x][y][z])) list_recb.push_back({x, y, z});
        if(isB && 2==states[x][y][z])                           list_AAtoAB.push_back({x, y, z}); // AA+B->A+AB
    }
    if(list_recb.size() !=0) goto letsRECB;
	
    for(int a=0; a<n2nbr; a ++){
        x= pbc(i+v2nbr[a][0], nx);
		y= pbc(j+v2nbr[a][1], ny);
		z= pbc(k+v2nbr[a][2], nz);
        if(! isB && (2==states[x][y][z] || 3==states[x][y][z])) list_recb.push_back({x, y, z});
        if(isB && 2==states[x][y][z])                           list_AAtoAB.push_back({x, y, z}); // AA+B->A+AB
    }
    if(list_recb.size() !=0) goto letsRECB;
	
    for(int a=0; a<n3nbr; a ++){
        x= pbc(i+v3nbr[a][0], nx);
		y= pbc(j+v3nbr[a][1], ny);
		z= pbc(k+v3nbr[a][2], nz);
        if(! isB && (2==states[x][y][z] || 3==states[x][y][z])) list_recb.push_back({x, y, z});
    }
    if(list_recb.size() !=0) goto letsRECB;
   
    if(list_AAtoAB.size() !=0){ // AA+B->A+AB
        int ran= (int) ran_generator()*list_AAtoAB.size();
        x= list_AAtoAB[ran].at(0);
        y= list_AAtoAB[ran].at(1);
        z= list_AAtoAB[ran].at(2);
        rules_recb(-1, x, y, z, -1, i, j, k); // AA+B->A+AB
        return true;
    }
    else return false;

letsRECB:
    int ran= (int) ran_generator()*list_recb.size();
    x= list_recb[ran].at(0);
    y= list_recb[ran].at(1);
    z= list_recb[ran].at(2);
    rules_recb(-1, x, y, z, id, i, j, k); // itlID is unknown, give -1

    return true;
}

void class_events::srf_check(int i, int j, int k){ // when vacuum changed, check if srf array changes
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

        srf[x][y][z]= is_srf; 
    }
}

// THE functions that have been deleted. To find them please go to ABVI_fixSINK //
// void class_events::recb_dir(int index){
// bool class_events::cal_dis(int d1, int d2, int d3){
// void class_events::recb_randomV(int index){
// bool class_events::recb_randomI(int index){
