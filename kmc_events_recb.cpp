#include <iostream>
#include <cmath>
#include <vector>
#include "kmc_global.h"
#include "kmc_events.h"
#include "kmc_par.h"

using namespace std;

void class_events::rules_recb(int ii, int xi, int yi, int zi, int iv, int xv, int yv, int zv){ // execute the recombination: vcc or vacuum
    if(-1==ii){ // if ii not know, search it
        for(int a=0; a<list_itl.size(); a ++) 
            if(list_itl[a].x==xi && list_itl[a].y==yi && list_itl[a].z==zi) ii= a;
        if(-1==ii) error(2, "(rules_recb) cant find the itl");
    }
    if(-1==iv && 0==states[xv][yv][zv]){
        for(int a=0; a<list_vcc.size(); a ++)
            if(list_vcc[a].x==xv && list_vcc[a].y==yv && list_vcc[a].z==zv) iv= a;
        if(-1==iv) error(2, "(rules_recb) cant find the vcc");
    }

    if(-1==states[xv][yv][zv]){ // AA+B->A+AB
        if(states[xi][yi][zi] != 2) error(2, "(rules_recb) a recb with B atom not AA itl (type)", 1, states[xi][yi][zi]);
        nAA --; nB --;
        states[xi][yi][zi]= 1;
        states[xv][yv][zv]= 3;
        nAB ++; nA ++;
        list_itl[ii].x= xv;
        list_itl[ii].y= yv;
        list_itl[ii].z= zv;
        recb_checki(ii);
    }
    else{ // recb: I+V or I+M
        if(4==states[xv][yv][zv]) nM --; 
        else if(0==states[xv][yv][zv]){
            nV --;
            list_vcc.erase(list_vcc.begin()+iv);
        }
        else error(2, "(rules_recb) the vcc isnt 0 or 4", 1, states[xv][yv][zv]);

        double e1, e2;
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
                e1= ecal_bond(xi, yi, zi, xv, yv, zv) + ecal_nonb(xi, yi, zi, xv, yv, zv);
	            states[xi][yi][zi]=-1;
	            states[xv][yv][zv]= 1;
                e2= ecal_bond(xi, yi, zi, xv, yv, zv) + ecal_nonb(xi, yi, zi, xv, yv, zv);
                if(e1<e2){ // comparing energy to decide whether A or B jumps to V
	                states[xi][yi][zi]= 1;
	                states[xv][yv][zv]=-1;
                }
                break;
		    default: error(2, "(rules_recb) an unknown itl type", 1, states[xi][yi][zi]);
	    }
			
	    list_itl.erase(list_itl.begin()+ii);
    
        if(nM !=0){
            srf_check(xi, yi, zi);
            srf_check(xv, yv, zv);
            cvcc_rates += update_ratesC(xi*ny*nz+yi*nz+zi, true);
            cvcc_rates += update_ratesC(xv*ny*nz+yv*nz+zv, true);
        }
    }
}

bool class_events::recb_checki(int id){
    int i= list_itl[id].x;
    int j= list_itl[id].y;
    int k= list_itl[id].z;
    if(i==x_sink){ // check if at sink
        sink(false, id);
        return true;
    }
    int stateI= states[i][j][k];
    
    vector<vector<int>> list_recb;
    vector<vector<int>> list_AAtoAB;
    double minE= 999999999; // if multiple cases, choose a minE case
    int x, y, z;
	for(int a=0; a<n123nbr; a ++){ // check recb
        x= pbc(i+v123nbr[a][0], nx);
		y= pbc(j+v123nbr[a][1], ny);
		z= pbc(k+v123nbr[a][2], nz);
        if(0==states[x][y][z] || 4==states[x][y][z]){
            int stateV= states[x][y][z];
            recb_check_ecal(true, minE, list_recb, i, j, k, stateI, x, y, z, stateV);
            // minE, list_recb are references, change in the func

            states[i][j][k]= stateI;
            states[x][y][z]= stateV;
        }
        if(a<n1nbr && ( -1==states[x][y][z] && 2==stateI )) list_AAtoAB.push_back({x, y, z}); // AA+B->A+AB
    }

    if(list_recb.size() !=0){ //recb
        double sro0;
        if(iscaldsro) sro0= cal_sro();
        
        int ran= (int) ( ran_generator()*list_recb.size() );
        x= list_recb[ran].at(0);
        y= list_recb[ran].at(1);
        z= list_recb[ran].at(2);
        rules_recb(id, i, j, k, -1, x, y, z); // vccID is unknown, give -1
        
        if(iscaldsro) acc_dsroRi += cal_sro()-sro0; 
        return true;
    }
    else if(list_AAtoAB.size() !=0){ // AA+B->A+AB
        double sro0;
        if(iscaldsro) sro0= cal_sro();
        
        int ran= (int) ( ran_generator()*list_AAtoAB.size() );
        x= list_AAtoAB[ran].at(0);
        y= list_AAtoAB[ran].at(1);
        z= list_AAtoAB[ran].at(2);
        rules_recb(id, i, j, k, -1, x, y, z); // vccID is unknown, give -1
        
        if(iscaldsro) acc_dsroRi += cal_sro()-sro0; 
        return true;
    }
    else return false;
}

bool class_events::recb_checkv(int id){
    int i= list_vcc[id].x;
    int j= list_vcc[id].y;
    int k= list_vcc[id].z;
    if(i==x_sink){ // check if at sink
        sink(true, id);
        return true;
    }
    int stateV= states[i][j][k];
    
    vector<vector<int>> list_recb;
    double minE= 999999999;
    int x, y, z;
	for(int a=0; a<n123nbr; a ++){ // check recb
        x= pbc(i+v123nbr[a][0], nx);
		y= pbc(j+v123nbr[a][1], ny);
		z= pbc(k+v123nbr[a][2], nz);

        if((2==states[x][y][z] || 3==states[x][y][z])){
            int stateI= states[x][y][z];
            recb_check_ecal(false, minE, list_recb, x, y, z, stateI, i, j, k, stateV);
            // minE, list_recb are references, change in the func

            states[x][y][z]= stateI;
            states[i][j][k]= stateV;
        }
    }

    if(list_recb.size() !=0){ //recb
        double sro0;
        if(iscaldsro) sro0= cal_sro(); 
        
        int ran= (int) ( ran_generator()*list_recb.size() );
        x= list_recb[ran].at(0);
        y= list_recb[ran].at(1);
        z= list_recb[ran].at(2);
        rules_recb(-1, x, y, z, id, i, j, k); // itlID is unknown, give -1
        
        if(iscaldsro) acc_dsroRv += cal_sro()-sro0; 
        return true;
    }
    else return false;
}

void class_events::recb_check_ecal(bool isitl, double& minE, vector<vector<int>>& list_recb, int xi, int yi, int zi, int stateI, int xv, int yv, int zv, int stateV){
    double e0, ediff;
    switch(stateI){ // perform image recb to cal ediff
        case 2:
            e0= ecal_bond(xi, yi, zi, xv, yv, zv);
            states[xi][yi][zi]= 1; 
            states[xv][yv][zv]= 1;
            ediff= ecal_bond(xi, yi, zi, xv, yv, zv) - e0;
            if(abs(ediff-minE)<1e-7){
                if(isitl)   list_recb.push_back({xv, yv, zv});
                else        list_recb.push_back({xi, yi, zi});
            }
            else if((ediff-minE)<0){
                list_recb.clear();
                if(isitl)   list_recb.push_back({xv, yv, zv});
                else        list_recb.push_back({xi, yi, zi});
                minE= ediff;
            }
                    
            break;
        case 3:
            e0= ecal_bond(xi, yi, zi, xv, yv, zv) + ecal_nonb(xi, yi, zi, xv, yv, zv);
            states[xi][yi][zi]=  1; 
            states[xv][yv][zv]= -1;
            ediff= ecal_bond(xi, yi, zi, xv, yv, zv) + ecal_nonb(xi, yi, zi, xv, yv, zv) - e0;
            if(abs(ediff-minE)<1e-7){
                if(isitl)   list_recb.push_back({xv, yv, zv});
                else        list_recb.push_back({xi, yi, zi});
            }
            else if((ediff-minE)<0){
                list_recb.clear();
                if(isitl)   list_recb.push_back({xv, yv, zv});
                else        list_recb.push_back({xi, yi, zi});
                minE= ediff;
            }
            
            states[xi][yi][zi]= -1; 
            states[xv][yv][zv]=  1;
            ediff= ecal_bond(xi, yi, zi, xv, yv, zv) + ecal_nonb(xi, yi, zi, xv, yv, zv) - e0;
            if(abs(ediff-minE)<1e-7){
                if(isitl)   list_recb.push_back({xv, yv, zv});
                else        list_recb.push_back({xi, yi, zi});
            }
            else if((ediff-minE)<0){
                list_recb.clear();
                if(isitl)   list_recb.push_back({xv, yv, zv});
                else        list_recb.push_back({xi, yi, zi});
                minE= ediff;
            }
                    
            break;
    }
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

void class_events::sink(bool isvcc, int index){ // execute the sink
	int x, y, z;
	double ran;
	
	if(isvcc){
		nV --;
		x= list_vcc[index].x;
		y= list_vcc[index].y;
		z= list_vcc[index].z;
		list_vcc.erase(list_vcc.begin()+index);

        ran= ran_generator();
		if(list_sink.size() != 0){ // if atoms in sink, use them
            int i= (int) ( ran*list_sink.size() );
			states[x][y][z]= list_sink[i];
            if(list_sink[i] != 1 && list_sink[i] != -1) error(2, "a atom from sink isnt atom", 1, list_sink[i]);
			list_sink.erase(list_sink.begin()+i);
		}
		else states[x][y][z]= (ran<par_compA) ? 1:-1;

		if(1==states[x][y][z])  nA ++;
		else                    nB ++;
	}
	else{
        x= list_itl[index].x;
        y= list_itl[index].y;
        z= list_itl[index].z;
		list_itl.erase(list_itl.begin()+index);
		
		switch(states[x][y][z]){
			case  2:
				nAA --;
				states[x][y][z]= 1; nA ++;
				list_sink.push_back(1);
				break;
			case  3:
				nAB --;
				states[x][y][z]= (ran_generator()>0.5) ? 1:-1;
                if(1==states[x][y][z]){ nA ++; list_sink.push_back(-1);}
                else                  { nB ++; list_sink.push_back( 1);}
				break;
			default: error(2, "(sink) an unknown type", 1, states[x][y][z]);
		}
	}
}

bool class_events::trap_check(int i, int j, int k){ // check if AB itl trapped
    if(states[i][j][k] != 3) error(2, "(trap_check) input not AB itl", 1, states[i][j][k]);
    if(! trap_included)      error(2, "(trap_check) trap_included not turn on");
    if(temp>1000)            error(2, "(trap_check) T >1000k, should not trapped", 1, temp);

    for(int a= 0; a<n1nbr; a ++){
        int x= pbc(i+v1nbr[a][0], nx);
        int y= pbc(j+v1nbr[a][1], ny);
        int z= pbc(k+v1nbr[a][2], nz);
        
        if(3==states[x][y][z] || -1==states[x][y][z]) return true; // if near by B or AB, trapped
    }

    return false; // not trapped
}

// THE functions that have been deleted. To find them please go to ABVI_fixSINK //
// void class_events::recb_dir(int index){
// bool class_events::cal_dis(int d1, int d2, int d3){
// void class_events::recb_randomV(int index){
// bool class_events::recb_randomI(int index){
