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
                if(e1<e2){
	                states[xi][yi][zi]= 1;
	                states[xv][yv][zv]=-1;
                }
                break;
		    default: error(2, "(rules_recb) an unknown itl type", 1, states[xi][yi][zi]);
	    }
			
	    list_itl.erase(list_itl.begin()+ii);
    
        if(nM !=0){
            srf_check(xi, yi, zi);
            srf_check(xv, yv, zi);
            cvcc_rates += update_ratesC(xi*ny*nz+yi*nz+zi, true);
            cvcc_rates += update_ratesC(xv*ny*nz+yv*nz+zv, true);
        }
    }
}

bool class_events::recb_checki(int id){
    int ltcp= list_itl[id].ltcp;
    int i= (int) (ltcp/nz)/ny;
	int j= (int) (ltcp/nz)%ny;
	int k= (int)  ltcp%nz;
    if(i==x_sink){ // check if at sink
        sink(false, id);
        return true;
    }
    int stateI= states[i][j][k];
    
    vector<vector<int>> list_recb;
    vector<vector<int>> list_AAtoAB;
    double minE= 999999999;
    int x, y, z;
	for(int a=0; a<n123nbr; a ++){ // check recb
        x= pbc(i+v123nbr[a][0], nx);
		y= pbc(j+v123nbr[a][1], ny);
		z= pbc(k+v123nbr[a][2], nz);
        if(0==states[x][y][z] || 4==states[x][y][z]){
            int stateV= states[x][y][z];
            double e0, ediff;
            switch(stateI){
                case 2:
                    e0= ecal_bond(x, y, z, i, j, k);
                    states[i][j][k]= 1; states[x][y][z]= 1;
                    ediff= ecal_bond(x, y, z, i, j, k) - e0;
                    if(abs(ediff-minE)<1e-7) list_recb.push_back({x, y, z});
                    else if((ediff-minE)<0){
                        list_recb.clear(); list_recb.push_back({x, y, z});
                        minE= ediff;
                    }
                    break;
                case 3:
                    e0= ecal_bond(x, y, z, i, j, k) + ecal_nonb(x, y, z, i, j, k);
                    
                    states[i][j][k]=  1; states[x][y][z]= -1;
                    ediff= ecal_bond(x, y, z, i, j, k) + ecal_nonb(x, y, z, i, j, k) - e0;
                    if(abs(ediff-minE)<1e-7) list_recb.push_back({x, y, z});
                    else if((ediff-minE)<0){
                        list_recb.clear(); list_recb.push_back({x, y, z});
                        minE= ediff;
                    }
                    
                    states[i][j][k]= -1; states[x][y][z]=  1;
                    ediff= ecal_bond(x, y, z, i, j, k) + ecal_nonb(x, y, z, i, j, k) - e0;
                    if(abs(ediff-minE)<1e-7) list_recb.push_back({x, y, z});
                    else if((ediff-minE)<0){
                        list_recb.clear(); list_recb.push_back({x, y, z});
                        minE= ediff;
                    }
                    break;
            }

            states[i][j][k]= stateI;
            states[x][y][z]= stateV;
        }
        if(a<n1nbr && ( -1==states[x][y][z] && 2==stateI )) list_AAtoAB.push_back({x, y, z}); // AA+B->A+AB
    }

    if(list_recb.size() !=0){ //recb
        int ran= (int) ran_generator()*list_recb.size();
        x= list_recb[ran].at(0);
        y= list_recb[ran].at(1);
        z= list_recb[ran].at(2);
        rules_recb(id, i, j, k, -1, x, y, z); // vccID is unknown, give -1
        return true;
    }
    else if(list_AAtoAB.size() !=0){ // AA+B->A+AB
        int ran= (int) ran_generator()*list_AAtoAB.size();
        x= list_AAtoAB[ran].at(0);
        y= list_AAtoAB[ran].at(1);
        z= list_AAtoAB[ran].at(2);
        rules_recb(id, i, j, k, -1, x, y, z); // vccID is unknown, give -1
        return true;
    }
    else return false;
}

bool class_events::recb_checkv(int id){
    int ltcp= list_vcc[id].ltcp;
    int i= (int) (ltcp/nz)/ny;
	int j= (int) (ltcp/nz)%ny;
	int k= (int)  ltcp%nz;
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
            double e0, ediff;
            switch(stateI){
                case 2:
                    e0= ecal_bond(x, y, z, i, j, k);
                    states[i][j][k]= 1; states[x][y][z]= 1;
                    ediff= ecal_bond(x, y, z, i, j, k) - e0;
                    if(abs(ediff-minE)<1e-7) list_recb.push_back({x, y, z});
                    else if((ediff-minE)<0){
                        list_recb.clear(); list_recb.push_back({x, y, z});
                        minE= ediff;
                    }
                    break;
                case 3:
                    e0= ecal_bond(x, y, z, i, j, k) + ecal_nonb(x, y, z, i, j, k);
                    
                    states[i][j][k]=  1; states[x][y][z]= -1;
                    ediff= ecal_bond(x, y, z, i, j, k) + ecal_nonb(x, y, z, i, j, k) - e0;
                    if(abs(ediff-minE)<1e-7) list_recb.push_back({x, y, z});
                    else if((ediff-minE)<0){
                        list_recb.clear(); list_recb.push_back({x, y, z});
                        minE= ediff;
                    }
                    
                    states[i][j][k]= -1; states[x][y][z]=  1;
                    ediff= ecal_bond(x, y, z, i, j, k) + ecal_nonb(x, y, z, i, j, k) - e0;
                    if(abs(ediff-minE)<1e-7) list_recb.push_back({x, y, z});
                    else if((ediff-minE)<0){
                        list_recb.clear(); list_recb.push_back({x, y, z});
                        minE= ediff;
                    }
                    break;
            }

            states[x][y][z]= stateI;
            states[i][j][k]= stateV;
        }
    }

    if(list_recb.size() !=0){ //recb
        int ran= (int) ran_generator()*list_recb.size();
        x= list_recb[ran].at(0);
        y= list_recb[ran].at(1);
        z= list_recb[ran].at(2);
        rules_recb(-1, x, y, z, id, i, j, k); // itlID is unknown, give -1
        return true;
    }
    else return false;
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
	int ltcp;
	double ran;
	
	if(isvcc){
		nV --;
		ltcp= list_vcc[index].ltcp;
		list_vcc.erase(list_vcc.begin()+index);

        ran= ran_generator();
		if(list_sink.size() != 0){
            int i= (int) ran*list_sink.size();
			*(&states[0][0][0]+ltcp)= list_sink[i];
            if(list_sink[i] != 1 && list_sink[i] != -1) error(2, "a atom from sink isnt atom", 1, list_sink[i]);
			list_sink.erase(list_sink.begin()+i);
		}
		else *(&states[0][0][0]+ltcp)= (ran<par_compA) ? 1:-1;

		if(1==*(&states[0][0][0]+ltcp)) nA ++;
		else                            nB ++;
	}
	else{
		ltcp= list_itl[index].ltcp;
		list_itl.erase(list_itl.begin()+index);
		
		switch(*(&states[0][0][0]+ltcp)){
			case  2:
				nAA --;
				list_sink.push_back(1);
				*(&states[0][0][0]+ltcp)= 1; nA ++;
				break;
			case  3:
				nAB --;
				ran= ran_generator();
				if(ran>0.5){
					list_sink.push_back(-1);
					*(&states[0][0][0]+ltcp)= 1; nA ++;
				}
				else{
					list_sink.push_back(1);
					*(&states[0][0][0]+ltcp)=-1; nB ++;
				}
				break;
			default: error(2, "(sink) an unknown type", 1, *(&states[0][0][0]+ltcp));
		}
	}
}
// THE functions that have been deleted. To find them please go to ABVI_fixSINK //
// void class_events::recb_dir(int index){
// bool class_events::cal_dis(int d1, int d2, int d3){
// void class_events::recb_randomV(int index){
// bool class_events::recb_randomI(int index){
