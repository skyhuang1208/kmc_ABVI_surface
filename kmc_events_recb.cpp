#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "kmc_global.h"
#include "kmc_events.h"
#include "kmc_par.h"

using namespace std;

void class_events::rules_recb(int ii, int xi, int yi, int zi, int iv, int xv, int yv, int zv){ // execute the recombination: vcc or vacuum
    int iltcp= xi*ny*nz + yi*nz + zi;
    int vltcp= xv*ny*nz + yv*nz + zv;
    if(-1==ii){
        for(int a=0; a<list_itl.size(); a ++) 
            if(list_itl[a].ltcp==iltcp) ii= a;
        if(-1==ii) error(2, "(rules_recb) cant find the itl");
    }
    if(-1==iv && (0==states[xv][yv][zv] || 5==states[xv][yv][zv]) ){
        for(int a=0; a<list_vcc.size(); a ++)
            if(list_vcc[a].ltcp==vltcp) iv= a;
        if(-1==iv) error(2, "(rules_recb) cant find the vcc");
    }

    int i; // ivoid
    switch (states[xv][yv][zv]){
        case 5:
            if(list_vcc[iv].ivoid==-1) error(2, "(rules_recb) ivoid=-1"); // check
            i= list_vcc[iv].ivoid;
            list_void[i].erase( find(list_void[i].begin(), list_void[i].end(), iv) );

            for(int a=list_void[i].size()-1; a>=0; a--){ // check if a void vcc disconnected
                int ivcc= list_void[i].at(a); // reverse order of a to void erase problem
                int n= 0;
                for(int b=0; b<n1nbr; b++){
                    int x= pbc( ( (int) (list_vcc[ivcc].ltcp/nz)/ny ) + v1nbr[b][0], nx);
	                int y= pbc( ( (int) (list_vcc[ivcc].ltcp/nz)%ny ) + v1nbr[b][1], ny);
	                int z= pbc( ( (int)  list_vcc[ivcc].ltcp%nz     ) + v1nbr[b][2], nz);

                    if(5==states[x][y][z] && x!=xv && y!=yv && z!=zv) n ++;
                } 
                if(n==0){
                    *(&states[0][0][0]+list_vcc[ivcc].ltcp)= 0;
                    list_vcc[ivcc].ivoid= -1;
                    list_void[i].erase(list_void[i].begin()+a);
                }
            }
            
            if(list_void[i].size()<=2){ // the void disappear
                for(int a=0; a<list_void[i].size(); a ++){
                    int iii= list_void[i].at(a); // get vid at the list
                    *(&states[0][0][0]+list_vcc[iii].ltcp)= 0;
                    list_vcc[iii].ivoid= -1;
                }
                list_void[i].clear();
                nVD --;
            } // no need break here

        case 0:
            nV --;
            
            erase_vcc(iv); // when erase in list_vcc, list_void needs correction (-- if >iv)
            break;
        case 4: 
            nM --; break;
        default: error(2, "(rules_recb) the vcc isnt 0 or 4 or 5", 1, states[xv][yv][zv]);
    }

    double e1, e2;
    nIrecb ++;
	nAA --; nA +=2;
	states[xi][yi][zi]= 1;
	states[xv][yv][zv]= 1;
			
	list_itl.erase(list_itl.begin()+ii);
}

bool class_events::rules_attach(int ii, int xi, int yi, int zi, int iv, int xv, int yv, int zv, bool attached){ // execute the recombination: vcc or vacuum
    if(states[xi][yi][zi] != 5) error(2, "(rules_attach) not a void", 1, states[xi][yi][zi]); // check
    if(states[xv][yv][zv] != 0) error(2, "(rules_attach) not a vcc",  1, states[xv][yv][zv]);

    int ltcp= xi*ny*nz + yi*nz + zi; // here i is not itl, but void
    for(int a=0; a<list_vcc.size(); a ++) 
        if(list_vcc[a].ltcp==ltcp) ii= a;
    if(-1==ii) error(2, "(rules_attach) cant find the void");

    if(list_vcc[ii].ivoid ==-1) error(2, "(rules_attach) ivoid==-1"); // check
    if(list_vcc[iv].ivoid !=-1) error(2, "(rules_attach) ivoid!=-1");

    if(attached) states[xv][yv][zv]= 5; // the vcc is right next to void
    else{
        int bkup=-1; // if cant find it after n=501; use backup value
	    int i, j, k;
        for(int n=1; n<=501; n ++){ // choose ltcp to attach 
            int ran= (int) ( ran_generator()*n1nbr );
            i= pbc(xi+v1nbr[ran][0], nx); // random 1st-nn of void
		    j= pbc(yi+v1nbr[ran][1], ny);
		    k= pbc(zi+v1nbr[ran][2], nz);

            if(n==501){ // if cant find any A that surrounded iv
                if(-1==bkup){ // no A atom surrounded ii (?!)
                    cout << " *** ATTACH FAIL *** ";
                    return false;
                }
                else{
                    cout << " *** FORCE ATTACH *** ";
                    i= pbc(xi+v1nbr[bkup][0], nx);
		            j= pbc(yi+v1nbr[bkup][1], ny);
		            k= pbc(zi+v1nbr[bkup][2], nz);
                    break;
                }
            }
            
            if(states[i][j][k] != 1) continue;
            bkup= ran;

            bool isFOUND= false;
	        for(int a=0; a<n1nbr; a ++){ // check if ltcp 1st-nn of vcc
                int x= pbc(xv+v1nbr[a][0], nx);
		        int y= pbc(yv+v1nbr[a][1], ny);
		        int z= pbc(zv+v1nbr[a][2], nz);
                if(i==x && j==y && k==z) isFOUND= true;
            }
            if(isFOUND) break;
        }
        
        states[xv][yv][zv]= 1; // switch position
        states[i][j][k]=    5;
        list_vcc[iv].ltcp= i*ny*nz + j*nz + k;
        xv= i; yv= j; zv= k;
    }

    int ivd= list_vcc[ii].ivoid;
    list_void[ivd].push_back(iv);
    list_vcc[iv].ivoid= ivd;
    
    int check= *(&states[0][0][0] + list_vcc[iv].ltcp);
    if(check != 5) error(2, "(rules_attach) the attached vcc isnt 5", 1, check);

    bool isRECB= recb_checkv(iv);
    if(! isRECB) coll_check(ivd, xv, yv, zv); // collision check
    return true;
}

void class_events::coll_check(int ivd, int i, int j, int k){ // collision check
    vector <int> list_coll; // list of ID of voids that mixed together
    for(int a=0; a<n1nbr; a++){
        int x= pbc(i+v1nbr[a][0], nx);
		int y= pbc(j+v1nbr[a][1], ny);
		int z= pbc(k+v1nbr[a][2], nz);
        if(5==states[x][y][z]){
            int iv;
            int ltcp= x*ny*nz + y*nz + z;
            for(int a=0; a<list_vcc.size(); a ++) 
                if(list_vcc[a].ltcp==ltcp && list_vcc[a].ivoid != ivd) list_coll.push_back(list_vcc[a].ivoid);
        }
    }

    for(int a=0; a<list_coll.size(); a++){ // do collision, mix voids
        int ivd2= list_coll[a]; // ivd2 will mix to ivd
        if(-1==ivd2) error(2, "(coll_check) found a -1 in ID of list_void");
        if(list_void[ivd2].size()!=0) nVD --;
        for(int b=0; b<list_void[ivd2].size(); b++){
            int iv= list_void[ivd2].at(b); // ID of list_vcc
            list_void[ivd].push_back(iv);
            list_vcc[iv].ivoid= ivd;
        }
        list_void[ivd2].clear();
    }
}

bool class_events::recb_checki(int id){
    int ltcp= list_itl[id].ltcp;
    int i= (int) (ltcp/nz)/ny;
	int j= (int) (ltcp/nz)%ny;
	int k= (int)  ltcp%nz;
    if(i==x_sink || j==y_sink || k==z_sink){ // check if at sink
        sink(false, id);
        return true;
    }
   
    vector<vector<int>> list_recb;
    int x, y, z;
	for(int a=0; a<n123nbr; a ++){ // check recb
        x= pbc(i+v123nbr[a][0], nx);
		y= pbc(j+v123nbr[a][1], ny);
		z= pbc(k+v123nbr[a][2], nz);
        if(0==states[x][y][z] || 4==states[x][y][z] || 5==states[x][y][z]) list_recb.push_back({x, y, z});
    }

    if(list_recb.size() !=0){ //recb
        int ran= (int) ( ran_generator()*list_recb.size() );
        x= list_recb[ran].at(0);
        y= list_recb[ran].at(1);
        z= list_recb[ran].at(2);
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
    if(i==x_sink || j==y_sink || k==z_sink){ // check if at sink
        sink(true, id);
        return true;
    }
    
    vector<vector<int>> list_recb;
    vector<vector<int>> list_attach; // single vcc attached onto void
    int x, y, z;
    bool attached= false; // != -1 if vcc is 1st-nn close to void
	for(int a=0; a<n123nbr; a ++){ // check recb
        x= pbc(i+v123nbr[a][0], nx);
		y= pbc(j+v123nbr[a][1], ny);
		z= pbc(k+v123nbr[a][2], nz);

        if(2==states[x][y][z]) list_recb.push_back({x, y, z});
        if(5==states[x][y][z]) list_attach.push_back({x, y, z});
        if(5==states[x][y][z] && a<n1nbr) attached= true;
    }

    if(list_recb.size() !=0){ //recb
        int ran= (int) ( ran_generator()*list_recb.size() );
        x= list_recb[ran].at(0);
        y= list_recb[ran].at(1);
        z= list_recb[ran].at(2);
        rules_recb(-1, x, y, z, id, i, j, k); // itlID is unknown, give -1
        return true;
    }
    else if(0==states[i][j][k] && list_attach.size() != 0){ // attach
        int ia; // index of attach
        if(attached) ia= 0;
        else         ia= (int) ( ran_generator()*list_attach.size() );
        x= list_attach[ia].at(0);
        y= list_attach[ia].at(1);
        z= list_attach[ia].at(2);
        bool is_attached= rules_attach(-1, x, y, z, id, i, j, k, attached); // voidID is unknown, give -1
        return is_attached; // if fail to attach, return false
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

		if(5==*(&states[0][0][0]+ltcp)) cout << "  *** void at sink ***";
        erase_vcc(index); // when erase in list_vcc, list_void needs correction (-- if >iv)

        ran= ran_generator();
		*(&states[0][0][0]+ltcp)= 1;

		nA ++;
	}
	else{
        nIsink ++;
		ltcp= list_itl[index].ltcp;
		list_itl.erase(list_itl.begin()+index);
		
		nAA --;
		*(&states[0][0][0]+ltcp)= 1;
        nA ++;
	}
}

void class_events::erase_vcc(int iv){ // when erase iv, list_void needs correct
    for(int a= 0; a<list_void.size(); a++){
        for(int b= 0; b<list_void[a].size(); b++){
            if     (list_void.at(a).at(b)==iv) list_void[a].erase(list_void[a].begin()+b);
            else if(list_void.at(a).at(b) >iv) list_void[a].at(b) --; // iv deleted, minus 1
    }}

	list_vcc.erase(list_vcc.begin()+iv);
}




// THE functions that have been deleted. To find them please go to ABVI_fixSINK //
// void class_events::recb_dir(int index){
// bool class_events::cal_dis(int d1, int d2, int d3){
// void class_events::recb_randomV(int index){
// bool class_events::recb_randomI(int index){

