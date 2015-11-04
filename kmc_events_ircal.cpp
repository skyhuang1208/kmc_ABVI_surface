#include <iostream>
#include <cmath>
#include <vector>
#include "kmc_global.h"
#include "kmc_events.h"

using namespace std;

double class_events::cal_ratesI(vector <int> &etype, vector <double> &rates, vector <int> &ilist, vector <int> &nltcp, vector <int> &jatom){
	double sum_rate= 0;
	if(nAA + nAB + nBB != list_itl.size()) error(2, "(cal_ratesI) itl number inconsistent");
	
	for(int ii=0; ii < list_itl.size(); ii ++){ // ii: index of interstitial
		int ltcp= list_itl[ii].ltcp;
		int i= (int) (ltcp/nz)/ny;
		int j= (int) (ltcp/nz)%ny;
		int k= (int)  ltcp%nz;
		const int stateI= states[i][j][k]; // state of the itl

		if(1==abs(states[i][j][k])) error(2, "(itl_jump) there's an non-interstitial in the itl list");
	
        int ja;
        if(2==stateI)   ja= 1;
        else            ja=-1; // for 0==stateI try -1 first and then 1

        do{
            vector< vector <int> > rlist(3); // collect xyz to recover marker
            
            double em, mu;
//		    if(1==ja) { em= emiA; mu= muiA;} // WARNING! the dir move (AB atoms move in 2 ways) and rot move removed
//		    else	  { em= emiB; mu= muiB;}
		    
            if(1==ja) mu= muiA; // !!! only for Dubey, CMS 2015 !!!
            else	  mu= muiB; // !!!
            switch(stateI){     // !!!
                case  2: em= emiAA; break;
                case  0: em= emiAB; break;
                case -2: em= emiBB; break;
                default: error(2, "(cal_ratesI) an unknown itl type", 1, stateI);
            }                   // !!!

		    for(int a=0; a<n1nbr; a ++){
			    int x= pbc(i+v1nbr[a][0], nx);
			    int y= pbc(j+v1nbr[a][1], ny);
			    int z= pbc(k+v1nbr[a][2], nz);
			    const int stateA= states[x][y][z]; // state of the atom
		
			    if(stateA != 1 && stateA != -1) continue; // altho recb jump checks 1st-nn of 1st-nn, but only possible if I next to A
                
                // check if it's a recb jump
                bool is_recb= false;
                for(int b=0; b<n1nbr; b ++){
                    int xb= pbc(x+v1nbr[b][0], nx);
		            int yb= pbc(y+v1nbr[b][1], ny);
		            int zb= pbc(z+v1nbr[b][2], nz);
                    const int stateB= states[xb][yb][zb];

                    if((stateB != 0 && stateB != 4) || itlAB[xb][yb][zb]) continue; // if it's not vcc or vacuum, don't do 
                    is_recb= true;
                    if(marker[xb][yb][zb]) continue; // if it's already counted, don't do (but is_recb should be turned true)
                    
                    marker[xb][yb][zb]= true; // mark this ltcp so won't do it again
                    rlist[0].push_back(xb);   // collect xyz for marker recovery 
                    rlist[1].push_back(yb);  
                    rlist[2].push_back(zb); 
                        
		            // calculate energy diff
		            double e0= cal_energy(true, i, j, k, xb, yb, zb);
                    
                    states[i][j][k]    -= ja;
                    states[xb][yb][zb]  = ja;
                    if(0==stateI) itlAB[i][j][k]= false;
				
                    double ediff= cal_energy(true, i, j, k, xb, yb, zb) - e0;
				
                    states[i][j][k]     = stateI; //transit back
		            states[xb][yb][zb]  = stateB;
                    if(0==stateI) itlAB[i][j][k]= true;
			
                    if((em+0.5*ediff)<0){
                        double e= em+0.5*ediff;
                        rates.push_back(e);
                        is_inf= true;
                        list_inf.push_back(rates.size()-1);
                        if(e<einf) einf= e;
                    }
                    else if(0==stateI)   rates.push_back(0.5 * mu * exp(-beta*(em+0.5*ediff))); // itlAB has 2 jumps(A, B), and hence divided by 2 here 
                    else                 rates.push_back(      mu * exp(-beta*(em+0.5*ediff)));

       			    etype.push_back(0);
   				    ilist.push_back(ii);
    			    nltcp.push_back(xb*ny*nz+yb*nz+zb);
    			    jatom.push_back(ja);
    				
                    sum_rate += rates.back();
                }

                if(is_recb)
                    marker[x][y][z]= true; // mark this ltcp. don't do in cal_vrates
                else{
		            // calculate energy diff
		            double e0= cal_energy(true, i, j, k, x, y, z);
                   
                    states[i][j][k] -= ja; // perform imag jump
                    states[x][y][z] += ja;
                    if(0==stateI)          itlAB[i][j][k]= false;
                    if(0==states[x][y][z]) itlAB[x][y][z]= true;
			
                    double ediff= cal_energy(true, i, j, k, x, y, z) - e0;

                    double ec;                      // !!! only for Dubey, CMS 2015 !!!
                    switch(states[x][y][z]-stateI){ // !!!
                        case  0: ec= 0; break;
                        case  2:
                                 if(0==stateI) ec= eciABtAA;
                                 else          ec= eciBBtAB;
                                 break;
                        case -2:
                                 if(0==stateI) ec= eciABtBB;
                                 else          ec= eciAAtAB;
                                 break;
                        default: error(2, "(cal_ratesI) an unknown conversion type: (diff)", 1, states[x][y][z]-stateI);
                    }                               // !!!
			
                    states[i][j][k]  = stateI; //transit back
		            states[x][y][z]  = stateA;
                    if(0==stateI)   itlAB[i][j][k]= true;
                                    itlAB[x][y][z]= false;
			
                    if((ec+em+0.5*ediff)<0){
                        double e= ec+em+0.5*ediff;
                        rates.push_back(e);
                        is_inf= true;
                        list_inf.push_back(rates.size()-1);
                        if(e<einf) einf= e;
                    }
                    else if(0==stateI) rates.push_back(0.5 * mu * exp(-beta*(ec+em+0.5*ediff))); // itlAB has 2 jumps: via A or B. so divided by 2 here 
                    else               rates.push_back(      mu * exp(-beta*(ec+em+0.5*ediff)));
					
    			    etype.push_back(0);
   				    ilist.push_back(ii);
    			    nltcp.push_back(x*ny*nz+y*nz+z);
    			    jatom.push_back(ja);
   				
                    sum_rate += rates.back();
                }
		    }

            for(int j= 0; j<rlist[0].size(); j ++)
                marker[rlist[0].at(j)][rlist[1].at(j)][rlist[2].at(j)]= false;
            
            ja += 2;
        }while(0==stateI && 1==ja);
    }

	return sum_rate;
}
