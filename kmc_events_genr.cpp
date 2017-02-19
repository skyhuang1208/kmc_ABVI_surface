#include <cstdio>
#include <iostream>
#include <vector>
#include "kmc_global.h"
#include "kmc_events.h"

using namespace std;

void class_events::genr(){
    // GENERATE VCC or VOID
    int vsize;
    double ran= ran_generator();
    for(int a=1; a<=1+n1nbr; a ++) if(ran > pVsize[a-1] && ran <= pVsize[a]) vsize= a;

    int x[vsize], y[vsize], z[vsize];
	for(int i=1; i<= 10000; i ++){
        bool isfound= true;
        x[0]= (int) (ran_generator()*nx); // central atom of the void
		y[0]= (int) (ran_generator()*ny);
		z[0]= (int) (ran_generator()*nz);

        bool picks[30]={}; // initial to false
        for(int j=1; j<vsize; j++){ // pick nbrs
            int ran;
            do{ ran= (int) (ran_generator()*n1nbr);} while(picks[ran]);
                
            picks[ran]= true;
            x[j]= pbc(x[0]+v1nbr[ran][0], nx); 
            y[j]= pbc(y[0]+v1nbr[ran][1], ny); 
            z[j]= pbc(z[0]+v1nbr[ran][2], nz);
        }
			
        for(int k=0; k<vsize; k++){ // check if void's on or next to other void
            if(states[x[k]][y[k]][z[k]] != 1)                   isfound= false;
            if(x[k]==x_sink || y[k]==y_sink || z[k]==z_sink)    isfound= false;

		    for(int a=0; a<n1nbr; a ++){ // check next
			    int x1= pbc(x[k]+v1nbr[a][0], nx);
			    int y1= pbc(y[k]+v1nbr[a][1], ny);
			    int z1= pbc(z[k]+v1nbr[a][2], nz);
                if(states[x1][y1][z1] != 1) isfound= false;
            }
        }

        if(isfound) break; // if found then exit
        if(i==10000) error(1, "(init_states_array) cant find a place to put void");
    }
            
    if(vsize>=3){ // found locations, now build void or vcc  
        nVD ++;
        list_void.push_back(vector<int>());
    }
    for(int n=0; n<vsize; n++){
        nV ++; nA --;

        list_vcc.push_back(vcc());
		list_vcc.back().ltcp= x[n]*ny*nz + y[n]*nz + z[n];

        if(vsize>=3){ // void
            states[x[n]][y[n]][z[n]]= 5;
			list_vcc.back().ivoid= list_void.size() -1;
            list_void.back().push_back(list_vcc.size()-1);
        }
        else{ // vcc
            states[x[n]][y[n]][z[n]]= 0;
            list_vcc.back().ivoid= -1;
        }
    }

    // GENERATE ITL
    for(int i=1; i<=vsize; i++){
        nAA ++; nA --;
        
        int xi, yi, zi;
        for(int n= 1; n<=5000; n ++){
            bool isfound= true;

            xi= (int) (ran_generator()*nx); // positions of itl
		    yi= (int) (ran_generator()*ny);
		    zi= (int) (ran_generator()*nz);
    
            if(states[xi][yi][zi] != 1)                 isfound= false;
            if(xi==x_sink || yi==y_sink || zi==z_sink)  isfound= false;

		    for(int a=0; a<n1nbr; a ++){ // check next
			    int x1= pbc(xi+v1nbr[a][0], nx);
			    int y1= pbc(yi+v1nbr[a][1], ny);
			    int z1= pbc(zi+v1nbr[a][2], nz);
                if(states[x1][y1][z1] != 1) isfound= false;
            }

            if(isfound) break;
            if(n==5000) error(2, "(genr) cant find a place to generate itl");
        }
        
	    list_itl.push_back(itl());
	    list_itl.back().ltcp= xi*ny*nz + yi*nz + zi; 
	    list_itl.back().dir= (int) (ran_generator()*n1nbr);
	    list_itl.back().head= 1;
        
	    states[xi][yi][zi]= 2;
    }
    
    // CHECK RECB
    for(int i=list_itl.size()-1; i>=0; i--) recb_checki(i);
	for(int v=0; v < nV; v ++){
        bool isRECB= recb_checkv(v);
        if(isRECB) v= 0;
    }
}
