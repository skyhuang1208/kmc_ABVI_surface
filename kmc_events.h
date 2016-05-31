#ifndef KMC_EVENTS_INCLUDED
#define KMC_EVENTS_INCLUDED
#include <vector>
#include <iostream>
#include <cmath>
#include <unordered_map>

using namespace std;

class class_events{
	public:
		class_events(){
			cout << "##Generation parameters (rate_genr) " << rate_genr << " (damage/s)" << endl;
			cout << "##Recombination parameters: 2nd nearest-neighbor distance (FIXED in SURFACE simulations) " << endl;
		
            cvcc_rates= init_ratesC();

            // create v12nbr
			for(int a=0; a<n1nbr; a ++) // 1st neighbors
			    v12nbr.push_back({ (*(v1nbr+a))[0], (*(v1nbr+a))[1], (*(v1nbr+a))[2] });
			for(int a=0; a<n2nbr; a ++) // 2nd neighbors
			    v12nbr.push_back({ (*(v2nbr+a))[0], (*(v2nbr+a))[1], (*(v2nbr+a))[2] });
        }
		
        // variables
        double cvcc_rates;
		
        // functions
		double main();
        double ecal_range(int xlo=0, int xhi=nx-1, int ylo=0, int yhi=ny-1, int zlo=0, int zhi=nz-1, bool is_corr= false);
		double ecal_whole() const; // DELETE IT: NO CONC DEP
	
	private:
        // ********************  VARIABLES  ******************** //
        // vcc creation at surface //
        #define RATIO_NOCVCC 0.5
        struct cvcc_info{
            long long int step; // modified step
            vector <int> mltcp;
            vector <double> rates;
        };
        unordered_map <int, struct cvcc_info> cvcc;

        // 1st-nn + 2nd-nn //
        int n12nbr= n1nbr + n2nbr;
        vector < vector<int> > v12nbr;
        // ********************  VARIABLES  ******************** //
        
		///// functions of energy calculation /////
		double cal_energy(bool is_itl, int x1, int y1, int z1, int x2, int y2, int z2) const; 
		int powc(int base, int index) const;
		double cal_c00(int type[], int ABA1, int ABB1, int ABA2, int ABB2) const;
		double ecal_bond(int x1, int y1, int z1, int x2, int y2, int z2) const; 
        double ecal_nonb(int x1, int y1, int z1, int x2, int y2, int z2) const; // cal non-broken A-B bonds 
		double ecal_otf (bool is_itl, int x1, int y1, int z1, int x2, int y2, int z2) const; // corrected H on the fly 
        double ecal_sp(int stateA1, int inbr, int i, int j, int k) const; // calculate saddle-point e 

		///// functions of events /////
		void actual_jumpI(int iid, int jatom);
        void actual_jumpV(int vid, int jatom);
		void genr();
        void create_vcc  (int altcp, int mltcp);
		
		///// functions for generation /////
        void genr_1strecb(int iid, int vid);
        
        ///// functions for rate calculations) /////
		double cal_ratesI(vector <int> &etype, vector <double> &rates, vector <int> &ilist, vector <int> &jatom);
		double cal_ratesV(vector <int> &etype, vector <double> &rates, vector <int> &ilist, vector <int> &jatom);
        double cal_ratesVsp(vector <int> &etype, vector <double> &rates, vector <int> &ilist, vector <int> &jatom);
        double update_ratesC(int ltcp_in); // update rates of vcc creation at srf
        double init_ratesC();
		
        ///// functions for recombination /////
        void rules_recb(int ii, int xi, int yi, int zi, int iv, int xv, int yv, int zv);
        bool recb_checki(int id);
        bool recb_checkv(int id, bool isB= false);
        void srf_check(int i, int j, int k);
};

#endif // KMC_EVENTS_INCLUDED
        
        // cvcc straight jump info (no flickering) //
//        int NF_id; // the newly created vcc that the straight jumps is under going
//        int NF_Nj; // the number of jumps remained
//        double NF_rates;
//        vector <int> list_nf; // a list of ids in rates vector that is used for NF
//        vector <int> list_update;

