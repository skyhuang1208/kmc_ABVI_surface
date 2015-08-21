#ifndef KMC_EVENTS_INCLUDED
#define KMC_EVENTS_INCLUDED
#include <vector>
#include <iostream>
#include <cmath>

using namespace std;

class class_events{
	public:
		class_events(){
			cout << "##Generation parameters: (dpa/s) " << 1.0/time_genr/nx/ny/nz << ", (time period)" << time_genr << endl;
			cout << "##Recombination parameters: 2nd nearest-neighbor distance (FIXED in SURFACE simulations) " << endl;
		}
		
		// functions
		void genr();
		double jump();
		double ecal_whole() const; 
	
	private:
		////// functions of energy calculation //////
		double cal_energy(bool is_itl, int x1, int y1, int z1, int x2, int y2, int z2) const; 
		int powc(int base, int index) const;
		double cal_c00(int type[], int ABA1, int ABB1, int ABA2, int ABB2) const;
		double ecal_bond(bool is_itl, int x1, int y1, int z1, int x2, int y2, int z2) const; 
		double ecal_otf (bool is_itl, int x1, int y1, int z1, int x2, int y2, int z2) const; // corrected H on the fly 

		////// functions for jumps(jumping rate calculations) //////
		double cal_ratesV(vector <bool> &isvcc, vector <double> &rates, vector <int> &ilist, vector <int> &nltcp, vector <int> &jatom);
		double cal_ratesI(vector <bool> &isvcc, vector <double> &rates, vector <int> &ilist, vector <int> &nltcp, vector <int> &jatom);
        void actual_jumpV(int vid, int nltcp, int jatom);
		void actual_jumpV(int vid, int nltcp);
		void actual_jumpI(int iid, int nltcp, int jatom);
	
		///// functions for recombination /////
        void genr_1strecb(int iid, int vid);
		
        ///// functions for recombination /////
        void rules_recb(bool is_vcm, int ii, int iv, int ja= 0);
        void srf_check(int vltcp);
};

#endif // KMC_EVENTS_INCLUDED
