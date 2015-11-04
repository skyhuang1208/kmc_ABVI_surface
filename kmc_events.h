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
		
            is_inf= false;
            einf= 0;
        
            crates= init_ratesC();
        }
		
        // variables
        double crates;
		
        // functions
		double main();
		double ecal_whole() const; 
	
	private:
        // infinite event //
        bool is_inf;
        double einf;
        vector <int> list_inf;

        // vcc creation at surface //
        #define RATIO_NOCVCC 0.5
        struct cvcc_info{
            long long int step; // modified step
            vector <int> mltcp;
            vector <double> rates;
        };
        unordered_map <int, struct cvcc_info> cvcc;
        double init_ratesC();

		///// functions of energy calculation /////
		double cal_energy(bool is_itl, int x1, int y1, int z1, int x2, int y2, int z2) const; 
		int powc(int base, int index) const;
		double cal_c00(int type[], int ABA1, int ABB1, int ABA2, int ABB2) const;
		double ecal_bond(bool is_itl, int x1, int y1, int z1, int x2, int y2, int z2) const; 
		double ecal_otf (bool is_itl, int x1, int y1, int z1, int x2, int y2, int z2) const; // corrected H on the fly 

		///// functions of events /////
		void actual_jumpI(int iid, int nltcp, int jatom);
        void actual_jumpV(int vid, int nltcp, int jatom);
		void genr();
        void create_vcc  (int altcp, int mltcp);
		
		///// functions for generation /////
        void genr_1strecb(int iid, int vid);
        
        ///// functions for rate calculations) /////
		double cal_ratesI(vector <int> &etype, vector <double> &rates, vector <int> &ilist, vector <int> &nltcp, vector <int> &jatom);
		double cal_ratesV(vector <int> &etype, vector <double> &rates, vector <int> &ilist, vector <int> &nltcp, vector <int> &jatom);
        double update_ratesC(int ltcp_in); // update rates of vcc creation at srf
        double cal_ratesC(vector <int> &etype, vector <double> &rates, vector <int> &ilist, vector <int> &nltcp, vector <int> &jatom);
		
        ///// functions for recombination /////
        void rules_recb(bool is_vcm, int ii, int iv, int ja= 0);
        void srf_check(int vltcp);
};

#endif // KMC_EVENTS_INCLUDED
