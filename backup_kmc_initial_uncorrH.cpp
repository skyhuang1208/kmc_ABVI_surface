#include <cstdlib>
#include <fstream>
#include <iostream>
#include "kmc_initial.h"
#include "kmc_par.h"

#define MAX_NNBR 20

using namespace std;

void class_initial::init_uncorrH(){
	// 1st-nn class 1
	unc1_44= ( (eAA1AA -8*eAA1A -8*eAA1B +12*eAA1AB +2*eAA1BB +12*eAB1BB -8*eA1BB -8*eB1BB +eBB1BB) +
	         (-48*eA1AB +36*eAB1AB -48*eAB1B) + (-48*eA1V +36*eV1V -48*eV1B) + (16*eA1A +32*eA1B +16*eB1B) )/576; 
	 
	unc1_42= ( (-eAA1AA +20*eAA1A +20*eAA1B -36*eAA1AB -2*eAA1BB -36*eAB1BB +20*eA1BB +20*eB1BB -eBB1BB) +
	         (216*eA1AB -180*eAB1AB +216*eAB1B) + (216*eA1V -180*eV1V +216*eV1B) + (-64*eA1A -128*eA1B -64*eB1B) )/576; 

	unc1_22= ( (eAA1AA -32*eAA1A -32*eAA1B +60*eAA1AB +2*eAA1BB +60*eAB1BB -32*eA1BB -32*eB1BB +eBB1BB) +
	         (-960*eA1AB +900*eAB1AB -960*eAB1B) + (-960*eA1V +900*eV1V -960*eV1B) + (256*eA1A +512*eA1B +256*eB1B) )/576; 
	
	// 1st-nn class 2
	unc1_43= ( (eAA1AA -6*eAA1A -2*eAA1B +6*eAA1AB -6*eAB1BB +2*eA1BB +6*eB1BB -eBB1BB) +
	         (-12*eA1AB +12*eAB1B) + (-12*eA1V +12*eV1B) + (8*eA1A -8*eB1B) )/288; 
	 
	unc1_41= ( (-eAA1AA +12*eAA1A -4*eAA1B -6*eAA1AB +6*eAB1BB +4*eA1BB -12*eB1BB +eBB1BB) +
	         (48*eA1AB -48*eAB1B) + (48*eA1V -48*eV1B) + (-32*eA1A +32*eB1B) )/288; 
	
	unc1_32= ( (-eAA1AA +18*eAA1A +14*eAA1B -30*eAA1AB +30*eAB1BB -14*eA1BB -18*eB1BB +eBB1BB) +
	         (60*eA1AB -60*eAB1B) + (60*eA1V -60*eV1B) + (-32*eA1A +32*eB1B) )/288; 
	
	unc1_21= ( (eAA1AA -24*eAA1A -8*eAA1B +30*eAA1AB -30*eAB1BB +8*eA1BB +24*eB1BB -eBB1BB) +
	         (-240*eA1AB +240*eAB1B) + (-240*eA1V +240*eV1B) + (128*eA1A -128*eB1B) )/288; 
	
	// 1st-nn class 3
	unc1_33= ( (eAA1AA -4*eAA1A +4*eAA1B -2*eAA1BB +4*eA1BB -4*eB1BB +eBB1BB) +
	         (4*eA1A -8*eA1B +4*eB1B) )/144; 
	
	unc1_31= ( (-eAA1AA +10*eAA1A -10*eAA1B +2*eAA1BB -10*eA1BB +10*eB1BB -eBB1BB) +
	         (-16*eA1A +32*eA1B -16*eB1B) )/144; 
	
	unc1_11= ( (eAA1AA -16*eAA1A +16*eAA1B -2*eAA1BB +16*eA1BB -16*eB1BB +eBB1BB) +
	         (64*eA1A -128*eA1B +64*eB1B) )/144; 
	
	// 2nd-nn class 1
	unc2_44= ( (eAA2AA -8*eAA2A -8*eAA2B +12*eAA2AB +2*eAA2BB +12*eAB2BB -8*eA2BB -8*eB2BB +eBB2BB) +
	         (-48*eA2AB +36*eAB2AB -48*eAB2B) + (-48*eA2V +36*eV2V -48*eV2B) + (16*eA2A +32*eA2B +16*eB2B) )/576; 
	 
	unc2_42= ( (-eAA2AA +20*eAA2A +20*eAA2B -36*eAA2AB -2*eAA2BB -36*eAB2BB +20*eA2BB +20*eB2BB -eBB2BB) +
	         (216*eA2AB -180*eAB2AB +216*eAB2B) + (216*eA2V -180*eV2V +216*eV2B) + (-64*eA2A -128*eA2B -64*eB2B) )/576; 

	unc2_22= ( (eAA2AA -32*eAA2A -32*eAA2B +60*eAA2AB +2*eAA2BB +60*eAB2BB -32*eA2BB -32*eB2BB +eBB2BB) +
	         (-960*eA2AB +900*eAB2AB -960*eAB2B) + (-960*eA2V +900*eV2V -960*eV2B) + (256*eA2A +512*eA2B +256*eB2B) )/576; 
	
	// 2nd-nn class 2
	unc2_43= ( (eAA2AA -6*eAA2A -2*eAA2B +6*eAA2AB -6*eAB2BB +2*eA2BB +6*eB2BB -eBB2BB) +
	         (-12*eA2AB +12*eAB2B) + (-12*eA2V +12*eV2B) + (8*eA2A -8*eB2B) )/288; 
	 
	unc2_41= ( (-eAA2AA +12*eAA2A -4*eAA2B -6*eAA2AB +6*eAB2BB +4*eA2BB -12*eB2BB +eBB2BB) +
	         (48*eA2AB -48*eAB2B) + (48*eA2V -48*eV2B) + (-32*eA2A +32*eB2B) )/288; 
	
	unc2_32= ( (-eAA2AA +18*eAA2A +14*eAA2B -30*eAA2AB +30*eAB2BB -14*eA2BB -18*eB2BB +eBB2BB) +
	         (60*eA2AB -60*eAB2B) + (60*eA2V -60*eV2B) + (-32*eA2A +32*eB2B) )/288; 
	
	unc2_21= ( (eAA2AA -24*eAA2A -8*eAA2B +30*eAA2AB -30*eAB2BB +8*eA2BB +24*eB2BB -eBB2BB) +
	         (-240*eA2AB +240*eAB2B) + (-240*eA2V +240*eV2B) + (128*eA2A -128*eB2B) )/288; 
	
	// 2nd-nn class 3
	unc2_33= ( (eAA2AA -4*eAA2A +4*eAA2B -2*eAA2BB +4*eA2BB -4*eB2BB +eBB2BB) +
	         (4*eA2A -8*eA2B +4*eB2B) )/144; 
	
	unc2_31= ( (-eAA2AA +10*eAA2A -10*eAA2B +2*eAA2BB -10*eA2BB +10*eB2BB -eBB2BB) +
	         (-16*eA2A +32*eA2B -16*eB2B) )/144; 
	
	unc2_11= ( (eAA2AA -16*eAA2A +16*eAA2B -2*eAA2BB +16*eA2BB -16*eB2BB +eBB2BB) +
	         (64*eA2A -128*eA2B +64*eB2B) )/144; 

	// 1st-nn constants (for itl jump only)
	unc1_40= ( ( eAA1AB +eAB1BB) + (-4*eA1AB  +6*eAB1AB  -4*eAB1B) + (-4*eA1V  +6*eV1V  -4*eV1B) )/24; 
	unc1_30= ( ( eAA1AB -eAB1BB) + (-2*eA1AB             +2*eAB1B) + (-2*eA1V           +2*eV1B) )/12; 
	unc1_20= ( (-eAA1AB -eAB1BB) + (16*eA1AB -30*eAB1AB +16*eAB1B) + (16*eA1V -30*eV1V +16*eV1B) )/24; 
	unc1_10= ( (-eAA1AB +eAB1BB) + ( 8*eA1AB             -8*eAB1B) + ( 8*eA1V           -8*eV1B) )/12; 
	unc1_00= eAB1AB + eV1V;

	// 2nd-nn constants (for itl jump only)
	unc2_40= ( ( eAA2AB +eAB2BB) + (-4*eA2AB  +6*eAB2AB  -4*eAB2B) + (-4*eA2V  +6*eV2V  -4*eV2B) )/24; 
	unc2_30= ( ( eAA2AB -eAB2BB) + (-2*eA2AB             +2*eAB2B) + (-2*eA2V           +2*eV2B) )/12; 
	unc2_20= ( (-eAA2AB -eAB2BB) + (16*eA2AB -30*eAB2AB +16*eAB2B) + (16*eA2V -30*eV2V +16*eV2B) )/24; 
	unc2_10= ( (-eAA2AB +eAB2BB) + ( 8*eA2AB             -8*eAB2B) + ( 8*eA2V           -8*eV2B) )/12; 
	unc2_00= eAB2AB + eV2V;

	// print out the parameters to log file
	cout << endl << "# uncorrected H: #" << endl;
	cout << "(1st neighbor)" << endl;
	printf("Class 1\nC44: %f, C42: %f, C22: %f\n", unc1_44, unc1_42, unc1_22);
	printf("Class 2\nC43: %f, C41: %f, C32: %f, C21: %f\n", unc1_43, unc1_41, unc1_32, unc1_21);
	printf("Class 3\nC33: %f, C31: %f, C11: %f\n", unc1_33, unc1_31, unc1_11);
	printf("Class 0\nC40: %f, C30: %f, C20: %f, C10: %f, C00: %f\n", unc1_40, unc1_30, unc1_20, unc1_10, unc1_00);
	cout << "(2nd neighbor)" << endl;
	printf("Class 1\nC44: %f, C42: %f, C22: %f\n", unc2_44, unc2_42, unc2_22);
	printf("Class 2\nC43: %f, C41: %f, C32: %f, C21: %f\n", unc2_43, unc2_41, unc2_32, unc2_21);
	printf("Class 3\nC33: %f, C31: %f, C11: %f\n", unc2_33, unc2_31, unc2_11);
	printf("Class 0\nC40: %f, C30: %f, C20: %f, C10: %f, C00: %f\n", unc2_40, unc2_30, unc2_20, unc2_10, unc2_00);
}

