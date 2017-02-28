## Kinetic Monte Carlo Simulator for ABVI Binary Alloy System ##

### Developed by ###

Sky Huang (Chen-Hsi Huang) (https://github.com/skyhuang1208)

A research project on "W-Re alloy evolution under irradiation" conducted by Marian's group in MSE department at UCLA  
(http://jmarian.bol.ucla.edu/)

Please send an email to the following address for more information:  
Sky Huang ( skyhuang@ucla.edu )  
Jaime Marian ( jmarian@ucla.edu )

### Introduction ###

Perform kinetic Monte Carlo simulations for a binary A-B alloy with defects of vacancies and interstitials. The atomic diffusions are assumed via defect mechanisms including: defect jumps, recombinations, annihilations at defect sinks. Frenkel-Pair generation is implemented for defect insertion. The simulation package can be used to study different senarios: radiation process, single-defect, with plane sinks, with surfaces as defect sinks and creation sites. Please see the following journal papers for details:  

[1] Chen-Hsi Huang, Jaime Marian, "A generalized Ising model for studying alloy evolution under irradiation and its use in kinetic Monte Carlo simulations", J. Phys.: Condens. Matter 28, 425201 (2016)  

[2] Chen-Hsi Huang, et al., "Mechanism of Re precipitation in irradiated W-Re alloys from kinetic Monte Carlo simulations", arXiv:1702.03019
(link: https://arxiv.org/abs/1702.03019 )  

For each step, the rates of all possible events (F-P genr, vcc jumps, itl jumps) are calculated first, and an event is randomly picked and performed. After that, recombination and annihilation are checked.   

### How to use ###

Before use: please install gcc/g++ for the newest version (c++11 or newer)

Set up parameters in "kmc_par.h", and then simply type "make". A executable file named exekmc will be generated. Or compile the source files by "g++ kmc_*.cpp -std=c++11". Run the simulations by using command "./exekmc".

### Physical mechanisms ###
1. Frenkel-pair generation  
    * Two atomic sites are randomly selected, which are then replaced by a vacancy and an interstitial with the type of combination of original atom typess
2. Defect jumping  
    * Vacancy: A vacancy jumps in 1st nearest neghbor. The jump rates are calculated using saddle-point model.
    * Self-interstitial atom(SIA): A SIA is assumed to move in 1-dimension; i.e., it sticks to move in one chosen 1st-nn direction, with certain probability to rotate.
    * mixed interstitial atom: It jumps in 2nd-nn, as an effective jumping mechanism of "bridge-mechanism".
3. Defect Sinks  
    Two different types of sinks:  
    * Plane sink: a plane is set to be the defect sink; when a defect jumps onto it, the defect is replaced by an atom selected in the atomic reservoir
    * Surface sink: a new species called vacuum site is introduced. The atoms between vacuum and bulk are surface atoms. When point defect jumps onto a surface atom, it disappears there lefting an intrusion or extrusion. The surfaces also act as vacancy creation sites.

### Program structure ###

* kmc_par.h: store parameters
* kmc_main.cpp: the main file of the program
* kmc_global.h: contains global variables and functions  
    kmc_global.cpp  
* kmc_initial.h: a class used to initialize the simulations  
    kmc_initial.cpp  
* kmc_events.h: a class with all event related variables and functions  
    kmc_events_main.cpp: calculate rates and perform one event  
        kmc_events_genr.cpp: perform F-P generation  
        kmc_events_crcal.cpp: update/calculate vcc creation rate at srf  
        kmc_events_ircal.cpp: calculate itl jumping rates  
        kmc_events_sp.cpp: calculate vcc jumping rates  
            kmc_events_bondecal.cpp: calculate energies using broken bond model  
            kmc_events_recb.cpp: check/perform recombination  


