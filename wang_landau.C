#include <stdio.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include "FreeEnergyHS.h"
#include "tools.h"
#include "cell.h"
#include "hashTable.h"
#include <time.h>


//adjusted parameters
#define ln_f_limit   1e-8                        //predefined value for adjustment factor to stop the simulation
#define ln_f_initial  2.7182818284590452         //initial value of the adjustment factor(e is chosen because we work with ln(density of states)).
#define flatness_criterion 0.85                  //in the original paper it is set to 0.8
                                                  //the condition to reset the histogram when 
                                                  //min(histogram) > average(histogram)*flatness
                                                  
#define ln_2 0.6931471806
#define reference_mc_steps  1000                  //minimum count in the refrence histogram to include an energy level

long *my_hist;
my_hist = new long [dtermine some size];
double ln_f;                                       //adjustment factor for density of states, as natural logarithm
double energy_list[2000];                          //results array holding the output at each sample point
void zero_dos (void);                              //initialize density of states by setting ln_g[i] = 0 (=> g[i] = 1)
long flat(void);                                   //check the flatness of the histogram agianst the flatness_criterion
void normalize_dos(void);                          //normalize te density of states so that ln[ground state] = ln(2)  *********not sure it is useful******
void num_sums(void);                               //perform the numerical sum over the density of states
void zero_hist_ref(void);                          //initialize the reference histogram
long energy_levels;






//check the flatness of the histogram agains the flatness_criterion

long flat(void){
    
    long i, h;
    long hist_max = 0;
    long hist_min 2147483647;
    double range;
    
    for(i=0; i< energy_levels; i++){
        
        if (hist_ref[i] >= 0){
            
            h = hist[i];
            if (h > hist_max) hist_max = h;
            if (h < hist_min) hist_min = h;
        }
    }
        
    range =((double) (hist_max - hist_min)) / ((double)(hist_max + hist_min));
    
    if (range < (1 - flatness_criterion))
        return -1;
    else
        return 0;
}

//initialize density of states by setting ln_g[i] = 0 (=>g[i] = 1)
void zero_dos(void){
    
    long i;
    
    for(i=0; i<energy_levels; i++){
        
        ln_g[i] = 0.0;
    } 
}

void zero_hist(void){
    
    long i;
    
    for(i=0; i<energy_levels; i++){
        
        hist[i] = 0;
    }
}
     
//initialize refrerence histogram     
void zero_hist_ref(void){
    
    long i;
    
    for(i=0; i<energy_levels; i++){
        
        hist_ref[i] = 0;
    }
}

//main function for wnag-landau method

void wang_landau(void){
    
    
    long flag = 0;
    double p, ln_g0, ln_g1;
    
    
    
}
        
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    
    
    
}
