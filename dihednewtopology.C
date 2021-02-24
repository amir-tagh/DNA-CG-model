#include <stdio.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <sstream>
#include <string>
#include "FreeEnergyHS.h"

using namespace std;


double hsc::dihed_new_topology(int clusterStart_1, int clusterStart_2, int clusterEnd_1, int clusterEnd_2, double *edihed1, double *edihed2){
    
    double neweDihed_1, neweDihed_2;
    
    
     if ( open_chain ) {
    if( clusterEnd_1 <= clusterStart_1 ){ //account for a break in the chain.
        neweDihed_1   = dihed_forward(&part[clusterStart_1],  &part[(N/2) - 1]);
        neweDihed_1  += dihed_forward(&part[0],               &part[clusterEnd_1 + 1]);
    }else{
        neweDihed_1 = dihed_forward(&part[clusterStart_1],    &part[clusterEnd_1]);
    }
    if( clusterEnd_2 <= clusterStart_2 ){ //account for a break in the chain.
        neweDihed_2   = dihed_forward(&part[clusterStart_2],  &part[N-1]);
        neweDihed_2  += dihed_forward(&part[N/2],             &part[clusterEnd_2 + 1]);
    }else{
        neweDihed_2 = dihed_forward(&part[clusterStart_2],    &part[clusterEnd_2 + 1]);
    }
  }else{
      
       if ( clusterEnd_1 < clusterStart_1 && clusterEnd_2 < clusterStart_2 && clusterEnd_1 == 0 ){
        neweDihed_1  = dihed_forward(&part[clusterStart_1],      &part[N/2 - 1]); 
        neweDihed_2  = dihed_forward(&part[clusterStart_2 + 2],  &part[clusterEnd_2 + 1]); 
     
       }else if ( clusterEnd_1 < clusterStart_1 && clusterStart_2 < clusterEnd_2 ){
        neweDihed_1  = dihed_forward(&part[clusterStart_1],  &part[(N/2) - 1]);
        neweDihed_1 += dihed_forward(&part[3],               &part[clusterEnd_1 + 1]);
        neweDihed_2  = dihed_forward(&part[clusterEnd_2],    &part[(N - 1)]); 
        neweDihed_2 += dihed_forward(&part[(N/2) + 3],       &part[clusterStart_2 + 1]); 
        
       }else if ( clusterEnd_2 == N/2){
        neweDihed_1 = dihed_forward(&part[clusterStart_1],   &part[clusterEnd_1 + 1]);    
        neweDihed_2  = dihed_forward(&part[clusterEnd_2],    &part[(N - 1)]);    
        neweDihed_2 += dihed_forward(&part[(N/2 + 3)],       &part[clusterStart_2 + 1]);  
        
       }else if (clusterStart_1 == 0 ){
         neweDihed_1 = dihed_forward(&part[clusterStart_1],   &part[N/2 - 1]);  
         neweDihed_1 += dihed_forward(&part[3],               &part[clusterEnd_1 + 1]); 
         neweDihed_2 = dihed_forward(&part[clusterEnd_2],     &part[clusterStart_2 + 1]);
         
       }else if (clusterStart_1 - 2 == 0 && clusterEnd_2 - 2 == (N/2)){ 
          neweDihed_1  = dihed_forward(&part[clusterStart_1],      &part[2]);   
          neweDihed_1 += dihed_forward(&part[clusterStart_1 + 1],  &part[clusterEnd_1 + 1]);
          neweDihed_2  = dihed_forward(&part[clusterEnd_2],            &part[(N/2)]);
          neweDihed_2  = dihed_forward(&part[(N/2) + 3],            &part[clusterStart_2 + 1]);
         
       }else if ((clusterStart_1) - 2 == 0 ){
          neweDihed_1  = dihed_forward(&part[clusterStart_1],      &part[2]);   
          neweDihed_1 += dihed_forward(&part[clusterStart_1 + 1],  &part[clusterEnd_1 + 1]); 
          neweDihed_2 = dihed_forward(&part[clusterEnd_2],         &part[clusterStart_2 + 1]);
         
       }else if (clusterEnd_2 < clusterStart_2){
         neweDihed_1  = dihed_forward(&part[clusterStart_1], &part[clusterEnd_1 + 1]);      
         neweDihed_2  = dihed_forward(&part[clusterEnd_2],   &part[(N/2)]);
         neweDihed_2 += dihed_forward(&part[(N/2) + 3],      &part[clusterStart_2 + 1]);
         
       } else {
        neweDihed_1 = dihed_forward(&part[clusterStart_1], &part[clusterEnd_1 + 1]);     
        neweDihed_2 = dihed_forward(&part[clusterEnd_2],   &part[clusterStart_2 + 1]);
        
       }
  }
  
  *edihed1 = neweDihed_1;
  *edihed2 = neweDihed_2;
  
  return 0;
}











