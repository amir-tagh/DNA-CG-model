#include <stdio.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <sstream>
#include <string>
#include "FreeEnergyHS.h"
#include "tools.h"

using namespace std;



double hsc::dihed_topology(int clusterStart_1, int clusterStart_2, int clusterEnd_1, int clusterEnd_2){
    
    double newEDihed;
    int    i, j;
    
    i = clusterStart_1;
    if( i < 0 )    i += N/2;
    j = clusterEnd_1+1;  //finish with the base-bead after the given phosphate
    if( j >= N/2 ) j -= N/2;
    
    //i = 0; j = N/2 - 1;
    newEDihed  = dihed_forward_topo(&part[i], &part[j]);
    
    
    i = clusterEnd_2;
    if( i < N/2   ) i += N/2;
    j = clusterStart_2 + 1;//finish with the base-bead after the given phosphate
    if( j >= N   )  j -= N/2;
   
    
    //i = N/2; j = N - 1;
    newEDihed += dihed_forward_topo(&part[i], &part[j]);
    
   return(newEDihed);
}











