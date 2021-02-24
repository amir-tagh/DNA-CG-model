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

double hsc::rigidbody_dihed(int clusterStart_1, int clusterStart_2, int clusterEnd_1, int clusterEnd_2, double *edihed1, double *edihed2, double *edihed3, double *edihed4){
    
    double oldeDihed_1, oldeDihed_2;
    double oldeDihed_3, oldeDihed_4;
    
    
   if (clusterStart_1 == 0 ){
      oldeDihed_1  = dihed_forward(&part[clusterStart_1],     &part[(N/2) - 1]);
      oldeDihed_1 += dihed_forward(&part[clusterStart_1 + 3], &part[clusterStart_1 + 1]);
      oldeDihed_2  = dihed_forward(&part[clusterEnd_1],       &part[clusterEnd_1 + 1]);
      oldeDihed_3  = dihed_forward(&part[clusterStart_2],     &part[N - 1]);
      oldeDihed_4  = dihed_forward(&part[clusterEnd_2],       &part[clusterEnd_2 + 1]);
      
  }else if (clusterEnd_1 == 0){
      oldeDihed_1  = dihed_forward(&part[clusterStart_1],   &part[clusterStart_1 + 1]);
      oldeDihed_2  = dihed_forward(&part[clusterEnd_1],     &part[(N/2) -1]);
      oldeDihed_2 += dihed_forward(&part[clusterEnd_1 + 3], &part[clusterEnd_1 + 1]);
      oldeDihed_3  = dihed_forward(&part[clusterStart_2],   &part[clusterStart_2 + 1]);
      oldeDihed_4  = dihed_forward(&part[clusterEnd_2],     &part[clusterEnd_2 + 1]);
      
  }else if (clusterEnd_1 - 2 == 0 && clusterStart_2 - 2 == (N/2)){
      oldeDihed_1  = dihed_forward(&part[clusterStart_1],       &part[clusterStart_1 + 1]);
      oldeDihed_2  = dihed_forward(&part[clusterEnd_1],         &part[0]);
      oldeDihed_2 += dihed_forward(&part[clusterEnd_1 + 1],     &part[clusterEnd_1 + 1]);
      oldeDihed_3  = dihed_forward(&part[clusterStart_2],       &part[clusterStart_2 + 1]);
      oldeDihed_3 += dihed_forward(&part[clusterStart_2 + 1],   &part[clusterStart_2 + 1]);
      oldeDihed_4  = dihed_forward(&part[clusterEnd_2],         &part[clusterEnd_2 + 1]);
      
  }else if (clusterEnd_1 - 2 == 0){
      oldeDihed_1  = dihed_forward(&part[clusterStart_1],   &part[clusterStart_1 + 1]);
      oldeDihed_2  = dihed_forward(&part[clusterEnd_1],     &part[0]);
      oldeDihed_2 += dihed_forward(&part[clusterEnd_1 + 1], &part[clusterEnd_1 + 1]);
      oldeDihed_3  = dihed_forward(&part[clusterStart_2],   &part[clusterStart_2 + 1]);
      oldeDihed_4  = dihed_forward(&part[clusterEnd_2],     &part[clusterEnd_2 + 1]);
      
  }else if (clusterStart_2 == (N/2)){
      oldeDihed_1  = dihed_forward(&part[clusterStart_1],     &part[clusterStart_1 + 1]);
      oldeDihed_2  = dihed_forward(&part[clusterEnd_1],       &part[clusterEnd_1 + 1]);
      oldeDihed_3  = dihed_forward(&part[clusterStart_2],     &part[(N - 1)]);
      oldeDihed_3 += dihed_forward(&part[clusterStart_2 + 3], &part[clusterStart_2 + 1]);
      oldeDihed_4  = dihed_forward(&part[clusterEnd_2],       &part[clusterEnd_2 + 1]);
      
  }else if (clusterStart_2 - 2 == (N/2)){
      oldeDihed_1  = dihed_forward(&part[clusterStart_1],     &part[clusterStart_1 + 1]);
      oldeDihed_2  = dihed_forward(&part[clusterEnd_1],       &part[clusterEnd_1 + 1]);
      oldeDihed_3  = dihed_forward(&part[clusterStart_2],     &part[N - 1]);
      oldeDihed_3 += dihed_forward(&part[clusterStart_2 + 1], &part[clusterStart_2 + 1]);
      oldeDihed_4  = dihed_forward(&part[clusterEnd_2],       &part[clusterEnd_2 + 1]);
      
  }else if (clusterStart_1 - 2 == 0){
      oldeDihed_1  = dihed_forward(&part[clusterStart_1],     &part[(N/2) - 1]);    
      oldeDihed_1 += dihed_forward(&part[clusterStart_1 + 1], &part[clusterStart_1 + 1]);
      oldeDihed_2  = dihed_forward(&part[clusterEnd_1],       &part[clusterEnd_1 + 1]);
      oldeDihed_3  = dihed_forward(&part[clusterStart_2],     &part[N - 1]);
      oldeDihed_4  = dihed_forward(&part[clusterEnd_2],       &part[clusterEnd_2 + 1]);
      
  }else if(clusterStart_2 - 2 == 0){
      oldeDihed_1  = dihed_forward(&part[clusterStart_1],     &part[clusterStart_1 + 1]);
      oldeDihed_2  = dihed_forward(&part[clusterEnd_1],       &part[clusterEnd_1 + 1]);
      oldeDihed_3  = dihed_forward(&part[clusterStart_2],     &part[N - 1]);
      oldeDihed_3 += dihed_forward(&part[clusterStart_2 + 1], &part[clusterStart_2 + 1]);
      oldeDihed_4  = dihed_forward(&part[clusterEnd_2],       &part[clusterEnd_2 + 1]);
      
  }else if (clusterEnd_2 == (N/2)){
      oldeDihed_1  = dihed_forward(&part[clusterStart_1], &part[clusterStart_1 + 1]);
      oldeDihed_2  = dihed_forward(&part[clusterEnd_1],   &part[clusterEnd_1 + 1]);
      oldeDihed_3  = dihed_forward(&part[clusterStart_2], &part[clusterStart_2 + 1]);
      oldeDihed_4  = dihed_forward(&part[clusterEnd_2],   &part[N - 1]);
      oldeDihed_4  += dihed_forward(&part[(N/2) + 3],     &part[clusterEnd_2 + 1]);
      
  }else if (clusterEnd_2 - 2 == 0){
      oldeDihed_1  = dihed_forward(&part[clusterStart_1], &part[clusterStart_1 + 1]);
      oldeDihed_2  = dihed_forward(&part[clusterEnd_1],   &part[clusterEnd_1 + 1]);
      oldeDihed_3  = dihed_forward(&part[clusterStart_2], &part[clusterStart_2 + 1]);
      oldeDihed_4  = dihed_forward(&part[clusterEnd_2],   &part[clusterEnd_2 + 1]);
      
  }else if (clusterEnd_2 - 2 == (N/2)){
      oldeDihed_1  = dihed_forward(&part[clusterStart_1],   &part[clusterStart_1 + 1]);
      oldeDihed_2  = dihed_forward(&part[clusterEnd_1],     &part[clusterEnd_1 + 1]);
      oldeDihed_3  = dihed_forward(&part[clusterStart_2],   &part[clusterStart_2 + 1]);
      oldeDihed_4  = dihed_forward(&part[clusterEnd_2],     &part[N - 1]);
      oldeDihed_4 += dihed_forward(&part[clusterEnd_2 + 1], &part[clusterEnd_2 + 1]);
  
  }else{
      oldeDihed_1  = dihed_forward(&part[clusterStart_1], &part[clusterStart_1 + 1]);
      oldeDihed_2  = dihed_forward(&part[clusterEnd_1],   &part[clusterEnd_1 + 1]);
      oldeDihed_3  = dihed_forward(&part[clusterStart_2], &part[clusterStart_2 + 1]);
      oldeDihed_4  = dihed_forward(&part[clusterEnd_2],   &part[clusterEnd_2 + 1]);
  }
  
  
  *edihed1 = oldeDihed_1;
  *edihed2 = oldeDihed_2;
  *edihed3 = oldeDihed_3;
  *edihed4 = oldeDihed_4;
  
  return 0;
}
  
  
  
  
  
  
  
  
  
  
  