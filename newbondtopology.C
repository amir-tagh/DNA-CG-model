//new bond for topology rotation
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


void hsc::bond_old_topology(){
 if ( open_chain ){
       if ( clusterEnd_1 <= clusterStart_1 ){
           oldEBond_1  = bond_bprotateEnergy(&part[clusterStart_1], &part[(N/2) - 1]);
           oldEBond_1 += bond_bprotateEnergy(&part[0],              &part[clusterEnd_1 + 1]); 
           
       }else{
           oldEBond_1 = bond_bprotateEnergy(&part[clusterStart_1],  &part[clusterEnd_1]);
       }
       if ( clusterEnd_2 <= clusterStart_2 ){
           oldEBond_2 = bond_bprotateEnergy(&part[clusterStart_2],  &part[N-1]);
           oldEBond_2 += bond_bprotateEnergy(&part[N/2],            &part[clusterEnd_2 + 1]);
       }else{
           oldEBond_2  = bond_bprotateEnergy(&part[clusterEnd_2],   &part[clusterStart_2 + 1]);
       }
   }else{
          
      if ( clusterEnd_1 < clusterStart_1 && clusterEnd_2 < clusterStart_2 && clusterEnd_1 == 0 ){
          oldEBond_1  = bond_bprotateEnergy(&part[clusterStart_1 - 2],  &part[(N/2) - 2]);
          oldEBond_2  = bond_bprotateEnergy(&part[N - 2],        &part[clusterEnd_2 + 2]);
          oldEBond_2  += bond_bprotateEnergy(&part[clusterStart_2],              &part[N - 2]);
          
      } else if ( clusterEnd_1 < clusterStart_1 && clusterStart_2 < clusterEnd_2){
          oldEBond_1  = bond_bprotateEnergy(&part[clusterStart_1 - 2],  &part[0]);
          oldEBond_1  += bond_bprotateEnergy(&part[0],                  &part[clusterEnd_1]);
          oldEBond_2  = bond_bprotateEnergy(&part[clusterEnd_2],        &part[N - 2]);
          oldEBond_2  += bond_bprotateEnergy(&part[N - 2],              &part[clusterStart_2 + 2]);
      
      }else if (clusterStart_1 == 0 ){
         oldEBond_1  = bond_bprotateEnergy(&part[clusterStart_1],  &part[clusterEnd_1]); 
         oldEBond_2  = bond_bprotateEnergy(&part[clusterStart_2],  &part[clusterEnd_2]);
            
      }else if ( clusterStart_2 < clusterEnd_2 ) {
          oldEBond_1  = bond_bprotateEnergy(&part[clusterStart_1 - 2],  &part[clusterEnd_1]);
          oldEBond_2  = bond_bprotateEnergy(&part[clusterStart_2],      &part[N - 2]);
          oldEBond_2  += bond_bprotateEnergy(&part[N/2],                &part[clusterEnd_2]);
               
      } else {
          oldEBond_1  = bond_bprotateEnergy(&part[clusterStart_1 - 2],  &part[clusterEnd_1]);
          oldEBond_2  = bond_bprotateEnergy(&part[clusterEnd_2],      &part[clusterStart_2 + 2]); 
      }
      
    }
}





void hsc::bond_new_topology(){

 if ( open_chain ){
       if ( clusterEnd_1 <= clusterStart_1 ){
           newEBond_1  = bond_bprotateEnergy(&part[clusterStart_1], &part[(N/2) - 1]);
           newEBond_1 += bond_bprotateEnergy(&part[0],              &part[clusterEnd_1 + 1]); 
           
       }else{
           newEBond_1 = bond_bprotateEnergy(&part[clusterStart_1],  &part[clusterEnd_1]);
       }
       if ( clusterEnd_2 <= clusterStart_2 ){
           newEBond_2 = bond_bprotateEnergy(&part[clusterStart_2],  &part[N-1]);
           newEBond_2 += bond_bprotateEnergy(&part[N/2],            &part[clusterEnd_2 + 1]);
       }else{
           newEBond_2  = bond_bprotateEnergy(&part[clusterEnd_2],   &part[clusterStart_2 + 1]);
       }
   }else{
          
     if ( clusterEnd_1 < clusterStart_1 && clusterEnd_2 < clusterStart_2 && clusterEnd_1 == 0 ){
          newEBond_1  = bond_bprotateEnergy(&part[clusterStart_1 - 2],  &part[(N/2) - 2]);
          newEBond_2  = bond_bprotateEnergy(&part[N - 2],        &part[clusterEnd_2 + 2]);
          newEBond_2  += bond_bprotateEnergy(&part[clusterStart_2],              &part[N - 2]);
          
      } else if ( clusterEnd_1 < clusterStart_1 && clusterStart_2 < clusterEnd_2){
          newEBond_1  = bond_bprotateEnergy(&part[clusterStart_1 - 2],  &part[0]);
          newEBond_1  += bond_bprotateEnergy(&part[0],                  &part[clusterEnd_1]);
          newEBond_2  = bond_bprotateEnergy(&part[clusterEnd_2],        &part[N - 2]);
          newEBond_2  += bond_bprotateEnergy(&part[N - 2],              &part[clusterStart_2 + 2]);
      
      }else if (clusterStart_1 == 0 ){
         newEBond_1  = bond_bprotateEnergy(&part[clusterStart_1],  &part[clusterEnd_1]); 
         newEBond_2  = bond_bprotateEnergy(&part[clusterStart_2],  &part[clusterEnd_2]);
            
      }else if ( clusterStart_2 < clusterEnd_2 ) {
          newEBond_1  = bond_bprotateEnergy(&part[clusterStart_1 - 2],  &part[clusterEnd_1]);
          newEBond_2  = bond_bprotateEnergy(&part[clusterStart_2],      &part[N - 2]);
          newEBond_2  += bond_bprotateEnergy(&part[N/2],                &part[clusterEnd_2]);
               
      } else {
          newEBond_1  = bond_bprotateEnergy(&part[clusterStart_1 - 2],  &part[clusterEnd_1]);
          newEBond_2  = bond_bprotateEnergy(&part[clusterEnd_2],      &part[clusterStart_2 + 2]); 
      }
      
    }
}