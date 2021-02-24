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
#include "hashTable.h"
#include "wangLandau.h"


using namespace std;

double hsc::LJ_topology(int clusterStart_1, int clusterSize){
    
    double neweLJ, LJ_pair;
    int    pNr_1_p, pNr_2_c, pNr_3_p, pNr_4_c ;
  
    neweLJ  = 0.0;
    LJ_pair = 0.0;
    
    for(int i=0; i<clusterSize; i++){
   
        Particle *parts[4];
        Particle *q;
      
        pNr_1_p  = (clusterStart_1 + (2*i))%(N/2);
        pNr_2_c  = pNr_1_p + 1;
        pNr_3_p  = N - pNr_1_p -2;
        pNr_4_c  = pNr_3_p +1; 
        
        parts[0] = &part[pNr_1_p];
        parts[1] = &part[pNr_2_c];
        parts[2] = &part[pNr_3_p];
        parts[3] = &part[pNr_4_c];
   
        
        
        for( int pp = 0; pp<4; pp++ ){
            
            q = parts[pp]->cell->firstParticle;
            while( q ){
                /** Do not double-count interactions */  
                if ( q->clusterId != 1 ||  q > parts[pp] ){
                    neweLJ += LJ_pairEnergy(parts[pp], q);
                    
                }
                q = q->next; 
            }
            /* Loop over particles in the 26 neighbour cells */
            for( int ii = 0; ii < 26; ii++){
                 
                 Cell *c = cellHash->getItemByKey(parts[pp]->cell->neighbours[ii]);
                 if ( c == NULL ) continue;
                
                q = c->firstParticle;
                while( q ){
                    
                     if ( q->clusterId != 1 ||  q > parts[pp] ){
                             LJ_pair = LJ_pairEnergy(parts[pp], q);
                             neweLJ += LJ_pair;
                     }
                    q  = q->next;
                }
            }
            
        }
    }
    return neweLJ;
}

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
 