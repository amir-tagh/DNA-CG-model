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
#include "LJ.h"
#include "cell.h"
#include "hashTable.h"

extern float ran2(long *idum);


#define _DEBUG_SWITCH


  
//extend by stretching a random single base-pair step
void hsc::SDNA_baseMove(){
 
  int     i, i_b1, i_b2;
  double  oldE, newE, flip;

#ifdef DEBUG_SWITCH
  double  oldET, newET;
  oldET = bond_totalEnergy();
#endif


  //pick a base-pair
  i = int(ran2(iran)*(N/4));

  //select the bases: switch to S in pairs.
  i_b1 = 2*i+1;
  i_b2 = part[N-i_b1].myId;
  
    
  oldE   = bond_onePartEnergy( &part[i_b1] );
  oldE  += bond_onePartEnergy( &part[i_b2] );


  //flip the status
  if( part[i_b1].S == 0 ){
    //cant be both S and interc.
    if ( part[i_b1].interCal == 1 || part[i_b1].bPrev->interCal == 1 ){
        return;
    }
    part[i_b1].S = 1;
    part[i_b2].S = 1;
    //cerr << " part 1: " << part[i_b1].S << " part 2: " << part[i_b2].S << endl;
    
  }else{
      
    part[i_b1].S = 0;
    part[i_b2].S = 0;
  }

  part[i_b1].S = 1;
  part[i_b2].S = 1;


  newE   = bond_onePartEnergy( &part[i_b1] );
  newE  += bond_onePartEnergy( &part[i_b2] );


#ifdef DEBUG_SWITCH
  newET = bond_totalEnergy();
  
  if( ( newET   - oldET ) - ( newE - oldE ) >   0.00001 ||
      ( newET   - oldET ) - ( newE - oldE ) <  -0.00001 ){
        
            cerr << "Bond error Switch.\n";
            cerr << "delta eBond    : " << newET  - oldET  << " versus: " << newE - oldE << endl;
            exit(1);
      }
#endif


  //accept or reject
  flip = ran2(iran);
  if( log(flip) < (oldE  - newE) ){
      //cerr << " S-tag acc " << endl;
    
    //global counter of current total energy
    eBondTot  += (newE - oldE); 
    ePotTot   += (newE - oldE);
    
    //counter of number bp steps intercalated
    if( part[i_b1].S == 1 )
        sdna_count += 1;
    else
        sdna_count -= 1;
    return;
    
  }
  //cerr << " S-tag rej " << endl;
  //rejected, so flip back
  
  //cerr << " delta E: " << (oldE  - newE) << " flip: " << log(flip) << endl;
  if( part[i_b1].S == 0 ){
    part[i_b1].S = 1;
    part[i_b2].S = 1;
    
  }else{
    part[i_b1].S = 0;
    part[i_b2].S = 0;
  }
  
  return;
  
}


  
  
 
 
 
 
 
 
 
 

 
 
 
