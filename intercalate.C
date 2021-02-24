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


#define _DEBUG_INTERC


  
//extend by stretching a random single base-pair step
void hsc::intercalate_bpMove(){
 
  int     i, i_b1, i_b2;
  double  oldE, newE, flip;

#ifdef DEBUG_INTERC
  double  oldET, newET;
  oldET = bond_totalEnergy();
#endif


  //pick a base-pair
  i = int(ran2(iran)*(N/4));

  //select the bases: flag acts forward, so we need a base on strand 1, and 
  //the base previous to the opposite base on strand 2.
  i_b1 = 2*i+1;
  i_b2 = part[N-i_b1].bPrev->myId;
  
  //cerr << " part1: " << part[i_b1].myId << " part2: " << part[i_b2].myId << endl;
    
  oldE   = bond_onePartEnergy( &part[i_b1] );
  oldE  += bond_onePartEnergy( &part[i_b2] );


  //flip the status
  if( part[i_b1].interCal == 0 ){
      
    //cant be both S and interc.
    //if ( part[i_b1].S == 1 || part[i_b1].bNext->S == 1 ){
      //  return;
   // }
    
    part[i_b1].interCal = 1;
    part[i_b2].interCal = 1;
    
    //cerr << "tag: " << part[i_b1].interCal << " part: " << part[i_b1].myId << " tag : " << part[i_b2].interCal <<  " part 2: " << part[i_b2].myId << endl; 
    
  }else{
      
    part[i_b1].interCal = 0;
    part[i_b2].interCal = 0;
    //cerr << " else tag part 1: " << part[i_b1].interCal << " else tag part 2: " << part[i_b2].interCal << endl;
    
  }
 
//set flag to zero for debug 
#if 0
#ifdef DEBUG_INTERC

part[i_b1].interCal = 0;
part[i_b2].interCal = 0;

#endif
#endif
  //cerr << " call one part 1 " << endl;
  newE   = bond_onePartEnergy( &part[i_b1] );
  //cerr << " call one part 2 " << endl;
  newE  += bond_onePartEnergy( &part[i_b2] );


#ifdef DEBUG_INTERC
  newET = bond_totalEnergy();
  
  if( ( newET - oldET ) - ( newE - oldE ) >   0.00001 ||
      ( newET - oldET ) - ( newE - oldE ) <  -0.00001 ){
        
            cerr << "Bond error intercalation.\n";
            cerr << "delta eBond    : " << newET  - oldET  << " versus: " << newE - oldE << endl;
            cerr << " newE: " << newE << " oldE: " << oldE << " newET: " << newET << "oldET: " << oldET << endl;
            //cerr << "tag part 1: " << part[i_b1].interCal << " tag part 2: " << part[i_b2].interCal << endl;
            exit(1);
      }
#endif

  //cerr << " old E: " << oldE << " new E: " << newE << " compared to: " << (oldE + mu_interCalator) - newE << endl; 
  //accept or reject
  flip = ran2(iran);
  //if( log(flip) < (oldE + (intercal_epsilon - mu_interCalator) - newE )){ //epsilon_2:binding energy of mu_interCalator //2
  if( log(flip) < (oldE + mu_interCalator) - newE ){  //1                                                                     //mu: chemical potential of intercalator 
      intercal_acc++;
      //cerr << " Inter acc " << intercal_acc << endl;
    
    //global counter of current total energy
    eBondTot  += (newE - oldE); 
    ePotTot   += (newE - oldE);
    
    //counter of number bp steps intercalated
    if( part[i_b1].interCal == 1 )
        intercal_count += 1;
    else
        intercal_count -= 1;
    return;
    
  }
  
  intercal_reject++;
  //cerr << " Inter rej: " << intercal_reject << endl;
  //rejected, so flip back
  if( part[i_b1].interCal == 0 ){
    part[i_b1].interCal = 1;
    part[i_b2].interCal = 1;
  }else{
    part[i_b1].interCal = 0;
    part[i_b2].interCal = 0;
  }
  return;
  
}
  
  
 
 
 
 
 
 
 
 

 
 
 
