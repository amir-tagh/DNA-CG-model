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


#define PN_A_TO_KBT 0.0241432434
#define EXT_ADAPT   1000
#define BOX_LIMIT   0.4


extern float ran2(long *idum);

#define NUM_DELTAS_EXT 1024
double delta_buf_ext[NUM_DELTAS_EXT];
int    delta_i_ext = 0;

void hsc::extension(ofstream *Ext_force){
 
 double    *Rold;
 double    oldEBond, oldeDihed,oldeSAK;
 double    oldE, newE;
 double    newEBond,neweDihed, neweSAK;
 double    oldeLJ, neweLJ; 
 Particle  *p;
  
 ExtensionMove++;
 //cerr << " ExtensionMove " << ExtensionMove << endl;
 
 
 double extension_scale;
 
 Rold      = new   double [3*N];
 
 //selecting the first particle (P) and save the old configuration
 for( int i = 0; i < N; i++){
     p = &part[i];
     Rold[3*i]     = p->R[0];
     Rold[3*i+1]   = p->R[1];
     Rold[3*i+2]   = p->R[2];
  }

  //get the old energies
  oldeDihed = totalDihedEnergy();
  oldeSAK   = SAK_totalEnergy();
  oldEBond  = bond_totalEnergy();
  oldeLJ    = LJ_totalEnergy();

    
  double box_oldLz = box->x[2];

  /*multiply the Z-coordinate by the FIXED step to the nearest other system in the transition matrix*/
  if( ran2(iran) < 0.5 ){
      extension_step *= -1;
  }

  extension_scale = (box->x[2] + extension_step)/box->x[2];
  box->x[2]      *= extension_scale;
  box->update();
  //cerr << " extension_scale: " << extension_scale << endl;
  
   
  // extend the DNA in Z direction
  for(int i = 0; i < N; i++){
         p = &part[i];
         p->R[2]  *=  extension_scale;
         intoBox(p->R, box);
  }
  part[0].R[2] = 1e-9 - box->halfx[2];
  
  Cell *c, *cNext;
  c = cellHash->allCells;
  while( c ){
          cNext = c->allCellsNext;
          cellHash->removeItem_byAddress( c );
          delete c;
          c = cNext;
  }
  
  //Quite expensive to init the hash 
  //delete cellHash;
  //cellHash = new HashTable();
  
  //put particles into cells
  for (int i=0; i<N; i++){
            int icell = cellIndex( part[i].R, box );
            insertToCell(&part[i], cellHash, icell, box);
  }
   
  //new energies
  neweDihed = totalDihedEnergy();
  neweSAK   = SAK_totalEnergy();
  newEBond  = bond_totalEnergy();
  neweLJ    = LJ_totalEnergy();
  
  oldE =  oldEBond +  oldeDihed + oldeSAK + oldeLJ;
  newE =  newEBond +  neweDihed + neweSAK + neweLJ;

  double genForce;
  genForce = (newE - oldE) / extension_step; 
    
#if 0
  delta_buf_ext[delta_i_ext++] = genForce;
  if( delta_i_ext >= NUM_DELTAS_EXT ){
       
       for( int i = 0; i < NUM_DELTAS_EXT - 1; i++ ){
             
            *Ext_force << delta_buf_ext[i] << "\n"; 
       } 
     *Ext_force << delta_buf_ext[NUM_DELTAS_EXT-1] << endl; 
     delta_i_ext = 0;
  }
#endif
    bool   accept = 0;
    double criterion   = log(ran2(iran));

    //include work done to extend/shrink.
    if (criterion < oldE - (newE - extension_step*tension_pn*PN_A_TO_KBT) ){
            accept = 1;
    }
    
    /*make compensations for overcompression of DNA*/
    double relative_length;
    relative_length = float(box->x[2]) / float(adaptive_ext->length_init); 
    if (relative_length < BOX_LIMIT ){
        accept = 0;
    }
    
   if( accept == 0 ){
    //cerr << "Not Accepted " << extension_step << " dE:" << (newE - oldE) << " dW: " << extension_step*tension_pn*PN_A_TO_KBT << endl;
        
    box->x[2] = box_oldLz;
    box->update();
    
    for (int i = 0; i < N; i++){
       part[i].R[0] = Rold[3*i];
       part[i].R[1] = Rold[3*i+1];
       part[i].R[2] = Rold[3*i+2];
    }
    
      //put particles into cells
      Cell *c, *cNext;
      c = cellHash->allCells;
      while( c ){
          cNext = c->allCellsNext;
          cellHash->removeItem_byAddress( c );
          delete c;
          c = cNext;
      }
      for (int i=0; i<N; i++){
            int icell = cellIndex( part[i].R, box );
            insertToCell(&part[i], cellHash, icell, box);
      }
    
   } else { //accepted
       //cerr << " Accepted " << endl;
        ePotTot   = newE;
        eDihedTot = neweDihed;
        eBondTot  = newEBond;
        eLJTot    = neweLJ;
        eSAKTot   = neweSAK;
   }
        
  delete [] Rold;
  
}
  
  
//extend by stretching a random single base-pair step
void hsc::extension_bpStep(){
 
 unsigned int  i_ext;
 double    *Rold;
 double    oldEBond, oldeDihed,oldeSAK;
 double    oldE, newE;
 double    newEBond,neweDihed, neweSAK;
 double    oldeLJ, neweLJ; 
 //Particle  *p;
  
 
 double    stepSize, maxStepSize = 1.;
 
 Rold      = new   double [3*N];
 
 //i_ext     = int(ran2(iran) * (N/4)) * 2;
 i_ext     = int( ran2(iran) * (N/4)) * 2;
 stepSize  = 2 * maxStepSize * ran2(iran) - maxStepSize;
 //cerr << " i_ext: " << i_ext << endl;
 
 
 //selecting the first particle (P) and save the old configuration
 for( int i = 0; i < N; i++){
     Rold[i*3]     = part[i].R[0];
     Rold[i*3+1]   = part[i].R[1];
     Rold[i*3+2]   = part[i].R[2];
  }

  
  //get the old energies
  oldeDihed = totalDihedEnergy();
  oldeSAK   = SAK_totalEnergy();
  oldEBond  = bond_totalEnergy();
  oldeLJ    = LJ_totalEnergy();

    
  double box_oldLz = box->x[2];
  
  box->x[2]       += stepSize;
  box->update();
  
  // extend the DNA in Z direction
  for(int i = i_ext; i < N/2; i++){
        part[i].R[2]     +=  stepSize;
        part[N-1-i].R[2] +=  stepSize;
  }
  part[0].R[2] = 1e-9 - box->halfx[2];
  
  for(int i = 0; i < N; i++ )
      intoBox(part[i].R, box);
  
  Cell *c, *cNext;
  c = cellHash->allCells;
  while( c ){
          cNext = c->allCellsNext;
          cellHash->removeItem_byAddress( c );
          delete c;
          c = cNext;
  }
  
  //put particles into cells
  for (int i=0; i<N; i++){
            int icell = cellIndex( part[i].R, box );
            insertToCell(&part[i], cellHash, icell, box);
  }
   
  //new energies
  neweDihed = totalDihedEnergy();
  neweSAK   = SAK_totalEnergy();
  newEBond  = bond_totalEnergy();
  neweLJ    = LJ_totalEnergy();
  
  oldE =  oldEBond +  oldeDihed + oldeSAK + oldeLJ;
  newE =  newEBond +  neweDihed + neweSAK + neweLJ;

  //double genForce;
  //genForce = (newE - oldE) / stepSize; 
    

    bool accept = 0;
    double criterion   = log(ran2(iran));

        //include work done to extend/shrink.
    if (criterion < oldE - (newE - stepSize*tension_pn*PN_A_TO_KBT) ){
            accept = 1;
    }
    
    /*make compensations for overcompression of DNA*/
    double relative_length;
    relative_length = float(box->x[2]) / float(adaptive_ext->length_init); 
    if (relative_length < BOX_LIMIT ){
        accept = 0;
    }
    
   if( accept == 0 ){
   
   // cerr << "Not Accepted " << extension_step << " dE:" << (newE - oldE) << " dW: " << extension_step*tension_pn*PN_A_TO_KBT << endl;
    //bp_extension_rej++;
    //cerr << " bp_extension_rej: " << bp_extension_rej << endl;
    box->x[2] = box_oldLz;
    box->update();
    
    for (int i = 0; i < N; i++ ){
          part[i].R[0] = Rold[3*i];
          part[i].R[1] = Rold[3*i+1];
          part[i].R[2] = Rold[3*i+2];
    }
    
      Cell *c, *cNext;
      c = cellHash->allCells;
      while( c ){
          cNext = c->allCellsNext;
          cellHash->removeItem_byAddress( c );
          delete c;
          c = cNext;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
      }
      //put particles into cells
      for (int i=0; i<N; i++){
            int icell = cellIndex( part[i].R, box );
            insertToCell(&part[i], cellHash, icell, box);
      }
    
   }
   else{ //accepted
        //bp_extension_acc++;
        //cerr << " bp_extension_acc: " << bp_extension_acc << endl;
        
        ePotTot   = newE;
        eDihedTot = neweDihed;
        eBondTot  = newEBond;
        eLJTot    = neweLJ;
        eSAKTot   = neweSAK;
   }
        
  delete [] Rold;
  
}
  
  
 
 
 
 
 
 
 
 

 
 
 
