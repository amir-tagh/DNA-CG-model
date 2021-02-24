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

void hsc::extension_shrink(){
 
 unsigned int  newCell_id_1;
 double    *Rold;
 double    oldEBond, oldeDihed,oldeSAK;
 double    oldE, newE;
 double    newEBond,neweDihed, neweSAK;
 double    oldeLJ, neweLJ; 
 Particle  *p, *q;
 Particle  **moveParts;
 Cell      **oldCell, **newCell;
 
 ExteShrinkMove++;
 
 int DNASize = N; //select the whole length of DNA
 
 moveParts = new Particle*[DNASize];
 oldCell   = new     Cell*[DNASize];
 newCell   = new     Cell*[DNASize]; 
 Rold      = new   double [(DNASize)*3];
 

 //selecting the first particle (P) and save the old configuration
 
 for( int i = 0; i < N; i++){
     
     p     = &part[i];
     moveParts[i] = p; 
     //cerr << " particle: " << p->myId << moveParts[i]->myId << endl;
     
     
     if(part[0].beadType != 'P'){
         cerr << "First particle should be a phosphor!" << endl;
         exit(8);
     }
   
     oldCell[i]    = p->cell; 
     
     Rold[i*3]     = p->R[0];
     Rold[i*3+1]   = p->R[1];
     Rold[i*3+2]   = p->R[2];
}

oldEBond  = 0;
oldeDihed = 0;
oldeSAK   = 0;
 
 //get the old energies
 for (int i = 0; i < DNASize; i++){
     
    p = &part[i];
 
    oldEBond  += bond_onePartEnergy(p);
    oldeDihed += dihed_onePartEnergy(p);
    oldeSAK   += SAK_pairEnergy(p, &part[SAK_myPair(p, part[i].myId, N)]);
  
 }
  
  oldeLJ    = 0.0;
 //LJ energy 
  for(int i = 0; i < N; i++){
     
     p = &part[i];
    
     
      q = oldCell[i]->firstParticle;
  do{
    if ( q != p){
      oldeLJ    += LJ_pairEnergy(p, q);
    }
    q = q->next;
  } while( q );
  
 
  /* Loop over particles in the 26 neighbour cells */
  for( int ii = 0; ii < 26; ii++){
    Cell *c = cellHash->getItemByKey(oldCell[i]->neighbours[ii]);
    if( c == NULL ) continue;

    q = c->firstParticle; //to update the cell if necessary
    while( q ){
        oldeLJ    += LJ_pairEnergy(p, q);
        q = q->next;
    }
  }
  }
  
 //cerr << " S dist befor: " << (part[2].R[2] - part[0].R[2]) << endl;
  
  //multiply the Z-coordinate by a random number to extend or shrink the DNA (X and Y are not changed) 
  double NewPos_1[3];
  double shrink = 0.2*(ran2(iran)) + 0.2; //a number less than 1 shrinks DNA
  //double extend = 1.2*(ran2(iran));       //a number greater than 1 extends DNA
  cerr << " extshr number: " << shrink << endl;
 
  double dist[3], Z_dist[3];
  
   for(int i = 0; i < DNASize; i++){
                
         p = &part[i];
         //cerr << " X1: " << p->R[0] << " Y1: " << p->R[1] << " Z1: " << p->R[2] << endl;
          
          NewPos_1[0]  = p->R[0];
          NewPos_1[1]  = p->R[1];
          NewPos_1[2]  = p->R[2] * shrink;
          
          p->R[0] = NewPos_1[0];
          p->R[1] = NewPos_1[1];
          p->R[2] = NewPos_1[2];
          //cerr << " X2: " << p->R[0] << " Y2: " << p->R[1] << " Z2: " << p->R[2] << endl;
   }
   
   ofstream myfile;
   myfile.open(extshrinkFileName.c_str(), ios::app);
   writeXYZ_frame( &myfile );
   myfile.close();
   
   //Z-distance between two base pairs after multiplication
  Z_dist[2] = (part[2].R[2] - part[0].R[2]);  
  cerr << " S dist after: " << (part[2].R[2] - part[0].R[2]) << endl; 
  //get the updated length of DNA 
  dist[2] = ((N/4)*(Z_dist[2]));
  cerr << " DNA length " << dist[2] << endl;
   
   
   //update the box information
   box->x[0] = dist[2];
   box->x[1] = dist[2];
   box->x[2] = dist[2];
   box->update();
   
   for(int i = 0; i < DNASize; i++){
       
       p = &part[i];
       intoBox(p->R, box);
   }
   cerr << " box X: " <<  box->x[0] << " box Y: " <<  box->x[1] << " box Z: " <<  box->x[2] << endl;
       
cerr << " get the new cells " << endl;

    for(int i = 0; i < DNASize; i++){  
        
        p = &part[i];
        
        newCell_id_1 =  cellIndex(p->R, box);
        cerr << " newCell Id: " <<  newCell_id_1 << endl;
        if(newCell_id_1 != oldCell[i*4]->myId ){
        newCell[i*4] = moveBetweenCells(p, oldCell[i*4], newCell_id_1, cellHash, box);
        } else {
                 newCell[i*4] = oldCell[i*4];
               }
   }
   
   
  //new energies
  newEBond  = 0;
  neweDihed = 0;
  neweSAK   = 0;
  
  for( int i = 0; i < DNASize; i++){
      
      p = &part[i];
      
  newEBond  += bond_onePartEnergy(p);
  neweDihed += dihed_onePartEnergy(p);
  neweSAK   += SAK_pairEnergy(p, &part[SAK_myPair(p, part[i].myId, N)]);
  
  }
  
  //LJ energy 
  neweLJ    = 0.0;

  for(int i = 0; i < N; i++){
     
      p = &part[i];
     
      q = newCell[i]->firstParticle;
  do{
    if ( q != p){
      neweLJ    += LJ_pairEnergy(p, q);
    }
    q = q->next;
  } while( q );
  
 
  /* Loop over particles in the 26 neighbour cells */
  for( int ii = 0; ii < 26; ii++){
    Cell *c = cellHash->getItemByKey(newCell[i]->neighbours[ii]);
    if( c == NULL ) continue;

    q = c->firstParticle; //to update the cell if necessary
    while( q ){
        neweLJ    += LJ_pairEnergy(p, q);
        q = q->next;
    }
  }
  }
  
  
  oldE =  oldEBond +  oldeDihed + oldeSAK + oldeLJ;
  newE =  newEBond +  neweDihed + neweSAK + neweLJ;
  
//Data collection for transition matrix regarding the extension and shrink  
double criteria = log(ran2(iran));

FILE *TransMatExt;
TransMatExt = fopen("TransMatExt", "a");

if (criteria < (oldE - newE) && (shrink < 1)){  //shrink accept
    
    Shrink++;
    
    if(ExteShrinkMove % Count_Twist == 0){
        fprintf(TransMatExt, "ShrinkAcc %f ExteShrinkMove %f ShrAcc/Tot %f", float(Shrink), float(ExteShrinkMove), float(Shrink)/float(ExteShrinkMove));
        
    }
}else if (criteria < (oldE - newE) && (shrink > 1)){//extension accept
    
    Extension++;
    
    if(ExteShrinkMove % Count_Twist == 0){
        fprintf(TransMatExt, "ExtensionAcc %f ExteShrinkMove %f ExtAcc/Tot %f", float(Extension), float(ExteShrinkMove), float(Extension)/ float(ExteShrinkMove));
    }
    
}else{
    
    ExtShrinkRej++;
    
    if(ExteShrinkMove % Count_Twist == 0){
        fprintf(TransMatExt, "ExtShrinkRej %f ExteShrinkMove %f ExtShrink/Tot %f", float(ExtShrinkRej), float(ExteShrinkMove), float(ExtShrinkRej)/ float(ExteShrinkMove));
        
    }
}
   
//Reject all movements and restore old config.   
  for (int i = 0; i < DNASize; i++){
      
          p = &part[i];
          
          p->R[0] =  Rold[i*3];
          p->R[1] =  Rold[i*3+1];
          p->R[2] =  Rold[i*3+2];
          
          if( newCell[i] != oldCell[i]) {
            moveBetweenCells(p, newCell[i], oldCell[i]->myId, cellHash, box);
          }
  }
  
  for(int i = 0; i < DNASize; i++){
              if( newCell[i]->firstParticle == NULL )
                  if( !cellHash->removeItem(newCell[i]->myId) )
                    newCell[i]       = NULL;
    }
    for(int i = 0; i < DNASize; i++){
              if(newCell[i] && newCell[i]->firstParticle == NULL)
                  delete newCell[i];
    }
  
  delete [] oldCell;
  delete [] newCell;
  delete [] Rold;
  delete [] moveParts;
  
}
  
  
 
 
 
 
 
 
 
 

 
 
 