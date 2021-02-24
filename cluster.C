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

//macro to get particle ID in circular chain.
//#define wrappedId(p, end, start) (p > end) ? ((p - end) + start) : (p)



extern float ran2(long *idum);


void hsc::clusterMove(){
 int       pNr_1, pNr_2, pNr_3, pNr_4;
 int       newCell_1, newCell_2, newCell_3, newCell_4;
 int       clusterSize, clusterStart_1; 
 double   *Rold;
 double    oldEBond, oldeDihed, oldeLJ, oldEBond_1, oldEBond_2, oldEBond_3, oldEBond_4, oldeDihed_1, oldeDihed_2, oldeDihed_3, oldeDihed_4,oldeLJ_1, oldeLJ_2, oldeLJ_3, oldeLJ_4;
 double    oldE, newE;
 double    newEBond, neweDihed, neweLJ, newEBond_1, newEBond_2, newEBond_3, newEBond_4, neweDihed_1, neweDihed_2, neweDihed_3, neweDihed_4, neweLJ_1, neweLJ_2, neweLJ_3, neweLJ_4 ;
 double    displace[3];
 Cell    **oldCell;
 Particle *p_1, *p_2, *p_3, *p_4, *q_1, *q_2, *q_3, *q_4;
 

  
 clusterSize    =  4 + 2*(int(ran2(iran)*((N/32))));
 if(clusterSize > 10)
     clusterSize = 10;
 cerr << " clustermove " << " clustersize: " <<  clusterSize << endl;
 
 
 //opposite bases
 clusterStart_1 = 2 * int(ran2(iran)*((N/4)));  
 /*clusterEnd_1   = (clusterStart_1 + clusterSize) % (N/2);
 if( clusterEnd_1 < clusterStart_1 ){
        int tmpi = clusterEnd_1;
        clusterEnd_1   = clusterStart_1;
        clusterStart_1 = tmpi;
 }*/
 
 //opposite bases
 //clusterEnd_2     = (N-1) - clusterStart_1;  //strand 2 cluster ends   opposite start of strand1
 //clusterStart_2   = (N-1) - clusterEnd_1;    //strand 2 cluster starts opposite end of strand1
 
 displace[0]  = (ran2(iran) - 0.5)*maxDisplace;
 displace[1]  = (ran2(iran) - 0.5)*maxDisplace;
 displace[2]  = (ran2(iran) - 0.5)*maxDisplace;
 
  oldCell = new   Cell*[4*(clusterSize+1)];
  Rold    = new double [4*(clusterSize+1)*3];
 
  
  transMoves ++;
  
  /*char fileName_1 [32];
  char fileName_2 [32];
  static int w = 0;
  
  sprintf(fileName_1, "clustertransition1_%03i.xyz", w);
  ofstream outfile_1;
  outfile_1.open(fileName_1, ios::out);
  
  sprintf(fileName_2, "clustertransition2_%03i.xyz",w);
  ofstream outfile_2;
  outfile_2.open(fileName_2, ios::out);
 
  outfile_1 << clusterSize << endl;
  outfile_1 << box->x[0] << " " << box->x[1] << " " << box->x[2] << endl;
 
  outfile_2 << clusterSize << endl;
  outfile_2 << box->x[0] << " " << box->x[1] << " " << box->x[2] << endl;*/
  
  
 for(int i=0; i <=clusterSize; i++){
     
     pNr_1 = (clusterStart_1 + (2*i))%(N/2); 
     p_1   = &part[pNr_1];
     p_1->clusterId  = 1; /* set this flag to mark that the particle is moved */
     
          
   /*first strand*/
   oldCell[i*4]    = p_1->cell; 
   Rold[i*4*3]     = p_1->R[0];
   Rold[i*4*3+1]   = p_1->R[1];
   Rold[i*4*3+2]   = p_1->R[2];
   
  
     pNr_2 = pNr_1 + 1;
     p_2   = &part[pNr_2];
     p_2->clusterId  = 1;
     
   /*second strand*/ 
   oldCell[i*4+1]      = p_2->cell; 
   Rold[(i*4+1)*3]     = p_2->R[0];
   Rold[(i*4+1)*3+1]   = p_2->R[1];
   Rold[(i*4+1)*3+2]   = p_2->R[2];
   
     pNr_3 = N - pNr_1 -2;
     p_3   = &part[pNr_3];
     p_3->clusterId  = 1;
     
    
     oldCell[i*4+2]      = p_3->cell; 
     Rold[(i*4+2)*3]     = p_3->R[0];
     Rold[(i*4+2)*3+1]   = p_3->R[1];
     Rold[(i*4+2)*3+2]   = p_3->R[2];
     
     
     pNr_4 = pNr_3 + 1;
     p_4 = &part[pNr_4];
     p_4->clusterId  = 1;
     
     oldCell[i*4+3]      = p_4->cell; 
     Rold[(i*4+3)*3]     = p_4->R[0];
     Rold[(i*4+3)*3+1]   = p_4->R[1];
     Rold[(i*4+3)*3+2]   = p_4->R[2];
    
         
 }
 //outfile_1.close();
 //outfile_2.close();
 //w++;
  //calc the energy for old config of particle p&q
  oldEBond_1  = 0.0;
  oldeDihed_1 = 0.0;
  //oldeSAK_1   = 0.0;
  oldeLJ_1    = 0.0;
 
  oldEBond_2  = 0.0;
  oldeDihed_2 = 0.0;
  //oldeSAK_2   = 0.0;
  oldeLJ_2    = 0.0; 
  
  oldEBond_3  = 0.0;
  oldeDihed_3 = 0.0;
  //oldeSAK_3   = 0.0;
  oldeLJ_3    = 0.0;
  
  oldEBond_4  = 0.0;
  oldeDihed_4 = 0.0;
  //oldeSAK_4   = 0.0;
  oldeLJ_4    = 0.0;
  
  
  
  //get energies 
 for(int i=0; i <=clusterSize; i++){
   
   
      pNr_1 = (clusterStart_1 + (2*i))%(N/2);    
      p_1   = &part[pNr_1];
     
   //energies for first strand cluster
   oldEBond_1  += bond_onePartEnergy(p_1);
   oldeDihed_1 += dihed_onePartEnergy(p_1);
   //oldeSAK_1   += SAK_pairEnergy(p_1, &part[SAK_myPair(p_1, pNr_1, N)]);
   
  
     pNr_2 = pNr_1 +1;
     p_2   = &part[pNr_2];
   
  
   //energies for second strand cluster
   oldEBond_2  += bond_onePartEnergy(p_2);
   oldeDihed_2 += dihed_onePartEnergy(p_2);
   //oldeSAK_2   += SAK_pairEnergy(p_2, &part[SAK_myPair(p_2, pNr_2, N)]);
   
   pNr_3 = N - pNr_1 -2;    
   p_3   = &part[pNr_3];
     
   //energies for first strand cluster
   oldEBond_1  += bond_onePartEnergy(p_3);
   oldeDihed_1 += dihed_onePartEnergy(p_3);
   //oldeSAK_1   += SAK_pairEnergy(p_1, &part[SAK_myPair(p_3, pNr_3, N)]);
   
   pNr_4 = pNr_3 + 1;   
   p_4   = &part[pNr_4];
     
   //energies for first strand cluster
   oldEBond_1  += bond_onePartEnergy(p_4);
   oldeDihed_1 += dihed_onePartEnergy(p_4);
   //oldeSAK_1   += SAK_pairEnergy(p_4, &part[SAK_myPair(p_4, pNr_4, N)]);
   
   
   /*LJ energy for first strand cluster*/
   /* Loop over particles in this cell to find LJ energies */
   q_1 = oldCell[i*4]->firstParticle;
   while( q_1 ){
     /** Do not consider self-interaction or interaction with other moved particle */  
     if (             q_1 != p_1
        && q_1->clusterId != 1 ) {
            oldeLJ_1 += LJ_pairEnergy(p_1, q_1);
     }
     q_1 = q_1->next; 
   } 
  
    /* LJ energy for second strand cluster*/ 
    q_2 = oldCell[i*4+1]->firstParticle;
   while( q_2 ){
     /** Do not consider self-interaction or interaction with other moved particle */  
     if (             q_2 != p_2
        && q_2->clusterId != 1 ) {
            oldeLJ_2 += LJ_pairEnergy(p_2, q_2);
     }
     q_2 = q_2->next; 
   }   
     
     q_3 = oldCell[i*4+2] ->firstParticle;
   while( q_2 ){
     /** Do not consider self-interaction or interaction with other moved particle */  
     if (             q_3 != p_3
        && q_3->clusterId != 1 ) {
            oldeLJ_3 += LJ_pairEnergy(p_3, q_3);
     }
     q_3 = q_3->next; 
   }   
   
   q_4 = oldCell[i*4+3]->firstParticle;
   while( q_4 ){
     /** Do not consider self-interaction or interaction with other moved particle */  
     if (             q_4 != p_4
        && q_4->clusterId != 1 ) {
            oldeLJ_4 += LJ_pairEnergy(p_4, q_4);
     }
     q_4 = q_4->next; 
   }   
     
     
     
  /* Loop over particles in the 26 neighbour cells */
  for( int ii = 0; ii < 26; ii++){
    q_1 = oldCell[i*4]->neighbours[ii]->firstParticle; //to update the cell if necessary
    while( q_1 ){
      if ( q_1->clusterId != 1 ){
        oldeLJ_1 += LJ_pairEnergy(p_1, q_1);
      }
      q_1 = q_1->next;
    }
  }
 
      /* Loop over particles in the 26 neighbour cells */
  for( int ii = 0; ii < 26; ii++){
    q_2 = oldCell[i*4+1]->neighbours[ii]->firstParticle; //to update the cell if necessary
    while( q_2 ){
      if ( q_2->clusterId != 1 ){
        oldeLJ_2 += LJ_pairEnergy(p_2, q_2);
      }
      q_2 = q_2->next;
    }
  }  
   
    for( int ii = 0; ii < 26; ii++){
    q_3 = oldCell[i*4+2]->neighbours[ii]->firstParticle; //to update the cell if necessary
    while( q_3 ){
      if ( q_3->clusterId != 1 ){
        oldeLJ_3 += LJ_pairEnergy(p_3, q_3);
      }
      q_3 = q_3->next;
    }
  }
    
    for( int ii = 0; ii < 26; ii++){
    q_4 = oldCell[i*4+3]->neighbours[ii]->firstParticle; //to update the cell if necessary
    while( q_4 ){
      if ( q_4->clusterId != 1 ){
        oldeLJ_4 += LJ_pairEnergy(p_4, q_4);
      }
      q_4 = q_4->next;
    }
  }
      
} 
 
  
   oldEBond  = oldEBond_1 + oldEBond_2 + oldEBond_3+ oldEBond_4;
   oldeDihed = oldeDihed_1 + oldeDihed_2 + oldeDihed_3 + oldeDihed_4;
   //oldeSAK   = oldeSAK_1 + oldeSAK_2 + oldeSAK_3 + oldeSAK_4;
   oldeLJ    = oldeLJ_1 + oldeLJ_2 + oldeLJ_3 + oldeLJ_4; 
   
   cerr << " clustermove " << " oldEBond: " << oldEBond << " oldEBond_1: :"<< oldEBond_1 << " oldEBond_2: " <<  oldEBond_2 << " oldeLJ1: " << oldeLJ_1 << endl;
   cerr << " clustermove " << " oldeDihed: " << oldeDihed << " oldeDihed_1: " << oldeDihed_1 << " oldeDihed_2: " << oldeDihed_2 << " oldeLJ2: " << oldeLJ_2 << endl;
   
    /* move the particles for the first & second strand*/ 
  for (int i=0; i<=clusterSize; i++){
    
   
        pNr_1 = (clusterStart_1 + (2*i))%(N/2);
        p_1   = &part[pNr_1];
        
    p_1->R[0] = displace[0];
    p_1->R[1] = displace[1];
    p_1->R[2] = displace[2];
    intoBox(p_1->R, box);
    
       
        pNr_2 = pNr_1 + 1;
        p_2   = &part[pNr_2];
        
        
    p_2->R[0] = displace[0];
    p_2->R[1] = displace[1];
    p_2->R[2] = displace[2];
    intoBox(p_2->R, box);
    
    pNr_3 = N - pNr_1 -2;
    p_3 = &part[pNr_3];
    
    p_3->R[0] = displace[0];
    p_3->R[1] = displace[1];
    p_3->R[2] = displace[2];
    intoBox(p_3->R, box);
    
    pNr_4 = pNr_3 + 1;
    p_4 = &part[pNr_4];
    
    p_4->R[0] = displace[0];
    p_4->R[1] = displace[1];
    p_4->R[2] = displace[2];
    intoBox(p_4->R, box);
     
    /*first strand*/ 
    newCell_1 =  int((p_1->R[0]+box->halfx[0])/box->xCell[0]) +                                
                 int((p_1->R[1]+box->halfx[1])/box->xCell[1])*box->nxCell[0] +
                 int((p_1->R[2]+box->halfx[2])/box->xCell[2])*box->nxCell[0]*box->nxCell[1];
    
    //if( newCell_1 < 0 || newCell_1 >= box->nxCell[0]*box->nxCell[1]*box->nxCell[2] ){
            //cerr << "new Cell_1: " << newCell_1 << endl;
            //exit(8);
                   
    if (&cells[newCell_1] != p_1->cell ){
       p_1->moveBetweenCells(*p_1->cell, cells[newCell_1]);
    }
  
     /*second strand*/ 
      newCell_2 =  int((p_2->R[0]+box->halfx[0])/box->xCell[0]) +                                
                   int((p_2->R[1]+box->halfx[1])/box->xCell[1])*box->nxCell[0] +
                   int((p_2->R[2]+box->halfx[2])/box->xCell[2])*box->nxCell[0]*box->nxCell[1];
    
    //if( newCell_2 < 0 || newCell_2 >= box->nxCell[0]*box->nxCell[1]*box->nxCell[2] ){
            //cerr << "new Cell_2: " << newCell_2 << endl;
            //exit(8);
                   
    if (&cells[newCell_2] != p_2->cell ){
       p_2->moveBetweenCells(*p_2->cell, cells[newCell_2]);
    }
      
      
      newCell_3 = int((p_3->R[0]+box->halfx[0])/box->xCell[0]) +                                
                  int((p_3->R[1]+box->halfx[1])/box->xCell[1])*box->nxCell[0] +
                  int((p_3->R[2]+box->halfx[2])/box->xCell[2])*box->nxCell[0]*box->nxCell[1];
     
                   
    if (&cells[newCell_3] != p_3->cell ){
       p_3->moveBetweenCells(*p_3->cell, cells[newCell_3]);
    }
    
    newCell_4 =  int((p_4->R[0]+box->halfx[0])/box->xCell[0]) +                                
                 int((p_4->R[1]+box->halfx[1])/box->xCell[1])*box->nxCell[0] +
                 int((p_4->R[2]+box->halfx[2])/box->xCell[2])*box->nxCell[0]*box->nxCell[1];
     
                   
    if (&cells[newCell_4] != p_4->cell ){
       p_4->moveBetweenCells(*p_4->cell, cells[newCell_4]);
          }
           
      
}
  
  newEBond_1  = 0.0;
  neweDihed_1 = 0.0;
  //neweSAK_1   = 0.0;
  neweLJ_1    = 0.0;
  
  newEBond_2  = 0.0;
  neweDihed_2 = 0.0;
  //neweSAK_2   = 0.0;
  neweLJ_2    = 0.0;
  
  newEBond_3  = 0.0;
  neweDihed_3 = 0.0;
  //neweSAK_3   = 0.0;
  neweLJ_3    = 0.0;
  
  newEBond_4  = 0.0;
  neweDihed_4 = 0.0;
 // neweSAK_4   = 0.0;
  neweLJ_4    = 0.0;
  
  /* get new energies for the first& second strand */
  for(int i=0; i<=clusterSize; i++){
  
      pNr_1 = (clusterStart_1 + (2*i))%(N/2);
      p_1   = &part[pNr_1];
   
   newEBond_1  += bond_onePartEnergy(p_1);
   neweDihed_1 += dihed_onePartEnergy(p_1);
   cerr << " neweDihed function is called: " << " neweDihed_1: " << neweDihed_1 << endl;
   //neweSAK_1   += SAK_pairEnergy(p_1, &part[SAK_myPair(p_1, pNr_1, N)]);
  
      pNr_2 = pNr_1 + 1;
      p_2   = &part[pNr_2];
   
   /*new energies for the second strand*/
   newEBond_2  += bond_onePartEnergy(p_2);
   neweDihed_2 += dihed_onePartEnergy(p_2);
   //neweSAK_2   += SAK_pairEnergy(p_2, &part[SAK_myPair(p_2, pNr_2, N)]);
   
    pNr_3 = N - pNr_1 -2;
    p_3   = &part[pNr_3];
          
   newEBond_3  += bond_onePartEnergy(p_3);
   neweDihed_3 += dihed_onePartEnergy(p_3);
   //neweSAK_3   += SAK_pairEnergy(p_3, &part[SAK_myPair(p_3, pNr_3_p, N)]);
   
          pNr_4 = pNr_3 + 1;
          p_4 = &part[pNr_4];
           
   newEBond_4  += bond_onePartEnergy(p_4);
   neweDihed_4 += dihed_onePartEnergy(p_4);
   //neweSAK_4   += SAK_pairEnergy(p_4, &part[SAK_myPair(p_4, pNr_4_c, N)]);
   
   
    /*LJ energy for the first strand cluster*/
  q_1 = p_1->cell->firstParticle;
   do{    
       if (             q_1 != p_1
          && q_1->clusterId != 1 ) {
            neweLJ_1 += LJ_pairEnergy(p_1, q_1);
        }
        q_1 = q_1->next;
   } while( q_1 );
   
   /* Loop over particles in the 26 neighbour cells */
   for( int ii = 0; ii < 26; ii++){
      q_1 = p_1->cell->neighbours[ii]->firstParticle;
      while( q_1 ){
        if (  q_1->clusterId != 1 ){
            neweLJ_1 += LJ_pairEnergy(p_1, q_1);
        } 
        q_1 = q_1->next;
      }
   }
  
   /*LJ energies for the second strand cluster*/
     q_2 = p_2->cell->firstParticle;
   do{    
       if (             q_2 != p_2
          && q_2->clusterId != 1 ) {
            neweLJ_2 += LJ_pairEnergy(p_2, q_2);
        }
        q_2 = q_2->next;
   } while( q_2 );
   
   /* Loop over particles in the 26 neighbour cells */
   for( int ii = 0; ii < 26; ii++){
      q_2 = p_2->cell->neighbours[ii]->firstParticle;
      while( q_2 ){
        if (  q_2->clusterId != 1 ){
            neweLJ_2 += LJ_pairEnergy(p_2, q_2);
        } 
        q_2 = q_2->next;
      }
   }   
     
      q_3 = p_3->cell->firstParticle;
   do{    
       if (             q_3 != p_3
          && q_3->clusterId != 1 ) {
            neweLJ_3 += LJ_pairEnergy(p_3, q_3);
        }
        q_3 = q_3->next;
   } while( q_3 );
   
   /* Loop over particles in the 26 neighbour cells */
   for( int ii = 0; ii < 26; ii++){
      q_3 = p_3->cell->neighbours[ii]->firstParticle;
      while( q_3 ){
        if (  q_3->clusterId != 1 ){
            neweLJ_3 += LJ_pairEnergy(p_3, q_3);
        } 
        q_3 = q_3->next;
      }
   }   
     
     
     q_4 = p_4->cell->firstParticle;
   do{    
       if (             q_4 != p_4
          && q_4->clusterId != 1 ) {
            neweLJ_4 += LJ_pairEnergy(p_4, q_4);
        }
        q_4 = q_4->next;
   } while( q_4 );
   
   /* Loop over particles in the 26 neighbour cells */
   for( int ii = 0; ii < 26; ii++){
      q_4 = p_4->cell->neighbours[ii]->firstParticle;
      while( q_4 ){
        if (  q_4->clusterId != 1 ){
            neweLJ_4 += LJ_pairEnergy(p_4, q_4);
        } 
        q_4 = q_4->next;
      }
   }   
     
 }
  
  
   newEBond  = newEBond_1 + newEBond_2 + newEBond_3 + newEBond_4;
   neweDihed = neweDihed_1 + neweDihed_2 + neweDihed_3 + neweDihed_4 ;
   //neweSAK   = neweSAK_1 + neweSAK_2;
   neweLJ = neweLJ_1 + neweLJ_2 + neweLJ_3 + neweLJ_4;
   cerr << " clustermove " << " newEBond: " << newEBond << "newEBond_1: " << newEBond_1 << "  newEBond_2: " <<  newEBond_2 << endl;
   cerr << " clustermove " << " neweDihed: " << neweDihed << " neweDihed_1: " << neweDihed_1 << " neweDihed_2: " << neweDihed_2 << endl;
   cerr << " clustermove " << " neweLJ: " << neweLJ << " neweLJ_1: " << neweLJ_1 << " neweLJ_2: " << neweLJ_2 << endl;
     
   oldE =  oldEBond +  oldeDihed + oldeLJ;
   newE =  newEBond +  neweDihed + neweLJ;
   cerr << " clustermove " << "oldE: " << oldE << " newE: " << newE << endl;
   
  
   
  if (log(ran2(iran)) < (oldE - newE)){
    
    ePotTot   += (newE - oldE);
    eDihedTot += (neweDihed - oldeDihed);
    eBondTot  += (newEBond - oldEBond);
    //eSAKTot   += (neweSAK - oldeSAK);
    eLJTot    += (neweLJ - oldeLJ);
    
    if (fabs(oldE - newE) > 1000.0){
        cerr << " Warning: large energy change " << " oldE: " << oldE << " newE: " << newE << endl;
    }
    
    
    transAccept ++;
      /*first strand*/
        for(int i=0; i<=clusterSize; i++) {
           
            pNr_1 = (clusterStart_1 + (2*i))%(N/2);
            p_1   = &part[pNr_1];
            p_1->clusterId  = 0; 
       
            
            pNr_2 = pNr_1 + 1; 
            p_2   = &part[pNr_2];
            p_2->clusterId  = 0; 
            
             pNr_3 = N - pNr_1 -2;
             p_3   = &part[pNr_3];
             p_3->clusterId  = 0; 
                
             pNr_4 = pNr_3 + 1;
             p_4 = &part[pNr_4];
             p_4->clusterId  = 0; 
            
            }
                 
     }else {
       
       for(int i=0; i<=clusterSize; i++) {
            
           
              pNr_1 = (clusterStart_1 + (2*i))%(N/2);
              p_1   = &part[pNr_1];
            
            
            p_1->R[0] =  Rold[i*4*3];
            p_1->R[1] =  Rold[i*4*3+1];
            p_1->R[2] =  Rold[i*4*3+2];
           
            p_1->clusterId  = 0;   
            
            if(p_1->cell != oldCell[i*4]) {
               p_1->moveBetweenCells(*(p_1->cell), *oldCell[i*4]);
            }
        
             pNr_2 = pNr_1 + 1;
             p_2   = &part[pNr_2];
            
            
             p_2->R[0] =  Rold[(i*4+1)*3];
             p_2->R[1] =  Rold[(i*4+1)*3+1];
             p_2->R[2] =  Rold[(i*4+1)*3+2];
           
             p_2->clusterId  = 0;   
            
            if(p_2->cell != oldCell[i*4+1]) {
               p_2->moveBetweenCells(*(p_2->cell), *oldCell[i*4+1]);
               
           }
           
            pNr_3 = N - pNr_1 -2;
            p_3   = &part[pNr_3];
              
             p_3->R[0] = Rold[(i*4+2)*3];  
             p_3->R[1] = Rold[(i*4+2)*3+1];   
             p_3->R[2] = Rold[(i*4+2)*3+2];
             p_3->clusterId  = 0;
            //cerr << " R'1P'2 " << Rold[(i*4+2)*3] << " R'2P'2 " << Rold[(i*4+2)*3+1] << " R'3P'2 " << Rold[(i*4+2)*3+2] << endl;
            //cerr << " i " << i << " P'3 " <<  pNr_2_p << " myId' " << p_3->myId << endl;
            
            if(p_3->cell != oldCell[i*4+2]) {
               p_3->moveBetweenCells(*(p_3->cell), *oldCell[i*4+2]);
          
           }
           
                pNr_4 = pNr_3 + 1;
                p_4 = &part[pNr_4];
                
                p_4->R[0] = Rold[(i*4+3)*3];
                p_4->R[1] = Rold[(i*4+3)*3+1];
                p_4->R[2] = Rold[(i*4+3)*3+2];
                p_4->clusterId  = 0; 
                 
                //cerr << " R'1C'2 " << Rold[(i*4+2)*3+2] << " R'2C'2 " << Rold[(i*4+3)*3+1] << " R'3C'2 " << Rold[(i*4+3)*3+2] << endl;
                //cerr << " i " << i << " P'4 " << pNr_2_c << " myId' " << p_4->myId << endl;
                
            if(p_4->cell != oldCell[i*4+3]) {
               p_4->moveBetweenCells(*(p_4->cell), *oldCell[i*4+3]);
               
               
     
          }
           
           
           
        }
   
     }
              
    delete [] oldCell;
    delete [] Rold;
              
}

   
   
   
      
      
      
      
  
  
  
  