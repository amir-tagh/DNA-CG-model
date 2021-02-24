#include <stdio.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <sstream>
#include <string>
#include <time.h>
#include "FreeEnergyHS.h"
#include "tools.h"
#include "hashTable.h"
#include "cell.h"
#include "wangLandau.h"

using namespace std;
  
inline void applyRot(double *m, double *R, double *to){
    
     to[0] = m[0] * R[0] + m[1] * R[1] + m[2] * R[2];
     to[1] = m[3] * R[0] + m[4] * R[1] + m[5] * R[2];
     to[2] = m[6] * R[0] + m[7] * R[1] + m[8] * R[2];  
}


struct Matrix {
    double m[9];
};
typedef struct Matrix Matrix; 


Matrix M, *m;

//rotation matrix corrected based on the orthonormalization method considering the round-off error
void matrix_rot(double *COM_axis, double Theta_Rigid) {
         
 
  //double  theta; 
  double  sTheta;
  double  cTheta;
  double  ux,uy,uz;
  double  length;
  
  Matrix *matrix;
  
  matrix = &M;
  
  length = sqrt((COM_axis[0] * COM_axis[0]) + (COM_axis[1] * COM_axis[1]) + (COM_axis[2] * COM_axis[2])); 
  
  ux = COM_axis[0]/length;
  uy = COM_axis[1]/length;
  uz = COM_axis[2]/length;
  
  sTheta = sin(Theta_Rigid);
  cTheta = cos(Theta_Rigid);
  
  double mcTheta = 1.0 - cTheta;
  
  matrix->m[0] = cTheta + (ux * ux * mcTheta);        matrix->m[1] = (ux * uy) * mcTheta - (uz * sTheta);    matrix->m[2] = ux * uz * mcTheta + (uy * sTheta);
  matrix->m[3] = uy * ux * mcTheta + (uz * sTheta);   matrix->m[4] = cTheta + uy * uy * mcTheta;             matrix->m[5] = uy * uz * mcTheta - (ux * sTheta);
  matrix->m[6] = uz * ux * mcTheta - (uy * sTheta);   matrix->m[7] = uz * uy * mcTheta + (ux * sTheta);      matrix->m[8] = cTheta + (uz * uz * mcTheta); 
  
}


#define _TEST_RB

extern float ran2(long *idum);


void hsc::RigidBodyRotation(){
 unsigned int   newCell_id_1;
 int       clusterSize, clusterStart_1, clusterStart_2, clusterEnd_1, clusterEnd_2; 
 double    *Rold;
 double    oldEBond, oldEBond_1, oldEBond_2, oldEBond_3, oldEBond_4,oldEBond_5, oldEBond_6, oldEBond_7, oldEBond_8;
 double    oldE, newE;
 double    newEBond, newEBond_1, newEBond_2, newEBond_3, newEBond_4,newEBond_5, newEBond_6, newEBond_7, newEBond_8;
 Cell      **oldCell;
 Cell      **newCell;
 Particle  *p_1;
 double    oldeLJ, neweLJ; 
 double    COM_str[3], COM_end[3],COM_axis[3];
 double    oldEDihed, newEDihed;
 
 unsigned int *moveThese;
 
 if( N/8 < 4 ){
     cerr <<"Error, system is too small for cluster move" << endl;
     exit(1);
 }
 
 //cluster always starts with a phosphate (so even number) on chain 1.
 //cluster size in base pairs: minimum 4, maximum N/8 = half of the system.
 //minimum 4 is so that dihedrals don't have to be counted backwards into the block.
 
#ifdef ROTORBOX
 clusterStart_1   = 2 + 2 * int(  ran2(iran)*((N/4)-16) ); //particle id
 
 //find max bps before wrap
 int clusterMax_bps;
 clusterMax_bps   = N/4 - clusterStart_1/2 - 1;
 
 clusterSize      = 4 + int( ran2(iran)*(clusterMax_bps-4) );  //number of base pairs: forbidding wrap.
 
#else
 clusterStart_1   = 2 * int(  ran2(iran)*(N/4) ); 
 clusterSize      = 4 + int( ran2(iran)*((N/8)-4) );  
#endif
 
 halfN = N/2;
    
 oldCell   = new   Cell*[4*(clusterSize)];
 newCell   = new   Cell*[4*(clusterSize)];
 Rold      = new double [4*(clusterSize)*3];
 moveThese = new unsigned int[4*clusterSize];
 
 rigidMoves++;
 
 //do the daptive rigid move
if(rigidMoves > 1000 ){
    
    adaptive_rigid->rigidMoves  = rigidMoves;
    adaptive_rigid->rigidAccept = rigidAccept; 
    adaptive_rigid->Adaptive_step_rigid();
}
 
 for(int i=0; i < clusterSize; i++){
     moveThese[4*i]   =  (clusterStart_1 + (2*i))%(N/2);
     moveThese[4*i+1] =  moveThese[4*i]   + 1;
     moveThese[4*i+2] =  N - moveThese[4*i] - 2;
     moveThese[4*i+3] =  moveThese[4*i+2] + 1;
 }
 clusterEnd_1     =  moveThese[clusterSize*4-4];
 clusterStart_2   =  moveThese[2];
 clusterEnd_2     =  moveThese[clusterSize*4-2];
 
  //select particles for two clusters and save the old configuration 
 for(int i=0; i < 4*clusterSize; i++){
     p_1     = &part[moveThese[i]];
     p_1->clusterId = 1;
     oldCell[i]     = p_1->cell; 
     Rold[i*3]      = p_1->R[0];
     Rold[i*3+1]    = p_1->R[1];
     Rold[i*3+2]    = p_1->R[2];
 } 
 
#ifdef TEST_RB
 //writeOutConfig("test_RB_start.xyz");
// print_bond_totalEnergy("test_RB_start_bondEnergy.dat");
// print_lj_totalEnergy("test_RB_start_ljEnergy.dat");
 
   double testOldEdihed = totalDihedEnergy();
   double testOldEBond  = bond_totalEnergy();
   double testOldELJ    = LJ_totalEnergy();
#endif
 
  oldEDihed = dihed_topology(clusterStart_1, clusterStart_2, clusterEnd_1, clusterEnd_2);
  
  oldEBond_1  = rigid_bond_onePartEnergy_backward(&part[clusterStart_1]);     
  oldEBond_2  = rigid_bond_onePartEnergy_backward(&part[clusterStart_1 + 1]); 
  oldEBond_3  = rigid_bond_onePartEnergy_forward(&part[clusterEnd_1]);       
  oldEBond_4  = rigid_bond_onePartEnergy_forward(&part[clusterEnd_1 + 1]);
  
  oldEBond_5  = rigid_bond_onePartEnergy_forward(&part[clusterStart_2]);
  oldEBond_6  = rigid_bond_onePartEnergy_forward(&part[clusterStart_2 + 1]);
  oldEBond_7  = rigid_bond_onePartEnergy_backward(&part[clusterEnd_2]);
  oldEBond_8  = rigid_bond_onePartEnergy_backward(&part[clusterEnd_2 + 1]);
  
  //LJ energies: using the function from the topology move
  //which includes interactions within the moved cluster
  //this is neccessary because of periodic imaging.
  oldeLJ = LJ_topology(clusterStart_1, clusterSize);
 
  oldEBond  = oldEBond_1  + oldEBond_2  + oldEBond_3  + oldEBond_4 + oldEBond_5 + oldEBond_6 + oldEBond_7 + oldEBond_8;
  
  
   /*select the first and last particles of first cluster*/
    COM_str[0] = 0.5 * (part[clusterStart_1].R[0] +
                        part[clusterStart_2].R[0]);
                         
    COM_str[1] = 0.5 * (part[clusterStart_1].R[1] + 
                        part[clusterStart_2].R[1]);
                             
    COM_str[2] = 0.5 * (part[clusterStart_1].R[2] + 
                        part[clusterStart_2].R[2]);   
    
     
    COM_end[0] = 0.5 * (part[clusterEnd_1].R[0] +
                        part[clusterEnd_2].R[0]);
                            
    COM_end[1] = 0.5 * (part[clusterEnd_1].R[1] +
                        part[clusterEnd_2].R[1]);
                            
    COM_end[2] = 0.5 * (part[clusterEnd_1].R[2] +
                        part[clusterEnd_2].R[2]);
    
    if (COM_end[2] < COM_str[2]){  //z-axis positive
           COM_end[2] += box->x[2];
    }
                    
    COM_axis[0] =  COM_end[0] - COM_str[0];
    COM_axis[1] =  COM_end[1] - COM_str[1];
    COM_axis[2] =  COM_end[2] - COM_str[2];
    
    
    double Theta_Rigid;
    double rotation_angle;
    rotation_angle = adaptive_rigid->theta_rigid;
    
    Theta_Rigid = (rotation_angle*(2*ran2(iran) -1));
    
    matrix_rot(COM_axis, Theta_Rigid);
    
    double NewPos_1[3];

    double  ref_p[3];
    ref_p[0] = part[clusterStart_1].R[0];
    ref_p[1] = part[clusterStart_1].R[1];
    ref_p[2] = part[clusterStart_1].R[2];
    
    
    for(int i=0; i < 4*clusterSize; i++){
           double dr[3];
      
           p_1 = &part[moveThese[i]];
           minVec(p_1->R, ref_p, dr, box->x, box->halfx);
           p_1->R[0] = ref_p[0] + dr[0];
           p_1->R[1] = ref_p[1] + dr[1];
           p_1->R[2] = ref_p[2] + dr[2];
           ref_p[0]  = p_1->R[0];
           ref_p[1]  = p_1->R[1];
           ref_p[2]  = p_1->R[2];
           applyRot(M.m, p_1->R, NewPos_1);
           
           p_1->R[0] = NewPos_1[0];
           p_1->R[1] = NewPos_1[1];
           p_1->R[2] = NewPos_1[2];
           intoBox(p_1->R, box);
           
           newCell_id_1 =  cellIndex(p_1->R, box);
           if ( newCell_id_1 != oldCell[i]->myId ){
                newCell[i] = moveBetweenCells(p_1, oldCell[i], newCell_id_1, cellHash, box);
           }else{
                newCell[i] = oldCell[i];
           }
    }

#ifndef ROTORBOX
  if( rigidMoves % 100000 == 0 ){ 
       
        ofstream myfile;
        myfile.open(rotationFileName.c_str(),ios::app);
        writeXYZ_frame( &myfile);
        myfile.close();
   
  }    
#endif
  
   newEDihed = dihed_topology(clusterStart_1, clusterStart_2, clusterEnd_1, clusterEnd_2);
  
   newEBond_1  = rigid_bond_onePartEnergy_backward(&part[clusterStart_1]); 
   newEBond_2  = rigid_bond_onePartEnergy_backward(&part[clusterStart_1 + 1]); 
   newEBond_3  = rigid_bond_onePartEnergy_forward(&part[clusterEnd_1]);       
   newEBond_4  = rigid_bond_onePartEnergy_forward(&part[clusterEnd_1 + 1]);
  
   newEBond_5  = rigid_bond_onePartEnergy_forward(&part[clusterStart_2]);
   newEBond_6  = rigid_bond_onePartEnergy_forward(&part[clusterStart_2 + 1]);
   newEBond_7  = rigid_bond_onePartEnergy_backward(&part[clusterEnd_2]);
   newEBond_8  = rigid_bond_onePartEnergy_backward(&part[clusterEnd_2 + 1]); 
   neweLJ = LJ_topology(clusterStart_1, clusterSize);
   
   newEBond  = newEBond_1  + newEBond_2  + newEBond_3  + newEBond_4 + newEBond_5 + newEBond_6 + newEBond_7 + newEBond_8;
   
   oldE      =  oldEBond +  oldEDihed + oldeLJ;
   newE      =  newEBond +  newEDihed + neweLJ;
  
#ifdef TEST_RB
   double testNewEdihed    = totalDihedEnergy();
   double testNewEBond     = bond_totalEnergy();
   double testNewELJ       = LJ_totalEnergy();
  
   int    status = 0;

   if(  fabs( ((testNewELJ-testOldELJ)/(neweLJ-oldeLJ)) - 1.) >  1e-6 
     && fabs(  (testNewELJ-testOldELJ)-(neweLJ-oldeLJ))       > 0.001 ) {
       cerr << "LJ ERROR Rigid" << endl;
       
       status = 1;
    }  
   if(( testNewEBond   - testOldEBond ) - (newEBond    - oldEBond) >   0.001 ||
      ( testNewEBond   - testOldEBond ) - (newEBond    - oldEBond) <  -0.001 ){
       
       cerr << "Bond ERROR Rigid" << endl;
       
       status = 1;
    }  

   if(( testNewEdihed   - testOldEdihed ) - (newEDihed    - oldEDihed) >   0.001 ||
      ( testNewEdihed   - testOldEdihed ) - (newEDihed    - oldEDihed) <  -0.001 ){
       
       cerr << "DH ERROR Rigid" << endl;
       status = 1;
    } 
    if( status ){
        cerr << "delta eDihed   : " << testNewEdihed - testOldEdihed << " versus: " << newEDihed - oldEDihed << endl;
        cerr << "delta eLJ      : " << testNewELJ    - testOldELJ    << " versus: " << neweLJ    - oldeLJ << endl;
        cerr << "delta eBond    : " << testNewEBond  - testOldEBond  << " versus: " << newEBond  - oldEBond << endl;
        
        writeOutConfig("test_RB_oops.xyz");
        cerr << "quitting" << endl;
        
       exit( 1 );
    }
   
#endif   
   
  if (log(ran2(iran)) < (oldE - newE)){
    
    ePotTot   += (newE - oldE);
    eDihedTot += (newEDihed - oldEDihed);
    eBondTot  += (newEBond - oldEBond);
    eLJTot    += (neweLJ - oldeLJ);
     
    rigidAccept++;
    
    //clean up the particle's old cell if it is now empty.this section causes a program crash
      /*for(int i=0; i<clusterSize; i++){
          for(int k=0; k<4; k++){
          
        if( oldCell[i*4+k]->firstParticle == NULL ){
           if (cellHash->removeItem( oldCell[i*4+k]->myId )) {
               delete oldCell[i*4+k];
               oldCell[i*4+k] = NULL;
           }
        }
          }
      }*/
   
           for(int i=0; i<4*clusterSize; i++) {
                p_1  = &part[moveThese[i]];
                p_1->clusterId = 0;
            }
     
    }else{
        
        rigidReject++;
        
        for(int i = 0; i < 4*clusterSize; i++) {
              
            p_1  = &part[moveThese[i]];
            p_1->R[0] =  Rold[i*3];
            p_1->R[1] =  Rold[i*3+1];
            p_1->R[2] =  Rold[i*3+2];
                        
            if( newCell[i] != oldCell[i]) {
                moveBetweenCells(p_1, newCell[i], oldCell[i]->myId, cellHash, box);
            }
            p_1->clusterId = 0;
            
          }
          for(int i = 0; i < 4*clusterSize; i++){
              if( newCell[i]->firstParticle == NULL )
                  if( !cellHash->removeItem(newCell[i]->myId) )
                    newCell[i]       = NULL;
          }
          for(int i = 0; i < 4*clusterSize; i++){
              if(newCell[i] && newCell[i]->firstParticle == NULL)
                  delete newCell[i];
          }
       }
  
    delete[] oldCell;
    delete[] newCell;
    delete[] Rold;  
    delete[] moveThese;
}

   
   
   
      
      
      
      
  
  
  
  
   
