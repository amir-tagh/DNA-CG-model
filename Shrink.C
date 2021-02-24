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
#include "wangLandau.h"


extern float ran2(long *idum);

#define NUM_DELTAS_SHRINK 1024
double delta_buf_shr_acc[NUM_DELTAS_SHRINK];
double delta_buf_shr_rej[NUM_DELTAS_SHRINK];
int    delta_i_shr_acc = 0;
int    delta_i_shr_rej = 0;

void hsc::shrink(ofstream *shr_name, ofstream *Shr_energy_acc, ofstream *Shr_energy_rej){
 
 unsigned int  newCell_id_1;
 double    *Rold;
 double    oldEBond, oldeDihed,oldeSAK;
 double    oldE, newE;
 double    newEBond,neweDihed, neweSAK;
 double    oldeLJ, neweLJ; 
 Particle  *p;
 Cell      **oldCell, **newCell;
 
 ShrinkMove++;
 
 
 int DNASize = N; //select the whole length of DNA
 
 
 oldCell   = new     Cell*[DNASize];
 newCell   = new     Cell*[DNASize]; 
 Rold      = new   double [(DNASize)*3];
 

 //selecting the first particle (P) and save the old configuration
 
 for( int i=0; i<N; i++){
     
     p = &part[i];
     
     oldCell[i]    = p->cell; 
     Rold[i*3]     = p->R[0];
     Rold[i*3+1]   = p->R[1];
     Rold[i*3+2]   = p->R[2];
 }

 oldEBond  = 0;
 oldeDihed = 0;
 oldeSAK   = 0;
 
 /*get the old energies*/
 oldeDihed = totalDihedEnergy();
 oldeSAK   = SAK_totalEnergy();
 oldEBond  = bond_totalEnergy();
 oldeLJ    = LJ_totalEnergy();
  
 double box_oldLz = box->x[2];
 
  
  //update the box information: shrnk step is that needed to bring to the neighbouring
  //extension value in the transition matrix
  box->x[2] = box->x[2] * shrink_step;
  box->update();
  double box_newLz = box->x[2] * shrink_step;
    
  for(int i=0; i<DNASize; i++){
                
         p = &part[i];
         p->R[2] = p->R[2]*shrink_step;
         //intoBox(p->R, box);
   }
   
   if(ShrinkMove % 1 == 0){
   writeOutConfig("shr_config.xyz");
   }
   
  

  for(int i=0; i<DNASize; i++){  
        
        p = &part[i];
        newCell_id_1 =  cellIndex(p->R, box);
        if(newCell_id_1 != oldCell[i]->myId ){
                 newCell[i] = moveBetweenCells(p, oldCell[i], newCell_id_1, cellHash, box);
        } else {
                 newCell[i] = oldCell[i];
        }
   }
   
  //new energies
  neweDihed = totalDihedEnergy();
  neweSAK   = SAK_totalEnergy();
  newEBond  = bond_totalEnergy();
  neweLJ    = LJ_totalEnergy();
  
  
  oldE = oldEBond +  oldeDihed + oldeSAK + oldeLJ;
  newE = newEBond +  neweDihed + neweSAK + neweLJ;
  
//Data collection for transition matrix regarding the extension and shrink  
double criteria   = log(ran2(iran));
double deltaE_shr = newE - oldE;


if (criteria < (oldE - newE)){
#if 0  
   delta_buf_ext_acc[delta_i_ext_acc++] = deltaE_ext;
   if( delta_i_ext_acc >= NUM_DELTAS_EXT ){
       
       for( int i = 0; i < NUM_DELTAS_EXT - 1; i++ ){
             
            *Ext_energy_acc << delta_buf_ext_acc[i] << "\n"; 
       } 
     *Ext_energy_acc << delta_buf_ext_acc[NUM_DELTAS_EXT-1] << endl; 
     delta_i_ext_acc = 0;
   }
#else
  *Shr_energy_acc << deltaE_shr << "\n"; //"\n" doesn't flush like endl does.
#endif
}else{   
    
#if 0  
   delta_buf_ext_rej[delta_i_ext_rej++] = deltaE_ext;
   if( delta_i_ext_rej >= NUM_DELTAS_EXT ){
       
       for( int i = 0; i < NUM_DELTAS_EXT - 1; i++ ){
             
            *Ext_energy_rej << delta_buf_ext_rej[i] << "\n"; 
       } 
     *Ext_energy_rej << delta_buf_ext_rej[NUM_DELTAS_EXT-1] << endl; 
     delta_i_ext_rej = 0;
   }
#else
  *Shr_energy_rej << deltaE_shr << "\n"; 
#endif
}

    /*flat histogram*/
    double x_twist, y_shrink;
    double y_new_shrink;
    double delta_hist, log_accept;
    bool   productionStage;
    int    i, j, ii, jj;
    
    x_twist   = wL->bin_min_twist_update;
    y_shrink  = wL->bin_min_stretch_update;
    //cerr << " SHR y_stretch: " << y_shrink << endl;
    productionStage = wL->productionStage;
    
    
    if( productionStage == true ){
        y_new_shrink  = y_shrink  -  wL->stepsize_stretch;
    }else{
        y_new_shrink  = y_shrink  - (wL->stepsize_stretch * wL->direction_stretch);
        //cerr << " SHR y_shrink: " << y_shrink << " SHR y_new_stretch: " << y_new_shrink << endl; 
    }
    
    if( y_new_shrink >= wL->bin_max_shrink ){
        wL->direction_stretch *= -1;
        wL->getBin_2dhist(x_twist, y_shrink, &i, &j);
        count_hist[i][j] += 1;
        //cerr << " HEREEEEEE " << endl;
    
    }else if( y_new_shrink < wL->bin_min_shrink ){ //change this to bin_min_shrink and set it to stretch bin size in wang class
        wL->direction_stretch *= -1;
        wL->getBin_2dhist(x_twist, y_shrink, &i, &j);
        count_hist[i][j] += 1;
        //cerr << " HEREEE " << " i: " << i << " j: " << j << endl;
        //cerr << " y_new_stretch: " << y_new_stretch << " bin_min_stretch: " << wL->bin_min_stretch << endl;
        
    }else{
        wL->getBin_2dhist(x_twist, y_shrink, &i, &j);
        count_hist[i][j] += 1;
        cerr << " count_hist[i][j]: " << count_hist[i][j] << " i: " << i << " j: " << j << endl;
    }
   
    wL->getBin_2dhist(x_twist, y_new_shrink, &ii, &jj);
    //cerr << " HERE SHR " << " SHR ii: " << ii << " SHR jj: " << jj << endl;
    delta_hist = bias_hist[ii][jj] - bias_hist[i][j];
    log_accept = delta_hist - deltaE_shr;


if (criteria < (log_accept)){        /*accepted*/
    
    ShrinkAcc++;
    cerr << " SHR Acc: " << ShrinkAcc << endl;
    cerr << " delta_hist: " << delta_hist << " deltaE_shr: " << deltaE_shr << " log_accept: " << log_accept << endl;
    
    wL->bin_min_stretch_update = y_new_shrink;
    //cerr << " i: " << i << " j: " << j << " ii: " << ii << " jj: " << jj << endl;
    //count_hist[i][j]    = count_hist[ii][jj];
     
     /*update the energy terms.*/
      ePotTot   += (newE - oldE);
      eDihedTot += (neweDihed - oldeDihed);
      eBondTot  += (newEBond - oldEBond);
      eLJTot    += (neweLJ - oldeLJ);
      
      /*update the box information*/
      box->x[2] = box_newLz;
      box->update(); 
    
    /*collect data*/
    if(ShrinkMove % 1000 == 0){ 
    *shr_name << " ShrinkAcc: " << ShrinkAcc << endl;
    }
    
}else{            /*move rejected*/
    
    ShrinkRej++;
    cerr << " SHR Rej: " << ShrinkRej << endl;
    //cerr << " SHR " << " ii: " << ii << " jj: " << jj << " i: " << i << " j: " << j << endl;
    
    if((ii != i || jj != j) && productionStage == false ){
         bias_hist[ii][jj] += wL->hist_increment;
         //cerr << " SHR bias hist : " << bias_hist[ii][jj] << " ii: " << ii << " jj: " << jj <<
         //" i: " << i << " j: " << j << endl;
    }
    
    if(ShrinkMove % 1000 == 0){ 
    *shr_name << " ShrinkRej: " << ShrinkRej << " TotalMove: " << ShrinkMove << endl;
    }

  /*update the box information*/
  box->x[2] = box_oldLz;
  box->update();  
  
   cerr << " CELL INDEX SHR " << endl;
   testCellIndices();
   
  for (int i=0; i<DNASize; i++){
      
          p = &part[i];
          p->R[0] =  Rold[i*3];
          p->R[1] =  Rold[i*3+1];
          p->R[2] =  Rold[i*3+2];
          
          if( newCell[i] != oldCell[i]) {
            moveBetweenCells(p, newCell[i], oldCell[i]->myId, cellHash, box);
          }
  }
  
  
  //If a cell is empty, take it out of the hash.
  //If you can't take it out of the hash, then it is out already.
  //If it is out already, then it will be freed (deleted) in the second pass.
  for(int i = 0; i < DNASize; i++){
      if( newCell[i]->firstParticle == NULL )
           if( !cellHash->removeItem(newCell[i]->myId) )
                newCell[i]       = NULL;
  }
  for(int i = 0; i < DNASize; i++){
           if(newCell[i] != NULL && newCell[i]->firstParticle == NULL)
                delete newCell[i];
  }
}
  
  delete [] oldCell;
  delete [] newCell;
  delete [] Rold;
   
}
  
  
 
 
 
 
 
 
 
 

 
 
 
