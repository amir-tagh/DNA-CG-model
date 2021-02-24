#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "FreeEnergyHS.h"
#include "tools.h"

#define TINY_LENGTH_VALUE 0.0001
#define TINY_SIN_VALUE    1e-10

using namespace std;


const double sigma2_1 = 10.7 * 10.7;
const double sigma2_2 = 5.35 * 5.35;

//main function to calculate LJ
double hsc::LJ_partEnergy_skipC( Particle *p, Cell *c ){
    
  Particle *q;  
  double    E = 0;
  
    
   /*LJ energy */
   q = c->firstParticle;
   do{    
        if ( q != p  && q->clusterId != 1 ){
            E += LJ_pairEnergy(p, q);
        }
     q = q->next;
   } while( q );
   
   
   /* Loop over particles in the 26 neighbour cells */
   for( int ii = 0; ii < 26; ii++){
       Cell *n = cellHash->getItemByKey( c->neighbours[ii] );
        if( n == NULL ) continue;
        
        q = n->firstParticle;
      
        while( q ){
            if( q->clusterId != 1 )
                E += LJ_pairEnergy(p, q);
            q  = q->next;
        }
   }

   return( E );
   
}
//debug function 
double hsc::LJ_partEnergy( Particle *p, Cell *c ){
    
  Particle *q;  
  double    E = 0;
  double e;
  
   /*LJ energy */
   q = c->firstParticle;
   do{    
       if ( q != p && q->clusterId != 1 ) {
            e  = LJ_pairEnergy(p, q);
            E += e;
        }
     q = q->next;
   } while( q );
   
   
   /* Loop over particles in the 26 neighbour cells */
   for( int ii = 0; ii < 26; ii++){
        Cell *n = cellHash->getItemByKey( c->neighbours[ii] );
        if( n == NULL ) continue;
        
        q = n->firstParticle;
      
        while( q ){
            if( q->clusterId != 1 ){
                
                e  = LJ_pairEnergy(p, q);
                E += e;
            }
            q  = q->next;
        }
   }

   return( E );
   
}

double hsc::LJ_pairEnergy(Particle *p, Particle *q ){
  
  double r2, rcut2;
  double e_lj;
  double sigma2;
  double dr[3];

      
  if(   p->bNext == q
    ||  p->bPrev == q
    ||  p->pNext == q
    ||  p->pPrev == q
    ||  ( SAK_myPair(p, p->myId, N) == q->myId ) ){
                 return(0.0);
  }
  if(p->beadType == 'P' && q->beadType == 'P'){

     if( (  p->pNext 
         && p->pNext->pNext == q )
     ||  (  q->pNext 
         && q->pNext->pNext == p ) ){
         
           return(0.0);
     }
     sigma2  = sigma2_1;
     rcut2   = 144;
  }else{
     sigma2  = sigma2_2;
     rcut2   = 36;    
  }
  
  minVec(p->R, q->R, dr, box->x, box->halfx);
  
  /* distance caculation */
  r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];  
  if( r2 >= rcut2 ){ 
        return( 0.0 );
  }
 
  //Lennard-Jones energy
  double s2r2 = sigma2/r2;
  double s6r6        = s2r2*s2r2*s2r2;
  double s12r12      = s6r6*s6r6;
  
  e_lj  = s12r12 - s6r6;
  e_lj += 0.25;  //correction for cut-shift
  e_lj *= 4.0;
    
  
  return( e_lj );
  
}

double hsc::LJ_pairEnergy_dbg(Particle *p, Particle *q ){
  
  double r2, rcut2;
  double e_lj;
  double sigma2;
  double dr[3];

      
  if(   p->bNext == q
    ||  p->bPrev == q
    ||  p->pNext == q
    ||  p->pPrev == q
    ||  ( SAK_myPair(p, p->myId, N) == q->myId ) ){
      
      
        cerr << p->myId << " " << q->myId << endl;
        cerr << "             "<< p->bNext->myId << endl;
        cerr << "             "<< p->bPrev->myId << endl;
        cerr << "             "<< p->pNext->myId << endl;
        cerr << "             "<< p->pPrev->myId << endl;
        cerr << "             "<< SAK_myPair(p, p->myId, N)<< endl;
        cerr << "Zero here" << endl;
      
                 return(0.0);
  }
    
  
  
  if(p->beadType == 'P' && q->beadType == 'P'){

     if( (  p->pNext 
         && p->pNext->pNext == q )
     ||  (  q->pNext 
         && q->pNext->pNext == p ) ){
         
        cerr << p->myId << " " << q->myId << endl;
        cerr << "Zero here now" << endl;
        
           return(0.0);
     }
     sigma2  = sigma2_1;
     rcut2   = 144;
  }else{
     sigma2  = sigma2_2;
     rcut2   = 36;    
  }
  
  minVec(p->R, q->R, dr, box->x, box->halfx);
  
  /* distance caculation */
  r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];   
  
  if( r2 >= rcut2 ){ 
      
        cerr << p->R[0] << " " << p->R[1] << " " << p->R[2] << endl;
        cerr << q->R[0] << " " << q->R[1] << " " << q->R[2] << endl;
        cerr << r2 << endl;
        cerr << rcut2 << endl;
        
        cerr << "El Zero hereo" << endl;
        return( 0.0 );
  }
 
  //Lennard-Jones energy
  double s2r2 = sigma2/r2;
  double s6r6        = s2r2*s2r2*s2r2;
  double s12r12      = s6r6*s6r6;
  
  e_lj  = s12r12 - s6r6;
  e_lj += 0.25;  //correction for cut-shift
  e_lj *= 4.0;
    
  return( e_lj );
  
}

#if 0
double hsc::LJ_totalEnergy(){
  
  double    Etot, e;
  Particle *p, *q;
  
  Etot = 0.0;
  
  for (int i=0; i<N - 1; i++){  
    p = &part[i];
    for (int j= i + 1; j<N; j++){
        q = &part[j];
        e = LJ_pairEnergy(p, q);
        Etot += e;
        
     }
  }  
 
   
  return(Etot);
}
#else 
double hsc::LJ_totalEnergy(){
  
  double    E;
  Particle *p, *q;
  
  E = 0.0;
  for (int i = 0; i < N; i++){  
      
    p = &part[i];
    Cell *c = p->cell;
    
    q = c->firstParticle;
    do{    
        if ( q > p ) {
            E += LJ_pairEnergy(p, q);
        }
        q = q->next;
    } while( q );
    //cerr << " p LJ: " << p->myId << " q: " << q->myId << " E LJ: " << E << endl;
    
    /* Loop over particles in the 26 neighbour cells */
    for( int ii = 0; ii < 26; ii++){
        Cell *n = cellHash->getItemByKey( c->neighbours[ii] );
        
        //only checking cells with a higher id number:
        //avoids double-counting without having to
        //check individual particle IDs.
        if( n == NULL )
            continue;
        if( c->myId > n->myId ) 
            continue;
        
        q = n->firstParticle;
        while( q ){
            E += LJ_pairEnergy(p, q);
            q  = q->next;
        }
    }
  }
  return( E );
}
#endif
////LJ debug function for rigid-rotation
//check the cluster particles before and after rotaton are the same
double hsc::LJ_pairEnergy_debug(int clusterStart_1, int clusterSize){
    
    
    int    pNr_1_p, pNr_2_c, pNr_3_p, pNr_4_c ;
    double r2, dr[3];
   
    
    for(int i=0; i<clusterSize; i++){
   
        Particle *parts[4];
        Particle *q;
        Particle *p;
      
        pNr_1_p  = (clusterStart_1 + (2*i))%(N/2);
        pNr_2_c  = pNr_1_p + 1;
        pNr_3_p  = N - pNr_1_p -2;
        pNr_4_c  = pNr_3_p +1; 
        
        parts[0] = &part[pNr_1_p];
        parts[1] = &part[pNr_2_c];
        parts[2] = &part[pNr_3_p];
        parts[3] = &part[pNr_4_c];

        for(int pp=0; pp<4; pp++){
         
            p = parts[pp];
            double debug_eLJ;
            int    pNeb;
            pNeb = 0;
            //debug_eLJ = 0.0;
            
            for (int j=0; j<N; j++){
                if( j == p->myId ) continue;
                q = &part[j];
                minVec(p->R, q->R, dr, box->x, box->halfx);
                r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
                double dist = sqrt(r2);
                
                
                if( r2 < rcut2 ){
                    
                    cerr << "\n    outer LJ test" << endl;
                    
                    debug_eLJ = LJ_pairEnergy_dbg(p, q);
                    
                    cerr << " p: " << p->myId << " q: " << q->myId << " dist: " << dist << " rcut_max: " << sqrt(rcut2) << " e_lj: " << debug_eLJ << endl;
                    
                    fprintf(stderr, "LJ %.16f\n\n", debug_eLJ);
                    
                    
                    if( debug_eLJ != 0.0 )
                        pNeb++;
                    
                }
            }  
            //debug_eLJ *= 0.5;
            //cerr << " debug eLJ: " << debug_eLJ << endl;
        }
    }
    
    return( 0. );
     
}
            
            
            
            
            
            
            
            
            



  
