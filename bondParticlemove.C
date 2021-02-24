#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "FreeEnergyHS.h" 
#include "LJ.h"
#include "Tools.h"
#include <fstream>
#include <sstream>
#define TINY_LENGTH_VALUE 0.0001
#define TINY_SIN_VALUE    1e-10
using namespace std;

extern float ran2(long *idum);





double hsc::bond_onePartEnergy(Particle *p){
  
  Particle  *q;
  double e_local; 
  double r, r0, k;
  double dr[3];
  
  e_local = 0.0;
  cerr << " p: " << p->myId << endl;
  
  /* test if particle is a phosphate or base bead*/
  if( p->beadType == 'P' ){
    
    q  = p->bNext;
    if( q )//null for free end DNA
    {                  
        r0 = 5.45;
        k  = 7.04;
        minVec(p->R, q->R, dr, box->x, box->halfx);
        r = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
        e_local  += k * (r-r0) * (r-r0) * 0.5; 
        //cerr << " P: " << p->myId << " B next: " << q->myId << " e: " << e_local << " r: " << r << endl;
    }
    /*q  = p->bPrev;
    if( q )
    {
        r0 = 6.09;
        k  = 16.14;
        minVec(p->R, q->R, dr);
        r = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
        e_local += k * (r-r0) * (r-r0) * 0.5; 
        //cerr << " P: " << p->myId << " B prev: " << q->myId << " e: " << e_local << " r: " << r <<  endl;
    }*/
    
    /*** Phos to phos bonds are symmetric */
    
    q  = p->pNext;
    if( q )
    {
        r0 = 6.14;
        k  = 20.36;
        minVec(p->R, q->R, dr, box->x, box->halfx);
        r = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
        e_local += k * (r-r0) * (r-r0) * 0.5;
        //cerr << " P: " << p->myId << " P next: " << q->myId << " e: " << e_local << " r: " << r << endl;
    }
    /*q  = p->pPrev;
    if( q )
    {
        minVec(p->R, q->R, dr);
        r = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
        e_local += k * (r-r0) * (r-r0) * 0.5;
        //cerr << " P: " << p->myId << " P prev: " << q->myId << " e: " << e_local << " r: " << r << endl;
    }*/
    
  } else {
    
    /*q  = p->pPrev;
    if( q )
    {
        r0 = 5.45;
        k  = 7.04;
        minVec(p->R, q->R, dr);
        r = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
        e_local += k * (r-r0) * (r-r0) * 0.5; //+ was missed
        //cerr << " B: " << p->myId << " P prev: " << q->myId << " e: " << e_local << " r: " << r << endl;
    }*/
        
    q  = p->pNext;
    if( q )
    {
        r0 = 6.09;
        k  = 16.14;
        minVec(p->R, q->R, dr, box->x, box->halfx);
        r = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
        e_local += k * (r-r0) * (r-r0) * 0.5; 
        //cerr << " B: " << p->myId << " P next: " << q->myId << " e: " << e_local << " r: " << r << endl;
    }
    
    /*** B to B bonds are symmetric */
   
    
    /*q  = p->bPrev;
    if( q )
    {
        minVec(p->R, q->R, dr);
        r = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
        e_local += k * (r-r0) * (r-r0) * 0.5; 
        //cerr << " B: " << p->myId << " B prev: " << q->myId << " e: " << e_local << " r: " << r << endl;
    }*/
    q = p->bNext;
    if( q )
    {
        r0 = 4.07;
        k  = 15.93;
        minVec(p->R, q->R, dr, box->x, box->halfx);
        r = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
        e_local += k * (r-r0) * (r-r0) * 0.5; 
        //cerr << " B: " << p->myId << " B next: " << q->myId << " e: " << e_local << " r: " << r << endl;
    }
  }
    
  return(e_local);
}