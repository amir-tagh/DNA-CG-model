#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "FreeEnergyHS.h" 
#include "tools.h"
#include "LJ.h"
#include <fstream>
#include <sstream>
#define TINY_LENGTH_VALUE 0.0001
#define TINY_SIN_VALUE    1e-10
using namespace std;

extern float ran2(long *idum);

//e_local += bond_s_ic(q, p,  r, k, r0 );

double hsc::rigid_bond_onePartEnergy_backward(Particle *p){
  
  Particle  *q;
  double e_local; 
  double r, r0, k;
  double dr[3];
  
  e_local = 0.0;
  
  /* test if particle is a phosphate or base*/
  if( p->beadType == 'P' ){
   
    q  = p->pPrev;
    {  
        r0 = 6.14;
        k  = 20.36;
        minVec(p->R, q->R, dr, box->x, box->halfx);
        r = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
#ifndef INTERCAL
        if( p->bPrev->interCal ){
            e_local += k * (r-(r0+r0_pp_intercal) ) * (r-(r0+r0_pp_intercal)) * 0.5;
        }else
#endif
        e_local += k * (r-r0) * (r-r0) * 0.5;
    }
    
    q  = p->bPrev;
    {
        r0 = 6.09;
        k  = 16.14;
        minVec(p->R, q->R, dr, box->x, box->halfx);
        r = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
        e_local += k * (r-r0) * (r-r0) * 0.5; 
    }
    
  } else {
    
    
    q  = p->bPrev;
    {   
        r0 = 4.07;
        k  = 15.93;
        minVec(p->R, q->R, dr, box->x, box->halfx);
        r = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
        

#ifdef INTERCAL 
        e_local += bond_s_ic(q, p,  r, k, r0 );
        //if( p->interCal ){
          //  e_local += Delta;
        //}
#else
	    e_local += k * (r-r0) * (r-r0) * 0.5; 
#endif
    }
  }
    
  return(e_local);
}


double hsc::rigid_bond_onePartEnergy_forward(Particle *p){
    
  Particle  *q;
  double e_local; 
  double r, r0, k;
  double dr[3];
  
  e_local = 0.0;
  
  /* test if particle is a phosphate or base bead*/
  if( p->beadType == 'P' ){
    
    q  = p->pNext;
    {   
        r0 = 6.14;
        k  = 20.36;
        minVec(p->R, q->R, dr, box->x, box->halfx);
        r = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
#ifndef INTERCAL
        if( p->bNext->interCal ){
            e_local += k * (r-(r0+r0_pp_intercal) ) * (r-(r0+r0_pp_intercal)) * 0.5;
        }else
#endif
        e_local += k * (r-r0) * (r-r0) * 0.5;
    }
    
    
  } else {
    
    q  = p->pNext;
    {
        r0 = 6.09;
        k  = 16.14;
        minVec(p->R, q->R, dr, box->x, box->halfx);
        r = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
        e_local += k * (r-r0) * (r-r0) * 0.5; 
    }
    
    q = p->bNext;
    {   
        r0 = 4.07;
        k  = 15.93;
        minVec(p->R, q->R, dr, box->x, box->halfx);
        r = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
        

#ifdef INTERCAL 
        e_local += bond_s_ic(p, q,  r, k, r0 );
        //if( q->interCal ){
            //e_local += Delta;
        //}
#else
	    e_local += k * (r-r0) * (r-r0) * 0.5; 
#endif
    }
  }
    
  return(e_local);
}










