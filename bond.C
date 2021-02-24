#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "FreeEnergyHS.h"
#include "tools.h"

using namespace std;

#define TINY_LENGTH_VALUE 0.0001
#define TINY_SIN_VALUE    1e-10
#define B_S_SHIFT         3.32


extern float ran2(long *idum);


//free energy gain to have a base loose on both sides
//#define HALF_K_THREEBODY 75.0  

//minimum range to begin seeing any gain
//#define RCUTON_THREEBODY    2.0
//#define RCUTON_THREEBODY_2  4.0

//well position
//#define R_THREEBODY         7.0

//#define OFFSET_THREEBODY   (HALF_K_THREEBODY * RCUTON_THREEBODY_2)

/*free energy penalties taken from Paul's paper*/
//#define LAMBDA         4.0  //cooperativity parameter, penalizes an interface between states 0 & 1.
//#define DELTA          3.0  //neighbor-exclusion penalty
//#define DELTA_X        3.2

//#define COST_S         2.0

#if 0
double hsc::bond_morse(Particle *p, Particle *q, double r, double k, double r0){
    
    double e_local, delta;
    double D_e; //well depth
    double a;   //width potential
    
    D_e = 0.04*38.94;      //D=0.04eV a=4.45 A*e-1 from PRE 47:1 1993,684 *****eV=38.94 KbT
    a = sqrt(k/(2*D_e));
    
    e_local += D_e*( 1 - exp(-a*(r-r0)))*( 1 - exp(-a*(r-r0)));
    
}
#endif

/*bond function for three state DNA*/
/*state '0' is B-DNA, state '1' is S-DNA and sate '2' is HS-DNA or sigma-DNA*/ 
double hsc::bond_s_ic(Particle *p, Particle *q, double r, double k_b, double r0 ){
 
    double e_local, delta, deltaU;
    
    
    if( p->interCal == 0 ){
      //if( p->S == 1 && q->S == 1 ){
      if( r > sStart_x ){
              e_local  =  k_b*( 0.5*(sStart_x - r0)*(sStart_x - r0) + 
                                    (sStart_x - r0)*(r - sStart_x) ); 
            //try a new format for S:
            //e_local  = exp((sStart_x-r0)) + exp((r0-sStart_x)*5) - 1;
            //threshold after which the bases are not interacting anymore - bases are separated
            if( e_local > U_max )  e_local = U_max;
	    e_local += lambda;
            
            
      //}else if( p->S == 1 || q->S == 1){
            //e_local  =  k_b*(r - r0)*(r - r0)*0.5;  //for exponential form
              //e_local  =  k_b*(sStart_x - r0)*(sStart_x - r0)*0.5;
              //if( e_local > U_max )  e_local = U_max;
            
              //e_local +=  lambda;  //free energy cost of B/S boundary 
            
      }else{
              e_local  =  k_b*(r - r0)*(r - r0)*0.5;
              if( e_local > U_max )  e_local = U_max;
      }
    //}
  
    }else{ //intercal == 1
    
                 delta   = r - intercal_r0;
           
                 //e_local  = intercal_epsilon + exp((delta)) + exp((delta)*5) - 1;   //epsilon:binding energy of an intercalator
                 e_local  = intercal_epsilon + (k_b * delta * delta * 0.5);//1
                 //e_local  = (k_b * delta * delta * 0.5); //2
                 
                 if( e_local > U_max )  e_local = U_max;//1
                 //e_local  = ((sStart_x-r0)) + exp((r0-sStart_x)*5) - 1;
                 //e_local += intercal_epsilon;
                 
                 //intercalated bp, neighboring overstretched bp.
                 //if( r > sStart_x ){//1
                     //e_local += ETA;
                 //}
    }
                 

    return(e_local);
    
}


double hsc::bond_onePartEnergy(Particle *p){

  Particle  *q;
  double e_local;
  double r, r0, k;
  double dr[3];

  e_local = 0.0;
  
  
  /*test if particle is a phosphate or base bead*/
  if( p->beadType == 'P' ){

    q  = p->bNext;
    {                  
        r0 = 5.45;
        k  = 7.04;
        minVec(p->R, q->R, dr, box->x, box->halfx);
        r = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
        e_local  += k * (r-r0) * (r-r0) * 0.5; 
    }
    q  = p->bPrev;
    {
        r0 = 6.09;
        k  = 16.14;
        minVec(p->R, q->R, dr, box->x, box->halfx);
        r = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
        e_local += k * (r-r0) * (r-r0) * 0.5; 
    }
    
    /*** Phos to phos bonds are symmetric */
    r0 = 6.14;
    k  = 20.36;
    q  = p->pNext;
    if( q )
    {
        minVec(p->R, q->R, dr, box->x, box->halfx);
        r = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
#ifndef INTERCAL
        if( p->bNext->interCal ){
            e_local += k * (r-(r0+r0_pp_intercal) ) * (r-(r0+r0_pp_intercal)) * 0.5;
        }else
#endif
        e_local += k * (r-r0) * (r-r0) * 0.5;
    }
    q  = p->pPrev;
    {
        minVec(p->R, q->R, dr, box->x, box->halfx);
        r = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
#ifndef INTERCAL
        if( p->bPrev->interCal ){
            e_local += k * (r-(r0+r0_pp_intercal) ) * (r-(r0+r0_pp_intercal)) * 0.5;
        }else
#endif
        e_local += k * (r-r0) * (r-r0) * 0.5;
    }

  } else {

    q  = p->pPrev;
    {
        r0 = 5.45;
        k  = 7.04;
        minVec(p->R, q->R, dr, box->x, box->halfx);
        r = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
        e_local += k * (r-r0) * (r-r0) * 0.5; 
    }

    q  = p->pNext;
    {
        r0 = 6.09;
        k  = 16.14;
        minVec(p->R, q->R, dr, box->x, box->halfx);
        r = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
        e_local += k * (r-r0) * (r-r0) * 0.5; 
    }

    /*** B to B bonds are symmetric */
    r0 = 4.07;
    k  = 15.93;
    q  = p->bPrev;
    {
        minVec(p->R, q->R, dr, box->x, box->halfx);
        r = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);


#ifdef INTERCAL
        e_local += bond_s_ic(q, p,  r, k, r0 );
        //cerr << " q: " << q->beadType << " p: " << p->beadType << endl;
        if ( p->interCal ){   //neighbor exclusion free energy penalty or 2->2 interface.this should shift the transition toward higher forces.
            e_local += Delta;
        }
        if ( q->interCal == 1 && r >= sStart_x ){
            e_local += ETA; //if the intercalated bp has a stretched neighbor.
        }
#else
        e_local += k * (r-r0) * (r-r0) * 0.5; 
#endif

    }
    /*** B to B bonds are symmetric */
    q = p->bNext;
    {
        minVec(p->R, q->R, dr, box->x, box->halfx);
        r = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);

#ifdef INTERCAL  
       e_local += bond_s_ic(p, q,  r, k, r0 );
       //cerr << " p: " << p->beadType << " q: " << q->beadType << endl;
       if( q->interCal ){
           e_local += Delta; //neighbor exclusion free energy penalty.
       }
       if( p->interCal == 1 &&  r >= sStart_x ){
           e_local += ETA;
       }
#else
	   e_local += k * (r-r0) * (r-r0) * 0.5; 
#endif
    }
  }
     
  return(e_local);
}
 
//Find total bond energy without double-counting.
//THIS FUNCTION ASSUMES CIRCULAR TOPOLOGY
double hsc::bond_totalEnergy(){
  
  Particle *p, *q;
  double e_tot; 
  
  const double r0_pb=5.45, r0_pp=6.14,  r0_bp=6.09,  r0_bb=4.07;
  const double k_pb =7.04, k_pp =20.36,  k_bp=16.14,  k_bb=15.93;
  
  double  r, dr[3];
  

  
  e_tot = 0.0;
  for (int i=0; i<N;){
      
      //phosphate
      p = &part[i++];
      q  = p->bNext;
        minVec(p->R, q->R, dr, box->x, box->halfx);
        r = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
        e_tot  += k_pb * (r-r0_pb) * (r-r0_pb) * 0.5; 
      
      q  = p->pNext;
        minVec(p->R, q->R, dr, box->x, box->halfx);
        r = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
        
#ifndef INTERCAL
        if( p->bNext->interCal ){
            e_tot += k_pp * (r-(r0_pp+r0_pp_intercal) ) * (r-(r0_pp+r0_pp_intercal)) * 0.5;
        }else
#endif
            e_tot += k_pp * (r-r0_pp) * (r-r0_pp) * 0.5;
            
      //base
      p = &part[i++];
      q  = p->pNext;
        minVec(p->R, q->R, dr, box->x, box->halfx);
        r = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
        e_tot += k_bp * (r-r0_bp) * (r-r0_bp) * 0.5;
        
      q = p->bNext;
        minVec(p->R, q->R, dr, box->x, box->halfx);
        r = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);

#ifdef INTERCAL 
        e_tot +=  bond_s_ic(p, q,  r, k_bb, r0_bb );
        if ( p->bPrev->interCal == 1 ){
            e_tot += Delta;
        }
        if( p->interCal == 1 && r >= sStart_x ){
           e_tot += ETA;
       }
#else
	    e_tot += k_bb * (r-r0_bb) * (r-r0_bb) * 0.5; 
#endif
  } 
  return(e_tot);
}
      
      
double hsc::bond_totalEnergy_noBaseStack(){
  
  Particle *q;
  double e_tot; 
  
  e_tot = 0.0;
  for (int i=0; i<N; i++){
      q = &part[i];
      
      e_tot       += rigid_bond_onePartEnergy_forward(q);
      
  }
    
  return(e_tot);
}
 
void hsc::stack_break(){
    
    Particle *p;
    double dr[3];
    double r_1;
    
    for(int i=1; i < (N/2); i+=2){
        
        p = &part[i> (N/2) ?i-(N/2) :i]; 
        
        minVec(p->R, p->bNext->R, dr, box->x, box->halfx); //axis vector
        r_1 = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);//fprintf(f, "bp %i B-domain %i\n", i, number_B);
        
            double BondThreshold = adaptive_ext->BondThreshold;
            if ( r_1 > BondThreshold  ){
                stack_break_count++;
            }
    }
    //cerr <<  "bond: " << stack_break_count << endl;
}
    
    
    
    
    
    
    
    
    
    
    
    
  
