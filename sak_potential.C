#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "FreeEnergyHS.h"
#include "tools.h"
#include "LJ.h"
#include "sak_potential.h"
#include "wangLandau.h"

using namespace std;


double hsc::SAK_totalEnergy(){
  
   double    SAK_totalEnergy;
   Particle *p;
  
   SAK_totalEnergy = 0.0;
  
   for (int i=0; i < N / 2; i++){
      p = &part[i];
      SAK_totalEnergy += SAK_pairEnergy(p, &part[SAK_myPair(p, i, N)]);
   }
    
   return( SAK_totalEnergy);
}


double hsc::SAK_pairEnergy(Particle *p, Particle *q){
  
  int                 i;
  double              r, d ,r_min, i_r, delta, E, *tab;
  double              dr[3]; 
  
 
 /* 0) check particle types (Phos or Base) and choose correct table */
  if (p->beadType == 'P' ){
        tab = (double *)sak_pp;
        d = 0.0065603603603;
        r_min = 11.000000 ;
        
  }else { 
        tab = (double *)sak_bb;
        d = 0.001801801802;
        r_min = 6.2000;
  }
  //r_min = tab[0] ;       
        
 
  minVec(p->R, q->R, dr, box->x, box->halfx);
  r = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);
  
  
  /* 2) find the index i of the table point below r */
  i_r = (r - r_min)/d;
  i = int(i_r);
  delta = i_r - i;
  
  
  //cerr << "SSAK: " << p->beadType << " " << q->beadType << " " << r << " " << r_min << " " << d << endl;
  
  /* 3) E = E[i]*(d2/(d1+d2)) + E[i+1]*(d1/(d1+d2)) */
  if( i >= 999 ){
#ifdef SAK_NOBREAK                                    //SAK_NOBREAK is controlled from the MAKEFILE with the same name flag
    if( p->beadType == 'B' ){
        //ARBITRARY weak attraction to prevent permanent bond-breaking: models non-specific attraction bases to opposite chain
        E = (i - 999) * 0.001;  
    }else{
        E = tab[999];
    }
#else      
    E = tab[999];
#endif
  }else if( i < 0 ){
    E = tab[0];
  }else{
    E = tab[i] * (1.-delta) + tab[i+1] * (delta);
  }
  
  //cerr << "SAK: " << p->beadType << " " << q->beadType << " " << r << " " << i << " " << E << endl;
  
  return(E);
  
}