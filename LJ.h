#include <iostream>
#include <fstream> 
#include <ostream>


using namespace std;
extern int N;


const double sigma2_1 = 10.7 * 10.7;
const double sigma2_2 = 5.35 * 5.35;

inline double LJ_pairEnergy(Particle *p, Particle *q, Box *box){
  
  double r2, rcut2;
  double e_lj;
  double sigma2;
  double dr[3];
 
   
      
  if(p->bNext == q)
                 return(0.0);
  if(p->bPrev == q)
                 return(0.0);
  if(p->pNext == q)
                 return(0.0);
  if(p->pPrev == q)
                 return(0.0);               
  if( SAK_myPair(p, p->myId, N) == q->myId )
                 return(0.0);    
    
  
  
  if(p->beadType == 'P' && q->beadType == 'P'){

     if( p->pNext 
     &&  p->pNext->pNext == q )
                  return(0.0);
     if( q->pNext 
     &&  q->pNext->pNext == p )
                  return(0.0);
     
     sigma2  = sigma2_1;
     rcut2   = 144;
  }else{
     sigma2  = sigma2_2;
     rcut2   = 36;    
  }
  
  minVec(p->R, q->R, dr, box->x, box->halfx);
  
  /* distance caculation */
  r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];   
  if( r2 >= rcut2 ) return( 0.0 );
 
  //Lennard-Jones energy
  double s2r2 = sigma2/r2;
  double s6r6        = s2r2*s2r2*s2r2;
  double s12r12      = s6r6*s6r6;
  
  e_lj  = s12r12 - s6r6;
  e_lj += 0.25;  //correction for cut-shift
  e_lj *= 4.0;
    
  return( e_lj );
  
}

inline double LJ_totalEnergy(Box *box){
  
  double    Etot;
  Particle *p, *q, *part;
  
  Etot = 0.0;
  
  for (int i=0; i<N - 1; i++){  
    p = &part[i];
    for (int j= i + 1; j<N; j++){
      q = &part[j];
      
      double E = LJ_pairEnergy(p, q, box);
      Etot += E;
     }
    }  
 
  return(Etot);
}


  
