#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "FreeEnergyHS.h"
#include "tools.h"
#include "wangLandau.h"

#define TINY_LENGTH_VALUE 0.0001
#define TINY_SIN_VALUE    1e-10

//estimate these based on simulation data with unmodified hamiltonians
//purpose is to ensure that chirality can drop as overstretching transition (already identified)
//takes place
//
//together with BB stack breaking this modification should allow non-chiral dna at large 
//extensions.
#define BONDIHED_BPBP_FADE  5.0
#define BONDIHED_PBPB_FADE  5.0


using namespace std;

extern float ran2(long *idum);
double hsc::dihed_forward_topo(Particle *p, Particle *q){


  double    eDihed, e;
  int       i, halfN, pp, qq, rr, ss, finish_i;
  int       ch_start, ch_end;
  int       dh_count, dh_count_max;

  
  i         = p->myId;
  finish_i  = q->myId;
  
  halfN  = N/2;
  eDihed = 0.0;

  pp = i - 3;
  qq = i - 2;
  rr = i - 1;
  ss = i;

  /* which chain */
  if( i < halfN ){
      ch_start = 0;
      ch_end   = halfN;
  }else{
      ch_start = halfN;
      ch_end   = N;
  }
  
  if ( pp < ch_start ) pp += halfN;
  if ( qq < ch_start ) qq += halfN;
  if ( rr < ch_start ) rr += halfN;
 
  //how many dihedrals do we need to evaluate
  dh_count_max = finish_i - i + 1;
  if( dh_count_max < 1 ){
        dh_count_max += halfN;
  }
  
  //have number of beads. number of dihedrals will normally include 
  //leading and trailing calculations, except for periodicity
  dh_count_max += 3;
  if( dh_count_max >= halfN ) dh_count_max = halfN;
  
  eDihed = 0.0;
  for(dh_count = 1; dh_count <= dh_count_max; dh_count++){
      
       e       = dihedEnergy(&part[pp], &part[qq], &part[rr], &part[ss]);
       eDihed += e;
 
                   
       pp++; if ( pp >= ch_end ) pp -= halfN;
       ss++; if ( ss >= ch_end ) ss -= halfN;
       qq++; if ( qq >= ch_end ) qq -= halfN;
       rr++; if ( rr >= ch_end ) rr -= halfN;   
  }
  
  return (eDihed);
}


//dihedral function for topological rotation
double hsc::dihed_forward(Particle *p, Particle *q){


  double    eDihed;
  int       i, halfN, pp, qq, rr, ss, finish_i;
  int       ch_start, ch_end;

  i         = p->myId;
  finish_i  = q->myId;

  halfN  = N/2;
  eDihed = 0.0;

  pp = i - 3;
  qq = i - 2;
  rr = i - 1;
  ss = i;

  /* which chain */
  if( i < halfN ){
      ch_start = 0;
      ch_end   = halfN;
  }else{
      ch_start = halfN;
      ch_end   = N;
  }
 
      if ( pp < ch_start ) pp += halfN;
      if ( qq < ch_start ) qq += halfN;
      if ( rr < ch_start ) rr += halfN;


  while( 1 ){
        double e;
        
        e       = dihedEnergy(&part[pp], &part[qq], &part[rr], &part[ss]);
        eDihed += e;
        

        //most common exit case: have done all four possible dihedral calcs
        pp++; //pp = (pp > ch_end ? (pp - halfN) : pp);
               
          if ( pp > finish_i ){
             return (eDihed);
          }else if ( pp >= ch_end ){
                 pp -= halfN;
          } 
        
        ss++;  if( ss >= ch_end ) ss -= halfN;
        qq++;  if( qq >= ch_end ) qq -= halfN;
        rr++;  if( rr >= ch_end ) rr -= halfN;
  }

  return (eDihed);

}


double hsc::dihed_boundaryEnergy_fromAbove(Particle *p){


  double    eDihed;
  int       i, halfN, pp, qq, rr, ss;
  int       ch_start, ch_end;

  i      = p->myId;
  halfN  = N / 2;
  eDihed = 0.0;

  pp = i - 3;
  qq = i - 2;
  rr = i - 1;
  ss = i;

  /* which chain */
  if( i < halfN ){
      ch_start = 0;
      ch_end   = halfN;
  }else{
      ch_start = halfN;
      ch_end   = N;
  }

  // if toplogy is circular, then wrap the dihedral
  // to the other end of the chain
  if( !open_chain  ){
        if( pp < ch_start ) pp += halfN;
        if( qq < ch_start ) qq += halfN;
        if( rr < ch_start ) rr += halfN;
  } else  {//if topology not circular, then start from beginning of the chain.
        while( pp < ch_start ){
            pp++;
            if ( pp == i ) return( 0. );
            qq++;
            rr++;
            ss++;
        }
  }


  while( 1 ){
        eDihed += dihedEnergy(&part[pp], &part[qq], &part[rr], &part[ss]);

        //most common exit case: have done all three possible dihedral calcs
        pp++;
        if ( pp > i ) {
            return( eDihed );
        }
        //else if ( pp >= ch_end && !initialize_chain_as_linear )
        else if ( pp >= ch_end )  //no consistently defined flag in the code at this time for linear chains.
            pp -= halfN;
        
        //other exit case: if linear DNA, then finish early if we run out of chain
        ss++;
        if( ss >= ch_end ){
            //no consistently defined flag in the code at this time for linear chains.
            //if( initialize_chain_as_linear ){
            //    return( eDihed );
            //}
            //else { 
                ss -= halfN;
            //}
        }
        
        qq++;  if( qq >= ch_end ) qq -= halfN;
        rr++;  if( rr >= ch_end ) rr -= halfN;
    }

  return ( eDihed );

}

double hsc::dihed_onePartEnergy(Particle *p){

  double    eDihed;
  int       i, halfN, pp, qq, rr, ss;
  int       ch_start, ch_end;


  i      = p->myId;
  halfN  = N / 2;
  eDihed = 0.0;
  

  pp = i - 3;
  qq = i - 2;
  rr = i - 1;
  ss = i;

  /* which chain */
  if( i < halfN ){
      ch_start = 0;
      ch_end   = halfN;
  }else{
      ch_start = halfN;
      ch_end   = N;
  }

  // if toplogy is circular, then wrap the dihedral
  // to the other ond of the chain
  if( pp < ch_start ) pp += halfN;
  if( qq < ch_start ) qq += halfN;
  if( rr < ch_start ) rr += halfN;
        
  
  while( 1 ){
        
        double this_d;
        
        this_d  = dihedEnergy(&part[pp], &part[qq], &part[rr], &part[ss]);
        eDihed += this_d;

        //most common exit case: have done all four possible dihedral calcs
        pp++;
        if ( pp == i + 1 ) {
            return( eDihed );
        }
        else if ( pp >= ch_end )
            pp -= halfN;

        //other exit case: wrap if we run out of chain
        ss++;  if( ss >= ch_end ) ss -= halfN;
        qq++;  if( qq >= ch_end ) qq -= halfN;
        rr++;  if( rr >= ch_end ) rr -= halfN;

        }
     
  return( eDihed );
}


double hsc::totalDihedEnergy(){


  double   eDihed;
  Particle *pp, *qq, *rr, *ss, *pStart;

  pp = &part[0];
  qq = &part[1];
  rr = &part[2];
  ss = &part[3];

  pStart = pp;
  eDihed = 0.0;
  
  //loop over the first chain
  do {

    double e = dihedEnergy(pp, qq, rr, ss);
    
    eDihed  += e;//dihedEnergy(pp, qq, rr, ss);

    pp = qq;
    qq = rr;
    rr = ss;
    if( ss->beadType == 'P' )
      ss = ss->bNext;
    else
      ss = ss->pNext;

  } while( pp != pStart && ss != NULL );

  //second chain
  pp = &part[N/2];
  qq = &part[N/2 + 1];
  rr = &part[N/2 + 2];
  ss = &part[N/2 + 3];
  pStart = pp;

  do {

    double e = dihedEnergy(pp, qq, rr, ss);
    
    eDihed  += e;

      pp = qq;
      qq = rr;
      rr = ss;
      if( ss->beadType == 'P' )
        ss = ss->bNext;
      else
        ss = ss->pNext;
      
  } while( pp != pStart  && ss != NULL);
  
  return( eDihed );
}

#define LEN(x)             sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])
#define DOT(x,y)           (x[0]*y[0]+x[1]*y[1]+x[2]*y[2])
#define INVLEN(x)          1.0/sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])
#define VECMUL(x, a)       x[0]*=a;x[1]*=a;x[2]*=a
#define VECMUL_TO(x, y, a) y[0]=a*x[0];y[1]=a*x[1];y[2]=a*x[2]
#define TINY_LENGTH_VALUE 0.0001
#define TINY_LENGTH_VALUE 0.0001

#define PRINT_VEC(x, a)  cerr << a << " " << x[0] << " " << x[1] << " " << x[2] << endl;


inline double dround(double x) {
    return floor(x+0.5);
}

double hsc::dihedEnergy(Particle *p, Particle *q, Particle *r, Particle *s){

  double a[3], b[3], c[3], n1[3],n2[3],m1[3];
  double x,y;
  //double len_aXb, len_bXc;
  double phi, k, deltaPhi, phi0;
  double e_dihed;
  
#ifdef ROTORBOX
  minVecForDihed_rotor(p->R, q->R, r->R, s->R, a, b, c, box->x, box->halfx );
#else
  minVec(q->R, p->R ,a, box->x, box->halfx);
  minVec(r->R, q->R ,b, box->x, box->halfx);
  minVec(s->R, r->R ,c, box->x, box->halfx);  
#endif  
  
 //normal vectors to the planes containing a,b,c;
  vector_product(a, b, n1);
  vector_product(b, c, n2);

  //orthonormal frame
  vector_product(n1, n2, m1 );


  //dot product of n1 & n2
  y = DOT(m1,b)/LEN(b);
  x = DOT(n1,n2);

  phi = atan2(y,x);

  if( p->beadType == 'P' ){
    phi0 =  3.62;
    k     = 25.40;
  } else {
    phi0 =  3.51;
    k     = 27.84;
  }

  deltaPhi = phi - phi0;
  e_dihed = k * ( 1 - cos(deltaPhi));

//#ifdef BONDIHED_PBPB
#ifdef INTERCAL  
  if( p->beadType == 'P' ){
    //reuse n1 to hold the distance between the two Bases
    n1[0] = b[0] + c[0];
    n1[1] = b[1] + c[1];
    n1[2] = b[2] + c[2];
    y     = sqrt(DOT(n1,n1));
    if( y > sStart_x )
        e_dihed /= ( BONDIHED_PBPB_FADE * (y - sStart_x) + 1.);
  }
#endif

//#ifdef BONDIHED_BPBB
#ifdef INTERCAL
  if( p->beadType == 'B' ){
    //reuse n1 to hold the distance between the two Bases
    n1[0] = a[0] + b[0];
    n1[1] = a[1] + b[1];
    n1[2] = a[2] + b[2];
    y     = sqrt(DOT(n1,n1));
    if( y > sStart_x )
        e_dihed /= ( BONDIHED_BPBP_FADE * (y - sStart_x) + 1.);
  }
#endif
  
  
#ifdef SATURATE_DIHED
  if( e_dihed > SATURATE_DIHED_MAX_E )
      e_dihed =  SATURATE_DIHED_MAX_E;
#endif
  
  return( e_dihed );
}


void hsc::dihed_break(){
    
  Particle *p, *q, *r, *s;
  double a[3], b[3], c[3], n1[3];
  double y;
  
  //int number_S = 1;
  //int number_B = 0;
  
   //FILE *f;
   //f = fopen(fname, "a");
   
    for(int i=0; i < (N/2); i++){
        
       p = &part[i> (N/2) ?i-(N/2) :i];
       q = &part[(i+1)> (N/2) ?(i+1)-(N/2) :i+1];
       r = &part[(i+2)> (N/2) ?(i+2)-(N/2) :i+2];
       s = &part[(i+3)> (N/2) ?(i+3)-(N/2) :i+3]; 
    
       minVecForDihed_rotor(p->R, q->R, r->R, s->R, a, b, c, box->x, box->halfx );
       
       if( p->beadType == 'P' ){
   
         n1[0] = b[0] + c[0];
         n1[1] = b[1] + c[1];
         n1[2] = b[2] + c[2];
         y     = sqrt(DOT(n1,n1));
         double BONDIHED_PBPB_RCUT;
         BONDIHED_PBPB_RCUT = adaptive_ext->BonDihed_PBPB_Rcut;
         if( y > BONDIHED_PBPB_RCUT ){
             BONDIHED_PBPB_RCUT_count++;
         }
       }
       
       p = &part[(i+1)> (N/2) ?(i+1)-(N/2) :i+1];
       q = &part[(i+2)> (N/2) ?(i+2)-(N/2) :i+2];
       r = &part[(i+3)> (N/2) ?(i+3)-(N/2) :i+3];
       s = &part[(i+4)> (N/2) ?(i+4)-(N/2) :i+4]; 
       
       minVecForDihed_rotor(p->R, q->R, r->R, s->R, a, b, c, box->x, box->halfx );
       
       if( p->beadType == 'B' ){
   
         n1[0] = a[0] + b[0];
         n1[1] = a[1] + b[1];
         n1[2] = a[2] + b[2];
         y     = sqrt(DOT(n1,n1));
         double BONDIHED_BPBP_RCUT;
         BONDIHED_BPBP_RCUT = adaptive_ext->BonDihed_BPBP_Rcut;
         if( y > BONDIHED_BPBP_RCUT ){
             BONDIHED_BPBP_RCUT_count++;
         }
       }
    }
    cerr << " P: " << BONDIHED_PBPB_RCUT_count << endl;
    cerr << " B: " << BONDIHED_BPBP_RCUT_count << endl;
}
  


 
