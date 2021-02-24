#ifndef HAVE_TOOLSH
#define HAVE_TOOLSH

#include <iostream>
#include <fstream>
#include <ostream>
#include "hashTable.h"
#include "FreeEnergyHS.h"
#include <assert.h>

// vector from pos2 to pos1.
//find minimum vector between two points given PBCs
//currently supports rectangular PBCs only

//WARNING ROTORBOX IS CONTROLLED FROM MAKEFILE


////globals to manage rotation of box with z-translation
extern double zrot, sin_zrot, cos_zrot;
extern int    box_turns_per_image, stepsize_box_turns_per_image;

void init_rotorBox(int input_turns_per_image);//function declaration, exposed from tools.C






inline void vector_product(double *a, double *b, double *c){

  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}


inline int SAK_myPair( Particle *p, int i, int N ){

    if( p->beadType == 'P' ){
       return( N - i - 2);
     }else{
       return( N - i );
    }
};

// image point back into box
// with fractional twist
inline void intoBox(double *dr, Box *box){
 
#ifdef ROTORBOX
    
   double hand;
    
   while( dr[2] >= box->halfx[2] ){ 
      dr[2] -= box->x[2];
      switch( box_turns_per_image ){
        case 3:                     /*90 degree counterclockwise*/               
            hand  = dr[1]; 
            dr[1] = dr[0];
            dr[0] = -1*hand;
            break;
        case 2:                     /*180 degree*/
            dr[1] *= -1;
            dr[0] *= -1;
            break;
        case 1:                     /*270 degree or 90 clockwise*/
            hand  = dr[1];
            dr[1] = -1*dr[0];
            dr[0] = hand;           
            break;
        default:
            continue;
      }
   }
   while( dr[2] < -box->halfx[2] ){ 
      dr[2] += box->x[2];
      switch( box_turns_per_image ){
        case 3:
            hand  = dr[1];
            dr[1] = -1*dr[0];
            dr[0] = hand;              
            break;
        case 2:
            dr[1] *= -1;
            dr[0] *= -1;
            break;
        case 1:
            hand  = dr[1];
            dr[1] = dr[0];
            dr[0] = -1*hand;
            break;
        default:
            continue;
      }
   }
    
#else
   while( dr[2] >= box->halfx[2] ){
            dr[2] -= box->x[2];
   }
   while( dr[2] < -1*box->halfx[2] ){
            dr[2] += box->x[2];
   }
#endif

   while( dr[0] >= box->halfx[0] ){
            dr[0] -= box->x[0];
   }
   while( dr[0] < -1*box->halfx[0] ){
            dr[0] += box->x[0];
   }
   while( dr[1] >=  box->halfx[1] ){
            dr[1] -= box->x[1];
   }
   while( dr[1] < -1*box->halfx[1] ){
            dr[1] += box->x[1];
   }
}


#ifdef ROTORBOX

//basic explanation: vectors to images are non-commutative with rotated boxes.
// x->y'  != -(y->x').
//
// solution to this is that imaging done in order to calculate a dihedral
// must be continuous, eg:  x->y' y'->z'  z'->w''  is an OK sequence but
// x->y' y->z z->w'  is not OK.
//
// need a new version of this function which handles 90 degree or 270 degree rotations.
//
inline void minVecForDihed_rotor(double *p, double *q, double *r, double *s, double *a, double *b, double *c, double *bx, double *b_halfx){
    
    double dr[3], hand;
    
    //first get vec p->q
    a[2] = q[2] - p[2]; 
    
    //find vector p->q OR p->q_imaged (if imaging is needed)
    int box_turns = 0;
    while( a[2] >=  b_halfx[2] ){
            a[2]      -= bx[2];
            box_turns -= box_turns_per_image;    
    }
    while( a[2] < -1*b_halfx[2] ){
            a[2]       += bx[2];  
            box_turns += box_turns_per_image;
    }
    a[0] = q[0]-p[0];
    a[1] = q[1]-p[1];
    box_turns = (box_turns + 444) % 4;
    switch( box_turns ){
        case 3:              
            hand  = a[1];
            a[1]  = -1*a[0];
            a[0]  = hand;
            break;
        case 2:
            a[0] *= -1;
            a[1] *= -1;
            break;
        case 1:
            hand  = a[1];
            a[1]  = a[0];
            a[0]  = -1*hand;
            break;
        default:
            break;
    }
    while( a[0] >= b_halfx[0] ) a[0] -= bx[0];
    while( a[0] < -b_halfx[0] ) a[0] += bx[0];
    while( a[1] >= b_halfx[1] ) a[1] -= bx[1];
    while( a[1] < -b_halfx[1] ) a[1] += bx[1];
    
    //dr now holds q or q_imaged.
    dr[0] = p[0] + a[0];
    dr[1] = p[1] + a[1];
    dr[2] = p[2] + a[2];
    
    //find vector q->r_imaged 
    b[2]      = r[2] - dr[2];
    box_turns = 0;
    while( b[2] >=  b_halfx[2] ){
            b[2]      -= bx[2];
            box_turns -= box_turns_per_image;  
    }
    while( b[2] < -1*b_halfx[2] ){
            b[2]       += bx[2];  
            box_turns += box_turns_per_image;  
    }
    b[0] = r[0]-dr[0];
    b[1] = r[1]-dr[1];
    box_turns = (box_turns + 444) % 4;
    switch( box_turns ){
        case 3:
            hand  = b[1];
            b[1]  = -1*b[0];
            b[0]  = hand;
            break;
        case 2:
            b[0] *= -1;
            b[1] *= -1;
            break;
        case 1:
            hand  = b[1];
            b[1]  = b[0];
            b[0]  = -1*hand;
            break;
        default:
            break;
    }
    
    //set dr to r_imaged
    dr[0] = dr[0] + b[0];
    dr[1] = dr[1] + b[1];
    dr[2] = dr[2] + b[2];
    
    //find vector r_imaged --> s_imaged
    c[2]      = s[2] - dr[2];
    box_turns = 0;
    while( c[2] >= b_halfx[2] ){
            c[2]      -= bx[2];
            box_turns -= box_turns_per_image;  
    }
    while( c[2] < -1*b_halfx[2] ){
            c[2]      += bx[2];  
            box_turns += box_turns_per_image;  
    }
    c[0] = s[0]-dr[0];
    c[1] = s[1]-dr[1];
    box_turns = (box_turns + 444) % 4;
    switch( box_turns ){
        case 3:
            hand  = c[1];
            c[1]  = -1*c[0];
            c[0]  = hand;
            break;
        case 2:
            c[0] *= -1;
            c[1] *= -1;
            break;
        case 1:
            hand  = c[1];
            c[1]  = c[0];
            c[0]  = -1*hand;
            break;
        default:
            break;
    }
    
    while( c[0] >= b_halfx[0] ) c[0] -= bx[0];
    while( c[0] < -b_halfx[0] ) c[0] += bx[0];
    while( c[1] >= b_halfx[1] ) c[1] -= bx[1];
    while( c[1] < -b_halfx[1] ) c[1] += bx[1];
    
}
//WARNING:NEVER TOUCH THIS BLOCK OF CODE


//basic explanation: vectors to images are non-commutative with rotated boxes.
// x->y'  != -(y->x').
//
// solution to this is that imaging done in order to calculate a dihedral
// must be continuous, eg:  x->y' y'->z'  z'->w''  is an OK sequence but
// x->y' y->z z->w'  is not OK.
//
// need a new version of this function which handles 90 degree or 270 degree rotations.
//
inline void minVecForDihed(double *p, double *q, double *r, double *s, double *a, double *b, double *c, double *bx, double *b_halfx){
    
    double dr[3];
    
    //first get vec p->q
    a[2] = q[2] - p[2]; 
    
    //find vector p->q OR p->q_imaged (if imaging is needed)
    int img_sign = 1;
    while( a[2] >= b_halfx[2] ){
            a[2]     -= bx[2];
            img_sign *= -1;    ///Multiply by minus 1: only works for 180 degree rotations.
    }
    while( a[2] < -1*b_halfx[2] ){
            a[2]     += bx[2];  
            img_sign *= -1;
    }
    a[0] = img_sign * q[0]-p[0];
    a[1] = img_sign * q[1]-p[1];
    while( a[0] >= b_halfx[0] ) a[0] -= bx[0];
    while( a[0] < -b_halfx[0] ) a[0] += bx[0];
    while( a[1] >= b_halfx[1] ) a[1] -= bx[1];
    while( a[1] < -b_halfx[1] ) a[1] += bx[1];
    
    //dr now holds q or q_imaged.
    dr[0] = p[0] + a[0];
    dr[1] = p[1] + a[1];
    dr[2] = p[2] + a[2];
    
    //find vector q->r_imaged 
    b[2]     = r[2] - dr[2];
    img_sign = 1;
    while( b[2] >= b_halfx[2] ){
            b[2]     -= bx[2];
            img_sign *= -1;
    }
    while( b[2] < -1*b_halfx[2] ){
            b[2]     += bx[2];  
            img_sign *= -1;
    }
    b[0] = img_sign * r[0]-dr[0];
    b[1] = img_sign * r[1]-dr[1];
    while( b[0] >= b_halfx[0] ) b[0] -= bx[0];
    while( b[0] < -b_halfx[0] ) b[0] += bx[0];
    while( b[1] >= b_halfx[1] ) b[1] -= bx[1];
    while( b[1] < -b_halfx[1] ) b[1] += bx[1];
    
    //set dr to r_imaged
    dr[0] = dr[0] + b[0];
    dr[1] = dr[1] + b[1];
    dr[2] = dr[2] + b[2];
    
    //find vector r_imaged --> s_imaged
    c[2]     = s[2] - dr[2];
    img_sign = 1;
    while( c[2] >= b_halfx[2] ){
            c[2]     -= bx[2];
            img_sign *= -1;
    }
    while( c[2] < -1*b_halfx[2] ){
            c[2]     += bx[2];  
            img_sign *= -1;
    }
    c[0] = img_sign * s[0]-dr[0];
    c[1] = img_sign * s[1]-dr[1];
    while( c[0] >= b_halfx[0] ) c[0] -= bx[0];
    while( c[0] < -b_halfx[0] ) c[0] += bx[0];
    while( c[1] >= b_halfx[1] ) c[1] -= bx[1];
    while( c[1] < -b_halfx[1] ) c[1] += bx[1];
    
}
#endif

inline void minVec(double *pos_1, double *pos_2, double *dr, double *bx, double *b_halfx){

        dr[2] = pos_1[2] - pos_2[2];

#ifdef ROTORBOX
            
        double hand;
        
        int box_turns = 0;
        while( dr[2] >= b_halfx[2] ){
            dr[2]     -= bx[2];
            box_turns -= box_turns_per_image;
        }
        while( dr[2] < -1*b_halfx[2] ){
            dr[2]     += bx[2];  
            box_turns += box_turns_per_image;
        }
        box_turns = (box_turns + 444) % 4;
        dr[0] = pos_1[0];
        dr[1] = pos_1[1];
        switch(box_turns){
            case 3:
                hand  = dr[1];
                dr[1] = -1*dr[0];
                dr[0] = hand;
            break;
            case 2:
                dr[0] *= -1;
                dr[1] *= -1;
            break;
            case 1:
                hand  = dr[1];
                dr[1] = dr[0];
                dr[0] = -1*hand;
            break;
            default:
            break;
        }
        
        dr[0] -= pos_2[0];
        dr[1] -= pos_2[1];
        
#else
        while( dr[2] >= b_halfx[2] ){
            dr[2] -= bx[2];
        }
        while( dr[2] < -1*b_halfx[2] ){
            dr[2] += bx[2];
        }
        dr[0] = pos_1[0]- pos_2[0];
        dr[1] = pos_1[1]- pos_2[1];
  
#endif
  
  
        while( dr[0] >= b_halfx[0] ){
            dr[0] -= bx[0];
        }
        while( dr[0] < -1*b_halfx[0] ){
            dr[0] += bx[0];
        }
        while( dr[1] >=  b_halfx[1] ){
            dr[1] -= bx[1];
        }
        while( dr[1] < -1*b_halfx[1] ){
            dr[1] += bx[1];
        }
}


inline int cellIndex( double *R, Box *box ){

    int cellIndex;
    
#define DEBUG_CELLINDEX
#ifdef  DEBUG_CELLINDEX
    if( ( R[0] >= box->halfx[0] || R[0] < -1*box->halfx[0] )
     || ( R[1] >= box->halfx[1] || R[1] < -1*box->halfx[1] )
     || ( R[2] >= box->halfx[2] || R[2] < -1*box->halfx[2] ) ){
         
        cerr << "Diverse alarms!!" << endl;
        cerr << " X: " << R[0] << " Y: " << R[1] << " Z: " << R[2] << " Box Half: " << box->halfx[0] << endl;
        exit( 1 );
    }
#endif
    
    cellIndex =  int((R[0]+box->halfx[0]) / box->xCell[0]) +
                 int((R[1]+box->halfx[1]) / box->xCell[1])*box->nxCell[0] +
                 int((R[2]+box->halfx[2]) / box->xCell[2])*box->nxCell[0]*box->nxCell[1];

#ifdef  DEBUG_CELLINDEX
              for( int i = 0; i < 3; i++ )
                 if( int( (R[i]+box->halfx[i]) / box->xCell[i] ) >= box->nxCell[i] 
                  || int( (R[i]+box->halfx[i]) / box->xCell[i] ) < 0 ) {
                   cerr << "Cell error" << " particle coord + half box: " << (R[i]+box->halfx[i])/(box->xCell[i]) << " is bigger than: " << 
                   box->nxCell[i] << " or smaller than a cell: " << box->xCell[i] << endl;  
                   
                    exit(8);
                 }
#endif


    return ( cellIndex );
}

inline int testCellIndex(double *R, Box *box, int cell_i){
    int check_i;
    
    check_i = cellIndex(R, box);
    
    if( check_i != cell_i ){
        cerr << "Cell index check failed" << endl;
        cerr << "have: " << cell_i << " should have: " << check_i << endl;
        exit(1);
    }
    return( 1 );
}



inline Cell* insertToCell( Particle *p, HashTable *cellHash, int icell, Box *b ){

    Cell *c;
    c = cellHash->getItemByKey(icell);

// p can not be NULL
    if( c == NULL ){
        c = new Cell(icell, b);
        c->firstParticle = p;
        p->next          = NULL;
        cellHash->insertItem( c );
    }else{
        p->next = c->firstParticle;
        if( c->firstParticle )
           c->firstParticle->prev = p;
        c->firstParticle = p;
    }
    p->cell = c;
    p->prev = NULL;

    return( c );

}

inline void removeFromCell( Particle *p, Cell *c ){

    //snip out of the list of particles
// the particle is in the middle
    //assert(p);
    //assert(c);

    //assert(p->cell == c);
    
// particle in the middle
    if( p->prev ){
        p->prev->next = p->next;
        if( p->next ) //if the next particle is not last particle
            p->next->prev = p->prev;		
    }
    else if( p->next ){ // the particle is first particle
        p->next->prev = NULL;
        c->firstParticle = p->next;
    }else // by elimination, is the only particle of the cell
        c->firstParticle = NULL;
    
    // if( p->next )
    //    p->next->prev = p->prev;

    //check if the cell's list needs to be amended.
    /* if( c->firstParticle == p ){
        c->firstParticle = p->next;
    } */
}

inline Cell *moveBetweenCells( Particle *p, Cell *oldCell, int jcell, HashTable *cellHash, Box *b ){

    Cell *c;

    removeFromCell(p, oldCell);
    c = insertToCell(p, cellHash, jcell, b);

    return( c );
}
#endif
