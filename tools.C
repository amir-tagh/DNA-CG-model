#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "FreeEnergyHS.h"
#include "tools.h"

using namespace std;

//this variable declared static, available globally via tools.h
int        box_turns_per_image, stepsize_box_turns_per_image;
double     zrot, sin_zrot, cos_zrot;


#ifdef ROTORBOX
void init_rotorBox(int input_turns_per_image){
    
    //test and set the static
    box_turns_per_image      = input_turns_per_image;
    
    if( box_turns_per_image < 0 || box_turns_per_image > 3){
        cerr << "number of box turns not 0,1,2 or 3" << endl; 
        cerr << "have: "          << box_turns_per_image << endl; 
        exit(1);
    }
    
    /*take care using angles: Pi is not a rational number*/
    switch( box_turns_per_image ){
        case 1:
            zrot = 0.5*M_PI;
            sin_zrot = 1.;
            cos_zrot = 0.;
        break;
        case 2:
            zrot = M_PI;
            sin_zrot =  0.;
            cos_zrot = -1.;
        break;
        case 3:
            zrot = 1.5*M_PI;
            sin_zrot = -1.;
            cos_zrot =  0.;
        break;
        default:
            zrot = 0.;
            sin_zrot = 0.;
            cos_zrot = 1.;
        break;
    }
    
}
#endif

//function to move the molecule such that the COM
//is at the origin.
void hsc::centreMolecule(){
  
   double          COM[3];
   unsigned int    icell;
   double          ref_p[3], dr[3];
   Particle       *p_1;
   
   ref_p[0] = part[0].R[0];
   ref_p[1] = part[0].R[1];
   ref_p[2] = part[0].R[2];
   
   COM[0] = 0.0;  
   COM[1] = 0.0;  
   COM[2] = 0.0; 
   
   for (int i=0; i<N; i++){
       
        p_1 = &part[i];
       
        minVec(p_1->R, ref_p, dr, box->x, box->halfx);
        ref_p[0] += dr[0];
        ref_p[1] += dr[1];
        ref_p[2] += dr[2];
       
        COM[0] += ref_p[0];
        COM[1] += ref_p[1];
        COM[2] += ref_p[2];
   }
   
   COM[0] /= N;
   COM[1] /= N;
   COM[2] /= N;
 
   for (int i=0; i<N; i++){
        part[i].R[0] -= COM[0];
        part[i].R[1] -= COM[1];
        part[i].R[2] -= COM[2];
        
        intoBox(part[i].R, box);
        icell = cellIndex( part[i].R, box );
        if( icell != part[i].cell->myId ){
            moveBetweenCells(&part[i], part[i].cell, icell, cellHash, box );
        }
    }
    return;
}
    
//function to move the molecule such that the COM
//of the last base-pair is at the centre of the z imaging plane.
void hsc::centre_zPlane(Box *box){
  
   double          COM[3];
   unsigned int    icell;
   double          ref_p[3], dr[3];
   Particle       *p_1;
   
   ref_p[0] = part[N/2].R[0];
   ref_p[1] = part[N/2].R[1];
   ref_p[2] = part[N/2].R[2];
   
   
   COM[0] = 0.0;  
   COM[1] = 0.0;  
   COM[2] = 0.0;  
   
   for (int i=(N/2)-2; i < (N/2)+2; i++){
       
       p_1 = &part[i];
       
       minVec(p_1->R, ref_p, dr, box->x, box->halfx);
       ref_p[0] += dr[0];
       ref_p[1] += dr[1];
       ref_p[2] += dr[2];
        
       COM[0] += ref_p[0];
       COM[1] += ref_p[1];
       COM[2] += ref_p[2];
       
    }
   //COM is of four particles: the "top" bp of the double helix.
   COM[0] *= 0.25;
   COM[1] *= 0.25;
   COM[2] *= 0.25;
   
   //centring in x-y plane, but z-plane needs to go to box->halfx[2] (top of the box).
   COM[2] -= box->halfx[2];
   COM[2] -= 0.0000000001;
   
   
   for (int i=0; i<N; i++){   
       
        part[i].R[0] -= COM[0];
        part[i].R[1] -= COM[1];
        part[i].R[2] -= COM[2];
        
        intoBox(part[i].R, box);
        icell = cellIndex( part[i].R, box );
        if( icell != part[i].cell->myId ){
            moveBetweenCells(&part[i], part[i].cell, icell, cellHash, box );
        } 
   }
     
   return;
}
    
  


