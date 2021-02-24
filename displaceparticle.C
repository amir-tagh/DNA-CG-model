#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "FreeEnergyHS.h"
#include "tools.h"
#define TINY_LENGTH_VALUE 0.0001
#define TINY_SIN_VALUE    1e-10

class Adaptivestepsize;

using namespace std;

extern float ran2(long *idum);


void hsc::displaceParticle(double newP[3], double oldP[3]) {

  double displace;
  for (int i=0; i<3; i++){
 
   //update the max displacement
    maxDisplace = adaptive->maxDisplace;
    displace = (ran2(iran) - 0.5)*maxDisplace;
    newP[i]  = oldP[i] + displace;
    
    //old imaging is both complex and inefficient
    //newP[i]  = -int((oldP[i]+displace)*box->halfxInv[i])*box->halfx[i] 
                       // + fmod((oldP[i]+displace),box->halfx[i]);
  }
   intoBox(newP, box);
}



// checks for overlap of particle pNr with all particles from list

bool hsc::checkOverlapAll(Particle &p) {

  bool overlap = false;
  Cell *curCell;
  Particle *curPart;

  curCell = p.cell;
  curPart = curCell->firstParticle;

  // Zelle um i
  while (curPart){
    overlap = checkOverlapTwo(p, *curPart);
    if (overlap) {
      if (curPart == &p) overlap = false;
      else break;
    }
    curPart = curPart->next;
  }

  // Nachbarzellen
  if (!overlap){
    for(int j=0; j<26; j++){

      Cell *curCell = cellHash->getItemByKey( p.cell->neighbours[j] );
      if( curCell == NULL ) continue;

      curPart = curCell->firstParticle;
      while (curPart){
	         overlap = checkOverlapTwo(p, *curPart);
	         if (overlap) {
	                if (curPart == &p) overlap = false;
	                else break;
             }
	         curPart = curPart->next;
      }
      if (overlap) break;
    }
  }

  return overlap;
}



bool hsc::checkOverlapTwo(Particle &p1, Particle &p2) {

    double rij[3], rijsq;

    for (int j=0; j<3; j++){
	rij[j] = p1.R[j] - p2.R[j];
	if (fabs(rij[j]) > box->halfx[j])
	    rij[j] -= copysign(box->x[j], rij[j]);
    }
    rijsq = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];

    if (rijsq < 1.0) return true;
    else return false;
}






