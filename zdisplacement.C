#include <stdio.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <sstream>
#include <string>
#include "FreeEnergyHS.h"
#include "tools.h"
#include "LJ.h"
#include "cell.h"
#include "hashTable.h"


void hsc::z_displacement(double newP[3], double oldP[3]){

 for( int i=0; i<N; i++){ 
  
     double z_displace;
     double ext_shrink = ((0.2*(2*ran2(iran) -1));
  
     newP[i]  = oldP[i] * ext_shrink;
   
  }

   intoBox(newP, box);
  
}