#include <stdio.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
using namespace std;


void writeXYZ_frame( ofstream *myfile, int i ){

  *myfile << "coordinates? " << i << endl;

  return();
}


int main(){
  
  /*declaratins here */
  int   i;
  std::ofstream myfile;

  /* commands: open the file*/
  myfile.open ("myTrajName.xyz", ios::out);

  for( i = 0; i < 100; i++){
  
  /* commands: al our function (defined above) with the arguments (&myFile, i) */
    writeXYZ_frame( &myfile, i );
    
  } 
  
  myfile.close();

  return( 0 );
}

