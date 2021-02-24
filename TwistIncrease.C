#include <stdio.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <sstream>
#include <string>
#include "FreeEnergyHS.h"
#include "Tools.h"
#include "LJ.h"
using namespace std;

extern float ran2(long *idum);

struct Matrix {
    double m[9];
};
typedef struct Matrix Matrix; 

  
inline void applyRot(double *m, double *R, double *to){
    
     to[0] = m[0] * R[0] + m[1] * R[1] + m[2] * R[2];
     to[1] = m[3] * R[0] + m[4] * R[1] + m[5] * R[2];
     to[2] = m[6] * R[0] + m[7] * R[1] + m[8] * R[2];      
     
}

Matrix *matrix_bprot(double *COM_axis, int clusterSize) {
         
  Matrix *matrix = new Matrix;     
  
  double  sTheta;
  double  cTheta;
  double  ux,uy,uz;
  double  length;
  
  length = sqrt ((COM_axis[0] * COM_axis[0]) + (COM_axis[1] * COM_axis[1]) + (COM_axis[2] * COM_axis[2])); 
  
  ux = COM_axis[0]/length;
  uy = COM_axis[1]/length;
  uz = COM_axis[2]/length;
  
  //cerr << " ux: " << ux << " uy: " << uy << " uz: " << uz << endl; 
  
  double theta = 2*M_PI/(clusterSize);       
  sTheta = sin(theta);
  cTheta = cos(theta);
  
  matrix->m[0] = cTheta + ux * ux * (1 - cTheta);
  matrix->m[1] = ux * uy * (1 - cTheta) - (uz * sTheta);
  matrix->m[2] = ux * uz * (1 - cTheta) + (uy * sTheta);
  matrix->m[3] = uy * ux * (1 - cTheta) + (uz * sTheta);
  matrix->m[4] = cTheta + uy * uy * (1 - cTheta);
  matrix->m[5] = uy * uz * (1 - cTheta) - (ux * sTheta);
  matrix->m[6] = uz * ux * (1 - cTheta) - (uy * sTheta);
  matrix->m[7] = uz * uy * (1 - cTheta) + (ux * sTheta);
  matrix->m[8] = cTheta + uz * uz * (1 - cTheta); 
  return matrix;
}