#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
//#include "FreeEnergyHS.h"
#include <string.h>
//#include "sak_potential.h"



//#define TINY_LENGTH_VALUE 0.0001
//#define TINY_SIN_VALUE    1e-10

//#define LEN(x)             sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])
//#define DOT(x,y)           (x[0]*y[0]+x[1]*y[1]+x[2]*y[2])
//#define INVLEN(x)          1.0/sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])


 using namespace std;
 
 
 int main(){
     
     double a;
     a = 1.002;
      if (a > 1){
          
          cerr << " wrong " << endl;
    } else {
        cerr << " right " << endl;
    }
     
     
}
 
#if 0
 Particle a, b, c, d;

//#if 0

void hsc::testLJ(){
    
    int    x;
    double e, f, g, h;
     
    
    
    memset(a.R, 0, 3*sizeof(double));
    memset(b.R, 0, 3*sizeof(double));
    memset(c.R, 0, 3*sizeof(double));
    memset(d.R, 0, 3*sizeof(double));
    x=0.05;
    
    a.beadType == 'P';
    b.beadType == 'B';
    c.beadType == 'P';
    d.beadType == 'B';
 
    
    for( int i =0; i < 1000; i++ ){
        
        
        b.R[0] = i * x;
        d.R[0] = i * -x;
        
        e = LJ_pairEnergy(a, b);
        f = LJ_pairEnergy(a, d);
        g = LJ_pairEnergy(c, d);
        h = LJ_pairEnergy(b, d);
           
    
        cerr << i*x << " " << e << " " << f << " " << g << " " << h << endl;
        
    }
    
}
#endif

#if 0
void hsc::testDihedral(){
 
  
  double n1[3],n2[3],m1[3];
  double x,y; 
  double x1[3], x2[3], x3[3], x4[3], theta_ref, theta_calculated_by_amir;
  double a[3], b[3], c[3];
  double len_a, len_b, len_c, b_norm;
  double phi, k, deltaPhi, phi0;
  double e_dihed;
  
 
    for(int i = 0; i <  1000; i++ ){
        
    theta_ref = double(i) * 2 *  M_PI / 1000.0;
        
        
    x1[0] =  1.0;
    x1[1] =  0.0;
    x1[2] =  0.0;
    
    x2[0] = 0.0;
    x2[1] = 0.0;
    x2[2] = 0.0;
    
    x3[0] = 0.0;
    x3[1] = 0.0;
    x3[2] = 1.0;
    
    x4[0] =  1.0 * cos(theta_ref);
    x4[1] =  1.0 * sin(theta_ref);
    x4[2] =  1.0;
    
    
    /*******These are the "bond vectors" */
    a[0] = x2[0] - x1[0];
    a[1] = x2[1] - x1[1];
    a[2] = x2[2] - x1[2];
   
    b[0] = x3[0] - x2[0];
    b[1] = x3[1] - x2[1];
    b[2] = x3[2] - x2[2];
    
    c[0] = x4[0] - x3[0];
    c[1] = x4[1] - x3[1];
    c[2] = x4[2] - x3[2];
    
    /**************Amir to insert trial dihedral calculation here */
  
 
 len_a  = 1.0 / LEN(a);
 len_b  = 1.0 / LEN(b);
 len_c  = 1.0 / LEN(c);

 
 for (int i=0; i<3; i++){
      a[i]  *= len_a;
      b[i]  *= len_b;
      c[i]  *= len_c;
    
  }
  
 
  
 //normal vectors to the planes containing a,b,c;
  vector_product(a, b, n1);
  
  vector_product(b, c, n2);
  
  //orthonormal frame
  vector_product(n1, n2, m1 );
  
 
  
  x = DOT(m1,b);
  y = DOT(n1,n2);
  
  
  //dot product of n1 & n2
  //x = DOT(n1,n2);
  //y = DOT(m1,n2);
  
  
  theta_calculated_by_amir = atan2(x,y); 
  
  if( theta_calculated_by_amir != theta_calculated_by_amir )
        theta_calculated_by_amir = 0.0;
  
  
   
  phi0 =  3.51;
  k     = 27.84;
  
  
  deltaPhi = theta_calculated_by_amir -  phi0;
  e_dihed   = k * ( 1 - cos(deltaPhi));
  
       
  cerr << theta_ref << " " << theta_calculated_by_amir << " " << e_dihed << endl;
  
    }
    
    
    
}

#endif
