#include <stdio.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include "FreeEnergyHS.h"

long int transMoves;
long int transAccept;
long int myAccFind;
double   maxDisplace;



void hsc::Adaptive_step_size(){
         
          if( transMoves > 1000 ){
                myAccFind = double(transAccept)/double(transMoves);
                if( fabs(myAccFind - 0.5) >= 0.01){
                if( myAccFind == myAccFind
                 && myAccFind != 0
                 && myAccFind != 1.0  ){
                    if( myAccFind > 0.5 + 0.01 ){
                     maxDisplace *= 1.001;
                    }else if (myAccFind < 0.5 - 0.01 ){
                     maxDisplace *= 0.999;
                    }
                    transAccept = 0;
                    transMoves  = 0;
                }
           }
       }
       cout << " step size: " << maxDisplace << " acc rate: " << myAccFind << endl;
}
    
 