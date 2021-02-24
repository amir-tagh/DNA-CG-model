#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

class wl2Dhist {
    
    public:
        double   **hist_2D;
        double     nhist_increment;
        double     flatness_criteria;
        double     hist_increment;
        int        nBins_stretch, nBins_twist;
        
        //constructor
        wl2Dhist(int nx, int ny){
            
            nBins_twist   = nx;
            nBins_stretch = ny;
        }
        
           /*function to make 2d arrays*/
            double **create_hist(int nx, int ny){
            
                int i,j;
                hist_2D = new double*[nx];
                for(i=0; i<nx; i++){
                    hist_2D[i] = new double[ny];
                }
                
                
            }
};


int main(){
    
    
    wl2Dhist *p1, *p2;
    
    p1 = new wl2Dhist(5,11);
    p2 = new wl2Dhist(8,24);
    
    double **hist;
    hist = p1->create_hist(5,11);
    cout << " array elem: " << hist[0][0] << endl;
    
    
    //cout << " x bins: " << p1->nBins_twist << " y bins: " << p1->nBins_stretch << endl;
    //cout << " x bins: " << p2->nBins_twist << " y bins: " << p2->nBins_stretch << endl;
    
    
}