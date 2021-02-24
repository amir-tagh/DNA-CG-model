//****************************************
//CLASS WANG LANDAU
//****************************************
#ifndef wangLandau_h
#define wangLandau_h


#include "FreeEnergyHS.h"

class wangLandau {
    
    public:
        
        double     binW_twist,       binW_stretch,     inv_binW_twist,     inv_binW_stretch; 
        double     bin_max_stretch,  bin_min_stretch,  bin_min_twist,      bin_max_twist;
        double     bin_max_shrink,   bin_min_shrink;
        double     centre_min_twist, centre_max_twist, centre_min_stretch, centre_max_stretch;
        double     bin_min_twist_update, bin_min_stretch_update;
        double     flatness_test; 
        double     hist_increment;
        double     stepsize_twist, stepsize_stretch, stepsize_shrink;
        double     **hist_2D;
        int        direction_twist,  direction_stretch;
        int        nBins_stretch, nBins_twist;
        int        n_twist, n_stretch;
        bool       productionStage;
       
        
        wangLandau(int nBins_twist, int nBins_stretch, double binW_twist, double binW_stretch){
         
            direction_twist      = 1;
            direction_stretch    = 1;
            flatness_test        = 0.6;
            hist_increment       = 1.0;
            
            n_twist   = nBins_twist;
            n_stretch = nBins_stretch;
            
            productionStage        = false;
            
            stepsize_twist         = 0.25;
            stepsize_stretch       = 0.028;
            
            bin_min_twist          =  -binW_twist*0.5;
            bin_min_twist_update   =  -binW_twist*0.5;
            bin_max_twist          =   nBins_twist*binW_twist - (binW_twist*0.5);
               
            centre_min_twist       =  bin_min_twist + (binW_twist*0.5);
            centre_max_twist       =  bin_max_twist - (binW_twist*0.5);
        
            bin_min_stretch        = -binW_stretch*0.5;
            bin_min_stretch_update = -binW_stretch*0.5;
            bin_max_stretch        =  nBins_stretch*binW_stretch - (binW_stretch*0.5);
                
            centre_min_stretch     =  bin_min_stretch + (binW_stretch*0.5);
            centre_max_stretch     =  bin_max_stretch - (binW_stretch*0.5);
            
            bin_min_shrink         = -binW_stretch*0.5;
            bin_max_shrink         =  nBins_stretch*binW_stretch - (binW_stretch*0.5);
        
            inv_binW_twist         = 1./binW_twist;
            inv_binW_stretch       = 1./binW_stretch;
        }
        
        
        double **create_hist(int nx, int ny){
                
            nBins_twist   = nx; 
            nBins_stretch = ny;
            hist_2D = new double*[nBins_stretch];
            for(int i = 0; i < nBins_stretch; i++){
                    hist_2D[i] = new double[nBins_twist];
            }
            zero_hist();
            return(hist_2D);
        }
        
        
        /*initialize 2darray to zero*/
        void zero_hist(){
            int i, j;
            for(i = 0; i < nBins_stretch; i++){
                for(j = 0; j < nBins_twist; j++){
                    hist_2D[i][j] = 0;
                }
            }
        } 
        
        
        /*get the bin index*/    
        void getBin_2dhist(double twist, double stretch, int *pt_i, int *pt_j){
    
            int    i, j;
            
            i = (int)((twist - bin_min_twist)*inv_binW_twist);
                if ( i < 0 ){
                    *pt_i = 0;
                }else if( i >= nBins_twist ){
                    *pt_i = (nBins_twist - 1);
                }else{
                    *pt_i = i;
                }
        
            j = (int)((stretch - bin_min_stretch)*inv_binW_stretch);
                if ( j < 0){
                    *pt_j = 0;
                }else if( j >= nBins_stretch ){
                    *pt_j  = (nBins_stretch - 1);
                }else{
                    *pt_j = j;
                }
         }
         
         
        double **Normalized_hist(double **hist_2D){
            
            double sum;
            int i, j;
            sum = 0;
            for( i = 0; i < nBins_stretch; i++){
                for( j = 0; j < nBins_twist; j++){
                    sum += hist_2D[i][j];
                }
            }
            
            for( i = 0; i < nBins_stretch; i++){
                for( j = 0; j < nBins_twist; j++){
                    hist_2D[i][j] = hist_2D[i][j]/sum;
                }
            }
            return(hist_2D);
        } 
       
        
        
        double Min_hist(double **hist_2D){
            
            double min, sum;
            int i, j;
            
            min = hist_2D[0][0];
            
            /*normalize the histogram*/
            sum = 0;
            for(i = 0; i < nBins_stretch; i++){
                for(j = 0; j < nBins_twist; j++){
                    sum += hist_2D[i][j];
                }
            }
            
            for(i = 0; i < nBins_stretch; i++){
                for(j = 0; j < nBins_twist; j++){
                    hist_2D[i][j] = hist_2D[i][j]/sum;
                }
            }
            /*get the minimum 0f the histogram*/
            for(i = 0; i < nBins_stretch; i++){
                for(j = 0; j < nBins_twist; j++){
                    if(hist_2D[i][j] < min)
                        min = hist_2D[i][j];
                }
            }
            return(min);
        }
        
        /*find the mean of normalized 2d array*/
        double getMean(double **hist_2D){
        
            double sum, mean;
            int    num_of_elements;
            int    i, j;
                        
            sum = 0;
            for(i = 0; i < nBins_stretch; i++){
                for(j = 0; j < nBins_twist; j++){
                    sum += hist_2D[i][j];
                }
            }
            /*normalize the histogram*/
            for(i = 0; i < nBins_stretch; i++){
                for(j = 0; j < nBins_twist; j++){
                    hist_2D[i][j] = hist_2D[i][j]/sum;
                }
            }
            /*take the sum of the normalized histogram*/
            sum = 0;
            for(i = 0; i < nBins_stretch; i++){
                for(j = 0; j < nBins_twist; j++){
                    sum += hist_2D[i][j];
                }
            }
            
            num_of_elements = (nBins_stretch*nBins_twist);
            mean = sum / num_of_elements;
    
            return(mean);
        }
    
        /*check the flatness of the histogram*/
        /*min(H(E)) > flatness_test*mean  source: J Comp. Chem. 32,816-821, 2011. */ 
        bool checkFlatness(double min, double mean){
            
                if( min > (flatness_test * mean)){
                    cout << " Histogram is flat! " << endl;
                    return true;
                }else{
                    return false;
                }
        }
        
        
        bool histogram_finished(double **hist_2D){
            
            double sum;
            double max_energy, min_energy;
            double min, max;
            double delta_PMF;
            int i, j;
            sum = 0;
            
            for(i = 0; i < nBins_stretch; i++){
                for(j = 0; j < nBins_twist; j++){
                    sum += hist_2D[i][j];
                }
            }
            
            /*normalize the histogram*/
            for(i = 0; i < nBins_stretch; i++){
                for(j = 0; j < nBins_twist; j++){
                    hist_2D[i][j] = hist_2D[i][j] / sum;
                }
            }
            
            min = hist_2D[0][0];
            max = hist_2D[0][0];
            for(i = 0; i < nBins_stretch; i++){
                for(j = 0; j < nBins_twist; j++){
                    if(hist_2D[i][j] < min)
                        min = hist_2D[i][j];
                    if(hist_2D[i][j] > max)
                        max = hist_2D[i][j];
                }
            }
            
            
            min_energy = -log(min);
            max_energy = -log(max);
            delta_PMF  = max_energy - min_energy;
            
            /*if the biggest spread in the PMF is less than 0.1 kBT then we can stop
               the nonequilibrium part of the calculation and start collecting final data*/
            cout << " Testing if histogram is complete, delta PMF is: " << delta_PMF << endl;
            
            if(delta_PMF < 0.1){
                return true;
            }else{
                return false;
            }
        }
        
        /*set the histogram to zero if flatness checks out*/
        double **make_zero(double **count_hist){
            
            int i, j;
            for( i = 0; i < nBins_stretch; i++){
                for( j = 0; j < nBins_twist; j++){
                    count_hist[i][j] = 0;
                }
            }
            return(count_hist);
        }
        
};

#endif





























