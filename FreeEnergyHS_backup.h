/*************************************************************************
 *  Monte Carlo code for Coarse Grained Model of DNA,                    *  
 * Linear as well as circular form of B-DNA                              *                           
 * DNA structtural parameters are based on the standard B_DNA model      *
 *                      
 *                                                                       *
 * -n: start a new Configuration (hcp)                                   *
 *                                                                       *
 * -f<filename>: start from a configuration file (output of an earlier   *
 *               simulation). Filename contains N, rho, box              *
 *               dimensions and center of mass vectors of all particles  *
 * -e<filename>: take centers of attractive wells from file, same        *
 *               structure as with -f                                    *
 *               if -e is not used, centers are put on hcp lattice       *
 *                                                                       *
 * -p<parameterfile>: parameterfile contains:                            *
 *                    number of particles, box dimensions,               *
 *                    #Blocks equilibration, #Blocks MC,                 *
 *                    #evaluations during MC,#evaluations during equil.  *
 *                    #Snapshot Configurations Eq and MC,                *
 *                    maximal displacement, a seed,                      *
 *                    a coupling strength epsilon, a cutoff radius rc    *
 *                                                                       *
 *                                                                       *
 *************************************************************************/
#ifndef HAVE_FEHS
#define HAVE_FEHS

#include <iostream>
#include <fstream>
#include <ostream>
#include <cstdlib>
#include <math.h>
#include "hashTable.h"

using namespace std;

#define LINEAR_FREE_ENDS  

#define RESET_ENERGY_TOTALS_EVERY 100 /*testfix numerical drift of total potential energy*/
#define CHECK_POINT 10000              /*trajectory collection in case of crash*/


/*defined values for wang landau*/
#define NUM_HIST             10000
#define LN_F_LIMIT           1e-8                        /*limit value of the adjustment factor */
#define INITIAL_LN_F         2.71828182845               /*initial value of the adjustment factor natural log of e*/
#define WEIGHT_FACTOR_LIMIT  1e-9                       /*weight factor limit for wangLandau method*/
#define INCREASED TWIST      0.25                       /*increased value of imposed twist on DNA*/    


class Particle;
class Box;
class Adaptivestepsize;
class AdaptiveStepRigid;
class wangLandau;
class wl2Dhist;




class hsc {

 public:

  bool newSetUp;
  bool readWells;
  const char *ParameterFile, *WellFile; // Input files
  
  std::string configFileName;
  std::string outputFileName;
  std::string restartFileName;
  std::string rotationFileName;
  std::string normalFileName;
  std::string rejectedFileName;
  std::string normalDNAFileName;
  std::string EquilFileName;
  std::string checkpointFileName;
  std::string toporotFileName;
  std::string extshrinkFileName;
  std::string ExtensionFileName;
  std::string ShrinkFileName;
  std::string TopologyFileName;
  std::string Topo_pot_ener_plusFileName;
  std::string Topo_pot_ener_minusFileName;
  std::string topo_ener_infoFileName;
  std::string Ext_energy_accFileName;
  std::string Ext_energy_rejFileName;
  std::string Shr_energy_accFileName;
  std::string Shr_energy_rejFileName;
  
  void usage();
  void run();

  Particle *part;
  Cell     *cells;           /*cells for overlap detection*/
  Box      *box;
  Adaptivestepsize *adaptive;
  AdaptiveStepRigid *adaptive_rigid;
  wangLandau *wL;
  wl2Dhist   *wLhist;

  HashTable *cellHash;

  void testCellIndices();
  void brute_debug_cells();
  int  reportNebs(int *nebs, int nnebs, int id1, int id2);
  int  testNeb(int *nebs, int nnebs, int id, int id2);

  /* Simulation Data*/
  long int transAccept;
  long int transMoves;
  
  /*topology move*/
  long int topoAccept_positive;  
  long int topoAccept_negative;
  long int topoReject_positive;
  long int topoReject_negative;
  long int topologyMoves;  
  long int topology_accept;
  long int topology_reject;

  /*rigid move*/
  long int rigidCount;
  long int rigidMoves;
  long int rigidAccept;
  long int rigidReject;
  
  /*extension and shrink move*/
  long int ShrinkMove;
  long int ShrinkAcc;
  long int ShrinkRej;
  long int ExtensionMove;
  long int ExtensionAcc;
  long int ExtensionRej;
  long int LineCount; //debug parameter
  
  /*MC moves to calculate trasition matrix entries.*/
  double   extension_step, shrink_step; 
  
  
  //debug counter
  long int new_cell_counter_up;
  long int new_cell_counter_down;
  
  double   maxDisplace;     /*translation move size*/
  double   thetaRigid;     /*rigid move size in radian*/
  
  
  int  runCount;
  int  measureCount;
  long seed;          //seed for ran3
  long *iran;         //pointer to seed

  int eqSteps;
  int eqEvalSteps;      // # trajectory averages taken during equilibration.
  int eqSnapSteps;      // # configuration snapshots taken during equil.
  int mcSteps;
  int mcEvalSteps;      /*measurements for averages*/
  int mcSnapSteps;      /*snapshots during mc*/
  //int steps;            //global variable to collect data   


  double bx;         /*box dimensions*/
  int N;            /*particle number*/
  double Ninv;
  double halfN;

  double ePotTot, eDihedTot, eBondTot, eLJTot, eSAKTot, phi;
  void   writeSystemEnergy();
  double resetSystemEnergy();
  
  

 //LJ ENERGY
  double LJ_pairEnergy(Particle *p, Particle *q);
  double LJ_pairEnergy_dbg(Particle *p, Particle *q);
  double LJ_totalEnergy();
  double LJ_partEnergy_skipC( Particle *p, Cell *c );
  double LJ_partEnergy( Particle *p, Cell *c );
  double LJ_totalEnergy_wClus(); //debug function
  double LJ_topology(int clusterStart_1, int clusterSize);
  double LJ_Debug_totalEnergy();
  double LJ_pairEnergy_debug(int clusterStart_1, int clusterSize); //LJ rigid debug function//

 //BONDING information
  double kBond, r0Bond;
  int    windingNumber;
  double bond_totalEnergy();
  double bond_bptotalEnergy();
  double bond_totalEnergy_noBaseStack();
  double bond_bprotate_totalEnergy();
  double bond_onePartEnergy(Particle *p);
  double bond_bprotateEnergy(Particle *p,Particle *q);
  double bond_bprotate(Particle *p);
  double bond_topology(int clusterStart_1, int clusterStart_2, int clusterEnd_1, int clusterEnd_2);
  double rigid_bond_onePartEnergy_backward(Particle *p);
  double rigid_bond_onePartEnergy_forward(Particle *p);

  //debug bonding
  double print_bond_totalEnergy(char *fname);
  double print_lj_totalEnergy(char *fname);
  double bond_total();

  /*DNA structural parameters*/
  double  R, h,  rBase, theta; 
  double x, y, z;
  
  //DNA topology: 1 for linear(topologically circular) and 0 for circular DNA
  int initialize_chain_as_linear;

  //flag for making free fluctuating linear form of DNA
  int open_chain;
  
  //flag for non-integer rotation of DNA in topology move
  double box_rotation_radians;

  double epsilon, eps2, eps3; //coupling strength
  double rcut, rcut2, rcut3, rcut2Inv, Dcut;  //cutoff of range of attraction

  /*dihedrals*/
  double dihed_onePartEnergy(Particle *p);
  double dihed_forward(Particle *p, Particle *q);
  double dihed_backward(Particle *p, Particle *q);
  double forward_backward_sum(Particle *p, Particle *q);
  double dihed_boundaryEnergy_fromAbove(Particle *p);
  double dihedEnergy(Particle *p, Particle *q, Particle *r, Particle *s);
  double espresso_dihedEnergy(Particle *p, Particle *q, Particle *r, Particle *s);
  double dihed_topology(int clusterStart_1, int clusterStart_2, int clusterEnd_1, int clusterEnd_2);
  double rigidbody_dihed(int clusterStart_1, int clusterStart_2, int clusterEnd_1, int clusterEnd_2, double *edihed1, double *edihed2, double *edihed3, double *edihed4);
  double dihed_forward_topo(Particle *p, Particle *q);

  //function for debug LJ
  void   NeighbourList(int clusterStart_1, int clusterSize);

  double totalDihedEnergy();

  //test harness
  //void testLJ();
  //void testDihedral();
  //void testDihedral_wikipedia();


  // observables
  // potential energy with respect to attractive wells
  double meanPhi;
  double meanDistanceFromWell;
  double meanRatioOutside;

  /* tabulated potential*/
  double SAK_pairEnergy(Particle *p, Particle *q);
  double SAK_totalEnergy();

  //routines
  //set up
  void setUpNew();
  void setUpFromFile();
  void setUpWells();
  void initialize();
  void readInWells();

  void packABCStretchAndExpand();

  void setupList();
  void setupListW();
  void setUpNeighbours();
  void setUpNeighboursW();

  /* simulation*/
  void equil();
  void MC();

  void particleMove();
  void ClusterTransitionMove();
  void relocMove();
  void centreMolecule();
  void centre_zPlane(Box *box);
  void topology_rotation(std::ofstream *of1, std::ofstream *of2, std::ofstream *of3);
  void RigidBodyRotation();
  void shrink(std::ofstream *of4, std::ofstream *of5, std::ofstream *of6);
  void extension(std::ofstream *of7, std::ofstream *of8, std::ofstream *of9);

  double computePhiOfConfig();
  double computePhiOfDist(double rsq);    // takes distance squared as input
  double computeFRef();
  double computeAverageDistanceFromWell();
  double computeRatioOutside();

  void displaceParticle(double newP[3], double oldP[3]);
  //void minVec(double *pos_1, double *pos_2, double *dr);
  void relocateParticle(Particle &p);
  bool checkOverlapAll(Particle &p);
  bool checkOverlapTwo(Particle &p1, Particle &p2);

  /*input & output*/
  void writeOutConfig(const char  name[16]);
  void writeOutSimData(const char name[16]);
  void readXYZ(const char* name);
  void readXYZ(std::string name);
  void writeXYZ(char name[16]);
  void writeXYZ_frame(std::ofstream *of);
  void writeOutWells(char name[16]);
  void writeOutBinding(char name[16]);
  void writeOutBindingToWells(char name[16]);
  
#ifdef ROTORBOX
  void changeBoxRotations(int new_box_turns_per_image);
#endif
  
};


/*wang landau algorithm*/
class wl2Dhist {
    
    public:
        double   **hist_2D;
        double     nhist_increment;
        double     flatness_test;
        double     hist_increment;
        int        nBins_stretch, nBins_twist;
        
        //constructor
        wl2Dhist(int nx, int ny){
            
            int i;
                
            nBins_twist   = nx;
            nBins_stretch = ny;
            hist_2D = new double*[nx];
            for(i = 0; i < nBins_twist; i++){
                    hist_2D[i] = new double[nBins_stretch];
            }
            zero_hist();
        }
        
        
        /*initialize 2darray to zero*/
        void zero_hist(){
            int i, j;
            for(i = 0; i < nBins_twist; i++){
                for(j = 0; j < nBins_stretch; j++){
                    hist_2D[i][j] = 0;
                }
            }
        } 
       
        
        
        double Min_hist(){
            
            double min; 
            double sum;
            int i,j;
            
            min = hist_2D[0][0];
            
            /*normalize the histogram*/
            sum = 0;
            for(i = 0; i < nBins_twist; i++){
                for(j = 0; j < nBins_stretch; j++){
                    sum += hist_2D[i][j];
                }
            }
            
            for(i = 0; i < nBins_twist; i++){
                for(j = 0; j < nBins_stretch; j++){
                    hist_2D[i][j] = hist_2D[i][j]/sum;
                }
            }
            /*get the minimum 0f the histogram*/
            for(i = 0; i < nBins_twist; i++){
                for(j = 0; j < nBins_stretch; j++){
                    if(hist_2D[i][j] < min)
                        min = hist_2D[i][j];
                }
            }
            return(min);
        }
        
        /*find the mean of normalized 2d array*/
        double getMean(){
        
            double sum, mean;
            int    num_of_elements;
            int    i, j;
                        
            sum = 0;
            for(i = 0; i < nBins_twist; i++){
                for(j = 0; j < nBins_stretch; j++){
                    sum += hist_2D[i][j];
                }
            }
            /*normalize the histogram*/
            for(i = 0; i < nBins_twist; i++){
                for(j = 0; j < nBins_stretch; j++){
                    hist_2D[i][j] = hist_2D[i][j]/sum;
                }
            }
            /*take the sum of the normalized histogram*/
            sum = 0;
            for(i = 0; i < nBins_twist; i++){
                for(j = 0; j < nBins_stretch; j++){
                    sum += hist_2D[i][j];
                }
            }
            
            num_of_elements = (nBins_stretch*nBins_twist);
            mean = sum / num_of_elements;
    
            return(mean);
        }
    
        /*check the flatness of the histogram*/
        /*min(H(E)) > flatness_test*mean  source: J Comp. Chem. 32,816-821, 2011. */ 
        bool checkFlatness(double mean, double min, double flatness_test){
            
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
            
            for(i = 0; i < nBins_twist; i++){
                for(j = 0; j < nBins_stretch; j++){
                    sum += hist_2D[i][j];
                }
            }
            
            /*normalize the histogram*/
            for(i = 0; i < nBins_twist; i++){
                for(j = 0; j < nBins_stretch; j++){
                    hist_2D[i][j] = hist_2D[i][j] / sum;
                }
            }
            
            min = hist_2D[0][0];
            max = hist_2D[0][0];
            for(i = 0; i < nBins_twist; i++){
                for(j = 0; j < nBins_stretch; j++){
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
        
};
         


class wangLandau {
    
    public:
        double     binW_twist,       binW_stretch,     inv_binW_twist,     inv_binW_stretch; 
        double     bin_max_stretch,  bin_min_stretch,  bin_min_twist,      bin_max_twist;
        double     centre_min_twist, centre_max_twist, centre_min_stretch, centre_max_stretch;
        double     flatness_test; 
        int        initial_twist,    initial_stretch,  directional_twist,  directional_stretch;
        int        hist_increment;
         
        bool       productionStage;
        
     
     wl2Dhist *count_wLhist, *bias_hist;
     
     /*this is a constructor*/
     wangLandau(int nBins_twist, int nBins_stretch, double binW_twist, double binW_stretch){
         
         
        initial_twist        = 0;
        initial_stretch      = 0;
        directional_twist    = 0;
        directional_stretch  = 0;
        flatness_test        = 0.6;
        hist_increment       = 0.05;
            
        productionStage     = false;
               
        centre_min_twist   = -binW_twist*(nBins_twist - 1)*0.5;
        centre_max_twist   =  binW_twist*(nBins_twist - 1)*0.5;
               
        centre_min_stretch = -binW_stretch*(nBins_stretch - 1)*0.5;
        centre_max_stretch =  binW_stretch*(nBins_stretch - 1)*0.5;
        
        
        bin_min_twist      = centre_min_twist - 0.5*binW_twist;
        bin_max_twist      = centre_max_twist + 0.5*binW_twist;
        
        bin_min_stretch    = centre_min_stretch - 0.5*binW_stretch;
        bin_max_stretch    = centre_max_stretch + 0.5*binW_stretch;
        
        inv_binW_twist     = 1./binW_twist;
        inv_binW_stretch   = 1./binW_stretch;
        
        count_wLhist = new wl2Dhist(nBins_twist, nBins_stretch);
        bias_hist    = new wl2Dhist(nBins_twist, nBins_stretch);
        
        
     }
     
     
        /*get the bin index*/    
        void getBin_2dhist(double twist_update, double stretch_update, int nBins_twist, int nBins_stretch, int *pt_i, int *pt_j){
    
            int    i, j;
            
            i = (int)((twist_update - bin_min_twist)*inv_binW_twist);
                if ( i < 0 ){
                    *pt_i = 0;
                }else if( i >= nBins_twist ){
                    *pt_i = (nBins_twist - 1);
                }else{
                    *pt_i = i;
                }
        
            j = (int)((stretch_update - bin_min_stretch)*inv_binW_stretch);
                if ( j < 0){
                    *pt_j = 0;
                }else if( j >= nBins_stretch ){
                    *pt_j  = (nBins_stretch - 1);
                }else{
                    *pt_j = j;
                }
         }
     

};


class Particle {

 public:
  int    myId;
  int    clusterId;
  double R[3];     // center of mass
  char   beadType; // beadType = 'P' or 'B'.

  Particle *next,  *prev;   // next, previous in cell list
  Particle *bNext, *bPrev;  // next, previous bases in DNA chain
  Particle *pNext, *pPrev;  // next, previous phosphates in DNA chain
  Cell     *cell;           // cell with respect to overlap radius


  Particle(){
      next      = 0;
      prev      = 0;
      cell      = 0;
      clusterId = 0;
      bNext   = NULL;
      bPrev   = NULL;
      pNext   = NULL;
      pPrev   = NULL;
  }
};

class Box {

 public:

  double x[3], halfx[3], halfxInv[3];  // dimensions
  double V;      // Volume
  
  double xCell[3]; //cell dimension
  int    nxCell[3]; // number of cells
  int    nCells;
  int    rcut; //cutoff distance for LJ 
  

  void update(){
    for (int i=0; i<3; i++){
	  halfx[i] = 0.5*x[i];
	  halfxInv[i] = 1.0/halfx[i];
    }
    V = x[0]*x[1]*x[2];
    
    for(int i=0; i<3; i++){
         nxCell[i] = int(x[i]/rcut); //number of cells per edge of a cube
        
         xCell[i]  = x[i]/nxCell[i]; //dimension of each cell
    }
     //cout << " number of cells: " << nxCell[0] << nxCell[1] << nxCell[2] << endl;
    nCells = nxCell[0]*nxCell[1]*nxCell[2];

  if ((nxCell[0] < 2)||(nxCell[1] < 2)||(nxCell[2] < 2)){
    cerr << "Error: Box dimension is less than 2" << endl;
    cerr << x[0] << " x box " << endl;
    exit(8);
  }
    //cerr << "Allocating : " << nCells << " cells, total size: " << nCells * sizeof(Cell) << " bytes " << endl;
  }
};
//dynamic step sizing for particle movement
class Adaptivestepsize {
    
public:
    
    long int transMoves;
    long int transAccept;
    long int myAccFind;
    double   maxDisplace;
    
    Adaptivestepsize(){
        
       transMoves  = 0;
       transAccept = 0;
       myAccFind   = 0;
       maxDisplace = 0;
        
    }
    
inline void Adaptive_step_size(){
         
   float(myAccFind) = float(transAccept)/float(transMoves);
   
   
    if( fabs(myAccFind - 0.5) >= 0.01){
        if(  myAccFind != 0
             && myAccFind != 1.0  ){
                if( myAccFind > 0.5 + 0.01 ){
                    maxDisplace *= 1.001;
                }else if (myAccFind < 0.5 - 0.01 ){
                     maxDisplace *= 0.999;
                }
        }
    }
}
};

//dynamic rotation angle for rigid body rotation
class AdaptiveStepRigid {
    
public:
    long int rigidMoves;
    long int rigidAccept;
    long int myAccFind_rigid;
    double   theta_rigid;
    
    AdaptiveStepRigid(){
        
        rigidMoves      = 0;
        rigidAccept     = 0;
        myAccFind_rigid = 0;
        theta_rigid     = 0;
    }
inline void Adaptive_step_rigid(){

   float(myAccFind_rigid) = float(rigidAccept)/float(rigidMoves); 

   if( fabs(myAccFind_rigid - 0.3) >= 0.01){
        if(  myAccFind_rigid != 0
             && myAccFind_rigid != 1.0  ){
                if( myAccFind_rigid > 0.3 + 0.01 ){
                    theta_rigid *= 1.001;
                }else if (myAccFind_rigid < 0.3 - 0.01 ){
                     theta_rigid *= 0.999;
                }
        }
    }
}
};

/*inline int SAK_myPair( Particle *p, int i, int N ){  
  
    if( p->beadType == 'P' ){
      return( N - i - 2);
    }else{
      return( N - i );
    }
  
};

// image point back into box
//currently recangular box only
inline void intoBox(double *dr, Box *box){
    
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
  while( dr[2] >= box->halfx[2] ){
    dr[2] -= box->x[2];
  }
  while( dr[2] < -1*box->halfx[2] ){
    dr[2] += box->x[2];
  } 
}*/



/*class Cell {

 public:

  Particle *firstParticle;
  Cell     *neighbours[26];      // neighbouring cells 
  Cell(void){
    firstParticle = 0;
    for (int i=0; i<26; i++){
      neighbours[i] = 0;
    }
  }

};
*/




#endif











