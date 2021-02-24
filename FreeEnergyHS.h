/*************************************************************************
 *  Monte Carlo code for Coarse Grained Model of DNA,                    *
 * Linear as well as circular form of B-DNA                              *
 * DNA structtural parameters are based on the standard B_DNA model      *
 *                                                                       * 
 *                                                                       *
 * -n: start a new Configuration (hcp)                                   *
 *                                                                       *
 * -f<filename>: start from a configuration file (output of an earlier   *
 *               simulation). Filename contains N, box                   *
 *               dimensions and coordinates of all particles             *
 *                                                                       *
 *                                                                       *
 *                                                                       *
 * -p<parameterfile>: parameterfile contains:                            *
 *                    number of particles, box dimensions,           q    *
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


#define RESET_ENERGY_TOTALS_EVERY 100000 /*test/fix numerical drift of total potential energy*/
#define CHECK_POINT 10000               /*trajectory collection in case of crash*/


/*defined values for wang landau*/
#define NUM_HIST             10000
#define LN_F_LIMIT           1e-8                        /*limit value of the adjustment factor */
#define INITIAL_LN_F         2.71828182845               /*initial value of the adjustment factor natural log of e*/
#define WEIGHT_FACTOR_LIMIT  1e-9                       /*weight factor limit for wangLandau method*/
#define INCREASED TWIST      0.25    /*increased value of imposed twist on DNA*/
#define WLCHECK              2000000



class Particle;
class Box;
class Adaptivestepsize;
class AdaptiveStepRigid;
class wangLandau;
class wl2Dhist;


/*adaptive extension to soften the bond and dihedral thresholds to get the second flat region of force-extension curve*/
class AdaptiveExtension{
    
public:
    
    double   length_init;
    double   length_ext;
    double   Relative_length; 
    double   force_imposed;
    double   BondThreshold;
    double   BonDihed_BPBP_Rcut, BonDihed_PBPB_Rcut;
    
    AdaptiveExtension(){
        
        BondThreshold      = 0;
        BonDihed_BPBP_Rcut = 0;
        BonDihed_PBPB_Rcut = 0;
        length_ext         = 0;
        length_init        = 0;
        length_ext         = 0;
        force_imposed      = 0;
        Relative_length    = 0;
    }
    
  inline void Adaptive_Extension(){
    
    float(Relative_length) = float(length_ext)/float(length_init);
    
    if( force_imposed > 60. && force_imposed <= 74. ){ 
        if( fabs(Relative_length + 0.001) > 1.01 ){
            BonDihed_BPBP_Rcut *= 0.988;
            BonDihed_PBPB_Rcut *= 0.988;
            BondThreshold      *= 0.988;
        
        }else if ( fabs(Relative_length + 0.001) > 1.72 ){
            BonDihed_BPBP_Rcut = BonDihed_BPBP_Rcut;
            BonDihed_PBPB_Rcut = BonDihed_PBPB_Rcut;
            BondThreshold      = BondThreshold;
        }
    }
  } 
};


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
  std::string Ext_forceFileName;
  std::string Shr_energy_accFileName;
  std::string Shr_energy_rejFileName;

  void usage();
  void run();

  Particle *part;
  Cell     *cells;           /*cells for overlap detection*/
  Box      *box;
  Adaptivestepsize *adaptive;
  AdaptiveStepRigid *adaptive_rigid;
  AdaptiveExtension *adaptive_ext;

  /*wang landau*/
  wangLandau *wL;
  double    **count_hist, **bias_hist;     /*make two main histograms for wang landau*/
  long int    count;

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
  long int intercal_reject;
  long int intercal_acc;

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
  long int stack_break_count;
  long int BONDIHED_PBPB_RCUT_count;
  long int BONDIHED_BPBP_RCUT_count;
  long int bp_extension_acc;
  long int bp_extension_rej;
  
  /*Intercalative move*/
  double intercal_r0;   //(2.0 * 3.2)  //pauls' paper: resting bp stack is about twice usual when intercalated
  double intercal_epsilon; //(-18.0) paul's paper
  double r0_pp_intercal;
  double intercal_delta;
  double sStart_x;
  double ETA;
  double Delta;         //free energy penalty for neighbour intercalation
  double Epsilon_one;   //free energy penalty from B to S-DNA(3.2)
  double sStretch_lim;
  double U_max;
  double lambda;
  /*MC moves to calculate trasition matrix entries.*/
  double   extension_step, shrink_step;

  //debug counter
  long int new_cell_counter_up;
  long int new_cell_counter_down;

  double   maxDisplace;     /*translation move size*/
  double   thetaRigid;     /*rigid move size in radian*/
  
#ifdef INTERCAL
  double BondThreshold;
  double BonDihed_BPBP_Rcut;
  double BonDihed_PBPB_Rcut;
#endif


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
  double bond_s_ic(Particle *p, Particle *q, double r, double k_b, double r0 );

  //debug bonding
  double print_bond_totalEnergy(char *fname);
  double print_lj_totalEnergy(char *fname);
  double bond_total();
  void   stack_break();
  void   dihed_break();


  //DNA topology: 1 for linear(topologically circular) and 0 for circular DNA
  int initialize_chain_as_linear;

  //flag for making free fluctuating linear form of DNA
  int open_chain;

#ifdef ROTORBOX
  float linking_number;
#else
  int   linking_number;
#endif

  double epsilon, eps2, eps3; //coupling strength
  double rcut, rcut2, rcut3, rcut2Inv, Dcut;  //cutoff of range of attraction
  
  ///////debug variable////////
  double part_init;

//flags for force-distance calculation
  int    use_constant_tension;
  double tension_pn;
  double mu_interCalator;
  int    intercal_count;
  int    sdna_count;
  
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
  void topology_rotation();
  void RigidBodyRotation();
  void shrink(std::ofstream *of4, std::ofstream *of5, std::ofstream *of6);
  void extension(std::ofstream *of7);
  void extension_bpStep();
  void intercalate_bpMove();
  void SDNA_baseMove();

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


class Particle {

 public:
  int       myId;
  int       clusterId;
  double    R[3];     // center of mass
  char      beadType; // beadType = 'P' or 'B'.
  short int interCal; // intercalation status
  short int S;        // S-DNA status

  Particle *next,  *prev;   // next, previous in cell list
  Particle *bNext, *bPrev;  // next, previous bases in DNA chain
  Particle *pNext, *pPrev;  // next, previous phosphates in DNA chain
  Cell     *cell;           // cell with respect to overlap radius

  Particle(){
      next      = 0;
      prev      = 0;
      cell      = 0;
      clusterId = 0;
      interCal  = 0;
      S         = 0;
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
     
    nCells = nxCell[0]*nxCell[1]*nxCell[2];

  if ((nxCell[0] < 2)||(nxCell[1] < 2)||(nxCell[2] < 2)){
    cerr << "Error: Box dimension is less than 2" << endl;
    cerr << x[0] << " x box " << endl;
    exit(8);
  }
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


#endif
