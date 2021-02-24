#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "FreeEnergyHS.h"
#include "cell.h"
#include "hashTable.h"
#include "tools.h"
#include "wangLandau.h"

class Adaptivestepsize;

using namespace std;

extern float ran2(long *idum);

//** external global variable for linear or circular DNA
//int initialize_chain_as_linear = 1;


void hsc::initialize(){

  if( rcut < 12.0 ){
    cerr << "interaction cutoff is lower than maximum LJ cutoff in PAK model" << endl;
    exit( 8 );
  }
  if( N % 2 != 0 ){
    cerr << "Require even number of particles (1 Base per Phosphate) in the PAK model" << endl;
    exit( 8 );
  }
  if( N % 4 != 0 ){
    cerr << "Probably better to have a multiple of 4 particles: 1bp = 2 bases+2 phosphates." << endl;
    exit( 8 );
  }
  if( N < 4 ){
    cerr << "Probably better to have at least 8 particles: 1bp = 2 bases+2 phosphates." << endl;
    exit( 8 );
  }

  Ninv  = 1.0/float(N);
  rcut2 = rcut*rcut;
  rcut3 = rcut2*rcut;
  rcut2Inv = 1.0/rcut2;
  Dcut = 2.0*rcut;
  eps2 = epsilon*epsilon;
  eps3 = eps2*epsilon;

  //set counters
  transAccept = 0;
  transMoves  = 0;

  //counter for tolpology move
  topologyMoves       = 0;
  topoAccept_positive = 0;
  topoAccept_negative = 0;
  topoReject_positive = 0;
  topoReject_negative = 0;
  topology_accept     = 0;
  topology_reject     = 0;
  intercal_reject     = 0;
  intercal_acc        = 0;


  //counter for rigid move
  rigidMoves  = 0;
  rigidAccept = 0;
  rigidReject = 0;

  //counter for extension & shrink move
  ShrinkMove    = 0;
  ShrinkAcc     = 0;
  ShrinkRej     = 0;
  ExtensionMove = 0;
  ExtensionAcc  = 0;
  ExtensionRej  = 0;
  
  stack_break_count = 0;
  BONDIHED_PBPB_RCUT_count = 0;
  BONDIHED_BPBP_RCUT_count = 0;
  
  bp_extension_acc = 0;
  bp_extension_rej = 0;

  /*counter foe wang landau*/
  count = 0;

  //debug counter
  LineCount = 0;

  setupList();


  for (int i=0; i < N; i++){
    part[i].myId = i;
    if( i % 2 == 0 ){
      part[i].beadType = 'P';
    }else{
      part[i].beadType = 'B';
    }
  }

  //setup bond pointers
  {

      for (int i=0; i< (N/2) - 2; i++){
        if( part[i].beadType == 'P' ){
          //setup bonds for phosphates
            part[i].bNext   = &part[i+1];
            part[i+1].pPrev = &part[i];
            part[i].pNext   = &part[i+2];
            part[i+2].pPrev = &part[i];

        }else{
           //setup bonds for bases
            part[i].bNext   = &part[i+2];
            part[i+2].bPrev = &part[i];
            part[i].pNext   = &part[i+1];
            part[i+1].bPrev = &part[i];
        }
      }

      if( part[(N/2)-1].beadType == 'P' ){
          cerr << "Probably better to have a multiple of 4 particles: 1bp = 2 bases+2 phosphates." << endl;
          exit( 8 );
      }else{
            //setup bonds for bases connect particles to their periodic image
            part[(N/2)-1].bNext   = &part[1];
            part[1].bPrev         = &part[(N/2)-1];
            part[(N/2)-1].pNext   = &part[0];
            part[0].bPrev         = &part[(N/2)-1];
      }

      if( part[(N/2)-2].beadType == 'P' ){
          //setup bonds for phosphates
            part[(N/2)-2].bNext  = &part[(N/2)-1];
            part[(N/2)-1].pPrev  = &part[(N/2)-2];
            part[(N/2)-2].pNext  = &part[0];
            part[0].pPrev        = &part[(N/2)-2];
      }else{
            cerr << "Probably better to have a multiple of 4 particles: 1bp = 2 bases+2 phosphates." << endl;
            exit( 8 );
      }


      for (int i=N/2; i< N - 2; i++){
          if( part[i].beadType == 'P' ){
              //setup bonds for phosphates
              part[i].bNext   = &part[i+1];
              part[i+1].pPrev = &part[i];
              part[i].pNext   = &part[i+2];
              part[i+2].pPrev = &part[i];
          }else{
              //setup bonds for bases
              part[i].bNext   = &part[i+2];
              part[i+2].bPrev = &part[i];
              part[i].pNext   = &part[i+1];
              part[i+1].bPrev = &part[i];
          }
      }


      if( part[N-1].beadType == 'P' ){
          cerr << "Probably better to have a multiple of 4 particles: 1bp = 2 bases+2 phosphates." << endl;
          exit( 8 );
      }else{
            //setup bonds for basesconnect particles to their periodic image
            part[N-1].bNext        = &part[(N/2) + 1];
            part[(N/2) + 1].bPrev  = &part[N-1];
            part[N-1].pNext        = &part[(N/2)];
            part[(N/2)].bPrev      = &part[N-1];
      }

      if( part[N-2].beadType == 'P' ){
          //setup bonds for phosphates connect particles to their periodic image
          part[N-2].bNext   = &part[N-1];
          part[N-1].pPrev   = &part[N-2];
          part[N-2].pNext   = &part[(N/2)];
          part[(N/2)].pPrev = &part[N-2];
       }else{
            cerr << "Probably better to have a multiple of 4 particles: 1bp = 2 bases+2 phosphates." << endl;
            exit( 8 );
       }

    }

    cout << "#read box: " << box->x[0]     << " " << box->x[1]     << " " << box->x[2]    << endl;
    cout << "#half box: " << box->halfx[0] << " " << box->halfx[1] << " " << box->halfx[2] << endl;


    //find the potential energy
    eBondTot  = bond_totalEnergy();
    cout << "#eBondTot " << eBondTot << endl;

    eSAKTot   = SAK_totalEnergy();
    cout << "#eSAKTot " << eSAKTot << endl;

    eDihedTot = totalDihedEnergy();
    cout << "#eDihedTot " << eDihedTot << endl;

    eLJTot    = LJ_totalEnergy();
    cout << "#eLJTot " << eLJTot << endl;

    ePotTot   = eLJTot + eBondTot + eDihedTot + eSAKTot;
    cout << "#ePotTot " << ePotTot << endl;
  }


void hsc::setupList(){

  int icell;

  //Maximum interaction length
  //double rc = 12.0;            ///WARNING LENNARD_JONES CUTOFF BASED ON THE MAXIMUM CUT OFF RANGE FOR PHOSPHATES //
  box->rcut = rcut;
  adaptive->maxDisplace = maxDisplace;
  adaptive_rigid->theta_rigid = thetaRigid;
  box->update();

  //cerr << "Allocating : " << box->nCells << " cells, total size: " << box->nCells * sizeof(Cell) << " bytes " << endl;

  cellHash = new HashTable();
  cout << "Rcut       : " << rcut << endl;
  cout << "Box        : " << box->x[0] << " " << box->x[1] << " " << box->x[2] << endl;

  //put particles into cells
  for (int i=0; i<N; i++){
    icell = cellIndex( part[i].R, box );
    insertToCell(&part[i], cellHash, icell, box);
  }


  //cellHash->getNumberOfItems();

  // find neighbouring cells
  //setUpNeighbours();
}


void hsc::setUpNew(){


  ifstream infile;
  char     line[200];
  char     quantity[200];
  float    number;
  int      check, seedSet;

  /*DNA structural parameters for setup*/
  double  R, h,  rBase, theta, Incl;
  double  x, y, z;

  open_chain    = 0;
  iran          = &seed;
  seedSet       = seed = -1;
  stepsize_box_turns_per_image = 1;
  int  turns_per_image_local = 0;

  bx = 0.;

  //default BP parameters
  R       =  9.49;   //outer distance between phosphates
  rBase   =  3.49;
  h       =  0.;     //3.16 is reasonable, but find this from box+N
  theta   =  0.598;  //(34.3 degree) rotation per base pair 0.598
  Incl    = -0.035;  //base pair inclination
  initialize_chain_as_linear = 1; //not a circle.

  use_constant_tension = 0;
  tension_pn           = 0.;
  mu_interCalator      = 0.;
  intercal_count       = 0;
  sdna_count           = 0;
  intercal_r0          = 0;
  intercal_epsilon     = 0;
  r0_pp_intercal       = 4.; //default
  Epsilon_one          = 0;
  ETA                  = 0;
  Delta                = 0;
  sStretch_lim         = 0;
  lambda               = 0;
  U_max                = 0;
  
  
  
  infile.open(ParameterFile);
  if (infile.bad()) {
    cerr << "Error " << ParameterFile << " not found\n";
    exit(8);
  }
  cerr << "Opened file, name: " << ParameterFile << endl;

  while (infile.peek() != EOF) {
    infile.getline(line,sizeof(line));
    cerr << "Read line: " << line << endl;
    check = sscanf(line, "%s%f", quantity, &number);
    if (check != 2) {
      cerr << "Error: wrong line format in " << ParameterFile << endl;
      cerr << "Line was: " << line << endl;
      exit(8);
    }

    cerr << "quant:" << quantity << " numb: " << number << endl;
    if (strcmp(quantity,"Number")==0) {
      N = int(number);
    }

	else if (strcmp(quantity, "MaxDisplacement")==0) {
        maxDisplace = double(number);
    }
    else if (strcmp(quantity,"Seed")==0) { seedSet = 0; seed=int(number); cout << "read seed " << seed <<  endl; }
    else if (strcmp(quantity,"McSteps")==0)     mcSteps=int(number);
    else if (strcmp(quantity,"McEvalSteps")==0) mcEvalSteps=int(number);
    else if (strcmp(quantity,"McSnapSteps")==0) mcSnapSteps=int(number);
    else if (strcmp(quantity,"EqSteps")==0)     eqSteps=int(number);
    else if (strcmp(quantity,"EqEvalSteps")==0) eqEvalSteps=int(number);
    else if (strcmp(quantity,"EqSnapSteps")==0) eqSnapSteps=int(number);
    else if (strcmp(quantity,"Rcutoff")==0)     rcut=number;
    else if (strcmp(quantity,"theta" )==0)      theta=double(number);
    else if (strcmp(quantity,"initialize_chain_as_linear" )==0) initialize_chain_as_linear=int(number);
    else if (strcmp(quantity,"openchainDNA" )==0) open_chain=int(number);
    else if (strcmp(quantity,"h" )==0)          h=double(number);
    else if (strcmp(quantity,"Boxside" )==0)    bx=double(number);
    else if (strcmp(quantity,"Extension_step" )==0)     extension_step               = double(number);
    else if (strcmp(quantity,"Use_constant_tension" )==0) use_constant_tension       = int(number);
    else if (strcmp(quantity,"mu_interCalator" )==0)      mu_interCalator            = double(number);
    else if (strcmp(quantity,"Tension_pn" )==0)           tension_pn                 = double(number);
    else if (strcmp(quantity,"r0_pp_intercal" )==0)       r0_pp_intercal             = double(number);
    else if (strcmp(quantity,"intercal_delta" )==0)       intercal_delta             = double(number);
    else if (strcmp(quantity,"Shrink_step" )==0)        shrink_step                  = double(number);
    else if (strcmp(quantity,"RigidTheta" )==0)         thetaRigid                   = double(number);
    else if (strcmp(quantity,"lambda" )==0)                  lambda                  = double(number); 
    else if (strcmp(quantity,"U_max" )==0)                   U_max                   = double(number);
    else if (strcmp(quantity,"ETA" )==0)                     ETA                     = double(number);
    else if (strcmp(quantity,"Delta" )==0)                   Delta                   = double(number);
    else if (strcmp(quantity,"sStart_x" )==0)                sStart_x                = double(number); 
#ifdef INTERCAL
    else if (strcmp(quantity,"BonDihed_BPBP_Rcut" )==0)      BonDihed_BPBP_Rcut      = double(number); 
    else if (strcmp(quantity,"BonDihed_PBPB_Rcut" )==0)      BonDihed_PBPB_Rcut      = double(number); 
    else if (strcmp(quantity,"Epsilon_one" )==0)             Epsilon_one             = double(number);
    else if (strcmp(quantity,"intercal_r0" )==0)             intercal_r0             = double(number); 
    else if (strcmp(quantity,"intercal_epsilon" )==0)        intercal_epsilon        = double(number); 
    else if (strcmp(quantity,"sStretch_lim" )==0)            sStretch_lim            = double(number);
#endif
    else if (strcmp(quantity,"start_box_turns_per_image" )==0)      turns_per_image_local = int(number);
    else if (strcmp(quantity,"stepsize_box_turns_per_image" )==0)   stepsize_box_turns_per_image = int(number); //global variable to reset the box rotation
  }

  infile.close();

  cerr << N << " Particle Numbers " << endl;
  if (h != 0.)  //if base rise imno passed instead of box.
    bx =  h  * (N/4);
  else
    h  =  bx / (N/4);
  part           = new Particle[N];
  box            = new Box;
  adaptive       = new Adaptivestepsize;
  adaptive_rigid = new AdaptiveStepRigid;
#ifdef INTERCAL
  adaptive_ext   = new AdaptiveExtension;
#endif

  box->x[0] = bx;
  box->x[1] = bx;
  box->x[2] = bx;

  box->rcut                   = rcut;
  adaptive->maxDisplace       = maxDisplace;
  adaptive_rigid->theta_rigid = thetaRigid;
  box->update();
#ifdef INTERCAL
  adaptive_ext->BondThreshold      = BondThreshold;
  adaptive_ext->BonDihed_BPBP_Rcut = BonDihed_BPBP_Rcut;
  adaptive_ext->BonDihed_PBPB_Rcut = BonDihed_PBPB_Rcut;
  adaptive_ext->length_init        = bx; //initial lengthg of DNA
  adaptive_ext->force_imposed      = tension_pn;
#endif

/*set the initial values for wanglandau histogram.bin number for twist,stretch and bin width for twist and stretch respectively.*/

#ifdef WANGLANDAU
  wL         = new wangLandau(5, 9, 0.25, 0.028);
  count_hist = wL->create_hist(5, 9);
  bias_hist  = wL->create_hist(5, 9);
#endif


#ifdef ROTORBOX
  init_rotorBox(turns_per_image_local);
#endif

  cout << "# read parameters" << endl;
  cout << "# N " << N << endl;
  cout << "# maxDisplacement " << maxDisplace << endl;
  cout << "# Boxside " << bx << " base-pair distance: " << h << endl;
  cout << "# Eq " << eqSteps << " EqEval " << eqEvalSteps << " EqSnap "
       << eqSnapSteps << endl;
  cout << "# Mc " << mcSteps << " McEval " << mcEvalSteps << " McSnap "
       << mcSnapSteps << endl;
  cout << "# seed " << seed << " Displ " << maxDisplace << endl;
  cout << "# theta " << theta << endl;
  cout << "# initialize_chain_as_linear " << initialize_chain_as_linear << endl;
  cout << "# open chain DNA " << open_chain << " Warning: this flag makes a conformationally and configurationally linear DNA " << endl;
  cout << "#Extension step: " << extension_step << endl;
  cout << "#theta rigid: " << thetaRigid << endl;
  cout << "#step size for box rotation: " << stepsize_box_turns_per_image << endl;
#ifdef INTERCAL
  cout << "#sStart_x: "        << sStart_x         << endl;
  cout << "#Umax: "            << U_max            << endl;
  cout << "#Mu intercal: "     << mu_interCalator  << endl;
  cout << "#lambda: "          << lambda           << endl;
  cout << "#Delta : "          << Delta            << endl;
  cout << "#ETA: "             << ETA              << endl;
  cout << "#sStretch_lim: "    << sStretch_lim     << endl;
  cout << "#intercal_r0: "    << intercal_r0      << endl;
  cout << "#sdna_count init: "     << sdna_count       << endl;
  cout << "#intercal_count init: " << intercal_count   << endl;
  cout << "#intercal_epsilon: " << intercal_epsilon << endl;
#endif

  if( seedSet == -1 ){
        cerr << "Did not seed RNG!" << endl;
        exit(8);
  }

  //** length of DNA
  double L = bx; //DNA length. Box length is adjusted with DNA length

  //** radius of Torus, parameter set for circular DNA
  double R_torus   = L / (2*M_PI);
  double del_theta = (2*M_PI)/(N/4);
  double del_phi   = theta;
  double u;


  if( L > box->x[2] ){
      cerr << "Warning... setup will place long double-helix in a short box. l_helix: " << h*(N/4.0) << " in: " << box->x[2]<< endl;
  }
  cerr << "DNA length: " << L << endl;

  //number of turns imposed on the DNA at setup.
  linking_number = int(0.5 + theta * (N/4.) / (2.*M_PI));
  cout << "##Initial linking number is round of: " << theta * (N/4.) / (2.*M_PI) <<"\n";
  cout << "##                                = " << linking_number <<"\n";

  for(int bpIndex=0; bpIndex<N/4; bpIndex++){

        //** if initialize_chain_as_linear = 1 then create linear
       //** otherwise circular
       //first base pair
  if(initialize_chain_as_linear){
    x = (R*cos(Incl)) * cos(theta*bpIndex);
    y = (R*cos(Incl)) * sin(theta*bpIndex);
    z = bpIndex*h + R*sin(Incl) - box->halfx[2];
    z += 0.35;
  }else{
    u = R_torus+(R/2)*cos(-del_phi*bpIndex);
    x = (R/2)*sin(-del_phi*bpIndex);
    y = u*cos(del_theta*bpIndex);
    z = u*sin(del_theta*bpIndex);
  }

  //First phosphate bead
  part[bpIndex * 2].R[0] = x;
  part[bpIndex * 2].R[1] = y;
  part[bpIndex * 2].R[2] = z;

  if(initialize_chain_as_linear){
    x = (rBase*cos(Incl)) * cos(theta*bpIndex);
	y = (rBase*cos(Incl)) * sin(theta*bpIndex);
    z =  bpIndex*h + rBase*sin(Incl) - box->halfx[2];
    z += 0.35;

  }else{
  	u = R_torus+(rBase/2)*cos(-del_phi*bpIndex);
  	x = (rBase/2)*sin(-del_phi*bpIndex);
  	y = u*cos(del_theta*bpIndex);
  	z = u*sin(del_theta*bpIndex);
  }
    //first base bead
    part[bpIndex * 2 + 1].R[0] = x;
    part[bpIndex * 2 + 1].R[1] = y;
    part[bpIndex * 2 + 1].R[2] = z;

  //second BP
  if(initialize_chain_as_linear){
	x = (R * cos(Incl)) * cos(theta * bpIndex + M_PI );
	y = (R * cos(Incl)) * sin(theta * bpIndex + M_PI );
	z = bpIndex * h - R*sin(Incl)  - box->halfx[2] ;
    z += 0.35;


  }else{
	u = R_torus+(R/2)*cos(-del_phi*bpIndex + M_PI);
	x = (R/2)*sin(-del_phi*bpIndex + M_PI);
	y = u*cos(del_theta*bpIndex);
	z = u*sin(del_theta*bpIndex);
 }

    //Second phosphate bead
    part[N - bpIndex*2 - 2].R[0] = x;
    part[N - bpIndex*2 - 2].R[1] = y;
    part[N - bpIndex*2 - 2].R[2] = z;

  if(initialize_chain_as_linear){
	x = (rBase*cos(Incl)) * cos(theta*bpIndex + M_PI );
	y = (rBase*cos(Incl)) * sin(theta*bpIndex + M_PI );
    z = bpIndex*h - rBase*sin(Incl) - box->halfx[2] ;
    z += 0.35;
  }else{
	u = R_torus+(rBase/2)*cos(-del_phi*bpIndex+M_PI);
	x = (rBase/2)*sin(-del_phi*bpIndex+M_PI);
	y = u*cos(del_theta*bpIndex);
	z = u*sin(del_theta*bpIndex);
  }
    //Second base bead
    part[N - bpIndex*2 - 1].R[0] = x;
    part[N - bpIndex*2 - 1].R[1] = y;
    part[N - bpIndex*2 - 1].R[2] = z;
  }



#ifdef ROTORBOX
//if rotating the box, then particle zero should always be fixed at the bottom plane of the system.
for(int i=0; i< N; i++){
        //part[i].R[0] -= (part[0].R[0]);
        //part[i].R[1] -= (part[0].R[1]);
    //cerr << " part 0 before change: " << part[0].R[2] << endl;
        part[i].R[2] -= (part[0].R[2] + box->halfx[2] - 1e-9);

        intoBox(part[i].R, box);
}
/////for debug//////
part_init = part[0].R[2];
#endif
       cerr << "fixed part 1 in z plane." << " part 0,Z: " << part[0].R[2] << " part 1,Z: " << part[1].R[2] << endl;

  /* image particles into the box... not needed, already checked that z was less than h */
  for(int i=0; i< N; i++){
    if ( fabs(part[i].R[2]) > box->halfx[2] ){
       cerr << "Error, did not initialize inside box" << endl;
       cerr << " particle's Z: " << part[i].R[2] << " Box_Z: " << box->halfx[2] << endl;


       cerr << " particle: " << i << " X: " << part[i].R[0]
                                  << " Y: " << part[i].R[1]
                                  << " Z: " << part[i].R[2] << endl;
       exit(1);
    }
  }

  initialize();

}


void hsc::setUpFromFile() {

  ifstream infile, infile2;
  char line[200];
  char quantity[200];
  float number;
  int check, seedSet;
  iran          = &seed;
  seedSet       = seed = -1;

  /*DNA structural parameters for setup*/
  double  h, theta;
  open_chain    = 0;
  stepsize_box_turns_per_image = 1;
  int  turns_per_image_local   = 0;

  use_constant_tension = 0;
  tension_pn           = 0.;
  mu_interCalator      = 0.;
  intercal_count       = 0;
  sdna_count           = 0;
  
  // read new simulation parameters from ParameterFile
  infile.open(ParameterFile);

  if (infile.bad()) {
      cerr << "Error " << ParameterFile << " not found\n";
      exit(8);
  }
  while (infile.peek() != EOF) {
    infile.getline(line,sizeof(line));
    check = sscanf(line, "%s%f", quantity, &number);
    if (check != 2) {
      cerr << "Error: wrong line format in " << ParameterFile << line << endl;
      exit(8);
    }
    else if (strcmp(quantity, "MaxDisplacement")==0) {
        maxDisplace = double(number);
    }

    else if (strcmp(quantity,"Seed")==0) { seedSet = 0; seed=int(number); cout << "read seed " << seed <<  endl; }
    else if (strcmp(quantity,"McSteps")==0)     mcSteps=int(number);
    else if (strcmp(quantity,"McEvalSteps")==0) mcEvalSteps=int(number);
    else if (strcmp(quantity,"McSnapSteps")==0) mcSnapSteps=int(number);
    else if (strcmp(quantity,"EqSteps")==0)     eqSteps=int(number);
    else if (strcmp(quantity,"EqEvalSteps")==0) eqEvalSteps=int(number);
    else if (strcmp(quantity,"EqSnapSteps")==0) eqSnapSteps=int(number);
    else if (strcmp(quantity,"Rcutoff")==0)     rcut=number;
    else if (strcmp(quantity,"theta" )==0)      theta=double(number);
    else if (strcmp(quantity,"openchainDNA" )==0) open_chain=int(number);
    else if (strcmp(quantity,"h" )==0)          h=double(number);
    else if (strcmp(quantity,"Boxside" )==0)    bx=double(number);
    else if (strcmp(quantity,"Extension_step" )==0)     extension_step               = double(number);
    else if (strcmp(quantity,"Use_constant_tension" )==0) use_constant_tension       = int(number);
    else if (strcmp(quantity,"Mu_intercalator" )==0)      mu_interCalator            = double(number);
    else if (strcmp(quantity,"Tension_pn" )==0)           tension_pn                 = double(number);
    else if (strcmp(quantity,"Shrink_step" )==0)        shrink_step                  = double(number);
    else if (strcmp(quantity,"RigidTheta" )==0)         thetaRigid                   = double(number);
    else if (strcmp(quantity,"Linking_num" )==0)        linking_number               = double(number);
    else if (strcmp(quantity,"start_box_turns_per_image" )==0)      turns_per_image_local = int(number);
    else if (strcmp(quantity,"stepsize_box_turns_per_image" )==0)   stepsize_box_turns_per_image = int(number); //global variable to reset the box rotation
  }

  infile.close();
  adaptive       = new Adaptivestepsize;
  adaptive_rigid = new AdaptiveStepRigid;
  box            = new Box;
  
  
  box->rcut                   = rcut;
  adaptive->maxDisplace       = maxDisplace;
  adaptive_rigid->theta_rigid = thetaRigid;
  

#ifdef ROTORBOX
  init_rotorBox(turns_per_image_local);
#endif

  cout << "#read "   << ParameterFile << endl;
  cout << "#hereX: " << configFileName << endl;
  cout << "#rcut: "  << rcut << endl;
  cout << "#max displace: " << maxDisplace << endl;
  cout << "#start_box_turns_per_image: " << turns_per_image_local << endl;
  cout << "#stepsize_box_turns_per_image: " << stepsize_box_turns_per_image << endl;
  cout << "#Linking number: " << linking_number << endl;
  cout << "#thetaRigid: " << thetaRigid << endl;


  readXYZ(configFileName);
  

  //box->update();
  cout << " reading config file" << endl;

  iran = &seed;
  cout << " seed: " << seed << endl;


  initialize();

  cout << "#read parameters and configuration" << endl;
  cout << "#N " << N << endl;
  cout << "#MaxDisplacement " << maxDisplace << endl;
  cout << "#Eq " << eqSteps << " EqEval " << eqEvalSteps << " EqSnap "
       << eqSnapSteps << endl;
  cout << "#Mc " << mcSteps << " McEval " << mcEvalSteps << " McSnap "
       << mcSnapSteps << endl;
  cout << "#seed " << seed << " Displ " << maxDisplace << endl;
  cout << "#Extension step: " << extension_step << "#Shrink step: " << shrink_step << " Rigid theta: " << thetaRigid << endl;


}
