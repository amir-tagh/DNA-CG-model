#include <stdio.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include "FreeEnergyHS.h"
#include "tools.h"
#include "cell.h"
#include "hashTable.h"
#include "wangLandau.h"
#include <time.h>
#include <typeinfo>
using namespace std;

//#define _DEBUG_BOND_ENERGY
const double Pi = M_PI;

extern float ran2(long *idum);

void hsc::run(){

  if (newSetUp)setUpNew();
  else setUpFromFile();
      //if( part[0].R[2] != 1e-9 - box->halfx[2] ){
        //    cerr << "Particle 0 is awry after setup. " << part[0].R[2] << " Part one: " <<
          //  part[1].R[2] + box->halfx[2] - 1e-9 << endl;  
            //exit(1);
      //}
  writeOutConfig("StartConfig.xyz");
  equil();
  MC();

  cout << "##Calculation reached max steps, saving stats and quitting." << endl;
  writeOutSimData("Stats");
}

void hsc::writeSystemEnergy(){

  cout << "#systemEnergy "
       << eBondTot/N  << " "
       << eSAKTot/N   << " "
       << eDihedTot/N << " "
       << eLJTot/N    << " "
       << ePotTot/N   << " "
       
#ifdef INTERCAL
       << intercal_count*4/float(N) << " "
#endif
#ifndef INTERCAL
       << sdna_count*4/float(N)     << " "
#endif
#ifdef ROTORBOX
       << linking_number << " ";
#else
       ;
#endif
    
  if( use_constant_tension ){
     cout << box->x[2] << endl; 
  }
  cout << endl;
       
}

//fix numerical drift of total potential energy
double hsc::resetSystemEnergy(){

  double delta_bond;
  double e_bond;
  double delta_dihed;
  double e_dihed;
  double delta_LJ;
  double e_LJ;
  double delta_SAK;
  double e_SAK;

  cout << endl;
  cout << "#resetSystemEnergy(before) "
        << eBondTot/N  <<     " "
        << eSAKTot/N   <<     " "
        << eDihedTot/N <<     " "
        << eLJTot/N    <<     " "
        << ePotTot/N   <<     " "
        << endl;


  e_bond      = bond_totalEnergy();
  delta_bond  = e_bond - eBondTot;
  eBondTot    = e_bond;

  e_dihed     = totalDihedEnergy();
  delta_dihed = e_dihed - eDihedTot;
  eDihedTot   = e_dihed;

  e_LJ        = LJ_totalEnergy();
  delta_LJ    = e_LJ - eLJTot;
  eLJTot      = e_LJ;

  e_SAK       = SAK_totalEnergy();
  delta_SAK   = e_SAK - eSAKTot;
  eSAKTot     = e_SAK;

  double deltaTot = delta_bond + delta_dihed + delta_LJ + delta_SAK;
  ePotTot = eBondTot + eSAKTot + eDihedTot + eLJTot;

  cout << "#resetSystemEnergy(after) "
        << eBondTot/N  <<     " "
        << eSAKTot/N   <<     " "
        << eDihedTot/N <<     " "
        << eLJTot/N    <<     " "
        << ePotTot/N   <<     " "
        << deltaTot/N  << endl;

   return ( deltaTot );

}

void hsc::equil(){

  int eqInterval;
  if (eqEvalSteps > 0 )
      eqInterval = eqSteps/eqEvalSteps;
  else eqInterval = 0;
  int measuring = eqInterval;
  int snapInterval;
  if (eqSnapSteps > 0)
      snapInterval = eqSteps/eqSnapSteps;
  else snapInterval = 0;
  int snapshot = snapInterval;

  cout << setw(20) << "Bond"
       << setw(8)  << "SAK"
       << setw(11) << "eDihed"
       << setw(8)  << "eLJ"
       << setw(8)  << "ePot"
#ifdef INTERCAL
       << setw(12)  << " Intercal "
#endif
#ifndef INTERCAL
       << setw(12)  << " S_tag "
#endif       
#ifdef ROTORBOX
       << setw(10)  << "Link"
#endif
       << setw(8)  << "box_z"
       <<endl;

  ofstream myfile;
  myfile.open (EquilFileName.c_str(), ios::out);

  for (int steps=1; steps<=eqSteps; steps++){

      measuring--;
      snapshot--;

      cout << steps - 1 << " ";
      writeSystemEnergy();
      
      RigidBodyRotation();

      for(int i=0; i<N; i++){
          particleMove();
      }

#ifdef INTERCAL
      //very cheap move, only bond terms.
      for(int i=0; i<N/4; i++){
        intercalate_bpMove();
      }
#endif


#ifndef INTERCAL
     //very cheap move, only bond terms.
      for(int i=0; i<N/4; i++){
        SDNA_baseMove();
      }
#endif

#ifndef WANGLANDAU
      topology_rotation();
#endif
       extension_bpStep();
//if ROTORBOX is defined, just fix one bead: no centring, which introduces artifacts.
#ifndef ROTORBOX
centreMolecule();
#endif

    //Do Rigid Rotation every RESET_ENERGY_TOTALS steps
    //in order to correct roundoff-drift
    if (steps % RESET_ENERGY_TOTALS_EVERY == 0) {

        double delta;
        delta = resetSystemEnergy();

        cerr << steps -1 << " reset energy, delta: " << delta << " of " << ePotTot << endl;

        if( fabs( delta ) > 0.01 && steps > 0){
            cerr << "Large drift in energy during equil, consider resetting more often: " << delta << endl;
            exit( 1 );
        }
    }
    
     if(measuring==0){
        measuring = eqInterval;
        runCount++;
     }

     if(snapshot == 0){
        snapshot = snapInterval;
     }
     //writeXYZ_frame( &myfile );
  }
  writeXYZ_frame( &myfile );
  cerr << " steps: " << eqSteps << endl;
  cout << "# Equilibration ended" << endl;

  myfile.close();
}

/*set the boolean for productionStage*/

void hsc::MC(){

  char outstring[16];
  char temp[16];
  int mcInterval;

  if (mcEvalSteps > 0 )
      mcInterval = mcEvalSteps;
  else mcInterval = 0;

  int measuring = mcInterval;
  int snapInterval;

  if (mcSnapSteps > 0)
      snapInterval  = mcSnapSteps;
  else snapInterval = 0;

  int snapshot = snapInterval;
  measureCount = 0;

  //!open trajectory file called "trajFileStream" here:
  ofstream myfile;
  myfile.open (outputFileName.c_str(), ios::out);
  
  ofstream Ext_force;
  Ext_force.open(Ext_forceFileName.c_str(), ios::out);

  cout << "#Running " << mcSteps << " production, printing snaps every " << snapInterval << " energies every " << mcInterval << endl;

  cout << setw(25) << "Bond"
       << setw(9)  << "SAK"
       << setw(9)  << "eDihed"
       << setw(9)  << "eLJ"
       << setw(11) << "ePot"
       //<< setw(12) << "deltaTot"
#ifdef INTERCAL
       << setw(12)  << "Intercal"
#endif
#ifndef INTERCAL
       << setw(12)  << "S_tag"
#endif
#ifdef ROTORBOX
       << setw(12)  << "Link"
#endif
       << setw(12)  << "box_z"
       << endl;

  for (int steps=1; steps<=mcSteps; steps++ ){
      measuring--;
      snapshot--;

//MC movements
       RigidBodyRotation();
       
       for(int i=0; i<N; i++){
           particleMove();
       }

#ifdef INTERCAL
      //very cheap move, only bond and dihedral terms.
      for(int i=0; i<N/4; i++){
        intercalate_bpMove();
      }
#endif

#ifndef INTERCAL
     //very cheap move, only bond terms.
      for(int i=0; i<N/4; i++){
        SDNA_baseMove();
      }
#endif

       topology_rotation();
       
       extension(&Ext_force);
       extension_bpStep();

#if 0       
       if ( steps % 10000 == 0 ){
           stack_break();
           dihed_break();
       }
#endif
#ifdef WANGLANDAU
       extension(&TransMatExt, &Ext_energy_acc, &Ext_energy_rej);
       shrink(&TransMatShrink, &Shr_energy_acc, &Shr_energy_rej);
#endif

#ifndef ROTORBOX //if doing rotorbox, atom 0 is never moved: this happens instead of centering.
      centreMolecule();
#endif

    //in order to correct roundoff-drift
    if (steps % RESET_ENERGY_TOTALS_EVERY == 0 ) {

        cout << steps -1 << " ";

        double delta;
        delta = resetSystemEnergy();

        if( fabs( delta ) > 0.1 ){
            cerr << "Large drift in energy during MC, consider resetting more often: " << delta << endl;
            exit( 1 );
        }
    }

/*flat histogram*/
/*collect data while running the algorithm*/
#ifdef WANGLANDAU
    ofstream out_visit_count("visit_count.txt", ios::app);
    ofstream out_hist_bias("histogram.txt", ios::app);
    ofstream out_PMF("histogram_plusPMF.txt", ios::app);
    if (steps % 100000 == 0){

        int twist_bin_num;
        int stretch_bin_num;;
        int i, j;
        double min;
        double **PMF;
        double **histogram_plusPMF;

        /*two new histograms to collect data*/
        /*WARNING: matrix sizes here should be the same as matrixes in setup*/
        PMF = wL->create_hist(5,9);
        histogram_plusPMF = wL->create_hist(5,9);

        twist_bin_num   = wL->n_twist;
        stretch_bin_num = wL->n_stretch;
        min = wL->Min_hist(bias_hist);

        for( i = 0; i < stretch_bin_num; i++){
            for( j = 0; j < twist_bin_num; j++){
                bias_hist[i][j] -= min;
            }
        }

        for( i = 0; i < stretch_bin_num; i++){
            for( j = 0; j < twist_bin_num; j++){
                out_visit_count << count_hist[i][j] << endl;
                out_hist_bias   << bias_hist[i][j]  << endl;
            }
        }

        /*PMF histogram*/
        for( i = 0; i < stretch_bin_num; i++){
            for( j = 0; j < twist_bin_num; j++){
                PMF[i][j] = -1*log(count_hist[i][j]);
            }
        }

        for( i = 0; i < stretch_bin_num; i++){
            for( j = 0; j < twist_bin_num; j++){
                histogram_plusPMF[i][j] = PMF[i][j] +
                                          bias_hist[i][j];

            }
        }

        for( i = 0; i < stretch_bin_num; i++){
            for( j = 0; j < twist_bin_num; j++){
                out_PMF << histogram_plusPMF << endl;
            }
        }
    }

    out_visit_count.close();
    out_hist_bias.close();
    out_PMF.close();

/*wang landau flatness check*/
/*WLCHECK = 1e6*/
    if ( steps > WLCHECK  && ( !wL->productionStage)){

        double mean, min;
        bool   productionStage;
        bool   isflat;

        /*see if the PMF is flat enough to not need any more histogram-bias building*/
        productionStage = wL->histogram_finished(count_hist);
        wL->productionStage = productionStage;                    /*if PMF is flat enough start collecting producticve data*/

        if ( productionStage == true )
            continue;

        mean = wL->getMean(count_hist);
        min  = wL->Min_hist(count_hist);

        isflat = wL->checkFlatness(min, mean);
        if (isflat == true){
            wL->hist_increment *= 0.9;
            wL->flatness_test = sqrt(wL->flatness_test);
            cout << " reset update factor to: " << wL->hist_increment << " and: "
                 << "flatness test to: "        << wL->flatness_test  << endl;
            wL->make_zero(count_hist);                         /*if flatness checks out set the histogram to zero*/
        }
    }
#endif

      if (measuring ==0){
          measureCount++;
          measuring = mcInterval;
#define _DEBUG_COLLISIONS
#ifdef  DEBUG_COLLISIONS
      cout << "#hash queries, collisions: " << queries        << " "        << collisions      << endl;
      cout << "#Item up, item down: "       << itemsInHashup  << "  "       << itemsInHashdown << endl;

#endif
        cout << steps - 1 << " ";
        writeSystemEnergy();
      }

     if(snapshot == 0) {
        snapshot = snapInterval;
        strcpy(outstring, "CoMC");
        sprintf(temp,"%.*d", 0, int(steps/snapInterval));
        strcat(outstring, temp);

        writeXYZ_frame( &myfile );
    }
  }
   myfile.close();
   Ext_force.close();

  //write the finishing state to a restart file
   myfile.open(restartFileName.c_str(),ios::out);
   writeXYZ_frame( &myfile );
   myfile.flush();
   myfile.close();

}
#define _TEST_RB

void hsc::particleMove(){

  int       pNr;
  Particle *p, *q;
  double    Rold[3];
  unsigned int  newCell_id;
  Cell     *oldCell, *newCell;
  double    oldEBond, oldeDihed, oldeSAK;
  double    newEBond, neweDihed, neweSAK;
  double    oldE, newE;
  double    oldeLJ = 0.0, neweLJ =0.0;

  transMoves ++;

  //do the adaptive step size
  if(transMoves > 1000 && transMoves < 50000){

    adaptive->transMoves  = transMoves;
    adaptive->transAccept = transAccept;
    adaptive->Adaptive_step_size();
}


#ifdef TEST_RB
 //writeOutConfig("test_particle_start.xyz");
// print_bond_totalEnergy("test_RB_start_bondEnergy.dat");
// print_lj_totalEnergy("test_RB_start_ljEnergy.dat");
 
   double testOldEdihed = totalDihedEnergy();
   double testOldEBond  = bond_totalEnergy();
   double testOldELJ    = LJ_totalEnergy();
#endif

#ifdef ROTORBOX
  //hold particle zero permanently fixed if allowing a rotor box.
  pNr =  int((N-1)*ran2(iran));
  pNr++; 
#else
  pNr =  int((N)*ran2(iran));
#endif

  p = &part[pNr];

  // save old configuration
  for (int i=0; i<3; i++){
    Rold[i] = p->R[i];
  }
  oldCell = p->cell;

  oldEBond  = bond_onePartEnergy(p);
  oldeDihed = dihed_onePartEnergy(p);
  oldeSAK   = SAK_pairEnergy(p, &part[SAK_myPair(p, pNr, N)]);
  oldeLJ    = 0.0;

 /* Loop over particles in this cell */
  q = oldCell->firstParticle;
  do{
    if ( q != p){
      oldeLJ    += LJ_pairEnergy(p, q);
    }
    q = q->next;
  } while( q );

  /* Loop over particles in the 26 neighbour cells */
  for( int i = 0; i < 26; i++){
    Cell *c = cellHash->getItemByKey(oldCell->neighbours[i]);
    if( c == NULL ) continue;

    q = c->firstParticle; //to update the cell if necessary
    while( q ){
        oldeLJ    += LJ_pairEnergy(p, q);
        q = q->next;
    }
  }
  /* move the particle */
  //double oldEDtot, newEDtot, oldED_1;

  //oldEDtot = totalDihedEnergy();
  //oldED_1  = dihed_onePartEnergy(p);


  displaceParticle(p->R, Rold); //displace and image particle.
  newCell_id = cellIndex(p->R, box);
  if ( newCell_id != oldCell->myId ){
      newCell = moveBetweenCells(p, oldCell, newCell_id, cellHash, box);
  }else{
      newCell = oldCell;
  }

  newEBond  = bond_onePartEnergy(p);
  neweDihed = dihed_onePartEnergy(p);
  neweSAK   = SAK_pairEnergy(p, &part[SAK_myPair(p, pNr, N)]);

  //debugging the dihed energy
  //newEDtot = totalDihedEnergy(NULL);
  //if( fabs( (newEDtot - oldEDtot)-(neweDihed - oldED_1) ) > 0.0001 ){
    //    cerr << "Dihedral panic, p = " << p->myId << " DE tot: " <<    fabs( (newEDtot - oldEDtot)-(neweDihed - oldED_1) ) << endl;
      //  exit(1);
  //}

  neweLJ = 0.0;
  q = newCell->firstParticle;
  do{
        if ( q != p ){
            neweLJ += LJ_pairEnergy(p, q);
        }
        q = q->next;
  } while( q );

  /* Loop over particles in the 26 neighbour cells */
  for( int i = 0; i < 26; i++){

    Cell *n = cellHash->getItemByKey( newCell->neighbours[i] );
    if( n == NULL ) continue;

    q = n->firstParticle;
    while( q ){
      neweLJ += LJ_pairEnergy(p, q);
      q = q->next;
    }
  }

  oldE =  oldEBond +  oldeDihed + oldeSAK + oldeLJ;
  newE =  newEBond +  neweDihed + neweSAK + neweLJ;
  
#ifdef TEST_RB
   double testNewEdihed    = totalDihedEnergy();
   double testNewEBond     = bond_totalEnergy();
   double testNewELJ       = LJ_totalEnergy();
  
   int    status = 0;
   
   if(  fabs( ((testNewELJ-testOldELJ)/(neweLJ-oldeLJ)) - 1.) >  1e-6 
     && fabs(  (testNewELJ-testOldELJ)-(neweLJ-oldeLJ))       > 0.001 ) {
       cerr << "LJ ERROR Particle" << endl;
       
       status = 1;
    }  
   if(( testNewEdihed   - testOldEdihed ) - (neweDihed    - oldeDihed) >   0.001 ||
      ( testNewEdihed   - testOldEdihed ) - (neweDihed    - oldeDihed) <  -0.001 ){
       
       cerr << "DH ERROR Particle" << endl;
       status = 1;
    } 
   if(( testNewEBond   - testOldEBond ) - (newEBond    - oldEBond) >   0.001 ||
      ( testNewEBond   - testOldEBond ) - (newEBond    - oldEBond) <  -0.001 ){
       
       cerr << "Bond ERROR Particle" << endl;
       
       status = 1;
    }  
    if( status ){
        cerr << "delta eDihed   : " << testNewEdihed - testOldEdihed << " versus: " << neweDihed - oldeDihed << endl;
        cerr << "delta eLJ      : " << testNewELJ    - testOldELJ    << " versus: " << neweLJ    - oldeLJ << endl;
        cerr << "delta eBond    : " << testNewEBond  - testOldEBond  << " versus: " << newEBond  - oldEBond << endl;
        
        //brute_debug_cells();
        //testCellIndices();
        
        writeOutConfig("test_particle_oops.xyz");
        
        for (int i=0; i<3; i++){
            p->R[i] = Rold[i];
        }
        writeOutConfig("test_particle_pre_oops.xyz");
        cerr << "quitting" << endl;
        
       exit( 1 );
    }
   
   
#endif     

  /* test for Metropolis Criterion */
  if (log(ran2(iran)) < (oldE - newE)){


    ePotTot   += (newE - oldE);
    eDihedTot += (neweDihed - oldeDihed);
    eBondTot  += (newEBond - oldEBond);
    eSAKTot   += (neweSAK - oldeSAK);
    eLJTot    += (neweLJ - oldeLJ);

    transAccept ++;

        //clean up the particle's old cell if it is now empty.
        if( oldCell->firstParticle == NULL ){
            cellHash->removeItem( oldCell->myId );
                delete oldCell;
        }

  }else{

       for (int i=0; i<3; i++){
            p->R[i] = Rold[i];
       }
       if ( newCell != oldCell ){

          moveBetweenCells(p, newCell, oldCell->myId, cellHash, box);
          if( newCell->firstParticle == NULL ){
              cellHash->removeItem( newCell->myId );
              delete newCell;
          }
       }
    }
}

void hsc::usage(){
    cerr << "./hsc [-n] [-f <input traj>] [-o <output traj>] -p <parameter file> -r <restart file>" << endl;
}


int main(int argc, char *argv[]){

  hsc Hsc;
  Hsc.newSetUp  = false;
  Hsc.readWells = false;
  //Hsc.testDihedral();
  //Hsc.testLJ();
  //exit(8);
  string iname;

  if (argc < 3) {
    cerr << "Too few arguments \n";
    Hsc.usage();
    exit( 8 );
  }

  //set some defaults
  Hsc.outputFileName       = string("DNA_MC.xyz");
  Hsc.EquilFileName        = string("DNA_Equil.xyz");
  Hsc.restartFileName      = string("DNA_MC_restart.xyz");
  Hsc.rotationFileName     = string("DNA_rigidrotate.xyz");
  Hsc.toporotFileName      = string("DNA_toporotate.xyz");
  Hsc.normalFileName       = string("DNA_notrotate.xyz");
  Hsc.rejectedFileName     = string("DNA_rejected.xyz");
  Hsc.normalDNAFileName    = string("DNA_normal.xyz");
  Hsc.checkpointFileName   = string("DNA_checkpoint.xyz");
  Hsc.extshrinkFileName    = string("DNA_extshrink.xyz");
  Hsc.ExtensionFileName    = string("Extension.dat");
  Hsc.ShrinkFileName       = string("Shrink.dat");
  Hsc.TopologyFileName     = string("Topology.dat");
  Hsc.Ext_forceFileName    = string("Ext_force.dat");
  Hsc.Topo_pot_ener_plusFileName   = string("Topology_Energy_diff_plus.dat");
  Hsc.Topo_pot_ener_minusFileName  = string("Topology_Energy_diff_minus.dat");
  
  //Hsc.topo_ener_infoFileName       = string("topo_ener_info_start.dat");
  Hsc.Shr_energy_accFileName               = string("Shrink_energy_acc.dat");
  Hsc.Shr_energy_rejFileName               = string("Shrink_energy_rej.dat");
    while ((argc > 1) && (argv[1][0] == '-')) {

          cout << argc << "<<" << argv[1][1] <<endl;

    switch (argv[1][1]) {
    case 'n':              // new Setup
      cerr << "Setting up a new system" << endl;
      Hsc.newSetUp=true;
      break;
    case 'f':
      Hsc.configFileName = string(&argv[1][2]);
      break;
    case 'o':
      Hsc.outputFileName = string(&argv[1][2]);
      break;
    case 'r':
      Hsc.restartFileName = string(&argv[1][2]);
      break;
    case 'p' :
      Hsc.ParameterFile= &argv[1][2];
      break;
    default:
      cerr << "Bad option " << argv[1] << '\n';
      Hsc.usage();
    }

    ++argv;
    --argc;
  }

  Hsc.run();

}
