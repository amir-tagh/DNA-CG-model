#include <stdio.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdio>
#include "FreeEnergyHS.h"
#include "tools.h"
#include "LJ.h"
#include "cell.h"
#include "hashTable.h"
#include "tools.h"
#include "wangLandau.h"
#include <assert.h>

using namespace std;

extern float ran2(long *idum);

#define NUM_DELTAS 1024
double delta_buf_plus[NUM_DELTAS];
double delta_buf_minus[NUM_DELTAS];
int    delta_i_p = 0;
int    delta_i_n = 0;

struct Matrix {
    double m[9];
};
typedef struct Matrix Matrix;


///after applying cluster move, need to rigid-body rotate particles
///which are "upstream" of the non-rigidly rotated section, providing that
///they are still inside the box.
/*
void rotateUpToBoundary(int clusterEnd_1, int clusterEnd_2, int axisDirection){

}
*/

/*change the box rotation*/
void hsc::changeBoxRotations(int new_box_turns_per_image){

    //set the global variables
    init_rotorBox(new_box_turns_per_image);
    //cerr << " new_box_turns_per_image: " << new_box_turns_per_image << endl;

    //loop over special cells on the z-boundary.
    Cell *c;
    c = cellHash->activeCells;
    //assert (c);
    while( c ){
        c->assignNeighbours(box);
        c = c->activeCellsNext;
    }
}


Matrix *matrix_myrot(double *COM_axis, double theta) {


  Matrix *matrix = new Matrix;
  double  sTheta;
  double  cTheta;
  double  ux,uy,uz;
  double  length;

  length = sqrt((COM_axis[0] * COM_axis[0]) + (COM_axis[1] * COM_axis[1]) + (COM_axis[2] * COM_axis[2]));

  ux = COM_axis[0]/length;
  uy = COM_axis[1]/length;
  uz = COM_axis[2]/length;

  sTheta = sin(theta);
  cTheta = cos(theta);

  double mcTheta = 1.0 - cTheta;

  matrix->m[0] = cTheta + (ux * ux * mcTheta);        matrix->m[1] = (ux * uy) * mcTheta - (uz * sTheta);    matrix->m[2] = ux * uz * mcTheta + (uy * sTheta);
  matrix->m[3] = uy * ux * mcTheta + (uz * sTheta);   matrix->m[4] = cTheta + uy * uy * mcTheta;             matrix->m[5] = uy * uz * mcTheta - (ux * sTheta);
  matrix->m[6] = uz * ux * mcTheta - (uy * sTheta);   matrix->m[7] = uz * uy * mcTheta + (ux * sTheta);      matrix->m[8] = cTheta + (uz * uz * mcTheta);

  return matrix;
}

inline void applyRot(double *m, double *R, double *to){

     to[0] = m[0] * R[0] + m[1] * R[1] + m[2] * R[2];
     to[1] = m[3] * R[0] + m[4] * R[1] + m[5] * R[2];
     to[2] = m[6] * R[0] + m[7] * R[1] + m[8] * R[2];
}

/*
Matrix *matrix_bprot(double *COM_axis, int clusterSize) {

  Matrix *matrix = &M;

  double  sTheta;
  double  cTheta;
  double  ux,uy,uz;
  double  length;

  length = sqrt ((COM_axis[0] * COM_axis[0]) + (COM_axis[1] * COM_axis[1]) + (COM_axis[2] * COM_axis[2]));

  ux = COM_axis[0]/length;
  uy = COM_axis[1]/length;
  uz = COM_axis[2]/length;


  double theta = 2.*M_PI/(clusterSize);

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
*/

extern float ran2(long *idum);

void hsc::topology_rotation(ofstream *topo_name, ofstream *topo_energy_plus,ofstream *topo_energy_minus){


 int           pNr_1_p, pNr_2_c, pNr_3_p, pNr_4_c ;
 unsigned int  newCell_id_1;
 int           clusterSize, clusterStart_1, clusterEnd_1, clusterStart_2, clusterEnd_2;
 double        *Rold;
 double        oldEBond,oldEDihed;
 double        oldE, newE;
 double        newEBond,newEDihed;
 Cell          **oldCell, **newCell;
 Particle      *p_1, *p_2, *p_3, *p_4, *bp[4];
 Particle      **moveParts;
 double        oldeLJ, neweLJ;
 double        COM_str[3], COM_end[3], COM_axis[3];


#define TEST_RB

#ifdef TEST_RB
   //writeOutConfig("test_TOPO_start.xyz");
   double testOldEdihed = totalDihedEnergy();
   double testOldEBond  = bond_totalEnergy();
   double testOldELJ    = LJ_totalEnergy();
#endif


 if ( N/4 < 11 ){
    cerr << "Warning, system is too small for topology-changing twist move" << endl;
    return;
 }
 /*min clustersize for bp rotation is one turn (42 particles) and max is whole DNA*/
 //clusterSize = 10 + int(ran2(iran)*((N/4)-10));
 clusterSize = (N/4);   //cluster size is set to include the whole DNA

#ifndef ROTORBOX
 /*cluster always starts with a phosphate */
 clusterStart_1   = 2 * int(ran2(iran)*((N/4)));
 clusterEnd_1     = (clusterStart_1 + 2*clusterSize - 2) % (N/2);
 clusterStart_2   = N - clusterStart_1 - 2;  //start on chain 2 has higher index than end... careful.
 clusterEnd_2     = N - clusterEnd_1   - 2;
#else
 clusterStart_1   = 0;
 clusterEnd_1     = N/2 - 2;
 clusterStart_2   = N - 2;  //start on chain 2 has higher index than end... careful.
 clusterEnd_2     = N - clusterEnd_1   - 2;
#endif


  moveParts = new Particle*[4*(clusterSize)];
  oldCell   = new     Cell*[4*(clusterSize)];
  newCell   = new     Cell*[4*(clusterSize)];
  Rold      = new   double [4*(clusterSize)*3];

  topologyMoves++;
  cerr << " topologyMoves: " << topologyMoves << endl;

  for(int i=0; i<N; i++){

     //first phosphate
     moveParts[i]  = &part[i];
     oldCell[i]    = part[i].cell;
     Rold[i*3]     = part[i].R[0];
     Rold[i*3+1]   = part[i].R[1];
     Rold[i*3+2]   = part[i].R[2];
     //part[i].clusterId = 1;
  }

  //ofstream out_oldCell_before("oldCell_before.txt",ios::app);
  //for(int i=0; i<N;i++){
      //out_oldCell_before << oldCell[i]->myId << endl;
  //}
  //out_oldCell_before.close();

  //bond energies
   oldEBond  = bond_totalEnergy();//bond_topology(clusterStart_1, clusterStart_2, clusterEnd_1, clusterEnd_2);

  //dihedral energies
   oldEDihed = dihed_topology(clusterStart_1, clusterStart_2, clusterEnd_1, clusterEnd_2);

  //LJ energies
   oldeLJ = LJ_totalEnergy();//LJ_topology(clusterStart_1, clusterSize);
   //NeighbourList(clusterStart_1, clusterSize);//debug function

   /*this section selects the center of mass and axis of rotation*/
    COM_str[0] = 0.5 * (part[clusterStart_1].R[0] +
                        part[clusterStart_2].R[0]);

    COM_str[1] = 0.5 * (part[clusterStart_1].R[1] +
                        part[clusterStart_2].R[1]);

    COM_str[2] = 0.5 * (part[clusterStart_1].R[2] +
                        part[clusterStart_2].R[2]);

    COM_end[0] = 0.5 * (part[clusterEnd_1].R[0] +
                        part[clusterEnd_2].R[0]);

    COM_end[1] = 0.5 * (part[clusterEnd_1].R[1] +
                        part[clusterEnd_2].R[1]);

    COM_end[2] = 0.5 * (part[clusterEnd_1].R[2] +
                        part[clusterEnd_2].R[2]);


     if (COM_end[2] < COM_str[2]){  //z-axis positive
         COM_end[2] += box->x[2];
     }

     /*define the axis direction to increase or decrease DNA twist*/
     //int axis_direction;
     //axis_direction  =  int(ran2(iran)*2.0);

     if(wL->direction_twist == 1){             /*positive rotation...increase twist*/
        COM_axis[0] =  COM_end[0] - COM_str[0];
        COM_axis[1] =  COM_end[1] - COM_str[1];
        COM_axis[2] =  COM_end[2] - COM_str[2];
     }else{                                    /*negative rotation...negative twist*/
        COM_axis[0] =  COM_str[0] - COM_end[0];
        COM_axis[1] =  COM_str[1] - COM_end[1];
        COM_axis[2] =  COM_str[2] - COM_end[2];
     }


#ifdef ROTORBOX
  int old_box_turns_per_image = box_turns_per_image;
  //cerr << " HERE old_box_turns_per_imag: " << old_box_turns_per_image << endl;
  //cerr << " TOPO CALL INIT BOX FUNCTION " << endl;
  if( wL->direction_twist == 1 )
    changeBoxRotations((box_turns_per_image + stepsize_box_turns_per_image + 444) % 4);   //adjust the box neighbouring cells to NEW value
  else
    changeBoxRotations((box_turns_per_image - stepsize_box_turns_per_image + 444) % 4);
#endif

  /*get the new position of the particles after rotation*/
     Matrix  *m;
     double  NewPos_1[3];
     double  ref_p[3];
     double  theta, dtheta;

     ref_p[0] = moveParts[0]->R[0];
     ref_p[1] = moveParts[0]->R[1];
     ref_p[2] = moveParts[0]->R[2];

     /*check the rotation angle for doing half-integer rotations in case the box is rotated as well*/
     /*for stepsize = 2 DNA is rotaed from -90 to +90 degree*/
     /*for stepsize = 4 DNA is rotated from -180 to +180 degree.*/
     dtheta = (stepsize_box_turns_per_image * 0.5 * M_PI);

#ifdef ROTORBOX

     for( int i=0; i<N; i++ ){

           m = matrix_myrot(COM_axis, (part[i].R[2]/box->x[2])*dtheta );
           //cerr << " part: " << i << " Z: " << part[i].R[2] << " theta: " << (part[i].R[2]/box->x[2])*dtheta*180/M_PI << endl;

           applyRot(m->m, part[i].R, NewPos_1);
           part[i].R[0] = NewPos_1[0];
           part[i].R[1] = NewPos_1[1];
           part[i].R[2] = NewPos_1[2];
           intoBox(part[i].R, box);

           newCell_id_1 = cellIndex(part[i].R, box);

           if(newCell_id_1 != oldCell[i]->myId ){
                   newCell[i] = moveBetweenCells(&part[i], oldCell[i], newCell_id_1, cellHash, box);
            }else{
                   newCell[i] = oldCell[i];
           }

           delete m;
     }

#endif

//cellHash->getNumberOfItems();
#if 0
if(topologyMoves % 10 == 0){
writeOutConfig("toplogy_config.xyz");
}
#endif


#ifdef ROTORBOX
if( topologyMoves % 100 == 0 ){
      ofstream myfile;
      myfile.open(toporotFileName.c_str(), ios::app);
      writeXYZ_frame( &myfile );
      myfile.close();
}
#endif

  //New Bond energies
  newEBond  = bond_totalEnergy();//bond_topology(clusterStart_1, clusterStart_2, clusterEnd_1, clusterEnd_2);
  newEDihed = dihed_topology(clusterStart_1, clusterStart_2, clusterEnd_1, clusterEnd_2);
  neweLJ    = LJ_totalEnergy();//LJ_topology(clusterStart_1, clusterSize);
  //NeighbourList(clusterStart_1, clusterSize);//debug function

  oldE =  oldEBond +  oldEDihed + oldeLJ;
  newE =  newEBond +  newEDihed + neweLJ;
  //cerr << " oldEBond " << oldEBond << " oldEDihed " << oldEDihed << " oldeLJ " << oldeLJ <<  endl;
  //cerr << " newEBond " << newEBond << " newEDihed " << newEDihed << " neweLJ " << neweLJ << endl;

#ifdef TEST_RB

   double testNewEdihed    = totalDihedEnergy();
   double testNewEBond     = bond_totalEnergy();
   double testNewELJ       = LJ_totalEnergy();
   int status = 0;

   if(  fabs( ((testNewELJ-testOldELJ)/(neweLJ-oldeLJ)) - 1.) >  1e-6
     && fabs(  (testNewELJ-testOldELJ)-(neweLJ-oldeLJ))       > 0.001 ) {
       cerr << "LJ ERROR TOPO" << endl;
       status = 1;


   }
   if(( testNewEdihed  - testOldEdihed ) - (newEDihed - oldEDihed) >   0.0001 ||
       ( testNewEdihed - testOldEdihed ) - (newEDihed - oldEDihed) <  -0.0001 ){

       cerr << "Dihedral ERROR TOPO" << endl;
       status = 1;
    }
   if(( testNewEBond - testOldEBond ) - (newEBond  - oldEBond) >   0.0001 ||
      ( testNewEBond - testOldEBond ) - (newEBond  - oldEBond) <  -0.0001 ){

       cerr << "Bond ERROR TOPO" << endl;
       status = 1;
    }
   if( status ){
        cerr << "delta eDihed   : " << testNewEdihed - testOldEdihed  << " versus: " << newEDihed - oldEDihed << endl;
        cerr << "delta eLJ      : " << testNewELJ    - testOldELJ     << " versus: " << neweLJ    - oldeLJ << endl;
        cerr << "delta eBond    : " << testNewEBond  - testOldEBond   << " versus: " << newEBond  - oldEBond << endl;

        //brute_debug_cells();
        testCellIndices();

        //writeOutConfig("test_TOPO_oops.xyz");
        exit( 1 );
    }
#endif

double deltaE = newE - oldE;
/*collecting data based on the potential energy difference*/
    if(wL->direction_twist == 1){ //increased twist

        delta_buf_plus[delta_i_p++] = deltaE;
        if( delta_i_p >= NUM_DELTAS ){

            for( int i = 0; i < NUM_DELTAS - 1; i++ ){

            *topo_energy_plus << delta_buf_plus[i] << "\n"; //newline does not flush to disk
        }
         *topo_energy_plus << delta_buf_plus[NUM_DELTAS-1] << endl; //endl flushes write to disk.
         delta_i_p = 0;
        }
    }else{ /*decreased twist*/

        delta_buf_minus[delta_i_n++] = deltaE;
        if( delta_i_n >= NUM_DELTAS ){

            for( int i = 0; i < NUM_DELTAS - 1; i++ ){

                *topo_energy_minus << delta_buf_minus[i] << "\n"; //newline does not flush to disk
             }
         *topo_energy_minus << delta_buf_minus[NUM_DELTAS-1] << endl; //endl flushes write to disk.
         delta_i_n = 0;
        }
    }

    /*wang landau*/
    double  x_twist, y_stretch;
    double  x_new_twist;
    double  delta_hist;
    double  log_accept;
    int     i, j, ii, jj;
    bool    productionStage;



    x_twist   = wL->bin_min_twist_update;
    y_stretch = wL->bin_min_stretch_update;
    //cerr << " TOPO " << " bin_min_stretch_update: " << wL->bin_min_stretch_update << endl;
    productionStage = wL->productionStage;


   if (productionStage == true){
       x_new_twist   = x_twist + wL->stepsize_twist;
   }else{
       x_new_twist   = x_twist + wL->stepsize_twist * wL->direction_twist;
       //cerr << " x_new_twist: " << x_new_twist << " x_twist: " << x_twist << " stepsize_twist: " <<  wL->stepsize_twist << endl;
   }


   if (x_new_twist >= wL->bin_max_twist){
        wL->direction_twist *= -1;
        wL->getBin_2dhist(x_twist, y_stretch, &i, &j);
        count_hist[i][j] += 1;

   }else if(x_new_twist < wL->bin_min_twist){
        wL->direction_twist *= -1;
        wL->getBin_2dhist(x_twist, y_stretch, &i, &j);
        count_hist[i][j] += 1;

   }else{
        wL->getBin_2dhist(x_twist, y_stretch, &i, &j);
        count_hist[i][j] += 1;
        //cerr << " TOPO count_hist[i][j]: " << count_hist[i][j] << " i: " << i << " j: " << j << endl;
        //cerr << " i: " << i << " j: " << j << endl;
   }

   //cerr << " row: " << wL->n_stretch << " column: " << wL->n_twist << endl;

   ofstream out("VISIT.TXT", ios::app);
   if(topologyMoves % 1000 == 0){

       for(int w = 0; w < wL->n_stretch; w++){
           for(int x = 0; x < wL->n_twist; x++){
               //cerr << " VISITTT " << count_hist[w][x] << endl;
               out << count_hist[w][x];
               out << " , ";
           }
           out << "\n";
       }
   }
     out.close();
     //exit(1);

   wL->getBin_2dhist(x_new_twist, y_stretch, &ii, &jj);
   //cerr << " HERE TOPO bias_hist[ii][jj]: " << bias_hist[ii][jj] << " ii: " << ii << " jj: " << jj << endl;
   delta_hist = bias_hist[ii][jj] - bias_hist[i][j];
   //cerr << " TOPOdelta_hist: " << delta_hist << endl;
   log_accept = delta_hist - deltaE;
   //cerr << " log_accept: " << log_accept << endl;

  /*metropolis criteria based on the modifiacation factor of Wang-Landau*/

  if (log(ran2(iran)) < (log_accept)){   /*move is accepted*/

        ePotTot   += (newE - oldE);
        eDihedTot += (newEDihed - oldEDihed);
        eBondTot  += (newEBond - oldEBond);
        eLJTot    += (neweLJ - oldeLJ);

        topology_accept++;
        cerr << " topology_accept (top_rotor): " << topology_accept << endl;

        /*flat histogram*/
        wL->bin_min_twist_update = x_new_twist;


    /*counter for collecting data(energy differences)*/
    if(wL->direction_twist == 1){
            topoAccept_positive++;
    }else{
            topoAccept_negative++;
    }

  }else{       /*move is rejected*/

        topology_reject++;
        cerr << " topology_reject " << topology_reject << endl;

       /*flat histogram*/
       if((ii != i || jj != j) && productionStage == false ){
             bias_hist[ii][jj] += wL->hist_increment;
             //cerr << " TOPO bias_hist[ii][jj] : " << bias_hist[ii][jj] << " ii: " << ii << " jj: " << jj << " i: " << i
             //<< " j: " << j << endl;
       }

       if(wL->direction_twist == 1){
            topoReject_positive ++;
       }else{
            topoReject_negative ++;
        }


    //collecting accept or reject counts for the transition matrix.
    if(topologyMoves % 1000 == 0){
      *topo_name << " topopositive accept: " << topoAccept_positive << " reject: " << topoReject_positive << " ";
      *topo_name << " toponegative accept: " << topoAccept_negative << " reject: " << topoReject_negative << endl;
    }

#ifdef ROTORBOX
  cerr << " old_box_turns_per_image: " << old_box_turns_per_image << endl;
    changeBoxRotations(old_box_turns_per_image); /*return back the cell neighbouring system after doing the rotation*/

#endif

    /*debug to test the cell indices*/
    cerr << " TOPO INC " << endl;
    testCellIndices();

    //ofstream out_oldCell_after("oldCell_after.txt");
    //for(int i=0; i<N;i++){
        //out_oldCell_after << oldCell[i]->myId << endl;
   //}
    //out_oldCell_after.close();


       for(int i=0; i<N; i++) {

           p_1 = moveParts[i];
           p_1->R[0] =  Rold[i*3];
           p_1->R[1] =  Rold[i*3+1];
           p_1->R[2] =  Rold[i*3+2];
           //p_1->clusterId  = 0;


           if( newCell[i] != oldCell[i]) {
             moveBetweenCells(p_1, newCell[i], oldCell[i]->myId, cellHash, box);
           }
       }

       for(int i = 0; i < N; i++){
              if( newCell[i]->firstParticle == NULL )
                  if( !cellHash->removeItem(newCell[i]->myId) )
                    newCell[i]       = NULL;
       }
       for(int i = 0; i < N; i++){
              if( newCell[i] != NULL && newCell[i]->firstParticle == NULL)
                  delete newCell[i];
       }



   }
   //cerr << " end of topo " << endl;

  //cellHash->getNumberOfItems();
  delete [] oldCell;
  delete [] newCell;
  delete [] Rold;
  delete [] moveParts;
  //cerr << " delete done " << endl;
}
