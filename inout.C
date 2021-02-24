#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "FreeEnergyHS.h"
#include "tools.h"


using namespace std;


void hsc::writeOutConfig(char const name[64]){

  //char temp[16];
  //char allname[50];

  //strcpy(allname, name);
  //strcat(allname, "_");
  //  sprintf(temp,"%.*lf", 3, N/box->V);
  //  strcat(allname, temp);
   //strcat(allname, "_");
  //sprintf(temp,"%.*lf", 2, epsilon);
  //strcat(allname, temp);

  ofstream outfile(name, ios::out);
  outfile.precision(10);

  outfile << N << endl;
  //outfile << "#Density " << N/box->V << endl;
  outfile << box->x[0] << " " <<
             box->x[1] << " " <<
             box->x[2] << endl;

  for (int i=0; i<N; i++){

    if( part[i].beadType == 'P' ){
       outfile << "P ";
    }else{
#ifdef INTERCAL
       if( part[i].interCal ){
          outfile << "O ";
       }
#else
       outfile << "C ";
#endif
    }
   outfile  << part[i].R[0] << " "
            << part[i].R[1] << " "
            << part[i].R[2] << "\n";
  }

  outfile.close();
}

void hsc::readXYZ(string name) {

  char lineIn[128], atName[8];


  ifstream infile(name.c_str(), ios::in);

  /* read N */
  infile.getline(lineIn,sizeof(lineIn));
  cout << lineIn << " " << endl;
  if( 1 != sscanf(lineIn, "%i", &N) ){
        cerr << "Error failed to read N atoms from XYZ file: " << name << endl;
        exit( 8 );
  }
  part = new Particle[N];


  /* read the comment line with box coordinates */
  infile.getline(lineIn,sizeof(lineIn));
  box = new Box;
  if( 3 != sscanf( lineIn, "%lf %lf %lf\n",
                   &box->x[0], &box->x[1], &box->x[2])){
        cerr << "Error failed to read coords from file: " << name << " box x: " << box->x[0] << endl;
        exit( 8 );
   }


  /* read the particles */
  for (int i=0; i<N; i++){
    if( infile.peek() == EOF )
      cerr << "Error: end of file at " << i << endl;

    infile.getline(lineIn,sizeof(lineIn));
    if( 4 != sscanf( lineIn, "%s %lf %lf %lf\n",
               atName, &part[i].R[0], &part[i].R[1], &part[i].R[2])){
            cerr << "Error failed to read coords from file: " << name << endl;
        exit( 8 );
    }
    if( atName[0] == 'P' ){
      part[i].beadType = 'P';
    }else if ( atName[0] == 'C' || atName[0] == 'O' ){
      part[i].beadType = 'B';
    }else{
        cerr << "Error did not recognise atom name in line: " << lineIn << endl;
        exit( 8 );
    }
   }

}
void hsc::readXYZ(const char* name) {

  char lineIn[128], atName[8];

  ifstream infile(name, ios::in);

  // read N
  infile.getline(lineIn,sizeof(lineIn));
  if( 1 != sscanf(lineIn, "%i", &N) ){
        cerr << "Error failed to read N atoms from XYZ file: " << name << endl;
        exit( 8 );
  }

  //This is corrected just to get x & y of the box as z is changed in Extension move
  // read the comment line with box coordinates
  infile.getline(lineIn,sizeof(lineIn));
  if( 3 != sscanf( lineIn, "%lf %lf %lf\n",
                   &box->x[0], &box->x[1], &box->x[2])){
        cerr << "Error failed to read coords from file: " << name << endl;
        exit( 8 );
   }

  // read the particles
  for (int i=0; i<N; i++){

    infile.getline(lineIn,sizeof(lineIn));

    if( 4 != sscanf( lineIn, "%s %lf %lf %lf\n",
               atName, &part[i].R[0], &part[i].R[1], &part[i].R[2])){
        cerr << "Error failed to read coords from file: " << name << endl;
        exit( 8 );
    }


    if( atName[0] == 'P' ){
      //part[i].beadType = 'P';
      cout << " part: " << part[i].beadType << endl;
    }else if ( atName[0] == 'C'  || atName[0] == 'O' ){
      part[i].beadType = 'B';
    }else{
        cerr << "Error did not recognise atom name in line: " << lineIn << endl;
        exit( 8 );
    }
  }
}

void hsc::writeXYZ_frame(ofstream *outfile){

   char outEnergy[64];
   outfile->precision(6);

  *outfile << N << endl;
  *outfile << box->x[0] << " " << box->x[1] << " " << box->x[2];

   sprintf(outEnergy," %.5f", ePotTot);  //energy time series
  *outfile << outEnergy << endl;

  for (int i=0; i<N; i++){

    if( part[i].beadType == 'P' ){
       *outfile << "P ";
    }else{
#ifdef INTERCAL
       if(part[i].interCal)
           *outfile << "O ";
       else
#endif
       *outfile << "C ";
    }

    if( i < N - 1 )
      *outfile        << part[i].R[0]
               << " " << part[i].R[1]
               << " " << part[i].R[2] <<  "\n";
    else
      *outfile        << part[i].R[0]
               << " " << part[i].R[1]
               << " " << part[i].R[2] <<  endl;  //last line flushes./

  }


}
void hsc::writeXYZ(char *name){

   ofstream outfile;
   outfile.open(name, ios::out);
   outfile.precision(6);

   outfile << N << endl;
   outfile << box->x[0] << " " << box->x[1] << " " << box->x[2];

  for (int i=0; i<N; i++){
    if( part[i].beadType == 'P' ){
       outfile << "P ";
    }else{
#ifdef INTERCAL
       if(part[i].interCal)
           outfile << "O ";
       else
#endif
       outfile << "C ";
    }
   outfile        << part[i].R[0]
                  << " " << part[i].R[1]
                  << " " << part[i].R[2] <<  "\n";
  }

  outfile.close();
}



void hsc::writeOutSimData(char const name[16]){

  //char temp[16];
  char allname[50];

  strcpy(allname, name);
  //strcat(allname, "_");
  //sprintf(temp,"%.*lf", 3, N/box->V);
  //strcat(allname, temp);
  //strcat(allname, "_");
  //sprintf(temp,"%.*lf", 2, epsilon);
  //strcat(allname, temp);

  ofstream outfile(allname, ios::out);

  outfile << "# N: " << N << endl;
  outfile << "# Box x, y, z: " << box->x[0] << " " << box->x[1]
	                           << " " << box->x[2] << endl;
  //outfile << "# Epsilon: " << epsilon << endl;
  outfile << "# Model: linear" << endl;
  outfile << "# Cutoff radius: " << rcut << endl;
  outfile << "# Translation acceptance rate: "
	      << float(transAccept)/float(transMoves) << endl;
  outfile << "# Rigid rotation acceptance rate: "
          << float(rigidAccept)/float(rigidMoves) << endl;
  outfile << "#Max Displacement: " <<    maxDisplace << endl;
  outfile << "#Rigid rotation: " << thetaRigid << endl;

  
  outfile.close();
}
