#include "FreeEnergyHS.h"
#include "cell.h"


class Particle;

#include "tools.h"

void Cell::assignNeighbours(Box *box){
    
      //get the indices of the neighbouring cells in 3D.
      int xx,yy,zz;
      int c = 0;
      
      for(int i = x-1; i<=x+1; i++){
          xx = (i + box->nxCell[0]) % box->nxCell[0];
          for(int j = y-1; j<=y+1; j++){
              yy = (j + box->nxCell[1]) % box->nxCell[1];
              for(int k = z-1; k<=z+1; k++){
                  if( k==z && j==y && i==x ) continue;
                  zz = (k + box->nxCell[2]) % box->nxCell[2];
#ifdef ROTORBOX
      int xx_rotated = xx; //default remains unchanged in most cases
      int yy_rotated = yy; //.
      
      if( z == 0 || z == box->nxCell[2] - 1 ){
        switch ( box_turns_per_image ){
            
            case 1:
                if( z == 0 && zz == box->nxCell[2] - 1){   //-90 degrees
                    xx_rotated = yy;
                    yy_rotated = box->nxCell[0] - xx - 1;
                }else if( z == box->nxCell[2] - 1 && zz == 0){
                    xx_rotated = box->nxCell[1] - yy - 1; //+90 degrees
                    yy_rotated = xx;
                }
            break;
            case 2:
                //for 180 degree rotation, up and down cases are the same.
                if( z == 0 && zz == box->nxCell[2] - 1){ 
                    xx_rotated = box->nxCell[0] - xx - 1;
                    yy_rotated = box->nxCell[1] - yy - 1;
                }else if( z == box->nxCell[2] - 1 && zz == 0){
                    xx_rotated = box->nxCell[0] - xx - 1;
                    yy_rotated = box->nxCell[1] - yy - 1;
                }
            break;
            case 3:
                if( z == 0 && zz == box->nxCell[2] - 1){   //+90 degrees
                    xx_rotated = box->nxCell[1] - yy - 1;;
                    yy_rotated = xx;
                }else if( z == box->nxCell[2] - 1 && zz == 0){
                    xx_rotated = yy; //-90 degrees
                    yy_rotated = box->nxCell[0] - xx - 1;
                }
            break;
            default:
            break;
        }
      }
      
         
      neighbours[c++] = (unsigned int)(xx_rotated + yy_rotated*box->nxCell[0] + zz*box->nxCell[0]*box->nxCell[1]); 

#else
      neighbours[c++] = (unsigned int)(xx + yy*box->nxCell[0] + zz*box->nxCell[0]*box->nxCell[1]);

#endif
              }
          }
      }
    
}


/*find neighbouring cells in the list.*/
Cell::Cell(unsigned int icell, Box *box){

      myId            = icell;
      firstParticle   = NULL;
      hashBucketPrev  = NULL;
      hashBucketNext  = NULL;
      activeCellsNext = NULL;
      activeCellsPrev = NULL;
      allCellsNext    = NULL;
      allCellsPrev    = NULL;
      zBound          = 0;
      
      
      x =  icell % box->nxCell[0];
      y = (icell / box->nxCell[0]) %  box->nxCell[1];
      z =  icell / (box->nxCell[0] *  box->nxCell[1]);
      
      if( z == 0 ){
            zBound = -1;
      }else if( z == box->nxCell[2] - 1 ){
            zBound =  1;
      }
      
      assignNeighbours(box);
}

/*debug function*/
void hsc::testCellIndices(){
    
    unsigned int cell_i;
    
    for( int i = 0; i < N; i++ ) {
        cell_i = cellIndex( part[i].R, box );
        
        if( cell_i != part[i].cell->myId ){
            cerr << "Cell Error Here!" << " cell index: " << cell_i << " part,Cell " << part[i].cell->myId << endl;
            exit(1);
        }
    }       
} 




