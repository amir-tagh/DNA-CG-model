//debug code for neighbors in LJ
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <sstream>
#include <string>
#include "FreeEnergyHS.h"

#include "tools.h"

int hsc::testNeb(int *nebs, int nnebs, int id, int id2){
               
    for( int nn = 0; nn < nnebs; nn++){
        if( nebs[nn] == id2 ){
               return( nn );
               break;
        }
    }
    
    double dr2, dr[3];
    minVec(part[id].R, part[id2].R, dr, box->x, box->halfx);
    dr2  = dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2];
    if( dr2 > rcut2 ){
        return( -2 );
    }
    
    return( -1 );
}


int hsc::reportNebs(int *nebs, int nnebs, int id1, int id2){
           
    double   dr2, dr[3];
               
    cerr << "nebs of " << id1 << endl;
    for( int nn = 0; nn < nnebs; nn++){
        cerr << "   " << nebs[nn] << endl; 
    }
    
    cerr << "distance id1 id2:" << endl;
    minVec(part[id1].R, part[id2].R, dr, box->x, box->halfx);
    
    dr2  = dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2];
    
    cerr << sqrt(dr2) << endl;
    
    return( -1 );
}

void hsc::brute_debug_cells(){
    
    cerr << "Expensive check of cell list" << endl;
    int     *nebs, *found, nnebs;
    double   dr2, dr[3];
    Particle *p,*q;
    
    rcut2 = rcut * rcut;
    nebs  = new int[512];
    found = new int[512];
    
    for(int i = 0; i < N; i++ ){
        
        nnebs = 0;
        for (int j = i + 1 ; j < N; j++ ){
            minVec(part[i].R, part[j].R, dr, box->x, box->halfx);
            
            dr2 = dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2];
            
            if( dr2 < rcut2 ){
                found[nnebs]  = 0;
                nebs[nnebs++] = j;
                
            }
        }
        
        p = &part[i];
        q = p->cell->firstParticle;
        do{
            if ( q->myId > p->myId){
                int n_found;
                n_found = testNeb(nebs, nnebs, p->myId, q->myId);
                if( -1 == n_found ) {
                    cerr << "Did not find neighbour. " << q->myId << endl;
                    reportNebs(nebs, nnebs, p->myId, q->myId);
                    exit( 1 );
                }else{
                    found[n_found] = 1;
                }
            }
            q = q->next;
        } while( q );
        
        /* Loop over particles in the 26 neighbour cells */
        for( int i = 0; i < 26; i++){
            Cell *c = cellHash->getItemByKey(p->cell->neighbours[i]);
            if( c == NULL ) continue;

            q = c->firstParticle; //to update the cell if necessary
            while( q ){
                if ( q->myId > p->myId){
                    int n_found;
                    n_found = testNeb(nebs, nnebs, p->myId, q->myId);
                    if( -1 == n_found ) {
                        cerr << "Did not find neighbour. " << q->myId << endl;
                        reportNebs(nebs, nnebs, p->myId, q->myId);
                        exit( 1 );
                    }else{
                        found[n_found] = 1;
                    }
                }
                q = q->next;
            }   
        }
        
        for(int nn = 0; nn < nnebs; nn++){
            if( found[nn] != 1 ){
                cerr << "missed one" << endl;
                exit(1);
            }
        }
    }
}