#include <stdio.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <sstream>
#include <string>
#include <assert.h>
#include "hashTable.h"
#include "cell.h"
#include "FreeEnergyHS.h"
#include "wangLandau.h"

class Cell;
//class Particle;

using namespace std;


#define _DEBUG_COLLISIONS
#ifdef  DEBUG_COLLISIONS
long int queries, collisions;
long int itemsInHashup, itemsInHashdown;
#endif


// Constructs the empty Hash Table object.
HashTable::HashTable()
{

    //pick a string of some decent length.
    char salt_string[] = "If you want to get a PhD you are going to have to dramatically increase the speed at which you work.  Try learning to type, for instance.";

    array  = new Cell*[ HASH_LENGTH ];

    itemsInHash = 0;
    
    for(unsigned int i = 0; i < HASH_LENGTH; i++)
        array[i] = NULL;

    //generate a random-seeming int from the salt phrase.
    salt = 5381;
    char *pt_c = &salt_string[0];
    while(*pt_c){
        unsigned int c = (unsigned int)(*pt_c++);
        salt = ((salt << 5) + salt) + c;
    }

    activeCells = NULL;
    allCells    = NULL;
    
//    printf("#Hash      salt: %ux\n", salt);
//    printf("#Effective salt: %ux\n", salt & HASH_LMASK);
//    printf("#Hash table size: %u\n", HASH_LENGTH);

}

/*Returns an array location for a given item key.*/
inline unsigned int HashTable::hashB( unsigned int itemKey )
{
    unsigned int value;

    //Scrambling the index to give pseudo-random mapping.
    //value = (itemKey ^ salt) % length; //for a prime number
    value = (itemKey ^ salt) & (HASH_LMASK); //for a number of 2^n

    assert( value < HASH_LENGTH );
    return(value);
}

 Cell *HashTable::getItemByKey( unsigned int itemKey )
{
        
    unsigned int index = hashB( itemKey );

#ifdef  DEBUG_COLLISIONS
    queries++;
#endif
    
    Cell *c = array[index];
    while( c != NULL ){
        if( c->myId == itemKey ) break;
#ifdef  DEBUG_COLLISIONS
        collisions++;
#endif
        c = c->hashBucketNext;
    }

    return( c );
        
    };


/*Adds an item to the Hash Table.*/
void HashTable::insertItem( Cell* newItem )
{
    unsigned int index;
    index = hashB( newItem->myId );
    
    //insert at head of list
    newItem->allCellsNext = allCells;
    if( allCells != NULL ){
        allCells->allCellsPrev = newItem;
    }
    allCells = newItem;
    
    if( !array[index] ){
     
        array[index] = newItem;
        newItem->hashBucketPrev = NULL;
        newItem->hashBucketNext = NULL;

    }else{
        
        array[index]->hashBucketPrev = newItem;
        newItem->hashBucketNext      = array[index];
        newItem->hashBucketPrev      = NULL;
        array[index]                 = newItem;
    }
    
    /*add to the head of the activeCells list *IF* this is a special cell on the z-boundary.*/
    if( newItem->zBound != 0 ){
        newItem->activeCellsNext     = activeCells;
        //assert( activeCells );
        if( activeCells )
            activeCells->activeCellsPrev = newItem;
        activeCells                  = newItem;
    }
        
#ifdef DEBUG_COLLISIONS    
    itemsInHashup++;
#endif  
}



// Deletes an Item from the Hash Table.
// Returns true if the operation is successful.
void HashTable::removeItem_byAddress( Cell *c )
{
    
    unsigned int index;
    
    index = hashB(c->myId);
    
    //remove from hash bucket list
    if( c->hashBucketPrev ) c->hashBucketPrev->hashBucketNext = c->hashBucketNext;
    else                    array[index]  = c->hashBucketNext;//will be null if bucket is now empty
    if( c->hashBucketNext ) c->hashBucketNext->hashBucketPrev = c->hashBucketPrev;
    
    //remove from the activeCells list of all cells
    if( c->zBound != 0 ){
        if( c->activeCellsPrev ) c->activeCellsPrev->activeCellsNext = c->activeCellsNext;
        else                     activeCells        = c->activeCellsNext;
        if( c->activeCellsNext ) c->activeCellsNext->activeCellsPrev = c->activeCellsPrev;
    }
    
    //remove from all cells list.
    if( c->allCellsPrev )        c->allCellsPrev->allCellsNext       = c->allCellsNext;
    else                         allCells           = c->allCellsNext;
    if( c->allCellsNext )        c->allCellsNext->allCellsPrev       = c->allCellsPrev;
    
    
#ifdef DEBUG_COLLISIONS
    itemsInHashdown--;
#endif
}

// Deletes an Item by key from the Hash Table.
// Returns true if the operation is successful.
bool HashTable::removeItem( unsigned int itemKey )
{
    unsigned int index = hashB(itemKey);

    Cell *c = array[index];
    while( c ){
        if( c->myId == itemKey ) break; 
       
        c = c->hashBucketNext;
    }
    
    if( !c ) return(false);
    
    //remove from hash bucket list
    if( c->hashBucketPrev ) c->hashBucketPrev->hashBucketNext = c->hashBucketNext;
    else                    array[index]  = c->hashBucketNext;//will be null if bucket is now empty
    if( c->hashBucketNext ) c->hashBucketNext->hashBucketPrev = c->hashBucketPrev;
    
    
    //remove from the activeCells list of all cells
    if( c->zBound != 0 ){
        if( c->activeCellsPrev ) c->activeCellsPrev->activeCellsNext = c->activeCellsNext;
        else                     activeCells        = c->activeCellsNext;
        if( c->activeCellsNext ) c->activeCellsNext->activeCellsPrev = c->activeCellsPrev;
    }
    
    //remove from all cells list.
    if( c->allCellsPrev )        c->allCellsPrev->allCellsNext       = c->allCellsNext;
    else                         allCells           = c->allCellsNext;
    if( c->allCellsNext )        c->allCellsNext->allCellsPrev       = c->allCellsPrev;
    
    
#ifdef DEBUG_COLLISIONS
    itemsInHashdown--;
    
#endif
 
    return(true);
}


// Display the contents of the Hash Table to console window.
void HashTable::printTable()
{
    cout << "\n\nHash Table:\n";
}

// Prints a histogram illustrating the Item distribution.
void HashTable::printHistogram()
{
    cout << "\n\nHash Table Contains errrr \n\n\n";
}

// Returns the number of locations in the Hash Table.
int HashTable::getLength()
{
    cerr << "errr" << endl;
    return( -1);
}

// Returns the number of Items (cells) in the Hash Table.

int HashTable::getNumberOfItems(){
    
    FILE *HashItems;
    HashItems = fopen("hashitems.dat", "a" );
    int cellCount = 0;
    int particleCount = 0;
    int cellPartCount = 0;
    Cell *c;
    Particle *p;
    
    
    //cerr << " loop start " << endl;
    for ( unsigned int i = 0; i < HASH_LENGTH; i++ )
    {
      c = array[i]; 
     
      if ( !c )
          continue;
       
        while( c ) {
           
            cellCount++;
            p = c->firstParticle;
            while(p){
               
                particleCount++;
                cellPartCount++;
                p = p->next;
            }
            
            fprintf(HashItems, "cells id %d cells %d particles %d\n", c->myId, cellCount, particleCount);
           
            c = c->hashBucketNext;
        }
    }
     fprintf(HashItems, "total cells %d total part %d\n", cellCount, cellPartCount);
     fclose(HashItems);
     
        return cellCount;
}
 

// De-allocates all memory used for the Hash Table.
HashTable::~HashTable()
{
    delete [] array;

}

//*****************************************************************
// End of File
//*****************************************************************
