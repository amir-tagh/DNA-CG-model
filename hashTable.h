//*****************************************************************
//  HashTable.h
//  HashTable
//
//  This header file contains the Hash Table class declaration.
//  Hash Table array elements consist of Linked List objects.
//  Hash Length is defined based on the size of DNA.
//*****************************************************************

#ifndef HashTable_h
#define HashTable_h

#include "cell.h"

#define HASH_LENGTH  0x4000000 //(hexadecimal notation of 2^26)
#define HASH_LMASK (HASH_LENGTH - 1) 
#define _DEBUG_COLLISIONS
#ifdef  DEBUG_COLLISIONS
extern long int queries, collisions;
extern long int itemsInHashup, itemsInHashdown;
#endif

//*****************************************************************
// Hash Table objects store a fixed number of Linked Lists.
//*****************************************************************
class HashTable
{
private:

    // Array is a reference to an array of Linked Lists.
    Cell **array;
    
    //salt is an int's worth of random bits to conflugulate
    //the keys.
    unsigned int salt;

    // Returns an array location for a given item key.
    unsigned int hashB( unsigned int itemKey );

    long int itemsInHash;
    
    
public:
    //list of cells on the special z boundary
    Cell  *activeCells;
    
    //list of inserted cells.
    Cell  *allCells;

    // Constructs the empty Hash Table object.
    HashTable();

    // Adds an item to the Hash Table.
    void insertItem( Cell* newItem );

    // Deletes an Item by key from the Hash Table.
    // Returns true if the operation is successful.
    bool removeItem( unsigned int itemKey );

    // Deletes an item from the hash.
    void removeItem_byAddress( Cell *c );
    
    // Returns an item from the Hash Table by key.
    // If the item isn't found, a null pointer is returned.
    Cell* getItemByKey( unsigned int itemKey );
    
    // Display the contents of the Hash Table to console window.
    void printTable();

    // Prints a histogram illustrating the Item distribution.
    void printHistogram();

    // Returns the number of locations in the Hash Table.
    int getLength();

    // Returns the number of Items in the Hash Table.
    int getNumberOfItems();
    
    //prints the number of particles in each cell
    int PrintParticleInCell();

    // De-allocates all memory used for the Hash Table.
    ~HashTable();
};

// Returns an item from the Hash Table by key.
// If the item isn't found, a null pointer is returned.
/*inline Cell *HashTable::getItemByKey( unsigned int itemKey )
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
}*/
#endif

//*****************************************************************
// End of File
//*****************************************************************










