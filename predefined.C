#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <typeinfo>

using namespace std;

#define smart 1   
#define plus -0.125   
/*
class foo
{
public:
    struct packet
    {
        int x;
        int y;
    };
    packet my_packet;
};

int main(void)
{
    foo myfoo;
    
    printf("value of x is: %d", myfoo.my_packet.x);
    printf("value of y is: %d", myfoo.my_packet.y);
}
*/
void getBinIndex(double twist_update, double stretch_update, int *pt_i, int *pt_j){
        
     int i,j;
     double bin_min_twist,bin_min_stretch;
     int nBins_twist, nBins_stretch;
     double binW_twist, binW_stretch;
     binW_twist = 0.25;
     binW_stretch = 0.5;
     nBins_twist = 5;
     nBins_stretch = 11;
     bin_min_twist = -0.125;
     bin_min_stretch = -0.1;
        
     i = (int)((twist_update - bin_min_twist)/binW_twist);
        if ( i < 0 ){
            *pt_i = 0;
        }else if( i >= nBins_twist ){
            *pt_i = (nBins_twist - 1);
        }else{
            *pt_i = i;
        }
        
     j = (int)((stretch_update - bin_min_stretch)/binW_stretch);
        if ( j < 0){
             *pt_j = 0;
        }else if( j >= nBins_stretch ){
             *pt_j  = (nBins_stretch - 1);
        }else{
             *pt_j = j;
        }
}



int main(){
    
    int i,j;
    int ii, jj;
    int k, l, m, n;
    int index_i, index_j;
    int index;
    int hhh[10][10];
    
    int **p;
    p = new int*[5];
    
    int **h;
    h = new int*[5];
    
    for(int i=0; i<11; i++){
        
        p[i] = new int[11];
    }
    
    for(int i=0; i<11; i++){
        h[i] = new int[11];
    }
    
    for(i=0; i<5; i++){
        for(j=0; j<11; j++){
            p[i][j] = 0;
        }
    }
    
    for(ii=0; ii<5; ii++){
        for(jj=0; jj<11; jj++){
            p[ii][jj] = 0;
        }
    }
    
    //cout << endl;
    getBinIndex(0.5, 1.5, &k, &l);
    cout << " h_k_l: " << h[k][l] << endl;
    //h[k][l] = 5;
    //cout << " h_k_l: " << h[k][l] << endl;
    
    //getBinIndex(1, 2.5, &m, &n);
    //p[m][n] = h[k][l];
    
    //cout << " p_m_n: " << p[m][n] << endl;
    //p[m][n] *= 0.5;
    //cout << " HERE p_m_n: " << p[m][n] << endl;
    
    
    
    /*
    cout << " count before: " << p[4][10] << endl;
    for(int k=0; k<4; k++){
    getBinIndex(6,11, &i, &j);
    p[i][j] +=1;
    }
    cout << " i: " << i << " j: " << j << endl;
    cout << " row 4: " << i << " column 10: " << j << endl;
    
    
    
    cout << " count after: " << p[4][10] << endl;*/
    
    
}







/*
typedef struct 
{
    int i;
    int j;
    
}coordinate;

coordinate myfunc(int a, int b)
{
    //coordinate retval = {a + b, a * b};
    //return retval;
    int i = a + b;
    int j = a * b;
    coordinate retval = {i,j};
    return retval;
}

int main (void)
{
    coordinate coord = myfunc(3,4);
    printf("X,Y:%d %d\n", coord.i, coord.j);
    //printf("Y: %d\n", coord.y);
}
*/

/*
int x,y; // globals
class enclose { // enclosing class
    int x; // note: private members
    static int s;
 public:
    struct inner { // nested class
        void f(int i) {
            x = i; // Error: can't write to non-static enclose::x without instance
            int a = sizeof x; // Error until C++11,
                              // OK in C++11: operand of sizeof is unevaluated,
                              // this use of the non-static enclose::x is allowed.
            s = i;   // OK: can assign to the static enclose::s
            ::x = i; // OK: can assign to global x
            y = i;   // OK: can assign to global y
        }
        void g(enclose* p, int i) {
            p->x = i; // OK: assign to enclose::x
        }
    };
};
*/

















/*
int main()
{
   int num_elem[3];
   double count[3];
    
   int arr_1[5] = {0,1,2,3,4};
   int arr_2[4] = {0,1,2,3};
   int arr_3[3] = {0,1,2};
   
   
   for(size_t li = 0, ui = 0, ni = 0;
       li < sizeof(arr_1) && ui < sizeof(arr_2) 
       && ni < sizeof(arr_3);
       ++li, ++ui, ++ni){

      count[0] = sizeof(arr_1)/sizeof(int);    
      count[1] = sizeof(arr_2)/sizeof(int); 
      count[2] = sizeof(arr_3)/sizeof(int); 
       
      }
    cout << " num_elem_arr_1 " << count[0] 
         << " num_elem_arr_2 " << count[1] 
         << " num_elem_arr_3 " << count[2] 
         << endl;
    cout << " int: " << int(0.5) << endl;
    cout << " num: " << double(smart+plus) << endl;
}
*/




/*
int main(int argc, char* argv[])
{
   char ls[] = {'a', 'b', 'c'};
   char us[] = {'A', 'B', 'C'};
   int ns[] = {1, 2, 3};
 
   for(size_t li = 0, ui = 0, ni = 0;
       li < sizeof(ls) && ui < sizeof(us) && ni 
< sizeof(ns) / sizeof(int);
       ++li, ++ui, ++ni)
   {
      std::cout << ls[li] << us[ui] << ns[ni] <<
 "\n";
   }
}
*/







 //clean up the particle's old cell if it is now empty.
   /*if((criteria) < (oldE - newE)){
         
       for(int i=0; i<clusterSize; i++){
           for(int k=0; k<4; k++){
             if( oldCell[i*4+k]->firstParticle == NULL ){
                if (cellHash->removeItem( oldCell[i*4+k]->myId )){
                  delete oldCell[i*4+k];
                  oldCell[i*4+k] = NULL;
                }
             }
           }
       }*/
       












