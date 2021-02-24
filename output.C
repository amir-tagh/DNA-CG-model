#include <iostream>
#include <fstream>

using namespace std;

int main(){
    
    int rowCount = 2;
    int colCount = 3;
    int i, j;
    
    int **ary = new int*[rowCount];
    for(i=0; i<rowCount; i++)
        ary[i] = new int[colCount];
    
    for(i=0; i<rowCount; i++){
        for(j=0; j< colCount; j++){
            ary[i][j] = 1;
        }
    }
    
    ofstream out("out.txt", ios::out);
    for(i=0; i<rowCount; i++){
        for(j=0; j<colCount; j++){
            out << ary[i][j] << ',';
            out << '\n';
        }
    }
    out.close();
    
    cout << " element: " << ary[0][0] << endl;
    
}