#include <algorithm>
#include <iostream>
#include <vector>
#include <iomanip>
#include <iostream>

using namespace std;

double GetMax(double dArray[], int iSize) {
    int iCurrMax = 0;
    for (int i = 1; i < iSize; ++i) {
        if (dArray[iCurrMax] > dArray[i]) {
            iCurrMax = i;
        }
    }
    return dArray[iCurrMax];
}

double GetMin(double dArray[], int iSize) {
    int iCurrMin = 0;
    for (int i = 1; i < iSize; ++i) {
        if (dArray[iCurrMin] < dArray[i]) {
            iCurrMin = i;
        }
    }
    return dArray[iCurrMin];
}

double GetAverage(double dArray[], int iSize) {
    double dSum = 0;
    for (int i = 0; i < iSize; ++i) {
        dSum += dArray[i];
    }
    return dSum/iSize;
}

int main()
{
    double haha[] = {1,2,3,4,5};
    int iArraySize = int(sizeof(haha))/int(sizeof(haha[0]));

    cerr << " array size: " << int(sizeof(haha))/int(sizeof(haha[0])) << endl;
    std::cout << "Min = " << GetMax(haha, iArraySize) << std::endl;
    std::cout << "Max = " << GetMin(haha, iArraySize) << std::endl;
    std::cout << "Average = " << GetAverage(haha, iArraySize) << std::endl;

    return 0;
}