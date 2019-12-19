#include <iostream>
#include <cmath>
#include "advection_forward.h"

using namespace std;
int main(int argc, char *argv[]){
    char flag[] = {"Gauss"};
    if(argc != 3){
        cout << "Usage: " << argv[0] << "input output" << endl;
        exit(1);
    }
    Advection2D ad = Advection2D(argv[1]);
    ad.Forward();
    ad.WriteSnapshot(argv[2]);
    return 0;
}
