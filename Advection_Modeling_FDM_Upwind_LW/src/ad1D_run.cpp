#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "advection_forward.h"

int main(int argc, char *argv[]){
    if(argc != 3){
        fprintf(stderr, "Usage: %s input output\n", argv[0]);
        exit(1);
    }
    Advection1D ad = Advection1D(argv[1]);
    ad.LW();
    ad.WriteSnapshot(20, argv[2]);
    return 0;
}
