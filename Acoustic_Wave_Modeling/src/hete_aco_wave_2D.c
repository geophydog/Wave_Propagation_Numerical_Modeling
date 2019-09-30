#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "seis_wave.h"

int main(int argc, char **argv){
	if(argc != 4){
		fprintf(stderr, "Usage: hete_aco_wave_2D velocity_file receiver_file wn\n");
		exit(1);
	}
	int wn = 0;
	wn = atoi(argv[3]);
	aco_2D_vel(argv[1], argv[2], wn); 
	return 0;
}
