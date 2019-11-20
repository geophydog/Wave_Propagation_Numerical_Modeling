/* --------------------------------------------------------------- */
/*   acoustic_wave_1D: Generate 1-D acoustic wavefield with        */
/*                        finite difference method.                */
/*   Initially coded by FENG Xuping @ SUSTech on Sept. 12th, 2019. */ 
/* --------------------------------------------------------------- */

#include <stdlib.h>
#include <stdio.h>
#include "seis_wave.h"

int main(int argc, char *argv[]){
    if(argc != 7) {
        fprintf(stderr, " Usage: acoustic_wave_1D nt nx sx dt f c\n");
        fprintf(stderr, "        nt: Number of time points;\n");
        fprintf(stderr, "        nx: Number of X axis point;\n");
        fprintf(stderr, "        sx: Position index of source on X axis;\n");
        fprintf(stderr, "        dt: Time interval;\n");
        fprintf(stderr, "        f:  Center frequency of source time function;\n");
        fprintf(stderr, "        c:  Constant velocity of acoustic wave propagation.\n");
        exit(1);
    }
    int nt, nx, sx;
    float dt, f, c, **u;
    // Get input parameters.
    nt = atoi(argv[1]); nx = atoi(argv[2]); sx = atoi(argv[3]);
    dt = atof(argv[4]); f = atof(argv[5]); c = atof(argv[6]);
    
    // Allocate dynamic memory of displacement matrix: u.
    u = (float **) malloc(sizeof(float *) * nt);
    for(int i = 0; i < nt; i ++)
        u[i] = (float *) malloc(sizeof(float) * nx);

    // Compute displacement.
    printf("Find solutions for wavefield ...\n");
    aco_1D(nt, nx, sx, dt, f, c, u);

    // Ouput displacement in a text file.
    printf("Write out wavefield data ...\n");
    FILE *fp;
    if((fp=fopen("WAVE_SNAPSHOT_1D.txt", "w"))){
        fprintf(fp, "%d %.6f %.6f %d %d\n", sx, dt, c*dt*1.01, nt, nx);
        for(int i = 0; i < nt; i ++)
            for(int j = 0; j < nx; j ++)
                if((j+1) == nx)
                    fprintf(fp, "%.6f\n", u[i][j]);
                else 
                    fprintf(fp, "%.6f ", u[i][j]);
        fclose(fp);
        printf("1-D wavefield data have been saved in file WAVE_SNAPSHOT_1D.txt!\n"); 
    } else
        fprintf(stderr, "Can not open file or write out wavefield data!!!\n");
    // Release dynamic memory of displacement maxtrix: u.
    for(int i = 0; i < nt; i ++)
        free(u[i]);
    free(u);

    return 0;
}
