/* ------------------------------------------------------------------- */
/*      Initially coded by FENG Xuping @ SUSTech on Sept. 12th, 2019.  */
/* ------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "seis_wave.h"

int main(int argc, char *argv[]) {
    if(argc != 10) {
        fprintf(stderr, "Usage: acoustic_wave_2D nt nx nz sx sz dt c f wn\n");
        fprintf(stderr, "       nt: Number of time points;\n");
        fprintf(stderr, "       nx: Number of distance points on X axis;\n");
        fprintf(stderr, "       nz: Number of distance points on Z axis;\n");
        fprintf(stderr, "       sx: Position index of source on X axis;\n");
        fprintf(stderr, "       sz: Pointion index of source on Z axis;\n");
        fprintf(stderr, "       dt: Time interval of wave propagation;\n");
        fprintf(stderr, "       c:  Constant velocity of wave propagation;\n");
        fprintf(stderr, "       f:  Center frequency of source time function;\n");
        fprintf(stderr, "       wn: Number of output wavefield snapshots interval;\n");
        exit(1);
    }

    int nt, nx, nz, sx, sz, wn;
    float c, f, dt, ***u;

    // Get input command line parameters.
    nt = atoi(argv[1]); nx = atoi(argv[2]); nz = atoi(argv[3]);
    sx = atof(argv[4]); sz = atof(argv[5]);
    dt = atof(argv[6]); c = atof(argv[7]); f = atof(argv[8]); wn = atoi(argv[9]);

    // Allocate dynamic memory of wavefield matrix: u.
    u = (float ***) malloc(sizeof(float **) * nt);
    for( int i = 0; i < nt; i ++) {
        u[i] = (float **) malloc(sizeof(float *) * nx);
        for(int j = 0; j < nx; j ++)
            u[i][j] = (float *) malloc(sizeof(float) * nz);
    }

    // Generate 2-D acoustic wavefield.
    aco_2D(nt, nx, nz, sx, sz, dt, c, f, u);

    // Ouput 2-D wavefield data.
    if(wn > 0 && wn < nt){
        system("mkdir -p Wavefield_2D_Snapshots");
        int index = nt / wn;
        for(int i = 0; i <= index; i ++) {
            char fname[128] = {"Wavefield_2D_Snapshots/WAVE_SNAPSHOT_2D_"};
            int wi = i * wn;
            if(wi<10)
                strcat(fname, "000");
            else if(wi>9 && wi < 100)
                strcat(fname, "00");
            else if(wi>99 && wi <1000)
                strcat(fname, "0");
            else
                ;
            char str[10];
            sprintf(str, "%d.txt", wi);
            strcat(fname, str);
            FILE *fp;
            if((fp=fopen(fname, "w"))){
                fprintf(fp, "%d %d %d %d %d %f %f %f %f %f\n", \
                        nt, nx, nz, sx, sz, dt, dt*c*2, wi*dt, f, c);
                for(int j = 0; j < nx; j ++){
                    for(int k = 0; k < nz; k ++)
                    if((k+1) == nz)
                        fprintf(fp, "%f\n", u[wi][j][k]);
                    else
                        fprintf(fp, "%f ", u[wi][j][k]);
                }

            } else
                fprintf(stderr, "Can not open file %s!!!", fname);
            fclose(fp);
        }
    } else
        printf("Number of output snapshot should be between 1 and %d!!!\n", nt);
    
    // Output source time function.
    FILE *fp;
    if(fp=fopen("Wavefield_2D_Snapshots/source_time_function.txt", "w")) {
        fprintf(fp, "%d %f\n", nt, dt);
        for(int i = 0; i < nt; i ++)
            fprintf(fp, "%f ", u[i][sx][sz]);
    } else
        fprintf(stderr, "Can not open file for writing source time function!\n");
    fclose(fp);

    // Release dynamic memory of matrices u and s.
    for(int i = 0; i < nt; i++) {
        for(int j = 0; j < nx; j ++)
            free(u[i][j]);
        free(u[i]);
    }
    free(u);
    printf("        Wavefield and source time function files\n\
            have been saved in directory \"Wavefield_2D_Snapshots\"!\n");

    return 0;
}
