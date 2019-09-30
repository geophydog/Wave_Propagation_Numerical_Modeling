/* -------------------------------------------------------------- */
/* txt2sac_1d: Convert 1-D wavefield data stored in a text file   */
/*                 generated by "acoustic_wave_1d" to SAC format  */
/*                 files.                                         */
/* Initially coded by FENG Xuping @ SUSTech on Sept. 12th, 2019.  */
/* -------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sacio.h"
#define SIZE 8092

int main(int argc, char **argv){
    if(argc < 2){
        fprintf(stderr, "Usage: txt2sac wave_snapshot.txt\n");
        fprintf(stderr, "       wave_snapshot.txt: Text file generated by \"acoustic_wave_1D\".\n");
        exit(1);
    }
    FILE *fp;
    float dt = 0., dx = 0., **amp;
    char ss[SIZE] = {};
    int count = 0, nt = 0, nx = 0, sx = 0, si = 0, xi = 0;
    fp = fopen(argv[1], "r");
    if(NULL == fp){
        fprintf(stderr, "Can not open file: %s!!!\n", argv[1]);
        exit(1);
    }
    // Read in data.
    while(fgets(ss, SIZE, fp)){
        // Get time interval, distance interval, numbers of time and distance. 
        if(0 == count) {
            sscanf(ss, "%d %f %f %d %d", &sx, &dt, &dx, &nt, &nx);
            // Matrix amp to store wavefield data.
            amp = (float **) malloc(sizeof(float*) * nt);
            for(int i = 0; i < nt; i ++)
                amp[i] = (float *) malloc(sizeof(float) * nx);
            for(int i = 0; i < nt; i ++)
                for(int j = 0; j < nx; j ++)
                    amp[i][j] = 0.;
        }
        // Get wavefield data from the 2nd line to the last line.
        else if(count > 0 && count < nt+1){
            si = 0;
            xi = 0;
            char num[16] = {};
            for(int i = 0; i < strlen(ss); i ++) {
                if( *(ss+i) != 32 && *(ss+i) != 10) {
                    *(num+si) = *(ss+i);
                    si ++;
                } else {
                    si = 0;
                    amp[count-1][xi] = atof(num);
                    xi ++;
                }
            }
            amp[count-1][xi] = atof(num);
            //printf("num: %s   amp: %f\n", num, atof(num));
        } else {
            fprintf(stderr, "There are too many lines in file %s!!!", argv[1]);
            break;
        }

        count ++;
    }
    fclose(fp);

    // Write out wavefield data and convert them to time series.
    printf("Now writing out SAC files ...\n");
    SACHEAD hd = new_sac_head(dt, nt, 0.);
    float *data = malloc(sizeof(float) * nt);
    for(int i = 0; i < nx; i ++) {
        hd.stlo = i*dx; hd.stla = 0.;
        hd.evlo = sx*dx; hd.evla = 0.;
        hd.npts = nt; hd.delta = dt; hd.b = 0.;
        char fname[32]="Acoustic_Wave_1D_", str[5];
        if(i < 10)
            strcat(fname, "000");
        else if(i >= 10 && i < 100)
            strcat(fname, "00");
        else if(i >= 100 && i < 1000) 
            strcat(fname, "0");

        else
            ;
        sprintf(str, "%d", i);
        strcat(fname, str);
        strcat(fname, ".SAC");
        for(int j = 0; j < nt; j ++)
            *(data+j) = *(amp[j]+i);
        //printf("SAC file: %s\n", fname);
        write_sac(fname, hd, data);
    }

    // Release dynamic memory of wavefield data matrix: amp and temporary matrix: data.
    for(int i = 0; i < nt; i ++)
        free(amp[i]);
    free(amp); free(data);

    // Move all converted SAC files into new directory.
    system("mkdir -p Time_Series_1D_SAC"); system("mv Acoustic_*1D*.SAC Time_Series_1D_SAC/");
    printf("SAC format files have been moved into \"Time_Series_1D_SAC\"!\n");
    return 0;

}
