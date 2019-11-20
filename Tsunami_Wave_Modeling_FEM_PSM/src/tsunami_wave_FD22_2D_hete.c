#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "wave_forward.h"

void tsunami_FD22_hete();

int main(int argc, char *argv[]){
	if(argc != 2){
		fprintf(stderr, "Usage: %s input_para.txt\n", argv[0]);
		exit(1);
	}
    tsunami_FD22_hete(argv[1]);
    return 0;
}
void tsunami_FD22_hete(char *input){
    int nt, nx, ny, sx, sy, i, j, k;
    float dt, dx, dy, **d, fc, ***h, g = 9.81;
    char ss[128], num[5];
	WAVEPARA wp;

// Read in parameters.
	wp = read_wave_para(input);
    nt = wp.nt; dt = wp.dt;
    nx = wp.nx; ny = wp.ny; sx = wp.sx; sy = wp.sy;
    fc = wp.fc; dx = wp.dx; dy = wp.dy;

// Allocate dynamic memories and source time function.
	printf("Alocate dynamic memories ...\n");
	h = (float ***) malloc(sizeof(float **) * 3);
	d = (float **) malloc(sizeof(float *) * ny);
	for(i = 0; i < ny; i ++)
		d[i] = (float *) malloc(sizeof(float) * nx);
	for(i = 0; i < 3; i ++){
		h[i] = (float **) malloc(sizeof(float *) * ny);
		for(j = 0; j < ny; j ++)
			h[i][j] = (float *) malloc(sizeof(float) * nx);
	}
// Read in sea depth data.
	read_field_2D(wp.field, nx, ny, d);

// Propagation modeling.
	float term1, term2;
    printf("Modeling ...\n");
    FILE *fs, *fp;
    float t1, t2;
    t1 = clock();
    fs = fopen("FDM22_HETE.tf", "w");
    fprintf(fs, "0.000000\n");
    for(i = 1; i < nt-1; i ++){
        for(j = 1; j < ny-1; j ++){
            h[1][j][0] = 0.; h[1][j][nx-1] = 0.;
            for(k = 1; k < nx-1; k ++){
                h[1][0][k] = 0.; h[1][ny-1][k] = 0.;
                if(j == sy && k == sx)
                    h[2][j][k] = ricker(i*dt, fc, 10*dt, 2.);
                else {
                    h[2][j][k] = 2*h[1][j][k] - h[0][j][k];
					term1 = g*d[j][k]*pow(dt/dx, 2.) * (h[1][j+1][k]+h[1][j-1][k]+h[1][j][k+1]+h[1][j][k-1]-4*h[1][j][k]);
					term2 = 1./4*g*pow(dt/dx, 2.)*( (d[j][k+1]-d[j][k-1])*(h[1][j][k+1]-h[1][j][k-1]) \
												 + (d[j+1][k]-d[j-1][k])*(h[1][j+1][k]-h[1][j-1][k]) );
					h[2][j][k] += (term1+term2);
				}
            }
		}
        fprintf(fs, "%f\n", h[2][300][600]);

        //Write out sea wave height snapshot.
        if(i%wp.wn == 0){
            sprintf(ss, "height_FD22_hete_%.3f.txt", i*dt);
            if(fp=fopen(ss, "w")){
                printf("Writing %s ...\n", ss);
                for(j = 0; j < ny; j ++)
                    for(k = 0; k < nx; k ++) {
                        if((k+1) == nx)
                            fprintf(fp, "%f\n", h[2][j][k]);
                        else
                            fprintf(fp, "%f ", h[2][j][k]);
                    }
            } else
                fprintf(stderr, "Can not open file %s!\n", ss);

        }
        // Exchange sea wave height values for iterations.
        for(j = 0; j < ny; j ++)
            for(k = 0; k < nx; k ++){
                h[0][j][k] = h[1][j][k];
                h[1][j][k] = h[2][j][k];
            }
    }
    fprintf(fs, "0.000000\n");
    fclose(fs); fclose(fp);
    t2 = clock();
	printf("Costing time: %.4f second!\n", (t2-t1)/CLOCKS_PER_SEC);

// Release dynamic memories.
	for(i = 0; i < 3; i ++){
		free(d[i]);
		for(j = 0; j < ny; j ++)
			free(h[i][j]);
		free(h[i]);
	}
	free(d); free(h);
}
