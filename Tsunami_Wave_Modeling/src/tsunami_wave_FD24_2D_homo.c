#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "wave_forward.h"

void tsunami_FD24_homo();

int main(int argc, char *argv[]){
	if(argc != 2){
		fprintf(stderr, "Usage: %s input_para.txt\n", argv[0]);
		exit(1);
	}
    tsunami_FD24_homo(argv[1]);
    return 0;
}
void tsunami_FD24_homo(char *input){
    int nt, nx, ny, sx, sy, i, j, k, index = 0;
    float dt, dx, dy, **d, fc, t0, A, ***h, g = 9.81;
	float fdm[] = {-1/12., 4/3., -5/2., 4/3., -1/12.}; 
	float t1, t2;
    char ss[128], num[5];
	WAVEPARA wp;

// Read in parameters.
	wp = read_wave_para(input);
    nt = wp.nt; dt = wp.dt;
    nx = wp.nx; ny = wp.ny; sx = wp.sx; sy = wp.sy;
    fc = wp.fc; t0 = 1./fc;

// Allocate dynamic memories and source time function.
	printf("Alocate dynamic memories ...\n");
	h = (float ***) malloc(sizeof(float **) * nt);
	d = (float **) malloc(sizeof(float *) * ny);
	for(i = 0; i < ny; i ++)
		d[i] = (float *) malloc(sizeof(float) * nx);
	for(i = 0; i < nt; i ++){
		h[i] = (float **) malloc(sizeof(float *) * ny);
		for(j = 0; j < ny; j ++)
			h[i][j] = (float *) malloc(sizeof(float) * nx);
	}
// Read in sea depth data.
	read_field_2D(wp.field, nx, ny, d);

    // Modeling.
	float c = 250.;
	dx = wp.dx;
	/*
	float dmax;
	dmax = max_or_min(d, nx, ny, 0);
	dx = 7. * sqrt(2.) * sqrt(dmax*g) * dt / 6.; dy = dx;
	*/
    printf("Modeling ...\n");
	float xterm, yterm;
    FILE *fs;
    char tfname[128];
    sprintf(tfname, "FDM24_HOMO_%.0f.tf", fc);
    fs = fopen(tfname, "w");
    fprintf(fs, "0.000000\n");
    t1 = clock();
    for(i = 1; i < nt-1; i ++){
        for(j = 2; j < ny-2; j ++) {
            h[i][j][0] = 0.; h[i][j][nx-1] = 0.;
			h[i][j][1] = 0; h[i][j][nx-2] = 0.;
            for(k = 2; k < nx-2; k ++){
                h[i][0][k] = 0.; h[i][ny-1][k] = 0.;
				h[i][1][k] = 0.; h[i][ny-2][k] = 0.;
				A = pow(dt*c/dx, 2);
                if(j == sy && k == sx)
                    h[i+1][j][k] = ricker(i*dt, fc, 10*dt, 2.);
                else {
                    h[i+1][j][k] = 2*h[i][j][k] - h[i-1][j][k];
                    yterm = A*(h[i][j+2][k]*fdm[0]+h[i][j+1][k]*fdm[1]+h[i][j][k]*fdm[2]+h[i][j-1][k]*fdm[3]+h[i][j-2][k]*fdm[4]);
					xterm = A*(h[i][j][k+2]*fdm[0]+h[i][j][k+1]*fdm[1]+h[i][j][k]*fdm[2]+h[i][j][k-1]*fdm[3]+h[i][j][k-2]*fdm[4]);
					h[i+1][j][k] += (xterm+yterm);
				}
            }
		}
        fprintf(fs, "%f\n", h[i][500][620]);
    }
    t2 = clock();
    fprintf(fs, "0.000000\n");
    fclose(fs);
	printf("Costing time: %.4f second!\n", (t2-t1)/CLOCKS_PER_SEC);

// Write out tsunami wave height data.
    FILE *fp;
	printf("Write out sea water height ...\n");
    if(fp=fopen("height_FD24.txt", "w")){
		for(i = 0; i < ny; i ++)
			for(j = 0; j <nx; j ++){
				if((j+1) == nx)
					fprintf(fp, "%.6f\n", h[wp.wn][i][j]);
				else
					fprintf(fp, "%.6f ", h[wp.wn][i][j]);
			}
	} else
		fprintf(stderr, "Can not open file!!!\n");
    fclose(fp);
// Release dynamic memories.
	for(i = 0; i < nt; i ++){
		free(d[i]);
		for(j = 0; j < ny; j ++)
			free(h[i][j]);
		free(h[i]);
	}
	free(d); free(h);
}
