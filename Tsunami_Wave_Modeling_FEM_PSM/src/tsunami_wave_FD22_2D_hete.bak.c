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
	float t1, t2;
	t1 = clock();
    tsunami_FD22_hete(argv[1]);
	t2 = clock();
	printf("Costing time: %.4f second!\n", (t2-t1)/CLOCKS_PER_SEC);
    return 0;
}
void tsunami_FD22_hete(char *input){
    int nt, nx, ny, sx, sy, i, j, k, index = 0;
    float dt, dx, dy, **d, fc, t0, A, ***h, *s, pi = 3.1415926535, g = 9.81;
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
	s = (float *) malloc(sizeof(float) * nt);
	for(i = 0; i < nt; i ++)
        s[i] = ricker(i*dt, fc, 20*dt, 2.);
	for(i = 0; i < nt; i ++){
		h[i] = (float **) malloc(sizeof(float *) * ny);
		for(j = 0; j < ny; j ++)
			h[i][j] = (float *) malloc(sizeof(float) * nx);
	}
// Read in sea depth data.
	read_field_2D(wp.field, nx, ny, d);

// Propagation modeling.
	float term1, term2;
	dx = wp.dx;
    printf("Modeling ...\n");
    FILE *fs;
    fs = fopen("FDM22_HETE.tf", "w");
    fprintf(fs, "0.000000\n");
    for(i = 1; i < nt-1; i ++){
        for(j = 1; j < ny-1; j ++){
            h[i][j][0] = 0.; h[i][j][nx-1] = 0.;
            for(k = 1; k < nx-1; k ++){
                h[i][0][k] = 0.; h[i][ny-1][k] = 0.;
                if(j == sy && k == sx)
                    h[i+1][j][k] = s[i];
                else {
                    h[i+1][j][k] = 2*h[i][j][k] - h[i-1][j][k];
					term1 = g*d[j][k]*pow(dt/dx, 2.) * (h[i][j+1][k]+h[i][j-1][k]+h[i][j][k+1]+h[i][j][k-1]-4*h[i][j][k]);
					term2 = 1./4*g*pow(dt/dx, 2.)*( (d[j][k+1]-d[j][k-1])*(h[i][j][k+1]-h[i][j][k-1]) \
												 + (d[j+1][k]-d[j-1][k])*(h[i][j+1][k]-h[i][j-1][k]) );
					h[i+1][j][k] += (term1+term2);
				}
            }
		}
        fprintf(fs, "%f\n", h[i][300][600]);
    }
    fprintf(fs, "0.000000\n");
    fclose(fs);

// Write out tsunami wave height data.
    FILE *fp;
	printf("Write out sea water height ...\n");
    if(fp=fopen("height_FD22_hete.txt", "w")){
		for(i = 0; i < ny; i ++)
			for(j = 0; j <nx; j ++){
				if((j+1) == nx)
					fprintf(fp, "%.6f\n", h[400][i][j]);
				else
					fprintf(fp, "%.6f ", h[400][i][j]);
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
	free(d); free(h); free(s);
}
