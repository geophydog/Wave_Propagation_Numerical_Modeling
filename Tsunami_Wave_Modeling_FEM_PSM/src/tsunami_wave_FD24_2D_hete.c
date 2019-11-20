#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "wave_forward.h"

#define G 9.81

void tsunami_FD24_hete();

int main(int argc, char *argv[]){
	if(argc != 2){
		fprintf(stderr, "Usage: %s input_para.txt\n", argv[0]);
		exit(1);
	}
    tsunami_FD24_hete(argv[1]);
    return 0;
}
void tsunami_FD24_hete(char *input){
    int nt, nx, ny, sx, sy, i, j, k, index = 0;
    float dt, dx, fc, **d, ***h;
    char ss[128], num[5];
	WAVEPARA wp;

// Read in parameters.
	wp = read_wave_para(input);
    nt = wp.nt; dt = wp.dt; dx = wp.dx;
    nx = wp.nx; ny = wp.ny; sx = wp.sx; sy = wp.sy;
    fc = wp.fc;

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

    float t1, t2;
    t1 = clock();
// Propagation modeling.
	float term1, term2, term3, term4;
    FILE *fp;
    printf("Modeling ...\n");
    for(i = 1; i < nt-1; i ++) {
        for(j = 2; j < ny-2; j ++) {
            h[1][j][0] = 0.; h[1][j][nx-1] = 0.;
			h[1][j][1] = 0.; h[1][j][nx-2] = 0.;
            for(k = 2; k < nx-2; k ++){
                h[1][0][k] = 0.; h[1][ny-1][k] = 0.;
				h[1][1][k] = 0.; h[1][ny-2][k] = 0.;
                if(j == sy && k == sx)
                    h[2][j][k] = ricker(i*dt, fc, 1./fc, 5.);
                else {
                    h[2][j][k] = -h[0][j][k] + (2-5*G*d[j][k]*pow(dt/dx,2))*h[1][j][k];
					term1 = 4/3.*G*d[j][k]*pow(dt/dx,2)*(h[1][j-1][k]+h[1][j+1][k]+h[1][j][k-1]+h[1][j][k+1]);
					term2 = -1/12.*G*d[j][k]*pow(dt/dx,2)*(h[1][j-2][k]+h[1][j+2][k]+h[1][j][k-2]+h[1][j][k+2]);
					term3 = 1/3.*G*pow(dt/dx,2)*( (d[j+1][k]-d[j-1][k])*(h[1][j+1][k]-h[1][j-1][k]) \
												 +(d[j][k+1]-d[j][k-1])*(h[1][j][k+1]-h[1][j][k-1]) );
					term4 = -1/6.*G*pow(dt/dx,2)*( (d[j+1][k]-d[j-1][k])*(h[1][j+2][k]-h[1][j-2][k]) \
												  +(d[j][k+1]-d[j][k-1])*(h[1][j][k+2]-h[1][j][k-2]) );
					h[2][j][k] += (term1+term2+term3+term4);
				}
            }
		}
        // Write out snapshot files.
        if( i%wp.wn == 0){
            if(i < 10)
                sprintf(ss, "height_FD24_hete_000%d_%.3f.txt", i, i*dt);
            else if(i >= 10 && i < 100)
                sprintf(ss, "height_FD24_hete_00%d_%.3f.txt", i, i*dt);
            else if(i >= 100 && i < 1000)
                sprintf(ss, "height_FD24_hete_0%d_%.3f.txt", i, i*dt);
            else
                sprintf(ss, "height_FD24_hete_%d_%.3f.txt", i, i*dt);

            if(fp=fopen(ss, "w")){
                printf("writing %s ...\n", ss);
                for(j = 0; j < ny; j ++)
                    for(k = 0; k < nx; k ++){
                        if((k+1)==nx)
                            fprintf(fp, "%f\n", h[2][j][k]);
                        else
                            fprintf(fp, "%f ", h[2][j][k]);
                    }
            } else
                fprintf(stderr, "Can not open file %s!!!\n", ss);
            fclose(fp);
        }
        // Exchange values for iterations.
        for(j = 2; j < ny-2; j ++)
            for(k = 2; k < nx-2; k ++){
                h[0][j][k] = h[1][j][k];
                h[1][j][k] = h[2][j][k];
            }
    }
    t2 = clock();
    printf("Time used: %f seconds!\n", (t2-t1)/CLOCKS_PER_SEC);
// Release dynamic memories.
	for(i = 0; i < 3; i ++){
		free(d[i]);
		for(j = 0; j < ny; j ++)
			free(h[i][j]);
		free(h[i]);
	}
	free(d); free(h);
}
