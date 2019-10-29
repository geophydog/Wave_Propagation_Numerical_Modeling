#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <math.h> 
#include <complex.h>
#include <fftw3.h>
#include <time.h>  
#include "wave_forward.h"

#define pi 3.14159265

void tsunami_PSM_HOMO(char *input);

int main(int argc, char *argv[]){
    if(argc != 2){
        fprintf(stderr, "Usage: %s input_file\n", argv[0]);
        exit(1);
    }
    tsunami_PSM_HOMO(argv[1]);
}

void tsunami_PSM_HOMO(char *input){
    float dx, dz, dt;
    float f0, c;  
    float **u1, **u2, **u3;
    float t1, t2;
    int nt, X, Z, x0,z0;
    int i,j,k;
    char snapshot[128];  
    float kz,kx;  
    FILE *fp;
    fftw_complex *in, *out;
    fftw_plan p;
    
    WAVEPARA wp = read_wave_para(input);
    nt = wp.nt; dt = wp.dt; c = 250.;
    X = wp.nx; Z = wp.ny;
    dx = wp.dx;
    dz = dx;
    f0 = wp.fc;  
    x0 = wp.sx;  
    z0 = wp.sy;      

    // Allocate dynamic memories.
    in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * Z * X);  
    out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * Z * X);
    u1 = (float **) malloc(sizeof(float *) * Z);
    u2 = (float **) malloc(sizeof(float *) * Z);
    u3 = (float **) malloc(sizeof(float *) * Z);
    for(i = 0; i < Z; i ++){
        u1[i] = (float *) malloc(sizeof(float) * X);
        u2[i] = (float *) malloc(sizeof(float) * X);
        u3[i] = (float *) malloc(sizeof(float) * X);
    }

    printf("Tsunami wave propagation modeling ...\n");
    t1 = clock();
    FILE *fs;
    char tfname[128];
    sprintf(tfname, "PSM_HOMO_%.0f.tf", f0);
    fs = fopen(tfname, "w");
    for(k = 0; k < nt; k ++){
        // Setting source with Ricker wavelet.
        u2[z0][x0] = ricker(k*dt, f0, 10*dt, 2.);  
        for(i = 0; i < Z; i ++)  
            for(j = 0; j < X;j ++)  
                in[i*X+j] = u2[i][j] + 0*I;  
  
        p = fftw_plan_dft_2d(Z, X, in, out, FFTW_FORWARD, FFTW_ESTIMATE);  
        fftw_execute(p);  
        fftw_destroy_plan(p);  
        fftw_cleanup();  
        // Generate 2nd order derivative operator in spatial domain with 2D FFT.
        for(i = 0; i < Z; i ++)  
        {  
            if (i < Z/2)  
                kz = pow(2*pi*i/dz/Z, 2.);  
            else  
                kz = pow(2*pi*(i-Z)/dz/Z, 2.);  
            for(j = 0; j < X; j ++)  
            {  
                if(j < X/2)  
                    kx = pow(2*pi*j/dx/X, 2.);  
                else  
                    kx = pow(2*pi*(j-X)/dx/X, 2.);  
                in[i*X+j] = -(kx+kz) * out[i*X+j];  
            }  
        }  
        p = fftw_plan_dft_2d(Z, X, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);  
        fftw_execute(p);  
        fftw_destroy_plan(p);  
        fftw_cleanup();
  
        // Shallow water wave propagation modeling
        for(i = 0; i < Z; i ++)  
            for(j = 0; j < X; j ++){  
                u3[i][j] = pow(c*dt, 2.) * creal(out[i*X+j]) / X / Z + 2 * u2[i][j] - u1[i][j];  
                u1[i][j] = u2[i][j];  
                u2[i][j] = u3[i][j];  
            }  
  
        // Output snapshots of wave propagation with "wn" interval.
        if((k%wp.wn==0) && k)  
        {  
            memset(snapshot, 0, 128);  
            sprintf(snapshot,"Snapshot_%.4f.txt", k*dt);  
            fp = fopen(snapshot, "w");  
            for(i = 0;i < Z; i ++)  
                for(j = 0; j < X; j ++) {
                    if((j+1) == X)
                        fprintf(fp, "%f\n", u3[i][j]);
                    else
                        fprintf(fp, "%f ", u3[i][j]);
                }
            fclose(fp);  
        }  
        fprintf(fs, "%f\n", u3[500][620]);
    }
    fclose(fs);
    t2 = clock();
    printf("Time used: %.4f seconds!\n", (t2-t1)/CLOCKS_PER_SEC);

    // Release dynamic memories.
    fftw_free(in); fftw_free(out);
    for(i = 0; i < Z; i ++){
        free(u1[i]); free(u2[i]); free(u3[i]);
    }
    free(u1); free(u2); free(u3);
}  
