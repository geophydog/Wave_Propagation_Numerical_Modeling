/* -------------------------------------------------------------------------------- */
/*           Initially coded by FENG Xuping @ SUSTech on Sept. 12th, 2019.          */
/* -------------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "wave_forward.h"

#define PI 3.1415926535
#define SIZE 32768

/* ------------------------ Default value of WAVEPARA struct --------------------------- */
static WAVEPARA wavepara_null = {
	-12345, -12345., -12345,
	-12345, -12345, -12345., -12345.,
	-12345.,
	-12345, -12345,
	{"-12345"},
	{"-12345"}
};

/* --------------------- Generate Ricker wavelet with shifted time --------------------- */
float ricker(float t,float f0, float t0, float amp){
    return amp*(1.-2.*pow(PI*f0, 2.)*pow((t-t0), 2.))*exp(-pow(PI*f0, 2.)*pow((t-t0), 2.));
}

/* ------------------- Find maximum or minimum of input 2D field ----------------------- */
float max_or_min(float **field, int nx, int ny, int m){
	int i, j;
	float v0 = field[0][0];
	for(i = 0; i < ny; i ++)
		for(j = 0; j < nx; j ++){
			if(0 == m){
				if(v0 < field[i][j])
					v0 = field[i][j];
			}
			else{
				if(v0 > field[i][j])
					v0 = field[i][j];
			}
		}
	return v0;
}

/* ------------------------- Read in input parameters --------------------------- */
WAVEPARA read_wave_para(char *input){
	FILE *fin;
	char ss[200];
	int index = 1;
	WAVEPARA wp = wavepara_null;
	if(fin=fopen(input, "r"))
		while(fgets(ss, 200, fin)){
			if(1 == index)
				sscanf(ss, "%*s %d %*s", &wp.nt);
			else if(2 == index)
				sscanf(ss, "%*s %f %*s", &wp.dt);
			else if(3 == index)
				sscanf(ss, "%*s %d %*s", &wp.wn);
			else if(4 == index)
				sscanf(ss, "%*s %d %*s", &wp.nx);
			else if(5 == index)
				sscanf(ss, "%*s %d %*s", &wp.ny);
			else if(6 == index)
				sscanf(ss, "%*s %f %*s", &wp.dx);
			else if(7 == index)
				sscanf(ss, "%*s %f %*s", &wp.dy);
			else if(8 == index)
				sscanf(ss, "%*s %f %*s", &wp.fc);
			else if(9 == index)
				sscanf(ss, "%*s %d %*s", &wp.sx);
			else if(10 == index)
				sscanf(ss, "%*s %d %*s", &wp.sy);
			else if(11 == index)
				sscanf(ss, "%*s %s %*s", wp.field);
			else if(12 == index)
				sscanf(ss, "%*s %s %*s", wp.receiver);
			else
				break;
			index ++;
		}
	fclose(fin);
	return wp;
}

/* ----------------- Read in 2D data (nx by ny, e.g. velocity) ------------------ */
void read_field_2D(char *fname, int nx, int ny, float **field){
	int index, i;
	char ss[SIZE];
	FILE *fin;

	printf("Read in data (%d by %d) from file \"%s\" ...\n", nx, ny, fname);
	if(fin=fopen(fname, "r")){
		index = 0;
		while(fgets(ss, SIZE, fin)){
			if(index < ny){
				char *num = strtok(ss, " ");
				i = 0;
				while(num != NULL){
					field[index][i] = atof(num);
					num = strtok(NULL, " ");
					i ++;
				}
			} else{
				fprintf(stderr, "Too many lines in file %s!!!\n", fname);
				break;
			}
			index ++;
		}
	} else{
		fprintf(stderr, "Can not open file %s!!!\n", fname);
	}
	fclose(fin);
}

/* ------- Generate 1-D acoustic wavefield with finite difference method -------- */
void aco_FD_1D(int nt, int nx, int sx, float dt, float f, float c, float **u) {
    // Calculate distance interval along X axis.
    float dx = c * dt * 1.01;
    float A = pow(dt*c/dx, 2);

    // Initialization of displacement matrix: u.
	int i, j;
    for(i = 0; i < nt; i ++)
        for(j = 0; j < nx; j ++)
            u[i][j] = 0.;

    // Generate source time function: s.
    float *s;
    s = (float *) malloc(sizeof(float) * nt );
    for(i = 0; i < nt; i ++)
        s[i] = sin(2.*PI*f*i*dt) * exp(-2.*PI*f*pow(i*dt, 2)*10);

    // Acoustic wave propagation modeling.
    for(i = 1; i < (nt-1); i ++)
        for(j = 1; j < (nx-1); j ++) {
            // Boundary condition: fixed boundary..
            u[i][0] = 0.; u[i][nx-1] = 0.;
            // Set source on sx position on X axis.
            if(j == sx)
                u[i+1][j] = s[i];
            // Solutions for wave equation.
            else
                u[i+1][j] = -u[i-1][j] + 2. * u[i][j] + \
                    A * (u[i][j+1] + u[i][j-1] - 2.*u[i][j]);
        }
}

/* ------- Generate 2-D acoustic wavefield with finite difference method ------ */
void aco_FD_2D(int nt, int nx, int nz, int sx, int sz, float dt, float c, float f, float ***u){

    float A, dx, *s, gamma = 3.;
    dx = dt * c * 8. / 3;
    A = pow(dt*c/dx, 2);
    s = (float *) malloc(sizeof(float) * nt);

    // Initialization of wavefield matrix: u and setting of source time function s.
	int i, j, k;
    for(i = 0; i < nt; i ++) {
        //s[i] = exp(-pow(2*PI*f*i*dt/gamma, 2)) * cos(2*PI*f*i*dt);
        s[i] = sin(2*PI*f*i*dt) * exp(-(2*PI*f*150*pow(i*dt, 2)));
        for(j = 0; j < nx; j ++)
            for(k = 0; k < nz; k ++)
                u[i][j][k] = 0.;
    }

    // Solutions for 2-D wavefield with finite difference method.
    for(i = 1; i < nt-1; i ++)
        for(j = 1; j < nx-1; j ++)
            for(k = 1; k < nz-1; k ++){
                // Boundary condition: fixed boundary.
                u[i][0][k] = 0.; u[i][nx-1][k] = 0.;
                u[i][j][0] = 0.; u[i][j][nz-1] = 0.;
                // Setting position of source.
                if( j == sx && k == sz)
                    u[i+1][j][k] = s[i];
                else
                    u[i+1][j][k] = A*(u[i][j+1][k]+u[i][j-1][k]+u[i][j][k+1]+u[i][j][k-1]-4*u[i][j][k])\
                                   + 2.*u[i][j][k] - u[i-1][j][k];
            }

    // Release dynamic memory of source time function s.
    free(s);
}

/* ------------ Generate 2D acoustic wavefield for given velocity model ------------------- */
void aco_FD_2D_vel(char *model, char *recei, int wn){
        int nt, nx, nz, sx, sz, index = 0, xi = 0, si = 0, sr = 0, b = 0;
		int i, j, k;
        float **vel, ***u, *s, f, dt, dx, v0= 0., A = 0.;
        FILE *fp;
        char ss[SIZE], num[32];

        printf("Reading in model file ...\n");
        if(fp=fopen(model, "r")) {
            while(fgets(ss, SIZE, fp)) {
                if(0 == index) {
                    // Scan input parameters.
					sscanf(ss, "%d %d %d %d %d %d %d %f %f", \
						&nt, &nx, &nz, &sx, &sz, &sr, &b, &f, &dt);
					//printf("nt: %d nz: %d nx: %d sx: %d, sz: %d sr: %d b: %d f: %f dt: %f\n", \
							nt, nz, nx, sx, sz, sr, b, f, dt);
					// Variables for Sine type source time function.
					int fi = (int)(1/f/dt), fh = (int)(1/f/dt/2);
                    // Allocate dynamic memory of velocity, source time function
                    // and displacement matrecs: vel, s and u.
                    vel = (float **) malloc(sizeof(float *) * nz);
                    u = (float ***) malloc(sizeof(float **) * nt);
                    s = (float *) malloc(sizeof(float) * nt);
                    for(i = 0; i < nz; i ++) {
                        vel[i] = (float *) malloc(sizeof(float) * nx);
                        u[i] = malloc(sizeof(float) * nx);
                    }
                    for(i = 0; i < nt; i ++) {
                    u[i] = (float **) malloc(sizeof(float *) * nz);
					// Initialization of source time function.
					s[i] = 0.;
                    for(j = 0; j < nz; j ++)
                        u[i][j] = (float *) malloc(sizeof(float) * nx);
                                }
                    // Generate source time function and initialization of displacement matrix: u.
                    FILE *fs;
                    if(fs=fopen("SOURCE_TIME_FUNCTION.txt", "w")){
                        fprintf(fs, "%d %f\n", nt, dt);
                        for(i = 0; i < nt; i ++) {					
						    // Ricker wavelet for source time function.
                            if(0 == sr)
							    s[i] = (1-2*pow(PI*f*(i*dt-1/f), 2)) * exp(-pow(PI*f*(i*dt-1/f), 2));
							// Like gaussian shape.
							else if(1 == sr)
								s[i] = exp(-pow(i*dt-1/f, 2)*PI*pow(f, 2.5));
							// Sine function squareed.
							else
                                if(i >= fh && i < fi+fh)
								    s[i] = pow(sin(PI*f*(i-fh)*dt), 2);
                                    fprintf(fs, "%f ", s[i]);
                                    for(j = 0; j < nz; j ++)
                                        for(k = 0; k < nx; k ++)
                                            u[i][j][k] = 0.;
                        }
                    } else
                        fprintf(stderr, "Can not open file to save source time function!!!\n");
                    fclose(fs);

                } else if(index > 0 && index < (nz+1)){
                    xi = 0; si = 0;
                    for(i = 0; i < strlen(ss); i ++) {
                        if(*(ss+i) != 32 && *(ss+i) != 10) {
                            *(num+si) = *(ss+i);
                            si ++;
                        } else {
                            si = 0;
                            vel[index-1][xi] = atof(num);
                            if(v0 < vel[index-1][xi])
                            v0 = vel[index-1][xi];
                            xi ++;
                        }
                    }
                    vel[index-1][xi] = atof(num);
                    if(v0 < vel[index-1][xi])
                        v0 = vel[index-1][xi];
                } else {
                    fprintf(stderr, "Too many lines in file %s!!!\n", model);
                    //break;
                }

                index ++;
                }
        } else {
                fprintf(stderr, "Can not open file %s!!!\n", model);
                fclose(fp);
                exit(1);
        }
        fclose(fp);

        // Acoustic wave propagation in given input velocity model.
        printf("Acoustic wave propagation Modeling ...\n");
        dx = dt * v0 * 7. * sqrt(2.) / 6.;
        for(i = 2; i < nt-2; i ++)
            for(j = 1; j < nz-1; j ++)
                for(k = 1; k < nx-1; k ++){
					A = dt * vel[j][k] / dx;
                    // Boundary condition.
					// Fixed boundary if b == 0..
					if(0 == b) {
                        u[i][0][k] = 0.; u[i][1][k] = 0.; u[i][nz-2][k] = 0.; u[i][nz-1][k] = 0.;
                        u[i][j][0] = 0.; u[i][j][1] = 0.; u[i][j][nx-2] = 0.; u[i][j][nx-1] = 0.;
					} 
					// Absorbing boundary condition.
					// 1. Clayton Engquist and Majda for 1st order.
					else if(1 == b){
					    // z = 0
					    u[i+1][0][k] = A*u[i][1][k] + (1-A)*u[i][0][k];
					    // z = nz-1
					    u[i+1][nz-1][k] = A*u[i][nz-2][k] + (1-A)*u[i][nz-1][k];
					    // x = 0
					    u[i+1][j][0] = A*u[i][j][1] + (1-A)*u[i][j][0];
					    // x= nz-1
					    u[i+1][j][nx-1] = A*u[i][j][nx-2] + (1-A)*u[i][j][nx-1];
					}
					// 2. Clayton Engquist and Majda for 2nd order and free surface on top of model..
					else if(2 == b){
						// Free surface boundary.
						u[i][0][k] = 0.; u[i][1][k] = 0.;
						// z = nz-1 on Z axis.
						u[i+1][nz-1][k] = (2-2*A-pow(A, 2))*u[i][nz-1][k] + 2*A*(1+A)*u[i][nz-2][k]\
								-pow(A, 2)*u[i][nz-3][k] + (2*A-1)*u[i-1][nz-1][k] - 2*A*u[i-1][nz-2][k];
						// x = 0 on X axis.
						u[i+1][j][0] = (2-2*A-pow(A, 2))*u[i][j][0] + 2*A*(1+A)*u[i][j][1]\
								-pow(A, 2)*u[i][j][2] + (2*A-1)*u[i-1][j][0] - 2*A*u[i-1][j][1];
						// x = nx-1 on X axis.
						u[i+1][j][nx-1] = (2-2*A-pow(A, 2))*u[i][j][nx-1] + 2*A*(1+A)*u[i][j][nx-2]\
								-pow(A, 2)*u[i][j][nx-3] + (2*A-1)*u[i-1][j][nx-1] - 2*A*u[i-1][j][nx-2];
					}
					// 2. Clayton Engquist and Majda for 2nd order.
					else{
						// z = 0 on Z axis.
						u[i+1][0][k] = (2-2*A-pow(A, 2))*u[i][0][k] + 2*A*(1+A)*u[i][1][k]\
								-pow(A, 2)*u[i][2][k] + (2*A-1)*u[i-1][0][k] - 2*A*u[i-1][1][k];
						// z = nz-1 on Z axis.
						u[i+1][nz-1][k] = (2-2*A-pow(A, 2))*u[i][nz-1][k] + 2*A*(1+A)*u[i][nz-2][k]\
								-pow(A, 2)*u[i][nz-3][k] + (2*A-1)*u[i-1][nz-1][k] - 2*A*u[i-1][nz-2][k];
						// x = 0 on X axis.
						u[i+1][j][0] = (2-2*A-pow(A, 2))*u[i][j][0] + 2*A*(1+A)*u[i][j][1]\
								-pow(A, 2)*u[i][j][2] + (2*A-1)*u[i-1][j][0] - 2*A*u[i-1][j][1];
							// x = nx-1 on X axis.
						u[i+1][j][nx-1] = (2-2*A-pow(A, 2))*u[i][j][nx-1] + 2*A*(1+A)*u[i][j][nx-2]\
								-pow(A, 2)*u[i][j][nx-3] + (2*A-1)*u[i-1][j][nx-1] - 2*A*u[i-1][j][nx-2];
					}
                    // Setting position of source.
                    if(j == sz && k == sx)
                        u[i+1][j][k] = s[i];
                    else
                        u[i+1][j][k] = pow(A, 2) * ( u[i][j+1][k] + u[i][j-1][k] \
                            + u[i][j][k+1] + u[i][j][k-1] - 4*u[i][j][k])\
                            + 2.*u[i][j][k] - u[i-1][j][k];
        }

        // Write out time series recorded by receivers.
        printf("Write out acoustic pressure recorded by receivers ...\n");
        FILE *fr, *fseis;
        int rx, rz;
        char sss[32];
        if((fr=fopen(recei, "r")) && (fseis=fopen("acoustic_pressure.txt", "w"))){
            while(fgets(ss, 32, fp)){
                sscanf(ss, "%d %d", &rx, &rz);
                if(rx > (nx-1) || rz > (nz-1))
                    continue;
                fprintf(fseis, "rx: %d rz: %d dx: %f\n", rx, rz, dx);
                for(i = 0; i < nt; i ++)
                    fprintf(fseis, "%f ", u[i][rz][rx]);
                fprintf(fseis, "\n");
            }
        }
        fclose(fr); fclose(fseis);

        // Write out acoustic wavefield snapshot files.
        printf("Writing out wavefield snapshot files ...\n");
        char str[16];
        int sn;
        sn = (int) nt/wn;
    if(sn*wn >= nt)
        sn --;
    system("rm -rf ACOU_2D_WAVE_SNAPSHOTS");
        system("mkdir -p ACOU_2D_WAVE_SNAPSHOTS");
        for(i = 0; i <= sn; i ++) {
            sprintf(str, "%.6f", i*wn*dt);
            FILE *fp;
            char fname[128] = {"ACOU_2D_WAVE_SNAPSHOTS/ACOUTIC_WAVE_2D_VEL_"};
            strcat(fname, str); strcat(fname, ".txt");
            if(fp=fopen(fname, "w")){
                fprintf(fp, "%d %d %d %d %d %f %f %f %d\n", nt, nx, nz, sx, sz, f, dt, dx, wn);
                for(j = 0; j < nz; j ++)
                    for(k = 0; k < nx; k++) {
                        if((k+1) == nx)
                            fprintf(fp, "%f\n", u[i*wn][j][k]);
                        else
                            fprintf(fp, "%f ", u[i*wn][j][k]);
                    }
                    printf("Write out wavefield snapshot file: %s  for time %f s ...\n", fname, i*wn*dt);
                } else {
                    fprintf(stderr, "Can not open file %s!!!\n", fname);
                    continue;
                }
            fclose(fp);
        }
        // Release dynamic memoried of velocity and displacement matrecs: vel and u.
        for(i = 0; i < nt; i ++) {
            for(j = 0; j < nz; j ++) {
                if(i == 0)
                    free(vel[j]);
                free(u[i][j]);
            }
            free(u[i]);
        }
        free(vel); free(u); free(s);
}
