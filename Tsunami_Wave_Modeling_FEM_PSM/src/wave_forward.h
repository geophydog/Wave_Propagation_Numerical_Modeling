/* ---------------------------------------------------------------- */
/*   Initially coded by FENG Xuping @ SUSTech on Sept. 12th, 2019.  */
/* ---------------------------------------------------------------- */

#ifndef _WAVE_FORWARD_H
#define _WAVE_FORWARD_H

/* --------------------------- STRUCT of input parameters ------------------------------ */
typedef struct wave_para{
	int nt;            // Number of time point;
	float dt;          // Time interval;
	int wn;            // Interval number of output snapshots;
	int nx;            // Number of sampling on X axis;
	int ny;            // Number of sampling on Y axis;
	float dx;          // X interval;
	float dy;          // Y interval;
	float fc;          // Center frequency of source;
	int sx;            // Position index of source on X axis;
	int sy;            // Position index of source on Y axis;
	char field[128];   // 2D field data (ny by nx) file name;
	char receiver[128];// Receiver (several "rx, ry" pairs) file name.
} WAVEPARA;

/* ++++++++++++++++ Function prototypes of solving acoustic wave equation ++++++++++++++ */
/* --------------------- Generate Ricker wavelet with shifted time --------------------- */
float ricker(float t,float f0, float t0, float amp);

/* ----------- Generate 1-D acoustic wavefield with finite difference method ----------- */
void aco_FD_1D(int nt, int nx, int sx, float dt, float f, float c, float **u);

/* ----------- Generate 1-D acoustic wavefield with finite difference method ----------- */
void aco_FD_2D(int nt, int nx, int nz, int sx, int sz, float dt, float c, float f, float ***u);

/* ----- Generate 2-D acoustic wavefield for given heterogeneous velocity model -------- */
void aco_FD_2D_vel(char *model, char *recei, int wn);

/* ------------------------------- Read in input parameters ---------------------------- */
WAVEPARA read_wave_para(char *input);

/* ---------------------------- Read in field data (nx by ny) -------------------------- */
void read_field_2D(char *input, int nx, int ny, float **field);

/* ------------------- Find maximum or minimum of input 2D field ----------------------- */
float max_or_min(float **field, int nx, int ny, int m);

#endif
