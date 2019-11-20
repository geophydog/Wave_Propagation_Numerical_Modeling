/* ---------------------------------------------------------------- */
/*   Initially coded by FENG Xuping @ SUSTech on Sept. 12th, 2019.  */
/* ---------------------------------------------------------------- */

#ifndef _SEIS_WAVE_H
#define _SEIS_WAVE_H

/* ++++++++++++++++ Function prototypes of solving acoustic wave equation ++++++++++++++ */
/* ------- Generate 1-D acoustic wavefield with finite difference method ------ */
void aco_1D(int nt, int nx, int sx, float dt, float f, float c, float **u);

/* ------- Generate 1-D acoustic wavefield with finite difference method ------ */
void aco_2D(int nt, int nx, int nz, int sx, int sz, float dt, float c, float f, float ***u);

/* ----- Generate 2-D acoustic wavefield for given heterogeneous velocity model ------- */
void aco_2D_vel(char *model, char *recei, int wn);

#endif
