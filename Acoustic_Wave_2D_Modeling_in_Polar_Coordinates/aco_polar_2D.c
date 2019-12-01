#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265

void aco_polar_2D_homo(char *input, char *output);

int main(int argc, char *argv[]){
    if(argc != 3) {
        fprintf(stderr, "Usage: aco_polar_2D input output\n");
        exit(1);
    }

    aco_polar_2D_homo(argv[1], argv[2]);
    return 0;
}

void aco_polar_2D_homo(char *input, char *output){
    int i=0, j=0, k=0, sr=0, sth=0, pn=0, nt, nr, nth, index=0;
    float dt, dr, dth, f, r1, r2, th1, th2, c=0;
    char ss[128];

    FILE *fin;
    if(fin=fopen(input, "r")){
        while(fgets(ss, 128, fin)){
            if(0==index){
                sscanf(ss, "%d %f %f %f", &nt, &dt, &f, &c);
            }else if(1==index){
                sscanf(ss, "%d %f %f %f %d %f %f %f", &nr, &r1, &r2, &dr, &nth, &th1, &th2, &dth);
            }else if(2==index){
                sscanf(ss, "%d %d %d", &sr, &sth, &pn);
            } else{
                break;
            }
            index ++;
        }
    } else{
        fprintf(stderr, "Cannot open file %s!!!\n", input);
        return;
    }
    fclose(fin);

    // Ricker wavelet source.
    float *s;
    s = (float *) malloc(sizeof(float) * nt);
    for(i = 0; i < nt; i ++)
        s[i] = (1.-2*pow(PI*f*(i*dt-1/f), 2.)) * exp(-pow(PI*f*(i*dt-1/f), 2.));

    // Allocate dynamic memory of acoustic pressure matrix p.
    float ***p;
    p = (float ***) malloc(sizeof(float **) * nt);
    for(i = 0; i < nt; i++){
        p[i] = (float **) malloc(sizeof(float *) * nr);
        for(j = 0; j < nr; j++)
            p[i][j] = (float *) malloc(sizeof(float) * nth);
    }
    // Initializtion of acoutic pressure matrix p.
    for(j = 0;j < nr; j++)
        for(k = 0; k < nth; k++){
            p[0][j][k] = 0.;
            p[nt-1][j][k] = 0.;
        }
    for(i = 0; i < nt; i ++)
        for(j = 0; j < nr; j++){
            p[i][j][0] = 0.;
            p[i][j][nth-1] = 0.;
        }
    for(i = 0; i < nt; i ++)
        for(k = 0; k < nth; k ++){
            p[i][0][k] = 0.;
            p[i][nr-1][k] = 0.;
        }

    // Modeling.
    float cc = 0.;
    for(i = 1; i < nt-1; i ++)
        for(j = 1; j < nr-1; j ++)
            for(k = 1; k < nth-1; k++){
                if(j < nr/3)
                    cc = c * 1.5;
                else if(j >= nr*2/3)
                    cc = c * 0.75;
                else
                    cc = c;

                if(sr==j && sth==k)
                    p[i+1][j][k] = s[i];
                else
                    p[i+1][j][k] = 2.*p[i][j][k] - p[i-1][j][k] + pow(c*dt, 2.) * ( \
                            (p[i][j+1][k]-2*p[i][j][k]+p[i][j-1][k])/pow(dr, 2.) + \
                            (p[i][j+1][k]-p[i][j-1][k])/(2.*dr*(r1+j*dr)) + \
                            (p[i][j][k+1]-2*p[i][j][k]+p[i][j][k-1])/pow((r1+j*dr)*dth, 2.) );
            }
    // Output acoustic pressure fieid.
    FILE *fout;
    if(fout=fopen(output, "w")){
        fprintf(fout, "%f %f %f %f %f %f\n", r1, r2, th1, th2, sr*dr+r1, sth*dth+th1);
        for(j = 0; j < nr; j ++)
            for(k = 0; k < nth; k ++){
                if((k+1)==nth)
                    fprintf(fout, "%f\n", p[pn][j][k]);
                else
                    fprintf(fout, "%f ", p[pn][j][k]);
            }
    } else{
        fprintf(stderr, "Cannot open file %s!!!\n", output);
        return;
    }
    fclose(fout);
}
