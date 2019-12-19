#include <iostream>
#include <fstream>
#include <cmath>
#include <string.h>
#include "advection_forward.h"

using namespace std;

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* --------------------- Advection Modeling for 1D Case --------------------------- */
Advection1D:: Advection1D(){}

Advection1D:: Advection1D(char *input){
    ifstream fin;
    char ss[128];
    int index = 1;
    try{
        fin.open(input, ios:: in);
        while(fin.getline(ss, 128)){
            switch(index){
                case 1:
                    sscanf(ss, "%*s %d %*s", &nt); break;
                case 2:
                    sscanf(ss, "%*s %lf %*s", &dt); break;
                case 3:
                    sscanf(ss, "%*s %d %*s", &nx); break;
                case 4:
                    sscanf(ss, "%*s %lf %*s", &dx); break;
                case 5:
                    sscanf(ss, "%*s %d %*s", &sx); break;
                case 6:
                    sscanf(ss, "%*s %lf %*s", &v); break;
                case 7:
                    sscanf(ss, "%*s %lf %*s", &sigma); break;
                default:
                    cout << "Ignore the other lines!" << endl; break;
            }
            index ++;
        }
    } catch( int e) {
        cout << "Canot open file: " << input << "!!!" << endl;
    }
    fin.close();
}

Advection1D:: Advection1D(char *input, bool direction){
    Advection1D{input};
    direction = direction;
}

Advection1D:: ~Advection1D(){
    cout << "Delete class Advection1D!" << endl;
}

// Initialization for array p.
void Advection1D:: Initialization(void){
    cout << "Initialization of array  ..." << endl;
    int i = 0, j = 0;

    p = new double *[nt];
    for(i = 0; i < nt; i ++)
        p[i] = new double[nx];
    // Initialization of array p.
    for(i = 0; i < nt; i ++)
        for(j = 0; j < nx; j ++)
            p[i][j] = 0.;
    // Initial condition.
    for(j = 0; j < nx; j ++)
        p[0][j] = exp(-pow((j*dx-sx*dx)/sigma, 2));
}

// Release dynamic memory for array p.
void Advection1D:: ReleaseMemory(void){
    cout << "Release dynamic memory of array ..." << endl;
    int i = 0;
    for(i = 0; i < nt; i ++)
        delete [] p[i];
    delete []p;
}

// Advection modeling with FDM with 2nd-order for space.
void Advection1D:: FDM(void){
    int i = 0, j = 0;

    // Initialization.
    Initialization();

    // Advection modeling.
    if(!direction)
        v = -v;
    for(i = 0; i < nt-1; i ++)
        for(j = 1; j < nx-1; j ++)
            p[i+1][j] = p[i][j] - v*dt/dx/2. * (p[i][j+1]-p[i][j-1]);
}

// Advection modeling with upwind FVM.
void Advection1D:: Upwind(void){
    int i = 0, j = 0;
    //Initialization.
    Initialization();
    if(!direction)
        v = -v;
    for(i = 0; i < nt-1; i ++)
        for(j = 1; j < nx; j++)
            p[i+1][j] = p[i][j] - v*dt/dx*(p[i][j]-p[i][j-1]);
}

// Advection modeling with Lax-Wendroff FVM.
void Advection1D:: LW(void){
    // Initialization.
    Initialization();
    cout << "Advection modeling for 1D case with Lax-Wendroff ..." << endl;
    int i = 0, j = 0;
    if(!direction)
        v = -v;
    for(i = 0; i < nt-1; i ++)
        for(j = 1; j < nx-1; j ++) {
            p[i+1][j] = p[i][j] - v*dt/dx*(p[i][j]-p[i][j-1]) + \
                        0.5*pow(v*dt/dx, 2.) * (p[i][j+1]-2.*p[i][j]+p[i][j-1]);
        }
}

void Advection1D:: WriteSnapshot(int tn, char *output){
    cout << "Write out snapshots ..." << endl;
    ofstream fout;
    int i = 0, j = 0;
    try{
        fout.open(output, ios:: out);
        for(i = 0; i < nt; i ++)
            for(j = 0; j < nx; j ++){
                if ((j+1) == nx)
                    fout << p[i][j] << endl;
                else
                    fout << p[i][j] << " ";
            }
    } catch(int e){
        cout << "Cannot open file: " <<  output << "!!!" << endl;
    }
    fout.close();
    // release memory.
    ReleaseMemory();
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* ------------------------ Advection Modeling for 2D Case ------------------------ */
Advection2D:: Advection2D(){}

Advection2D:: Advection2D(char *input){
    ifstream fin;
    int index = 1;
    char ss[128]={0};
    try{
        fin.open(input, ios:: in);
        while(!fin.eof()){
            fin.getline(ss, 128, '\n');
            switch(index){
                case 1:
                sscanf(ss, "%*s %d %*s", &nt); break;
                case 2:
                sscanf(ss, "%*s %lf %*s", &dt); break;
                case 3:
                sscanf(ss, "%*s %d %*s", &nx); break;
                case 4:
                sscanf(ss, "%*s %lf %*s", &dx); break;
                case 5:
                sscanf(ss, "%*s %d %*s", &ny); break;
                case 6:
                sscanf(ss, "%*s %lf %*s", &dy); break;
                case 7:
                sscanf(ss, "%*s %d %*s", &sx); break;
                case 8:
                sscanf(ss, "%*s %d %*s", &sy); break;
                case 9:
                sscanf(ss, "%*s %lf %*s", &vx); break;
                case 10:
                sscanf(ss, "%*s %lf %*s", &vy); break;
                case 11:
                sscanf(ss, "%*s %lf %*s", &sigma_x); break;
                case 12:
                sscanf(ss, "%*s %lf %*s", &sigma_y); break;
                case 13:
                sscanf(ss, "%*s %d %*s", &directions[0]); break;
                case 14:
                sscanf(ss, "%*s %d %*s", &directions[1]); break;
                case 15:
                sscanf(ss, "%*s %d %*s", &ns); break;
                case 16:
                sscanf(ss, "%*s %s %*s", method); break;
                case 17:
                sscanf(ss, "%*s %s %*s", SrcType); break;
                default:
                cout << "Ignore the other line(s)!" << endl; break;
            }
            index ++;
        }
    } catch( int e){
        cout << "Cannot open file: " << input << "!!!" << endl;
    }
    fin.close();
}

Advection2D:: ~Advection2D(){
    cout << "Delete class of Advection2D!" << endl;
}


void Advection2D:: Gauss(void){
    int j = 0, k = 0;
    // Initial condition.
    double tmp = 0.;
    for(j = 0; j < ny; j ++)
        for(k = 0; k < nx; k ++) {
            tmp = -pow((k*dx-sx*dx)/sigma_x, 2.) - pow((j*dy-sy*dy)/sigma_y, 2.);
            p[0][j][k] = exp(tmp);
        }
}

void Advection2D:: AttenuationSine(void){
    int j = 0, k = 0;
    double tmp1 = 0., tmp2 = 0.;
    // Initial condition.
    for(j = 0; j < ny; j ++)
        for(k = 0; k < nx; k ++){
            tmp1 = pow((k*dx-sx*dx)/sigma_x, 2.) + pow((j*dy-sy*dy)/sigma_y, 2.);
            tmp2 = -pow((k*dx-sx*dx)/sigma_x, 2.) - pow((j*dy-sy*dy)/sigma_y, 2.);
            p[0][j][k] = sin(tmp1) * exp(tmp2); 
        }
}

void Advection2D:: Initialization(void){
    cout << "Initialization of array p in Advection2D ..." << endl;
    int i = 0, j = 0, k = 0;
    // Allocate dynamic memory.
    p = new double **[nt];
    for(i = 0; i < nt; i ++){
        p[i] = new double *[ny];
        for(j = 0; j < ny; j ++)
            p[i][j] = new double [nx];
    }
    // Setting default values.
    for(i = 1; i < nt; i ++)
        for(j = 0; j < ny; j ++)
            for(k = 0; k < nx; k ++)
                p[i][j][k] = 0.;

    // Initial condition.
    if(strcmp("Gauss", SrcType))
        Gauss();
    else
        AttenuationSine();
}

void Advection2D:: ReleaseMemory(void){
    cout << "Release dynamic memory of array p in Advection2D ..." << endl;
    int i = 0, j = 0;
    for(i = 0; i < nt; i ++){
        for(j = 0; j < ny; j ++)
            delete []p[i][j];
        delete []p[i];
    }
    delete []p;
}

void Advection2D:: FDM(void){
    double xterm = 0., yterm = 0.;
    int i = 0, j = 0, k = 0;

    // Initialization.
    Initialization();
    // Advection directions.
    if(!directions[0])
        vx = -vx;
    if(!directions[1])
        vy = -vy;
    cout << "Advection modeling for 2D case with FDM ..." << endl;
    // Advection modeling.
    for(i = 0; i < nt-1; i ++)
        for(j = 1; j < ny-1; j ++)
            for(k = 1; k < nx-1; k ++){
                xterm = (p[i][j][k+1]-p[i][j][k-1])/dx/2.;
                yterm = (p[i][j+1][k]-p[i][j-1][k])/dy/2.;
                p[i+1][j][k] = p[i][j][k] - vx*dt*xterm - vy*dt*yterm;
            }
}

void Advection2D:: Upwind(void){
    double xterm = 0., yterm = 0.;
    int i = 0, j = 0, k = 0;
    
    // Initialization.
    Initialization();
    // Advection directions.
    if(!directions[0])
        vx = -vx;
    if(!directions[1])
        vy = -vy;
    cout << "Advection modeling for 2D case with upwind ..." << endl;
    // Advection modeling.
    for(i = 0; i < nt-1; i ++)
        for(j = 1; j < ny; j ++)
            for(k = 1; k < nx; k ++){
                xterm = (p[i][j][k]-p[i][j][k-1])/dx;
                yterm = (p[i][j][k]-p[i][j-1][k])/dy;
                p[i+1][j][k] = p[i][j][k] - vx*dt*xterm - vy*dt*yterm;
            }
}

void Advection2D:: LW(void){
    double xterm1 = 0., xterm2 = 0., yterm1 = 0., yterm2 = 0.;
    int i = 0, j = 0, k = 0;

    // Initialization.
    Initialization();
    // Advection directions.
    if(!directions[0])
        vx = -vx;
    if(!directions[1])
        vy = -vy;
    // Advection modeling.
    cout << "Advection modeling with L-W for 2D case ..." << endl;
    for(i = 0; i < nt-1; i ++)
        for(j = 1; j < ny-1; j ++)
            for(k = 1; k < nx-1; k ++){
                xterm1 = (p[i][j][k+1]-p[i][j][k-1])/dx/2.;
                yterm1 = (p[i][j+1][k]-p[i][j-1][k])/dy/2.;
                xterm2 = (p[i][j][k+1]-2*p[i][j][k]+p[i][j][k-1])/dx/dx;
                xterm2 = (p[i][j+1][k]-2*p[i][j][k]+p[i][j-1][k])/dy/dy;
                p[i+1][j][k] = p[i][j][k] - vx*dt*xterm1 - vy*dt*yterm1 +\
                               0.5*pow(vx*dt, 2.)*xterm2 + 0.5*pow(vy*dt, 2.)*yterm2;
            }
}

void Advection2D:: Forward(void){
    char m1[] = "FDM", m2[]="Upwind", m3[]="LW";
    if(strcmp(m1, method))
        FDM();
    else if(strcmp(m2, method))
        Upwind();
    else if(strcmp(m3, method))
        LW();
    else
        Upwind();
}

void Advection2D:: WriteSnapshot(char *output){
    int i = 0, j = 0, k = 0;
    char ss[256];
    int tt = nt / ns;
    ofstream fout;
    system("mkdir -p Advection2D_Snapshots");
    system("rm Advection2D_Snapshots/*");
    for(i = 0; i <= ns; i ++){
        sprintf(ss, "Advection2D_Snapshots/advection2D_%lf_%s", i*tt*dt, output);
        try{
            cout << "Write out file: " << ss << "..." << endl;
            fout.open(ss, ios:: out);
            for(j = 0; j < ny; j ++)
                for(k = 0; k < nx; k ++){
                    if((k+1) == nx)
                        fout << p[i*tt][j][k] << endl;
                    else
                        fout << p[i*tt][j][k] << " ";
                }
        } catch( int e){
            cout << "Cannot open file: " << ss << "!!!" << endl;
        }
        fout.close();
    }

    // Release dynamic memory of array p.
    ReleaseMemory();
}
