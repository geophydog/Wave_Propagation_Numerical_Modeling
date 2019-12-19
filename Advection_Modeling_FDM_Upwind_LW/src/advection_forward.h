#ifndef _ADVECTION_FORWARD_H
#define _ADVECTION_FORWARD_H

/* --------------- Class of advection modeling for 1D case ------------------ */
class Advection1D{
    public:
        Advection1D();
        Advection1D(char *input);
        Advection1D(char *input, bool direction);
        ~Advection1D();

        void Initialization(void);
        void ReleaseMemory(void);
        void FDM(void);
        void Upwind(void);
        void LW(void);
        void WriteSnapshot(int tn, char *output);

    private:
        int nt = 501;
        double dt = 1.e-3;
        int nx = 201;
        double dx = 1.;
        int sx = 100;
        double v = 500.;
        double sigma = 15.;
        bool direction = true;
        double **p;
};

/* --------------- Class of advection modeling for 2D case ------------------ */
class Advection2D{
    public:
        Advection2D();
        Advection2D(char *input);
        ~Advection2D();

        void Initialization(void);
        void ReleaseMemory(void);
        void Gauss(void);
        void AttenuationSine(void);
        void FDM(void);
        void Upwind(void);
        void LW(void);
        void Forward(void);
        void WriteSnapshot(char *output);
    private:
        int nt = 501;
        double dt = 1.e-3;
        int nx = 101, ny = 101;
        double dx = 1., dy = 1.;
        int sx = 50, sy = 50;
        double vx = 500., vy = 500;
        double sigma_x = 15, sigma_y = 10.;
        int directions[2] = {true, true};
        int ns = 50;
        char method[64]={0};
        char SrcType[64] = {0};
        double ***p;
};
#endif
