/**
 * File: FDTD2Dsolver.c
 * 
 * Direct convertion from MATLAB script file linearadvectionFOU2D.m,
 * which is also uploaded to the project.
 * 
 * Following is the original file document from the MATLAB script:
 * 
 * Description: Solves the 2D linear advection equation
 * dU/dt + vx dU/dx + vy dU/dy = 0,
 * using first order forward difference in time
 * and first order backward differences in space.
 * The solution is calculated over [p, q] x [r, s] using NX, NY points
 * in the x and y directions respectively and plotted
 * after ntimesteps time steps, i.e. the final
 * solution is at time, ntimesteps*dt seconds.
 *
 * Boundary conditions: Dirichlet boundary conditions are used everywhere.
 *
 * Subfunction: gaussian2D
 *
 * Note:
 * This 2D problem has an analytical solution.
 * If initial condition is U(x,0) = f(x, y) then it can be shown
 * that the exact solution is U(x, y, t) = f(x - vx t, y - vy t).
 * i.e. the initial state is translated (advected) in the x and y directions
 * with speeds vx and vy respectively.
 *
 * Stability analysis for the timestep is very complicated so a heuristic
 * formula has been used based on the 1D case and a safety factor, F.
 */

#include <stdio.h>  /* file writing */
#include <string.h> /* memcpy */
#include <math.h>   /* exp */
#include <stdlib.h> /* malloc */

#define MIN(A, B) ((A) < (B) ? (A) : (B))

/* p, q, r, s specify size of domain */
const double p = 0;
const double q = 100;
const double r = 0;
const double s = 100;
const double vx = 0.5;      /* water speed in x direction */
const double vy = 0.5;      /* water speed in y direction */

/**
 * 2 dimensional Gaussian function over [p, q]x[r, s], height 1, centred on
 * (cenx, ceny) which becomes zero rad away from the centre.
 */
double gaussian2D(double x, double y)
{
    /* cenx, ceny centre for Gaussian */
    double cenx = (p + q) / 2;
    double ceny = (r + s) / 2;
    
    /* rad is Gaussian 'radius' */
    double rad = 20;

    /* square of distance from centre */
    double d2 = (x - cenx) * (x - cenx) + (y - ceny) * (y - ceny);
    
    /* value of Gaussian function at (x, y) */
    double z = (d2 < rad * rad) ? exp(-0.01 * d2) : 0;

    return z;
}

double* linearadvectionFOU2D(int NX, int NY, double* pt)
{
    /* spatial step size in x */
    double dx = (q - p) / (NX - 1);
    /* spatial step size in x */
    double dy = (s - r) / (NY - 1);
    /* Increase NX and NY to create an extra row and col with zero
     * as boundary conditions. */
    NX++; NY++;
    /* vector of grid points in x */
    double *x = (double*)malloc(NX * sizeof(double));
    for(int i = 0; i < NX - 1; i++)
        x[i + 1] = p + i * dx;
    /* vector of grid points in y */
    double *y = (double*)malloc(NY * sizeof(double));
    for(int i = 0; i < NY - 1; i++)
        y[i + 1] = r + i * dy;
    /* initial u vector, extended array for 'ghost values' */
    double *u0 = (double*)calloc(NX * NY, sizeof(double));
    for(int i = 1; i < NX; i++)
        for(int j = 1; j < NY; j++)
            u0[j * NX + i] = gaussian2D(x[i], y[j]);

    /* initial time */
    double t = 0;
    /* number of time steps */
    int Ntimesteps = 500;
    /* safety factor */
    double F = 0.4;
    /* heuristic time step calc */
    double dt  = F * MIN(dx / fabs(vx), dy / fabs(vy)); 
    
    /* Courant number in the x direction */
    double Cx = dt * vx / dx;
    /* Courant number in the y direction */
    double Cy = dt * vy / dy;
    
    /* define correct sized numerical solution array */
    double *u = (double*)calloc(NX * NY, sizeof(double));
    /* extended array for 'ghost values' mentioned below */
    for(int timecount = 0; timecount < Ntimesteps; timecount++)
    {//P
        t += dt;
        for(int i = 1; i < NX; i++)
            for(int j = 1; j < NY; j++)
                /*  FOU scheme in 2D */
                u[j * NX + i] = u0[j * NX + i]
                    - Cx * (u0[j * NX + i] - u0[j * NX + i - 1])
                    - Cy * (u0[j * NX + i] - u0[(j - 1) * NX + i]);
        /* copy solution to initial conditions for next iteration */
        memcpy(u0, u, NX * NY * sizeof(double));
    }

    *pt = t;
    free(x);
    free(y);
    free(u);

    return u0;
}

void writetofile(double *u0, double t, int NX, int NY)
{
    /* spatial step size in x */
    double dx = (q - p) / (NX - 1);
    /* spatial step size in x */
    double dy = (s - r) / (NY - 1);
    /* Increase NX and NY to create an extra row and col with zero
     * as boundary conditions. */
    NX++; NY++;
    /* Output */
    FILE* u0file = fopen("u0.txt", "w");
    FILE* exactfile = fopen("exact.txt", "w");
    FILE* difffile = fopen("diff.txt", "w");
    double Eu0, Eexact, x, y;
    for(int i = 1; i < NX; i++)
    {
        for(int j = 1; j < NY; j++)
        {
            x = p + (i - 1) * dx;
            y = r + (j - 1) * dy;
            Eu0 = u0[j * NX + i];
            Eexact = gaussian2D(x - vx * t, y - vy * t);
            fprintf(u0file, "%10.6f ", Eu0); 
            fprintf(exactfile, "%10.6f ", Eexact);
            fprintf(difffile, "%10.6f ", Eu0 - Eexact);
        }
        fprintf(u0file, "\n");
        fprintf(exactfile, "\n");
        fprintf(difffile, "\n");
    }
    fclose(u0file);
    fclose(exactfile);
    fclose(difffile);
}

int main()
{
    int NX = 2048;        /* number of grid points in x direction */
    int NY = 2048;        /* number of grid points in y direction */

    double *u0, t;

    u0 = linearadvectionFOU2D(NX, NY, &t);
    printf("Writing to files ...\n");
    writetofile(u0, t, NX, NY);
    printf("Done.\n");

    free(u0);

    return 0;
}
