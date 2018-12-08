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
#include <time.h>   /* timing */

#define MIN(A, B) ((A) < (B) ? (A) : (B))

/* p, q, r, s specify size of domain */
#define p  0.0f
#define q  100.0f
#define r  0.0f
#define s  100.0f
/* water speed in x direction */
#define vx  0.5f
/* water speed in y direction */
#define vy  0.5f
/* cenx, ceny centre for Gaussian */
#define cenx ((p) + (q)) / 2
#define ceny ((r) + (s)) / 2
/* rad is Gaussian 'radius' */
#define rad 20.0f


/**
 * 2 dimensional Gaussian function over [p, q]x[r, s], height 1, centred on
 * (cenx, ceny) which becomes zero rad away from the centre.
 */
__device__ float gaussian2D(float x, float y)
{
    /* square of distance from centre */
    float d2 = (x - cenx) * (x - cenx) + (y - ceny) * (y - ceny);
    /* value of Gaussian function at (x, y) */
    return (d2 < rad * rad) ? exp(-0.01 * d2) : 0;
}

__global__ void iteration(float *u, float *u0, float Cx, float Cy, int NX, int NY)
{
    int ix0 = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int iy0 = blockIdx.y * blockDim.y + threadIdx.y + 1;

    int stride_x = blockDim.x * gridDim.x;
    int stride_y = blockDim.y * gridDim.y;

    for(int ix = ix0; ix < NX; ix += stride_x)
        for(int iy = iy0; iy < NY; iy += stride_y)
        {
            u[iy * NX + ix] = u0[iy * NX + ix]
                - Cx * (u0[iy * NX + ix] - u0[iy * NX + ix - 1])
                - Cy * (u0[iy * NX + ix] - u0[(iy - 1) * NX + ix]);
        }
}

__global__ void initArray(float *u0, int NX, int NY, float dx, float dy)
{
    int ix0 = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int iy0 = blockIdx.y * blockDim.y + threadIdx.y + 1;

    int stride_x = blockDim.x * gridDim.x;
    int stride_y = blockDim.y * gridDim.y;

    for(int ix = ix0; ix < NX; ix += stride_x)
        for(int iy = iy0; iy < NY; iy += stride_y)
            u0[iy * NX + ix] = gaussian2D(p + (ix - 1) * dx,
                                          r + (iy - 1) * dy);
}

__global__ void exactSolution(float *u, int NX, int NY, float dx, float dy, float t)
{
    int ix0 = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int iy0 = blockIdx.y * blockDim.y + threadIdx.y + 1;

    int stride_x = blockDim.x * gridDim.x;
    int stride_y = blockDim.y * gridDim.y;

    for(int ix = ix0; ix < NX; ix += stride_x)
        for(int iy = iy0; iy < NY; iy += stride_y)
            u[iy * NX + ix] = gaussian2D(p + (ix - 1) * dx - vx * t,
                                         r + (iy - 1) * dy - vy * t);
}

void linearadvectionFOU2D(int NX, int NY, float *u0, float *u)
{
    /* spatial step size in x */
    float dx = (q - p) / (NX - 1);
    /* spatial step size in x */
    float dy = (s - r) / (NY - 1);
    /* initial time */
    float t = 0;
    /* safety factor */
    float F = 0.4;
    /* heuristic time step calc */
    float dt  = F * MIN(dx / fabs(vx), dy / fabs(vy));
    /* number of time steps */
    int Ntimesteps = MIN(500, 20 / dt);
    /* Courant number in the x direction */
    float Cx = dt * vx / dx;
    /* Courant number in the y direction */
    float Cy = dt * vy / dy;

    dim3 number_of_blocks(15, 16, 1);
    dim3 threads_per_block(16, 16, 1);

    /* Increase NX and NY to create an extra row and col with zero
     * as boundary conditions. */
    NX++; NY++;
    /* initial u vector, extended array for 'ghost values' */
    float *temp;
    initArray<<<number_of_blocks, threads_per_block>>>(u0, NX, NY, dx, dy);
    cudaDeviceSynchronize();

    /* extended array for 'ghost values' mentioned below */
    for(int Ntimestep = 0; Ntimestep < Ntimesteps; Ntimestep++)
    {
        iteration<<<number_of_blocks, threads_per_block>>>(u, u0, Cx, Cy, NX, NY);
        cudaDeviceSynchronize();
        temp = u0, u0 = u, u = temp;
        t += dt;
    }
    /* store the analytic solution in u */
    exactSolution<<<number_of_blocks, threads_per_block>>>(u, NX, NY, dx, dy, t);

    printf("Final t = %f\n", t);
}

void writetofile(float *u0, float *exact, int NX, int NY)
{
    NX++; NY++;
    /* Output */
    FILE* u0file = fopen("u0.txt", "w");
    FILE* exactfile = fopen("exact.txt", "w");
    FILE* difffile = fopen("diff.txt", "w");
    float Eu0, Eexact;
    for(int i = 1; i < NX; i++)
    {
        for(int j = 1; j < NY; j++)
        {
            Eu0 = u0[j * NX + i];
            Eexact = exact[j * NX + i];
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
    time_t start = clock();

    int NX = 2048;        /* number of grid points in x direction */
    int NY = 2048;        /* number of grid points in y direction */

    float *ud0, *ud, *uh0, *uh;
    /*cudaMallocManaged(&ud0, NX * NY, sizeof(float));*/
    /*cudaMallocManaged(&ud, NX * NY, sizeof(float));*/
    cudaMalloc(&ud0, (NX + 1) * (NY + 1) * sizeof(float));
    cudaMalloc(&ud, (NX + 1) * (NY + 1) * sizeof(float));

    linearadvectionFOU2D(NX, NY, ud0, ud);

    /* copy back to host memory */
    uh = (float*)malloc((NX + 1) * (NY + 1) * sizeof(float));
    uh0 = (float*)malloc((NX + 1) * (NY + 1) * sizeof(float));
    cudaMemcpy(uh, ud, (NX + 1) * (NY + 1) * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(uh0, ud0, (NX + 1) * (NY + 1) * sizeof(float), cudaMemcpyDeviceToHost);

    int testx = (int)(NX * 0.6);
    int testy = (int)(NY * 0.6);
    printf("u[%d][%d] = %f\n", testx, testy,
            uh0[(testx + 1) * (NX + 1) + testy + 1]);

    time_t end = clock();
    printf("Time: %f ms.\n", 1000.0 * (end - start) / CLOCKS_PER_SEC);

    printf("Writing to files ...\n");
    writetofile(uh0, uh, NX, NY);
    printf("Done.\n");

    cudaFree(ud0);
    cudaFree(ud);
    free(uh0);
    free(uh);

    return 0;
}
