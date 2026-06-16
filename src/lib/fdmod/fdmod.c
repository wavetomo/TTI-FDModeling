/********** **********/
#include <fdmod.h>
/**/
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/**/
#include <omp.h>
/**/
#include <cwp.h>
#include <par.h>
/********** **********/

/**
 * the 8th order accuracy finite difference to calculate the first derivative
 *  the 8th order accuracy finite difference to calculate the first
 * derivative
 * @param float *dat_fd: out 1D array
 * @param float *dat: input 1D array
 * @param float d: the samping interval of data
 * @param int n: number of samples in 1st dimension
 * @return {*}
 * @LastEditors: LaurentCheung
 * @Version: 1.0
 */
void derivate1_fd1d(float *dat_fd, float *dat, float d, int n)
{
  int i;
  float invd, C81, C82, C83, C84;
  invd = 1.0f / d;
  C81 = 0.839725908387582f * invd;
  C82 = -0.243216342043984f * invd;
  C83 = 0.059588302317624f * invd;
  C84 = -0.008074531266178f * invd;
#pragma omp parallel for private(i) firstprivate(C81, C82, C83, C84, invd)
  for (i = 4; i < n - 4; i++)
  {
    dat_fd[i] = C81 * (dat[i + 1] - dat[i - 1]) + C82 * (dat[i + 2] - dat[i - 2]) +
                C83 * (dat[i + 3] - dat[i - 3]) + C84 * (dat[i + 4] - dat[i - 4]);
  }
#pragma omp parallel for private(i) firstprivate(invd)
  for (i = 0; i < 4; i++)
  {
    dat_fd[i] = (dat[i + 1] - dat[i]);
    dat_fd[i] *= invd;
  }
#pragma omp parallel for private(i) firstprivate(invd)
  for (i = n - 4; i < n; i++)
  {
    dat_fd[i] = (dat[i] - dat[i - 1]);
    dat_fd[i] *= invd;
  }
}
/**
 * the 8th order accuracy finite difference to calculate the first derivative
 *  the 8th order accuracy finite difference to calculate the first
 * derivative
 * @param float **dat_fd: out 2D array
 * @param float **dat: input 2D array
 * @param float d: the interval of the calculation direction
 * @param int nz: number of samples in 1st dimension
 * @param int nx: number of samples in 2nd dimension
 * @param int dim: the direction for calculating directive, =1: calcualte
 * derivative along z (vertical), =2: calculate derivative along x (horizontal);
 * @return {*}
 * @LastEditors: LaurentCheung
 * @Version: 1.0
 */
void derivate1_fd2d(float **dat_fd, float **dat, float d, int nz, int nx, int dim)
{
  int i, j;
  float invd, C81, C82, C83, C84;
  invd = 1.0f / d;
  C81 = 0.839725908387582f * invd;
  C82 = -0.243216342043984f * invd;
  C83 = 0.059588302317624f * invd;
  C84 = -0.008074531266178f * invd;
  switch (dim)
  {
  case 1:
#pragma omp parallel for private(i, j) firstprivate(C81, C82, C83, C84, invd)
    for (i = 0; i < nx; i++)
    {
#pragma ivdep
      for (j = 4; j < nz - 4; j++)
      {
        dat_fd[i][j] =
            C81 * (dat[i][j + 1] - dat[i][j - 1]) + C82 * (dat[i][j + 2] - dat[i][j - 2]) +
            C83 * (dat[i][j + 3] - dat[i][j - 3]) + C84 * (dat[i][j + 4] - dat[i][j - 4]);
      }
#pragma ivdep
      for (j = 0; j < 4; j++)
      {
        dat_fd[i][j] = (dat[i][j + 1] - dat[i][j]);
        dat_fd[i][j] *= invd;
      }
#pragma ivdep
      for (j = nz - 4; j < nz; j++)
      {
        dat_fd[i][j] = (dat[i][j] - dat[i][j - 1]);
        dat_fd[i][j] *= invd;
      }
    }
    break;
  case 2:
#pragma omp parallel for private(i, j) firstprivate(C81, C82, C83, C84, invd)
    for (j = 0; j < nz; j++)
    {
#pragma ivdep
      for (i = 4; i < nx - 4; i++)
      {
        dat_fd[i][j] =
            C81 * (dat[i + 1][j] - dat[i - 1][j]) + C82 * (dat[i + 2][j] - dat[i - 2][j]) +
            C83 * (dat[i + 3][j] - dat[i - 3][j]) + C84 * (dat[i + 4][j] - dat[i - 4][j]);
      }
#pragma ivdep
      for (i = 0; i < 4; i++)
      {
        dat_fd[i][j] = (dat[i + 1][j] - dat[i][j]);
        dat_fd[i][j] *= invd;
      }
#pragma ivdep
      for (i = nx - 4; i < nx; i++)
      {
        dat_fd[i][j] = (dat[i][j] - dat[i - 1][j]);
        dat_fd[i][j] *= invd;
      }
    }
    break;
  }
}
/**
 * the 8th order accuracy finite difference to calculate the first derivative
 *  the 8th order accuracy finite difference to calculate the first
 * derivative
 * @param float ***dat_fd: out 2D array
 * @param float ***dat: input 2D array
 * @param float d: the interval of the calculation direction
 * @param int nz: number of samples in 1st dimension
 * @param int nx: number of samples in 2nd dimension
 * @param int ny: number of samples in 3rd dimension
 * @param int dim: the direction for calculating directive, =1: calcualte
 * derivative along z (vertical), =2: calculate derivative along x
 * (horizontal),=3: calculate derivative along y (horizontal);
 * @return {*}
 * @LastEditors: LaurentCheung
 * @Version: 1.0
 */
void derivate1_fd3d(float ***dat_fd, float ***dat, float d, int nz, int nx, int ny, int dim)
{
  int i, j, k;
  float invd, C81, C82, C83, C84;
  invd = 1.0f / d;
  C81 = 0.839725908387582f * invd;
  C82 = -0.243216342043984f * invd;
  C83 = 0.059588302317624f * invd;
  C84 = -0.008074531266178f * invd;
  switch (dim)
  {
  case 1:
#pragma omp parallel for private(i, j, k) firstprivate(C81, C82, C83, C84, invd)
    for (i = 0; i < ny; i++)
    {
      for (j = 0; j < nx; j++)
      {
#pragma ivdep
        for (k = 4; k < nz - 4; k++)
        {
          dat_fd[i][j][k] = C81 * (dat[i][j][k + 1] - dat[i][j][k - 1]) +
                            C82 * (dat[i][j][k + 2] - dat[i][j][k - 2]) +
                            C83 * (dat[i][j][k + 3] - dat[i][j][k - 3]) +
                            C84 * (dat[i][j][k + 4] - dat[i][j][k - 4]);
        }
#pragma ivdep
        for (k = 0; k < 4; k++)
        {
          dat_fd[i][j][k] = (dat[i][j][k + 1] - dat[i][j][k]);
          dat_fd[i][j][k] *= invd;
        }
#pragma ivdep
        for (k = nz - 4; k < nz; k++)
        {
          dat_fd[i][j][k] = (dat[i][j][k] - dat[i][j][k - 1]);
          dat_fd[i][j][k] *= invd;
        }
      }
    }
    break;
  case 2:
#pragma omp parallel for private(i, j, k) firstprivate(C81, C82, C83, C84, invd)
    for (i = 0; i < ny; i++)
    {
      for (k = 0; k < nz; k++)
      {
#pragma ivdep
        for (j = 4; j < nx - 4; j++)
        {
          dat_fd[i][j][k] = C81 * (dat[i][j + 1][k] - dat[i][j - 1][k]) +
                            C82 * (dat[i][j + 2][k] - dat[i][j - 2][k]) +
                            C83 * (dat[i][j + 3][k] - dat[i][j - 3][k]) +
                            C84 * (dat[i][j + 4][k] - dat[i][j - 4][k]);
        }
#pragma ivdep
        for (j = 0; j < 4; j++)
        {
          dat_fd[i][j][k] = (dat[i][j + 1][k] - dat[i][j][k]);
          dat_fd[i][j][k] *= invd;
        }
#pragma ivdep
        for (j = nx - 4; j < nx; j++)
        {
          dat_fd[i][j][k] = (dat[i][j][k] - dat[i][j - 1][k]);
          dat_fd[i][j][k] *= invd;
        }
      }
    }
    break;
  case 3:
#pragma omp parallel for private(i, j, k) firstprivate(C81, C82, C83, C84, invd)
    for (j = 0; j < nx; j++)
    {
      for (k = 0; k < nz; k++)
      {
#pragma ivdep
        for (i = 4; i < ny - 4; i++)
        {
          dat_fd[i][j][k] = C81 * (dat[i + 1][j][k] - dat[i - 1][j][k]) +
                            C82 * (dat[i + 2][j][k] - dat[i - 2][j][k]) +
                            C83 * (dat[i + 3][j][k] - dat[i - 3][j][k]) +
                            C84 * (dat[i + 4][j][k] - dat[i - 4][j][k]);
        }
#pragma ivdep
        for (i = 0; i < 4; i++)
        {
          dat_fd[i][j][k] = (dat[i + 1][j][k] - dat[i][j][k]);
          dat_fd[i][j][k] *= invd;
        }
#pragma ivdep
        for (i = ny - 4; i < ny; i++)
        {
          dat_fd[i][j][k] = (dat[i][j][k] - dat[i - 1][j][k]);
          dat_fd[i][j][k] *= invd;
        }
      }
    }
    break;
  }
}
/**
 *
 * @param float **gx
 * @param float **gz
 * @param float **dat
 * @param float dz
 * @param float dx
 * @param int nz
 * @param int nx
 * @return {*}
 * @LastEditors: LaurentCheung
 * @Version: 1.0
 */
void gradient2d(float **gz, float **gx, float **dat, float dz, float dx, int nz, int nx)
{
  int i, j;
  float invdx, invdz, C81, C82, C83, C84;
  invdx = 1.0f / dx;
  invdz = 1.0f / dz;
  C81 = 0.839725908387582f;
  C82 = -0.243216342043984f;
  C83 = 0.059588302317624f;
  C84 = -0.008074531266178f;
#pragma omp parallel for private(i, j) firstprivate(C81, C82, C83, C84, invdz)
  for (i = 0; i < nx; i++)
  {
#pragma ivdep
    for (j = 4; j < nz - 4; j++)
    {
      gz[i][j] = C81 * (dat[i][j + 1] - dat[i][j - 1]) + C82 * (dat[i][j + 2] - dat[i][j - 2]) +
                 C83 * (dat[i][j + 3] - dat[i][j - 3]) + C84 * (dat[i][j + 4] - dat[i][j - 4]);
      gz[i][j] *= invdz;
    }
#pragma ivdep
    for (j = 0; j < 4; j++)
    {
      gz[i][j] = (dat[i][j + 1] - dat[i][j]);
      gz[i][j] *= invdz;
    }
#pragma ivdep
    for (j = nz - 4; j < nz; j++)
    {
      gz[i][j] = (dat[i][j] - dat[i][j - 1]);
      gz[i][j] *= invdz;
    }
  }
#pragma omp parallel for private(i, j) firstprivate(C81, C82, C83, C84, invdx)
  for (j = 0; j < nz; j++)
  {
#pragma ivdep
    for (i = 4; i < nx - 4; i++)
    {
      gx[i][j] = C81 * (dat[i + 1][j] - dat[i - 1][j]) + C82 * (dat[i + 2][j] - dat[i - 2][j]) +
                 C83 * (dat[i + 3][j] - dat[i - 3][j]) + C84 * (dat[i + 4][j] - dat[i - 4][j]);
      gx[i][j] *= invdx;
    }
#pragma ivdep
    for (i = 0; i < 4; i++)
    {
      gx[i][j] = (dat[i + 1][j] - dat[i][j]);
      gx[i][j] *= invdx;
    }
#pragma ivdep
    for (i = nx - 4; i < nx; i++)
    {
      gx[i][j] = (dat[i][j] - dat[i - 1][j]);
      gx[i][j] *= invdx;
    }
  }
}
/**
 *
 * @param float ***gz:
 * @param float ***gx:
 * @param float ***gy:
 * @param float ***dat:
 * @param float dz:
 * @param float dx:
 * @param float dy:
 * @param int nz:
 * @param int nx:
 * @param int ny:
 * @return {*}
 * @LastEditors: LaurentCheung
 * @Version: 1.0
 */
void gradient3d(float ***gz, float ***gx, float ***gy, float ***dat, float dz, float dx, float dy,
                int nz, int nx, int ny)
{
  int i, j, k;
  float invdx, invdy, invdz, C81, C82, C83, C84;
  invdx = 1.0f / dx;
  invdy = 1.0f / dy;
  invdz = 1.0f / dz;
  C81 = 0.839725908387582f;
  C82 = -0.243216342043984f;
  C83 = 0.059588302317624f;
  C84 = -0.008074531266178f;
#pragma omp parallel for private(i, j, k) firstprivate(C81, C82, C83, C84, invdz)
  for (i = 0; i < ny; i++)
  {
    for (j = 0; j < nx; j++)
    {
#pragma ivdep
      for (k = 4; k < nz - 4; k++)
      {
        gx[i][j][k] = C81 * (dat[i][j][k + 1] - dat[i][j][k - 1]) +
                      C82 * (dat[i][j][k + 2] - dat[i][j][k - 2]) +
                      C83 * (dat[i][j][k + 3] - dat[i][j][k - 3]) +
                      C84 * (dat[i][j][k + 4] - dat[i][j][k - 4]);
        gx[i][j][k] *= invdz;
      }
#pragma ivdep
      for (k = 0; k < 4; k++)
      {
        gx[i][j][k] = (dat[i][j][k + 1] - dat[i][j][k]);
        gx[i][j][k] *= invdz;
      }
#pragma ivdep
      for (k = nz - 4; k < nz; k++)
      {
        gx[i][j][k] = (dat[i][j][k] - dat[i][j][k - 1]);
        gx[i][j][k] *= invdz;
      }
    }
  }
#pragma omp parallel for private(i, j, k) firstprivate(C81, C82, C83, C84, invdx)
  for (i = 0; i < ny; i++)
  {
    for (k = 0; k < nz; k++)
    {
#pragma ivdep
      for (j = 4; j < nx - 4; j++)
      {
        gx[i][j][k] = C81 * (dat[i][j + 1][k] - dat[i][j - 1][k]) +
                      C82 * (dat[i][j + 2][k] - dat[i][j - 2][k]) +
                      C83 * (dat[i][j + 3][k] - dat[i][j - 3][k]) +
                      C84 * (dat[i][j + 4][k] - dat[i][j - 4][k]);
        gx[i][j][k] *= invdx;
      }
#pragma ivdep
      for (j = 0; j < 4; j++)
      {
        gx[i][j][k] = (dat[i][j + 1][k] - dat[i][j][k]);
        gx[i][j][k] *= invdx;
      }
#pragma ivdep
      for (j = nx - 4; j < nx; j++)
      {
        gx[i][j][k] = (dat[i][j][k] - dat[i][j - 1][k]);
        gx[i][j][k] *= invdx;
      }
    }
  }
#pragma omp parallel for private(i, j, k) firstprivate(C81, C82, C83, C84, invdy)
  for (j = 0; j < nx; j++)
  {
    for (k = 0; k < nz; k++)
    {
#pragma ivdep
      for (i = 4; i < ny - 4; i++)
      {
        gy[i][j][k] = C81 * (dat[i + 1][j][k] - dat[i - 1][j][k]) +
                      C82 * (dat[i + 2][j][k] - dat[i - 2][j][k]) +
                      C83 * (dat[i + 3][j][k] - dat[i - 3][j][k]) +
                      C84 * (dat[i + 4][j][k] - dat[i - 4][j][k]);
        gy[i][j][k] *= invdy;
      }
#pragma ivdep
      for (i = 0; i < 4; i++)
      {
        gy[i][j][k] = (dat[i + 1][j][k] - dat[i][j][k]);
        gy[i][j][k] *= invdy;
      }
#pragma ivdep
      for (i = ny - 4; i < ny; i++)
      {
        gy[i][j][k] = (dat[i][j][k] - dat[i - 1][j][k]);
        gy[i][j][k] *= invdy;
      }
    }
  }
}
/**
 *
 * @param float **div:
 * @param float **Fx:
 * @param float **Fz:
 * @param float dz:
 * @param float dx:
 * @param int nz:
 * @param int nx:
 * @return {*}
 * @LastEditors: LaurentCheung
 * @Version: 1.0
 */
void divergence2d(float **div, float **Fx, float **Fz, float dz, float dx, int nz, int nx)
{
  int i, j;
  float invdx, invdz, fz, C81, C82, C83, C84;
  invdx = 1.0f / dx;
  invdz = 1.0f / dz;
  C81 = 0.839725908387582f;
  C82 = -0.243216342043984f;
  C83 = 0.059588302317624f;
  C84 = -0.008074531266178f;
#pragma omp parallel for private(i, j) firstprivate(C81, C82, C83, C84, invdx)
  for (i = 0; i < nx; i++)
  {
#pragma ivdep
    for (j = 4; j < nz - 4; j++)
    {
      div[i][j] = C81 * (Fz[i][j + 1] - Fz[i][j - 1]) + C82 * (Fz[i][j + 2] - Fz[i][j - 2]) +
                  C83 * (Fz[i][j + 3] - Fz[i][j - 3]) + C84 * (Fz[i][j + 4] - Fz[i][j - 4]);
      div[i][j] *= invdz;
    }
    for (j = 0; j < 4; j++)
    {
      div[i][j] = (Fz[i][j + 1] - Fz[i][j]);
      div[i][j] *= invdz;
    }
    for (j = nz - 4; j < nz; j++)
    {
      div[i][j] = (Fz[i][j] - Fz[i][j - 1]);
      div[i][j] *= invdz;
    }
  }
#pragma omp parallel for private(i, j, fz) firstprivate(C81, C82, C83, C84, invdz)
  for (j = 0; j < nz; j++)
  {
#pragma ivdep
    for (i = 4; i < nx - 4; i++)
    {
      fz = C81 * (Fx[i + 1][j] - Fx[i - 1][j]) + C82 * (Fx[i + 2][j] - Fx[i - 2][j]) +
           C83 * (Fx[i + 3][j] - Fx[i - 3][j]) + C84 * (Fx[i + 4][j] - Fx[i - 4][j]);
      fz *= invdx;
      div[i][j] += fz;
    }
    for (i = 0; i < 4; i++)
    {
      fz = (Fx[i + 1][j] - Fx[i][j]);
      fz *= invdx;
      div[i][j] += fz;
    }
    for (i = nx - 4; i < nx; i++)
    {
      fz = (Fx[i][j] - Fx[i - 1][j]);
      fz *= invdx;
      div[i][j] += fz;
    }
  }
}
/**
 *
 * @param float **laplace:
 * @param float **dat:
 * @param float dz:
 * @param float dx:
 * @param int nz:
 * @param int nx:
 * @return {*}
 * @LastEditors: LaurentCheung
 * @Version: 1.0
 */
void laplace2d(float **laplace, float **dat, float dz, float dx, int nz, int nx)
{
  float c0, c1, c2, c3, c4, c5;
  float invdxdx, invdzdz, dxx, dzz;
  int i, j;
  invdxdx = 1.0f / dx / dx;
  invdzdz = 1.0f / dz / dz;
  c0 = -1.46361111111123f;
  c1 = 1.66666666666659f;
  c2 = -0.238095238095227f;
  c3 = 0.0396825396825378f;
  c4 = -0.00496031746031723f;
  c5 = 0.000317460317460303f;
#pragma omp parallel for private(i, j, dxx, dzz)                                                   \
    firstprivate(invdxdx, invdzdz, c0, c1, c2, c3, c4, c5)
  for (i = 5; i < nx - 5; i++)
  {
#pragma ivdep
    for (j = 5; j < nz - 5; j++)
    {
      dxx = c0 * (dat[i][j] + dat[i][j]) + c1 * (dat[i + 1][j] + dat[i - 1][j]) +
            c2 * (dat[i + 2][j] + dat[i - 2][j]) + c3 * (dat[i + 3][j] + dat[i - 3][j]) +
            c4 * (dat[i + 4][j] + dat[i - 4][j]) + c5 * (dat[i + 5][j] + dat[i - 5][j]);
      dzz = c0 * (dat[i][j] + dat[i][j]) + c1 * (dat[i][j + 1] + dat[i][j - 1]) +
            c2 * (dat[i][j + 2] + dat[i][j - 2]) + c3 * (dat[i][j + 3] + dat[i][j - 3]) +
            c4 * (dat[i][j + 4] + dat[i][j - 4]) + c5 * (dat[i][j + 5] + dat[i][j - 5]);
      laplace[i][j] = dxx * invdxdx + dzz * invdzdz;
    }
  }
}
/**
 *
 * @param float **p
 * @param float **f
 * @param float D80
 * @param float D81
 * @param float D82
 * @param float D83
 * @param float D84
 * @param float invdxdx
 * @param float invdzdz
 * @param float beta
 * @param int Nz
 * @param int Nx
 * @param int N
 * @param int maxit
 * @return {*}
 * @LastEditors: LaurentCheung
 * @Version: 1.0
 */
void poissonsolver2d(float **p, float **f, float D80, float D81, float D82, float D83, float D84,
                     float invdxdx, float invdzdz, float beta, int Nz, int Nx, int N, int maxit)
{
  float px, pz;
  int i, j, k;
  for (k = 0; k < maxit; k++)
  {
#pragma omp parallel for private(i, j, px, pz)                                                     \
    firstprivate(invdxdx, invdzdz, D80, D81, D82, D83, D84)
    for (i = N; i <= Nx - N - 1; i++)
    {
#pragma ivdep
      for (j = N; j <= Nz - N - 1; j++)
      {
        pz = D81 * (p[i + 1][j] + p[i - 1][j]) + D82 * (p[i + 2][j] + p[i - 2][j]) +
             D83 * (p[i + 3][j] + p[i - 3][j]) + D84 * (p[i + 4][j] + p[i - 4][j]);
        px = D81 * (p[i][j + 1] + p[i][j - 1]) + D82 * (p[i][j + 2] + p[i][j - 2]) +
             D83 * (p[i][j + 3] + p[i][j - 3]) + D84 * (p[i][j + 4] + p[i][j - 4]);
        p[i][j] =
            beta * (px * invdxdx + pz * invdzdz - f[i][j]) / (-2.0f * D80 * (invdxdx + invdzdz)) +
            (1.0f - beta) * p[i][j];
      }
    }
#pragma omp parallel for private(i, j, px, pz)                                                     \
    firstprivate(invdxdx, invdzdz, D80, D81, D82, D83, D84)
    for (i = Nx - N - 1; i >= N; i--)
    {
#pragma ivdep
      for (j = N; j <= Nz - N - 1; j++)
      {
        pz = D81 * (p[i + 1][j] + p[i - 1][j]) + D82 * (p[i + 2][j] + p[i - 2][j]) +
             D83 * (p[i + 3][j] + p[i - 3][j]) + D84 * (p[i + 4][j] + p[i - 4][j]);
        px = D81 * (p[i][j + 1] + p[i][j - 1]) + D82 * (p[i][j + 2] + p[i][j - 2]) +
             D83 * (p[i][j + 3] + p[i][j - 3]) + D84 * (p[i][j + 4] + p[i][j - 4]);
        p[i][j] =
            beta * (px * invdxdx + pz * invdzdz - f[i][j]) / (-2.0f * D80 * (invdxdx + invdzdz)) +
            (1.0f - beta) * p[i][j];
      }
    }
#pragma omp parallel for private(i, j, px, pz)                                                     \
    firstprivate(invdxdx, invdzdz, D80, D81, D82, D83, D84)
    for (i = N; i <= Nx - N - 1; i++)
    {
#pragma ivdep
      for (j = Nz - N - 1; j >= 0; j--)
      {
        pz = D81 * (p[i + 1][j] + p[i - 1][j]) + D82 * (p[i + 2][j] + p[i - 2][j]) +
             D83 * (p[i + 3][j] + p[i - 3][j]) + D84 * (p[i + 4][j] + p[i - 4][j]);
        px = D81 * (p[i][j + 1] + p[i][j - 1]) + D82 * (p[i][j + 2] + p[i][j - 2]) +
             D83 * (p[i][j + 3] + p[i][j - 3]) + D84 * (p[i][j + 4] + p[i][j - 4]);
        p[i][j] =
            beta * (px * invdxdx + pz * invdzdz - f[i][j]) / (-2.0f * D80 * (invdxdx + invdzdz)) +
            (1.0f - beta) * p[i][j];
      }
    }
#pragma omp parallel for private(i, j, px, pz)                                                     \
    firstprivate(invdxdx, invdzdz, D80, D81, D82, D83, D84)
    for (i = Nx - N - 1; i >= N; i--)
    {
#pragma ivdep
      for (j = Nz - N - 1; j >= 0; j--)
      {
        pz = D81 * (p[i + 1][j] + p[i - 1][j]) + D82 * (p[i + 2][j] + p[i - 2][j]) +
             D83 * (p[i + 3][j] + p[i - 3][j]) + D84 * (p[i + 4][j] + p[i - 4][j]);
        px = D81 * (p[i][j + 1] + p[i][j - 1]) + D82 * (p[i][j + 2] + p[i][j - 2]) +
             D83 * (p[i][j + 3] + p[i][j - 3]) + D84 * (p[i][j + 4] + p[i][j - 4]);
        p[i][j] =
            beta * (px * invdxdx + pz * invdzdz - f[i][j]) / (-2.0f * D80 * (invdxdx + invdzdz)) +
            (1.0f - beta) * p[i][j];
      }
    }
  }
}