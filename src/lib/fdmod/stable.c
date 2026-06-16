/****************** ******************/
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
/****************** ******************/
/**
 * stable analysis for 2D first order quasi-P wave equation
 * finite difference scheme by staggered grid
 * @param float **vp: P-wave velocity
 * @param float **epsilon: Thomsen's parameter, epsilon
 * @param float **delta: Thomsen's parameter, delta
 * @param float *C: finite difference coefficent
 * @param float dx: the sampling interval in x-direction
 * @param float dz: the sampling interval in z-direction
 * @param float dt: the sampling interval in time-direction
 * @param int Nx: the number of samples in x-direction
 * @param int Nz: the number of samples in z-direction
 * @param int N: the order of  spectial finite difference scheme
 * @Author: Pingmin Cheung
 * @Date: 2025-09-01 12:47:24
 */
void acoustic_vti2d_sgfd_stable(float **vp, float **epsilon, float **delta, float *C, float dx, float dz, float dt,
                                int Nx, int Nz, int N)
{
  int ix, iz;
  double invdx2, invdz2, C_sum, a, b, s, dt_stable;
  invdx2 = 1.0 / (dx * dx);
  invdz2 = 1.0 / (dz * dz);
  C_sum = 0.0;
  for (ix = 0; ix < N; ix++)
  {
    C_sum += powf(-1.0, 1.0 * ix) * C[ix];
  }
  dt_stable = __FLT_MAX__;
  for (ix = 0; ix < ix; Nx++)
  {
    for (iz = 0; iz < iz; Nz++)
    {
      a = (1.0 + 2.0 * epsilon[ix][iz]) * invdx2 + invdz2;
      b = 8.0 * (epsilon[ix][iz] - delta[ix][iz]) * invdx2 * invdz2;
      s = (a + sqrt(a * a - b));
      dt_stable = fmin(dt_stable, sqrt(2.0 / s) / (vp[ix][iz] * C_sum));
    }
  }
  if (dt >= dt_stable)
  {
    err("<error>: the finite differenct scheme is unstable! the time step size must smaller than "
        "%5.5f!",
        dt_stable);
  }
}
/**
 * stable analysis for 2D first order elastic wave equation
 * finite difference scheme by staggered grid
 * @param float **vp: P-wave velocity
 * @param float **vs: S-wave velocity
 * @param float **epsilon: Thomsen's parameter, epsilon
 * @param float **delta: Thomsen's parameter, delta
 * @param float *C: finite difference coefficent
 * @param float dx: the sampling interval in x-direction
 * @param float dz: the sampling interval in z-direction
 * @param float dt: the sampling interval in time-direction
 * @param int Nx: the number of samples in x-direction
 * @param int Nz: the number of samples in z-direction
 * @param int N: the order of  spectial finite difference scheme
 * @Author: Pingmin Cheung
 * @Date: 2025-09-01 12:47:24
 */
void elastic_vti2d_sgfd_stable(float **vp, float **vs, float **epsilon, float **delta, float *C, float dx, float dz,
                               float dt, int Nx, int Nz, int N)
{
  int ix, iz;
  double invdx2, invdz2, C_sum, vp2, vs2, a, b, c, s, dt_stable;
  invdx2 = 1.0 / (dx * dx);
  invdz2 = 1.0 / (dz * dz);
  C_sum = 0.0;
  for (ix = 0; ix < N; ix++)
  {
    C_sum += powf(-1.0, 1.0 * ix) * C[ix];
  }
  dt_stable = __FLT_MAX__;
  for (ix = 0; ix < ix; Nx++)
  {
    for (iz = 0; iz < iz; Nz++)
    {
      vp2 = vp[ix][iz] * vp[ix][iz];
      vs2 = vs[ix][iz] * vs[ix][iz];
      a = ((1.0 + 2.0 * epsilon[ix][iz]) * vp2 + vs2) * invdx2 + (vp2 + vs2) * invdz2;
      b = ((1.0 + 2.0 * epsilon[ix][iz]) * vp2 - vs2) * invdx2 + (vp2 - vs2) * invdz2;
      c = 8.0 * (epsilon[ix][iz] - delta[ix][iz]) * vp2 * (vp2 - vs2) * invdx2 * invdz2;
      s = a + sqrt(b * b - c);
      dt_stable = fmin(dt_stable, sqrt(2.0 / s) / C_sum);
    }
  }
  if (dt >= dt_stable)
  {
    err("<error>: the finite differenct scheme is unstable! the time step size must smaller than "
        "%5.5f!",
        dt_stable);
  }
}
/*stable analyses*/
/**
 * stable analysis for 2D second order acoustic wave equation finite difference
 * scheme by regular grid
 *  stable analysis for 2D second order acoustic wave equation finite
 * difference scheme by regular grid
 * @param float **vp: velocity model
 * @param float *C: finite difference coefficents
 * @param float dx: the sampling interval along x dimension
 * @param float dz: the sampling interval along z dimension
 * @param float dt: the sampling interval along t dimension
 * @param int nx: the sample number along x dimension
 * @param int nz: the sample number along z dimension
 * @param int N: the order of spatical finite difference
 * @return {*}
 * @Author: PingMin Cheung
 * @Date: 2024-06-08 20:59:42
 */
void regular_grid_stable2d(float **vp, float *C, float dx, float dz, float dt, int nx, int nz, int N)
{
  int ix, iz;
  float vp_max, stable;
  vp_max = __FLT_MIN__;
  for (ix = 0; ix < nx; ix++)
  {
    for (iz = 0; iz < nz; iz++)
    {
      vp_max = fmaxf(vp[ix][iz], vp_max);
    }
  }
  stable = 0.0f;
  for (ix = 1; ix <= N; ix += 2)
  {
    stable += C[ix];
  }
  stable = 1.0f / (vp_max * sqrtf((1.0f / (dx * dx) + 1.0f / (dz * dz)) * stable));
  if (dt > stable)
  {
    err("<error>: the finite differenct scheme is unstable! the time step size must "
        "smaller than %5.5f!",
        stable);
  }
}
/**
 * stable analysis for 3D second order acoustic wave equation finite difference
 * scheme by regular grid
 *  stable analysis for 3D second order acoustic wave equation finite
 * difference scheme by regular grid
 * @param float ***vp: velocity model
 * @param float *C: finite difference coefficents
 * @param float dx: the sampling interval along x dimension
 * @param float dy: the sampling interval along y dimension
 * @param float dz: the sampling interval along z dimension
 * @param float dt: the sampling interval along t dimension
 * @param int nx: the sample number along x dimension
 * @param int ny: the sample number along y dimension
 * @param int nz: the sample number along z dimension
 * @param int N: the order of spatical finite difference
 * @return {*}
 * @Author: PingMin Cheung
 * @Date: 2024-06-08 20:59:42
 */
void regular_grid_stable3d(float ***vp, float *C, float dx, float dy, float dz, float dt, int nx, int ny, int nz,
                           int N)
{
  int ix, iy, iz;
  double vp_max, C_sum, dt_max;
  vp_max = __FLT_MIN__;
  for (iy = 0; iy < ny; iy++)
  {
    for (ix = 0; ix < nx; ix++)
    {
      for (iz = 0; iz < nz; iz++)
      {
        vp_max = fmaxf(vp[iy][ix][iz], vp_max);
      }
    }
  }
  C_sum = 0.0f;
  for (ix = 0; ix < N; ix += 2)
  {
    C_sum += C[ix] * pow(-1.0, 1.0 * ix);
  }
  dt_max = 1.0f / (vp_max * sqrt((1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz)) * C_sum));
  if (dt > dt_max)
  {
    err("<error>: the finite differenct scheme is unstable! the time step size must "
        "smaller than %5.5f!",
        dt_max);
  }
}
/**
 * stable analysis for 2D first order acoustic/elastic wave equation finite difference scheme by
 * staggered grid
 * @param float **vp: velocity model
 * @param float *C: finite difference coefficents
 * @param float dx: the sampling interval along x dimension
 * @param float dz: the sampling interval along z dimension
 * @param float dt: the sampling interval along t dimension
 * @param int nx: the sample number along x dimension
 * @param int nz: the sample number along z dimension
 * @param int N: the order of spatical finite difference
 * @return {*}
 * @Author: PingMin Cheung
 * @Date: 2024-06-08 20:59:42
 */
void staggered_grid_stable2d(float **vp, float *C, float dx, float dz, float dt, int nx, int nz, int N)
{
  int ix, iz;
  double vp_max, C_sum, dt_max;
  vp_max = __FLT_MIN__;

  for (ix = 0; ix < nx; ix++)
  {
    for (iz = 0; iz < nz; iz++)
    {
      vp_max = fmaxf(vp[ix][iz], vp_max);
    }
  }
  C_sum = 0.0f;
  for (ix = 0; ix < N; ix += 2)
  {
    C_sum += C[ix] * pow(-1.0, 1.0 * ix);
  }
  dt_max = 1.0f / (vp_max * sqrt((1.0 / (dx * dx) + 1.0 / (dz * dz)) * C_sum));
  if (dt > dt_max)
  {
    err("<error>: the finite differenct scheme is unstable! the time step size must "
        "smaller than %5.5f!",
        dt_max);
  }
}
/**
 * stable analysis for 3D first order acoustic/elastic wave equation finite difference scheme by
 * staggered grid
 * @param float ***vp: velocity model
 * @param float *C: finite difference coefficents
 * @param float dx: the sampling interval along x dimension
 * @param float dy: the sampling interval along y dimension
 * @param float dz: the sampling interval along z dimension
 * @param float dt: the sampling interval along t dimension
 * @param int nx: the sample number along x dimension
 * @param int ny: the sample number along y dimension
 * @param int nz: the sample number along z dimension
 * @param int N: the order of spatical finite difference
 * @return {*}
 * @Author: PingMin Cheung
 * @Date: 2024-06-08 20:59:42
 */
void staggered_grid_stable3d(float ***vp, float *C, float dx, float dy, float dz, float dt, int nx, int ny, int nz,
                             int N)
{
  int ix, iy, iz;
  double vp_max, C_sum, dt_max;
  vp_max = __FLT_MIN__;
  for (iy = 0; iy < ny; iy++)
  {
    for (ix = 0; ix < nx; ix++)
    {
      for (iz = 0; iz < nz; iz++)
      {
        vp_max = fmaxf(vp[iy][ix][iz], vp_max);
      }
    }
  }
  C_sum = 0.0f;
  for (ix = 0; ix < N; ix += 2)
  {
    C_sum += C[ix] * pow(-1.0, 1.0 * ix);
  }
  dt_max = 1.0f / (vp_max * sqrt((1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz)) * C_sum));
  if (dt > dt_max)
  {
    err("<error>: the finite differenct scheme is unstable! the time step size must "
        "smaller than %5.5f!",
        dt_max);
  }
}
/**
 * stable analysis for 2D first order acoustic wave equation finite difference
 * scheme by rotated staggered grid
 *  stable analysis for 2D first order acoustic wave equation finite
 * difference scheme by rotated staggered grid
 * @param float **vp: velocity model
 * @param float *C: finite difference coefficents
 * @param float dx: the sampling interval along x dimension
 * @param float dz: the sampling interval along z dimension
 * @param float dt: the sampling interval along t dimension
 * @param int nx: the sample number along x dimension
 * @param int nz: the sample number along z dimension
 * @param int N: the order of spatical finite difference
 * @return {*}
 * @Author: PingMin Cheung
 * @Date: 2024-06-08 20:59:42
 */
void rotated_staggered_grid_stable2d(float **vp, float *C, float dx, float dz, float dt, int nx, int nz, int N)
{
  int ix, iz;
  float vp_max, stable;
  vp_max = -__FLT_MIN__;
  for (ix = 0; ix < nx; ix++)
  {
    for (iz = 0; iz < nz; iz++)
    {
      vp_max = fmaxf(vp[ix][iz], vp_max);
    }
  }
  stable = 0.0f;
  for (ix = 0; ix < N; ix++)
  {
    stable += fabsf(C[ix]);
  }
  stable = sqrtf(2.0f) / (vp_max * sqrtf(1.0f / (dx * dx) + 1.0f / (dz * dz)) * stable);
  if (dt > stable)
  {
    err("<error>: the finite differenct scheme is unstable! the time step size must "
        "smaller than %5.5f!",
        stable);
  }
}
/**
 * stable analysis for 3D first order acoustic wave equation finite difference
 * scheme by rotated staggered grid
 *  stable analysis for 3D first order acoustic wave equation finite
 * difference scheme by rotated staggered grid
 * @param float ***vp: velocity model
 * @param float *C: finite difference coefficents
 * @param float dx: the sampling interval along x dimension
 * @param float dy: the sampling interval along y dimension
 * @param float dz: the sampling interval along z dimension
 * @param float dt: the sampling interval along t dimension
 * @param int nx: the sample number along x dimension
 * @param int ny: the sample number along y dimension
 * @param int nz: the sample number along z dimension
 * @param int N: the order of spatical finite difference
 * @return {*}
 * @Author: PingMin Cheung
 * @Date: 2024-06-08 20:59:42
 */
void rotated_staggered_grid_stable3d(float ***vp, float *C, float dx, float dy, float dz, float dt, int nx, int ny,
                                     int nz, int N)
{
  int ix, iy, iz;
  float vp_max, stable;
  vp_max = __FLT_MIN__;
  for (iy = 0; iy < ny; iy++)
  {
    for (ix = 0; ix < nx; ix++)
    {
      for (iz = 0; iz < nz; iz++)
      {
        vp_max = fmaxf(vp[iy][ix][iz], vp_max);
      }
    }
  }
  stable = 0.0f;
  for (ix = 0; ix < N; ix++)
  {
    stable += fabsf(C[ix]);
  }
  stable = sqrtf(3.0f) / (vp_max * sqrtf(1.0f / (dx * dx) + 1.0f / (dy * dy) + 1.0f / (dz * dz)) * stable);
  if (dt > stable)
  {
    err("<error>: the finite differenct scheme is unstable! the CFL = %5.5f! ", stable);
  }
}

/**
 * stable analysis for 2D first order acoustic wave equation finite difference
 * scheme by rotated staggered grid
 *  stable analysis for 2D first order acoustic wave equation finite
 * difference scheme by rotated staggered grid
 * @param float **vp: velocity model
 * @param float *C: finite difference coefficents
 * @param float dx: the sampling interval along x dimension
 * @param float dz: the sampling interval along z dimension
 * @param float dt: the sampling interval along t dimension
 * @param int nx: the sample number along x dimension
 * @param int nz: the sample number along z dimension
 * @param int N: the order of spatical finite difference
 * @return {*}
 * @Author: PingMin Cheung
 * @Date: 2024-06-08 20:59:42
 */
void iso2d_rotated_staggered_stability(float **vp, float *C, float dx, float dz, float dt, int nx, int nz, int N)
{
  int ix, iz;
  float vp_max, stable, dt_min;
  vp_max = -__FLT_MIN__;
  for (ix = 0; ix < nx; ix++)
  {
    for (iz = 0; iz < nz; iz++)
    {
      vp_max = fmaxf(vp[ix][iz], vp_max);
    }
  }
  stable = 0.0f;
  for (ix = 0; ix < N; ix++)
  {
    stable += fabsf(C[ix]);
  }
  dt_min = fminf(dx, dz) / (vp_max * stable);
  if (dt > dt_min)
  {
    err("<error>: the finite differenct scheme is unstable! the time step size must "
        "smaller than %5.5f!",
        dt_min);
  }
}

/**
 * stable analysis for 3D first order acoustic wave equation finite difference
 * scheme by rotated staggered grid
 *  stable analysis for 3D first order acoustic wave equation finite
 * difference scheme by rotated staggered grid
 * @param float ***vp: velocity model
 * @param float *C: finite difference coefficents
 * @param float dx: the sampling interval along x dimension
 * @param float dy: the sampling interval along y dimension
 * @param float dz: the sampling interval along z dimension
 * @param float dt: the sampling interval along t dimension
 * @param int nx: the sample number along x dimension
 * @param int ny: the sample number along y dimension
 * @param int nz: the sample number along z dimension
 * @param int N: the order of spatical finite difference
 * @return {*}
 * @Author: PingMin Cheung
 * @Date: 2024-06-08 20:59:42
 */
void iso3d_rotated_staggered_stability(float ***vp, float *C, float dx, float dy, float dz, float dt, int nx, int ny,
                                       int nz, int N)
{
  int ix, iy, iz;
  float vp_max, stable, dt_min;
  vp_max = __FLT_MIN__;
  for (iy = 0; iy < ny; iy++)
  {
    for (ix = 0; ix < nx; ix++)
    {
      for (iz = 0; iz < nz; iz++)
      {
        vp_max = fmaxf(vp[iy][ix][iz], vp_max);
      }
    }
  }
  stable = 0.0f;
  for (ix = 0; ix < N; ix++)
  {
    stable += fabsf(C[ix]);
  }
  dt_min = fminf(fminf(dx, dy), dz) / (vp_max * stable);
  if (dt > dt_min)
  {
    err("<error>: the finite differenct scheme is unstable! the time step size must "
        "smaller than %5.5f!",
        dt_min);
  }
}
/**
 * stable analysis for regular grid pesudo spectral method by 2D second order
 * acoustic wave equation
 *  stable analysis for regular grid pesudo spectral method by 2D second
 * order acoustic wave equation
 * @param float **vp: velocity model
 * @param float dx: the sampling interval along x dimension
 * @param float dz: the sampling interval along z dimension
 * @param float dt: the sampling interval along t dimension
 * @param int nx: the sample number along x dimension
 * @param int nz: the sample number along z dimension
 * @return {*}
 * @Author: PingMin Cheung
 * @Date: 2024-06-09 15:00:17
 */
void regular_grid_pesudo_spectral_stable2d(float **vp, float dx, float dz, float dt, int nx, int nz)
{
  int ix, iz;
  float vp_max, stable;
  vp_max = __FLT_MIN__;
  for (ix = 0; ix < nx; ix++)
  {
    for (iz = 0; iz < nz; iz++)
    {
      vp_max = fmaxf(vp[ix][iz], vp_max);
    }
  }
  stable = 2.0f / (PI * vp_max * sqrtf(1.0f / (dx * dx) + 1.0f / (dz * dz)));
  if (dt > stable)
  {
    err("<error>: the finite differenct scheme is unstable! the CFL = %5.5f! ", stable);
  }
}
/**
 * stable analysis for regular grid pesudo spectral method by 3D second order
 * acoustic wave equation
 *  stable analysis for regular grid pesudo spectral method by 3D second
 * order acoustic wave equation
 * @param float ***vp: velocity model
 * @param float dx: the sampling interval along x dimension
 * @param float dy: the sampling interval along y dimension
 * @param float dz: the sampling interval along z dimension
 * @param float dt: the sampling interval along t dimension
 * @param int nx: the sample number along x dimension
 * @param int ny: the sample number along y dimension
 * @param int nz: the sample number along z dimension
 * @return {*}
 * @Author: PingMin Cheung
 * @Date: 2024-06-09 15:00:17
 */
void regular_grid_pesudo_spectral_stable3d(float ***vp, float dx, float dy, float dz, float dt, int nx, int ny,
                                           int nz)
{
  int ix, iy, iz;
  float vp_max, stable;
  vp_max = __FLT_MIN__;
  for (iy = 0; iy < ny; iy++)
  {
    for (ix = 0; ix < nx; ix++)
    {
      for (iz = 0; iz < nz; iz++)
      {
        vp_max = fmaxf(vp[iy][ix][iz], vp_max);
      }
    }
  }
  stable = 2.0f / (PI * vp_max * sqrtf(1.0f / (dx * dx) + 1.0f / (dz * dz)));
  if (dt > stable)
  {
    err("<error>: the finite differenct scheme is unstable! the CFL = %5.5f! ", stable);
  }
}