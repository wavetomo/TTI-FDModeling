
/************************************/
#include <rockmod.h>
/**/
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/**/
#include <omp.h>
/************************************/

/**
 * convert velocity to slowness
 * @param float *vp: velocity in m/s
 * @param float *slowness: slowness in s/m
 * @param int n: number of points
 * @return float *slowness: slowness in s/m
 * @Author: PingMin Cheung
 * @Date: 2023-06-11 13:47:56
 */
void vel2slowness(float *vp, float *slowness, int n)
{
  int i;
#pragma omp parallel for private(i)
  for (i = 0; i < n; i++)
  {
    slowness[i] = vp[i] == 0.0f ? 0.0f : 1.0 / vp[i];
  }
}

/**
 * convert P-wave velocity to S-wave velocity by Vs/Vp ratios of Possion media
 * @param float *vp: P-wave velocity in m/s
 * @param float *vs: S-wave velocity in m/s
 * @param int n: number of points
 * @return float *vs: S-wave velocity in m/s
 * @Author: PingMin Cheung
 * @Date: 2023-06-11 13:47:56
 */
void vp2vs(float *vp, float *vs, int n)
{
  int i;
#pragma omp parallel for private(i)
  for (i = 0; i < n; i++)
  {
    if (vp[i] <= 1500.0)
    {
      vs[i] = 0.0;
    }
    // else if (vp[i] > 1500.0 && vp[i] < 3500.0)
    // {
    //   vs[i] = 0.6 * vp[i];
    // }
    else
    {
      vs[i] = 0.577350269189626 * vp[i];
    }
  }
}

/**
 * convert P-wave velocity to S-wave velocity by Fliedner's method
 * @param float *vp: P-wave velocity in m/s
 * @param float *vs: S-wave velocity in m/s
 * @param int n: number of points
 * @return float *vs: S-wave velocity in m/s
 * @Author: PingMin Cheung
 * @Date: 2023-06-11 13:47:56
 */
void vp2vsFliedner(float *vp, float *vs, int n)
{
  int i;
#pragma omp parallel for private(i)
  for (i = 0; i < n; i++)
  {
    if (vp[i] <= 1500.0)
    {
      vs[i] = 0.0;
    }
    else if (vp[i] > 1500.0 && vp[i] < 3500.0)
    {
      vs[i] = (vp[i] - 1360.0) / 1.16;
    }
    else
    {
      vs[i] = 0.53 * vp[i];
    }
  }
}

/**
 * convert P-wave velocity to density
 * @param float *vp: P-wave velocity in m/s
 * @param float *rho: density in kg/m^3
 * @param int n: number of points
 * @return float *rho: density in kg/m^3
 * @Author: PingMin Cheung
 * @Date: 2023-06-11 13:47:56
 */
void vp2density(float *vp, float *rho, int n)
{
  int i;
#pragma omp parallel for private(i)
  for (i = 0; i < n; i++)
  {
    if (vp[i] <= 1500.0)
    {
      rho[i] = 1000.0;
    }
    else
    {
      rho[i] = 310.0 * pow(vp[i], 0.25);
    }
  }
}

/**
 * convert P-wave velocity to density by Gardner relation
 * @param float *vp: P-wave velocity in m/s
 * @param float *rho: density in kg/m^3
 * @param int n: number of points
 * @return float *rho: density in kg/m^3
 * @Author: PingMin Cheung
 * @Date: 2023-06-11 13:47:56
 */
void vp2densityGardner(float *vp, float *rho, int n)
{
  int i;
#pragma omp parallel for private(i)
  for (i = 0; i < n; i++)
  {
    if (vp[i] <= 1500.0)
    {
      rho[i] = 1000.0;
    }
    else
    {
      rho[i] = 1741.0 * pow(vp[i] / 1000.0, 0.25); /* Gardner relation */
    }
  }
}

/**
 * convert P-wave velocity to density by Ludwig relation
 * @param float *vp: P-wave velocity in m/s
 * @param float *rho: density in kg/m^3
 * @param int n: number of points
 * @return float *rho: density in kg/m^3
 * @Author: PingMin Cheung
 * @Date: 2023-06-11 13:47:56
 */
void vp2densityLudwig(float *vp, float *rho, int n)
{
  int i;
#pragma omp parallel for private(i)
  for (i = 0; i < n; i++)
  {
    rho[i] = 1000.0 + 0.2 * (vp[i] - 1500.0); /* Ludwig relation */
  }
}

/**
 * convert P-wave velocity to density by Hamilton's method
 * @param float *vp: P-wave velocity in m/s
 * @param float *rho: density in kg/m^3
 * @param int n: number of points
 * @return float *rho: density in kg/m^3
 * @Author: PingMin Cheung
 * @Date: 2023-06-11 13:47:56
 */
void vp2densityHamilton(float *vp, float *rho, int n)
{
  int i;
#pragma omp parallel for private(i)
  for (i = 0; i < n; i++)
  {
    if (vp[i] <= 1500.0)
    {
      rho[i] = 1000.0;
    }
    else if (vp[i] > 1500.0 && vp[i] < 2200.0)
    {
      rho[i] = 2351.0 - 7497.0 * pow(vp[i] / 1000.0, -4.656);
    }
    else
    {
      rho[i] = 1740.0 * pow(vp[i] / 1000.0, 0.25);
    }
  }
}

/**
 * fix the water layer's parameters in the top of 2D model
 * @param float **vp: P-wave velocity in m/s
 * @param float **vs: S-wave velocity in m/s
 * @param float **rho: density in kg/m^3
 * @param float vp_water: P-Wave velocity of water
 * @param float rho_water: density of water
 * @param int *interface: the interface of the bottom of water layer or sea bed
 * @param int nx: the number of samples in the x dimension
 * @param int nz: the number of samples in the z dimension
 * @Author: PingMin Cheung
 * @Date: 2024-10-13 17:06:38
 */
void waterlayer2dfix(float **vp, float **vs, float **rho, float vp_water, float rho_water, int *interface, int nx,
                     int nz)
{
  int ix, iz;
  if (vp != NULL)
  {
#pragma omp parallel for private(ix, iz)
    for (ix = 0; ix < nx; ix++)
    {
#pragma ivdep
      for (iz = 0; iz < interface[ix]; iz++)
      {
        vp[ix][iz] = vp_water;
      }
    }
  }
  if (vs != NULL)
  {
#pragma omp parallel for private(ix, iz)
    for (ix = 0; ix < nx; ix++)
    {
#pragma ivdep
      for (iz = 0; iz < interface[ix]; iz++)
      {
        vs[ix][iz] = 0.0;
      }
    }
  }
  if (rho != NULL)
  {
#pragma omp parallel for private(ix, iz)
    for (ix = 0; ix < nx; ix++)
    {
#pragma ivdep
      for (iz = 0; iz < interface[ix]; iz++)
      {
        rho[ix][iz] = rho_water;
      }
    }
  }
}

/**
 * fix the water layer's parameters in the top of 3D model
 * @param float ***vp: P-wave velocity in m/s
 * @param float ***vs: S-wave velocity in m/s
 * @param float ***rho: density in kg/m^3
 * @param float vp_water: P-Wave velocity of water
 * @param float rho_water: density of water
 * @param int *interface: the interface of the bottom of water layer or sea bed
 * @param int nx: the number of samples in the x dimension
 * @param int ny: the number of samples in the y dimension
 * @param int nz: the number of samples in the z dimension
 * @Author: PingMin Cheung
 * @Date: 2024-10-13 17:06:38
 */
void waterlayer3dfix(float ***vp, float ***vs, float ***rho, float vp_water, float rho_water, int **interface, int nx,
                     int ny, int nz)
{
  int ix, iy, iz;
  if (vp != NULL)
  {
#pragma omp parallel for private(ix, iy, iz)
    for (iy = 0; iy < ny; iy++)
    {
      for (ix = 0; ix < nx; ix++)
      {
#pragma ivdep
        for (iz = 0; iz < interface[iy][ix]; iz++)
        {
          vp[iy][ix][iz] = vp_water;
        }
      }
    }
  }
  if (vs != NULL)
  {
#pragma omp parallel for private(ix, iy, iz)
    for (iy = 0; iy < ny; iy++)
    {
      for (ix = 0; ix < nx; ix++)
      {
#pragma ivdep
        for (iz = 0; iz < interface[iy][ix]; iz++)
        {
          vs[iy][ix][iz] = 0.0;
        }
      }
    }
  }
  if (rho != NULL)
  {
#pragma omp parallel for private(ix, iy, iz)
    for (iy = 0; iy < ny; iy++)
    {
      for (ix = 0; ix < nx; ix++)
      {
#pragma ivdep
        for (iz = 0; iz < interface[iy][ix]; iz++)
        {
          rho[iy][ix][iz] = rho_water;
        }
      }
    }
  }
}