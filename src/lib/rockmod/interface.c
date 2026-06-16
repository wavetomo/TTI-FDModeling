
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
 * gridding of 2D interfaces
 * @param float *interface: coordinates of interfaces at z-axis
 * @param float fz: first coordinates at z-axis
 * @param float dz: sample interval at z-axis
 * @param int *fixpos: interface position after gridding
 * @param int nx: the number of samples at x-axis
 * @param int nz: the number of samples at z-axis
 * @return int *fixpos: interface position after gridding
 * @Author: PingMin Cheung
 * @Date: 2024-10-06 15:20:17
 */
void interface2d_gridding(float *interface, float fz, float dz, int *fixpos, int nx, int nz)
{
  int ix;
#pragma omp parallel for private(ix)
  for (ix = 0; ix < nx; ix++)
  {
    interface[ix] = (interface[ix] - fz) / dz;
    if (interface[ix] - (int)interface[ix] == 0)
    {
      fixpos[ix] = fminf(fmaxf((int)interface[ix], 0), nz - 1);
    }
    else
    {
      fixpos[ix] = fminf(fmaxf((int)interface[ix] + 1, 0), nz - 1);
    }
  }
}
/**
 * gridding of 3D interfaces
 * @param float **interface: coordinates of interfaces at z-axis
 * @param float fz: first coordinates at z-axis
 * @param float dz: sample interval at z-axis
 * @param int **fixpos: interface position after gridding
 * @param int nx: the number of samples at x-axis
 * @param int ny: the number of samples at y-axis
 * @param int nz: the number of samples at z-axis
 * @return int *fixpos: interface position after gridding
 * @Author: PingMin Cheung
 * @Date: 2024-10-06 15:20:17
 */
void interface3d_gridding(float **interface, float fz, float dz, int **fixpos, int nx, int ny, int nz)
{
  int ix, iy;
#pragma omp parallel for private(ix)
  for (iy = 0; iy < ny; iy++)
  {
#pragma ivdep
    for (ix = 0; ix < nx; ix++)
    {
      interface[iy][ix] = (interface[iy][ix] - fz) / dz;
      if (interface[iy][ix] - (int)interface[iy][ix] == 0)
      {
        fixpos[iy][ix] = fminf(fmaxf((int)interface[iy][ix], 0), nz - 1);
      }
      else
      {
        fixpos[iy][ix] = fminf(fmaxf((int)interface[iy][ix] + 1, 0), nz - 1);
      }
    }
  }
}