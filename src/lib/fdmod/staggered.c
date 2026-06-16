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
 * calculate the model's parameters on the center of 2d staggered grid
 * @param float **mod_in : input model parameters
 * @param float **mod_out : output model parameters
 * @param int nz: the number of grid points in the z direction
 * @param int nx: the number of grid points in the x direction
 * @param int type: type of grid
 *                  1 for inverse of rho in x direction,
 *                  2 for inverse of rho in z direction,
 *                  3 for miu in the center of staggered grid
 * @return float **mod_out : output model parameters
 * @Author: PingMin Cheung
 * @Date: 2024-07-31 15:11:49
 */
void staggered_grid2d_parameter(float **mod_in, float **mod_out, int nz, int nx, int type)
{
  int ix, iz;
  float **temp = NULL;
  if (mod_in == mod_out)
  {
    temp = alloc2float(nz, nx);
    memcpy(temp[0], mod_in[0], nz * nx * FSIZE);
  }
  else
  {
    temp = mod_in;
  }
  if (type == 1)
  {
#pragma omp parallel for private(ix, iz)
    for (ix = 0; ix < nx; ix++)
    {
#pragma ivdep
      for (iz = 0; iz < nz; iz++)
      {
        if (ix < nx - 1)
        {
          mod_out[ix][iz] = 2.0f / (temp[ix][iz] + temp[ix + 1][iz]);
        }
        else
        {
          mod_out[ix][iz] = 1.0f / temp[ix][iz];
        }
      }
    }
  }
  if (type == 2)
  {
#pragma omp parallel for private(ix, iz)
    for (ix = 0; ix < nx; ix++)
    {
#pragma ivdep
      for (iz = 0; iz < nz; iz++)
      {
        if (iz < nz - 1)
        {
          mod_out[ix][iz] = 2.0f / (temp[ix][iz] + temp[ix][iz + 1]);
        }
        else
        {
          mod_out[ix][iz] = 1.0f / temp[ix][iz];
        }
      }
    }
  }
  if (type == 3)
  {
#pragma omp parallel for private(ix, iz)
    for (ix = 0; ix < nx; ix++)
    {
#pragma ivdep
      for (iz = 0; iz < nz; iz++)
      {
        if (ix < nx - 1 && iz < nz - 1)
        {
          mod_out[ix][iz] = 4.0f / (1.0f / temp[ix][iz] + 1.0f / temp[ix + 1][iz] + 1.0f / temp[ix][iz + 1] +
                                    1.0f / temp[ix + 1][iz + 1]);
          if (isnan(mod_out[ix][iz]) || isinf(mod_out[ix][iz]))
          {
            mod_out[ix][iz] = 0.25 * (temp[ix][iz] + temp[ix + 1][iz] + temp[ix][iz + 1] + temp[ix + 1][iz + 1]);
          }
        }
        else
        {
          mod_out[ix][iz] = temp[ix][iz];
        }
      }
    }
  }
  if (mod_in == mod_out)
  {
    free2float(temp);
  }
}
/**
 * calculate the model's parameters on the center of 3d staggered grid
 * @param float ***mod_in : input model parameters
 * @param float ***mod_out : output model parameters
 * @param int nz: the number of grid points in the z direction
 * @param int nx: the number of grid points in the x direction
 * @param int ny: the number of grid points in the y direction
 * @return float ***mod_out : output model parameters
 * @param int type: type of grid
 *                  1 for inverse of rho in x direction,
 *                  2 for inverse of rho in y direction,
 *                  3 for inverse of rho in z direction,
 *                  4 for miu in the center of staggered grid at XOY plane
 *                  5 for miu in the center of staggered grid at XOZ plane
 *                  6 for miu in the center of staggered grid at YOZ plane
 * @Author: PingMin Cheung
 * @Date: 2024-07-31 15:11:49
 */
void staggered_grid3d_parameter(float ***mod_in, float ***mod_out, int nz, int nx, int ny, int type)
{
  int ix, iy, iz;
  float ***temp = NULL;
  if (mod_in == mod_out)
  {
    temp = alloc3float(nz, nx, ny);
    memcpy(temp[0][0], mod_in[0][0], nz * nx * ny * FSIZE);
  }
  else
  {
    temp = mod_in;
  }
  if (type == 1)
  {
#pragma omp parallel for private(ix, iy, iz)
    for (ix = 0; ix < nx; ix++)
    {
      for (iy = 0; iy < ny; iy++)
      {
#pragma ivdep
        for (iz = 0; iz < nz; iz++)
        {
          if (ix < nx - 1)
          {
            mod_out[iy][ix][iz] = 2.0f / (temp[iy][ix][iz] + temp[iy][ix + 1][iz]);
          }
          else
          {
            mod_out[iy][ix][iz] = 1.0f / temp[iy][ix][iz];
          }
        }
      }
    }
  }
  if (type == 2)
  {
#pragma omp parallel for private(ix, iy, iz)
    for (ix = 0; ix < nx; ix++)
    {
      for (iy = 0; iy < ny; iy++)
      {
#pragma ivdep
        for (iz = 0; iz < nz; iz++)
        {
          if (iy < ny - 1)
          {
            mod_out[iy][ix][iz] = 2.0f / (temp[iy][ix][iz] + temp[iy + 1][ix][iz]);
          }
          else
          {
            mod_out[iy][ix][iz] = 1.0f / temp[iy][ix][iz];
          }
        }
      }
    }
  }
  if (type == 3)
  {
#pragma omp parallel for private(ix, iz)
    for (ix = 0; ix < nx; ix++)
    {
#pragma ivdep
      for (iz = 0; iz < nz; iz++)
      {
        if (iz < nz - 1)
        {
          mod_out[iy][ix][iz] = 2.0f / (temp[iy][ix][iz] + temp[iy][ix][iz + 1]);
        }
        else
        {
          mod_out[iy][ix][iz] = 1.0f / temp[iy][ix][iz];
        }
      }
    }
  }
  if (type == 4)
  {
#pragma omp parallel for private(ix, iy, iz)
    for (ix = 0; ix < nx; ix++)
    {
      for (iy = 0; iy < ny; iy++)
      {
#pragma ivdep
        for (iz = 0; iz < nz; iz++)
        {
          if (ix < nx - 1 && iy < ny - 1)
          {
            mod_out[iy][ix][iz] = 4.0f / (1.0f / temp[iy][ix][iz] + 1.0f / temp[iy][ix + 1][iz] +
                                          1.0f / temp[iy + 1][ix][iz] + 1.0f / temp[iy + 1][ix + 1][iz]);
            if (isnan(mod_out[iy][ix][iz]) || isinf(mod_out[iy][ix][iz]))
            {
              mod_out[iy][ix][iz] =
                  0.25 * (temp[iy][ix][iz] + temp[iy][ix + 1][iz] + temp[iy + 1][ix][iz] + temp[iy + 1][ix + 1][iz]);
            }
          }
          else
          {
            mod_out[iy][ix][iz] = temp[iy][ix][iz];
          }
        }
      }
    }
  }
  if (type == 5)
  {
#pragma omp parallel for private(ix, iy, iz)
    for (ix = 0; ix < nx; ix++)
    {
      for (iy = 0; iy < ny; iy++)
      {
#pragma ivdep
        for (iz = 0; iz < nz; iz++)
        {
          if (ix < nx - 1 && iz < nz - 1)
          {
            mod_out[iy][ix][iz] = 4.0f / (1.0f / temp[iy][ix][iz] + 1.0f / temp[iy][ix + 1][iz] +
                                          1.0f / temp[iy][ix][iz + 1] + 1.0f / temp[iy][ix + 1][iz + 1]);
            if (isnan(mod_out[iy][ix][iz]) || isinf(mod_out[iy][ix][iz]))
            {
              mod_out[iy][ix][iz] =
                  0.25 * (temp[iy][ix][iz] + temp[iy][ix + 1][iz] + temp[iy][ix][iz + 1] + temp[iy][ix + 1][iz + 1]);
            }
          }
          else
          {
            mod_out[iy][ix][iz] = temp[iy][ix][iz];
          }
        }
      }
    }
  }
  if (type == 6)
  {
#pragma omp parallel for private(ix, iy, iz)
    for (ix = 0; ix < nx; ix++)
    {
      for (iy = 0; iy < ny; iy++)
      {
#pragma ivdep
        for (iz = 0; iz < nz; iz++)
        {
          if (iy < ny - 1 && iz < nz - 1)
          {
            mod_out[iy][ix][iz] = 4.0f / (1.0f / temp[iy][ix][iz] + 1.0f / temp[iy + 1][ix][iz] +
                                          1.0f / temp[iy][ix][iz + 1] + 1.0f / temp[iy + 1][ix][iz + 1]);
            if (isnan(mod_out[iy][ix][iz]) || isinf(mod_out[iy][ix][iz]))
            {
              mod_out[iy][ix][iz] =
                  0.25 * (temp[iy][ix][iz] + temp[iy + 1][ix][iz] + temp[iy][ix][iz + 1] + temp[iy + 1][ix][iz + 1]);
            }
          }
          else
          {
            mod_out[iy][ix][iz] = temp[iy][ix][iz];
          }
        }
      }
    }
  }
  if (mod_in == mod_out)
  {
    free3float(temp);
  }
}
/**
 * calculate the model's parameters on the center of 2d rotated staggered grid
 * @param float **mod_in : input model parameters
 * @param float **mod_out : output model parameters
 * @param int nz: the number of grid points in the z direction
 * @param int nx: the number of grid points in the x direction
 * @return float **mod_out : output model parameters
 * @Author: PingMin Cheung
 * @Date: 2024-07-31 15:11:49
 */
void rotated_staggered_grid2d_parameter(float **mod_in, float **mod_out, int nz, int nx)
{
  int ix, iz;
  float **temp = NULL;
  if (mod_in == mod_out)
  {
    temp = alloc2float(nz, nx);
    memcpy(temp[0], mod_in[0], nz * nx * FSIZE);
  }
  else
  {
    temp = mod_in;
  }
#pragma omp parallel for private(ix, iz)
  for (ix = 0; ix < nx; ix++)
  {
#pragma ivdep
    for (iz = 0; iz < nz; iz++)
    {
      if (ix < nx - 1 && iz < nz - 1)
      {
        mod_out[ix][iz] = 0.25 * (temp[ix][iz] + temp[ix + 1][iz] + temp[ix][iz + 1] + temp[ix + 1][iz + 1]);
      }
      else
      {
        mod_out[ix][iz] = temp[ix][iz];
      }
    }
  }
  if (mod_in == mod_out)
  {
    free2float(temp);
  }
}
/**
 * calculate the model's parameters on the center of 3d rotated staggered grid
 * @param float ***mod_in : input model parameters
 * @param float ***mod_out : output model parameters
 * @param int nz: the number of grid points in the z direction
 * @param int nx: the number of grid points in the x direction
 * @param int ny: the number of grid points in the y direction
 * @return float ***mod_out : output model parameters
 * @Author: PingMin Cheung
 * @Date: 2024-07-31 15:11:49
 */
void rotated_staggered_grid3d_parameter(float ***mod_in, float ***mod_out, int nz, int nx, int ny)
{
  int ix, iy, iz;
  float ***temp = NULL;
  if (mod_in == mod_out)
  {
    temp = alloc3float(nz, nx, ny);
    memcpy(temp[0][0], mod_in[0][0], nz * nx * ny * FSIZE);
  }
  else
  {
    temp = mod_in;
  }
#pragma omp parallel for private(ix, iy, iz)
  for (ix = 0; ix < nx; ix++)
  {
    for (iy = 0; iy < ny; iy++)
    {
#pragma ivdep
      for (iz = 0; iz < nz; iz++)
      {
        if (ix < nx - 1 && iy < ny - 1 && iz < nz - 1)
        {
          mod_out[ix][iy][iz] = 0.125 * (temp[ix][iy][iz] + temp[ix + 1][iy][iz] + temp[ix][iy + 1][iz] +
                                         temp[ix + 1][iy + 1][iz] + temp[ix][iy][iz + 1] + temp[ix + 1][iy][iz + 1] +
                                         temp[ix][iy + 1][iz + 1] + temp[ix + 1][iy + 1][iz + 1]);
        }
        else
        {
          mod_out[ix][iy][iz] = temp[ix][iy][iz];
        }
      }
    }
  }
  if (mod_in == mod_out)
  {
    free3float(temp);
  }
}