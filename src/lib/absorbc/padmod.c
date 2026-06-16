/*
 * @description:
 * @Author: PingMin Cheung
 * @Date: 2024-03-17 14:54:57
 * @LastEditTime: 2024-05-19 21:04:06
 * @FilePath: /src/fdmod/lib/absorbc/padmod.c
 *
 * Copyright (c) 2024 by WaveTomo, All Rights Reserved.
 */
/************************************/
#include <absorbc.h>
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
/************************************/
/**
 * padding 2D model
 * @description: padding 2D model
 * @param {float**} mod: model array
 * @param int Nx: the sample number along x dimension
 * @param int Nz: the sample number along z dimension
 * @param int nlayer: extend layers
 * @LastEditors: LaurentCheung
 * @Version: 1.0
 */
void padModel2d(float **mod, int Nx, int Nz, int nlayer)
{
  int ix, iz, ilayer;
  /* padding top and bottom of model */
#pragma omp parallel for private(ix, ilayer) firstprivate(nlayer)
  for (ix = nlayer; ix < Nx - nlayer; ix++)
  {
#pragma ivdep
    for (ilayer = 1; ilayer <= nlayer; ilayer++)
    {
      mod[ix][nlayer - ilayer] = mod[ix][nlayer];
      mod[ix][Nz - nlayer - 1 + ilayer] = mod[ix][Nz - nlayer - 1];
    }
  }
  /* padding left and right of model */
#pragma omp parallel for private(iz, ilayer) firstprivate(nlayer)
  for (ilayer = 1; ilayer <= nlayer; ilayer++)
  {
#pragma ivdep
    for (iz = 0; iz < Nz; iz++)
    {
      mod[nlayer - ilayer][iz] = mod[nlayer][iz];
      mod[Nx - nlayer - 1 + ilayer][iz] = mod[Nx - nlayer - 1][iz];
    }
  }
}
/**
 * extend 2D model
 * @description: extend 2D model
 * @param {float**} mod: model array
 * @param int Nx: the sample number along x dimension
 * @param int Nz: the sample number along z dimension
 * @param int nlayer: extend layers
 * @LastEditors: LaurentCheung
 * @Version: 1.0
 */
void extendModel2d(float **modEx, float **mod, int Nx, int Nz, int nlayer)
{
  int ix, iz;
  /* center */
#pragma omp parallel for private(ix, iz) firstprivate(nlayer)
  for (ix = 0; ix < Nx - 2 * nlayer; ix++)
  {
#pragma ivdep
    for (iz = 0; iz < Nz - 2 * nlayer; iz++)
    {
      modEx[ix + nlayer][iz + nlayer] = mod[ix][iz];
    }
  }
  /* pading */
  padModel2d(modEx, Nx, Nz, nlayer);
}
/**
 * padding 3D model
 * @description: padding 2D model
 * @param {float**} mod: model array
 * @param int Nz: the sample number along z dimension
 * @param int Nx: the sample number along x dimension
 * @param int Ny: the sample number along y dimension
 * @param int nlayer: extend layers
 * @LastEditors: LaurentCheung
 * @Version: 1.0
 */
void padModel3d(float ***mod, int Nx, int Ny, int Nz, int nlayer)
{
  int ix, iy, iz, ilayer;
  /* padding top and bottom of model */
#pragma omp parallel for private(ix, iy, ilayer) firstprivate(nlayer)
  for (iy = nlayer; iy < Ny - nlayer; iy++)
  {
    for (ix = nlayer; ix < Nx - nlayer; ix++)
    {
#pragma ivdep
      for (ilayer = 1; ilayer <= nlayer; ilayer++)
      {
        mod[iy][ix][nlayer - ilayer] = mod[iy][ix][nlayer];
        mod[iy][ix][Nz - nlayer - 1 + ilayer] = mod[iy][ix][Nz - nlayer - 1];
      }
    }
  }
  /* padding left and right of model */
#pragma omp parallel for private(iy, iz, ilayer) firstprivate(nlayer)
  for (iy = nlayer; iy < Ny - nlayer; iy++)
  {
    for (ilayer = 1; ilayer <= nlayer; ilayer++)
    {
#pragma ivdep
      for (iz = 0; iz < Nz; iz++)
      {
        mod[iy][nlayer - ilayer][iz] = mod[iy][nlayer][iz];
        mod[iy][Nx - nlayer - 1 + ilayer][iz] = mod[iy][Nx - nlayer - 1][iz];
      }
    }
  }
  /* padding front and back of model */
#pragma omp parallel for private(ix, iz, ilayer) firstprivate(nlayer)
  for (ilayer = 1; ilayer <= nlayer; ilayer++)
  {
    for (ix = 0; ix < Nx; ix++)
    {
#pragma ivdep
      for (iz = 0; iz < Nz; iz++)
      {
        mod[nlayer - ilayer][ix][iz] = mod[nlayer][ix][iz];
        mod[Ny - nlayer - 1 + ilayer][ix][iz] = mod[Ny - nlayer - 1][ix][iz];
      }
    }
  }
}
/**
 * extend 3D model
 * @description: extend 3D model
 * @param {float**} mod: model array
 * @param int Nx: the sample number along x dimension
 * @param int Ny: the sample number along y dimension
 * @param int Nz: the sample number along z dimension
 * @param int nlayer: extend layers
 * @LastEditors: LaurentCheung
 * @Version: 1.0
 */
void extendModel3d(float ***modEx, float ***mod, int Nx, int Ny, int Nz, int nlayer)
{
  int ix, iy, iz;
  /* center */
#pragma omp parallel for private(ix, iy, iz) firstprivate(nlayer)
  for (iy = 0; iy < Ny - 2 * nlayer; iy++)
  {
    for (ix = 0; ix < Nx - 2 * nlayer; ix++)
    {
#pragma ivdep
      for (iz = 0; iz < Nz - 2 * nlayer; iz++)
      {
        modEx[iy + nlayer][ix + nlayer][iz + nlayer] = mod[iy][ix][iz];
      }
    }
  }
  /* pading */
  padModel3d(modEx, Nx, Ny, Nz, nlayer);
}