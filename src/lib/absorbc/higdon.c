/*
 * @description:
 * @Author: PingMin Cheung
 * @Date: 2024-03-26 10:13:53
 * @LastEditTime: 2024-05-11 21:07:32
 * @FilePath: /src/fdmod/lib/absorbc/higdon.c
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
 * @description:
 * @param float **p
 * @param float **pn
 * @param float **pn_1
 * @param float **vp
 * @param float dx
 * @param float dz
 * @param float dt
 * @param int Nx
 * @param int Nz
 * @param int nlayer
 * @param int freesurface
 * @return {*}
 * @Author: PingMin Cheung
 * @Date: 2024-03-26 10:59:27
 */
void acoustic2dClayton(float **p, float **pn, float **pn_1, float **vp, float dx, float dz, float dt, int Nx, int Nz,
                       int nlayer, int freesurface)
{

  int ix, iz;
  float qx, qz;
  float lambda, miu;
  qx = dt / dx;
  qz = dt / dz;
  /* top */
  if (freesurface == 0)
  {
#pragma omp parallel for private(ix, iz, lambda, miu) firstprivate(qz)
    for (ix = 0; ix < Nx; ix++)
    {
#pragma ivdep
      for (iz = nlayer - 1; iz >= 0; iz--)
      {
        lambda = 2.0f * vp[ix][iz] * qz;
        miu = vp[ix][iz] * vp[ix][iz] * qz * qz;

        p[ix][iz] = (2.0f - lambda - miu) * pn[ix][iz] + (lambda + 2.0f * miu) * pn[ix][iz + 1] -
                    miu * pn[ix][iz + 2] + (lambda - 1.0f) * pn_1[ix][iz] - lambda * pn_1[ix][iz + 1];
      }
    }
  }
  /* bottom */
#pragma omp parallel for private(ix, iz, lambda, miu) firstprivate(qz)
  for (ix = 0; ix < Nx; ix++)
  {
#pragma ivdep
    for (iz = Nz - nlayer; iz < Nz; iz++)
    {
      lambda = 2.0f * vp[ix][iz] * qz;
      miu = vp[ix][iz] * vp[ix][iz] * qz * qz;

      p[ix][iz] = (2.0f - lambda - miu) * pn[ix][iz] + (lambda + 2.0f * miu) * pn[ix][iz - 1] - miu * pn[ix][iz - 2] +
                  (lambda - 1.0f) * pn_1[ix][iz] - lambda * pn_1[ix][iz - 1];
    }
  }
  if (freesurface == 0)
  {
    /* left-top */
#pragma omp parallel for private(ix, iz, lambda, miu) firstprivate(qx)
    for (ix = nlayer - 1; ix >= 0; ix--)
    {
#pragma ivdep
      for (iz = 0; iz < nlayer; iz++)
      {
        lambda = 2.0f * vp[ix][iz] * qx;
        miu = vp[ix][iz] * vp[ix][iz] * qx * qx;

        p[ix][iz] =
            0.5f * ((2.0f - lambda - miu) * pn[ix][iz] + (lambda + 2.0f * miu) * pn[ix + 1][iz] -
                    miu * pn[ix + 2][iz] + (lambda - 1.0f) * pn_1[ix][iz] - lambda * pn_1[ix + 1][iz] + p[ix][iz]);
      }
    }
    /* right-top */
#pragma omp parallel for private(ix, iz, lambda, miu) firstprivate(qx)
    for (ix = Nx - nlayer; ix < Nx; ix++)
    {
#pragma ivdep
      for (iz = 0; iz < nlayer; iz++)
      {
        lambda = 2.0f * vp[ix][iz] * qx;
        miu = vp[ix][iz] * vp[ix][iz] * qx * qx;

        p[ix][iz] =
            0.5f * ((2.0f - lambda - miu) * pn[ix][iz] + (lambda + 2.0f * miu) * pn[ix - 1][iz] -
                    miu * pn[ix - 2][iz] + (lambda - 1.0f) * pn_1[ix][iz] - lambda * pn_1[ix - 1][iz] + p[ix][iz]);
      }
    }
  }
#pragma omp parallel for private(ix, iz, lambda, miu) firstprivate(qx)
  for (ix = nlayer - 1; ix >= 0; ix--)
  {
    /* left */
#pragma ivdep
    for (iz = nlayer; iz < Nz - nlayer; iz++)
    {
      lambda = 2.0f * vp[ix][iz] * qx;
      miu = vp[ix][iz] * vp[ix][iz] * qx * qx;

      p[ix][iz] = (2.0f - lambda - miu) * pn[ix][iz] + (lambda + 2.0f * miu) * pn[ix + 1][iz] - miu * pn[ix + 2][iz] +
                  (lambda - 1.0f) * pn_1[ix][iz] - lambda * pn_1[ix + 1][iz];
    }
    /* left-bottom */
#pragma ivdep
    for (iz = Nz - nlayer; iz < Nz; iz++)
    {
      lambda = 2.0f * vp[ix][iz] * qx;
      miu = vp[ix][iz] * vp[ix][iz] * qx * qx;

      p[ix][iz] =
          0.5f * ((2.0f - lambda - miu) * pn[ix][iz] + (lambda + 2.0f * miu) * pn[ix + 1][iz] - miu * pn[ix + 2][iz] +
                  (lambda - 1.0f) * pn_1[ix][iz] - lambda * pn_1[ix + 1][iz] + p[ix][iz]);
    }
  }
#pragma omp parallel for private(ix, iz, lambda, miu) firstprivate(qx)
  for (ix = Nx - nlayer; ix < Nx; ix++)
  {
    /* right */
#pragma ivdep
    for (iz = nlayer; iz < Nz - nlayer; iz++)
    {
      lambda = 2.0f * vp[ix][iz] * qx;
      miu = vp[ix][iz] * vp[ix][iz] * qx * qx;

      p[ix][iz] = (2.0f - lambda - miu) * pn[ix][iz] + (lambda + 2.0f * miu) * pn[ix - 1][iz] - miu * pn[ix - 2][iz] +
                  (lambda - 1.0f) * pn_1[ix][iz] - lambda * pn_1[ix - 1][iz];
    }
    /* right-bottom */
#pragma ivdep
    for (iz = Nz - nlayer; iz < Nz; iz++)
    {
      lambda = 2.0f * vp[ix][iz] * qx;
      miu = vp[ix][iz] * vp[ix][iz] * qx * qx;

      p[ix][iz] =
          0.5f * ((2.0f - lambda - miu) * pn[ix][iz] + (lambda + 2.0f * miu) * pn[ix - 1][iz] - miu * pn[ix - 2][iz] +
                  (lambda - 1.0f) * pn_1[ix][iz] - lambda * pn_1[ix - 1][iz] + p[ix][iz]);
    }
  }
}
/**
 * @description:
 * @param float **p
 * @param float **px
 * @param float **pz
 * @param float **vp
 * @param float dx
 * @param float dz
 * @param float dt
 * @param int Nx
 * @param int Nz
 * @param int nlayer
 * @param int freesurface
 * @return {*}
 * @Author: PingMin Cheung
 * @Date: 2024-03-26 10:59:27
 */
void Higdon2d(float **p, float **px, float **pz, float **vp, float dx, float dz, float dt, int Nx, int Nz, int nlayer,
              int freesurface)
{

  int ix, iz;
  float b, beta, bx, bz, qt;
  float rx, rz, qx, qz, qxt, qzt;
  beta = 1.0f;
  b = 0.4f;
  bx = dt / dx;
  bz = dt / dz;
  qt = b / (b - 1.0f);
  /* initial */
#pragma omp parallel for private(ix) firstprivate(Nz, nlayer)
  for (ix = 0; ix < Nx; ix++)
  {
    pz[ix][nlayer] = p[ix][nlayer];
    pz[ix][Nz - nlayer - 1] = p[ix][Nz - nlayer - 1];
  }
#pragma omp parallel for private(iz) firstprivate(Nx, nlayer)
  for (iz = 0; iz < Nz; iz++)
  {
    px[nlayer][iz] = p[nlayer][iz];
    px[Nx - nlayer - 1][iz] = p[Nx - nlayer - 1][iz];
  }
  /* left */
  if (freesurface == 0)
  {
    for (ix = nlayer - 1; ix >= 0; ix--)
    {
#pragma omp parallel for private(iz, rx, qx, qxt) firstprivate(ix, b, beta, bx, qt)
      for (iz = 0; iz < nlayer; iz++)
      {
        rx = bx * vp[ix][iz];
        qx = (b * (beta + rx) - rx) / ((beta + rx) * (1.0f - b));
        qxt = (b * (beta + rx) - beta) / ((beta + rx) * (1.0f - b));
        px[ix][iz] = -qx * px[ix + 1][iz] - qxt * p[ix][iz] - qt * p[ix + 1][iz];
      }
    }
  }
  for (ix = nlayer - 1; ix >= 0; ix--)
  {
#pragma omp parallel for private(iz, rx, qx, qxt) firstprivate(ix, b, beta, bx, qt)
    for (iz = nlayer; iz < Nz; iz++)
    {
      rx = bx * vp[ix][iz];
      qx = (b * (beta + rx) - rx) / ((beta + rx) * (1.0f - b));
      qxt = (b * (beta + rx) - beta) / ((beta + rx) * (1.0f - b));
      px[ix][iz] = -qx * px[ix + 1][iz] - qxt * p[ix][iz] - qt * p[ix + 1][iz];
    }
  }
  /* right */
  if (freesurface == 0)
  {
    for (ix = Nx - nlayer; ix < Nx; ix++)
    {
#pragma omp parallel for private(iz, rx, qx, qxt) firstprivate(ix, b, beta, bx, qt)
      for (iz = 0; iz < nlayer; iz++)
      {
        rx = bx * vp[ix][iz];
        qx = (b * (beta + rx) - rx) / ((beta + rx) * (1.0f - b));
        qxt = (b * (beta + rx) - beta) / ((beta + rx) * (1.0f - b));
        px[ix][iz] = -qx * px[ix - 1][iz] - qxt * p[ix][iz] - qt * p[ix - 1][iz];
      }
    }
  }
  for (ix = Nx - nlayer; ix < Nx; ix++)
  {
#pragma omp parallel for private(iz, rx, qx, qxt) firstprivate(ix, b, beta, bx, qt)
    for (iz = nlayer; iz < Nz; iz++)
    {
      rx = bx * vp[ix][iz];
      qx = (b * (beta + rx) - rx) / ((beta + rx) * (1.0f - b));
      qxt = (b * (beta + rx) - beta) / ((beta + rx) * (1.0f - b));
      px[ix][iz] = -qx * px[ix - 1][iz] - qxt * p[ix][iz] - qt * p[ix - 1][iz];
    }
  }
  /* top */
  if (freesurface == 0)
  {
    for (iz = nlayer - 1; iz >= 0; iz--)
    {
#pragma omp parallel for private(ix, rz, qz, qzt) firstprivate(iz, b, beta, bz, qt)
      for (ix = 0; ix < Nx; ix++)
      {
        rz = bz * vp[ix][iz];
        qz = (b * (beta + rz) - rz) / ((beta + rz) * (1.0f - b));
        qzt = (b * (beta + rz) - beta) / ((beta + rz) * (1.0f - b));
        pz[ix][iz] = -qz * pz[ix][iz + 1] - qzt * p[ix][iz] - qt * p[ix][iz + 1];
      }
    }
  }
  /* bottom */
  for (iz = Nz - nlayer; iz < Nz; iz++)
  {
#pragma omp parallel for private(ix, rz, qz, qzt) firstprivate(iz, b, beta, bz, qt)
    for (ix = 0; ix < Nx; ix++)
    {
      rz = bz * vp[ix][iz];
      qz = (b * (beta + rz) - rz) / ((beta + rz) * (1.0f - b));
      qzt = (b * (beta + rz) - beta) / ((beta + rz) * (1.0f - b));
      pz[ix][iz] = -qz * pz[ix][iz - 1] - qzt * p[ix][iz] - qt * p[ix][iz - 1];
    }
  }
  /* setting boundary */
  if (freesurface == 0)
  {
#pragma omp parallel for private(ix, iz)
    for (ix = 0; ix < nlayer; ix++)
    {
#pragma ivdep
      for (iz = 0; iz < nlayer; iz++)
      {
        p[ix][iz] = 0.5f * (px[ix][iz] + pz[ix][iz]);
      }
    }
#pragma omp parallel for private(ix, iz)
    for (ix = nlayer; ix < Nx - nlayer; ix++)
    {
#pragma ivdep
      for (iz = 0; iz < nlayer; iz++)
      {
        p[ix][iz] = pz[ix][iz];
      }
    }
#pragma omp parallel for private(ix, iz)
    for (ix = Nx - nlayer; ix < Nx; ix++)
    {
#pragma ivdep
      for (iz = 0; iz < nlayer; iz++)
      {
        p[ix][iz] = 0.5f * (px[ix][iz] + pz[ix][iz]);
      }
    }
  }
#pragma omp parallel for private(ix, iz)
  for (ix = nlayer; ix < Nx - nlayer; ix++)
  {
#pragma ivdep
    for (iz = Nz - nlayer; iz < Nz; iz++)
    {
      p[ix][iz] = pz[ix][iz];
    }
  }
#pragma omp parallel for private(ix, iz)
  for (ix = 0; ix < nlayer; ix++)
  {
#pragma ivdep
    for (iz = nlayer; iz < Nz - nlayer; iz++)
    {
      p[ix][iz] = px[ix][iz];
    }
  }
#pragma omp parallel for private(ix, iz)
  for (ix = Nx - nlayer; ix < Nx; ix++)
  {
#pragma ivdep
    for (iz = nlayer; iz < Nz - nlayer; iz++)
    {
      p[ix][iz] = px[ix][iz];
    }
  }
#pragma omp parallel for private(ix, iz)
  for (ix = 0; ix < nlayer; ix++)
  {
#pragma ivdep
    for (iz = Nz - nlayer; iz < Nz; iz++)
    {
      p[ix][iz] = 0.5f * (px[ix][iz] + pz[ix][iz]);
    }
  }
#pragma omp parallel for private(ix, iz)
  for (ix = Nx - nlayer; ix < Nx; ix++)
  {
#pragma ivdep
    for (iz = Nz - nlayer; iz < Nz; iz++)
    {
      p[ix][iz] = 0.5f * (px[ix][iz] + pz[ix][iz]);
    }
  }
}