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
 * @param float alpha
 * @param int l
 * @param int nlayer
 * @return {*
 * @Author: PingMin Cheung
 * @Date: 2023-06-11 14:05:37
 */
static float absorbWeight(float alpha, int l, int nlayer)
{
  float w, d;
  d = (0.1f * l) / (0.1f * nlayer);
  w = expf(logf(2.0f) * powf(d, alpha)) - 1.0f;
  return w;
}
/**
 * 2D MEAL boundary condition
 * @description: 2D MEAL boundary condition
 * @param float** D: the value of damping function
 * @param float** Beta: scaling factor
 * @param float** vp: velocity of P-wave
 * @param float** vs: velocity of S-wave
 * @param float fpeak: the dominant frequency of the source
 * @param float alpha: the order of Gaussian weight function
 * @param float dz: the sampling interval along z, the 1st dimension
 * @param float dx: the sampling interval along x, the 2nd dimension
 * @param int Nz: The sample number along 1st dimension
 * @param int Nx: The sample number along 2nd dimension
 * @param int nlayer: number of EAL
 * @param int freesurface: is free-surface boundary condition used.
 *                         0 for not, 1 for yes
 * @return float** D: the value of damping function
 * @return float** Beta: scaling factor
 * @LastEditors: LaurentCheung
 * @Version: 1.0
 */
void mealCoefs2d(float **D, float **Beta, float **vp, float **vs, float fpeak, float alpha, float dz, float dx,
                 int Nz, int Nx, int nlayer, int freesurface)
{
  int ix, iz, izStart, l;
  float weight;
  float pd, Rcoef, PPWfd, PPWfcs;
  float thicknessx, thicknessz;
  float d0x, d0z;
  float betax, betaz;
  float **Vs = NULL;
  float **Dx = NULL, **Dz = NULL;
  float **Betax = NULL, **Betaz = NULL;
  izStart = freesurface == 0 ? 0 : nlayer;
  pd = 2;
  Rcoef = powf(10.0f, -1.0f * ((log10f(nlayer) - 1.0f) / (log10f(2.0f))) - 3.0f - 0.5f * pd);
  PPWfd = 10.0f;
  thicknessx = nlayer * dx;
  thicknessz = nlayer * dz;
  Vs = alloc2float(Nz, Nx);
  Dx = alloc2float(Nz, Nx);
  Dz = alloc2float(Nz, Nx);
  Betax = alloc2float(Nz, Nx);
  Betaz = alloc2float(Nz, Nx);
  if (alpha == -1.0f)
  {
    alpha = 2.5f;
  }
  if (vs == NULL)
  {
#pragma omp parallel for private(ix, iz)
    for (ix = 0; ix < Nx; ix++)
    {
#pragma ivdep
      for (iz = 0; iz < Nz; iz++)
      {
        Vs[ix][iz] = vp[ix][iz] / sqrtf(3.0f);
      }
    }
  }
  else
  {
#pragma omp parallel for private(ix, iz)
    for (ix = 0; ix < Nx; ix++)
    {
#pragma ivdep
      for (iz = 0; iz < Nz; iz++)
      {
        if (vs[ix][iz] == 0.0f)
        {
          Vs[ix][iz] = vp[ix][iz] / sqrtf(3.0f);
        }
        else
        {
          Vs[ix][iz] = vs[ix][iz];
        }
      }
    }
  }
  /* initialization */
#pragma omp parallel for private(ix, iz)
  for (ix = 0; ix < Nx; ix++)
  {
#pragma ivdep
    for (iz = 0; iz < Nz; iz++)
    {
      D[ix][iz] = 0.0f;
      Dx[ix][iz] = 0.0f;
      Dz[ix][iz] = 0.0f;
      Beta[ix][iz] = 1.0f;
      Betax[ix][iz] = 1.0f;
      Betaz[ix][iz] = 1.0f;
    }
  }
  /* top and bottom boundary */
  if (freesurface == 0)
  {
#pragma omp parallel for private(ix, iz, l, weight, d0z, PPWfcs, betaz) firstprivate(PPWfd, Rcoef, thicknessz, dz)
    for (ix = 0; ix < Nx; ix++)
    {
#pragma ivdep
      for (l = 1; l <= nlayer; l++)
      {
        weight = absorbWeight(alpha, l, nlayer);
        /* top */
        iz = nlayer - l;
        d0z = -1.0f * (1.0f + pd) * vp[ix][iz] * logf(Rcoef) / (2.0f * thicknessz);
        PPWfcs = Vs[ix][iz] / (dz * fpeak);
        betaz = 2.0f * PPWfcs / PPWfd;
        Dz[ix][iz] = d0z * weight;
        Betaz[ix][iz] = 1.0f + (betaz - 1.0f) * weight;
      }
    }
  }
#pragma omp parallel for private(ix, iz, l, weight, d0z, PPWfcs, betaz) firstprivate(PPWfd, Rcoef, thicknessz, dz)
  for (ix = 0; ix < Nx; ix++)
  {
#pragma ivdep
    for (l = 1; l <= nlayer; l++)
    {
      weight = absorbWeight(alpha, l, nlayer);
      /* bottom */
      iz = Nz - nlayer - 1 + l;
      d0z = -1.0f * (1.0f + pd) * vp[ix][iz] * logf(Rcoef) / (2.0f * thicknessz);
      PPWfcs = Vs[ix][iz] / (dz * fpeak);
      betaz = 2.0f * PPWfcs / PPWfd;
      Dz[ix][iz] = d0z * weight;
      Betaz[ix][iz] = 1.0f + (betaz - 1.0f) * weight;
    }
  }
  /* left and right boundary */
#pragma omp parallel for private(ix, iz, l, weight, d0x, PPWfcs, betax) firstprivate(PPWfd, Rcoef, thicknessx, dx)
  for (l = 1; l <= nlayer; l++)
  {
    weight = absorbWeight(alpha, l, nlayer);
#pragma ivdep
    for (iz = izStart; iz < Nz; iz++)
    {
      /* left */
      ix = nlayer - l;
      d0x = -1.0f * (1.0f + pd) * vp[ix][iz] * logf(Rcoef) / (2.0f * thicknessx);
      PPWfcs = Vs[ix][iz] / (dx * fpeak);
      betax = 2.0f * PPWfcs / PPWfd;
      Dx[ix][iz] = d0x * weight;
      Betax[ix][iz] = 1.0f + (betax - 1.0f) * weight;
      /* right */
      ix = Nx - nlayer - 1 + l;
      d0x = -1.0f * (1.0f + pd) * vp[ix][iz] * logf(Rcoef) / (2.0f * thicknessx);
      PPWfcs = Vs[ix][iz] / (dx * fpeak);
      betax = 2.0f * PPWfcs / PPWfd;
      Dx[ix][iz] = d0x * weight;
      Betax[ix][iz] = 1.0f + (betax - 1.0f) * weight;
    }
  }
  /* absorbing coefficent */
#pragma omp parallel for private(ix, iz)
  for (ix = 0; ix < Nx; ix++)
  {
#pragma ivdep
    for (iz = 0; iz < Nz; iz++)
    {
      D[ix][iz] = sqrtf(powf(Dx[ix][iz], 2.0f) + powf(Dz[ix][iz], 2.0f));
      Beta[ix][iz] = sqrtf(powf(Betax[ix][iz] - 1.0f, 2.0f) + powf(Betaz[ix][iz] - 1.0f, 2.0f)) + 1.0f;
    }
  }
  /* free memory */
  free2float(Vs);
  free2float(Dx);
  free2float(Dz);
  free2float(Betax);
  free2float(Betaz);
}
/**
 * 3D MEAL boundary condition
 * @description: 2D MEAL boundary condition
 * @param float*** D: the value of damping function
 * @param float*** Beta: scaling factor
 * @param float*** vp: velocity of P-wave
 * @param float*** vs: velocity of S-wave
 * @param float fpeak: the dominant frequency of the source
 * @param float alpha: the order of Gaussian weight function
 * @param float dz: the sampling interval along z, the 1st dimension
 * @param float dx: the sampling interval along x, the 2nd dimension
 * @param float dy: the sampling interval along y, the 3rd dimension
 * @param int Nz: The sample number along 1st dimension
 * @param int Nx: The sample number along 2nd dimension
 * @param int Ny: The sample number along 3rd dimension
 * @param int nlayer: number of EAL
 * @param int freesurface: is free-surface boundary condition used.
 *                         0 for not, 1 for yes
 * @return float*** D: the value of damping function
 * @return float*** Beta: scaling factor
 * @LastEditors: LaurentCheung
 * @Version: 1.0
 */
void mealCoefs3d(float ***D, float ***Beta, float ***vp, float ***vs, float fpeak, float alpha, float dz, float dx,
                 float dy, int Nz, int Nx, int Ny, int nlayer, int freesurface)
{
  int ix, iy, iz, izStart, l;
  float weight;
  float pd, Rcoef, PPWfd, PPWfcs;
  float thicknessx, thicknessy, thicknessz;
  float d0x, d0y, d0z;
  float betax, betay, betaz;
  float ***Vs = NULL;
  float ***Dx = NULL, ***Dy = NULL, ***Dz = NULL;
  float ***Betax = NULL, ***Betay = NULL, ***Betaz = NULL;
  izStart = freesurface == 0 ? 0 : nlayer;
  pd = 2;
  Rcoef = powf(10.0f, -1.0f * ((log10f(nlayer) - 1.0f) / (log10f(2.0f))) - 3.0f - 0.5f * pd);
  PPWfd = 10.0f;
  thicknessx = nlayer * dx;
  thicknessy = nlayer * dy;
  thicknessz = nlayer * dz;
  Vs = alloc3float(Nz, Nx, Ny);
  Dx = alloc3float(Nz, Nx, Ny);
  Dy = alloc3float(Nz, Nx, Ny);
  Dz = alloc3float(Nz, Nx, Ny);
  Betax = alloc3float(Nz, Nx, Ny);
  Betay = alloc3float(Nz, Nx, Ny);
  Betaz = alloc3float(Nz, Nx, Ny);
  if (alpha == -1.0f)
  {
    alpha = 2.5f;
  }
  if (vs == NULL)
  {
#pragma omp parallel for private(ix, iy, iz)
    for (iy = 0; iy < Ny; iy++)
    {
      for (ix = 0; ix < Nx; ix++)
      {
#pragma ivdep
        for (iz = 0; iz < Nz; iz++)
        {
          Vs[iy][ix][iz] = vp[iy][ix][iz] / sqrtf(3.0f);
        }
      }
    }
  }
  else
  {
#pragma omp parallel for private(ix, iy, iz)
    for (iy = 0; iy < Ny; iy++)
    {
      for (ix = 0; ix < Nx; ix++)
      {
#pragma ivdep
        for (iz = 0; iz < Nz; iz++)
        {
          if (vs[iy][ix][iz] == 0.0f)
          {
            Vs[iy][ix][iz] = vp[iy][ix][iz] / sqrtf(3.0f);
          }
          else
          {
            Vs[iy][ix][iz] = vs[iy][ix][iz];
          }
        }
      }
    }
  }
  /* initialization */
#pragma omp parallel for private(ix, iy, iz)
  for (iy = 0; iy < Ny; iy++)
  {
    for (ix = 0; ix < Nx; ix++)
    {
#pragma ivdep
      for (iz = 0; iz < Nz; iz++)
      {
        D[iy][ix][iz] = 0.0f;
        Dx[iy][ix][iz] = 0.0f;
        Dy[iy][ix][iz] = 0.0f;
        Dz[iy][ix][iz] = 0.0f;
        Beta[iy][ix][iz] = 1.0f;
        Betax[iy][ix][iz] = 1.0f;
        Betay[iy][ix][iz] = 1.0f;
        Betaz[iy][ix][iz] = 1.0f;
      }
    }
  }
  /* top and bottom boundary */
  if (freesurface == 0)
  {
#pragma omp parallel for private(ix, iy, iz, l, weight, d0z, PPWfcs, betaz) firstprivate(PPWfd, Rcoef, thicknessz, dz)
    for (iy = 0; iy < Ny; iy++)
    {
      for (ix = 0; ix < Nx; ix++)
      {
#pragma ivdep
        for (l = 1; l <= nlayer; l++)
        {
          weight = absorbWeight(alpha, l, nlayer);
          /* top */
          iz = nlayer - l;
          d0z = -1.0f * (1.0f + pd) * vp[iy][ix][iz] * logf(Rcoef) / (2.0f * thicknessz);
          PPWfcs = Vs[iy][ix][iz] / (dz * fpeak);
          betaz = 2.0f * PPWfcs / PPWfd;
          Dz[iy][ix][iz] = d0z * weight;
          Betaz[iy][ix][iz] = 1.0f + (betaz - 1.0f) * weight;
        }
      }
    }
  }
#pragma omp parallel for private(ix, iy, iz, l, weight, d0z, PPWfcs, betaz) firstprivate(PPWfd, Rcoef, thicknessz, dz)
  for (iy = 0; iy < Ny; iy++)
  {
    for (ix = 0; ix < Nx; ix++)
    {
#pragma ivdep
      for (l = 1; l <= nlayer; l++)
      {
        weight = absorbWeight(alpha, l, nlayer);
        /* bottom */
        iz = Nz - nlayer - 1 + l;
        d0z = -1.0f * (1.0f + pd) * vp[iy][ix][iz] * logf(Rcoef) / (2.0f * thicknessz);
        PPWfcs = Vs[iy][ix][iz] / (dz * fpeak);
        betaz = 2.0f * PPWfcs / PPWfd;
        Dz[iy][ix][iz] = d0z * weight;
        Betaz[iy][ix][iz] = 1.0f + (betaz - 1.0f) * weight;
      }
    }
  }
  /* left and right boundary */
#pragma omp parallel for private(ix, iy, iz, l, weight, d0x, PPWfcs, betax) firstprivate(PPWfd, Rcoef, thicknessx, dx)
  for (iy = 0; iy < Ny; iy++)
  {
    for (l = 1; l <= nlayer; l++)
    {
      weight = absorbWeight(alpha, l, nlayer);
#pragma ivdep
      for (iz = izStart; iz < Nz; iz++)
      {
        /* left */
        ix = nlayer - l;
        d0x = -1.0f * (1.0f + pd) * vp[iy][ix][iz] * logf(Rcoef) / (2.0f * thicknessx);
        PPWfcs = Vs[iy][ix][iz] / (dx * fpeak);
        betax = 2.0f * PPWfcs / PPWfd;
        Dx[iy][ix][iz] = d0x * weight;
        Betax[iy][ix][iz] = 1.0f + (betax - 1.0f) * weight;
        /* right */
        ix = Nx - nlayer - 1 + l;
        d0x = -1.0f * (1.0f + pd) * vp[iy][ix][iz] * logf(Rcoef) / (2.0f * thicknessx);
        PPWfcs = Vs[iy][ix][iz] / (dx * fpeak);
        betax = 2.0f * PPWfcs / PPWfd;
        Dx[iy][ix][iz] = d0x * weight;
        Betax[iy][ix][iz] = 1.0f + (betax - 1.0f) * weight;
      }
    }
  }
  /* front and back boundary */
#pragma omp parallel for private(ix, iy, iz, l, weight, d0y, PPWfcs, betay) firstprivate(PPWfd, Rcoef, thicknessy, dy)
  for (l = 1; l <= nlayer; l++)
  {
    weight = absorbWeight(alpha, l, nlayer);
    for (ix = izStart; ix < Nx; ix++)
    {
#pragma ivdep
      for (iz = izStart; iz < Nz; iz++)
      {
        /* front */
        iy = nlayer - l;
        d0y = -1.0f * (1.0f + pd) * vp[iy][ix][iz] * logf(Rcoef) / (2.0f * thicknessy);
        PPWfcs = Vs[iy][ix][iz] / (dy * fpeak);
        betay = 2.0f * PPWfcs / PPWfd;
        Dy[iy][ix][iz] = d0y * weight;
        Betay[iy][ix][iz] = 1.0f + (betay - 1.0f) * weight;
        /* back */
        iy = Ny - nlayer - 1 + l;
        d0y = -1.0f * (1.0f + pd) * vp[iy][ix][iz] * logf(Rcoef) / (2.0f * thicknessy);
        PPWfcs = Vs[iy][ix][iz] / (dy * fpeak);
        betay = 2.0f * PPWfcs / PPWfd;
        Dy[iy][ix][iz] = d0y * weight;
        Betay[iy][ix][iz] = 1.0f + (betay - 1.0f) * weight;
      }
    }
  }
  /* absorbing coefficent */
#pragma omp parallel for private(ix, iy, iz)
  for (iy = 0; iy < Ny; iy++)
  {
    for (ix = 0; ix < Nx; ix++)
    {
#pragma ivdep
      for (iz = 0; iz < Nz; iz++)
      {
        D[iy][ix][iz] = sqrtf(powf(Dx[iy][ix][iz], 2.0f) + powf(Dy[iy][ix][iz], 2.0f) + powf(Dz[iy][ix][iz], 2.0f));
        Beta[iy][ix][iz] = sqrtf(powf(Betax[iy][ix][iz] - 1.0f, 2.0f) + powf(Betay[iy][ix][iz] - 1.0f, 2.0f) +
                                 powf(Betaz[iy][ix][iz] - 1.0f, 2.0f)) +
                           1.0f;
      }
    }
  }
  /* free memory */
  free3float(Vs);
  free3float(Dx);
  free3float(Dy);
  free3float(Dz);
  free3float(Betax);
  free3float(Betay);
  free3float(Betaz);
}

/**
 * calculate the coefficent of finite difference scheme for first or second order equation
 * @param float *coef1: first coefficent
 * @param float *coef2: second coefficent
 * @param float *coef3: third coefficent
 * @param float *D: the value of damping function
 * @param float *Beta: scaling factor
 * @param float *vp: velocity of P-wave
 * @param float *vs: velocity of S-wave
 * @param float *rho: density
 * @param float dt: time stepsize
 * @param int Np: the number of points
 * @param int order: the order of equation
 * @Author: Pingmin Cheung
 * @Date: 2025-08-26 10:22:59
 */
void mealFDCoefs(float *coef1, float *coef2, float *coef3, float *D, float *Beta, float *vp, float *vs, float *rho,
                 float dt, int Np, int order)
{
  int i;
  if (vp != NULL)
  {
#pragma omp parallel for private(i)
    for (i = 0; i < Np; i++)
    {
      vp[i] /= Beta[i];
    }
  }
  if (vs != NULL)
  {
#pragma omp parallel for private(i)
    for (i = 0; i < Np; i++)
    {
      vs[i] /= Beta[i];
    }
  }
  if (rho != NULL)
  {
#pragma omp parallel for private(i)
    for (i = 0; i < Np; i++)
    {
      rho[i] *= Beta[i];
    }
  }
#pragma omp parallel for private(i)
  for (i = 0; i < Np; i++)
  {
    D[i] /= Beta[i];
  }
  if (order == 1 && coef1 != NULL && coef2 != NULL)
  {
#pragma omp parallel for private(i)
    for (i = 0; i < Np; i++)
    {
      coef1[i] = (1.0 - 0.5 * dt * D[i]) / (1.0 + 0.5 * dt * D[i]);
      coef2[i] = dt / (1.0 + 0.5 * dt * D[i]);
    }
  }
  else if (order == 2 && coef1 != NULL && coef2 != NULL && coef3 != NULL)
  {
#pragma omp parallel for private(i)
    for (i = 0; i < Np; i++)
    {
      coef1[i] = (2.0 - D[i] * D[i] * dt * dt) / (1.0 + dt * D[i]);
      coef2[i] = -(1.0 - dt * D[i]) / (1.0 + dt * D[i]);
      coef3[i] = dt * dt / (1.0 + dt * D[i]);
    }
  }
}