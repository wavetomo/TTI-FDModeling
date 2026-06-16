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
/********** elastic parameters **********/
/**
 * calculate the Lame constants from P- and S-wave velocities and density
 * @param float *lambda: Lame constants, lambda
 * @param float *miu: Lame constants, miu
 * @param float *vp: P-wave velocity
 * @param float *vs: S-wave velocity
 * @param float *rho: density
 * @param int np: the length of pointer
 * @return float *lambda: Lame constants, lambda, length: np
 * @return float *miu: Lame constants, miu, length: np
 * @Author: PingMin Cheung
 * @Date: 2024-08-01 10:29:15
 */
void lameConstants(float *lambda, float *miu, float *vp, float *vs, float *rho, int np)
{
  int i;
  if (vs == NULL)
  {
#pragma omp parallel for private(i)
    for (i = 0; i < np; i++)
    {
      lambda[i] = rho[i] * vp[i] * vp[i];
    }
  }
  else
  {
#pragma omp parallel for private(i)
    for (i = 0; i < np; i++)
    {
      lambda[i] = rho[i] * (vp[i] * vp[i] - 2.0f * vs[i] * vs[i]);
      miu[i] = rho[i] * vs[i] * vs[i];
    }
  }
}
/**
 * calculate 2D elastic coefficients from
 * P- and S-wave velocities, density
 * @param float *C11: elastic coefficient
 * @param float *C13: elastic coefficient
 * @param float *C55: elastic coefficient
 * @param float *vp: P-wave velocity
 * @param float *vs: S-wave velocity
 * @param float *rho: density
 * @param int np: the length of pointers
 * @return float *Cij: elastic coefficient
 * @Author: PingMin Cheung
 * @Date: 2024-08-01 10:37:34
 */
void elasticCoefficients2d(float *C11, float *C13, float *C55, float *vp, float *vs, float *rho, int np)
{
  int i;
  if (vs == NULL)
  {
#pragma omp parallel for private(i)
    for (i = 0; i < np; i++)
    {
      C11[i] = rho[i] * vp[i] * vp[i];
    }
  }
  else
  {
#pragma omp parallel for private(i)
    for (i = 0; i < np; i++)
    {
      C11[i] = rho[i] * vp[i] * vp[i];
      C13[i] = rho[i] * (vp[i] * vp[i] - 2.0 * vs[i] * vs[i]);
      C55[i] = rho[i] * vs[i] * vs[i];
    }
  }
}
/********** receive wavefield **********/
/**
 *
 * @param {Receiver2D} receiver_vx
 * @param {Receiver2D} receiver_vz
 * @param {Receiver2D} receiver_p
 * @param float **vx
 * @param float **vz
 * @param float **sigma_xx
 * @param float **sigma_zz
 * @param int it
 * @param int mt
 * @param int ealThick
 * @param int seistype
 * @Author: PingMin Cheung
 * @Date: 2024-05-16 21:33:26
 */
void elastic2dReceiver(Receiver2D receiver_vx, Receiver2D receiver_vz, Receiver2D receiver_p, float **vx, float **vz,
                       float **sigma_xx, float **sigma_zz, int it, int mt, int ealThick, int seistype)
{
  int itrace, irecv;
  int igx, igz;
  if (it % mt == 0)
  {
    irecv = it / mt;
    switch (seistype)
    {
    case 1:
#pragma omp parallel for private(itrace, igx, igz) firstprivate(irecv, ealThick)
      for (itrace = 0; itrace < receiver_vx.ntrace; itrace++)
      {
        igx = receiver_vx.posx[itrace] + ealThick;
        igz = receiver_vx.posz[itrace] + ealThick;
        receiver_vx.recv[itrace][irecv] = vx[igx][igz];
      }
#pragma omp parallel for private(itrace, igx, igz) firstprivate(irecv, ealThick)
      for (itrace = 0; itrace < receiver_vz.ntrace; itrace++)
      {
        igx = receiver_vz.posx[itrace] + ealThick;
        igz = receiver_vz.posz[itrace] + ealThick;
        receiver_vz.recv[itrace][irecv] = vz[igx][igz];
      }
      break;
    case 2:
#pragma omp parallel for private(itrace, igx, igz) firstprivate(irecv, ealThick)
      for (itrace = 0; itrace < receiver_p.ntrace; itrace++)
      {
        igx = receiver_p.posx[itrace] + ealThick;
        igz = receiver_p.posz[itrace] + ealThick;
        receiver_p.recv[itrace][irecv] = 0.5f * (sigma_xx[igx][igz] + sigma_zz[igx][igz]);
      }
      break;
    case 3:
      elastic2dReceiver(receiver_vx, receiver_vz, receiver_p, vx, vz, sigma_xx, sigma_zz, it, mt, ealThick, 1);
      elastic2dReceiver(receiver_vx, receiver_vz, receiver_p, vx, vz, sigma_xx, sigma_zz, it, mt, ealThick, 2);
      break;
    }
  }
}
/**
 *
 * @param {Receiver3D} receiver_vx
 * @param {Receiver3D} receiver_vy
 * @param {Receiver3D} receiver_vz
 * @param {Receiver3D} receiver_p
 * @param float ***vx
 * @param float ***vy
 * @param float ***vz
 * @param float **sigma_xx
 * @param float **sigma_yy
 * @param float **sigma_zz
 * @param int it
 * @param int mt
 * @param int ealThick
 * @param int seistype
 * @return {*}
 * @Author: PingMin Cheung
 * @Date: 2024-05-17 14:30:57
 */
void elastic3dReceiver(Receiver3D receiver_p, Receiver3D receiver_vx, Receiver3D receiver_vy, Receiver3D receiver_vz,
                       float ***vx, float ***vy, float ***vz, float ***sigma_xx, float ***sigma_yy, float ***sigma_zz,
                       int it, int mt, int ealThick, int seistype)
{
  int itrace, irecv;
  int igx, igy, igz;
  if (it % mt == 0)
  {
    irecv = it / mt;
    switch (seistype)
    {
    case 1:
#pragma omp parallel for private(itrace, igx, igy, igz) firstprivate(irecv, ealThick)
      for (itrace = 0; itrace < receiver_vy.ntrace; itrace++)
      {
        igy = receiver_vx.posy[itrace] + ealThick;
        igx = receiver_vx.posx[itrace] + ealThick;
        igz = receiver_vx.posz[itrace] + ealThick;
        receiver_vx.recv[itrace][irecv] = vx[igy][igx][igz];
      }
#pragma omp parallel for private(itrace, igx, igy, igz) firstprivate(irecv, ealThick)
      for (itrace = 0; itrace < receiver_vz.ntrace; itrace++)
      {
        igy = receiver_vy.posy[itrace] + ealThick;
        igx = receiver_vy.posx[itrace] + ealThick;
        igz = receiver_vy.posz[itrace] + ealThick;
        receiver_vy.recv[itrace][irecv] = vy[igy][igx][igz];
      }
#pragma omp parallel for private(itrace, igx, igy, igz) firstprivate(irecv, ealThick)
      for (itrace = 0; itrace < receiver_vz.ntrace; itrace++)
      {
        igy = receiver_vz.posy[itrace] + ealThick;
        igx = receiver_vz.posx[itrace] + ealThick;
        igz = receiver_vz.posz[itrace] + ealThick;
        receiver_vz.recv[itrace][irecv] = vz[igy][igx][igz];
      }
      break;
    case 2:
#pragma omp parallel for private(itrace, igx, igy, igz) firstprivate(irecv, ealThick)
      for (itrace = 0; itrace < receiver_p.ntrace; itrace++)
      {
        igy = receiver_p.posy[itrace] + ealThick;
        igx = receiver_p.posx[itrace] + ealThick;
        igz = receiver_p.posz[itrace] + ealThick;
        receiver_p.recv[itrace][irecv] =
            0.333333333f * (sigma_xx[igy][igx][igz] + sigma_yy[igy][igx][igz] + sigma_zz[igy][igx][igz]);
      }
      break;
    case 3:
      elastic3dReceiver(receiver_p, receiver_vx, receiver_vy, receiver_vz, vx, vy, vz, sigma_xx, sigma_yy, sigma_zz,
                        it, mt, ealThick, 1);
      elastic3dReceiver(receiver_p, receiver_vx, receiver_vy, receiver_vz, vx, vy, vz, sigma_xx, sigma_yy, sigma_zz,
                        it, mt, ealThick, 2);
      break;
    }
  }
}
/********** snapshot save **********/
/**
 *
 * @param {char} *prefix
 * @param float **p
 * @param float **vx
 * @param float **vz
 * @param float dt
 * @param int it
 * @param int Nz
 * @param int Nx
 * @param int ealThick
 * @param int savew
 * @param int wsave
 * @return {*}
 * @LastEditors: LaurentCheung
 * @Version: 1.0
 */
void elastic2dSnapshot(char *prefix, float **vx, float **vz, float **sigma_xx, float **sigma_zz, float dt, int it,
                       int Nz, int Nx, int ealThick, int savew, int wsave)
{
  char fileName[256];
  int nz, nx;
  int ix, iz;
  float **p = NULL;
  nz = savew == 0 ? Nz - 2 * ealThick : Nz;
  nx = savew == 0 ? Nx - 2 * ealThick : Nx;
  switch (wsave)
  {
  case 1:
    sprintf(fileName, "%s_vx_%05d_%05d_%05dms.bin", prefix, nz, nx, (int)(1000 * it * dt));
    write_snapshot2d(fileName, vx, Nz, Nx, ealThick, savew);
    sprintf(fileName, "%s_vz_%05d_%05d_%05dms.bin", prefix, nz, nx, (int)(1000 * it * dt));
    write_snapshot2d(fileName, vz, Nz, Nx, ealThick, savew);
    break;
  case 2:
    p = alloc2float(Nz, Nx);
#pragma omp parallel for private(ix, iz) firstprivate(savew, ealThick, Nx, Nz)
    for (ix = (savew == 0 ? ealThick : 0); ix < (savew == 0 ? Nx - ealThick : Nx); ix++)
    {
#pragma ivdep
      for (iz = (savew == 0 ? ealThick : 0); iz < (savew == 0 ? Nz - ealThick : Nz); iz++)
      {
        p[ix][iz] = 0.5f * (sigma_xx[ix][iz] + sigma_zz[ix][iz]);
      }
    }
    sprintf(fileName, "%s_p_%05d_%05d_%05dms.bin", prefix, nz, nx, (int)(1000 * it * dt));
    write_snapshot2d(fileName, p, Nz, Nx, ealThick, savew);
    free2float(p);
    break;
  case 3:
    elastic2dSnapshot(prefix, vx, vz, sigma_xx, sigma_zz, dt, it, Nz, Nx, ealThick, savew, 1);
    elastic2dSnapshot(prefix, vx, vz, sigma_xx, sigma_zz, dt, it, Nz, Nx, ealThick, savew, 2);
    break;
  case 4:
    sprintf(fileName, "%s_vxp_%05d_%05d_%05dms.bin", prefix, nz, nx, (int)(1000 * it * dt));
    write_snapshot2d(fileName, vx, Nz, Nx, ealThick, savew);
    sprintf(fileName, "%s_vzp_%05d_%05d_%05dms.bin", prefix, nz, nx, (int)(1000 * it * dt));
    write_snapshot2d(fileName, vz, Nz, Nx, ealThick, savew);
    break;
  case 5:
    sprintf(fileName, "%s_vxs_%05d_%05d_%05dms.bin", prefix, nz, nx, (int)(1000 * it * dt));
    write_snapshot2d(fileName, vx, Nz, Nx, ealThick, savew);
    sprintf(fileName, "%s_vzs_%05d_%05d_%05dms.bin", prefix, nz, nx, (int)(1000 * it * dt));
    write_snapshot2d(fileName, vz, Nz, Nx, ealThick, savew);
    break;
  }
}
/**
 *
 * @param {char} *prefix
 * @param float ***vx
 * @param float ***vy
 * @param float ***vz
 * @param float ***sigma_xx
 * @param float ***sigma_yy
 * @param float ***sigma_zz
 * @param float dt
 * @param int it
 * @param int Nz
 * @param int Nx
 * @param int Ny
 * @param int ealThick
 * @param int savew
 * @param int wsave
 * @return {*}
 * @Author: PingMin Cheung
 * @Date: 2024-05-23 20:39:57
 */
void elastic3dSnapshot(char *prefix, float ***vx, float ***vy, float ***vz, float ***sigma_xx, float ***sigma_yy,
                       float ***sigma_zz, float dt, int it, int Nz, int Nx, int Ny, int ealThick, int savew,
                       int wsave)
{
  char fileName[256];
  int nz, nx, ny;
  int iz, ix, iy;
  float ***p = NULL;
  nz = savew == 0 ? Nz - 2 * ealThick : Nz;
  nx = savew == 0 ? Nx - 2 * ealThick : Nx;
  ny = savew == 0 ? Ny - 2 * ealThick : Ny;
  switch (wsave)
  {
  case 1:
    sprintf(fileName, "%s_vx_%05d_%05d_%05d_%05dms.bin", prefix, nz, nx, ny, (int)(1000 * it * dt));
    write_snapshot3d(fileName, vx, Nz, Nx, Ny, ealThick, savew);
    sprintf(fileName, "%s_vy_%05d_%05d_%05d_%05dms.bin", prefix, nz, nx, ny, (int)(1000 * it * dt));
    write_snapshot3d(fileName, vy, Nz, Nx, Ny, ealThick, savew);
    sprintf(fileName, "%s_vz_%05d_%05d_%05d_%05dms.bin", prefix, nz, nx, ny, (int)(1000 * it * dt));
    write_snapshot3d(fileName, vz, Nz, Nx, Ny, ealThick, savew);
    break;
  case 2:
    p = alloc3float(Nz, Nx, Ny);
#pragma omp parallel for private(ix, iy, iz) firstprivate(savew, ealThick, Nx, Ny, Nz)
    for (iy = (savew == 0 ? ealThick : 0); iy < (savew == 0 ? Ny - ealThick : Ny); iy++)
    {
      for (ix = (savew == 0 ? ealThick : 0); ix < (savew == 0 ? Nx - ealThick : Nx); ix++)
      {
#pragma ivdep
        for (iz = (savew == 0 ? ealThick : 0); iz < (savew == 0 ? Nz - ealThick : Nz); iz++)
        {
          p[iy][ix][iz] = 0.333333333f * (sigma_xx[iy][ix][iz] + sigma_yy[iy][ix][iz] + sigma_zz[iy][ix][iz]);
        }
      }
    }
    sprintf(fileName, "%s_p_%05d_%05d_%05d_%05dms.bin", prefix, nz, nx, ny, (int)(1000 * it * dt));
    write_snapshot3d(fileName, p, Nz, Nx, Ny, ealThick, savew);
    free3float(p);
    break;
  case 3:
    elastic3dSnapshot(prefix, vx, vy, vz, sigma_xx, sigma_yy, sigma_zz, dt, it, Nz, Nx, Ny, ealThick, savew, 1);
    elastic3dSnapshot(prefix, vx, vy, vz, sigma_xx, sigma_yy, sigma_zz, dt, it, Nz, Nx, Ny, ealThick, savew, 2);
    break;
  }
}