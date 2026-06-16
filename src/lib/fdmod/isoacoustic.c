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

/***** receive wavefield *****/
/**
 *
 * @param {Receiver2D} receiver_p
 * @param {Receiver2D} receiver_vx
 * @param {Receiver2D} receiver_vz
 * @param float **p
 * @param float **vx
 * @param float **vz
 * @param int it
 * @param int mt
 * @param int ealThick
 * @param int seistype
 * @Author: PingMin Cheung
 * @Date: 2024-05-16 21:33:26
 */
void acoustic2dReceiver(Receiver2D receiver_p, Receiver2D receiver_vx, Receiver2D receiver_vz, float **p, float **vx,
                        float **vz, int it, int mt, int ealThick, int seistype)
{
  int itrace, irec;
  int igx, igz;
  if (it % mt == 0)
  {
    irec = it / mt;
    switch (seistype)
    {
    case 1:
#pragma omp parallel for private(itrace, igx, igz) firstprivate(irec, ealThick)
      for (itrace = 0; itrace < receiver_p.ntrace; itrace++)
      {
        igx = receiver_p.posx[itrace] + ealThick;
        igz = receiver_p.posz[itrace] + ealThick;
        receiver_p.recv[itrace][irec] = p[igx][igz];
      }
      break;
    case 2:
#pragma omp parallel for private(itrace, igx, igz) firstprivate(irec, ealThick)
      for (itrace = 0; itrace < receiver_vx.ntrace; itrace++)
      {
        igx = receiver_vx.posx[itrace] + ealThick;
        igz = receiver_vx.posz[itrace] + ealThick;
        receiver_vx.recv[itrace][irec] = vx[igx][igz];
      }
#pragma omp parallel for private(itrace, igx, igz) firstprivate(irec, ealThick)
      for (itrace = 0; itrace < receiver_vz.ntrace; itrace++)
      {
        igx = receiver_vz.posx[itrace] + ealThick;
        igz = receiver_vz.posz[itrace] + ealThick;
        receiver_vz.recv[itrace][irec] = vz[igx][igz];
      }
      break;
    case 3:
      acoustic2dReceiver(receiver_p, receiver_vx, receiver_vz, p, vx, vz, it, mt, ealThick, 1);
      acoustic2dReceiver(receiver_p, receiver_vx, receiver_vz, p, vx, vz, it, mt, ealThick, 2);
      break;
    }
  }
}
/**
 *
 * @param {Receiver3D} receiver_p
 * @param {Receiver3D} receiver_vx
 * @param {Receiver3D} receiver_vy
 * @param {Receiver3D} receiver_vz
 * @param float ***p
 * @param float ***vx
 * @param float ***vy
 * @param float ***vz
 * @param int it
 * @param int mt
 * @param int ealThick
 * @param int seistype
 * @return {*}
 * @Author: PingMin Cheung
 * @Date: 2024-05-17 14:30:57
 */
void acoustic3dReceiver(Receiver3D receiver_p, Receiver3D receiver_vx, Receiver3D receiver_vy, Receiver3D receiver_vz,
                        float ***p, float ***vx, float ***vy, float ***vz, int it, int mt, int ealThick, int seistype)
{
  int itrace, irec;
  int igx, igy, igz;
  if (it % mt == 0)
  {
    irec = it / mt;
    switch (seistype)
    {
    case 1:
#pragma omp parallel for private(itrace, igx, igy, igz) firstprivate(irec, ealThick)
      for (itrace = 0; itrace < receiver_p.ntrace; itrace++)
      {
        igy = receiver_p.posy[itrace] + ealThick;
        igx = receiver_p.posx[itrace] + ealThick;
        igz = receiver_p.posz[itrace] + ealThick;
        receiver_p.recv[itrace][irec] = p[igy][igx][igz];
      }
      break;
    case 2:
#pragma omp parallel for private(itrace, igx, igy, igz) firstprivate(irec, ealThick)
      for (itrace = 0; itrace < receiver_vy.ntrace; itrace++)
      {
        igy = receiver_vx.posy[itrace] + ealThick;
        igx = receiver_vx.posx[itrace] + ealThick;
        igz = receiver_vx.posz[itrace] + ealThick;
        receiver_vx.recv[itrace][irec] = vx[igy][igx][igz];
      }
#pragma omp parallel for private(itrace, igx, igy, igz) firstprivate(irec, ealThick)
      for (itrace = 0; itrace < receiver_vz.ntrace; itrace++)
      {
        igy = receiver_vy.posy[itrace] + ealThick;
        igx = receiver_vy.posx[itrace] + ealThick;
        igz = receiver_vy.posz[itrace] + ealThick;
        receiver_vy.recv[itrace][irec] = vy[igy][igx][igz];
      }
#pragma omp parallel for private(itrace, igx, igy, igz) firstprivate(irec, ealThick)
      for (itrace = 0; itrace < receiver_vz.ntrace; itrace++)
      {
        igy = receiver_vz.posy[itrace] + ealThick;
        igx = receiver_vz.posx[itrace] + ealThick;
        igz = receiver_vz.posz[itrace] + ealThick;
        receiver_vz.recv[itrace][irec] = vz[igy][igx][igz];
      }
      break;
    case 3:
      acoustic3dReceiver(receiver_p, receiver_vx, receiver_vy, receiver_vz, p, vx, vy, vz, it, mt, ealThick, 1);
      acoustic3dReceiver(receiver_p, receiver_vx, receiver_vy, receiver_vz, p, vx, vy, vz, it, mt, ealThick, 2);
      break;
    }
  }
}
/***** snapshot save *****/
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
void acoustic2dSnapshot(char *prefix, float **p, float **vx, float **vz, float dt, int it, int Nz, int Nx,
                        int ealThick, int savew, int wsave)
{
  char fileName[256];
  int nz, nx;
  nz = savew == 0 ? Nz - 2 * ealThick : Nz;
  nx = savew == 0 ? Nx - 2 * ealThick : Nx;
  switch (wsave)
  {
  case 1:
    sprintf(fileName, "%s_p_%05d_%05d_%05dms.bin", prefix, nz, nx, (int)(1000 * it * dt));
    write_snapshot2d(fileName, p, Nz, Nx, ealThick, savew);
    break;
  case 2:
    sprintf(fileName, "%s_vx_%05d_%05d_%05dms.bin", prefix, nz, nx, (int)(1000 * it * dt));
    write_snapshot2d(fileName, vx, Nz, Nx, ealThick, savew);
    sprintf(fileName, "%s_vz_%05d_%05d_%05dms.bin", prefix, nz, nx, (int)(1000 * it * dt));
    write_snapshot2d(fileName, vz, Nz, Nx, ealThick, savew);
    break;
  case 3:
    acoustic2dSnapshot(prefix, p, vx, vz, dt, it, Nz, Nx, ealThick, savew, 1);
    acoustic2dSnapshot(prefix, p, vx, vz, dt, it, Nz, Nx, ealThick, savew, 2);
    break;
  case 4:
    sprintf(fileName, "%s_vy_%05d_%05d_%05dms.bin", prefix, nz, nx, (int)(1000 * it * dt));
    write_snapshot2d(fileName, p, Nz, Nx, ealThick, savew);
    break;
  }
}
/**
 *
 * @param {char} *prefix
 * @param float ***p
 * @param float ***vx
 * @param float ***vy
 * @param float ***vz
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
 * @Date: 2024-05-17 15:43:08
 */
void acoustic3dSnapshot(char *prefix, float ***p, float ***vx, float ***vy, float ***vz, float dt, int it, int Nz,
                        int Nx, int Ny, int ealThick, int savew, int wsave)
{
  char fileName[256];
  int nz, nx, ny;
  nz = savew == 0 ? Nz - 2 * ealThick : Nz;
  nx = savew == 0 ? Nx - 2 * ealThick : Nx;
  ny = savew == 0 ? Ny - 2 * ealThick : Ny;
  switch (wsave)
  {
  case 1:
    sprintf(fileName, "%s_p_%05d_%05d_%05d_%05dms.bin", prefix, nz, nx, ny, (int)(1000 * it * dt));
    write_snapshot3d(fileName, p, Nz, Nx, Ny, ealThick, savew);
    break;
  case 2:
    sprintf(fileName, "%s_vx_%05d_%05d_%05d_%05dms.bin", prefix, nz, nx, ny, (int)(1000 * it * dt));
    write_snapshot3d(fileName, vx, Nz, Nx, Ny, ealThick, savew);
    sprintf(fileName, "%s_vy_%05d_%05d_%05d_%05dms.bin", prefix, nz, nx, ny, (int)(1000 * it * dt));
    write_snapshot3d(fileName, vy, Nz, Nx, Ny, ealThick, savew);
    sprintf(fileName, "%s_vz_%05d_%05d_%05d_%05dms.bin", prefix, nz, nx, ny, (int)(1000 * it * dt));
    write_snapshot3d(fileName, vz, Nz, Nx, Ny, ealThick, savew);
    break;
  case 3:
    acoustic3dSnapshot(prefix, p, vx, vy, vz, dt, it, Nz, Nx, Ny, ealThick, savew, 1);
    acoustic3dSnapshot(prefix, p, vx, vy, vz, dt, it, Nz, Nx, Ny, ealThick, savew, 2);
    break;
  }
}