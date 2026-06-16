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
#include <su.h>
/**/
#include <header.h>
#include <segy.h>
/**/
#include <supro.h>
/********** **********/
/** observasion **/
/**
 * read source wavelet and position from file
 * @param Source2D *source: 2D source structure pointer
 * @param char *sourcefile: file name of source
 * @param float dt: sample interval in time dimension
 * @param int nx: number of samples in x-space dimension
 * @param int nz: number of samples in z-space dimension
 * @param int nt: number of samples in time dimension
 * @return Source2D *source: 2D source structure
 * @Author: PingMin Cheung
 * @Date: 2024-08-14 14:13:19
 */
void read2dsource(Source2D *source, char *sourcefile, float dt, int nx, int nz, int nt)
{
  FILE *fp = NULL;
  segy tr;
  int it, itrace;
  float *time = NULL;
  source2d_init(source);
  time = alloc1float(nt);
  for (it = 0; it < nt; it++)
  {
    time[it] = it * dt;
  }
  fp = efopen(sourcefile, "rb");
  if (fp == NULL)
  {
    err("<error>: failed to open sourcefile (%s)!", sourcefile);
  }
  else
  {
    fgettr(fp, &tr);
    source->ntrace = tr.ntr;
    source->ntsample = nt;
    source->posx = alloc1int(source->ntrace);
    source->posz = alloc1int(source->ntrace);
    source->src = alloc2float(source->ntsample, source->ntrace);
    itrace = -1;
    warn("<warning>: the number of trace is %d in sourcefile (%s)!", source->ntrace, sourcefile);
    do
    {
      itrace++;
      if (itrace >= source->ntrace)
      {
        err("<error>: the number of traces in sourcefile (%s) is out of %d!", sourcefile, source->ntrace);
      }
      source->posz[itrace] = head_scalar_read(tr.sdepth, tr.scalel);
      if (source->posz[itrace] < 0 || source->posz[itrace] > nz)
      {
        err("<error>: the position of source in 1st dimension (%d) must be "
            "between 0 and %d!",
            source->posz[itrace], nz);
      }
      source->posx[itrace] = head_scalar_read(tr.sx, tr.scalco);
      if (source->posx[itrace] < 0 || source->posx[itrace] > nx)
      {
        err("<error>: the position of source in 2nd dimension (%d) must be "
            "between 0 and %d!",
            source->posx[itrace], nx);
      }
      ints8r(tr.ns, tr.dt / 1000000.0f, tr.f1, tr.data, 0.0f, 0.0f, nt, time, source->src[itrace]);
    } while (fgettr(fp, &tr));
    efclose(fp);
  }
  free1float(time);
}
/**
 * read source wavelet and position from file
 * @param Source2D *source: 2D source structure pointer
 * @param char *sourcefile: file name of source
 * @param float dt: sample interval in time dimension
 * @param int nx: number of samples in x-space dimension
 * @param int ny: number of samples in y-space dimension
 * @param int nz: number of samples in z-space dimension
 * @param int nt: number of samples in time dimension
 * @return Source2D *source: 2D source structure
 * @Author: PingMin Cheung
 * @Date: 2024-08-14 14:13:19
 */
void read3dsource(Source3D *source, char *sourcefile, float dt, int nx, int ny, int nz, int nt)
{
  FILE *fp = NULL;
  segy tr;
  int it, itrace;
  float *time = NULL;
  source3d_init(source);
  time = alloc1float(nt);
  for (it = 0; it < nt; it++)
  {
    time[it] = it * dt;
  }
  fp = efopen(sourcefile, "rb");
  if (fp == NULL)
  {
    err("<error>: failed to open sourcefile (%s)!", sourcefile);
  }
  else
  {
    fgettr(fp, &tr);
    source->ntrace = tr.ntr;
    source->ntsample = nt;
    source->posx = alloc1int(source->ntrace);
    source->posy = alloc1int(source->ntrace);
    source->posz = alloc1int(source->ntrace);
    source->src = alloc2float(source->ntsample, source->ntrace);
    itrace = -1;
    warn("<warning>: the number of trace is %d in sourcefile (%s)!", source->ntrace, sourcefile);
    do
    {
      itrace++;
      if (itrace >= source->ntrace)
      {
        err("<error>: the number of traces in sourcefile (%s) is out of %d!", sourcefile, source->ntrace);
      }
      source->posz[itrace] = head_scalar_read(tr.sdepth, tr.scalel);
      if (source->posz[itrace] < 0 || source->posz[itrace] > nz)
      {
        err("<error>: the position of source in 1st dimension (%d) must be "
            "between 0 and %d!",
            source->posz[itrace], nz);
      }
      source->posx[itrace] = head_scalar_read(tr.sx, tr.scalco);
      if (source->posx[itrace] < 0 || source->posx[itrace] > nx)
      {
        err("<error>: the position of source in 2nd dimension (%d) must be "
            "between 0 and %d!",
            source->posx[itrace], nx);
      }
      source->posy[itrace] = head_scalar_read(tr.sy, tr.scalco);
      if (source->posy[itrace] < 0 || source->posy[itrace] > ny)
      {
        err("<error>: the position of source in 2nd dimension (%d) must be "
            "between 0 and %d!",
            source->posy[itrace], ny);
      }
      ints8r(tr.ns, tr.dt / 1000000.0f, tr.f1, tr.data, 0.0f, 0.0f, nt, time, source->src[itrace]);
    } while (fgettr(fp, &tr));
  }
  efclose(fp);
  free1float(time);
}
/**
 * read 2D receiver's information from file
 * @param Receiver2D *receiver: 2D receiver structure pointer
 * @param char *recposfile: file name of receiver
 * @param int nx: number of samples in x-space dimension
 * @param int nz: number of samples in z-space dimension
 * @param int nt: number of samples in time dimension
 * @param int mt: receiver's sample interval in time dimension
 * @return Receiver2D *receiver: 2D receiver structure pointer
 * @Author: PingMin Cheung
 * @Date: 2024-08-17 14:26:13
 */
void read2dreceiver(Receiver2D *receiver, char *recposfile, int nx, int nz, int nt, int mt)
{
  FILE *fp = NULL;
  segy trx, trz;
  int it, itrace;
  receiver2d_init(receiver);
  fp = fopen(recposfile, "rb");
  if (fp == NULL)
  {
    err("<error>: failed to open recposfile (%s)!", recposfile);
  }
  else
  {
    if (fvgettr(fp, &trx) && fvgettr(fp, &trz))
    {
      if (trx.ns != trz.ns)
      {
        err("<error>: the number of Group coordinate's samples are not match, Group "
            "X coordinate : %d, Group Z coordinate : %d.",
            trx.ns, trz.ns);
      }
      else
      {
        receiver->ntrace = trx.ns;
        receiver->ntsample = (nt - 1) / mt + 1;
        receiver->posx = alloc1int(receiver->ntrace);
        receiver->posz = alloc1int(receiver->ntrace);
        receiver->recv = alloc2float(receiver->ntsample, receiver->ntrace);
        for (itrace = 0; itrace < receiver->ntrace; itrace++)
        {
          receiver->posz[itrace] = trz.data[itrace];
          if (receiver->posz[itrace] < 0 || receiver->posz[itrace] > nz)
          {
            err("<error>: the position of receiver in 1st dimension (%d) must be "
                "between 0 and %d!",
                receiver->posz[itrace], nz);
          }
          receiver->posx[itrace] = trx.data[itrace];
          if (receiver->posx[itrace] < 0 || receiver->posx[itrace] > nx)
          {
            err("<error>: the position of receiver in 2nd dimension (%d) must be "
                "between 0 and %d!",
                receiver->posx[itrace], nx);
          }
        }
      }
      warn("the number of traces is %d!", receiver->ntrace);
      warn("the number of samples is %d in traces!", receiver->ntsample);
    }
    else
    {
      err("<error>: failed to get traces in recposfile (%s)!", recposfile);
    }
    efclose(fp);
  }
}
/**
 * read 3D receiver's information from file
 * @param Receiver3D *receiver: 3D receiver structure pointer
 * @param char *recposfile: file name of receiver
 * @param int nx: number of samples in x-space dimension
 * @param int ny: number of samples in y-space dimension
 * @param int nz: number of samples in z-space dimension
 * @param int nt: number of samples in time dimension
 * @param int mt: receiver's sample interval in time dimension
 * @return Receiver3D *receiver: 3D receiver structure pointer
 * @Author: PingMin Cheung
 * @Date: 2024-08-17 14:26:13
 */
void read3dreceiver(Receiver3D *receiver, char *recposfile, int nx, int ny, int nz, int nt, int mt)
{
  FILE *fp = NULL;
  segy trx, try, trz;
  int it, itrace;
  receiver3d_init(receiver);
  fp = fopen(recposfile, "rb");
  if (fp == NULL)
  {
    err("<error>: failed to open recposfile (%s)!", recposfile);
  }
  else
  {
    if (fvgettr(fp, &trx) && fvgettr(fp, &try) && fvgettr(fp, &trz))
    {
      if (trx.ns == try.ns && trx.ns == trz.ns)
      {
        receiver->ntrace = trx.ns;
        receiver->ntsample = (nt - 1) / mt + 1;
        receiver->posx = alloc1int(receiver->ntrace);
        receiver->posy = alloc1int(receiver->ntrace);
        receiver->posz = alloc1int(receiver->ntrace);
        receiver->recv = alloc2float(receiver->ntsample, receiver->ntrace);
        for (itrace = 0; itrace < receiver->ntrace; itrace++)
        {
          receiver->posx[itrace] = trx.data[itrace];
          if (receiver->posx[itrace] < 0 || receiver->posx[itrace] > nx)
          {
            err("<error>: the position of receiver in 2nd dimension (%d) must be "
                "between 0 and %d!",
                receiver->posx[itrace], nx);
          }
          receiver->posy[itrace] = try.data[itrace];
          if (receiver->posy[itrace] < 0 || receiver->posy[itrace] > ny)
          {
            err("<error>: the position of receiver in 3rd dimension (%d) must be "
                "between 0 and %d!",
                receiver->posy[itrace], ny);
          }
          receiver->posz[itrace] = trz.data[itrace];
          if (receiver->posz[itrace] < 0 || receiver->posz[itrace] > nz)
          {
            err("<error>: the position of receiver in 1st dimension (%d) must be "
                "between 0 and %d!",
                receiver->posz[itrace], nz);
          }
        }
      }
      else
      {
        err("<error>: the number of Group coordinate's samples are not match, Group "
            "X coordinate : %d, Group Y coordinate : %d, Group Z coordinate : %d.",
            trx.ns, try.ns, trz.ns);
      }
    }
    else
    {
      err("<error>: failed to get traces in recposfile (%s)!", recposfile);
    }
    efclose(fp);
    warn("the number of traces is %d!", receiver->ntrace);
    warn("the number of samples is %d in traces!", receiver->ntsample);
  }
}