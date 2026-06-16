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
/***** traces save *****/
/**
 *
 * @param {char} *fileName
 * @param {Source2D} source
 * @param {Receiver2D} receiver
 * @param float fx
 * @param float fz
 * @param float dx
 * @param float dz
 * @param float dt
 * @return {*}
 * @LastEditors: LaurentCheung
 * @Version: 1.0
 */
void putraces_2d(char *fileName, Source2D source, Receiver2D receiver, float fx, float fz, float dx, float dz,
                 float dt, int shotID)
{
  FILE *fp;
  segy tr;
  int itrace;
  fp = efopen(fileName, "wb");
  memset(&tr, 0, sizeof(segy));
  tr.tracl = 1;
  tr.trid = 1;
  tr.ns = receiver.ntsample;
  tr.dt = dt * 1000000; /* microsecond */
  tr.ntr = receiver.ntrace;
  tr.scalel = -1000;
  tr.scalco = -1000;
  tr.f1 = 0.0f;
  tr.d1 = dt;
  tr.fldr = shotID;
  /**/
  tr.sx = head_scalar_write(fx + source.posx[0] * dx, tr.scalco);
  tr.sdel = head_scalar_write(fz + source.posz[0] * dz, tr.scalel);
  for (itrace = 0; itrace < receiver.ntrace; itrace++)
  {
    tr.tracr = itrace + 1;
    tr.trid = receiver.posx[itrace];
    /**/
    tr.gx = head_scalar_write(fx + receiver.posx[itrace] * dx, tr.scalco);
    tr.gdel = head_scalar_write(fz + receiver.posz[itrace] * dz, tr.scalel);
    /**/
    tr.offset = tr.gx - tr.sx;
    tr.cdp = 0.5f * (tr.gx + tr.sx);
    /**/
    memcpy(tr.data, receiver.recv[itrace], tr.ns * FSIZE);
    fputtr(fp, &tr);
  }
  efclose(fp);
}
/**
 *
 * @param {char} *fileName
 * @param {Source2D} source
 * @param {Receiver2D} receiver
 * @param float fx
 * @param float fy
 * @param float fz
 * @param float dx
 * @param float dy
 * @param float dz
 * @param float dt
 * @return {*}
 * @LastEditors: LaurentCheung
 * @Version: 1.0
 */
void putraces_3d(char *fileName, Source3D source, Receiver3D receiver, float fx, float fy, float fz, float dx,
                 float dy, float dz, float dt, int shotID)
{
  FILE *fp;
  segy tr;
  int itrace;
  fp = efopen(fileName, "wb");
  memset(&tr, 0, sizeof(segy));
  tr.tracl = 1;
  tr.trid = 1;
  tr.ns = receiver.ntsample;
  tr.dt = dt * 1000000; /* microsecond */
  tr.ntr = receiver.ntrace;
  tr.scalel = -1000;
  tr.scalco = -1000;
  tr.f1 = 0.0f;
  tr.d1 = dt;
  tr.fldr = shotID;
  /**/
  tr.sx = head_scalar_write(fx + source.posx[0] * dx, tr.scalco);
  tr.sy = head_scalar_write(fy + source.posy[0] * dy, tr.scalco);
  tr.sdel = head_scalar_write(fz + source.posz[0] * dz, tr.scalel);
  for (itrace = 0; itrace < receiver.ntrace; itrace++)
  {
    tr.tracl = receiver.posy[itrace];
    tr.tracr = itrace + 1;
    tr.trid = itrace + 1;
    /**/
    tr.gx = head_scalar_write(fx + receiver.posx[itrace] * dx, tr.scalco);
    tr.gy = head_scalar_write(fy + receiver.posy[itrace] * dy, tr.scalco);
    tr.gdel = head_scalar_write(fz + receiver.posz[itrace] * dz, tr.scalel);
    /**/
    tr.offset = tr.gx - tr.sx;
    tr.cdp = 0.5f * (tr.gx + tr.sx);
    /**/
    memcpy(tr.data, receiver.recv[itrace], tr.ns * FSIZE);
    fputtr(fp, &tr);
  }
  efclose(fp);
}
/**
 *
 * @param {char} *trace_prefix
 * @param {Source2D} source
 * @param {Receiver2D} receiver
 * @param float fx
 * @param float fz
 * @param float dx
 * @param float dz
 * @param float dt
 * @param int mt
 * @param int shotID
 * @param int flag
 * @return {*}
 * @Author: PingMin Cheung
 * @Date: 2023-07-16 22:10:54
 */
void traces2d_output(char *trace_prefix, Source2D source, Receiver2D receiver, float fx, float fz, float dx, float dz,
                     float dt, int mt, int shotID, int flag)
{
  dt = mt * dt;
  char fileName[256];
  switch (flag)
  {
  case 1:
    sprintf(fileName, "%s_vx_%04d_%05d_%05d_%05.2fms.su", trace_prefix, shotID, receiver.ntsample, receiver.ntrace,
            dt * 1000.0f);
    break;
  case 2:
    sprintf(fileName, "%s_vy_%04d_%05d_%05d_%05.2fms.su", trace_prefix, shotID, receiver.ntsample, receiver.ntrace,
            dt * 1000.0f);
    break;
  case 3:
    sprintf(fileName, "%s_vz_%04d_%05d_%05d_%05.2fms.su", trace_prefix, shotID, receiver.ntsample, receiver.ntrace,
            dt * 1000.0f);
    break;
  case 4:
    sprintf(fileName, "%s_p_%04d_%05d_%05d_%05.2fms.su", trace_prefix, shotID, receiver.ntsample, receiver.ntrace,
            dt * 1000.0f);
    break;
  case 5:
    sprintf(fileName, "%s_vx_p_%04d_%05d_%05d_%05.2fms.su", trace_prefix, shotID, receiver.ntsample, receiver.ntrace,
            dt * 1000.0f);
    break;
  case 6:
    sprintf(fileName, "%s_vy_p_%04d_%05d_%05d_%05.2fms.su", trace_prefix, shotID, receiver.ntsample, receiver.ntrace,
            dt * 1000.0f);
    break;
  case 7:
    sprintf(fileName, "%s_vz_p_%04d_%05d_%05d_%05.2fms.su", trace_prefix, shotID, receiver.ntsample, receiver.ntrace,
            dt * 1000.0f);
    break;
  case 8:
    sprintf(fileName, "%s_vx_s_%04d_%05d_%05d_%05.2fms.su", trace_prefix, shotID, receiver.ntsample, receiver.ntrace,
            dt * 1000.0f);
    break;
  case 9:
    sprintf(fileName, "%s_vy_s_%04d_%05d_%05d_%05.2fms.su", trace_prefix, shotID, receiver.ntsample, receiver.ntrace,
            dt * 1000.0f);
    break;
  case 10:
    sprintf(fileName, "%s_vz_s_%04d_%05d_%05d_%05.2fms.su", trace_prefix, shotID, receiver.ntsample, receiver.ntrace,
            dt * 1000.0f);
    break;
  }
  /**/
  putraces_2d(fileName, source, receiver, fx, fz, dx, dz, dt, shotID);
}
/**
 *
 * @param {char} *trace_prefix:
 * @param {Source3D} source:
 * @param {Receiver3D} receiver:
 * @param float fx:
 * @param float fy:
 * @param float fz:
 * @param float dx:
 * @param float dy:
 * @param float dz:
 * @param float dt:
 * @param int mt:
 * @param int shotID:
 * @param int flag:
 * @LastEditors: LaurentCheung
 * @Version: 1.0
 */
void traces3d_output(char *trace_prefix, Source3D source, Receiver3D receiver, float fx, float fy, float fz, float dx,
                     float dy, float dz, float dt, int mt, int shotID, int flag)
{
  char fileName[256];
  dt = mt * dt;
  switch (flag)
  {
  case 1:
    sprintf(fileName, "%s_vx_%04d_%05d_%05d_%05.2fms.su", trace_prefix, shotID, receiver.ntsample, receiver.ntrace,
            dt * 1000.0f);
    break;
  case 2:
    sprintf(fileName, "%s_vy_%04d_%05d_%05d_%05.2fms.su", trace_prefix, shotID, receiver.ntsample, receiver.ntrace,
            dt * 1000.0f);
    break;
  case 3:
    sprintf(fileName, "%s_vz_%04d_%05d_%05d_%05.2fms.su", trace_prefix, shotID, receiver.ntsample, receiver.ntrace,
            dt * 1000.0f);
    break;
  case 4:
    sprintf(fileName, "%s_p_%04d_%05d_%05d_%05.2fms.su", trace_prefix, shotID, receiver.ntsample, receiver.ntrace,
            dt * 1000.0f);
    break;
  case 5:
    sprintf(fileName, "%s_vx_p_%04d_%05d_%05d_%05.2fms.su", trace_prefix, shotID, receiver.ntsample, receiver.ntrace,
            dt * 1000.0f);
    break;
  case 6:
    sprintf(fileName, "%s_vy_p_%04d_%05d_%05d_%05.2fms.su", trace_prefix, shotID, receiver.ntsample, receiver.ntrace,
            dt * 1000.0f);
    break;
  case 7:
    sprintf(fileName, "%s_vz_p_%04d_%05d_%05d_%05.2fms.su", trace_prefix, shotID, receiver.ntsample, receiver.ntrace,
            dt * 1000.0f);
    break;
  case 8:
    sprintf(fileName, "%s_vx_s_%04d_%05d_%05d_%05.2fms.su", trace_prefix, shotID, receiver.ntsample, receiver.ntrace,
            dt * 1000.0f);
    break;
  case 9:
    sprintf(fileName, "%s_vy_s_%04d_%05d_%05d_%05.2fms.su", trace_prefix, shotID, receiver.ntsample, receiver.ntrace,
            dt * 1000.0f);
    break;
  case 10:
    sprintf(fileName, "%s_vz_s_%04d_%05d_%05d_%05.2fms.su", trace_prefix, shotID, receiver.ntsample, receiver.ntrace,
            dt * 1000.0f);
    break;
  }
  /**/
  putraces_3d(fileName, source, receiver, fx, fy, fz, dx, dy, dz, dt, shotID);
}
/**
 *
 * @param {char} *trace_prefix:
 * @param {Source2D} source:
 * @param {Receiver2D} receiver:
 * @param float fx:
 * @param float fz:
 * @param float dx:
 * @param float dz:
 * @param float dt:
 * @param int shotID:
 * @param int flag:
 * @LastEditors: LaurentCheung
 * @Version: 1.0
 */
void TItraces2d_output(char *trace_prefix, Source2D source, Receiver2D receiver, float fx, float fz, float dx,
                       float dz, float dt, int shotID, int flag)
{
  char trace_suffix[256], fileName[256];
  switch (flag)
  {
  case 1:
    sprintf(trace_suffix, "_vx_%04d_%05d_%05d_%05.2fms.su", shotID, receiver.ntsample, receiver.ntrace, dt * 1000.0f);
    break;
  case 2:
    sprintf(trace_suffix, "_vz_%04d_%05d_%05d_%05.2fms.su", shotID, receiver.ntsample, receiver.ntrace, dt * 1000.0f);
    break;
  case 3:
    sprintf(trace_suffix, "_p_%04d_%05d_%05d_%05.2fms.su", shotID, receiver.ntsample, receiver.ntrace, dt * 1000.0f);
    break;
  case 4:
    sprintf(trace_suffix, "_q_%04d_%05d_%05d_%05.2fms.su", shotID, receiver.ntsample, receiver.ntrace, dt * 1000.0f);
    break;
  }
  /**/
  sprintf(fileName, "%s%s", trace_prefix, trace_suffix);
  putraces_2d(fileName, source, receiver, fx, fz, dx, dz, dt, shotID);
}
/**
 *
 * @param {char} *trace_prefix
 * @param {Source3D} source
 * @param {Receiver3D} receiver
 * @param float fx
 * @param float fy
 * @param float fz
 * @param float dx
 * @param float dy
 * @param float dz
 * @param float dt
 * @param int shotID
 * @param int flag
 * @return {*}
 * @LastEditors: LaurentCheung
 * @Version: 1.0
 */
void TItraces3d_output(char *trace_prefix, Source3D source, Receiver3D receiver, float fx, float fy, float fz,
                       float dx, float dy, float dz, float dt, int shotID, int flag)
{
  char trace_suffix[256], fileName[256];
  switch (flag)
  {
  case 1:
    sprintf(trace_suffix, "_vx_%04d_%05d_%05d_%05.2fms.su", shotID, receiver.ntsample, receiver.ntrace, dt * 1000.0f);
    break;
  case 2:
    sprintf(trace_suffix, "_vy_%04d_%05d_%05d_%05.2fms.su", shotID, receiver.ntsample, receiver.ntrace, dt * 1000.0f);
    break;
  case 3:
    sprintf(trace_suffix, "_vz_%04d_%05d_%05d_%05.2fms.su", shotID, receiver.ntsample, receiver.ntrace, dt * 1000.0f);
    break;
  case 4:
    sprintf(trace_suffix, "_p_%04d_%05d_%05d_%05.2fms.su", shotID, receiver.ntsample, receiver.ntrace, dt * 1000.0f);
    break;
  case 5:
    sprintf(trace_suffix, "_q_%04d_%05d_%05d_%05.2fms.su", shotID, receiver.ntsample, receiver.ntrace, dt * 1000.0f);
    break;
  }
  /**/
  sprintf(fileName, "%s%s", trace_prefix, trace_suffix);
  putraces_3d(fileName, source, receiver, fx, fy, fz, dx, dy, dz, dt, shotID);
}