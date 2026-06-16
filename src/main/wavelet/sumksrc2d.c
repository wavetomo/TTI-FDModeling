/*
 * @Description:
 * @Version: 1.0
 * @Autor: LaurentCheung
 * @Date: 2022-06-07 16:00:57
 * @LastEditTime: 2026-06-16 14:48:05
 */
/********************/
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
#include <wavelet.h>
/********************/
char *sdoc[] = {
    /* */
    " SUMKSRC2D - write source file                           ",
    " sumksrc2d sdepth sx dt nt wtype A fpeak delay > stdout  ",
    " Required parameter:                                     ",
    " sdepth  source position in 1st dimension (unit: cell)   ",
    " sx      source position in 2nd dimension (unit: cell)   ",
    " dt      sampling interval in time dimension             ",
    " nt      number of samples in time dimension             ",
    " Optional Parameters:                                    ",
    " wfile   a biniary file storeing wavelet                 ",
    " wtype   source type(default: 1):                        ",
    "         =1: Ricker wavelet                              ",
    "         =2: First derivative of Gaussian wavelet        ",
    "         =3: First derivative of Ricker wavelet          ",
    "         =4: Gaussian wavelet                            ",
    " A       peak amplitude of wavelet                       ",
    " fpeak   peak frequency of wavelet (unit: Hz)            ",
    " delay   source delay time (unit: s) can avoid the       ",
    "         first sample of the source be a big value       ",
    "         default: delay = 1.5 / fpeak                    ",
    " order   order of derivative (default = 0)               ",
    " shift   shift of wavelet(default = 0)                   ",
    NULL};
/**/
int main(int argc, char *argv[])
{
  /* Initialize */
  initargs(argc, argv);
  requestdoc(0);
  /* Get parameters */
  int nt, wtype, sdepth, sx;
  float dt, A, fpeak, delay;
  char *wfile = NULL;
  if (!getparstring("wfile", &wfile))
  {
    wfile = NULL;
    warn("<warning>: didn't specify a wavelet file! The program will use one of wavelet in the "
         "document!");
  }
  if (!getparint("nt", &nt))
  {
    err("<error>: must specify number of samples in time dimension!");
  }
  if (!getparint("wtype", &wtype))
  {
    wtype = 1;
    warn("<warning>: didn't specify type of wavelet! The program will set it %d!", wtype);
  }
  if (!getparint("sdepth", &sdepth))
  {
    err("<error>: must specify source position in 2nd dimension!");
  }
  if (!getparint("sx", &sx))
  {
    err("<error>: must specify source position in 1st dimension!");
  }
  if (!getparfloat("dt", &dt))
  {
    err("<error>: must specify sampleing interval in time dimension!");
  }
  if (!getparfloat("A", &A))
  {
    A = 1.0f;
    warn("<warning>: didn't specify peak amplitude of wavelet! The program will set it as %10.5f.", A);
  }
  if (!getparfloat("fpeak", &fpeak))
  {
    fpeak = 10.0f;
    warn("<warning>: didn't specify fpeak frequency of wavelet! The program will set it as %10.5f "
         "Hz.",
         fpeak);
  }
  if (!getparfloat("delay", &delay))
  {
    delay = 1.5f / fpeak;
    warn("<warning>: didn't specify delay time of wavelet! The program will set delay = 1.5 / "
         "fpeak =%10.5f Hz!",
         delay);
  }
  /**/
  FILE *fp = NULL;
  float *src = NULL;
  src = alloc1float(nt);
  if (wfile != NULL)
  {
    fp = efopen(wfile, "rb");
    if (fp == NULL)
    {
      err("<error>: failed to open wavelet file (%s)!", wfile);
    }
    else
    {
      efread(src, FSIZE, nt, fp);
    }
  }
  else
  {
    switch (wtype)
    {
    case 1:
      riker_wlt(src, A, fpeak, delay, dt, nt);
      break;
    case 2:
      fdgauss_wlt(src, A, fpeak, delay, dt, nt);
      break;
    case 3:
      driker_wlt(src, A, fpeak, delay, dt, nt);
      break;
    case 4:
      gauss_wlt(src, A, fpeak, delay, dt, nt);
      break;
    }
  }
  /**/
  segy tr;
  memset(&tr, 0, sizeof(segy));
  tr.tracl = 1;
  tr.trid = 1;
  tr.scalel = -1000;
  tr.scalco = -1000;
  tr.sdepth = head_scalar_write(sdepth, tr.scalco);
  tr.sx = head_scalar_write(sx, tr.scalco);
  tr.ns = nt;
  tr.dt = (dt * 1000000); /* micro-second */
  tr.f1 = 0.0f;
  tr.d1 = dt;
  tr.ntr = 1;
  memcpy(tr.data, src, tr.ns * FSIZE);
  puttr(&tr);
  /**/
  free1float(src);
  return (CWP_Exit());
}