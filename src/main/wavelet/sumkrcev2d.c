/*
 * @Description:
 * @Version: 1.0
 * @Autor: LaurentCheung
 * @Date: 2022-06-07 16:00:57
 * @LastEditTime: 2024-05-19 17:19:23
 */
/********************/
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/**/
#include <cwp.h>
#include <par.h>
#include <su.h>
/**/
#include <header.h>
#include <segy.h>
/**/
#include <supro.h>
/********************/
char *sdoc[] = {
    /**/
    " SUMKRCEV2D - write receiver position file                    ",
    " sumkrcev2d sx gdel offsetmin offsetmax dtr nx > stdout       ",
    " Required parameter:                                          ",
    " gdel       receiver position in 1st dimension (unit: cell)   ",
    " sx         source position in 2nd dimension (unit: cell)     ",
    " offsetmin  maximum offset (unit: cell)                       ",
    " offsetmax  minimum offset (unit: cell)                       ",
    " dtr        interval of trace (unit: cell)                    ",
    " nx         number of samples in 2nd dimension (unit: cell)   ",
    " mode       the mode of observation                           ",
    "            = 0: generate observation at both side of source  ",
    "            = 1: generate observation at left side of source  ",
    "            = 2: generate observation at right side of source ",
    NULL};
/**/
int main(int argc, char *argv[])
{
  /* Initialize */
  initargs(argc, argv);
  requestdoc(0);
  /* Get parameters */
  int sx, gdel, offsetmin, dtr, offsetmax, nx, mode;
  if (!getparint("sx", &sx))
  {
    err("<error>: didn't specify source position in 2nd dimension (unit: "
        "cell)!");
  }
  if (!getparint("gdel", &gdel))
  {
    gdel = 0;
    warn("<warning>: didn't specify receiver position in 1st dimension (unit: "
         "cell)! The program will set it as %d.",
         gdel);
  }
  if (!getparint("offsetmin", &offsetmin))
  {
    offsetmin = 0;
    warn("<warning>: didn't specify minimum offset (unit: cell). The program will "
         "set it as %d",
         offsetmin);
  }
  if (!getparint("dtr", &dtr))
  {
    dtr = 1;
    warn("<warning>: didn't specify interval of trace (unit: cell). The program "
         "will set it as %d.",
         dtr);
  }
  if (!getparint("offsetmax", &offsetmax))
  {
    err("<error>: didn't specify maximum offset(unit: cell).");
  }
  if (!getparint("nx", &nx))
  {
    err("<error>: didn't specify number of samples in 2nd dimension (unit: cell).");
  }
  if (!getparint("mode", &mode))
  {
    mode = 0;
    warn("<warning>: didn't specify the mode of observation. The program will set "
         "it as %d.",
         mode);
  }
  /**/
  int ntr_left, ntr_right, ntr;
  int gx_left, gx_right;
  ntr_left = 0;
  ntr_right = 0;
  if (mode == 1 || mode == 0)
  {
    gx_left = sx - offsetmin;
    while (gx_left >= 0 && sx - gx_left <= offsetmax)
    {
      ntr_left++;
      gx_left -= dtr;
    }
    gx_left += dtr;
  }
  if (mode == 2 || mode == 0)
  {
    gx_right = sx + offsetmin;
    while (gx_right < nx && gx_right - sx <= offsetmax)
    {
      ntr_right++;
      gx_right += dtr;
    }
    gx_right -= dtr;
  }
  ntr = ntr_left + ntr_right;
  if (mode == 0 && offsetmin == 0)
  {
    ntr--;
    ntr_left--;
  }
  warn("<warning>: the number of traces is %d!", ntr);
  if (ntr > SU_NFLTS)
  {
    err("<error>: the number of traces is too much.");
  }
  /**/
  segy tr;
  int itrace;
  memset(&tr, 0, sizeof(segy));
  tr.tracl = 1;
  tr.trid = 1;
  tr.scalel = -1000;
  tr.scalco = -1000;
  tr.gdel = head_scalar_write(gdel, tr.scalel);
  tr.sx = head_scalar_write(sx, tr.scalco);
  tr.ns = ntr;
  tr.ntr = 2;
  tr.dt = dtr * 10000000; /* micro-second */
  tr.f1 = 0.0f;
  tr.d1 = dtr;
  if (mode == 1 || mode == 0)
  {
    for (itrace = 0; itrace < ntr_left; itrace++)
    {
      tr.data[itrace] = gx_left + itrace * dtr;
    }
  }
  if (mode == 2 || mode == 0)
  {
    for (itrace = 0; itrace < ntr_right; itrace++)
    {
      tr.data[ntr - 1 - itrace] = gx_right - itrace * dtr;
    }
  }
  puttr(&tr);
  /**/
  for (itrace = 0; itrace < ntr; itrace++)
  {
    tr.data[itrace] = gdel;
  }
  puttr(&tr);
  /**/
  return (CWP_Exit());
}
