/************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/* */
#include <omp.h>
/* */
#include <cwp.h>
#include <par.h>
#include <su.h>
/* */
#include <segy.h>
/* */
#include <absorbc.h>
#include <fdmod.h>
#include <rockmod.h>
/************************************/
char *sdoc[] = {
    /**/
    "                                                                                     ",
    " SUFDMOD2DTTIELASTIC.V1.0 nz= nx= nt= dx= dz= dt= fpeak= ...                         ",
    " Rotated staggered grid finite difference simulation of two-dimensional              ",
    " first-order elastic wave equation (rotated operator) in TTI media                   ",
    " illustration of stagger grid:                                                       ",
    "                                                                                     ",
    "  sigma_xx,sigma_zz,sigma_xz------------------------------sigma_xx,sigma_zz,sigma_xz ",
    "               |(ix,iz)                     |                          |(ix+1,iz)    ",
    "               |                            |                          |             ",
    "               |                            |                          |             ",
    "               |                            |                          |             ",
    "               |                            |                          |             ",
    "               |                            |                          |             ",
    "               |                            |                          |             ",
    "               |                            |                          |             ",
    "               ---------------------------vz,vx-------------------------             ",
    "               |                            |(ix,iz)                   |             ",
    "               |                            |                          |             ",
    "               |                            |                          |             ",
    "               |                            |                          |             ",
    "               |                            |                          |             ",
    "               |                            |                          |             ",
    "               |                            |                          |             ",
    "               |(ix,iz+1)                   |                          |(ix+1,iz+1)  ",
    "  sigma_xx,sigma_zz,sigma_xz------------------------------sigma_xx,sigma_zz,sigma_xz ",
    "                                                                                     ",
    " Required Parameters:                                                                ",
    " dz            sampling interval in 1st dimension.                                   ",
    " dx            sampling interval in 2nd dimension.                                   ",
    " dt            sampling interval in time dimension.                                  ",
    " nz            number of samples in 1st dimension.                                   ",
    " nx            number of samples in 2nd dimension.                                   ",
    " nt            number of samples in time dimension.                                  ",
    " vpfile        P-wave velocity file.                                                 ",
    " epsilonfile   Thomsen parameter (epsilon) file.                                     ",
    " deltafile     Thomsen parameter (delta) file.                                       ",
    " thetafile     dip angle file.                                                       ",
    " sourcefile    source file (format: su).                                             ",
    " Optional Parameters:                                                                ",
    " vsfile        S-wave velocity file.                                                 ",
    " rhofile       density file. default: rho =  310.0 * pow(vp, 0.25).                  ",
    " scalar        scalar to be applied to dip angle (default: 1.0).                     ",
    " fpeak         peak frequency of source wavelet.                                     ",
    " stype         which type of source will be load (default: 1).                       ",
    "               stype = 1: load P-wave source (divergence force which load on         ",
    "                          sigma_xx and sigma_zz).                                    ",
    "               stype = 2: load P-wave source (divergence force which load on         ",
    "                          vx and vz).                                                ",
    "               stype = 3: load S-wave source (curl force which load on               ",
    "                          vx and vz).                                                ",
    "               stype = 4: load equal energy P- and S-wave source (force which        ",
    "                          load on sigma_xx, sigma_zz and sigma_xz).                  ",
    "               stype = 5: load horizonal concentrated force source which load on vx. ",
    "               stype = 6: load vertical concentrated force source which load on vz.  ",
    "               stype = 7: load tilted concentrated force source which load on        ",
    "                          vx and vz and need to be improved.                         ",
    "               stype = 8: quasi curl force source which load on vx and vz.           ",
    " recposfile    receiver_p position file.                                             ",
    " seistype      which wavefiled will be recorded (default: 1).                        ",
    "               seistype = 0: don't record.                                           ",
    "               seistype = 1: record the particle velocity wavefiled                  ",
    "               seistype = 2: record the perssure.                                    ",
    "               seistype = 3: record the all wavefiled.                               ",
    " mt            sampling interval of seismic receiver (default: 1).                   ",
    " ealThick      thickness of effective absorbing layers (EAL), default: 30.           ",
    " eal_alpha     the order of gaussian weighting function in effective absorbing layer ",
    "               default: eal_alpha = 0.05 * ealThick + 1.5.                          ",
    " freesurface   whether to implement free surface boundary condition                  ",
    "               freesurface = 0: don't implement                                      ",
    "               freesurface = 1: implement free surface boundary condition by         ",
    "                                stress image method.                                 ",
    " mode          mode of modelng (default: 1).                                         ",
    "               mode = 1: execute the forward modeling.                               ",
    "               mode = 2: execute the backward modeling.                              ",
    " energy        whether to calculate the energy of wavefiled                          ",
    " savew         whether to save the wavefiled (default: 0).                           ",
    "               savew = 0: save the wavefiled without the absorbing                   ",
    "                          boundary area.                                             ",
    "               savew = 1: save the wavefiled with the absorbing boundary             ",
    "                          area.                                                      ",
    " wsave         which wavefiled will be saved the wavefiled (default: 1).             ",
    "               wsave = 0: don't save the wavefiled.                                  ",
    "               wsave = 1: save the particle velocity wavefiled.                      ",
    "               wsave = 2: save the pressure wavefiled.                               ",
    "               wsave = 3: save the all wavefiled.                                    ",
    " snap_sample   sampling interval of snapshot.                                        ",
    "               default: snap_sample =  500 / (dt * 1000)                             ",
    " fx            first sample in 2nd dimension.                                        ",
    " fz            first sample in 1st dimension.                                        ",
    " shotID        ID number of shot.                                                    ",
    " snap_prefix   prefix of snapshot file.                                              ",
    " trace_prefix  prefix of traces file.                                                ",
    NULL};
/******************************/
int main(int argc, char *argv[])
{
  /********* Hook up getpar *********/
  initargs(argc, argv);
  requestdoc(0);
  /* Start and end times only used for recording the cost*/
  double timeStart, timeEnd, timeCost, tStart, tEnd, tCost;
  timeStart = omp_get_wtime();
  /********* get required parameters *********/
  int nz, nx, nt;
  int freesurface, ealThick;
  int mode, savew, wsave, stype, seistype, mt, snap_sample, shotID;
  int energy;
  float dx, dz, dt, fx, fz, scalar, fpeak, eal_alpha;
  char *vpfile = NULL, *vsfile = NULL, *rhofile = NULL, *thetafile = NULL;
  char *epsilonfile = NULL, *deltafile = NULL;
  char *sourcefile = NULL, *recposfile = NULL;
  char *trace_prefix = NULL, *snap_prefix = NULL;
  /*** float ***/
  if (!getparfloat("dx", &dx))
  {
    err("<error>: must specify sampling interval in 1st dimension!");
  }
  if (!getparfloat("dz", &dz))
  {
    err("<error>: must specify sampling interval in 2nd dimension!");
  }
  if (!getparfloat("dt", &dt))
  {
    err("<error>: must specify sampling interval in time dimension!");
  }
  if (!getparfloat("fx", &fx))
  {
    fx = 0.0f;
    warn("<warning>: didn't specify first sample in 2nd dimension. The program will set it as %f.", fx);
  }
  if (!getparfloat("fz", &fz))
  {
    fz = 0.0f;
    warn("<warning>: didn't specify first sample in 1st dimension. The program will set it as %f.", fz);
  }
  if (!getparfloat("scalar", &scalar))
  {
    scalar = 1.0f;
    warn("<warning>: didn't specify the scalar to be apllied to dip angle. The program will set it as %f.", scalar);
  }
  if (!getparfloat("fpeak", &fpeak))
  {
    err("<error>: must specify peak frequency of source wavelet!");
  }
  if (!getparfloat("eal_alpha", &eal_alpha))
  {
    eal_alpha = -1.0f;
    warn("<warning>: didn't specify the order of gaussian weighting function! The program will set it as 0.05 * "
         "ealThick + 1.5.");
  }
  /*** int ***/
  if (!getparint("nz", &nz))
  {
    err("<error>: must specify the number of samples in 1st dimension!");
  }
  if (!getparint("nx", &nx))
  {
    err("<error>: must specify the number of samples in 2nd dimension!!");
  }
  if (!getparint("nt", &nt))
  {
    err("<error>: must specify the number of samples in time dimension!!");
  }
  if (!getparint("freesurface", &freesurface))
  {
    freesurface = 0;
    warn("<warning>: didn't specify whether to use free surface boundary conditions. The program will not use it.");
  }
  if (!getparint("ealThick", &ealThick))
  {
    ealThick = 30;
    warn("<warning>: didn't specify the number of absorbing layers! The program will set it %d.", ealThick);
  }
  if (!getparint("mode", &mode))
  {
    mode = 1;
    warn("<warning>: didn't specify the mode of modelng! The program will it as %d which means the program will "
         "execute the forward modeling.",
         mode);
  }
  if (!getparint("savew", &savew))
  {
    savew = 0;
    warn("<warning>: didn't specify whether to save the wavefiled. The program will set it as %d which means the "
         "program will not save the wavefiled of absoring boundary area.",
         savew);
  }
  if (!getparint("wsave", &wsave))
  {
    wsave = 1;
    warn("<warning>: didn't specify which wavefiled will be saved. The program will set it as %d which means the "
         "program will not save the wavefiled.",
         wsave);
  }
  if (!getparint("stype", &stype))
  {
    stype = 1;
    warn("<warning>: didn't specify the type of source,  and will set stype as %d!", stype);
  }
  if (!getparint("seistype", &seistype))
  {
    seistype = 0;
    warn("<warning>: didn't specify which wavefiled will be recorded. The program will set it as %d which means the "
         "program will not save the seismic record.",
         savew);
  }
  if (!getparint("mt", &mt))
  {
    mt = 1;
    warn("<warning>: didn't specify the sampling interval of seismic receiver_p. The program will set it as %d.", mt);
  }
  if (!getparint("snap_sample", &snap_sample))
  {
    snap_sample = 500 / (dt * 1000);
    warn("<warning>: didn't specify the sampling interval of snapshot. The program will set it as 500 / (dt * 1000) "
         "= %d.",
         snap_sample);
  }
  if (!getparint("shotID", &shotID))
  {
    shotID = 1;
    warn("<warning>: didn't specify the ID number of shot. The program will set it as %d.", shotID);
  }
  if (!getparint("energy", &energy))
  {
    energy = 0;
    warn("<warning>: didn't specify whether to calculate the energy of wavefiled. The program will set it as %d.",
         energy);
  }
  /*** char ***/
  if (!getparstring("vpfile", &vpfile))
  {
    err("<error>: must specify P-wave velocity!");
  }
  if (!getparstring("vsfile", &vsfile))
  {
    warn("<warning>: didn't specify S-wave velocity! The program will set S-wave velocity as vp / sqrt(3)!");
  }
  if (!getparstring("rhofile", &rhofile))
  {
    warn("<warning>: didn't specify density file! The program will set density as 310 * pow(vp, 0.25)!");
  }
  if (!getparstring("epsilonfile", &epsilonfile))
  {
    err("<error>: must specify Thomsen parameter (epsilon) file!");
  }
  if (!getparstring("deltafile", &deltafile))
  {
    err("<error>: must specify Thomsen parameter (delta) file!");
  }
  if (!getparstring("thetafile", &thetafile))
  {
    err("<error>: must specify dip angle (theta) file!");
  }
  if (!getparstring("sourcefile", &sourcefile))
  {
    err("<error>: must specify sourcefile!");
  }
  if (!getparstring("recposfile", &recposfile))
  {
    warn("<warning>: didn't specify receiver_p position file!");
  }
  if (!getparstring("snap_prefix", &snap_prefix))
  {
    snap_prefix = (char *)malloc(256 * sizeof(char));
    sprintf(snap_prefix, "snap");
    warn("<warning>: didn't specify prefix of snapshot file! The program will set it as %s.", snap_prefix);
  }
  if (!getparstring("trace_prefix", &trace_prefix))
  {
    trace_prefix = (char *)malloc(256 * sizeof(char));
    sprintf(trace_prefix, "traces");
    warn("<warning>: didn't specify the prefix of traces file. The program will set it as %s", trace_prefix);
  }
  /********* allocate memery *********/
  int N;
  int np;
  int Nz, Nx, Np;
  float **vp = NULL, **vs = NULL, **rho = NULL;
  float **epsilon = NULL, **delta = NULL, **theta = NULL;
  float **invrho = NULL, **c11 = NULL, **c13 = NULL, **c33 = NULL, **c55 = NULL;
  float **vx = NULL, **vz = NULL;
  float **sigma_xx = NULL, **sigma_zz = NULL, **sigma_xz = NULL;
  float **sin_theta = NULL, **cos_theta = NULL;
  float **sin_theta_v = NULL, **cos_theta_v = NULL;
  float **D = NULL, **Beta = NULL;
  float **coef1 = NULL, **coef2 = NULL;
  float *E = NULL, *Ek = NULL, *Ep = NULL;
  N = 5;
  ealThick += N;
  Nz = nz + 2 * ealThick;
  Nx = nx + 2 * ealThick;
  Np = Nz * Nx;
  np = nz * nx;
  vp = alloc2float(Nz, Nx);
  vs = alloc2float(Nz, Nx);
  rho = alloc2float(Nz, Nx);
  epsilon = alloc2float(Nz, Nx);
  delta = alloc2float(Nz, Nx);
  theta = alloc2float(Nz, Nx);
  invrho = alloc2float(Nz, Nx);
  c11 = alloc2float(Nz, Nx);
  c13 = alloc2float(Nz, Nx);
  c33 = alloc2float(Nz, Nx);
  c55 = alloc2float(Nz, Nx);
  vx = alloc2float(Nz, Nx);
  vz = alloc2float(Nz, Nx);
  sigma_xx = alloc2float(Nz, Nx);
  sigma_zz = alloc2float(Nz, Nx);
  sigma_xz = alloc2float(Nz, Nx);
  sin_theta = alloc2float(Nz, Nx);
  cos_theta = alloc2float(Nz, Nx);
  sin_theta_v = alloc2float(Nz, Nx);
  cos_theta_v = alloc2float(Nz, Nx);
  D = alloc2float(Nz, Nx);
  Beta = alloc2float(Nz, Nx);
  coef1 = alloc2float(Nz, Nx);
  coef2 = alloc2float(Nz, Nx);
  if (energy == 1)
  {
    E = alloc1float(nt);
    Ek = alloc1float(nt);
    Ep = alloc1float(nt);
  }
  /********* read and set velocity and density *********/
  FILE *fp = NULL;
  char fileName[256];
  int ix, iz, it;
  fp = fopen(vpfile, "rb");
  if (fp == NULL)
  {
    err("<error>: fail to open P-wave velocity file: %s!", vpfile);
  }
  else
  {
    for (ix = 0; ix < nx; ix++)
    {
      efread(vp[ix + ealThick] + ealThick, FSIZE, nz, fp);
    }
    efclose(fp);
  }
  if (vsfile == NULL)
  {
    vp2vs(vp[0], vs[0], Np);
    sprintf(fileName, "vs_%05d_%05d.bin", nz, nx);
    fp = fopen(fileName, "wb");
    for (ix = 0; ix < nx; ix++)
    {
      efwrite(vs[ix + ealThick] + ealThick, FSIZE, nz, fp);
    }
    efclose(fp);
  }
  else
  {
    fp = fopen(vsfile, "rb");
    if (fp == NULL)
    {
      err("<error>: fail to open S-wave velocity file: %s!", vsfile);
    }
    else
    {
      for (ix = 0; ix < nx; ix++)
      {
        efread(vs[ix + ealThick] + ealThick, FSIZE, nz, fp);
      }
      efclose(fp);
    }
  }
  if (rhofile == NULL)
  {
    vp2densityGardner(vp[0], rho[0], Np);
    sprintf(fileName, "rho_%05d_%05d.bin", nz, nx);
    fp = fopen(fileName, "wb");
    for (ix = 0; ix < nx; ix++)
    {
      efwrite(rho[ix + ealThick] + ealThick, FSIZE, nz, fp);
    }
    efclose(fp);
  }
  else
  {
    fp = fopen(rhofile, "rb");
    if (fp == NULL)
    {
      err("<error>: fail to open density file: %s!", rhofile);
    }
    else
    {
      for (ix = 0; ix < nx; ix++)
      {
        efread(rho[ix + ealThick] + ealThick, FSIZE, nz, fp);
      }
      efclose(fp);
    }
  }
  fp = fopen(epsilonfile, "rb");
  if (fp == NULL)
  {
    err("<error>: fail to open Thomsen parameter (epsilon) file: %s!", epsilonfile);
  }
  else
  {
    for (ix = 0; ix < nx; ix++)
    {
      efread(epsilon[ix + ealThick] + ealThick, FSIZE, nz, fp);
    }
    efclose(fp);
  }
  fp = fopen(deltafile, "rb");
  if (fp == NULL)
  {
    err("<error>: fail to open Thomsen parameter (delta) file: %s!", deltafile);
  }
  else
  {
    for (ix = 0; ix < nx; ix++)
    {
      efread(delta[ix + ealThick] + ealThick, FSIZE, nz, fp);
    }
    efclose(fp);
  }
  fp = fopen(thetafile, "rb");
  if (fp == NULL)
  {
    err("<error>: fail to open dip angle file: %s!", thetafile);
  }
  else
  {
    for (ix = 0; ix < nx; ix++)
    {
      efread(theta[ix + ealThick] + ealThick, FSIZE, nz, fp);
    }
    efclose(fp);
  }
  /********* padding modeling *********/
  padModel2d(vp, Nx, Nz, ealThick);
  padModel2d(vs, Nx, Nz, ealThick);
  padModel2d(rho, Nx, Nz, ealThick);
  padModel2d(epsilon, Nx, Nz, ealThick);
  padModel2d(delta, Nx, Nz, ealThick);
  padModel2d(theta, Nx, Nz, ealThick);
  /********* FD coefficients *********/
  float C1, C2, C3, C4, C5;
  /* TE FD coef1 */
  // float C[5] = {1.21124267580382f, -0.0897216796891598f, 0.0138427734377559f,
  // -0.00176565987726478f, 0.000118679470488305f};
  /* LS FD coef1 */
  float C[5] = {1.24362739770907f, -0.114021164233934f, 0.0272965993810630f, -0.00699199493927320f,
                0.00132022981661381f};
  C1 = C[0];
  C2 = C[1];
  C3 = C[2];
  C4 = C[3];
  C5 = C[4];
  /********* stability analysis *********/
  staggered_grid_stable2d(vp, C, dx, dz, dt, Nx, Nz, N);
  /********* MEAL *********/
  mealCoefs2d(D, Beta, vp, vs, fpeak, eal_alpha, dz, dx, Nz, Nx, ealThick, freesurface);
  mealFDCoefs(coef1[0], coef2[0], NULL, D[0], Beta[0], vp[0], vs[0], rho[0], dt, Np, 1);
  sprintf(fileName, "D_%05d_%05d.bin", Nz, Nx);
  fp = fopen(fileName, "wb");
  efwrite(D[0], FSIZE, Np, fp);
  efclose(fp);
  sprintf(fileName, "Beta_%05d_%05d.bin", Nz, Nx);
  fp = fopen(fileName, "wb");
  efwrite(Beta[0], FSIZE, Np, fp);
  efclose(fp);
  sprintf(fileName, "coef1_%05d_%05d.bin", Nz, Nx);
  fp = fopen(fileName, "wb");
  efwrite(coef1[0], FSIZE, Np, fp);
  efclose(fp);
  sprintf(fileName, "coef2_%05d_%05d.bin", Nz, Nx);
  fp = fopen(fileName, "wb");
  efwrite(coef2[0], FSIZE, Np, fp);
  efclose(fp);
  /********* model parameters *********/
  double f;
#pragma omp parallel for private(ix, iz, f)
  for (ix = 0; ix < Nx; ix++)
  {
#pragma ivdep
    for (iz = 0; iz < Nz; iz++)
    {
      f = 1.0 - pow(vs[ix][iz] / vp[ix][iz], 2.0);
      c11[ix][iz] = rho[ix][iz] * vp[ix][iz] * vp[ix][iz] * (1.0 + 2.0 * epsilon[ix][iz]);
      c13[ix][iz] =
          rho[ix][iz] * (vp[ix][iz] * vp[ix][iz] * sqrt(f * (f + 2.0 * delta[ix][iz])) - vs[ix][iz] * vs[ix][iz]);
      c33[ix][iz] = rho[ix][iz] * vp[ix][iz] * vp[ix][iz];
      c55[ix][iz] = rho[ix][iz] * vs[ix][iz] * vs[ix][iz];

      theta[ix][iz] = scalar * theta[ix][iz] / 180.0 * PI;
      sin_theta[ix][iz] = sin(theta[ix][iz]);
      cos_theta[ix][iz] = cos(theta[ix][iz]);
    }
  }
#pragma omp parallel for private(ix, iz)
  for (ix = 0; ix < Nx; ix++)
  {
#pragma ivdep
    for (iz = 0; iz < Nz; iz++)
    {
      if (ix < Nx - 1 && iz < Nz - 1)
      {
        sin_theta_v[ix][iz] =
            0.25 * (sin_theta[ix][iz] + sin_theta[ix + 1][iz] + sin_theta[ix][iz + 1] + sin_theta[ix + 1][iz + 1]);
        cos_theta_v[ix][iz] =
            0.25 * (cos_theta[ix][iz] + cos_theta[ix + 1][iz] + cos_theta[ix][iz + 1] + cos_theta[ix + 1][iz + 1]);
        invrho[ix][iz] = 4.0 / (rho[ix][iz] + rho[ix + 1][iz] + rho[ix][iz + 1] + rho[ix + 1][iz + 1]);
      }
      else
      {
        sin_theta_v[ix][iz] = sin_theta[ix][iz];
        cos_theta_v[ix][iz] = cos_theta[ix][iz];
        invrho[ix][iz] = 1.0 / rho[ix][iz];
      }
    }
  }
  sprintf(fileName, "invrho_%05d_%05d.bin", Nz, Nx);
  fp = fopen(fileName, "wb");
  efwrite(invrho[0], FSIZE, Np, fp);
  efclose(fp);
  sprintf(fileName, "c11_%05d_%05d.bin", Nz, Nx);
  fp = fopen(fileName, "wb");
  efwrite(c11[0], FSIZE, Np, fp);
  efclose(fp);
  sprintf(fileName, "c13_%05d_%05d.bin", Nz, Nx);
  fp = fopen(fileName, "wb");
  efwrite(c13[0], FSIZE, Np, fp);
  efclose(fp);
  sprintf(fileName, "c33_%05d_%05d.bin", Nz, Nx);
  fp = fopen(fileName, "wb");
  efwrite(c33[0], FSIZE, Np, fp);
  efclose(fp);
  sprintf(fileName, "c55_%05d_%05d.bin", Nz, Nx);
  fp = fopen(fileName, "wb");
  efwrite(c55[0], FSIZE, Np, fp);
  efclose(fp);
  /********* Source and Reciver *********/
  Source2D source;
  Receiver2D receiver_p, receiver_vx, receiver_vz;
  receiver2d_init(&receiver_p);
  receiver2d_init(&receiver_vx);
  receiver2d_init(&receiver_vz);
  read2dsource(&source, sourcefile, dt, nx, nz, nt);
  if (seistype == 1 || seistype == 3)
  {
    read2dreceiver(&receiver_vx, recposfile, nx, nz, nt, mt);
    read2dreceiver(&receiver_vz, recposfile, nx, nz, nt, mt);
  }
  if (seistype == 2 || seistype == 3)
  {
    read2dreceiver(&receiver_p, recposfile, nx, nz, nt, mt);
  }
  /********* show parameter *********/
  fprintf(stderr, "%s: Input Parameters:\n", argv[0]);
  fprintf(stderr, "%s: nx = %d, nz = %d, nt = %d\n", argv[0], nx, nz, nt);
  fprintf(stderr, "%s: Nx = %d, Nz = %d\n", argv[0], Nx, Nz);
  fprintf(stderr, "%s: dx = %6.5f, dz = %6.5f, dt = %6.5f\n", argv[0], dx, dz, dt);
  /********* initial these arrays with zero *********/
  memset(vx[0], 0, Np * FSIZE);
  memset(vz[0], 0, Np * FSIZE);
  memset(sigma_xx[0], 0, Np * FSIZE);
  memset(sigma_zz[0], 0, Np * FSIZE);
  memset(sigma_xz[0], 0, Np * FSIZE);
  if (energy == 1)
  {
    memset(E, 0, nt * FSIZE);
    memset(Ek, 0, nt * FSIZE);
    memset(Ep, 0, nt * FSIZE);
  }
  /********* FD modeling *********/
  int irec, isrc, isx, isz;
  int ixx, izz;
  int ixStart, ixEnd, izStart, izEnd;
  float d1, d2;
  float dvxdx, dvxdz, dvzdx, dvzdz;
  float dsigma_xxdx, dsigma_xxdz, dsigma_zzdx, dsigma_zzdz, dsigma_xzdx, dsigma_xzdz;
  float invdx, invdz;
  float alpha, weight_src, rx, rz;
  float ek, ep;
  ixStart = N;
  ixEnd = Nx - N;
  izStart = (freesurface == 0 ? N : ealThick);
  izEnd = Nz - N;
  invdx = 1.0 / dx;
  invdz = 1.0 / dz;
  alpha = logf(0.8f) / (dx * dx + dz * dz);
  fprintf(stderr, "%s: ixStart = %d, ixEnd = %d, izStart = %d, izEnd = %d\n", argv[0], ixStart, ixEnd, izStart,
          izEnd);
  tStart = omp_get_wtime();
  for (it = (mode == 1 ? 0 : nt - 1); (mode == 1 ? it < nt : it > -1); (mode == 1 ? it++ : it--))
  {
    /********* finite difference scheme *********/
#pragma omp parallel for private(ix, iz, dsigma_xxdx, dsigma_xxdz, dsigma_zzdx, dsigma_zzdz, dsigma_xzdx,            \
                                     dsigma_xzdz, d1, d2) firstprivate(invdx, invdz, C1, C2, C3, C4, C5)
    for (ix = ixStart; ix < ixEnd; ix++)
    {
#pragma ivdep
      for (iz = izStart; iz < izEnd; iz++)
      {
        d1 = C1 * (sigma_xx[ix + 1][iz + 1] - sigma_xx[ix][iz]) +
             C2 * (sigma_xx[ix + 2][iz + 2] - sigma_xx[ix - 1][iz - 1]) +
             C3 * (sigma_xx[ix + 3][iz + 3] - sigma_xx[ix - 2][iz - 2]) +
             C4 * (sigma_xx[ix + 4][iz + 4] - sigma_xx[ix - 3][iz - 3]) +
             C5 * (sigma_xx[ix + 5][iz + 5] - sigma_xx[ix - 4][iz - 4]);
        d2 = C1 * (sigma_xx[ix + 1][iz] - sigma_xx[ix][iz + 1]) +
             C2 * (sigma_xx[ix + 2][iz - 1] - sigma_xx[ix - 1][iz + 2]) +
             C3 * (sigma_xx[ix + 3][iz - 2] - sigma_xx[ix - 2][iz + 3]) +
             C4 * (sigma_xx[ix + 4][iz - 3] - sigma_xx[ix - 3][iz + 4]) +
             C5 * (sigma_xx[ix + 5][iz - 4] - sigma_xx[ix - 4][iz + 5]);

        dsigma_xxdx = 0.5f * (d1 + d2) * invdx;
        dsigma_xxdz = 0.5f * (d1 - d2) * invdx;

        d1 = C1 * (sigma_zz[ix + 1][iz + 1] - sigma_zz[ix][iz]) +
             C2 * (sigma_zz[ix + 2][iz + 2] - sigma_zz[ix - 1][iz - 1]) +
             C3 * (sigma_zz[ix + 3][iz + 3] - sigma_zz[ix - 2][iz - 2]) +
             C4 * (sigma_zz[ix + 4][iz + 4] - sigma_zz[ix - 3][iz - 3]) +
             C5 * (sigma_zz[ix + 5][iz + 5] - sigma_zz[ix - 4][iz - 4]);
        d2 = C1 * (sigma_zz[ix + 1][iz] - sigma_zz[ix][iz + 1]) +
             C2 * (sigma_zz[ix + 2][iz - 1] - sigma_zz[ix - 1][iz + 2]) +
             C3 * (sigma_zz[ix + 3][iz - 2] - sigma_zz[ix - 2][iz + 3]) +
             C4 * (sigma_zz[ix + 4][iz - 3] - sigma_zz[ix - 3][iz + 4]) +
             C5 * (sigma_zz[ix + 5][iz - 4] - sigma_zz[ix - 4][iz + 5]);

        dsigma_zzdx = 0.5f * (d1 + d2) * invdx;
        dsigma_zzdz = 0.5f * (d1 - d2) * invdz;

        d1 = C1 * (sigma_xz[ix + 1][iz + 1] - sigma_xz[ix][iz]) +
             C2 * (sigma_xz[ix + 2][iz + 2] - sigma_xz[ix - 1][iz - 1]) +
             C3 * (sigma_xz[ix + 3][iz + 3] - sigma_xz[ix - 2][iz - 2]) +
             C4 * (sigma_xz[ix + 4][iz + 4] - sigma_xz[ix - 3][iz - 3]) +
             C5 * (sigma_xz[ix + 5][iz + 5] - sigma_xz[ix - 4][iz - 4]);
        d2 = C1 * (sigma_xz[ix + 1][iz] - sigma_xz[ix][iz + 1]) +
             C2 * (sigma_xz[ix + 2][iz - 1] - sigma_xz[ix - 1][iz + 2]) +
             C3 * (sigma_xz[ix + 3][iz - 2] - sigma_xz[ix - 2][iz + 3]) +
             C4 * (sigma_xz[ix + 4][iz - 3] - sigma_xz[ix - 3][iz + 4]) +
             C5 * (sigma_xz[ix + 5][iz - 4] - sigma_xz[ix - 4][iz + 5]);

        dsigma_xzdx = 0.5f * (d1 + d2) * invdx;
        dsigma_xzdz = 0.5f * (d1 - d2) * invdz;

        vx[ix][iz] =
            coef1[ix][iz] * vx[ix][iz] + coef2[ix][iz] * invrho[ix][iz] *
                                             (cos_theta_v[ix][iz] * dsigma_xxdx - sin_theta_v[ix][iz] * dsigma_xxdz +
                                              sin_theta_v[ix][iz] * dsigma_xzdx + cos_theta_v[ix][iz] * dsigma_xzdz);
        vz[ix][iz] =
            coef1[ix][iz] * vz[ix][iz] + coef2[ix][iz] * invrho[ix][iz] *
                                             (cos_theta_v[ix][iz] * dsigma_xzdx - sin_theta_v[ix][iz] * dsigma_xzdz +
                                              sin_theta_v[ix][iz] * dsigma_zzdx + cos_theta_v[ix][iz] * dsigma_zzdz);
      }
    }
    /********* Load Source *********/
    if (stype == 2 || stype == 4)
    { /* P-wave source (divergence force which load on vx and vz) */
      for (isrc = 0; isrc < source.ntrace; isrc++)
      {
        isx = source.posx[isrc] + ealThick;
        isz = source.posz[isrc] + ealThick;
#pragma omp parallel for private(ix, iz, weight_src, rx, rz)                                                         \
    firstprivate(dt, dx, dz, invdx, invdz, alpha, isx, isz, N)
        for (ix = MAX(isx - N, 0); ix <= MIN(isx + N, Nx - 1); ix++)
        {
          rx = (ix - isx) * dx;
#pragma ivdep
          for (iz = MAX(isz - N, 0); iz <= MIN(isz + N, Nz - 1); iz++)
          {
            rz = (iz - isz) * dz;
            weight_src = expf(alpha * (rx * rx + rz * rz)) * invdx * invdz * dt;
            vx[ix][iz] += 2.0f * alpha * rx * weight_src * source.src[isrc][it];
            vz[ix][iz] += 2.0f * alpha * rz * weight_src * source.src[isrc][it];
          }
        }
      }
    }
    if (stype == 3 || stype == 4)
    { /* S-wave source (curl force source which load on vx and vz) */
      for (isrc = 0; isrc < source.ntrace; isrc++)
      {
        isx = source.posx[isrc] + ealThick;
        isz = source.posz[isrc] + ealThick;
#pragma omp parallel for private(ix, iz, weight_src, rx, rz)                                                         \
    firstprivate(dt, dx, dz, invdx, invdz, alpha, isx, isz, N)
        for (ix = MAX(isx - N, 0); ix <= MIN(isx + N, Nx - 1); ix++)
        {
          rx = (ix - isx) * dx;
#pragma ivdep
          for (iz = MAX(isz - N, 0); iz <= MIN(isz + N, Nz - 1); iz++)
          {
            rz = (iz - isz) * dz;
            weight_src = expf(alpha * (rx * rx + rz * rz)) * invdx * invdz * dt;
            vx[ix][iz] -= 2.0f * alpha * rz * weight_src * source.src[isrc][it];
            vz[ix][iz] += 2.0f * alpha * rx * weight_src * source.src[isrc][it];
          }
        }
      }
    }
    if (stype == 5)
    { /* horizonal concentrated force source which load on vx */
      for (isrc = 0; isrc < source.ntrace; isrc++)
      {
        isx = source.posx[isrc] + ealThick;
        isz = source.posz[isrc] + ealThick;
#pragma omp parallel for private(ix, iz, weight_src, rx, rz)                                                         \
    firstprivate(dt, dx, dz, invdx, invdz, alpha, isx, isz, N)
        for (ix = MAX(isx - N, 0); ix <= MIN(isx + N, Nx - 1); ix++)
        {
          rx = (ix - isx) * dx;
#pragma ivdep
          for (iz = MAX(isz - N, 0); iz <= MIN(isz + N, Nz - 1); iz++)
          {
            rz = (iz - isz) * dz;
            weight_src = expf(alpha * (rx * rx + rz * rz)) * invdx * invdz * dt;
            vx[ix][iz] += weight_src * source.src[isrc][it];
          }
        }
      }
    }
    if (stype == 6)
    { /* vertical concentrated force source which load on vz */
      for (isrc = 0; isrc < source.ntrace; isrc++)
      {
        isx = source.posx[isrc] + ealThick;
        isz = source.posz[isrc] + ealThick;
#pragma omp parallel for private(ix, iz, weight_src, rx, rz)                                                         \
    firstprivate(dt, dx, dz, invdx, invdz, alpha, isx, isz, N)
        for (ix = MAX(isx - N, 0); ix <= MIN(isx + N, Nx - 1); ix++)
        {
          rx = (ix - isx) * dx;
#pragma ivdep
          for (iz = MAX(isz - N, 0); iz <= MIN(isz + N, Nz - 1); iz++)
          {
            rz = (iz - isz) * dz;
            weight_src = expf(alpha * (rx * rx + rz * rz)) * invdx * invdz * dt;
            vz[ix][iz] += weight_src * source.src[isrc][it];
          }
        }
      }
    }
    if (stype == 7)
    { /* tilted concentrated force source which load on vx and vz */
      for (isrc = 0; isrc < source.ntrace; isrc++)
      {
        isx = source.posx[isrc] + ealThick;
        isz = source.posz[isrc] + ealThick;
#pragma omp parallel for private(ix, iz, weight_src, rx, rz)                                                         \
    firstprivate(dt, dx, dz, invdx, invdz, alpha, isx, isz, N)
        for (ix = MAX(isx - N, 0); ix <= MIN(isx + N, Nx - 1); ix++)
        {
          rx = (ix - isx) * dx;
#pragma ivdep
          for (iz = MAX(isz - N, 0); iz <= MIN(isz + N, Nz - 1); iz++)
          {
            rz = (iz - isz) * dz;
            weight_src = expf(alpha * (rx * rx + rz * rz)) * invdx * invdz * dt;
            vx[ix][iz] += 0.5f * sqrtf(2.0f) * weight_src * source.src[isrc][it];
            vz[ix][iz] += 0.5f * sqrtf(2.0f) * weight_src * source.src[isrc][it];
          }
        }
      }
    }
    if (stype == 8)
    { /* quasi curl force source which load on vx and vz and need to be improved */
#pragma omp parallel for private(isrc, isx, isz) firstprivate(dt, invdx, invdz, ealThick)
      for (isrc = 0; isrc < source.ntrace; isrc++)
      {
        isx = source.posx[isrc] + ealThick;
        isz = source.posz[isrc] + ealThick;

        vx[isx][isz] -= 0.125f * sqrtf(2.0f) * source.src[isrc][it];
        vz[isx][isz] += 0.125f * sqrtf(2.0f) * source.src[isrc][it];

        vx[isx - 1][isz] -= 0.125f * sqrtf(2.0f) * source.src[isrc][it];
        vz[isx - 1][isz] -= 0.125f * sqrtf(2.0f) * source.src[isrc][it];

        vx[isx - 1][isz - 1] += 0.125f * sqrtf(2.0f) * source.src[isrc][it];
        vz[isx - 1][isz - 1] -= 0.125f * sqrtf(2.0f) * source.src[isrc][it];

        vx[isx][isz - 1] += 0.125f * sqrtf(2.0f) * source.src[isrc][it];
        vz[isx][isz - 1] += 0.125f * sqrtf(2.0f) * source.src[isrc][it];
      }
    }
    /********* finite difference scheme *********/
#pragma omp parallel for private(ix, iz, dvxdx, dvxdz, dvzdx, dvzdz, d1, d2)                                         \
    firstprivate(invdx, invdz, C1, C2, C3, C4, C5)
    for (ix = ixStart; ix < ixEnd; ix++)
    {
#pragma ivdep
      for (iz = izStart; iz < izEnd; iz++)
      {
        d1 = C1 * (vx[ix][iz] - vx[ix - 1][iz - 1]) + C2 * (vx[ix + 1][iz + 1] - vx[ix - 2][iz - 2]) +
             C3 * (vx[ix + 2][iz + 2] - vx[ix - 3][iz - 3]) + C4 * (vx[ix + 3][iz + 3] - vx[ix - 4][iz - 4]) +
             C5 * (vx[ix + 4][iz + 4] - vx[ix - 5][iz - 5]);
        d2 = C1 * (vx[ix][iz - 1] - vx[ix - 1][iz]) + C2 * (vx[ix + 1][iz - 2] - vx[ix - 2][iz + 1]) +
             C3 * (vx[ix + 2][iz - 3] - vx[ix - 3][iz + 2]) + C4 * (vx[ix + 3][iz - 4] - vx[ix - 4][iz + 3]) +
             C5 * (vx[ix + 4][iz - 5] - vx[ix - 5][iz + 4]);

        dvxdx = 0.5f * (d1 + d2) * invdx;
        dvxdz = 0.5f * (d1 - d2) * invdz;

        d1 = sin_theta[ix][iz] * dvxdx + cos_theta[ix][iz] * dvxdz;
        dvxdx = cos_theta[ix][iz] * dvxdx - sin_theta[ix][iz] * dvxdz;
        dvxdz = d1;

        d1 = C1 * (vz[ix][iz] - vz[ix - 1][iz - 1]) + C2 * (vz[ix + 1][iz + 1] - vz[ix - 2][iz - 2]) +
             C3 * (vz[ix + 2][iz + 2] - vz[ix - 3][iz - 3]) + C4 * (vz[ix + 3][iz + 3] - vz[ix - 4][iz - 4]) +
             C5 * (vz[ix + 4][iz + 4] - vz[ix - 5][iz - 5]);
        d2 = C1 * (vz[ix][iz - 1] - vz[ix - 1][iz]) + C2 * (vz[ix + 1][iz - 2] - vz[ix - 2][iz + 1]) +
             C3 * (vz[ix + 2][iz - 3] - vz[ix - 3][iz + 2]) + C4 * (vz[ix + 3][iz - 4] - vz[ix - 4][iz + 3]) +
             C5 * (vz[ix + 4][iz - 5] - vz[ix - 5][iz + 4]);

        dvzdx = 0.5f * (d1 + d2) * invdx;
        dvzdz = 0.5f * (d1 - d2) * invdz;

        d1 = sin_theta[ix][iz] * dvzdx + cos_theta[ix][iz] * dvzdz;
        dvzdx = cos_theta[ix][iz] * dvzdx - sin_theta[ix][iz] * dvzdz;
        dvzdz = d1;

        sigma_xx[ix][iz] =
            coef1[ix][iz] * sigma_xx[ix][iz] + coef2[ix][iz] * (c11[ix][iz] * dvxdx + c13[ix][iz] * dvzdz);
        sigma_zz[ix][iz] =
            coef1[ix][iz] * sigma_zz[ix][iz] + coef2[ix][iz] * (c13[ix][iz] * dvxdx + c33[ix][iz] * dvzdz);
        sigma_xz[ix][iz] = coef1[ix][iz] * sigma_xz[ix][iz] + coef2[ix][iz] * c55[ix][iz] * (dvzdx + dvxdz);
      }
    }
    /********* free surface boundary condition *********/
    if (freesurface == 1)
    { /* Stress Imaged Method */
      iz = izStart;
#pragma omp parallel for private(ix, iz) firstprivate(izStart)
      for (ix = ixStart; ix < ixEnd; ix++)
      {
        sigma_xz[ix][izStart] = 0.0f;
        sigma_zz[ix][izStart] = 0.0f;
#pragma ivdep
        for (iz = 1; iz <= N; iz++)
        {
          sigma_xx[ix][izStart - iz] = 2.0f * sigma_xx[ix][izStart] - sigma_xx[ix][izStart + iz];
          sigma_zz[ix][izStart - iz] = -sigma_zz[ix][izStart + iz];
          sigma_xz[ix][izStart - iz] = -sigma_xz[ix][izStart + iz];
        }
      }
    }
    /********* Load Source *********/
    if (stype == 1)
    { /* P-wave source (divergence force which load on sigma_xx and sigma_zz) */
      for (isrc = 0; isrc < source.ntrace; isrc++)
      {
        isx = source.posx[isrc] + ealThick;
        isz = source.posz[isrc] + ealThick;
#pragma omp parallel for private(ix, iz, weight_src, rx, rz)                                                         \
    firstprivate(dt, dx, dz, invdx, invdz, alpha, isx, isz, N)
        for (ix = MAX(isx - N, 0); ix <= MIN(isx + N, Nx - 1); ix++)
        {
          rx = (ix - isx) * dx;
#pragma ivdep
          for (iz = MAX(isz - N, 0); iz <= MIN(isz + N, Nz - 1); iz++)
          {
            rz = (iz - isz) * dz;
            weight_src = expf(alpha * (rx * rx + rz * rz)) * invdx * invdz * dt;
            sigma_xx[ix][iz] += 0.5f * weight_src * source.src[isrc][it];
            sigma_zz[ix][iz] += 0.5f * weight_src * source.src[isrc][it];
          }
        }
      }
    }
    /********* snap shot *********/
    if (it % snap_sample == 0 && wsave > 0)
    {
      elastic2dSnapshot(snap_prefix, vx, vz, sigma_xx, sigma_zz, dt, it, Nz, Nx, ealThick, savew, wsave);
    }
    /********* traces *********/
    if (seistype > 0)
    {
      elastic2dReceiver(receiver_vx, receiver_vz, receiver_p, vx, vz, sigma_xx, sigma_zz, it, mt, ealThick, seistype);
    }
    /********* energy *********/
    if (energy == 1)
    {
      ep = 0.0f;
      ek = 0.0f;
#pragma omp parallel for reduction(+ : ek, ep) private(ix, iz, ixx, izz) firstprivate(ealThick)
      for (ix = 0; ix < nx; ix++)
      {
        ixx = ix + ealThick;
        for (iz = 0; iz < nz; iz++)
        {
          izz = iz + ealThick;
          ek += powf(vx[ixx][izz], 2.0f) / invrho[ixx][izz] + powf(vz[ixx][izz], 2.0f) / invrho[ixx][izz];
          ep += 0.25f *
                    (((c11[ixx][izz] + 2.0f * c55[ixx][izz]) * sigma_xx[ixx][izz] -
                      c11[ixx][izz] * sigma_zz[ixx][izz]) *
                         sigma_xx[ixx][izz] +
                     (c11[ixx][izz] * sigma_xx[ixx][izz] -
                      (c11[ixx][izz] + 2.0f * c55[ixx][izz]) * sigma_zz[ixx][izz]) *
                         sigma_zz[ixx][izz]) /
                    (c55[ixx][izz] * (c11[ixx][izz] + c55[ixx][izz])) +
                sigma_xz[ixx][izz] * sigma_xz[ixx][izz] / c55[ixx][izz];
        }
      }
      Ek[it] = 0.5f * ek * dx * dz;
      Ep[it] = 0.5f * ep * dx * dz;
      E[it] = Ek[it] + Ep[it];
    }
    // fprintf(stderr, "%s: Time step: %d\r", argv[0], it);
  }
  tEnd = omp_get_wtime();
  fprintf(stderr, "\n");
  tCost = tEnd - tStart;
  fprintf(stdout, "Finite difference modeling takes %lf s.\n", tCost);
  /********* write traces *********/
  if (seistype == 1 || seistype == 3)
  {
    traces2d_output(trace_prefix, source, receiver_vx, fx, fz, dx, dz, dt, mt, shotID, 1);
    traces2d_output(trace_prefix, source, receiver_vz, fx, fz, dx, dz, dt, mt, shotID, 3);
  }
  if (seistype == 2 || seistype == 3)
  {
    traces2d_output(trace_prefix, source, receiver_p, fx, fz, dx, dz, dt, mt, shotID, 4);
  }
  /********* energy output *********/
  if (energy == 1)
  {
    fp = fopen("src.txt", "w");
    for (it = 0; it < nt; it++)
    {
      fprintf(fp, "%10.5e\n", source.src[0][it]);
    }
    efclose(fp);
    fp = fopen("Eng.txt", "w");
    for (it = 0; it < nt; it++)
    {
      fprintf(fp, "%10.5e\n", E[it]);
    }
    efclose(fp);
    fp = fopen("Ek.txt", "w");
    for (it = 0; it < nt; it++)
    {
      fprintf(fp, "%10.5e\n", Ek[it]);
    }
    efclose(fp);
    fp = fopen("Ep.txt", "w");
    for (it = 0; it < nt; it++)
    {
      fprintf(fp, "%10.5e\n", Ep[it]);
    }
    efclose(fp);
  }
  /********* free memery *********/
  free2float(vp);
  free2float(vs);
  free2float(rho);
  free2float(epsilon);
  free2float(delta);
  free2float(theta);
  free2float(invrho);
  free2float(c11);
  free2float(c13);
  free2float(c33);
  free2float(c55);
  free2float(vx);
  free2float(vz);
  free2float(sigma_xx);
  free2float(sigma_zz);
  free2float(sigma_xz);
  free2float(sin_theta);
  free2float(cos_theta);
  free2float(sin_theta_v);
  free2float(cos_theta_v);
  free2float(D);
  free2float(Beta);
  free2float(coef1);
  free2float(coef2);
  if (energy == 1)
  {
    free1float(E);
    free1float(Ek);
    free1float(Ep);
  }
  freeSource2d(&source);
  freeReceiver2d(&receiver_p);
  freeReceiver2d(&receiver_vx);
  freeReceiver2d(&receiver_vz);
  /********* runing time *********/
  timeEnd = omp_get_wtime();
  timeCost = timeEnd - timeStart;
  fprintf(stdout, "Total time takes %lf s.\n", timeCost);
  return (CWP_Exit());
}