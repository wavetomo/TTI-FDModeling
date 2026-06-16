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

/**
 *
 * @param {char} *fileName
 * @param float *
 * @param int Nz
 * @param int Nx
 * @param int ealThick
 * @param int savew
 * @return {*}
 * @LastEditors: LaurentCheung
 * @Version: 1.0
 */
void write_snapshot2d(char *fileName, float **snap, int Nz, int Nx, int ealThick, int savew)
{
  FILE *fp = NULL;
  int ix;
  fp = efopen(fileName, "wb");
  for (ix = (savew == 0 ? ealThick : 0); ix < (savew == 0 ? Nx - ealThick : Nx); ix++)
  {
    efwrite(snap[ix] + (savew == 0 ? ealThick : 0), FSIZE, (savew == 0 ? Nz - 2 * ealThick : Nz),
            fp);
  }
  efclose(fp);
}
/**
 *
 * @param {char} *fileName
 * @param float *
 * @param int Nz
 * @param int Nx
 * @param int Ny
 * @param int ealThick
 * @param int savew
 * @return {*}
 * @LastEditors: LaurentCheung
 * @Version: 1.0
 */
void write_snapshot3d(char *fileName, float ***snap, int Nz, int Nx, int Ny, int ealThick,
                      int savew)
{
  FILE *fp = NULL;
  int iy, ix;
  fp = efopen(fileName, "wb");
  for (iy = (savew == 0 ? ealThick : 0); iy < (savew == 0 ? Ny - ealThick : Ny); iy++)
  {
    for (ix = (savew == 0 ? ealThick : 0); ix < (savew == 0 ? Nx - ealThick : Nx); ix++)
    {
      efwrite(snap[iy][ix] + (savew == 0 ? ealThick : 0), FSIZE,
              (savew == 0 ? Nz - 2 * ealThick : Nz), fp);
    }
  }
  efclose(fp);
}