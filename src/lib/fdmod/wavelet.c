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
/** wavelet **/
/**
 * Ricker wavelet
 *  Ricker wavelet
 * @param float *src: Ricker wavelet
 * @param float A: peak amplitude
 * @param float fpeak: peak frequency
 * @param float delay: delay time
 * @param float dt: the time step
 * @param int nt: umber of samples along t
 * @LastEditors: LaurentCheung
 * @Version: 1.0
 */
void riker(float *src, float A, float fpeak, float delay, float dt, int nt)
{
  for (int i = 0; i < nt; i++)
  {
    if (i * dt <= 2 * delay)
    {
      src[i] = A * (1 - 2 * powf(PI * fpeak * (i * dt - delay), 2)) *
               expf(-1 * powf(PI * fpeak * (i * dt - delay), 2));
    }
    else
    {
      src[i] = 0;
    }
  }
}
/**
 * primitive function of Riker wavelet
 *  primitive function of Riker wavelet
 * @param float *src: primitive function of Riker wavelet
 * @param float A: peak amplitude
 * @param float fpeak: peak frequency
 * @param float delay: delay time
 * @param float dt: the time step
 * @param int nt: umber of samples along t
 * @LastEditors: LaurentCheung
 * @Version: 1.0
 */
void iriker(float *src, float A, float fpeak, float delay, float dt, int nt)
{
  for (int i = 0; i < nt; i++)
  {
    if (i * dt <= 2 * delay)
    {
      src[i] = A * (i * dt - delay) * expf(-1 * powf(fpeak * PI * (i * dt - delay), 2));
    }
    else
    {
      src[i] = 0;
    }
  }
}
/**
 * Gaussian function
 *  Gaussian function
 * @param float *src: Gaussian function
 * @param float A: peak amplitude
 * @param float fpeak: peak frequency
 * @param float delay: delay time
 * @param float dt: the time step
 * @param int nt: umber of samples along t
 * @Version: 1.0
 */
void gauss(float *src, float A, float fpeak, float delay, float dt, int nt)
{
  for (int i = 0; i < nt; i++)
  {
    if (i * dt <= 2 * delay)
    {
      src[i] = A * expf(-2 * powf(PI * fpeak * (i * dt - delay), 2));
    }
    else
    {
      src[i] = 0;
    }
  }
}
/**
 *  First derivative of Gaussian function
 *   First derivative of Gaussian function
 * @param float *src:  First derivative of Gaussian function
 * @param float A: peak amplitude
 * @param float fpeak: peak frequency
 * @param float delay: delay time
 * @param float dt: the time step
 * @param int nt: umber of samples along t
 * @LastEditors: LaurentCheung
 * @Version: 1.0
 */
void fdgauss(float *src, float A, float fpeak, float delay, float dt, int nt)
{
  for (int i = 0; i < nt; i++)
  {
    if (i * dt <= 2 * delay)
    {
      src[i] = -4 * A * (i * dt - delay) * powf(PI * fpeak, 2) *
               expf(-2 * powf(PI * fpeak * (i * dt - delay), 2));
    }
    else
    {
      src[i] = 0;
    }
  }
}