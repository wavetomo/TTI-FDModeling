/************************************/
#include <wavelet.h>
/**/
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/************************************/
/**
 * Gaussian function
 * @description: Gaussian function
 * @param float *src: Gaussian function
 * @param float A: peak amplitude
 * @param float fpeak: peak frequency
 * @param float delay: delay time
 * @param float dt: the time step
 * @param int nt: umber of samples along t
 * @Version: 1.0
 */
void gauss_wlt(float *src, float A, float fpeak, float delay, float dt, int nt)
{
  float t, pi;
  pi = 3.141592654f;
  for (int it = 0; it < nt; it++)
  {
    t = it * dt - delay;
    if (t <= delay)
    {
      src[it] = 0.5f / powf(pi * fpeak, 2.0f) * expf(-powf(pi * fpeak * t, 2.0f));
    }
    else
    {
      src[it] = 0.0f;
    }
  }
  src[0] = 0.0f;
}
/**
 *  First derivative of Gaussian function
 * @description:  First derivative of Gaussian function
 * @param float *src:  First derivative of Gaussian function
 * @param float A: peak amplitude
 * @param float fpeak: peak frequency
 * @param float delay: delay time
 * @param float dt: the time step
 * @param int nt: umber of samples along t
 * @LastEditors: LaurentCheung
 * @Version: 1.0
 */
void fdgauss_wlt(float *src, float A, float fpeak, float delay, float dt, int nt)
{
  float t, pi;
  pi = 3.141592654f;
  for (int it = 0; it < nt; it++)
  {
    t = it * dt - delay;
    if (t <= delay)
    {
      src[it] = A * t * expf(-powf(pi * fpeak * t, 2.0f));
    }
    else
    {
      src[it] = 0.0f;
    }
  }
  src[0] = 0.0f;
}
/**
 * Ricker wavelet
 * @description: Ricker wavelet
 * @param float *src: Ricker wavelet
 * @param float A: peak amplitude
 * @param float fpeak: peak frequency
 * @param float delay: delay time
 * @param float dt: the time step
 * @param int nt: umber of samples along t
 * @LastEditors: LaurentCheung
 * @Version: 1.0
 */
void riker_wlt(float *src, float A, float fpeak, float delay, float dt, int nt)
{
  float t, pi;
  pi = 3.141592654f;
  for (int it = 0; it < nt; it++)
  {
    t = it * dt - delay;
    if (t <= delay)
    {
      src[it] =
          A * (1.0f - 2.0f * powf(pi * fpeak * t, 2.0f)) * expf(-1.0f * powf(pi * fpeak * t, 2.0f));
    }
    else
    {
      src[it] = 0.0f;
    }
  }
  src[0] = 0.0f;
}
/**
 * Ricker wavelet
 * @description: Ricker wavelet
 * @param float *src: Ricker wavelet
 * @param float A: peak amplitude
 * @param float fpeak: peak frequency
 * @param float delay: delay time
 * @param float dt: the time step
 * @param int nt: umber of samples along t
 * @LastEditors: LaurentCheung
 * @Version: 1.0
 */
void driker_wlt(float *src, float A, float fpeak, float delay, float dt, int nt)
{
  float t, pi;
  pi = 3.141592654f;
  for (int it = 0; it < nt; it++)
  {
    t = it * dt - delay;
    if (t <= delay)
    {
      src[it] = A * (2.0f * powf(pi * fpeak, 2.0f) * t) *
                (2.0f * powf(pi * fpeak * t, 2.0f) - 3.0f) *
                expf(-1.0f * powf(pi * fpeak * t, 2.0f));
    }
    else
    {
      src[it] = 0.0f;
    }
  }
  src[0] = 0.0f;
}