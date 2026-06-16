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
/********** anisotropic elastic parameters **********/
/**
 * calculate 2D VTI elastic coefficients from
 * P- and S-wave velocities, density, Thomsen's parameters
 * @param float *C11: elastic coefficient
 * @param float *C13: elastic coefficient
 * @param float *C33: elastic coefficient
 * @param float *C55: elastic coefficient
 * @param float *vp: P-wave velocity
 * @param float *vs: S-wave velocity
 * @param float *rho: density
 * @param float *epsilon: Thomsen's parameter
 * @param float *delta: Thomsen's parameter
 * @param int np: the length of pointers
 * @return float *Cij: elastic coefficient
 * @Author: PingMin Cheung
 * @Date: 2024-08-01 10:37:34
 */
void vti_elasticCoefficients2d(float *C11, float *C13, float *C33, float *C55, float *vp, float *vs, float *rho,
                               float *epsilon, float *delta, int np)
{
  int i;
  double f;
  if (vs == NULL)
  {
#pragma omp parallel for private(i)
    for (i = 0; i < np; i++)
    {
      C11[i] = rho[i] * vp[i] * vp[i] * (1.0 + 2.0 * epsilon[i]);
      C13[i] = rho[i] * vp[i] * vp[i] * sqrt(1.0 + 2.0 * delta[i]);
      C33[i] = rho[i] * vp[i] * vp[i];
    }
  }
  else
  {
#pragma omp parallel for private(i, f)
    for (i = 0; i < np; i++)
    {
      f = 1.0 - powf(vs[i] / vp[i], 2.0);
      C11[i] = rho[i] * vp[i] * vp[i] * (1.0 + 2.0 * epsilon[i]);
      C13[i] = rho[i] * (vp[i] * vp[i] * sqrtf(f + (f + 2.0 * delta[i])) - vs[i] * vs[i]);
      C33[i] = rho[i] * vp[i] * vp[i];
      C55[i] = rho[i] * vs[i] * vs[i];
    }
  }
}
/**
 * calculate 2D TTI elastic coefficients from
 * P- and S-wave velocities, density, Thomsen's parameters
 * @param float *C11: elastic coefficient
 * @param float *C13: elastic coefficient
 * @param float *C15: elastic coefficient
 * @param float *C33: elastic coefficient
 * @param float *C35: elastic coefficient
 * @param float *C55: elastic coefficient
 * @param float *vp: P-wave velocity
 * @param float *vs: S-wave velocity
 * @param float *rho: density
 * @param float *epsilon: Thomsen's parameter
 * @param float *delta: Thomsen's parameter
 * @param int np: the length of pointers
 * @return float *Cij: elastic coefficient
 * @Author: PingMin Cheung
 * @Date: 2024-08-01 10:37:34
 */
void tti_elasticCoefficients2d(float *C11, float *C13, float *C15, float *C33, float *C35, float *C55, float *vp,
                               float *vs, float *rho, float *epsilon, float *delta, float *theta, int np)
{
  int i;
  double c11, c13, c33, c55, f;
  double sin2_theta, cos2_theta, sin4_theta, cos4_theta, sin_2theta, cos_2theta, sin2_2theta, cos2_2theta;
  if (vs == NULL)
  {
#pragma omp parallel for private(i, c11, c13, c33, sin2_theta, cos2_theta, sin4_theta, cos4_theta, sin_2theta,       \
                                     cos_2theta, sin2_2theta, cos2_2theta)
    for (i = 0; i < np; i++)
    {
      c11 = rho[i] * vp[i] * vp[i] * (1.0 + 2.0 * epsilon[i]);
      c13 = rho[i] * vp[i] * vp[i] * sqrt(1.0 + 2.0 * delta[i]);
      c33 = rho[i] * vp[i] * vp[i];

      sin2_theta = sin(theta[i]) * sin(theta[i]);
      cos2_theta = cos(theta[i]) * cos(theta[i]);
      sin4_theta = sin2_theta * sin2_theta;
      cos4_theta = cos2_theta * cos2_theta;
      sin_2theta = sin(2.0 * theta[i]);
      cos_2theta = cos(2.0 * theta[i]);
      sin2_2theta = sin_2theta * sin_2theta;
      cos2_2theta = cos_2theta * cos_2theta;

      C11[i] = c11 * cos4_theta + c33 * sin4_theta + 0.5 * c13 * sin2_2theta;
      C13[i] = 0.25 * (c11 + c33) * sin2_2theta + c13 * (sin4_theta + cos4_theta);
      C15[i] = 0.5 * (-1.0 * c11 * cos2_theta + c33 * sin2_theta + c13 * cos_2theta) * sin_2theta;
      C33[i] = c11 * sin4_theta + c33 * cos4_theta + 0.5 * c13 * sin2_2theta;
      C35[i] = 0.5 * (-1.0 * c11 * sin2_theta + c33 * cos2_theta - c13 * cos_2theta) * sin_2theta;
      C55[i] = 0.25 * (c11 + c33 - 2.0 * c13) * sin2_2theta;
    }
  }
  else
  {
#pragma omp parallel for private(i, c11, c13, c33, c55, f, sin2_theta, cos2_theta, sin4_theta, cos4_theta,           \
                                     sin_2theta, cos_2theta, sin2_2theta, cos2_2theta)
    for (i = 0; i < np; i++)
    {
      f = 1.0 - pow(vs[i] / vp[i], 2.0);
      c11 = rho[i] * vp[i] * vp[i] * (1.0 + 2.0 * epsilon[i]);
      c13 = rho[i] * (vp[i] * vp[i] * sqrt(f * (f + 2.0 * delta[i])) - vs[i] * vs[i]);
      c33 = rho[i] * vp[i] * vp[i];
      c55 = rho[i] * vs[i] * vs[i];

      sin2_theta = sin(theta[i]) * sin(theta[i]);
      cos2_theta = cos(theta[i]) * cos(theta[i]);
      sin4_theta = sin2_theta * sin2_theta;
      cos4_theta = cos2_theta * cos2_theta;
      sin_2theta = sin(2.0 * theta[i]);
      cos_2theta = cos(2.0 * theta[i]);
      sin2_2theta = sin_2theta * sin_2theta;
      cos2_2theta = cos_2theta * cos_2theta;

      C11[i] = c11 * cos4_theta + c33 * sin4_theta + (0.5 * c13 + c55) * sin2_2theta;
      C13[i] = (0.25 * (c11 + c33) - c55) * sin2_2theta + c13 * (sin4_theta + cos4_theta);
      C15[i] = 0.5 * (-1.0 * c11 * cos2_theta + c33 * sin2_theta + (c13 + 2.0 * c55) * cos_2theta) * sin_2theta;
      C33[i] = c11 * sin4_theta + c33 * cos4_theta + (0.5 * c13 + c55) * sin2_2theta;
      C35[i] = 0.5 * (-1.0 * c11 * sin2_theta + c33 * cos2_theta - (c13 + 2.0 * c55) * cos_2theta) * sin_2theta;
      C55[i] = 0.25 * (c11 + c33 - 2.0 * c13) * sin2_2theta + c55 * cos2_2theta;
    }
  }
}
/**
 * calculate 3D VTI elastic coefficients from
 * P- and S-wave velocities, density, Thomsen's parameters
 * @param float *C11: elastic coefficient
 * @param float *C12: elastic coefficient
 * @param float *C13: elastic coefficient
 * @param float *C33: elastic coefficient
 * @param float *C55: elastic coefficient
 * @param float *C66: elastic coefficient
 * @param float *vp: P-wave velocity
 * @param float *vs: S-wave velocity
 * @param float *rho: density
 * @param float *epsilon: Thomsen's parameter
 * @param float *delta: Thomsen's parameter
 * @param int np: the length of pointers
 * @return float *Cij: elastic coefficient
 * @Author: PingMin Cheung
 * @Date: 2024-08-01 10:37:34
 */
void vti_elasticCoefficients3d(float *C11, float *C12, float *C13, float *C33, float *C55, float *C66, float *vp,
                               float *vs, float *rho, float *epsilon, float *delta, float *gama, int np)
{
  int i;
  float f;
  if (vs == NULL)
  {
#pragma omp parallel for private(i)
    for (i = 0; i < np; i++)
    {
      C11[i] = rho[i] * vp[i] * vp[i] * (1.0f + 2.0f * epsilon[i]);
      C13[i] = rho[i] * vp[i] * vp[i] * sqrtf(1.0f + 2.0f * delta[i]);
      C33[i] = rho[i] * vp[i] * vp[i];
    }
  }
  else
  {
#pragma omp parallel for private(i, f)
    for (i = 0; i < np; i++)
    {
      f = 1.0f - powf(vs[i] / vp[i], 2.0f);
      C11[i] = rho[i] * vp[i] * vp[i] * (1.0f + 2.0f * epsilon[i]);
      C12[i] = rho[i] * (vp[i] * vp[i] * (1.0f + 2.0f * epsilon[i]) - 2.0f * vs[i] * vs[i] * (1.0f + 2.0f * gama[i]));
      C13[i] = rho[i] * (vp[i] * vp[i] * sqrtf(f + (f + 2.0f * delta[i])) - vs[i] * vs[i]);
      C33[i] = rho[i] * vp[i] * vp[i];
      C55[i] = rho[i] * vs[i] * vs[i];
      C66[i] = rho[i] * vs[i] * vs[i] * (1.0f + 2.0f * gama[i]);
    }
  }
}