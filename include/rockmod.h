#pragma once

#ifdef __cplusplus
extern "C"
{
#endif
  /**
   * gridding of 2D interfaces
   * @param float *interface: coordinates of interfaces at z-axis
   * @param float fz: first coordinates at z-axis
   * @param float dz: sample interval at z-axis
   * @param int *fixpos: interface position after gridding
   * @param int nx: the number of samples at x-axis
   * @param int nz: the number of samples at z-axis
   * @return int *fixpos: interface position after gridding
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-10-06 15:20:17
   */
  void interface2d_gridding(float *interface, float fz, float dz, int *fixpos, int nx, int nz);
  /**
   * gridding of 2D interfaces
   * @param float **interface: coordinates of interfaces at z-axis
   * @param float fz: first coordinates at z-axis
   * @param float dz: sample interval at z-axis
   * @param int **fixpos: interface position after gridding
   * @param int nx: the number of samples at x-axis
   * @param int ny: the number of samples at y-axis
   * @param int nz: the number of samples at z-axis
   * @return int *fixpos: interface position after gridding
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-10-06 15:20:17
   */
  void interface3d_gridding(float **interface, float fz, float dz, int **fixpos, int nx, int ny,
                            int nz);
  /**
   * convert velocity to slowness
   * @param float *vp: velocity in m/s
   * @param float *slowness: slowness in s/m
   * @param int n: number of points
   * @return float *slowness: slowness in s/m
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2023-06-11 13:47:56
   */
  void vel2slowness(float *vp, float *slowness, int n);

  /**
   * convert P-wave velocity to S-wave velocity by Vs/Vp ratios of Possion media
   * @param float *vp: P-wave velocity in m/s
   * @param float *vs: S-wave velocity in m/s
   * @param int n: number of points
   * @return float *vs: S-wave velocity in m/s
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2023-06-11 13:47:56
   */
  void vp2vs(float *vp, float *vs, int n);

  /**
   * convert P-wave velocity to S-wave velocity by Fliedner's method
   * @param float *vp: P-wave velocity in m/s
   * @param float *vs: S-wave velocity in m/s
   * @param int n: number of points
   * @return float *vs: S-wave velocity in m/s
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2023-06-11 13:47:56
   */
  void vp2vsFliedner(float *vp, float *vs, int n);

  /**
   * convert P-wave velocity to density
   * @param float *vp: P-wave velocity in m/s
   * @param float *rho: density in kg/m^3
   * @param int n: number of points
   * @return float *rho: density in kg/m^3
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2023-06-11 13:47:56
   */
  void vp2density(float *vp, float *rho, int n);

  /**
   * convert P-wave velocity to density by Gardner relation
   * @param float *vp: P-wave velocity in m/s
   * @param float *rho: density in kg/m^3
   * @param int n: number of points
   * @return float *rho: density in kg/m^3
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2023-06-11 13:47:56
   */
  void vp2densityGardner(float *vp, float *rho, int n);

  /**
   * convert P-wave velocity to density by Ludwig relation
   * @param float *vp: P-wave velocity in m/s
   * @param float *rho: density in kg/m^3
   * @param int n: number of points
   * @return float *rho: density in kg/m^3
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2023-06-11 13:47:56
   */
  void vp2densityLudwig(float *vp, float *rho, int n);

  /**
   * convert P-wave velocity to density by Hamilton's method
   * @param float *vp: P-wave velocity in m/s
   * @param float *rho: density in kg/m^3
   * @param int n: number of points
   * @return float *rho: density in kg/m^3
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2023-06-11 13:47:56
   */
  void vp2densityHamilton(float *vp, float *rho, int n);

  /**
   * fix the water layer's parameters in the top of 2D model
   * @param float **vp: P-wave velocity in m/s
   * @param float **vs: S-wave velocity in m/s
   * @param float **rho: density in kg/m^3
   * @param float vp_water: P-Wave velocity of water
   * @param float rho_water: density of water
   * @param int *interface: the interface of the bottom of water layer or sea bed
   * @param int nx: the number of samples in the x dimension
   * @param int nz: the number of samples in the z dimension
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-10-13 17:06:38
   */
  void waterlayer2dfix(float **vp, float **vs, float **rho, float vp_water, float rho_water,
                       int *interface, int nx, int nz);

  /**
   * fix the water layer's parameters in the top of 3D model
   * @param float ***vp: P-wave velocity in m/s
   * @param float ***vs: S-wave velocity in m/s
   * @param float ***rho: density in kg/m^3
   * @param float vp_water: P-Wave velocity of water
   * @param float rho_water: density of water
   * @param int *interface: the interface of the bottom of water layer or sea bed
   * @param int nx: the number of samples in the x dimension
   * @param int ny: the number of samples in the y dimension
   * @param int nz: the number of samples in the z dimension
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-10-13 17:06:38
   */
  void waterlayer3dfix(float ***vp, float ***vs, float ***rho, float vp_water, float rho_water,
                       int **interface, int nx, int ny, int nz);
#ifdef __cplusplus
}
#endif