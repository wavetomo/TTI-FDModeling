#pragma once
#ifdef __cplusplus
extern "C"
{
#endif
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
  void gauss_wlt(float *src, float A, float fpeak, float delay, float dt, int nt);
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
  void fdgauss_wlt(float *src, float A, float fpeak, float delay, float dt, int nt);
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
  void riker_wlt(float *src, float A, float fpeak, float delay, float dt, int nt);
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
  void driker_wlt(float *src, float A, float fpeak, float delay, float dt, int nt);
#ifdef __cplusplus
}
#endif