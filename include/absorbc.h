#pragma once
#ifdef __cplusplus
extern "C"
{
#endif
  /********* model pad and extend *********/
  /**
   * padding 2D model
   * @description: padding 2D model
   * @param float** mod: model array
   * @param int Nx: the sample number along x dimension
   * @param int Nz: the sample number along z dimension
   * @param int nlayer: extend layers
   * @LastEditors: LaurentCheung
   * @Version: 1.0
   */
  void padModel2d(float **mod, int Nx, int Nz, int nlayer);
  /**
   * extend 2D model
   * @description: extend 2D model
   * @param float** mod: model array
   * @param int Nx: the sample number along x dimension
   * @param int Nz: the sample number along z dimension
   * @param int nlayer: extend layers
   * @LastEditors: LaurentCheung
   * @Version: 1.0
   */
  void extendModel2d(float **modEx, float **mod, int Nx, int Nz, int nlayer);
  /**
   * padding 3D model
   * @description: padding 2D model
   * @param float** mod: model array
   * @param int Nx: the sample number along x dimension
   * @param int Ny: the sample number along y dimension
   * @param int Nz: the sample number along z dimension
   * @param int nlayer: extend layers
   * @LastEditors: LaurentCheung
   * @Version: 1.0
   */
  void padModel3d(float ***mod, int Nx, int Ny, int Nz, int nlayer);
  /**
   * extend 3D model
   * @description: extend 3D model
   * @param float** mod: model array
   * @param int Nx: the sample number along x dimension
   * @param int Ny: the sample number along y dimension
   * @param int Nz: the sample number along z dimension
   * @param int nlayer: extend layers
   * @LastEditors: LaurentCheung
   * @Version: 1.0
   */
  void extendModel3d(float ***modEx, float ***mod, int Nx, int Ny, int Nz, int nlayer);
  /********* Clayton absorbing boundary condition *********/
  /**
   * @description:
   * @param float **p
   * @param float **pn
   * @param float **pn_1
   * @param float **vp
   * @param float dx
   * @param float dz
   * @param float dt
   * @param int Nx
   * @param int Nz
   * @param int nlayer
   * @param int freesurface
   * @return {*
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-03-26 10:59:27
   */
  void acoustic2dClayton(float **p, float **pn, float **pn_1, float **vp, float dx, float dz,
                         float dt, int Nx, int Nz, int nlayer, int freesurface);
  /**
   * @description:
   * @param float **p
   * @param float **px
   * @param float **pz
   * @param float **vp
   * @param float dx
   * @param float dz
   * @param float dt
   * @param int Nx
   * @param int Nz
   * @param int nlayer
   * @param int freesurface
   * @return {*
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-03-26 10:59:27
   */
  void Higdon2d(float **p, float **px, float **pz, float **vp, float dx, float dz, float dt, int Nx,
                int Nz, int nlayer, int freesurface);
  /********* perfectly matched layer absorbing boundary condition *********/

  /********* modified effective absorbing boundary condition *********/
  /**
   * 2D MEAL boundary condition
   * @description: 2D MEAL boundary condition
   * @param float** D: the value of damping function
   * @param float** Beta: scaling factor
   * @param float** vp: velocity of P-wave
   * @param float** vs: velocity of S-wave
   * @param float fpeak: the dominant frequency of the source
   * @param float alpha: the order of Gaussian weight function
   * @param float dz: the sampling interval along z, the 1st dimension
   * @param float dx: the sampling interval along x, the 2nd dimension
   * @param int Nz: The sample number along 1st dimension
   * @param int Nx: The sample number along 2nd dimension
   * @param int nlayer: number of EAL
   * @param int freesurface: is free-surface boundary condition used.
   *                         0 for not, 1 for yes
   * @return float** D: the value of damping function
   * @return float** Beta: scaling factor
   * @LastEditors: LaurentCheung
   * @Version: 1.0
   */
  void mealCoefs2d(float **D, float **Beta, float **vp, float **vs, float fpeak, float alpha,
                   float dz, float dx, int Nz, int Nx, int nlayer, int freesurface);
  /**
   * 3D MEAL boundary condition
   * @description: 2D MEAL boundary condition
   * @param float*** D: the value of damping function
   * @param float*** Beta: scaling factor
   * @param float*** vp: velocity of P-wave
   * @param float*** vs: velocity of S-wave
   * @param float fpeak: the dominant frequency of the source
   * @param float alpha: the order of Gaussian weight function
   * @param float dz: the sampling interval along z, the 1st dimension
   * @param float dx: the sampling interval along x, the 2nd dimension
   * @param float dy: the sampling interval along y, the 3rd dimension
   * @param int Nz: The sample number along 1st dimension
   * @param int Nx: The sample number along 2nd dimension
   * @param int Ny: The sample number along 3rd dimension
   * @param int nlayer: number of EAL
   * @param int freesurface: is free-surface boundary condition used.
   *                         0 for not, 1 for yes
   * @return float*** D: the value of damping function
   * @return float*** Beta: scaling factor
   * @LastEditors: LaurentCheung
   * @Version: 1.0
   */
  void mealCoefs3d(float ***D, float ***Beta, float ***vp, float ***vs, float fpeak, float alpha,
                   float dz, float dx, float dy, int Nz, int Nx, int Ny, int nlayer,
                   int freesurface);
  /**
   * calculate the coefficent of finite difference scheme for first or second order equation
   * @param float *coef1: first coefficent
   * @param float *coef2: second coefficent
   * @param float *coef3: third coefficent
   * @param float *D: the value of damping function
   * @param float *Beta: scaling factor
   * @param float *vp: velocity of P-wave
   * @param float *vs: velocity of S-wave
   * @param float *rho: density
   * @param float dt: time stepsize
   * @param int Np: the number of points
   * @param int order: the order of equation
   * @Author: Pingmin Cheung
   * @Date: 2025-08-26 10:22:59
   */
  void mealFDCoefs(float *coef1, float *coef2, float *coef3, float *D, float *Beta, float *vp,
                   float *vs, float *rho, float dt, int Np, int order);
#ifdef __cplusplus
}
#endif