
#pragma once
#ifndef swaptr
#define swaptr(ptr1, ptr2)                                                                         \
  {                                                                                                \
    void *temp = (ptr1);                                                                           \
    (ptr1) = (ptr2);                                                                               \
    (ptr2) = temp;                                                                                 \
  }
#endif
/******************** TYPEDEFS ********************/
#pragma pack(8)
// line source
// int ntrace: number of trace
// int ntsample: time sample number
// float dt: sampling interval
// int *posx: store the x position of each trace: unit: cell
// int *posz: store the z position of each trace: unit: cell
// float **src: 2d source array, 1st dimension is time, 2nd dimension is x,y
typedef struct Source2D
{
  int ntrace;
  int ntsample;
  int *posx;
  int *posz;
  float **src;
} Source2D;
// line receiver
// int ntrace: number of trace
// int ntsample: time sample number
// int *posx: store the x position of each trace: unit: cell
// int *posz: store the z position of each trace: unit: cell
// float **recv: 2d record array, 1st dimension is time, 2nd dimension is x,y
typedef struct Receiver2D
{
  int ntrace;
  int ntsample;
  int *posx;
  int *posz;
  float **recv;
} Receiver2D;
// line source
// int ntrace: trace number
// int ntsample: time sample number
// float dt: sampling interval
// int *posx: store the x position of each trace, unit: cell
// int *posz: store the z position of each trace, unit: cell
// int *posy: store the y position of each trace, unit: cell
// float **src: 2d source array, 1st dimension is time, 2nd dimension is x,y
typedef struct Source3D
{
  int ntrace;
  int ntsample;
  int *posx;
  int *posz;
  int *posy;
  float **src;
} Source3D;
// line receiver
// int ntrace: trace number
// int ntsample: time sample number
// int *posx: store the x position of each trace, unit: cell
// int *posz: store the z position of each trace, unit: cell
// int *posy: store the y position of each trace, unit: cell
// float **recv: 2d record array, 1st dimension is time, 2nd dimension is x,y
typedef struct
{
  int ntrace;
  int ntsample;
  int *posx;
  int *posz;
  int *posy;
  float **recv;
} Receiver3D;
#pragma pack()
/******************** FUNCTIONS ********************/
#ifdef __cplusplus
extern "C"
{
#endif
  /******************** types initialize ********************/
  /**
   * initialize 2D source structure
   * @param Source2D *source: 2D source structure pointer
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-08-14 14:18:10
   */
  void source2d_init(Source2D *source);
  /**
   * initialize 3D source structure
   * @param Source2D *source: 3D source structure pointer
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-08-14 14:18:10
   */
  void source3d_init(Source3D *source);
  /**
   * initialize 2D receiver structure
   * @param Receiver2D *receiver: 2D receiver structure pointer
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-08-14 14:18:10
   */
  void receiver2d_init(Receiver2D *receiver);
  /**
   * initialize 3D receiver structure
   * @param Receiver3D *receiver: 3D receiver structure pointer
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-08-14 14:18:10
   */
  void receiver3d_init(Receiver3D *receiver);
  /******************** types free ********************/
  /**
   * free memery of Source2D
   *  free memery of Source2D
   * @param {Source2D} source: 2D source structure pointer
   *                               contain the source information
   * @LastEditors: LaurentCheung
   * @Version: 1.0
   */
  void freeSource2d(Source2D *source);
  /**
   * free memery of Source3D
   *  free memery of Source3D
   * @param {Source3D} source: 3D source structure pointer
   *                               contain the source information
   * @LastEditors: LaurentCheung
   * @Version: 1.0
   */
  void freeSource3d(Source3D *source);
  /**
   * free memery of Receiver2D
   *  free memery of Receiver2D
   * @param {Receiver2D} receiver: 2D receiver structure pointer
   *                               contain the receiver information
   * @LastEditors: LaurentCheung
   * @Version: 1.0
   */
  void freeReceiver2d(Receiver2D *receiver);
  /**
   * free memery of Receiver3D
   *  free memery of Receiver3D
   * @param {Receiver3D} receiver: 3D receiver structure pointer
   *                               contain the receiver information
   * @LastEditors: LaurentCheung
   * @Version: 1.0
   */
  void freeReceiver3d(Receiver3D *receiver);
  /******************** types copy ********************/
  /**
   * copy 2D receiver structure
   * @param Receiver2D *dest: 2D receiver structure pointer
   * @param Receiver2D *src: 2D receiver structure pointer
   * @return *
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-09-01 18:03:07
   */
  void receiver2d_copy(Receiver2D *dest, Receiver2D *src);
  /**
   * copy 3D receiver structure
   * @param Receiver2D *dest: 2D receiver structure pointer
   * @param Receiver2D *src: 2D receiver structure pointer
   * @return *
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-09-01 18:03:07
   */
  void receiver3d_copy(Receiver3D *dest, Receiver3D *src);
  /******************** observation ********************/
  /**
   * read source wavelet and position from file
   * @param Source2D *source: 2D source structure pointer
   * @param char *sourcefile: file name of source
   * @param float dt: sample interval in time dimension
   * @param int nx: number of samples in x-space dimension
   * @param int nz: number of samples in z-space dimension
   * @param int nt: number of samples in time dimension
   * @return Source2D *source: 2D source structure
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-08-14 14:13:19
   */
  void read2dsource(Source2D *source, char *sourcefile, float dt, int nx, int nz, int nt);
  /**
   * read source wavelet and position from file
   * @param Source2D *source: 2D source structure pointer
   * @param char *sourcefile: file name of source
   * @param float dt: sample interval in time dimension
   * @param int nx: number of samples in x-space dimension
   * @param int ny: number of samples in y-space dimension
   * @param int nz: number of samples in z-space dimension
   * @param int nt: number of samples in time dimension
   * @return Source2D *source: 2D source structure
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-08-14 14:13:19
   */
  void read3dsource(Source3D *source, char *sourcefile, float dt, int nx, int ny, int nz, int nt);
  /**
   * read 2D receiver's information from file
   * @param Receiver2D *receiver: 2D receiver structure pointer
   * @param char *recposfile: file name of receiver
   * @param int nx: number of samples in x-space dimension
   * @param int nz: number of samples in z-space dimension
   * @param int nt: number of samples in time dimension
   * @param int mt: receiver's sample interval in time dimension
   * @return Receiver2D *receiver: 2D receiver structure pointer
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-08-17 14:26:13
   */
  void read2dreceiver(Receiver2D *receiver, char *recposfile, int nx, int nz, int nt, int mt);
  /**
   * read 3D receiver's information from file
   * @param Receiver3D *receiver: 3D receiver structure pointer
   * @param char *recposfile: file name of receiver
   * @param int nx: number of samples in x-space dimension
   * @param int ny: number of samples in y-space dimension
   * @param int nz: number of samples in z-space dimension
   * @param int nt: number of samples in time dimension
   * @param int mt: receiver's sample interval in time dimension
   * @return Receiver3D *receiver: 3D receiver structure pointer
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-08-17 14:26:13
   */
  void read3dreceiver(Receiver3D *receiver, char *recposfile, int nx, int ny, int nz, int nt,
                      int mt);
  /******************** receive wavefeild ********************/
  /**
   *
   * @param {Receiver2D} receiver_p
   * @param {Receiver2D} receiver_vx
   * @param {Receiver2D} receiver_vz
   * @param float **p
   * @param float **vx
   * @param float **vz
   * @param int it
   * @param int mt
   * @param int ealThick
   * @param int seistype
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-05-16 21:33:26
   */
  void acoustic2dReceiver(Receiver2D receiver_p, Receiver2D receiver_vx, Receiver2D receiver_vz,
                          float **p, float **vx, float **vz, int it, int mt, int ealThick,
                          int seistype);
  /**
   *
   * @param {Receiver3D} receiver_p
   * @param {Receiver3D} receiver_vx
   * @param {Receiver3D} receiver_vy
   * @param {Receiver3D} receiver_vz
   * @param float ***p
   * @param float ***vx
   * @param float ***vy
   * @param float ***vz
   * @param int it
   * @param int mt
   * @param int ealThick
   * @param int seistype
   * @return {*}
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-05-17 14:30:57
   */
  void acoustic3dReceiver(Receiver3D receiver_p, Receiver3D receiver_vx, Receiver3D receiver_vy,
                          Receiver3D receiver_vz, float ***p, float ***vx, float ***vy, float ***vz,
                          int it, int mt, int ealThick, int seistype);
  /**
   *
   * @param {Receiver2D} receiver_vx
   * @param {Receiver2D} receiver_vz
   * @param {Receiver2D} receiver_p
   * @param float **vx
   * @param float **vz
   * @param float **sigma_xx
   * @param float **sigma_zz
   * @param int it
   * @param int mt
   * @param int ealThick
   * @param int seistype
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-05-16 21:33:26
   */
  void elastic2dReceiver(Receiver2D receiver_vx, Receiver2D receiver_vz, Receiver2D receiver_p,
                         float **vx, float **vz, float **sigma_xx, float **sigma_zz, int it, int mt,
                         int ealThick, int seistype);
  /**
   *
   * @param {Receiver3D} receiver_vx
   * @param {Receiver3D} receiver_vy
   * @param {Receiver3D} receiver_vz
   * @param {Receiver3D} receiver_p
   * @param float ***vx
   * @param float ***vy
   * @param float ***vz
   * @param float **sigma_xx
   * @param float **sigma_yy
   * @param float **sigma_zz
   * @param int it
   * @param int mt
   * @param int ealThick
   * @param int seistype
   * @return {*}
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-05-17 14:30:57
   */
  void elastic3dReceiver(Receiver3D receiver_p, Receiver3D receiver_vx, Receiver3D receiver_vy,
                         Receiver3D receiver_vz, float ***vx, float ***vy, float ***vz,
                         float ***sigma_xx, float ***sigma_yy, float ***sigma_zz, int it, int mt,
                         int ealThick, int seistype);
  /******************** trace output ********************/
  /**
   *
   * @param {char} *fileName
   * @param {Source2D} source
   * @param {Receiver2D} receiver
   * @param float fx
   * @param float fz
   * @param float dx
   * @param float dz
   * @param float dt
   * @return {*}
   * @LastEditors: LaurentCheung
   * @Version: 1.0
   */
  void putraces_2d(char *fileName, Source2D source, Receiver2D receiver, float fx, float fz,
                   float dx, float dz, float dt, int shotID);
  /**
   *
   * @param {char} *fileName
   * @param {Source2D} source
   * @param {Receiver2D} receiver
   * @param float fx
   * @param float fy
   * @param float fz
   * @param float dx
   * @param float dy
   * @param float dz
   * @param float dt
   * @return {*}
   * @LastEditors: LaurentCheung
   * @Version: 1.0
   */
  void putraces_3d(char *fileName, Source3D source, Receiver3D receiver, float fx, float fy,
                   float fz, float dx, float dy, float dz, float dt, int shotID);
  /**
   *
   * @param {char} *trace_prefix
   * @param {Source2D} source
   * @param {Receiver2D} receiver
   * @param float fx
   * @param float fz
   * @param float dx
   * @param float dz
   * @param float dt
   * @param int mt
   * @param int shotnumber
   * @param int flag
   * @return {*}
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2023-07-16 22:10:54
   */
  void traces2d_output(char *trace_prefix, Source2D source, Receiver2D receiver, float fx, float fz,
                       float dx, float dz, float dt, int mt, int shotnumber, int flag);
  /**
   *
   * @param {char} *trace_prefix:
   * @param {Source3D} source:
   * @param {Receiver3D} receiver:
   * @param float fx:
   * @param float fy:
   * @param float fz:
   * @param float dx:
   * @param float dy:
   * @param float dz:
   * @param float dt:
   * @param int mt:
   * @param int shotnumber:
   * @param int flag:
   * @LastEditors: LaurentCheung
   * @Version: 1.0
   */
  void traces3d_output(char *trace_prefix, Source3D source, Receiver3D receiver, float fx, float fy,
                       float fz, float dx, float dy, float dz, float dt, int mt, int shotnumber,
                       int flag);
  /******************** snapshot output ********************/
  /**
   *
   * @param {char} *snap_filename
   * @param float *
   * @param int Nz
   * @param int Nx
   * @param int ealThick
   * @param int savew
   * @return {*}
   * @LastEditors: LaurentCheung
   * @Version: 1.0
   */
  void write_snapshot2d(char *snap_filename, float **snap, int Nz, int Nx, int ealThick, int savew);
  /**
   *
   * @param {char} *snap_filename
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
  void write_snapshot3d(char *snap_filename, float ***snap, int Nz, int Nx, int Ny, int ealThick,
                        int savew);
  /**
   *
   * @param {char} *prefix
   * @param float **p
   * @param float **vx
   * @param float **vz
   * @param float dt
   * @param int it
   * @param int Nz
   * @param int Nx
   * @param int ealThick
   * @param int savew
   * @param int wsave
   * @return {*}
   * @LastEditors: LaurentCheung
   * @Version: 1.0
   */
  void acoustic2dSnapshot(char *prefix, float **p, float **vx, float **vz, float dt, int it, int Nz,
                          int Nx, int ealThick, int savew, int wsave);
  /**
   *
   * @param {char} *prefix
   * @param float ***p
   * @param float ***vx
   * @param float ***vy
   * @param float ***vz
   * @param float dt
   * @param int it
   * @param int Nz
   * @param int Nx
   * @param int Ny
   * @param int ealThick
   * @param int savew
   * @param int wsave
   * @return {*}
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-05-17 15:43:08
   */
  void acoustic3dSnapshot(char *prefix, float ***p, float ***vx, float ***vy, float ***vz, float dt,
                          int it, int Nz, int Nx, int Ny, int ealThick, int savew, int wsave);
  /**
   *
   * @param {char} *prefix
   * @param float **p
   * @param float **vx
   * @param float **vz
   * @param float dt
   * @param int it
   * @param int Nz
   * @param int Nx
   * @param int ealThick
   * @param int savew
   * @param int wsave
   * @return {*}
   * @LastEditors: LaurentCheung
   * @Version: 1.0
   */
  void elastic2dSnapshot(char *prefix, float **vx, float **vz, float **sigma_xx, float **sigma_zz,
                         float dt, int it, int Nz, int Nx, int ealThick, int savew, int wsave);
  /**
   *
   * @param {char} *prefix
   * @param float ***vx
   * @param float ***vy
   * @param float ***vz
   * @param float ***sigma_xx
   * @param float ***sigma_yy
   * @param float ***sigma_zz
   * @param float dt
   * @param int it
   * @param int Nz
   * @param int Nx
   * @param int Ny
   * @param int ealThick
   * @param int savew
   * @param int wsave
   * @return {*}
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-05-23 20:39:57
   */
  void elastic3dSnapshot(char *prefix, float ***vx, float ***vy, float ***vz, float ***sigma_xx,
                         float ***sigma_yy, float ***sigma_zz, float dt, int it, int Nz, int Nx,
                         int Ny, int ealThick, int savew, int wsave);
  /****** stable analyses ********************/
  /**
   * stable analysis for 2D first order quasi-P wave equation
   * finite difference scheme by staggered grid
   * @param float **vp: P-wave velocity
   * @param float **epsilon: Thomsen's parameter, epsilon
   * @param float **delta: Thomsen's parameter, delta
   * @param float *C: finite difference coefficent
   * @param float dx: the sampling interval in x-direction
   * @param float dz: the sampling interval in z-direction
   * @param float dt: the sampling interval in time-direction
   * @param int Nx: the number of samples in x-direction
   * @param int Nz: the number of samples in z-direction
   * @param int N: the order of  spectial finite difference scheme
   * @Author: Pingmin Cheung
   * @Date: 2025-09-01 12:47:24
   */
  void acoustic_vti2d_sgfd_stable(float **vp, float **epsilon, float **delta, float *C, float dx,
                                  float dz, float dt, int Nx, int Nz, int N);
  /**
   * stable analysis for 2D first order elastic wave equation
   * finite difference scheme by staggered grid
   * @param float **vp: P-wave velocity
   * @param float **vs: S-wave velocity
   * @param float **epsilon: Thomsen's parameter, epsilon
   * @param float **delta: Thomsen's parameter, delta
   * @param float *C: finite difference coefficent
   * @param float dx: the sampling interval in x-direction
   * @param float dz: the sampling interval in z-direction
   * @param float dt: the sampling interval in time-direction
   * @param int Nx: the number of samples in x-direction
   * @param int Nz: the number of samples in z-direction
   * @param int N: the order of  spectial finite difference scheme
   * @Author: Pingmin Cheung
   * @Date: 2025-09-01 12:47:24
   */
  void elastic_vti2d_sgfd_stable(float **vp, float **vs, float **epsilon, float **delta, float *C,
                                 float dx, float dz, float dt, int Nx, int Nz, int N);
  /**
   * stable analysis for 2D second order acoustic wave equation finite difference
   * scheme by regular grid
   *  stable analysis for 2D second order acoustic wave equation finite
   * difference scheme by regular grid
   * @param float **vp: velocity model
   * @param float *C: finite difference coefficents
   * @param float dx: the sampling interval along x dimension
   * @param float dz: the sampling interval along z dimension
   * @param float dt: the sampling interval along t dimension
   * @param int nx: the sample number along x dimension
   * @param int nz: the sample number along z dimension
   * @param int N: the order of spatical finite difference
   * @return {*}
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-06-08 20:59:42
   */
  void regular_grid_stable2d(float **vp, float *C, float dx, float dz, float dt, int nx, int nz,
                             int N);
  /**
   * stable analysis for 3D second order acoustic wave equation finite difference
   * scheme by regular grid
   *  stable analysis for 3D second order acoustic wave equation finite
   * difference scheme by regular grid
   * @param float ***vp: velocity model
   * @param float *C: finite difference coefficents
   * @param float dx: the sampling interval along x dimension
   * @param float dy: the sampling interval along y dimension
   * @param float dz: the sampling interval along z dimension
   * @param float dt: the sampling interval along t dimension
   * @param int nx: the sample number along x dimension
   * @param int ny: the sample number along y dimension
   * @param int nz: the sample number along z dimension
   * @param int N: the order of spatical finite difference
   * @return {*}
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-06-08 20:59:42
   */
  void regular_grid_stable3d(float ***vp, float *C, float dx, float dy, float dz, float dt, int nx,
                             int ny, int nz, int N);
  /**
   * stable analysis for 2D first order acoustic/elastic wave equation finite difference scheme by
   * staggered grid
   * @param float **vp: velocity model
   * @param float *C: finite difference coefficents
   * @param float dx: the sampling interval along x dimension
   * @param float dz: the sampling interval along z dimension
   * @param float dt: the sampling interval along t dimension
   * @param int nx: the sample number along x dimension
   * @param int nz: the sample number along z dimension
   * @param int N: the order of spatical finite difference
   * @return {*}
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-06-08 20:59:42
   */
  void staggered_grid_stable2d(float **vp, float *C, float dx, float dz, float dt, int nx, int nz,
                               int N);
  /**
   * stable analysis for 3D first order acoustic/elastic wave equation finite difference scheme by
   * staggered grid
   * @param float ***vp: velocity model
   * @param float *C: finite difference coefficents
   * @param float dx: the sampling interval along x dimension
   * @param float dy: the sampling interval along y dimension
   * @param float dz: the sampling interval along z dimension
   * @param float dt: the sampling interval along t dimension
   * @param int nx: the sample number along x dimension
   * @param int ny: the sample number along y dimension
   * @param int nz: the sample number along z dimension
   * @param int N: the order of spatical finite difference
   * @return {*}
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-06-08 20:59:42
   */
  void staggered_grid_stable3d(float ***vp, float *C, float dx, float dy, float dz, float dt,
                               int nx, int ny, int nz, int N);
  /**
   * stable analysis for 2D first order acoustic wave equation finite difference
   * scheme by rotated staggered grid
   *  stable analysis for 2D first order acoustic wave equation finite
   * difference scheme by rotated staggered grid
   * @param float **vp: velocity model
   * @param float *C: finite difference coefficents
   * @param float dx: the sampling interval along x dimension
   * @param float dz: the sampling interval along z dimension
   * @param float dt: the sampling interval along t dimension
   * @param int nx: the sample number along x dimension
   * @param int nz: the sample number along z dimension
   * @param int N: the order of spatical finite difference
   * @return {*}
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-06-08 20:59:42
   */
  void rotated_staggered_grid_stable2d(float **vp, float *C, float dx, float dz, float dt, int nx,
                                       int nz, int N);
  /**
   * stable analysis for 3D first order acoustic wave equation finite difference
   * scheme by rotated staggered grid
   *  stable analysis for 3D first order acoustic wave equation finite
   * difference scheme by rotated staggered grid
   * @param float ***vp: velocity model
   * @param float *C: finite difference coefficents
   * @param float dx: the sampling interval along x dimension
   * @param float dy: the sampling interval along y dimension
   * @param float dz: the sampling interval along z dimension
   * @param float dt: the sampling interval along t dimension
   * @param int nx: the sample number along x dimension
   * @param int ny: the sample number along y dimension
   * @param int nz: the sample number along z dimension
   * @param int N: the order of spatical finite difference
   * @return {*}
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-06-08 20:59:42
   */
  void rotated_staggered_grid_stable3d(float ***vp, float *C, float dx, float dy, float dz,
                                       float dt, int nx, int ny, int nz, int N);
  /**
   * stable analysis for regular grid pesudo spectral method by 2D second order
   * acoustic wave equation
   *  stable analysis for regular grid pesudo spectral method by 2D
   * second order acoustic wave equation
   * @param float **vp: velocity model
   * @param float dx: the sampling interval along x dimension
   * @param float dz: the sampling interval along z dimension
   * @param float dt: the sampling interval along t dimension
   * @param int nx: the sample number along x dimension
   * @param int nz: the sample number along z dimension
   * @return {*}
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-06-09 15:00:17
   */
  void regular_grid_pesudo_spectral_stable2d(float **vp, float dx, float dz, float dt, int nx,
                                             int nz);
  /**
   * stable analysis for regular grid pesudo spectral method by 3D second order
   * acoustic wave equation
   *  stable analysis for regular grid pesudo spectral method by 3D
   * second order acoustic wave equation
   * @param float ***vp: velocity model
   * @param float dx: the sampling interval along x dimension
   * @param float dy: the sampling interval along y dimension
   * @param float dz: the sampling interval along z dimension
   * @param float dt: the sampling interval along t dimension
   * @param int nx: the sample number along x dimension
   * @param int ny: the sample number along y dimension
   * @param int nz: the sample number along z dimension
   * @return {*}
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-06-09 15:00:17
   */
  void regular_grid_pesudo_spectral_stable3d(float ***vp, float dx, float dy, float dz, float dt,
                                             int nx, int ny, int nz);
  /******************** wavelet ********************/
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
  void riker(float *src, float A, float fpeak, float delay, float dt, int nt);
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
  void iriker(float *src, float A, float fpeak, float delay, float dt, int nt);
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
  void gauss(float *src, float A, float fpeak, float delay, float dt, int nt);
  /**
   *  First derivative of Gaussian function
   *  First derivative of Gaussian function
   * @param float *src: First derivative of Gaussian function
   * @param float A: peak amplitude
   * @param float fpeak: peak frequency
   * @param float delay: delay time
   * @param float dt: the time step
   * @param int nt: umber of samples along t
   * @LastEditors: LaurentCheung
   * @Version: 1.0
   */
  void fdgauss(float *src, float A, float fpeak, float delay, float dt, int nt);
  /******************** grid parameters ********************/
  /**
   * calculate the Lame constants from P- and S-wave velocities and density
   * @param float *lambda: Lame constants, lambda
   * @param float *miu: Lame constants, miu
   * @param float *vp: P-wave velocity
   * @param float *vs: S-wave velocity
   * @param float *rho: density
   * @param int np: the length of pointer
   * @return float *lambda: Lame constants, lambda, length: np
   * @return float *miu: Lame constants, miu, length: np
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-08-01 10:29:15
   */
  void lameConstants(float *lambda, float *miu, float *vp, float *vs, float *rho, int np);
  /**
   * calculate 2D elastic coefficients from
   * P- and S-wave velocities, density
   * @param float *C11: elastic coefficient
   * @param float *C13: elastic coefficient
   * @param float *C55: elastic coefficient
   * @param float *vp: P-wave velocity
   * @param float *vs: S-wave velocity
   * @param float *rho: density
   * @param int np: the length of pointers
   * @return float *Cij: elastic coefficient
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-08-01 10:37:34
   */
  void elasticCoefficients2d(float *C11, float *C13, float *C55, float *vp, float *vs, float *rho,
                             int np);
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
   * @LastEditors: PingMin Cheung
   * @Date: 2024-08-01 10:37:34
   */
  void vti_elasticCoefficients2d(float *C11, float *C13, float *C33, float *C55, float *vp,
                                 float *vs, float *rho, float *epsilon, float *delta, int np);
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
   * @LastEditors: PingMin Cheung
   * @Date: 2024-08-01 10:37:34
   */
  void tti_elasticCoefficients2d(float *C11, float *C13, float *C15, float *C33, float *C35,
                                 float *C55, float *vp, float *vs, float *rho, float *epsilon,
                                 float *delta, float *theta, int np);
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
   * @LastEditors: PingMin Cheung
   * @Date: 2024-08-01 10:37:34
   */
  void vti_elasticCoefficients3d(float *C11, float *C12, float *C13, float *C33, float *C55,
                                 float *C66, float *vp, float *vs, float *rho, float *epsilon,
                                 float *delta, float *gama, int np);
  /**
   * calculate the model's parameters on the center of 2d staggered grid
   * @param float **mod_in : input model parameters
   * @param float **mod_out : output model parameters
   * @param int nz: the number of grid points in the z direction
   * @param int nx: the number of grid points in the x direction
   * @param int type: type of grid
   *                  1 for inverse of rho in x direction,
   *                  2 for inverse of rho in z direction,
   *                  3 for miu in the center of staggered grid
   * @return float **mod_out : output model parameters
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-07-31 15:11:49
   */
  void staggered_grid2d_parameter(float **mod_in, float **mod_out, int nz, int nx, int type);
  /**
   * calculate the model's parameters on the center of 3d staggered grid
   * @param float ***mod_in : input model parameters
   * @param float ***mod_out : output model parameters
   * @param int nz: the number of grid points in the z direction
   * @param int nx: the number of grid points in the x direction
   * @param int ny: the number of grid points in the y direction
   * @return float ***mod_out : output model parameters
   * @param int type: type of grid
   *                  1 for inverse of rho in x direction,
   *                  2 for inverse of rho in y direction,
   *                  3 for inverse of rho in z direction,
   *                  4 for miu in the center of staggered grid at XOY plane
   *                  5 for miu in the center of staggered grid at XOZ plane
   *                  6 for miu in the center of staggered grid at YOZ plane
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-07-31 15:11:49
   */
  void staggered_grid3d_parameter(float ***mod_in, float ***mod_out, int nz, int nx, int ny,
                                  int type);
  /**
   * calculate the model's parameters on the center of 2d rotated staggered grid
   * @param float **mod_in : input model parameters
   * @param float **mod_out : output model parameters
   * @param int nz: the number of grid points in the z direction
   * @param int nx: the number of grid points in the x direction
   * @return float **mod_out : output model parameters
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-07-31 15:11:49
   */
  void rotated_staggered_grid2d_parameter(float **mod_in, float **mod_out, int nz, int nx);
  /**
   * calculate the model's parameters on the center of 3d rotated staggered grid
   * @param float ***mod_in : input model parameters
   * @param float ***mod_out : output model parameters
   * @param int nz: the number of grid points in the z direction
   * @param int nx: the number of grid points in the x direction
   * @param int ny: the number of grid points in the y direction
   * @return float ***mod_out : output model parameters
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-07-31 15:11:49
   */
  void rotated_staggered_grid3d_parameter(float ***mod_in, float ***mod_out, int nz, int nx,
                                          int ny);
  /******************** finite difference scheme ********************/
  /**
   * the 8th order accuracy finite difference to calculate the first derivative
   *  the 8th order accuracy finite difference to calculate the first
   * derivative
   * @param float *dat_fd: out 1D array
   * @param float *dat: input 1D array
   * @param float d: the samping interval of data
   * @param int n: number of samples in 1st dimension
   * @return {*}
   * @LastEditors: LaurentCheung
   * @Version: 1.0
   */
  void derivate1_fd1d(float *dat_fd, float *dat, float d, int n);
  /**
   * the 8th order accuracy finite difference to calculate the first derivative
   *  the 8th order accuracy finite difference to calculate the first
   * derivative
   * @param float **dat_fd: out 2D array
   * @param float **dat: input 2D array
   * @param float d: the interval of the calculation direction
   * @param int nz: number of samples in 1st dimension
   * @param int nx: number of samples in 2nd dimension
   * @param int dim: the direction for calculating directive, =1: calcualte
   * derivative along z (vertical), =2: calculate derivative along x (horizontal);
   * @return {*}
   * @LastEditors: LaurentCheung
   * @Version: 1.0
   */
  void derivate1_fd2d(float **dat_fd, float **dat, float d, int nz, int nx, int dim);
  /**
   * the 8th order accuracy finite difference to calculate the first derivative
   *  the 8th order accuracy finite difference to calculate the first
   * derivative
   * @param float ***dat_fd: out 2D array
   * @param float ***dat: input 2D array
   * @param float d: the interval of the calculation direction
   * @param int nz: number of samples in 1st dimension
   * @param int nx: number of samples in 2nd dimension
   * @param int ny: number of samples in 3rd dimension
   * @param int dim: the direction for calculating directive, =1: calcualte
   * derivative along z (vertical), =2: calculate derivative along x
   * (horizontal),=3: calculate derivative along y (horizontal);
   * @return {*}
   * @LastEditors: LaurentCheung
   * @Version: 1.0
   */
  void derivate1_fd3d(float ***dat_fd, float ***dat, float d, int nz, int nx, int ny, int dim);
  /**
   *
   * @param float **gx
   * @param float **gz
   * @param float **dat
   * @param float dz
   * @param float dx
   * @param int nz
   * @param int nx
   * @return {*}
   * @LastEditors: LaurentCheung
   * @Version: 1.0
   */
  void gradient2d(float **gz, float **gx, float **dat, float dz, float dx, int nz, int nx);
  /**
   *
   * @param float ***gz:
   * @param float ***gx:
   * @param float ***gy:
   * @param float ***dat:
   * @param float dz:
   * @param float dx:
   * @param float dy:
   * @param int nz:
   * @param int nx:
   * @param int ny:
   * @return {*}
   * @LastEditors: LaurentCheung
   * @Version: 1.0
   */
  void gradient3d(float ***gz, float ***gx, float ***gy, float ***dat, float dz, float dx, float dy,
                  int nz, int nx, int ny);
  /**
   *
   * @param float **div:
   * @param float **Fx:
   * @param float **Fz:
   * @param float dz:
   * @param float dx:
   * @param int nz:
   * @param int nx:
   * @return {*}
   * @LastEditors: LaurentCheung
   * @Version: 1.0
   */
  void divergence2d(float **div, float **Fx, float **Fz, float dz, float dx, int nz, int nx);
  /**
   *
   * @param float **laplace:
   * @param float **dat:
   * @param float dz:
   * @param float dx:
   * @param int nz:
   * @param int nx:
   * @return {*}
   * @LastEditors: LaurentCheung
   * @Version: 1.0
   */
  void laplace2d(float **laplace, float **dat, float dz, float dx, int nz, int nx);
  /**
   *
   * @param float **q
   * @param float **p
   * @param float D80
   * @param float D81
   * @param float D82
   * @param float D83
   * @param float D84
   * @param float invdxdx
   * @param float invdzdz
   * @param float beta
   * @param int Nz
   * @param int Nx
   * @param int N
   * @param int maxit
   * @return {*}
   * @LastEditors: LaurentCheung
   * @Version: 1.0
   */
  void poissonsolver2d(float **p, float **f, float D80, float D81, float D82, float D83, float D84,
                       float invdxdx, float invdzdz, float beta, int Nz, int Nx, int N, int maxit);
#ifdef __cplusplus
}
#endif