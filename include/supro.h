#pragma once
#include <su.h>
/**/
#include <segy.h>
/**/
#ifndef TRUNCATE
#define TRUNCATE 100
#endif

#ifdef __cplusplus
extern "C"
{
#endif
  /*********** BINARY IO **********/
  /**
   * open a file with error checking
   * @param char *fileName: file name
   * @param char *mode: open mode
   * @return FILE *fp: file pointer
   *                   NULL for failure else for success
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-06-19 16:27:56
   */
  FILE *cfopen(char *fileName, const char *mode);
  /**
   * close a file with error checking
   * @param FILE *fp: file pointer
   * @param char *fileName: file name
   * @return int value: 1 for success, 0 for failure
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-06-19 16:29:05
   */
  int cfclose(FILE *fp, char *fileName);
  /**
   * get the size of a binary file
   * @param char *fileName: file name
   * @return long int: the size of the file
   *                   0 for failure else for success
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-06-19 16:30:40
   */
  long int sizeofile(char *fileName);
  /**
   * read data from a binary file
   * @param char *fileName: file name
   * @param void *ptr: pointer to the data
   * @param size_t size: size of each element
   * @param size_t nmemb: number of elements
   * @return int value: 1 for success, 0 for failure
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-06-19 16:31:57
   */
  int readBinaryFile(char *fileName, void *ptr, size_t size, size_t nmemb);
  /**
   * write data to a binary file
   * @param char *fileName: file name
   * @param void *ptr: pointer to the data
   * @param size_t size: size of each element
   * @param size_t nmemb: number of elements
   * @return int value: 1 for success, 0 for failure
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-06-19 16:44:39
   */
  int writeBinaryFile(char *fileName, void *ptr, size_t size, size_t nmemb);
  /**
   * write a double percision data to a binary file with float percision
   * @param char *fileName: file name
   * @param double *ptr: pointer to the data
   * @param size_t nmemb: number of elements
   * @return int value: 1 for success, 0 for failure
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-06-19 16:44:39
   */
  int writeBinaryFile_d2f(char *fileName, double *ptr, size_t nmemb);
  /**
   * @param char *prefix: prefix of output file name
   * @param float **ptr: pointer to the wavefield data
   * @param float dt: time step
   * @param int Nx: the number of samples in x direction
   * @param int Nz: the number of samples in z direction
   * @param int nlayer: the thickness of absorbing boundary layers
   * @param int it: time index of wavefield
   * @param int savew: whether to save wavefield with absorbing boundary
   *                   1 for yes, 0 for no
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-10-16 14:49:25
   */
  void snapshot2doutput(char *prefix, float **ptr, float dt, int Nx, int Nz, int nlayer, int it,
                        int savew);
  /*********** SU IO **********/
  segy suzeros(void);
  /**
   * get the number of traces from a su format file
   * @param FILE *fp: FILE pointer to a su format file for reading
   * @return unsigned long int: the number of traces
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-10-20 15:56:06
   */
  unsigned long int fgetntr(FILE *fp);
  /**
   * get the number of shots from a su format file
   * @param FILE *fp: FILE pointer to a su format file for reading
   * @return int: the number of shots
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-10-20 15:56:06
   */
  int checkShotsGather(FILE *fp);
  /**
   * get the number of traces from a su format file
   * @param char *fileName: file name
   * @return int ntr: the number of traces
   *                  0 for failure else for success
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2023-03-29 20:04:05
   */
  int numberTraces(char *fileName);
  /**
   * get the number of shots from a su format file
   * @param char *fileName: the fileName of su format file
   * @return int numShots: he number of shots
   *                  0 for failure else for success
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2023-04-12 19:47:37
   */
  int numberShots(char *fileName);
  /*********** SU HEAD IO **********/
  /**
   * apply scalar to the value of head for reading su format file
   * @param float key: the value of head
   * @param int scalar: scalar to be applied
   * @return int val: the value of head after applying scalar
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-04-09 15:50:26
   */
  int head_scalar_read(float key, int scalar);
  /**
   * apply scalar to the value of head for writting su format file
   * @param float key: the value of head
   * @param int scalar: scalar to be applied
   * @return int val: the value of head after applying scalar
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-04-09 15:50:26
   */
  int head_scalar_write(float key, int scalar);
  /**
   * remainder after applying scalar to the value of head for reading su
   * format file
   * @param float key: the value of head
   * @param int scalar: scalar to be applied
   * @return int val: the value of remainder after applying scalar
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2024-04-09 15:50:26
   */
  int head_scalar_read_mod(float key, int scalar);
  /**
   *
   * @param int keyword
   * @param int scale
   * @return {*}
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2023-04-15 16:36:22
   */
  int header_scale(int keyword, int scale);
  /**
   *
   * @param int keyword
   * @param int scale
   * @return {*}
   * @Author: PingMin Cheung
   * @LastEditors: PingMin Cheung
   * @Date: 2023-04-15 16:36:18
   */
  int header_scale_mod(int keyword, int scale);
#ifdef __cplusplus
}
#endif