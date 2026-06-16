/************************************/
#include <supro.h>
/**/
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/**/
#include <cwp.h>
#include <par.h>
/************************************/
/**
 * open a file with error checking
 * @param char *fileName: file name
 * @param char *mode: open mode
 * @return FILE *fp: file pointer
 *                   NULL for failure else for success
 * @Author: PingMin Cheung
 * @Date: 2024-06-19 16:27:56
 */
FILE *cfopen(char *fileName, const char *mode)
{
  FILE *fp = NULL;
  fp = fopen(fileName, mode);
  if (fp == NULL)
  {
    err("<error>: failed to close file: %s!", fileName);
  }
  return fp;
}
/**
 * close a file with error checking
 * @param FILE *fp: file pointer
 * @param char *fileName: file name
 * @return int value: 1 for success, 0 for failure
 * @Author: PingMin Cheung
 * @Date: 2024-06-19 16:29:05
 */
int cfclose(FILE *fp, char *fileName)
{
  if (fclose(fp) == EOF)
  {
    warn("<warning>: failed to close file: %s!", fileName);
    return 0;
  }
  else
  {
    return 1;
  }
}
/**
 * get the size of a binary file
 * @param char *fileName: file name
 * @return long int: the size of the file
 *                   0 for failure else for success
 * @Author: PingMin Cheung
 * @Date: 2024-06-19 16:30:40
 */
long int sizeofile(char *fileName)
{
  FILE *fp = NULL;
  long int len;
  len = 0;
  fp = cfopen(fileName, "r");
  if (fp != NULL)
  {
    if (fseek(fp, 0, SEEK_END) != 0)
    {
      err("<error>: failed to seek file: %s!", fileName);
    }
    else
    {
      len = ftell(fp);
      if (len == -1L)
      {
        len = 0;
        err("<error>: failed to get file size: %s!", fileName);
      }
      else
      {
        warn("<warning>: the size of the file: %s is %f Mb!", fileName, 1.0 * len / (1024.0 * 1024.0));
      }
    }
    cfclose(fp, fileName);
  }
  return len;
}
/**
 * read data from a binary file
 * @param char *fileName: file name
 * @param void *ptr: pointer to the data
 * @param size_t size: size of each element
 * @param size_t nmemb: number of elements
 * @return int value: 1 for success, 0 for failure
 * @Author: PingMin Cheung
 * @Date: 2024-06-19 16:31:57
 */
int readBinaryFile(char *fileName, void *ptr, size_t size, size_t nmemb)
{
  FILE *fp = NULL;
  fp = cfopen(fileName, "r");
  if (fp != NULL)
  {
    if (fread(ptr, size, nmemb, fp) != nmemb)
    {
      warn("<warning>: failed to read file: %s!", fileName);
    }
    cfclose(fp, fileName);
  }
  else
  {
    return 0;
  }
  return 1;
}
/**
 * write data to a binary file
 * @param char *fileName: file name
 * @param void *ptr: pointer to the data
 * @param size_t size: size of each element
 * @param size_t nmemb: number of elements
 * @return int value: 1 for success, 0 for failure
 * @Author: PingMin Cheung
 * @Date: 2024-06-19 16:44:39
 */
int writeBinaryFile(char *fileName, void *ptr, size_t size, size_t nmemb)
{
  FILE *fp = NULL;
  fp = cfopen(fileName, "w");
  if (fp != NULL)
  {
    if (fwrite(ptr, size, nmemb, fp) != nmemb)
    {
      cfclose(fp, fileName);
      err("<error>: failed to write file: %s!", fileName);
      return 0;
    }
    cfclose(fp, fileName);
    return 1;
  }
  else
  {
    return 0;
  }
}
/**
 * write a double percision data to a binary file with float percision
 * @param char *fileName: file name
 * @param double *ptr: pointer to the data
 * @param size_t nmemb: number of elements
 * @return int value: 1 for success, 0 for failure
 * @Author: PingMin Cheung
 * @Date: 2024-06-19 16:44:39
 */
int writeBinaryFile_d2f(char *fileName, double *ptr, size_t nmemb)
{
  int i;
  float *dat = NULL;
  dat = alloc1float(nmemb);
  for (i = 0; i < nmemb; i++)
  {
    dat[i] = (float)ptr[i];
  }
  i = writeBinaryFile(fileName, dat, sizeof(float), nmemb);
  free1float(dat);
  return i;
}
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
 * @Date: 2024-10-16 14:49:25
 */
void snapshot2doutput(char *prefix, float **ptr, float dt, int Nx, int Nz, int nlayer, int it, int savew)
{
  FILE *fp = NULL;
  char fileName[256];
  int ix, iz;
  int nx, nz;
  nx = (savew == 0 ? Nx - 2 * nlayer : Nx);
  nz = (savew == 0 ? Nz - 2 * nlayer : Nz);
  sprintf(fileName, "%s_%05d_%05d_%05dms.bin", prefix, nz, nx, (int)(it * dt * 1000));
  fp = cfopen(fileName, "wb");
  if (fp == NULL)
  {
    warn("<warning>: faild to open file %s\n", fileName);
  }
  else
  {
    if (savew == 0)
    {
      for (ix = 0; ix < nx; ix++)
      {
        if (fwrite(ptr[nlayer + ix] + nlayer, FSIZE, nz, fp) != nz)
        {
          warn("<warning>: faild to write data to file: %s!", fileName);
        }
      }
    }
    else
    {
      if (fwrite(ptr[0], FSIZE, nx * nz, fp) != nx * nz)
      {
        warn("<warning>: faild to write data to file: %s!", fileName);
      }
    }
    cfclose(fp, fileName);
  }
}