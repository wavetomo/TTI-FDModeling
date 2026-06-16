
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
#include <su.h>
/**/
#include <header.h>
#include <segy.h>
#include <tapehdr.h>
/************************************/
/**
 * get the number of traces from a su format file
 * @param FILE *fp: FILE pointer to a su format file for reading
 * @return unsigned long int: the number of traces
 * @Author: PingMin Cheung
 * @Date: 2024-10-20 15:56:06
 */
unsigned long int fgetntr(FILE *fp)
{
  unsigned long int ntr;
  int index, offset;
  unsigned short ns, temp;
  size_t size_ns;
  size_ns = sizeof(unsigned short);
  index = getindex("ns");
  offset = tapehdr[index].offs;
  if (fp != NULL)
  {
    if (fseek(fp, 0, SEEK_SET) == 0)
    {
      if (efseek(fp, offset, SEEK_CUR) == 0)
      {
        ntr = 0;
        temp = 0;
        while (efread(&ns, 1, size_ns, fp) == size_ns)
        {
          if (ns == 0)
          {
            ntr = 0;
            warn("<warning>: the number of samples is an illegal value (%d) in trace (%d).", ns, ntr + 1);
            break;
          }
          else
          {
            ntr++;
            if (temp != 0 && temp != ns)
            {
              warn("<warning>: the number of samples  (%d) in trace (%d) is different to previous "
                   "trace (%d).",
                   ns, ntr, temp);
            }
            temp = ns;
            if (efseek(fp, HDRBYTES + ns * FSIZE - size_ns, SEEK_CUR) != 0)
            {
              warn("<warning>: failed to seek the next trace header position.");
            }
          }
        }
      }
      else
      {
        warn("<warning>: failed to seek the first trace header 'ns' position.");
      }
    }
    else
    {
      ntr = 0;
      warn("<warning>: failed to seek the beginning position of file.");
    }
  }
  else
  {
    ntr = 0;
    warn("<warning>: the file is not opened for reading.");
  }
  return ntr;
}

/**
 * get the number of shots from a su format file
 * @param FILE *fp: FILE pointer to a su format file for reading
 * @return int: the number of shots
 * @Author: PingMin Cheung
 * @Date: 2024-10-20 15:56:06
 */
int checkShotsGather(FILE *fp)
{
  segy tr, tp;
  int itr, nShot;
  memset(&tr, 0, sizeof(segy));
  memset(&tp, 0, sizeof(segy));
  itr = 0;
  nShot = 1;
  if (fp != NULL)
  {
    if (fseek(fp, 0, SEEK_SET) == 0)
    {
      while (efread(&tr, 1, HDRBYTES, fp) == HDRBYTES)
      {
        itr++;
        if (tr.ns == 0)
        {
          nShot = 0;
          warn("<warning>: the number of samples is an illegal value (%d) in shot "
               "(%d) trace (%d).",
               tr.ns, nShot, itr);
          break;
        }
        else
        {
          if (tp.ns != 0)
          {
            if (tr.sx == tp.sx && tr.sy == tp.sy)
            {
              if (tr.fldr != tp.fldr)
              {
                nShot = 0;
                warn("<warning>: the original field record number (%d) of trace "
                     "(%d) is different to previous traces in shot (%d).",
                     tr.fldr, itr, nShot);
                break;
              }
              if (tr.ntr != tp.ntr)
              {
                nShot = 0;
                warn("<warning>: the number of traces (%d) of trace (%d) is "
                     "different to previous traces in shot (%d).",
                     tr.ntr, itr, nShot);
                break;
              }
              if (tr.ns != tp.ns)
              {
                warn("<warning>: the number of samples (%d) of trace (%d) is "
                     "different to previous traces in shot (%d).",
                     tr.ns, itr, nShot);
              }
              if (tr.dt != tp.dt)
              {
                warn("<warning>: the sampling interval (%d) of trace (%d) is "
                     "different to previous traces in shot (%d).",
                     tr.dt, itr, nShot);
              }
            }
            else
            {
              itr = 0;
              nShot++;
            }
          }
          memcpy(&tp, &tr, HDRBYTES);
        }
        if (efseek(fp, tr.ns * FSIZE, SEEK_CUR) != 0)
        {
          nShot = 0;
          warn("<warning>: failed to seek the end position of traces (%d).", itr);
          break;
        }
      }
    }
    else
    {
      nShot = 0;
      warn("<warning>: failed to seek the beginning position.");
    }
  }
  else
  {
    nShot = 0;
    warn("<warning>: the file is not opened for reading.");
  }
  return nShot;
}
/**
 * get the number of traces from a su format file
 * @param char *fileName: file name
 * @return int ntr: the number of traces
 *                  0 for failure else for success
 * @Author: PingMin Cheung
 * @Date: 2023-03-29 20:04:05
 */
int numberTraces(char *fileName)
{
  FILE *fp = NULL;
  int ntr;
  fp = cfopen(fileName, "rb");
  ntr = fgetntr(fp);
  cfclose(fp, fileName);
  return ntr;
}
/**
 * get the number of shots from a su format file
 * @param char *fileName: the fileName of su format file
 * @return int numShots: the number of shots
 *                       0 for failure else for success
 * @Author: PingMin Cheung
 * @Date: 2023-04-12 19:47:37
 */
int numberShots(char *fileName)
{
  FILE *fp = NULL;
  int numShots;
  fp = cfopen(fileName, "rb");
  numShots = checkShotsGather(fp);
  cfclose(fp, fileName);
  return numShots;
}