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
/* types initialize */
/**
 * initialize 2D source structure
 * @param Source2D *source: 2D source structure pointer
 * @Author: PingMin Cheung
 * @Date: 2024-08-14 14:18:10
 */
void source2d_init(Source2D *source)
{
  source->ntrace = 0;
  source->ntsample = 0;
  source->posx = NULL;
  source->posz = NULL;
  source->src = NULL;
}
/**
 * initialize 3D source structure
 * @param Source2D *source: 3D source structure pointer
 * @Author: PingMin Cheung
 * @Date: 2024-08-14 14:18:10
 */
void source3d_init(Source3D *source)
{
  source->ntrace = 0;
  source->ntsample = 0;
  source->posx = NULL;
  source->posy = NULL;
  source->posz = NULL;
}
/**
 * initialize 2D receiver structure
 * @param Receiver2D *receiver: 2D receiver structure pointer
 * @Author: PingMin Cheung
 * @Date: 2024-08-14 14:18:10
 */
void receiver2d_init(Receiver2D *receiver)
{
  receiver->ntrace = 0;
  receiver->ntsample = 0;
  receiver->posx = NULL;
  receiver->posz = NULL;
  receiver->recv = NULL;
}
/**
 * initialize 3D receiver structure
 * @param Receiver3D *receiver: 3D receiver structure pointer
 * @Author: PingMin Cheung
 * @Date: 2024-08-14 14:18:10
 */
void receiver3d_init(Receiver3D *receiver)
{
  receiver->ntrace = 0;
  receiver->ntsample = 0;
  receiver->posx = NULL;
  receiver->posy = NULL;
  receiver->recv = NULL;
}
/* types free */
/**
 * free memery of Source2D
 *  free memery of Source2D
 * @param {Source2D} source: 2D source structure pointer
 *                               contain the source information
 * @LastEditors: LaurentCheung
 * @Version: 1.0
 */
void freeSource2d(Source2D *source)
{
  if (source->posx != NULL)
  {
    free1int(source->posx);
    source->posx = NULL;
  }
  if (source->posz != NULL)
  {
    free1int(source->posz);
    source->posz = NULL;
  }
  if (source->src != NULL)
  {
    free2float(source->src);
    source->src = NULL;
  }
}
/**
 * free memery of Source3D
 *  free memery of Source3D
 * @param {Source3D} source: 3D source structure pointer
 *                               contain the source information
 * @LastEditors: LaurentCheung
 * @Version: 1.0
 */
void freeSource3d(Source3D *source)
{
  if (source->posx != NULL)
  {
    free1int(source->posx);
    source->posx = NULL;
  }
  if (source->posy != NULL)
  {
    free1int(source->posy);
    source->posy = NULL;
  }
  if (source->posz != NULL)
  {
    free1int(source->posz);
    source->posz = NULL;
  }
  if (source->src != NULL)
  {
    free2float(source->src);
    source->src = NULL;
  }
}
/**
 * free memery of Receiver2D
 *  free memery of Receiver2D
 * @param {Receiver2D} receiver: 2D receiver structure pointer
 *                               contain the receiver information
 * @LastEditors: LaurentCheung
 * @Version: 1.0
 */
void freeReceiver2d(Receiver2D *receiver)
{
  if (receiver->posx != NULL)
  {
    free1int(receiver->posx);
    receiver->posx = NULL;
  }
  if (receiver->posz != NULL)
  {
    free1int(receiver->posz);
    receiver->posz = NULL;
  }
  if (receiver->recv != NULL)
  {
    free2float(receiver->recv);
    receiver->recv = NULL;
  }
}
/**
 * free memery of Receiver3D
 *  free memery of Receiver3D
 * @param {Receiver3D} receiver: 3D receiver structure pointer
 *                               contain the receiver information
 * @LastEditors: LaurentCheung
 * @Version: 1.0
 */
void freeReceiver3d(Receiver3D *receiver)
{
  if (receiver->posx != NULL)
  {
    free1int(receiver->posx);
    receiver->posx = NULL;
  }
  if (receiver->posy != NULL)
  {
    free1int(receiver->posy);
    receiver->posy = NULL;
  }
  if (receiver->posz != NULL)
  {
    free1int(receiver->posz);
    receiver->posz = NULL;
  }
  if (receiver->recv != NULL)
  {
    free2float(receiver->recv);
    receiver->recv = NULL;
  }
}
/* types copy */
/**
 * copy 2D receiver structure
 * @param Receiver2D *dest: 2D receiver structure pointer
 * @param Receiver2D *src: 2D receiver structure pointer
 * @return *
 * @Author: PingMin Cheung
 * @Date: 2024-09-01 18:03:07
 */
void receiver2d_copy(Receiver2D *dest, Receiver2D *src)
{ // copy src to dest
  dest->ntrace = src->ntrace;
  dest->ntsample = src->ntsample;
  dest->posx = alloc1int(dest->ntrace);
  dest->posz = alloc1int(dest->ntrace);
  dest->recv = alloc2float(dest->ntsample, dest->ntrace);
  memcpy(dest->posx, src->posx, ISIZE * dest->ntrace);
  memcpy(dest->posz, src->posz, ISIZE * dest->ntrace);
  memcpy(dest->recv[0], src->recv[0], FSIZE * dest->ntrace * dest->ntsample);
}
/**
 * copy 3D receiver structure
 * @param Receiver2D *dest: 2D receiver structure pointer
 * @param Receiver2D *src: 2D receiver structure pointer
 * @return *
 * @Author: PingMin Cheung
 * @Date: 2024-09-01 18:03:07
 */
void receiver3d_copy(Receiver3D *dest, Receiver3D *src)
{ // copy src to dest
  dest->ntrace = src->ntrace;
  dest->ntsample = src->ntsample;
  dest->posx = alloc1int(dest->ntrace);
  dest->posz = alloc1int(dest->ntrace);
  dest->posy = alloc1int(dest->ntrace);
  dest->recv = alloc2float(dest->ntsample, dest->ntrace);
  memcpy(dest->posx, src->posx, ISIZE * dest->ntrace);
  memcpy(dest->posz, src->posz, ISIZE * dest->ntrace);
  memcpy(dest->posy, src->posy, ISIZE * dest->ntrace);
  memcpy(dest->recv[0], src->recv[0], FSIZE * dest->ntrace * dest->ntsample);
}