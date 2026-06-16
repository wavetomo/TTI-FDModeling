/*
 * @description:
 * @Author: PingMin Cheung
 * @Date: 2024-02-26 15:49:28
 * @LastEditTime: 2024-03-17 15:08:51
 * @FilePath: /src/supro/lib/supro.c
 *
 * Copyright (c) 2024 by WaveTomo, All Rights Reserved.
 */

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
/************************************/

segy suzeros(void)
{
  segy tr;
  memset(&tr, 0, sizeof(segy));
  return tr;
}