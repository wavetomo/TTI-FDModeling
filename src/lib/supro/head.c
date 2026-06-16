/*
 * @description:
 * @Author: PingMin Cheung
 * @Date: 2024-04-09 15:22:32
 * @LastEditors: Pingmin Cheung
 * git config user.email & please set dead value or install git & please set dead value or install
 * git
 * @LastEditTime: 2026-05-09 10:40:15
 * @FilePath: /src/supro/lib/head.c
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
/**
 * apply scalar to the value of head for reading su format file
 * @param float key: the value of head
 * @param int scalar: scalar to be applied
 * @return int val: the value of head after applying scalar
 * @Author: PingMin Cheung
 * @Date: 2024-04-09 15:50:26
 */
int head_scalar_read(float key, int scalar)
{
  int val, temp;
  if (scalar == 0)
  {
    err("<error> : scalar must not be %d!", scalar);
    val = key;
  }
  else
  {
    temp = SGN(scalar) * log10((double)abs(scalar));
    val = key * pow(10.0, temp);
  }
  return val;
}
/**
 * apply scalar to the value of head for writting su format file
 * @param float key: the value of head
 * @param int scalar: scalar to be applied
 * @return int val: the value of head after applying scalar
 * @Author: PingMin Cheung
 * @Date: 2024-04-09 15:50:26
 */
int head_scalar_write(float key, int scalar) { return head_scalar_read(key, -1 * scalar); }
/**
 * remainder after applying scalar to the value of head for reading su
 * format file
 * @param float key: the value of head
 * @param int scalar: scalar to be applied
 * @return int val: the value of remainder after applying scalar
 * @Author: PingMin Cheung
 * @Date: 2024-04-09 15:50:26
 */
int head_scalar_read_mod(float key, int scalar)
{
  int val, temp;
  float temp_val;
  if (scalar == 0)
  {
    err("<error> : scalar must not be %d!", scalar);
    val = key;
  }
  else
  {
    temp = SGN(scalar) * log10((double)abs(scalar));
    temp_val = key * pow(10.0, temp);
    val = TRUNCATE * (temp_val - (int)(temp_val));
  }
  return val;
}