/*
 * @Description:
 * @Version: 1.0
 * @Autor: LaurentCheung
 * @Date: 2022-06-06 00:03:15
 * @LastEditTime: 2024-06-11 16:07:24
 */
/**/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/**/
#include <cwp.h>
#include <par.h>
//
char *sdoc[] = {
    /**/
    "CUTBIN2D - cut a binary file with a constant value               ",
    "cutbin2d <stdin n1= n2= fx1= fx2= nx1= nx2= >stdout              ",
    "Required parameter:                                              ",
    "n1       number of samples in 1st (fast) dimension               ",
    "n2       number of samples in 2nd (slow) dimension               ",
    "fx1      first number of outputs samples in 1st (fast) dimension ",
    "fx2      second number of output samples in 2nd (slow) dimension ",
    "nx1      number of output samples in 1st (fast) dimension        ",
    "nx2      number of output samples in 2nd (slow) dimension        ",
    "Example:                                                         ",
    NULL};
int main(int argc, char *argv[])
{
  /* Initialize */
  initargs(argc, argv);
  requestdoc(0); /* stdin not used */
  /* Make sure stdout is a file or pipe */
  int n1, n2, fx1, fx2, nx1, nx2;
  if (!getparint("n1", &n1))
  {
    err("<error>: must specify number of samples in 1st (fast) dimension!");
  }
  if (!getparint("n2", &n2))
  {
    err("<error>: must specify number of samples in 2nd (slow) dimension!");
  }
  if (!getparint("fx1", &fx1))
  {
    fx1 = 0;
    warn("<warning>: didn't specify first number of outputs samples in 1st "
         "(fast) dimension! The program will set it as %d!",
         fx1);
  }
  if (!getparint("fx2", &fx2))
  {
    fx2 = 0;
    warn("<warning>: didn't specify second number of outputs samples in 2nd "
         "(slow) dimension! The program will set it as %d!",
         fx2);
  }
  if (!getparint("nx1", &nx1))
  {
    nx1 = n1;
    warn("<warning>: didn't specify number of output samples in 1st (fast) "
         "dimension! The program will set it as %d!",
         nx1);
  }
  if (!getparint("nx2", &nx2))
  {
    nx2 = n2;
    warn("<warning>: didn't specify number of output samples in 2nd (slow) "
         "dimension! The program will set it as %d!",
         nx2);
  }
  int i1, i2;
  float **dat = NULL, **dat_cut = NULL;
  dat = alloc2float(n1, n2);
  dat_cut = alloc2float(nx1, nx2);
  efread(dat[0], FSIZE, n1 * n2, stdin);
  for (i1 = 0; i1 < nx2; i1++)
  {
    for (i2 = 0; i2 < nx1; i2++)
    {
      dat_cut[i1][i2] = dat[i1 + fx2][i2 + fx1];
    }
  }
  stdout->_IO_write_ptr;
  efwrite(dat_cut[0], FSIZE, nx1 * nx2, stdout);
  fflush(stdin);
  fclose(stdin);
  fflush(stdout);
  fclose(stdout);
  free2float(dat);
  free2float(dat_cut);
  return (CWP_Exit());
}