#ifndef __CSMINWEL_H__
#define __CSMINWEL_H__

   #define VERBOSE_WARNINGS   // display warnings
   #define VERBOSE_DETOUPUT   // display detailed output

   #include <math.h>
   #include <string.h>
   #include <stdlib.h>
   #include <stdio.h>

   #include "tzmatlab.h"

   void csminwel(double (*fcn)(double *x, int n, double **args, int *dims),
               double *x, int n, double *H, double *gh,
               int (*grad)(double *x, int n, double *g, double **args,
               int *dims), double *fh, double crit, int *itct, int nit,
               int *fcount, int *retcodeh, double **args, int *dims);
   // Alternative but less clear way:  ... (double (*fcn)(double *, int, double **, int *), ...

   int numgrad(double *g, double *x, int n,
               double (*fcn)(double *x, int n, double **args, int *dims),
               double **args, int *dims);

   void csminit(double *fhat, double *xhat, int *fcount, int *retcode,
               double *x0, double f0, double *g, int badg, double *H0, int n,
               double (*fcn)(double *x, int n, double **args, int *dims),
               double **args, int *dims);

   void bfgsi(double *H, double *dg, double *dx, int n, int nn);

   int peakwall(double *g, int retcode, double *x, int n,
               int (*gfcn)(double *x, int n, double *g, double **args, int *dims),
               double (*fcn)(double *x, int n, double **args, int *dims),
               double **args, int *dims);

   double times(double *x, double *y, int n);

   double *mtimes(double *x, double *y, int n, int nn);

   double *mminus(double *x, double *y, int n);

#endif
