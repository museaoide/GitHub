#ifndef __CONGRADMIN_H__
#define __CONGRADMIN_H__

   #include <math.h>

   #include "tzmatlab.h"

   void frprmn(double p[], int n, int *iter, double *fret,
               double (*func)(double [], int), void (*dfunc)(double [], double [], int, double (*func)(double [], int), double, double),
               double ftol, double grdh);

//   #include <string.h>
//   #include <stdlib.h>




//   void csminwel(double (*fcn)(double *x, int n, double **args, int *dims),
//               double *x, int n, double *H, double *gh,
//               int (*grad)(double *x, int n, double *g, double **args,
//               int *dims), double *fh, double crit, int *itct, int nit,
//               int *fcount, int *retcodeh, double **args, int *dims);
//   // Alternative but less clear way:  ... (double (*fcn)(double *, int, double **, int *), ...

//   int numgrad(double *g, double *x, int n,
//               double (*fcn)(double *x, int n, double **args, int *dims),
//               double **args, int *dims);

//   void csminit(double *fhat, double *xhat, int *fcount, int *retcode,
//               double *x0, double f0, double *g, int badg, double *H0, int n,
//               double (*fcn)(double *x, int n, double **args, int *dims),
//               double **args, int *dims);

//   void bfgsi(double *H, double *dg, double *dx, int n, int nn);

//   int peakwall(double *g, int retcode, double *x, int n,
//               int (*gfcn)(double *x, int n, double *g, double **args, int *dims),
//               double (*fcn)(double *x, int n, double **args, int *dims),
//               double **args, int *dims);

//   double times(double *x, double *y, int n);

//   double *mtimes(double *x, double *y, int n, int nn);

//   double *mminus(double *x, double *y, int n);

#endif
