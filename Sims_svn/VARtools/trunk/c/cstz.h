#ifndef __CSTZ_H__
#define __CSTZ_H__
   #include <math.h>
   #include "tzmatlab.h"
   #include "mathlib.h"

   #if defined ( CSMINWEL_OPTIMIZATION )
      void fn_gradcd(double *g, double *x, int n, double grdh,
                     double (*fcn)(double *x, int n, double **args, int *dims),
                     double **args, int *dims);

      void fn_hesscd(double *H, double *x, int n, double grdh,
                     double (*fcn)(double *x, int n, double **args, int *dims),
                     double **args, int *dims);
   #elif defined ( IMSL_OPTIMIZATION )
      void fn_gradcd(double *g, double *x, int n, double grdh,
                  double fcn(int n, double *x));
      void fn_hesscd(double *H, double *x, int n, double grdh,
                  double fcn(int n, double *x));
   #endif

   //=== For the conjugate gradient method I or II
   void gradcd_gen(double *g, double *x, int n, double (*fcn)(double *x, int n), double *grdh, double f0);
   void gradfd_gen(double *g, double *x, int n, double (*fcn)(double *x, int n), double *grdh, double f0);


   int next_permutation(int *first, int *last);


   //void fn_ergodp(double **aop, int *aod, mxArray *cp);
   void fn_cumsum(double **aos_v, int *aods_v, double *v, int d_v);
   int fn_cumsum_int(int *x_v, const int d_x_v);
   double fn_cumsum_lf(double *x_v, const int d_x_v);
   double fn_mean(const double *a_v, const int _n);


   //=== For sorting according to x_dv.
   void fn_SetBaseArrayForComp(TSdvector *x_dv);
   int fn_compare(const void *i1, const void *i2);
   int fn_compare2(const void *i1, const void *i2);
   //=== Normalization for VARs.
   void fn_wznormalization(TSdvector *wznmlz_dv, TSdmatrix *A0draw_dm, TSdmatrix *A0peak_dm);


   //=== Handling under or over flows with log values.
   int fn_update_logofsum(int N, double ynew, double *Y_N_dp, double *y_Nmax_dp);
   double fn_replace_logofsumsbt(double *yold, double _a, double ynew, double _b);


   //=== Special functions.
   double fn_chi2inv(double p, double df);
        //p = int_{0}^{\infty} chi2pdf(t, df) dt
   double fn_betainv(double p, double _alpha, double _beta);
        //p = int_{0}^{\infty} betapdf(t, _alpha, _beta) dt where betapdf(t,_alpha,_beta) \propt t^{_alpha-1}*(1-t)^(_beta-1}.
   double fn_gammalog(double x);
        //log gamma (x) where gamma(n+1) = n! and gamma(x) = int_0^{\infty} e^{-t} t^{x-1} dt.
   double fn_betalog(double x, double y);
        //log beta(x, y) where beta(x, y) = gamma(x)*gamm(y)/gamma(x+y).


   //---------------------------- Old Interface. ---------------------
   double gammalog(double x);
        //log gamma (x) where gamma(n+1) = n! and gamma(x) = int_0^{\infty} e^{-t} t^{x-1} dt.


   //---------------------------- New functions with my own matrix notations. ---------------------
   double MinVector_lf(TSdvector *x_dv);
#endif
