#include "cstz.h"


#if defined( CSMINWEL_OPTIMIZATION )
   #define STPS 6.0554544523933391e-6    /* step size = pow(DBL_EPSILON,1.0/3) */
   void fn_gradcd(double *g, double *x, int n, double grdh,
                  double (*fcn)(double *x, int n, double **args, int *dims),
                  double **args, int *dims) {
      //Outputs:
      //  g: the gradient n-by-1 g (no need to be initialized).
      //Inputs:
      //  grdh: step size.  If ==0.0, then dh is set automatically; otherwise, grdh is taken as a step size, often set as 1.0e-004.
      //  x:  no change in the end although will be added or substracted by dh during the function (but in the end the original value will be put back).

      double dh, fp, fm, tmp, *xp;
      int i;
      for (i=0, xp=x; i<n; i++, xp++, g++) {
         dh = grdh?grdh:fabs(*xp)<1?STPS:STPS*(*xp);
         tmp = *xp;
         *xp += dh;
         dh = *xp - tmp;                   // This increases the precision slightly.
         fp = fcn(x,n,args,dims);
         *xp = tmp - dh;
         fm = fcn(x,n,args,dims);
         *g = (fp-fm)/(2*dh);
         *xp = tmp;                        // Put the original value of x[i] back to x[i] so that the content x[i] is still unaltered.
      }
   }
   #undef STPS

   #define STPS 6.0554544523933391e-6    /* step size = pow(DBL_EPSILON,1.0/3) */
   void fn_hesscd(double *H, double *x, int n, double grdh,
                  double (*fcn)(double *x, int n, double **args, int *dims),
                  double **args, int *dims) {
      double dhi, dhj, f1, f2, f3, f4, tmpi, tmpj, *xpi, *xpj;
      int i, j;
      for (i=0, xpi=x; i<n; i++, xpi++) {
         dhi = grdh?grdh:fabs(*xpi)<1?STPS:STPS*(*xpi);
         tmpi = *xpi;
         for (j=i, xpj=x+i; j<n; j++, xpj++)
            if (i==j) {
               /* f2 = f3 when i = j */
               f2 = fcn(x,n,args,dims);

               /* this increases precision slightly */
               *xpi += dhi;
               dhi = *xpi - tmpi;

               /* calculate f1 and f4 */
               *xpi = tmpi + 2*dhi;
               f1 = fcn(x,n,args,dims);
               *xpi = tmpi - 2*dhi;
               f4 = fcn(x,n,args,dims);

               /* diagonal element */
               H[i*(n+1)] = (f1-2*f2+f4)/(4*dhi*dhi);

               /* reset to intial value */
               *xpi = tmpi;
            } else {
               dhj = grdh?grdh:fabs(*xpj)<1?STPS:STPS*(*xpj);
               tmpj = *xpj;

               /* this increases precision slightly */
               *xpi += dhi;
               dhi = *xpi - tmpi;
               *xpj += dhj;
               dhj = *xpj - tmpj;

               /* calculate f1, f2, f3 and f4 */
               *xpi = tmpi + dhi;
               *xpj = tmpj + dhj;
               f1 = fcn(x,n,args,dims);
               *xpi = tmpi - dhi;
               f2 = fcn(x,n,args,dims);
               *xpi = tmpi + dhi;
               *xpj = tmpj - dhj;
               f3 = fcn(x,n,args,dims);
               *xpi = tmpi - dhi;
               f4 = fcn(x,n,args,dims);

               /* symmetric elements */
               H[i+j*n] = H[j+i*n] = (f1-f2-f3+f4)/(4*dhi*dhj);

               /* reset to intial values */
               *xpi = tmpi;
               *xpj = tmpj;
            }
      }
   }
   #undef STPS
#elif defined( IMSL_OPTIMIZATION )
   #define STPS 6.0554544523933391e-6    /* step size = pow(DBL_EPSILON,1.0/3) */
   void fn_gradcd(double *g, double *x, int n, double grdh,
                  double fcn(int n, double *x) // IMSL
                  //void NAG_CALL fcn(Integer n,double x[],double *f,double g[],Nag_Comm *comm)
                  ) {
      //Outputs:
      //  g: the gradient n-by-1 g (no need to be initialized).
      //Inputs:
      //  grdh: step size.  If ==0.0, then dh is set automatically; otherwise, grdh is taken as a step size, often set as 1.0e-004.
      //  x:  no change in the end although will be added or substracted by dh during the function (but in the end the original value will be put back).

      double dh, fp, fm, tmp, *xp;
      int i;
      for (i=0, xp=x; i<n; i++, xp++, g++) {
         dh = grdh?grdh:fabs(*xp)<1?STPS:STPS*(*xp);
         tmp = *xp;
         *xp += dh;
         dh = *xp - tmp;                   // This increases the precision slightly.
         fp = fcn(n,x); // IMSL
         //fcn(n,x,&fp,NULL,NULL); /* NAG */
         *xp = tmp - dh;
         fm = fcn(n,x); // IMSL
         //fcn(n,x,&fm,NULL,NULL);
         *g = (fp-fm)/(2*dh);
         *xp = tmp;                        // Put the original value of x[i] back to x[i] so that the content x[i] is still unaltered.
      }
   }
   #undef STPS

   #define STPS 6.0554544523933391e-6    /* step size = pow(DBL_EPSILON,1.0/3) */
   void fn_hesscd(double *H, double *x, int n, double grdh,
                  double fcn(int n, double *x) // IMSL
                  //void NAG_CALL fcn(Integer n,double x[],double *f,double g[],Nag_Comm *comm)
                  ) {
      double dhi, dhj, f1, f2, f3, f4, tmpi, tmpj, *xpi, *xpj;
      int i, j;
      for (i=0, xpi=x; i<n; i++, xpi++) {
         dhi = grdh?grdh:fabs(*xpi)<1?STPS:STPS*(*xpi);
         tmpi = *xpi;
         for (j=i, xpj=x+i; j<n; j++, xpj++)
            if (i==j) {
               /* f2 = f3 when i = j */
               f2 = fcn(n,x); // IMSL
               //fcn(n,x,&f2,NULL,NULL);

               /* this increases precision slightly */
               *xpi += dhi;
               dhi = *xpi - tmpi;

               /* calculate f1 and f4 */
               *xpi = tmpi + 2*dhi;
               f1 = fcn(n,x); // IMSL
               //fcn(n,x,&f1,NULL,NULL);
               *xpi = tmpi - 2*dhi;
               f4 = fcn(n,x); /* IMSL */
               //fcn(n,x,&f4,NULL,NULL);

               /* diagonal element */
               H[i*(n+1)] = (f1-2*f2+f4)/(4*dhi*dhi);

               /* reset to intial value */
               *xpi = tmpi;
            } else {
               dhj = grdh?grdh:fabs(*xpj)<1?STPS:STPS*(*xpj);
               tmpj = *xpj;

               /* this increases precision slightly */
               *xpi += dhi;
               dhi = *xpi - tmpi;
               *xpj += dhj;
               dhj = *xpj - tmpj;

               /* calculate f1, f2, f3 and f4 */
               *xpi = tmpi + dhi;
               *xpj = tmpj + dhj;
               f1 = fcn(n,x); // IMSL
               //fcn(n,x,&f1,NULL,NULL);
               *xpi = tmpi - dhi;
               f2 = fcn(n,x); // IMSL
               //fcn(n,x,&f2,NULL,NULL);
               *xpi = tmpi + dhi;
               *xpj = tmpj - dhj;
               f3 = fcn(n,x); // IMSL
               //fcn(n,x,&f3,NULL,NULL);
               *xpi = tmpi - dhi;
               f4 = fcn(n,x); // IMSL
               //fcn(n,x,&f4,NULL,NULL);

               /* symmetric elements */
               H[i+j*n] = H[j+i*n] = (f1-f2-f3+f4)/(4*dhi*dhj);

               /* reset to intial values */
               *xpi = tmpi;
               *xpj = tmpj;
            }
      }
   }
   #undef STPS
#endif



//-------------------------------
//Modified from fn_gradcd() in cstz.c for the minimization problem for  the conjugate gradient method I or II
//-------------------------------
#define STPS 1.0e-04    // 6.0554544523933391e-6 step size = pow(DBL_EPSILON,1.0/3)
void gradcd_gen(double *g, double *x, int n, double (*fcn)(double *x, int n), double *grdh, double f0) {
   //Outputs:
   //  g: the gradient n-by-1 g (no need to be initialized).
   //Inputs:
   //  x: the vector point at which the gradient is evaluated.  No change in the end although will be added or substracted by dh during the function (but in the end the original value will be put back).
   //  n: the dimension of g or x.
   //  fcn(): the function for which the gradient is evaluated
   //  grdh: step size.  If NULL, then dh is set automatically; otherwise, grdh is taken as a step size, often set as 1.0e-004.
   //  f0: the value of (*fcn)(x).   NOT used in this function except dealing with the boundary (NEARINFINITY) for the
   //    minimization problem, but to be compatible with a genral function call where, say, gradfw_gen() and cubic
   //    interpolation of central difference method will use f0.

   double dh, dhi, dh2i, fp, fm, tmp, *xp;
   int i;
   if (grdh) {
      dh2i = (dhi=1.0/(dh=*grdh))/2.0;
      for (i=0, xp=x; i<n; i++, xp++, g++) {
         tmp = *xp;
         *xp += dh;
         //The following statement is bad because dh does not get reset at the beginning of the loop and thus may get changed continually within the loop.
         //  dh = *xp - tmp;                   // This increases the precision slightly.
         fp = fcn(x, n); //For frprmn() CGI_OPTIMIZATION
         //fp = fcn(n,x); // IMSL
         //fcn(n,x,&fp,NULL,NULL); /* NAG */
         *xp = tmp - dh;
         fm = fcn(x, n); //For frprmn() CGI_OPTIMIZATION
         //fm = fcn(n,x); // IMSL
         //fcn(n,x,&fm,NULL,NULL);

         //=== Checking the boundary condition for the minimization problem.
         if (fp >= NEARINFINITY)  *g = (f0-fm)*dhi;
         else if (fm >= NEARINFINITY)  *g = (fp-f0)*dhi;
         else  *g = (fp-fm)*dh2i;

         *xp = tmp;                        // Put the original value of x[i] back to x[i] so that the content x[i] is still unaltered.
      }

   }
   else {
      for (i=0, xp=x; i<n; i++, xp++, g++) {
         dh = fabs(*xp)<=1 ? STPS : STPS*(*xp);
         tmp = *xp;
         *xp += dh;
         dh = *xp - tmp;                   // This increases the precision slightly.
         fp = fcn(x, n);   //For frprmn() CGI_OPTIMIZATION
         //fp = fcn(n,x); // IMSL
         //fcn(n,x,&fp,NULL,NULL); /* NAG */
         *xp = tmp - dh;
         fm = fcn(x, n); //For frprmn() CGI_OPTIMIZATION
         //fm = fcn(n,x); // IMSL
         //fcn(n,x,&fm,NULL,NULL);

         //=== Checking the boundary condition for the minimization problem.
         if (fp >= NEARINFINITY)  *g = (f0-fm)/dh;
         else if (fm >= NEARINFINITY)  *g = (fp-f0)/dh;
         else  *g = (fp-fm)/(2.0*dh);

         *xp = tmp;                        // Put the original value of x[i] back to x[i] so that the content x[i] is still unaltered.
      }
   }
}
#undef STPS
//-------------------------------
//Forward difference gradient: much faster than gradcd_gen() when the objective function is very expensive to evaluate.
//-------------------------------
#define STPS 1.0e-04    // 6.0554544523933391e-6 step size = pow(DBL_EPSILON,1.0/3)
void gradfd_gen(double *g, double *x, int n, double (*fcn)(double *x, int n), double *grdh, double f0) {
   //Outputs:
   //  g: the gradient n-by-1 g (no need to be initialized).
   //Inputs:
   //  x: the vector point at which the gradient is evaluated.  No change in the end although will be added or substracted by dh during the function (but in the end the original value will be put back).
   //  n: the dimension of g or x.
   //  fcn(): the function for which the gradient is evaluated
   //  grdh: step size.  If NULL, then dh is set automatically; otherwise, grdh is taken as a step size, often set as 1.0e-004.
   //  f0: the value of (*fcn)(x).   NOT used in this function except dealing with the boundary (NEARINFINITY) for the
   //    minimization problem, but to be compatible with a genral function call where, say, gradfw_gen() and cubic
   //    interpolation of central difference method will use f0.

   double dh, dhi, fp, tmp, *xp;
   int i;
   if (grdh) {
      dhi = 1.0/(dh=*grdh);
      for (i=0, xp=x; i<n; i++, xp++, g++) {
         dh = fabs(*xp)<=1 ? STPS : STPS*(*xp);
         tmp = *xp;
         *xp += dh;
         if ( (fp=fcn(x, n)) < NEARINFINITY )  *g = (fp-f0)*dhi;   //For frprmn() CGI_OPTIMIZATION
         else {
            //Switches to the other side of the boundary.
            *xp = tmp - dh;
            *g = (f0-fcn(x,n))*dhi;
         }
         *xp = tmp;                        // Put the original value of x[i] back to x[i] so that the content x[i] is still unaltered.
      }

   }
   else {
      for (i=0, xp=x; i<n; i++, xp++, g++) {
         dh = fabs(*xp)<=1 ? STPS : STPS*(*xp);
         tmp = *xp;
         *xp += dh;
         dh = *xp - tmp;                   // This increases the precision slightly.
         if ( (fp=fcn(x, n)) < NEARINFINITY )  *g = (fp-f0)/dh;   //For frprmn() CGI_OPTIMIZATION
         else {
            //Switches to the other side of the boundary.
            *xp = tmp - dh;
            *g = (f0-fcn(x,n))/dh;
         }

         *xp = tmp;                        // Put the original value of x[i] back to x[i] so that the content x[i] is still unaltered.
      }
   }
}
#undef STPS




int next_permutation(int *first, int *last)
{
   // Given the permulation, say, [3 2 1 0], the ouput is the next permulation [0 1 2 3], and so on.
   // Note that last is simply a pointer.  Because it is not allocated to a memory, it cannot be accessed.
   //   So last is used for (1) gauging the dimension size of the array first;
   //                       (2) being accssed but with --last (which points to a valid memory place), NOT last.
   //
   // first: n-by-1 vector of integers filled with 0, 1, 2, ..., n.
   // last:  simply a pointer to the address after the last element of first.  Note that no memory is allocated.

   int *i = last, *ii, *j, tmp;
   if (first == last || first == --i)
      return 0;

   for(; ; ) {
      ii = i;
      if (*--i < *ii) {
         j = last;
         while (!(*i < *--j));
         tmp = *i; *i = *j; *j = tmp;
         for (; ii != last && ii != --last; ++ii) {
            tmp = *ii; *ii = *last; *last = tmp;
         }
         return 1;
      }
      if (i == first) {
         for (; first != last && first != --last; ++first) {
            tmp = *first; *first = *last; *last = tmp;
         }
         return 0;
      }
   }
}



/**
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void permute_matrix(double *a, int n, int *indx) {
   double *b;
   int nn=n*n;
   register int i;
   b = calloc(nn,sizeof(double));
   memcpy(b, a, nn*sizeof(double));
   for (i=0; i<nn; i++, a++)
      *a = b[indx[i%n]+indx[i/n]*n];
}

int main() {
   double a[9]={1,2,3,4,5,6,7,8,9};
   int indx[3]={1,2,0};
   permute_matrix(a,3,indx);
   return 0;
}
/**/


int fn_cumsum_int(int *x_v, const int d_x_v) {
   //Outputs:
   //  x_v: an int vector of cumulative sums over an input int vector.
   //  return: the sum of an input int vector.
   //Inputs:
   //  x_v: a vector of ints.
   //  d_x_v: dimension of x_v.
   //
   // Compute cumulative sums of a vector of ints.
   int _i;

   if (x_v==NULL) fn_DisplayError(".../cstz/fn_cumsum_lf:  x_v must be allocated with memory");

   for (_i=1; _i<d_x_v; _i++) {
      x_v[_i] = x_v[_i-1] + x_v[_i];
   }

   return (x_v[d_x_v-1]);
}


double fn_cumsum_lf(double *x_v, const int d_x_v) {
   //Outputs:
   //  x_v: a double vector of cumulative sums over an input double vector.
   //  return: the sum of an input double vector.
   //Inputs:
   //  x_v: a vector of doubles.
   //  d_x_v: dimension of x_v.
   //
   // Compute cumulative sums of a vector of doubles.
   int _i;

   if (!x_v) fn_DisplayError(".../cstz/fn_cumsum_lf:  x_v must be allocated with memory");

   for (_i=1; _i<d_x_v; _i++) {
      x_v[_i] = x_v[_i-1] + x_v[_i];
   }

   return (x_v[d_x_v-1]);
}


double fn_mean(const double *a_v, const int _n) {
   int _i;
   double x=0.0;

   for (_i=0; _i<_n; _i++)  x += a_v[_i];
   x /= (double)_n;

   return x;
}

//<<---------------
static double *tz_BaseForComp;          // This base variable is to be sorted and thus made global for this source file.
void fn_SetBaseArrayForComp(TSdvector *x_dv)
{
   if ( !x_dv->flag )   fn_DisplayError(".../cstz.c/ftd_SetBaseArrayForComp(): input vector used for comparison must be given legal values");
   else  tz_BaseForComp = x_dv->v;
}
int fn_compare(const void *i1, const void *i2)
{
   // Ascending order according to tz_BaseForComp.
   return ( (tz_BaseForComp[*((int*)i1)]<tz_BaseForComp[*((int*)i2)]) ? -1 : (tz_BaseForComp[*((int*)i1)]>tz_BaseForComp[*((int*)i2)]) ? 1 : 0 );
}
int fn_compare2(const void *i1, const void *i2)
{
   // Descending order according to tz_BaseForComp.
   return ( (tz_BaseForComp[*((int*)i1)]<tz_BaseForComp[*((int*)i2)]) ? 1 : (tz_BaseForComp[*((int*)i1)]>tz_BaseForComp[*((int*)i2)]) ? -1 : 0);
}
//--------------->>



//<<---------------
// WZ normalization on VARs.
//--------------->>
void fn_wznormalization(TSdvector *wznmlz_dv, TSdmatrix *A0draw_dm, TSdmatrix *A0peak_dm)
{
   //Outputs:
   //  wznmlz_dv (n-by-1):  If negative, the sign of the equation must switch; if positive: no action needs be taken.
   //  A0draw_dm (n-by-n):  replaced by wz-normalized draw.
   //Inputs:
   //  A0draw_dm (n-by-n):  a draw of A0.
   //  A0peak_dm (n-by-n):  reference point to which normalized A0draw_dm is closest.
   int _j, _n,
       errflag = -2;
   double *v;
   TSdmatrix *X_dm = NULL;

   if ( !A0peak_dm )  fn_DisplayError(".../cstz.c/fn_wznormalization():  input matrix for ML estimates must be created (memory allocated) and have legal values");
        //This is a minimum check to prevent crash without error messages.  More robust checks are done in BdivA_rgens().

   _n = A0peak_dm->nrows;
   X_dm = CreateMatrix_lf(_n, _n);

   if ( errflag=BdivA_rgens(X_dm, A0peak_dm, '\\', A0draw_dm) ) {
      printf(".../cstz.c/fn_wznormalization(): errors when calling BdivA_rgens() with error flag %d", errflag);
      exit(EXIT_FAILURE);
   }

   diagdv(wznmlz_dv, X_dm);
   v = wznmlz_dv->v;

   for (_j=_n-1; _j>=0; _j--)
      if (v[_j]<0)  ScalarTimesColofMatrix((TSdvector *)NULL, -1.0, A0draw_dm, _j);


   //=== Destroys memory allocated for this function only.
   X_dm = DestroyMatrix_lf(X_dm);
}




//<<---------------
// Handling under or over flows with log values.
//--------------->>
int fn_update_logofsum(int N, double ynew, double *Y_N_dp, double *y_Nmax_dp)
{
   //Recursive algorithm to update Y_N (=log(sum of x_i)) for i=1, ..., N with the new value ynew = log(x_{N+1}).
   //See TVBVAR Notes p.81a.
   if (*y_Nmax_dp>=ynew)  *Y_N_dp = log( exp(*Y_N_dp - *y_Nmax_dp) + exp(ynew - *y_Nmax_dp) ) + *y_Nmax_dp;
   else {
      *y_Nmax_dp = ynew;
      *Y_N_dp = log( exp(*Y_N_dp - ynew) + 1.0 ) + ynew;
   }
   return (N+1);
}
double fn_replace_logofsumsbt(double *yold, double _a, double ynew, double _b)
{
   //Outputs:
   //  *yold is replaced by log abs(a*xold + b*xnew).
   //  1.0 or -1.0: sign of a*xold + b*xnew.
   //
   //Given yold=log(xold) and ynew=log(xnew), it updates and returns yold = log abs(a*xold + b*xnew).
   //sbt: subtraction or subtract.
   //See TVBVAR Notes p.81a.
   double tmpd;
   //*yold = (*yold > ynew) ? (log( _a + _b*exp(ynew - *yold)) + *yold) : (log( _a*exp(*yold - ynew) + _b) + ynew);

   if (*yold > ynew) {
      if ((tmpd=_a + _b*exp(ynew - *yold) ) < 0.0) {
         // printf("WARNING! .../cstz.c/fn_replace_logofsumsbt(): Expression inside log is negative and the function returns the negative sign!\n");
         *yold += log(fabs(tmpd));
         return (-1.0);
      }
      else {
         *yold += log(tmpd);
         return (1.0);
      }
   }
   else {
      if ((tmpd=_a*exp(*yold - ynew) + _b) < 0.0 ) {
         // printf("WARNING! .../cstz.c/fn_replace_logofsumsbt(): Expression inside log is negative and the function returns the negative sign!\n");
         *yold = log(fabs(tmpd)) + ynew;
         return (-1.0);
      }
      else {
         *yold = log(tmpd) + ynew;
         return (1.0);
      }
   }
}





//<<---------------
// Evaluating the inverse of the chi-square cumulative distribution function.
//--------------->>
#if defined( IMSL_STATISTICSTOOLBOX )
double fn_chi2inv(double p, double df)
{
   //p = int_{0}^{\infty} chi2pdf(t, df) dt
   return (imsls_d_chi_squared_inverse_cdf(p, df));
}
#else
//No default routine yet.
#endif


//<<---------------
// Evaluating the inverse of the beta cumulative distribution function.
//--------------->>
#if defined( IMSL_STATISTICSTOOLBOX )
double fn_betainv(double p, double _alpha, double _beta)
{
   //p = int_{0}^{\infty} betapdf(t, _alpha, _beta) dt where betapdf(t,_alpha,_beta) \propt t^{_alpha-1}*(1-t)^(_beta-1}.
   return (imsls_d_beta_inverse_cdf(p, _alpha, _beta));
}
#else
//No default routine yet.
#endif


//<<---------------
// Computes log gamma (x) where gamma(n+1) = n! and gamma(x) = int_0^{\infty} e^{-t} t^{x-1} dt.
//--------------->>
#if defined( IMSL_STATISTICSTOOLBOX )
double fn_gammalog(double x)
{
   return (imsl_d_log_gamma(x));
}
#else
//No default routine yet.
#endif


//<<---------------
// Computes log beta(x, y) where beta(x, y) = gamma(x)*gamm(y)/gamma(x+y).
//--------------->>
#if defined( IMSL_STATISTICSTOOLBOX )
double fn_betalog(double x, double y)
{
   return (imsl_d_log_beta(x, y));
}
#else
//No default routine yet.
#endif



//<<---------------
// Computes log gamma (x) where gamma(n+1) = n! and gamma(x) = int_0^{\infty} e^{-t} t^{x-1} dt.
//--------------->>
#if defined( IMSL_STATISTICSTOOLBOX )
double gammalog(double x)
{
   return (imsl_d_log_gamma(x));
}
#else
//No default routine yet.
#endif




//<<---------------
// P2 algorithm ???????
//--------------->>
void psqr(double *q, int *m, double x, const double *p, int n)
{
   //Outputs:
   //  q: n-by-1 vector of
   //  m: n-by-1 vector of
   //  x: a random draw.
   //------
   //Inputs:
   //  p:  n-by-1 vector of cumulative cut-off probabilties for the error bands.
   static double qm, dq;
   static int i, dm, dn;

   for (i=0; q[i]<=x && i<n; i++) ;
   if (i==0) { q[0]=x; i++; }
   if (i==n) { q[n-1]=x; i--; }
   for (; i<n; i++) m[i]++;
   for (i=1; i<n-1; i++) {
      dq = p[i]*m[n-1];
      if (m[i]+1<=dq && (dm=m[i+1]-m[i])>1) {
         dn = m[i]-m[i-1];
         dq = ((dn+1)*(qm=q[i+1]-q[i])/dm+
            (dm-1)*(q[i]-q[i-1])/dn)/(dm+dn);
         if (qm<dq) dq = qm/dm;
         q[i] += dq;
         m[i]++;
      } else
      if (m[i]-1>=dq && (dm=m[i]-m[i-1])>1) {
         dn = m[i+1]-m[i];
         dq = ((dn+1)*(qm=q[i]-q[i-1])/dm+
            (dm-1)*(q[i+1]-q[i])/dn)/(dm+dn);
         if (qm<dq) dq = qm/dm;
         q[i] -= dq;
         m[i]--;
      }
   }
}
void piksrt(double *arr, int n)
{
   //Outputs:
   //  arr: replaced by new values.
   //Inputs:
   //  arr: n-by-1 vector ??????
   int i, j;
   double a;

   for (j=1; j<n; j++) {
      a = arr[j];
      for (i=j-1; i>=0 && arr[i]>a; i--)
         arr[i+1] = arr[i];
      arr[i+1]=a;
   }
}





//---------------------------- New functions with my own matrix notations. ---------------------
double MinVector_lf(TSdvector *x_dv) {
   //Input: no change for x_dv in this function.
   int _i, n;
   double minvalue;
   double *v;

   if (!x_dv || !x_dv->flag) fn_DisplayError(".../cstz.c/MinVector_lf():  Input vector x_dv must be (1) allocated memory and (2) assigned legal values");
   n = x_dv->n;
   v = x_dv->v;

   minvalue = v[0];
   for (_i=n-1; _i>0; _i--)
      if (v[_i]<minvalue)  minvalue = v[_i];

   return( minvalue );
}




//---------------------------- Not used often ---------------------
void fn_cumsum(double **aos_v, int *aods_v, double *v, int d_v) {
   // Compute a cumulative sum of a vector.
   //
   // v: an n-by-1 vector.
   // d_v: n -- size of the vector v to be used for a cumulative sum.
   // aos_v: address of the pointer to the n-by-1 vector s_v.
   // aods_v: address of the size of the dimension of s_v.
   //----------
   // *aos_v: An n-by-1 vector of cumulative sum s_v.
   // *aods_v: n -- size of the dimension for s_v.

   int ki;

   *aos_v = tzMalloc(d_v, double);
   (*aods_v) = d_v;                               // n for the n-by-1 vector s_v.
   *(*aos_v) = *v;
   if (d_v>1) {
      for (ki=1; ki<d_v; ki++) (*aos_v)[ki] = (*aos_v)[ki-1] + v[ki];
   }
}



/**
void fn_ergodp(double **aop, int *aod, mxArray *cp) {
   // Compute the ergodic probabilities.  See Hamilton p.681.
   //
   // cp: n-by-n Markovian transition matrix.
   // aop: address of the pointer to the n-by-1 vector p.
   // aod: address of the size of the dimension of p.
   //----------
   // *aop: n-by-1 vector of ergodic probabilities p.  @@Must be freed outside this function.@@
   // *aod: n -- size of the dimension for p (automatically supplied within this function).

   mxArray *gpim=NULL, *gpid=NULL;   // m: n-by-n eigvector matrix; d: n-by-n eigvalue diagonal.
   double *gpim_p, *gpid_p;             // _p:  a pointer to the corresponding mxArray whose name occurs before _p.
         //------- Note the following two lines will cause Matlab or C to crash because gpim has not been initialized so it points to garbage.
         //   double *gpim_p = mxGetPr(gpim);
         //   double *gpid_p = mxGetPr(gpid);
   int eigmaxindx,                      // Index of the column corresponding to the max eigenvalue.
       n, ki;
   double gpisum=0.0,
          eigmax, tmpd0;

   n=mxGetM(cp);                        // Get n for the n-by-n mxArray cp.
   (*aod)=n;

   *aop = tzMalloc(n, double);

   gpim = mlfEig(&gpid,cp,NULL,NULL);
   gpim_p = mxGetPr(gpim);
   gpid_p = mxGetPr(gpid);

   eigmax = *gpid_p;
   eigmaxindx = 0;
   if (n>1) {
      for (ki=1;ki<n;ki++) {
         if (gpid_p[n*ki+ki] > eigmax) {
            eigmax=gpid_p[n*ki+ki];
            eigmaxindx=ki;
         }                           // Note that n*ki+ki refers to a diagonal location in the n-by-n matrix.
      }
   }
   for (ki=0;ki<n;ki++) {
      gpisum += gpim_p[n*eigmaxindx+ki];                 // Sum over the eigmaxindx_th column.
   }
   tmpd0 = 1.0/gpisum;
   for (ki=0;ki<n;ki++) {
      (*aop)[ki] = gpim_p[n*eigmaxindx+ki]*tmpd0;                 // Normalized eigmaxindx_th column as ergodic probabilities.
   }

   mxDestroyArray(gpim);                // ????? free(gpim_p)
   mxDestroyArray(gpid);
}
/**/
