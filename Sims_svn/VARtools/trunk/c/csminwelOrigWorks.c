/*=========================================================
 * csminwel.c
 *
 * Unconstrained minimization.  Uses a quasi-Newton method with BFGS update of
 * the estimated inverse hessian.  It is robust against certain pathologies
 * common on likelihood functions.  It attempts to be robust against "cliffs",
 * i.e. hyperplane discontinuities, though it is not really clear whether what
 * it does in such cases succeeds reliably.
 *
 * function [fhat,xhat,ghat,Hhat,itct,fcount,retcodehat] = csminwelmex(fcn,x0,H0,grad,crit,nit,varargin)
 * fcn:   string naming the objective function to be minimized
 * x0:    initial value of the parameter vector
 * H0:    initial value for the inverse Hessian.  Must be positive definite.
 * grad:  Either a string naming a function that calculates the gradient, or the null matrix.
 *        If it's null, the program calculates a numerical gradient.  In this case fcn must
 *        be written so that it can take a matrix argument and produce a row vector of values.
 * crit:  Convergence criterion.  Iteration will cease when it proves impossible to improve the
 *        function value by more than crit.
 * nit:   Maximum number of iterations.
 * varargin: A list of optional length of additional parameters that get handed off to fcn each
 *        time it is called.
 *        Note that if the program ends abnormally, it is possible to retrieve the current x,
 *        f, and H from the files g1.mat and H.mat that are written at each iteration and at each
 *        hessian update, respectively.  (When the routine hits certain kinds of difficulty, it
 *        write g2.mat and g3.mat as well.  If all were written at about the same time, any of them
 *        may be a decent starting point.  One can also start from the one with best function value.)
 *
 * Note: to set the level of display output, change preprocessor definitions in csminwel.h
 *       to display all output, uncomment both VERBOSE_WARNINGS and VERBOSE_DETOUTPUT
 *       to display only warnings without output, uncomment VERBOSE_WARNINGS
 *       to display no ouput, comment both VERBOSE_DETOUTPUT and VERBOSE_WARNINGS
 *
 * MATLAB algorithm by Christopher Sims
 * C implementation by Iskander Karibzhanov
 * Modified by Dan Waggoner and Tao Zha
 *
 * Copyright(c) 1996 Christopher Sims
 * Copyright(c) 2003 Karibzhanov, Waggoner, and Zha
 *=======================================================
 * Revision history:
 *
 *  10/3/2002  -  1. corrected problem with memory corruption in C-MEX-file (csminwelmex.c)
 *                   (needed to switch fcnRhs[0] back to x[0] before destroying it.
 *                   If we don't do this, we will later clear previously destroyed array
 *                   (one of x[1], x[2] or x[3]) which causes memory fault.
 *                   The reason why fcnRhs[0] pointed to array other than x[0] is
 *                   because we use mxSetPr in feval and gfeval.
 *                   This was not a problem in C-file (csminwel.c).
 *
 * 10/11/2002  -  1. changed csminit function to avoid using fPeak without first initializing it
 *                2. added two switches in csminit function to assign retcode to 7 for lambda>=4
 *                3. added one more verbose level to display only warnings or all output *
 *=======================================================*/

#include "csminwel.h"
#define STRLEN 192
static int numgrad(double *g, double *x, int n,
                   double (*fcn)(double *x, int n, double **args, int *dims),
                   double **args, int *dims);
static void csminit(double *fhat, double *xhat, int *fcount, int *retcode,
                    double *x0, double f0, double *g, int badg, double *H0, int n,
                    double (*fcn)(double *x, int n, double **args, int *dims),
                    double **args, int *dims);
static void bfgsi(double *H, double *dg, double *dx, int n, int nn);
static int peakwall(double *g, int retcode, double *x, int n,
                    int (*gfcn)(double *x, int n, double *g, double **args, int *dims),
                    double (*fcn)(double *x, int n, double **args, int *dims),
                    double **args, int *dims);
static double times(double *x, double *y, int n);
static double *mtimes(double *x, double *y, int n, int nn);
static double *mminus(double *x, double *y, int n);


static FILE *fptr_interesults = (FILE *)NULL;   //Printing intermediate results to a file.
static char filename_sp2vecs[STRLEN];  //Two vectors.  1st row: numerical gradient; 2nd row: vectorized parameters.


void csminwel(double (*fcn)(double *x, int n, double **args, int *dims),
              double *xh, int n, double *H, double *gh,
              int (*gfcn)(double *x, int n, double *g, double **args,
              int *dims), double *fh, double crit, int *itct, int nit,
              int *fcount, int *retcodeh, double **args, int *dims) {

   int done=0, badg[4], badgh, nogh=1, stuck=0;
   double *x[4], *g[4], f[4], *dg, *dx;
   int retcode[3], fc=0, ih, nn, i;
   #ifdef VERBOSE_DETOUTPUT
   time_t begtime, currentime;
   #endif


   nn = n*n;     /* n: dimension size of x or xh */
   *itct = -1;    /* itct: number of actual iterations */
   *fcount = -1;  /* fcount: number of evaluations of the function */

   for (i=0; i<4; i++)
      x[i] = tzMalloc(n, double);  //x[i] = calloc(n, sizeof(double)); Commented out by TZ.
   memcpy(x[0],xh,n*sizeof(double));

   for (i=0; i<4; i++)
      g[i] = tzMalloc(n, double);    //calloc(n, sizeof(double));  Commented out by TZ.

   f[0] = fcn(x[0],n,args,dims);

   if (f[0] > 1.0e+50) {
      printf("Bad initial parameter.\n");
      return;
   }

   if (gfcn)
      /* if grad is a string, compute it */
      badg[0] = gfcn(x[0],n,g[0],args,dims);
   else
      /* if grad is not string, compute it */
      badg[0] = numgrad(g[0],x[0],n,fcn,args,dims);
   retcode[2] = 101;
   /* iterate until done is false */
   while (!done) {
      #ifdef VERBOSE_DETOUTPUT
      time(&begtime);
      #endif

      for (i=0; i<n; i++)
         g[1][i] = g[2][i] = g[3][i] = 0;

//         #ifdef VERBOSE_DETOUTPUT
//         printf("-----------------\n-----------------\n");
//         printf("f at the beginning of new iteration, %.10f\nx = ",f[0]);
//         for (i=0; i<n; i++) {
//            printf("%15.8g ",x[0][i]);
//            if (i%4==3) printf("\n");
//         }
//         if (i%4>0) printf("\n");
//         #endif

         (*itct)++;
         csminit(&f[1],x[1],&fc,&retcode[0],x[0],f[0],g[0],badg[0],H,n,fcn,args,dims);
         *fcount += fc;
         /* if retcode1=1 gradient is zero and you are at the peak */
         if (retcode[0]!=1) {
            badg[1] = peakwall(g[1],retcode[0],x[1],n,gfcn,fcn,args,dims);
            /* Bad gradient or back and forth on step length.
               Possibly at cliff edge. Try perturbing search direction. */
            if (badg[1]) {
               double *Hcliff = tzMalloc(nn, double);   //calloc(nn,sizeof(double));  Commented out by TZ.
               double randmax=1/RAND_MAX;
               /* if stuck, give it another try by perturbing Hessian */
               memcpy(Hcliff,H,nn*sizeof(double));
               for (i=0; i<nn; i+=n+1)
                  Hcliff[i] *= 1+rand()*randmax;

               #ifdef VERBOSE_WARNINGS
               printf("Cliff.  Perturbing search direction.\n");
               #endif

               csminit(&f[2],x[2],&fc,&retcode[1],x[0],f[0],g[0],badg[0],Hcliff,n,fcn,args,dims);
               *fcount += fc;
               if (f[2] < f[0]) {
                  badg[2] = peakwall(g[2],retcode[1],x[2],n,gfcn,fcn,args,dims);
                  if (badg[2]) {
                     double *xx = tzMalloc(n, double), nx;   //calloc(n,sizeof(double)), nx;  Commented out by TZ.

                  #ifdef VERBOSE_WARNINGS
                  printf("Cliff again.  Try traversing.\n");
                  #endif

                  for (i=0; i<n; i++)
                     xx[i] = x[2][i]-x[1][i];
                  nx = times(xx,xx,n);
                  if (sqrt(nx) < 1e-13) {
                     f[3] = f[0];
                     memcpy(x[3],x[0],n*sizeof(double));
                     badg[3] = 1;
                     retcode[2] = 101;
                  } else {
                     double *gcliff = tzMalloc(n, double),  //calloc(n,sizeof(double)),  Commented out by TZ.
                            *eye = tzMalloc(n, double);  //calloc(n,sizeof(double));  Commented out by TZ.
                     double dfnx = (f[2]-f[1])/nx;
                     for (i=0; i<n; i++) {
                        gcliff[i] = dfnx*xx[i];
                        eye[i*(n+1)] = 1;
                     }
                     csminit(&f[3],x[3],&fc,&retcode[2],x[0],f[0],gcliff,0,eye,n,fcn,args,dims);
                     *fcount += fc;
                     badg[3] = peakwall(g[3],retcode[2],x[3],n,gfcn,fcn,args,dims);
                     free(eye);
                     free(gcliff);
                  }
                  free(xx);
               } else {
                  f[3] = f[0];
                  memcpy(x[3],x[0],n*sizeof(double));
                  badg[3] = 1;
                  retcode[2] = 101;
               }
            } else {
               f[3] = f[0];
               memcpy(x[3],x[0],n*sizeof(double));
               badg[3] = 1;
               retcode[2] = 101;
            }
            free(Hcliff);
         } else {
            /* normal iteration, no walls, or else we're finished here. */
            f[2] = f[0];
            f[3] = f[0];
            badg[2] = 1;
            badg[3] = 1;
            retcode[1] = 101;
            retcode[2] = 101;
         }
      }


      /* how to pick gh and xh */
      if (f[3]<f[0] && badg[3]==0) {
         /* if 3 (transversing) was needed, it improved and gradient is good, take that */
         ih = 3;
         *fh = f[3];
         memcpy(xh,x[3],n*sizeof(double));
         memcpy(gh,g[3],n*sizeof(double));
         badgh = badg[3];
         *retcodeh = retcode[2];
      }
      else if (f[2]<f[0] && badg[2]==0) {
         /* if 2 (perturbig) was needed, it improved and gradient is good, take that */
         ih = 2;
         *fh = f[2];
         memcpy(xh,x[2],n*sizeof(double));
         memcpy(gh,g[2],n*sizeof(double));
         badgh = badg[2];
         *retcodeh = retcode[1];
      }
      else if (f[1]<f[0] && badg[1]==0) {
         /* if first try went fine, take that */
         ih = 1;
         *fh = f[1];
         memcpy(xh,x[1],n*sizeof(double));
         memcpy(gh,g[1],n*sizeof(double));
         badgh = badg[1];
         *retcodeh = retcode[0];
      }
      else {
         /* if nothing worked, just take the min of your attempts and compute the gradient */
         if (f[1] <= f[2])
            if (f[1] <= f[3]) ih = 1;
            else ih = 3;
         else if (f[2] <= f[3]) ih = 2;
            else ih = 3;
         *fh = f[i];
         memcpy(xh,x[ih],n*sizeof(double));
         *retcodeh = retcode[ih-1];
         if (nogh) {
            nogh = 0;
            if (gfcn)
               badgh = gfcn(xh,n,gh,args,dims);
            else
               badgh = numgrad(gh,xh,n,fcn,args,dims);
         }
         badgh = 1;
      }
      /* end of picking */
      stuck = fabs(*fh-f[0]) < crit;
      /* if nothing REALLY worked, too bad, you're stuck */
      if (!badg[0] && !badgh && !stuck) {
         /* if you are not stuck, update H0 matrix */
         dg = mminus(gh,g[0],n);
         dx = mminus(xh,x[0],n);
         bfgsi(H,dg,dx,n,nn);
         free(dx);
         free(dg);
      }

      #ifdef VERBOSE_DETOUTPUT
      //=== Prints out intermediate results.
      printf("========================================\n");
      printf(" (1) New value of the obj. func. on iteration %d: %.9f\n (2) Old value: %.9f\n (3) Downhill improvement: %.9f\n",
             (int)*itct, *fh, f[0], f[0]-(*fh));

      time(&currentime);
      //=== Times the iterative progress.
      printf(" (4) Seconds to complete one iteration: %0.4f\n (5) Current time of day: %s\n\n", difftime(currentime, begtime), ctime(&currentime));
      fflush(stdout);                // Flush the buffer to get out this message without delay.
      #endif

      //--------- Prints outputs to a file. ---------
      if ( !(fptr_interesults = fopen(filename_sp2vecs,"w")) ) {
         printf("\n\nUnable to create the starting point data file %s in csminwel.c!\n", filename_sp2vecs);
         getchar();
         exit(EXIT_FAILURE);
      }
      fprintf(fptr_interesults, "--------Numerical gradient---------\n");
      for (i=0; i<n; i++)  fprintf(fptr_interesults, " %0.16e ", gh[i]);
      fprintf(fptr_interesults, "\n");
      fprintf(fptr_interesults, "--------Restarting point---------\n");
      for (i=0; i<n; i++)  fprintf(fptr_interesults, " %0.16e ", xh[i]);
      fprintf(fptr_interesults, "\n\n");
      tzFclose(fptr_interesults);

      if ((int)*itct > nit) {
         #ifdef VERBOSE_WARNINGS
         printf("\nWarning: termination as the maximum number of iterations is reached.\n");
         #endif
         done = 1;
      }
      else if (stuck) {
         #ifdef VERBOSE_DETOUTPUT
         printf("improvement < crit termination\n");
         #endif

         done = 1;
      }

      #ifdef VERBOSE_WARNINGS
      switch ((int)*retcodeh) {
         case 1:
            printf("\nWarning: Zero gradient.\n"); break;
         case 2:
            printf("\nWarning: Back adjustment of stepsize didn't finish.\n"); break;
         case 3:
            printf("\nWarning: Smallest stepsize still improving too slow.\n"); break;
         case 4:
            printf("\nWarning: Forth adjustment of stepsize didn't finish.\n"); break;
         case 6:
            printf("\nWarning: Smallest step still improving too slow, reversed gradient.\n"); break;
         case 5:
            printf("\nWarning: Largest stepsize still improving too fast.\n"); break;
         case 7:
            printf("\nWarning: Possible inaccuracy in Hessian matrix.\n"); break;
      }
      #endif

      f[0] = *fh;
      memcpy(x[0],xh,n*sizeof(double));
      memcpy(g[0],gh,n*sizeof(double));
      badg[0] = badgh;
   }
   for (i=0; i<4; i++) {
      free(g[i]);
      free(x[i]);
   }
}

/**/
#define SCALE 1.0
static int numgrad(double *g, double *x, int n,
                   double (*fcn)(double *x, int n, double **args, int *dims),
                   double **args, int *dims) {
   //Forward difference gradient method.
   double delta, deltai;
   double f0, g0, ff, tmp, *xp;
   int i;
   int badg;
   f0 = fcn(x,n,args,dims);
   badg = 0;
   for (i=0, xp=x; i<n; i++, xp++, g++) {
      delta=SCALE*1.0e-5, deltai=1.0e+5/SCALE;

      tmp = *xp;
      *xp += delta;
      delta = *xp - tmp;                   // This increases the precision slightly.  Added by TZ.
      if ( (ff=fcn(x,n,args,dims)) < NEARINFINITY )   g0 = (ff-f0)*deltai;   //Not over the boundary.
      else {
         //Switches to the other side of the boundary.
         *xp = tmp - delta;
         g0 = (f0-fcn(x,n,args,dims))*deltai;
      }

      *xp = tmp;       //Puts back to the original place.  TZ, 9/03.
      if (fabs(g0)<1.0e+15)
         *g = g0;
      else {
         #ifdef VERBOSE_WARNINGS
         printf("Bad gradient.\n");
         #endif

         *g = 0;
         badg = 1;
      }
   }
   return badg;
}
#undef SCALE

/**
#define STPS 1.0e-04    // 6.0554544523933391e-6 step size = pow(DBL_EPSILON,1.0/3)
static int numgrad(double *g, double *x, int n,
                   double (*fcn)(double *x, int n, double **args, int *dims),
                   double **args, int *dims) {
   //Central difference gradient method.  Added by TZ.
   double dh;
   double f0, fp, fm, tmp, *xp;
   int i;
   int badg;

   badg = 0;
   for (i=0, xp=x; i<n; i++, xp++, g++) {
      dh = fabs(*xp)<=1 ? STPS : STPS*(*xp);

      tmp = *xp;
      *xp += dh;
      dh = *xp - tmp;                   // This increases the precision slightly.
      fp = fcn(x,n,args,dims);
      *xp = tmp - dh;
      fm = fcn(x,n,args,dims);

      //=== Checking the boundary condition for the minimization problem.
      if (fp >= NEARINFINITY) {
         *xp = tmp;       //Puts back to the original place.  TZ, 9/03.
         f0 = fcn(x,n,args,dims);
         *g = (f0-fm)/dh;
      }
      else if (fm >= NEARINFINITY) {
         //Switches to the other side of the boundary.
         *xp = tmp;       //Puts back to the original place.  TZ, 9/03.
         f0 = fcn(x,n,args,dims);
         *g = (fp-f0)/dh;
      }
      else {
         *g = (fp-fm)/(2.0*dh);
         *xp = tmp;       //Puts back to the original place.  TZ, 9/03.
      }

      if (fabs(*g)>1.0e+15) {
         #ifdef VERBOSE_WARNINGS
         printf("Bad gradient.\n");
         #endif
         *g = 0.0;
         badg = 1;
      }
   }
   return badg;
}
#undef STPS
/**/


#define ANGLE 0.05  //When output of this variable becomes negative, we have a wrong analytical graident.
                    //.005 works for identified VARs and OLS.
#define THETA .4    //(0<THETA<.5) THETA near .5 makes long line searches, possibly fewer iterations.
                    //.1 works for OLS or other nonlinear functions.
                    //.3 works for identified VARs.
#define FCHANGE 1000
#define MINLAMB 1e-9
#define MINDFAC .01
static void csminit(double *fhat, double *xhat, int *fcount, int *retcode,
                    double *x0, double f0, double *g, int badg, double *H0, int n,
                    double (*fcn)(double *x, int n, double **args, int *dims),
                    double **args, int *dims) {
   double lambda=1, gnorm=0, dxnorm=0, factor=3, lambdaPeak=0;
   double f, dfhat, a, tmp, fPeak=f0, lambdaMax=DBL_MAX;
   double *dx, *dxtest;
   int done=0, shrink=1, shrinkSignal, growSignal;
   int i;

   *fhat = f0;
   *fcount = 0;
   *retcode = 0;
   gnorm = sqrt(times(g,g,n));
   if ((gnorm < 1.e-12) && !badg)
      *retcode = 1;  /* gradient convergence */
   else {
      /* with badg 1, we don't try to match rate of improvement to directional
         derivative.  We're satisfied just to get some improvement in f. */
      dx = tzMalloc(n, double);   //dx = calloc(n, sizeof(double));  Commented out by TZ.
      //if (!dx) printf("Dynamic memory allocation error.\n");  Commnted out by TZ.
      for (i=0; i<n; i++)
         dx[i] = -times(&H0[i*n],g,n);
      dxnorm = sqrt(times(dx,dx,n));
      if (dxnorm > 1e12) {
         #ifdef VERBOSE_WARNINGS
         printf("Near-singular H problem.\n");
         #endif

         for (i=0; i<n; i++)
            dx[i] *= FCHANGE/dxnorm;
      }
      dfhat = times(dx,g,n);
      if (!badg) {
         /* test for alignment of dx with gradient and fix if necessary */
         a = -dfhat/(gnorm*dxnorm);
         if (a<ANGLE) {
            tmp = (ANGLE*dxnorm+dfhat/gnorm)/gnorm;
            for (i=0; i<n; i++)
               dx[i] -= tmp*g[i];
            dfhat = times(dx,g,n);
            dxnorm = sqrt(times(dx,dx,n));

            #ifdef VERBOSE_DETOUTPUT
            printf("Correct for low angle: %g\n",a);
            #endif
         }
      }

      #ifdef VERBOSE_DETOUTPUT
      printf("Predicted improvement: %18.9f, Norm of gradient: %18.9f\n",-dfhat/2,gnorm);
      #endif

      dxtest = tzMalloc(n, double);  //calloc(n, sizeof(double));  Commented out by TZ.
      while (!done) {
         for (i=0; i<n; i++)
            dxtest[i] = x0[i]+dx[i]*lambda;
         f = fcn(dxtest,n,args,dims);

         #ifdef VERBOSE_DETOUTPUT
         printf("lambda = %10.5g; f = %20.7f\n",lambda,f);
         #endif

         if (f<*fhat) {
            *fhat = f;
            memcpy(xhat,dxtest,n*sizeof(double));
         }
         (*fcount)++;
         tmp = -THETA*dfhat*lambda;

         /* the optimal lambda should be such that f0-f > -THETA*dfhat*lambda (see Berndt et al.)
            If that's not the case and grad is good, OR
            if grad is bad and f is not going down, shrinkSignal = 1 */
         shrinkSignal = ( !badg && (f0-f <= (tmp>0?tmp:0)) ) ||
                         ( badg && (f0-f < 0 ) );

         /* the optimal lambda should also be such that f0-f<-(1-THETA)*dfhat*lambda
            If that's not the case with lambda>0, AND grad is good, growthSignal = 1 */
         growSignal = !badg && ( (lambda > 0)  &&  (f0-f >= -(1-THETA)*dfhat*lambda) );

         /* If shrinkSignal=1 AND ( lambda>lambdaPeak or lambda negative )
            (note when lambdaPeak=0 the second part only excludes lambda=0)
            try shrinking lambda */
         if ( shrinkSignal && ( (lambda>lambdaPeak) || (lambda<0) ) ) {
            /* if shrink=0 OR lambda/factor is already smaller than lambdaPeak, increase factor */
            if ( (lambda>0) && ((!shrink) || (lambda/factor <= lambdaPeak)) ) {
               shrink = 1;
               factor = pow(factor,.6);
               while (lambda/factor <= lambdaPeak)
                  factor = pow(factor,.6);
               if (fabs(factor-1)<MINDFAC) {
                  if (fabs(lambda) < 4)
                     *retcode = 2;
                  else
                     *retcode = 7;
                  done = 1;
               }
            }
            if ((lambda<lambdaMax) && (lambda>lambdaPeak))
               lambdaMax=lambda;
            /* shrink lambda */
            lambda /= factor;
            /* if lambda has already been shrunk as much as possible */
            if (fabs(lambda) < MINLAMB)
               /* if lambda is positive AND you have not made any improvement
                  try going against gradient, which may be inaccurate */
               if ((lambda > 0) && (f0 <= *fhat))
                  lambda = -lambda*pow(factor,6);
               else {
                  /* if lambda is negative: let it be known and quit trying */
                  if (lambda < 0)
                     *retcode = 6;
                  /* if you have not made any imporvement:
                     let it be known and quit trying */
                  else
                     *retcode = 3;
                  done = 1;
               }
         }
         /* If growSignal=1 and lambda positive OR ( lambda>lambdaPeak or lambda negative )
            (note when lambdaPeak=0 the second part only excludes lambda=0)
            try increase lambda */
         else
            if ( (growSignal && (lambda > 0) ) ||
               ( shrinkSignal && (lambda <= lambdaPeak) && (lambda > 0) ) ) {
               if (shrink) {
                  shrink = 0;
                  factor = pow(factor,.6);
                  if (fabs(factor-1) < MINDFAC) {
                     if (fabs(lambda) < 4)
                        *retcode = 4;
                     else
                        *retcode = 7;
                     done = 1;
                  }
               }
               if ( (f<fPeak) && (lambda>0) ) {
                  fPeak = f;
                  lambdaPeak = lambda;
                  if (lambdaMax <= lambdaPeak)
                     lambdaMax = lambdaPeak*factor*factor;
               }
               /* increase lambda (up to 1e20) */
               lambda *= factor;
               /* if lambda has been increased up to the limit and
                  you have not made any imporvement:
                  let it be known and quit trying */
               if (fabs(lambda) > 1e20) {
                  *retcode = 5;
                  done = 1;
               }
            }
            /* If growthSignal=shrinkSignal=0 you found a good lambda, you are done */
            else {
               done = 1;
               *retcode = factor<1.2 ? 7 : 0;
            }
      }
      free(dxtest);
      free(dx);
   }
   #ifdef VERBOSE_DETOUTPUT
   printf("Norm of dx %10.5g\n", dxnorm);
   #endif
}
#undef ANGLE
#undef THETA
#undef FCHANGE
#undef MINLAMB
#undef MINDFAC


static double times(double *x, double *y, int n) {
   double z = 0;
   int i;
   for (i=0; i<n; i++, x++, y++)
      z += (*x)*(*y);
   return z;
}

static int peakwall(double *g, int retcode, double *x, int n,
                    int (*gfcn)(double *x, int n, double *g, double **args, int *dims),
                    double (*fcn)(double *x, int n, double **args, int *dims),
                    double **args, int *dims) {
   /* if retcode=2 or 4 you have shrunk or increased lambda as much as you could:
      exhausted search possibilities the csminit step has failed */
   if (retcode==2 || retcode==4)
      return 1;
   else
      /* if you are not at the peak but the csminit has improved,
         compute the gradient again to update H0 */
      if (gfcn)
         return gfcn(x,n,g,args,dims);
      else
         return numgrad(g,x,n,fcn,args,dims);
}

static void bfgsi(double *H, double *dg, double *dx, int n, int nn) {
   double *Hdg, *dxdx, *dxHdg, *Hdgdx;
   double dgdx, m;
   int i;

   Hdg = tzMalloc(n, double);  //calloc(n, sizeof(double));  Commented out by TZ.
   //if (!Hdg) printf("Dynamic memory allocation error.\n");  Commented out by TZ.

   /* Hdg = H0*dg; */
   for (i=0; i<n; i++)
      Hdg[i] = times(&H[i*n],dg,n);
   /* dgdx = dg'*dx; */
   dgdx = 1/times(dg,dx,n);
   if (fabs(dgdx)<1e12) {
      dxdx = mtimes(dx,dx,n,nn);
      dxHdg = mtimes(dx,Hdg,n,nn);
      Hdgdx = mtimes(Hdg,dx,n,nn);
      m = 1+times(dg,Hdg,n)*dgdx;
      for (i=0; i<nn; i++, H++, dxdx++, dxHdg++, Hdgdx++)
         *H += (m*(*dxdx)-(*dxHdg)-(*Hdgdx))*dgdx;
      free(Hdgdx-nn);
      free(dxHdg-nn);
      free(dxdx-nn);
   }
   else {
   #ifdef VERBOSE_WARNINGS
   printf("BFGS update failed.\n");
   printf("|dg| = %f  |dx| = %f\n",sqrt(times(dg,dg,n)),sqrt(times(dx,dx,n)));
   printf("dg\'*dx = %f\n",dgdx);
   printf("|H*dg| = %f\n",sqrt(times(Hdg,Hdg,n)));
   #endif
   }
   free(Hdg);
}

static double *mtimes(double *x, double *y, int n, int nn) {
   double *x0;
   double *z;
   int i, j;
   z = tzMalloc(nn, double);  //calloc(nn, sizeof(double));  Commented out by TZ.
   for (i=0, x0=x; i<n; i++, y++)
      for (j=0, x=x0; j<n; j++, x++, z++)
         *z = (*x)*(*y);
   return z-nn;
}

static double *mminus(double *x, double *y, int n) {
   double *z;
   int i;
   z = tzMalloc(n, double);  //calloc(n, sizeof(double));  Commented out by TZ.
   for (i=0; i<n; i++, x++, y++, z++)
      *z = (*x)-(*y);
   return z-n;
}


//=== Extern function to be accessed by other C files.
void csminwel_SetPrintFile(char *filename) {
   if (!filename)   sprintf(filename_sp2vecs, "outdata5csminwel.prn");  //Default filename.
   else if (STRLEN-1 < strlen(filename))  fn_DisplayError(".../csminwel.c:  the allocated length STRLEN for filename_sp2vecs is too short.  Must increase the string length");
   else  strcpy(filename_sp2vecs, filename);
}

#undef STRLEN
