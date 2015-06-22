/*************************************************************
 *  Conjugate Gradient Minimization Methods.  See Numerical Recipes in C by Press, Flannery, Teukolsky, and Vetterling.
 *  (I)  frprmn():  Plolak-Ribiere method with the line minimization without using the derivative information.
 *  (II) dlinmin():  Fletcher-Reeves method with the line minimization using the derivative information.
 *
 * Modified by Tao Zha, 27 August 2003.
*************************************************************/

#include "congradmin.h"
//#include "nrutil.h"


static void linmin(double p[], double xi[], int n, double *fret, double (*func)(double [], int));
static double f1dim(double x);
static void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, double (*func)(double));
static double brent(double ax, double bx, double cx, double (*f)(double), double tol, double *xmin);

#define ITMAX 200          //Maximum number of iterations.
#define EPS 1.0e-10        //Small number to rectify special case of converging to exactly zero function value.
#define FREEALL {tzDestroy(xi); tzDestroy(h); tzDestroy(g);}
void frprmn(double p[], int n, int *iter, double *fret,
            double (*func)(double [], int), void (*dfunc)(double [], double [], int, double (*func)(double [], int), double, double),
            double ftol, double grdh) {
   //Outputs:
   //  p[0, ..., n-1]:  the location of the minimum if it converges, which replaces the starting value.
   //  iter:  the number of iterations that were performed.
   //  fret:  the minimum value of the function.
   //Inputs:
   //  p[0, ..., n-1]:  a starting point for the minimization.
   //  n:  the dimension of p.
   //  ftol:  the convergence tolerance on the objective function value.
   //  grdh:  the user's specified step size for a numerical gradient.  If 0.0, dfunc() (i.e., gradcd_gen()) will select it automatically.
   //  func():  the objective function.
   //  dfunc(): the gradient function computing the numerical gradient.  In the form of gradcd_gen() in cstz.c.
   int j, its;
   double gg, gam, fp, dgg;
   double *g=NULL, *h=NULL, *xi=NULL;

   g=tzMalloc(n, double);
   h=tzMalloc(n, double);
   xi=tzMalloc(n, double);
   fp=(*func)(p, n);
   (*dfunc)(xi, p, n, func, grdh, fp);
   for (j=n-1;j>=0;j--) {
		g[j] = -xi[j];
		xi[j]=h[j]=g[j];
	}
   for (its=0;its<ITMAX;its++) {
		*iter=its;
      #if defined (CGI_OPTIMIZATION)
         linmin(p,xi,n,fret,func);       //Next statement is the normal return.
      #elif defined (CGII_OPTIMIZATION)
         fn_DisplayError("I have not got time to put dlinmin() in frprmn() for CGII_OPTIMIZATION");
      #else
         fn_DisplayError("The minimization routine frprmn() requires activating CGI_OPTIMIZATION or CGII_OPTIMIZATION in tzmatlab.h")
      #endif
		if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)) {
			FREEALL
			return;
		}
      fp=(*func)(p, n);
      (*dfunc)(xi, p, n, func, grdh, fp);
		dgg=gg=0.0;
      for (j=n-1;j>=0;j--) {
			gg += g[j]*g[j];
			dgg += (xi[j]+g[j])*xi[j];
		}
		if (gg == 0.0) {
			FREEALL
			return;
		}
		gam=dgg/gg;
      for (j=n-1;j>=0;j--) {
			g[j] = -xi[j];
			xi[j]=h[j]=g[j]+gam*h[j];
		}
	}
   fn_DisplayError("The maximum number of iterations is reached before convergence in frprmn()");
}
#undef ITMAX
#undef EPS
#undef FREEALL



static int ncom;
static double *pcom=NULL, *xicom=NULL, (*nrfunc)(double [], int);   //nrfunc(), pcom, ncom, and xicom will be used by f1dim().
#define TOL 2.0e-4
static void linmin(double p[], double xi[], int n, double *fret, double (*func)(double [], int)) {
   //Outputs:
   //  p[0, ..., n-1]:  a returned and reset value.
   //  xi[0, ..., n-1]:  a value repaced by the actual vector displacement that p was moved.
   //  fret:  the value of func at the returned location p.
   //Inputs:
   //  p[0, ..., n-1]:  a given point.
   //  xi[0, ..., n-1]:  a given multidimensional direction.
   //  n:  the dimension of p and xi.
   //  func():  the objective function.
	int j;
   double xx,xmin,fx,fb,fa,bx,ax;

	ncom=n;
   pcom = tzMalloc(n, double);
   xicom = tzMalloc(n, double);
	nrfunc=func;
   for (j=n-1;j>=0;j--) {
		pcom[j]=p[j];
		xicom[j]=xi[j];
	}
	ax=0.0;
	xx=1.0;
	mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
	*fret=brent(ax,xx,bx,f1dim,TOL,&xmin);
   for (j=n-1;j>=0;j--) {
		xi[j] *= xmin;
		p[j] += xi[j];
	}
   tzDestroy(xicom);
   tzDestroy(pcom);
}
#undef TOL


//extern int ncom;
//extern float *pcom,*xicom,(*nrfunc)(float []);

//===========  Must accompany limin().
static double f1dim(double x) {
	int j;
   double f,*xt=NULL;

   xt = tzMalloc(ncom, double);
   for (j=ncom-1;j>=0;j--) xt[j]=pcom[j]+x*xicom[j];
   f=(*nrfunc)(xt, ncom);
   tzDestroy(xt);
	return f;
}


#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d)  {(a)=(b);(b)=(c);(c)=(d);}
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, double (*func)(double)) {
   double ulim,u,r,q,fu,dum, tmpd;

	*fa=(*func)(*ax);
	*fb=(*func)(*bx);
	if (*fb > *fa) {
		SHFT(dum,*ax,*bx,dum)
		SHFT(dum,*fb,*fa,dum)
	}
	*cx=(*bx)+GOLD*(*bx-*ax);
	*fc=(*func)(*cx);
	while (*fb > *fc) {
		r=(*bx-*ax)*(*fb-*fc);
		q=(*bx-*cx)*(*fb-*fa);
		u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
           (2.0*SIGN((tmpd=fabs(q-r))>TINY ? tmpd : TINY,q-r));   //Original: (2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
		ulim=(*bx)+GLIMIT*(*cx-*bx);
		if ((*bx-u)*(u-*cx) > 0.0) {
			fu=(*func)(u);
			if (fu < *fc) {
				*ax=(*bx);
				*bx=u;
				*fa=(*fb);
				*fb=fu;
				return;
			} else if (fu > *fb) {
				*cx=u;
				*fc=fu;
				return;
			}
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u);
		} else if ((*cx-u)*(u-ulim) > 0.0) {
			fu=(*func)(u);
			if (fu < *fc) {
				SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
				SHFT(*fb,*fc,fu,(*func)(u))
			}
		} else if ((u-ulim)*(ulim-*cx) >= 0.0) {
			u=ulim;
			fu=(*func)(u);
		} else {
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u);
		}
		SHFT(*ax,*bx,*cx,u)
		SHFT(*fa,*fb,*fc,fu)
	}
}
#undef GOLD
#undef GLIMIT
#undef TINY
#undef SHFT
#undef SIGN

#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d)  {(a)=(b);(b)=(c);(c)=(d);}
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
static double brent(double ax, double bx, double cx, double (*f)(double), double tol, double *xmin) {
	int iter;
   double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
   double e=0.0;

	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=(*f)(x);
   for (iter=0;iter<ITMAX;iter++) {
		xm=0.5*(a+b);
		tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
			*xmin=x;
			return fx;
		}
		if (fabs(e) > tol1) {
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0) p = -p;
			q=fabs(q);
			etemp=e;
			e=d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
			else {
				d=p/q;
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
					d=SIGN(tol1,xm-x);
			}
		} else {
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		}
		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		fu=(*f)(u);
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			SHFT(v,w,x,u)
			SHFT(fv,fw,fx,fu)
		} else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			} else if (fu <= fv || v == x || v == w) {
				v=u;
				fv=fu;
			}
		}
	}
   fn_DisplayError("The maximum number of iterations is reached before convergence in brent()");
	*xmin=x;
	return fx;
}
#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT
#undef SIGN



