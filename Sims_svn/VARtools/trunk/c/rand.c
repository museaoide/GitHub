
#include "rand.h"


//==========================================================
// I. My own library
//==========================================================
static long idum=-1;

long initialize_generator(long init) {
 //=== init can be a positive or negative integer.
 if (init > 0)
   idum=-init;
 else
   if (init < 0)
     idum=init;
    else                                // when init=0.
     idum=-(long)time((time_t *)NULL);
 return idum;
}



/*
   Returns a uniform variate, adapted from Numerical Recipes in C
   (C) Copr. 1986-92 Numerical Recipes Software $=|''1`2.
   Ran2() in Numerical Recipes in C.  The cycle period is
   about 2.3*10^18  -- practically infinite.
*/
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double unirnd(void)
{
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
   double temp;

   if (idum <= 0) {
      if (-(idum) < 1) idum=1;
      else idum = -(idum);
      idum2=(idum);
		for (j=NTAB+7;j>=0;j--) {
         k=(idum)/IQ1;
         idum=IA1*(idum-k*IQ1)-k*IR1;
         if (idum < 0) idum += IM1;
         if (j < NTAB) iv[j] = idum;
		}
		iy=iv[0];
	}
   k=(idum)/IQ1;
   idum=IA1*(idum-k*IQ1)-k*IR1;
   if (idum < 0) idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
   iv[j] = idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX


/*
   Returns a standard gaussian variate.  The density function for the
   standard gaussian is

                          1
                     ----------- exp(-0.5*x^2)
                      sqrt(2*Pi)

*/
double gaussrnd(void)
{
 static long iset=0;
 static double gset;
 double fac,r,v1,v2;

 if  (iset == 0)
   {
    do
     {
      v1=(double)(2.0*unirnd()-1.0);
      v2=(double)(2.0*unirnd()-1.0);
      r=v1*v1+v2*v2;
     }
    while (r >= 1.0);
    fac=(double)sqrt(-2.0*log(r)/r);
    gset=v1*fac;
    iset=1;
    return v2*fac;
	}
  else
   {
    iset=0;
    return gset;
	}
}

/*
   Returns a standard gamma variate.  The density function for a standard gamma
   distribution is

                                     1
                         p(x) =  --------- x^(a-1) exp(-x).
                                  gamma(a)

   where a>0 and gamma denotes the gamma function, which is defined as the
   integral with respect to t from 0 to infinity of exp(-t)*t^(z-1).

   When a==1.0, then gamma is exponential. (Devroye, page 405).
   When a<1.0, Johnk's generator (Devroye, page 418).
   When a>1.0, a rejection method or Best's algorithm (Devroye, page 410).

   A general gamma variate can be obtained as follows.  Let z=b*x.  Then,
   z is drawn from a general gamma distribution whose density is

                                      1
                         p(z) = -------------- z^(a-1) exp(-z/b).
                                 gamma(a) b^a

   Uses algorithm translated by Iskander Karibzhanov from the Matlab function
   gamrnd.m, which follows Johnk's generator in Devroye ("Non-Uniform Random
   Variate Generation", Springer-Verlag, 1986, page 418).
*/
double gammrnd(double a) {
   double b = a-1.0,
          u, v, w, x, y, z;

   if (a <= 0.0) {
      //** When used with a C MEX-file for MATLAB.
      printf("\nThe input argument x for gammrnd(x) in rand.c must be a positive double.");
      #ifdef WIN_MATLABAPI
         mexErrMsgTxt("Error! ???.");       // This terminates the entire Matlab program.
      #else
         //@ When used entirely with C.
         getchar();
         exit(1);                          // This exits the entire C program.
      #endif
   }
   if (a==1.0)
      return -log(unirnd());
   if (a < 1.0) {
      for (;;) {
         x = pow(unirnd(), 1.0/a),
         y = pow(unirnd(), 1.0/(1.0-a));
         if (x+y <= 1.0)
            return -log(unirnd())*x/(x+y);
      }
   }
   else {
      for (;;) {
         u = unirnd();
         v = unirnd();
         w = u*(1.0-u);
         y = sqrt((3.0*a-0.75)/w)*(u-0.5);
         x = b+y;
         if (x >= 0.0) {
            z = 64.0*w*w*w*v*v;
            if (z <= (1.0-2.0*y*y/x))
               return x;
            else if (log(z) <= (2.0*(b*log(x/b)-y)))
               return x;
         }
      }
   }
}



/*
   Returns the integral from -infinity to x of 1/sqrt(2*PI)*exp(-y^2/2).
   Routine adapted from Numerical Recipes in C.
*/
double cumulative_normal(double x)
{
 double z=fabs(0.7071067811865*x), t=2.0/(2.0+z);

 return (x > 0) ?
           1.0-0.5*t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+
             t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*(-1.13520398+
               t*(1.48851587+t*(-0.82215223+t*0.17087277)))))))))
                :
           0.5*t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+
             t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*(-1.13520398+
               t*(1.48851587+t*(-0.82215223+t*0.17087277)))))))));

}



//==========================================================
// II. IMSL Library
//==========================================================
void InitializeGlobleSeed(int seednumber) {
   #ifdef IMSL_RANDOMNUMBERGENERATOR             // IMSL optimization routines.
      imsls_random_option( 7 );    //Selects a random number generator from 1 to 7.  When 7, the GFSR algorithm is used.
      imsls_random_seed_set( seednumber );    //Initializes a random seed from 0 to 2147483646 inclusive.
                   //If seednumber==0, a value is computed using the system clock.
   #else CASE2_RANDOMNUMBERGENERATOR   //Imported from the C recipe book -- Case 2 and my own (Iskander) code for generating a gamma distribution.
       initialize_generator(seednumber);      //If seednumber==0, a value is computed using the system clock.
   #endif
}

double UniformDouble(void) {
   #ifdef IMSL_RANDOMNUMBERGENERATOR             // IMSL optimization routines.
      //The GFSR algorithm to generate random numbers.
      double x;
      imsls_d_random_uniform(1, IMSLS_RETURN_USER, &x, 0);
      return x;
   #else CASE2_RANDOMNUMBERGENERATOR   //Imported from the C recipe book -- Case 2 and my own (Iskander) code for generating a gamma distribution.
      return unirnd();
   #endif
}

double StandardNormalDouble(void) {
   #ifdef IMSL_RANDOMNUMBERGENERATOR             // IMSL optimization routines.
      //The GFSR algorithm to generate random numbers.
      double x;
      imsls_d_random_normal(1, IMSLS_RETURN_USER, &x, 0);
      return x;
   #else CASE2_RANDOMNUMBERGENERATOR   //Imported from the C recipe book -- Case 2 and my own (Iskander) code for generating a gamma distribution.
      return gaussrnd();
   #endif
}

double GammaDouble(const double a, const double b) {
   //The probability density of x is of the form: ( x^(a-1) exp(-x/b) )/( Gamma(a) * b^a ).
   //The same form as in the MATLAB (not Gelman et al's) notation.
   #ifdef IMSL_RANDOMNUMBERGENERATOR             // IMSL optimization routines.
      //The GFSR algorithm to generate random numbers.
      double x;
      imsls_d_random_gamma(1, a, IMSLS_RETURN_USER, &x, 0);
      if (b != 1.0)  x *= b;
      return x;
   #else CASE2_RANDOMNUMBERGENERATOR   //Imported from the C recipe book -- Case 2 and my own (Iskander) code for generating a gamma distribution.
      return ( (b==1.0) ? gammrnd(a) : b*gammrnd(a) );
   #endif
}

void StandardNormalVector(TSdvector *x_dv) {
   #ifdef CASE2_RANDOMNUMBERGENERATOR   //Imported from the C recipe book -- Case 2 and my own (Iskander) code for generating a gamma distribution.
      int _i;
   #endif
   double *v;

   if ( !x_dv ) fn_DisplayError(".../rand.c/StandardNormalVector(): Input vector must be created (memory-allocated)");
   else  v = x_dv->v;

   #if defined (IMSL_RANDOMNUMBERGENERATOR)             // IMSL optimization routines.
      //The GFSR algorithm to generate random numbers.
      imsls_d_random_normal(x_dv->n, IMSLS_RETURN_USER, v, 0);
      x_dv->flag = V_DEF;
   #else   //Default to CASE2_RANDOMNUMBERGENERATOR, imported from the C recipe book -- Case 2 and my own (Iskander) code for generating a gamma distribution.
      for (_i=x_dv->n-1; _i>=0; _i--)  v[_i] = gaussrnd();
      x_dv->flag = V_DEF;
   #endif
}


void StandardNormalMatrix(TSdmatrix *X_dm) {
   #ifndef IMSL_RANDOMNUMBERGENERATOR   //Default to the C recipe book -- Case 2 and my own (Iskander) code for generating a gamma distribution.
   int _i;
   #endif
   double *M;

   if ( !X_dm ) fn_DisplayError(".../rand.c/StandardNormalMatrix(): Input matrix must be created (memory-allocated)");
   M = X_dm->M;

   #if defined (IMSL_RANDOMNUMBERGENERATOR)             // IMSL optimization routines.
      //The GFSR algorithm to generate random numbers.
      imsls_d_random_normal(X_dm->nrows*X_dm->ncols, IMSLS_RETURN_USER, M, 0);
      X_dm->flag = M_GE;
   #else    //Default to CASE2_RANDOMNUMBERGENERATOR, imported from the C recipe book -- Case 2 and my own (Iskander) code for generating a gamma distribution.
      for (_i=X_dm->nrows*X_dm->ncols-1; _i>=0; _i--)  M[_i] = gaussrnd();
      X_dm->flag = M_GE;
   #endif
}


void GammaMatrix(TSdmatrix *X_dm, TSdmatrix *A_dm, TSdmatrix *B_dm) {
   //The probability density of each element of X_dm is of the form: ( x^(a-1) exp(-x/b) )/( Gamma(a) * b^a ).
   //The same form as in the MATLAB (not Gelman et al's) notation.
   #ifdef IMSL_RANDOMNUMBERGENERATOR             // IMSL optimization routines.
      //The GFSR algorithm to generate random numbers.
      int _i, nels, nrows, ncols;
      double *X, *A, *B, a, b;

      if ( !X_dm || !A_dm || !B_dm )  fn_DisplayError(".../rand.c/GammaMatrix():  All input matrices must be created (memory-allocated)");
      else if ( !A_dm->flag || !B_dm->flag )  fn_DisplayError(".../rand.c/GammaMatrix():  Two R input matrices must be given values");
      else {
         nrows = X_dm->nrows;
         ncols = X_dm->ncols;
         nels = nrows*ncols;
         X = X_dm->M;
         A = A_dm->M;
         B = B_dm->M;
      }

      if ( (nrows != A_dm->nrows) || (nrows != B_dm->nrows) || (ncols != A_dm->ncols) || (ncols != B_dm->ncols) )
         fn_DisplayError(".../rand.c/GammaMatrix():  Dimensions of all input matrices must match");


      if ( A_dm->flag & M_CN ) {
         imsls_d_random_gamma(nels, a=A[0], IMSLS_RETURN_USER, X, 0);  //Same for all elements of A_dm->M.
         if ( B_dm->flag & M_CN ) {
            if ( !((b = B[0])==1.0) )  cblas_dscal(nels, b, X, 1);   //Same for all elements of B_dm->M.
         }
         else {
            for (_i=nels-1; _i>=0; _i--)  X[_i] *= B[_i];
         }
      }
      else {
         for (_i=nels-1; _i>=0; _i--) {
            imsls_d_random_gamma(1, A[_i], IMSLS_RETURN_USER, &X[_i], 0);
            X[_i] *= B[_i];
         }
      }

      X_dm->flag = M_GE;
   #else CASE2_RANDOMNUMBERGENERATOR   //Imported from the C recipe book -- Case 2 and my own (Iskander) code for generating a gamma distribution.
      int _i, nels, nrows, ncols;
      double *X, *A, *B;

      if ( !X_dm || !A_dm || !B_dm )  fn_DisplayError(".../rand.c/GammaMatrix():  All input matrices must be created (memory-allocated)");
      else if ( !A_dm->flag || !B_dm->flag )  fn_DisplayError(".../rand.c/GammaMatrix():  Two R input matrices must be given values");
      else {
         nrows = X_dm->nrows;
         ncols = X_dm->ncols;
         nels = nrows*ncols;
         X = X_dm->M;
         A = A_dm->M;
         B = B_dm->M;
      }

      if ( (nrows != A_dm->nrows) || (nrows != B_dm->nrows) || (ncols != A_dm->ncols) || (ncols != B_dm->ncols) )
         fn_DisplayError(".../rand.c/GammaMatrix():  Dimensions of all input matrices must match");

      for (_i=nels-1; _i>=0; _i--)  X[_i] = B[_i]*gammrnd(A[_i]);
      X_dm->flag = M_GE;
   #endif
}


double ChisquareDouble(const double v) {
   //The probability density of x is of the form: ( x^(v/2-1) exp(-x/2) )/( Gamma(v/2) * 2^(v/2) ).
   //= GammaDouble(v/2.0, 2.0).
   #ifdef IMSL_RANDOMNUMBERGENERATOR             // IMSL optimization routines.
      //The GFSR algorithm to generate random numbers.
      double x;
      imsls_d_random_chi_squared(1, v, IMSLS_RETURN_USER, &x, 0);
      return x;
   #else CASE2_RANDOMNUMBERGENERATOR   //Imported from the C recipe book -- Case 2 and my own (Iskander) code for generating a gamma distribution.
      return ( 2.0*gammrnd(0.5*v) );
   #endif
}


/**
//void UniformDouble(double *x) {
//   #ifdef IMSL_RANDOMNUMBERGENERATOR             // IMSL optimization routines.
//      //The GFSR algorithm to generate random numbers.
//      imsls_d_random_uniform(1, IMSLS_RETURN_USER, x, 0);
//   #endif
//   #ifdef CASE2_RANDOMNUMBERGENERATOR   //Imported from the C recipe book -- Case 2 and my own (Iskander) code for generating a gamma distribution.
//      *x = unirnd();
//   #endif
//}
/**/
