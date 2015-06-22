
#ifndef __RANDOM__
#define __RANDOM__
   #include <math.h>
   #include <time.h>
   #include "tzmatlab.h"                        // Only when used with the MATLAB interface in gammrnd().


   //==========================================================
   // I. My own library
   //==========================================================
//   long initialize_generator(long init);
   long initialize_generator(long init);
   double unirnd(void);
   double gaussrnd(void);
   double gammrnd(double a);
   double cumulative_normal(double x);



   //==========================================================
   // II. IMSL Library
   //==========================================================
//   #define SetGlobleSeed(seed)     imsls_random_seed_set( seed )
//      //Initializes a random seed for the IMSL random number generator.
//      //seed: an integer number from 0 to 2147483646 inclusive.  If seed==0, a value is computed using the system clock.
//   #define SetGenerator(indicator)     imsls_random_option( indicator )
//      //Select a random number generator from the IMSL.
//      //indicator: an integer number from 1 to 7.  If indicator==7, the GFSR algorithm is used.

   void InitializeGlobleSeed(int seednumber);
   double UniformDouble(void);
   double StandardNormalDouble(void);
   double GammaDouble(const double a, const double b);
   double ChisquareDouble(const double v);
   void StandardNormalVector(TSdvector *x_dv);
   void StandardNormalMatrix(TSdmatrix *X_dm);
   void GammaMatrix(TSdmatrix *X_dm, TSdmatrix *A_dm, TSdmatrix *B_dm);
#endif
