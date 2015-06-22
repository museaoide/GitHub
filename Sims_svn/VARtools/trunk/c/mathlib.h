#ifndef __MATHLIB_H__
#define __MATHLIB_H__
   #include <math.h>
   #include <memory.h>

   #include "tzmatlab.h"

   //------------------------------------------------------
   // LAPACK routines -- all based on Intel MKL (or IMSL C Math library).
   //------------------------------------------------------
   int lurgen(TSdmatrix *lu_dm, TSivector *pivot_dv, const TSdmatrix *x_dm);
   int eigrsym(TSdvector *eval_dv, TSdmatrix *eVec_dm, const TSdmatrix *S_dm);
   int invrtri(TSdmatrix *X_dm, TSdmatrix *A_dm, const char un);
   int eigrgen(TSdzvector *vals_dzv, TSdzmatrix *rights_dzm, TSdzmatrix *lefts_dzm, const TSdmatrix *x_dm);
   int chol(TSdmatrix *D_dm, TSdmatrix *S_dm, const char ul);
   int invrgen(TSdmatrix *X_dm, TSdmatrix *A_dm);
      //Inverse of a general real matrix A.
      //If X=A, A (and X) will be replaced by inv(A).
   int BdivA_rrect(TSdmatrix *X_dm, const TSdmatrix *B_dm, const char lr, const TSdmatrix *A_dm);
   int BdivA_rgens(TSdmatrix *X_dm, const TSdmatrix *B_dm, const char lr, const TSdmatrix *A_dm);
   void Aldivb_spd(TSdvector *x_dv, TSdmatrix *A_dm, TSdvector *b_dv, char an);
   double logdeterminant(TSdmatrix *A_dm);
   //
   //void eig_rgen_all(double *eval_v, double *evec_m, const double *x_m, const int _n);
   int chol_decomp(double *D, const double *x_m, const int _n, const char ul);
   int eigrgen_decomp(double *evalr_v, double *evali_v, double *revecr_m, double *reveci_m,  double *levecr_m, double *leveci_m, const double *x_m, const int _n);
   int eigrsym_decomp(double *eval_v, double *evec_m, const double *s_m, const int _n, const char ul);
   int inv_spd(double *D, const double *s_m, const int _n, const char ul);



   //------------------------------------------------------
   // BLAS routines -- all based on Intel MKL (or IMSL C Math library).
   //------------------------------------------------------
   double VectorDotVector(TSdvector *x1_dv, TSdvector *x2_dv);
      //Output: Return sum(x1[i] * x2[i]) over i=1, ..., n.
      //Inputs:
      //  x1_dv: _n-by-1 double vector.
      //  x2_dv: _n-by-1 double vector.
   void ScalarTimesVectorUpdate(TSdvector *x2_dv, const double _alpha, TSdvector *x1_dv);
      //Output:  x2 = alpha * x1 + x2;
      //Inputs:
      //  alpha:  a double scalar;
      //  x1:  n-by-1 double vector.
   void ScalarTimesVector(TSdvector *x_dv, const double _alpha, TSdvector *a_dv, const double _beta);
      //Output: x_dv = alpha*a_dv + beta*x_dv where x_dv is n-by-1.
      //  When beta=0.0 and x_dv->v = a_dv->v, x_dv->v will be replaced by new values.
      //Inputs:
      //  a_dv: n-by-1.
      //  _alpha: a double scalar.
      //  _beta: a double scalar.
   void VectorPlusMinusVectorUpdate(TSdvector *x_dv, const TSdvector *b_dv, double _alpha);
      //Output: x_dv = _alpha * b_dv + x_dv where x_dv is _n-by-1.
      //Inputs:
      //  b_dv: _n-by-1 double vector.
      //  _alpha: double scalar.
   void VectorPlusMinusVector(TSdvector *x_dv, const TSdvector *a_dv, const TSdvector *b_dv, double _alpha);
      //???????? NOT finished yet.
      //????????Must add _beta for x_dv = alpha*a_dv + beta*b_dv.
      //??????????? NOT fully tested yet.
      //Output: x_dv = a_dv + _alpha * b_dv where x_dv is _n-by-1.
      //Inputs:
      //  a_dv: _n-by-1 double vector.
      //  b_dv: _n-by-1 double vector.
      //  _alpha: double scalar.
   void VectorTimesSelf(TSdmatrix *C_dm, const TSdvector *a_dv, const double _alpha, const double _beta, const char ul);
      //Switch between MKL and my own C code.
      //Output is the matrix C and all other arguments are inputs.
      //Computes C = alpah*a*a' + beta*C where
      //  a is m-by-1,
      //  C is m-by-m symmetric matrix,
      //  alpha: a double scalar,
      //  beta: a double scalar.
      //  ul: if=='u' or 'U', only the upper triangular part of C is to be referenced; otherwise, only the lower triangular part of C is to be referenced;
   void VectorTimesVector(TSdmatrix *C_dm, const TSdvector *a_dv, const TSdvector *b_dv, const double _alpha, const double _beta);
      //Output is C and all other arguments are inputs.
      //Computes C = alpah*a*b + beta*C where
      //  a is m-by-1,
      //  b is 1-by-n,
      //  C is m-by-n symmetric matrix,
      //  alpha: a double scalar,
      //  beta: a double scalar.
   void MatrixPlusMinusMatrixUpdate(TSdmatrix *X_dm, TSdmatrix *A_dm, double _alpha);
      //$$$$$ If A_dm or X_dm is only upper or lower symmetric, it will be always converted to a general (and symmetric) matrix.  $$$$$$
      //Output: X =_alpha * A + X where X_dm is an m-by-n general (and possibly symmetric) matrix.
      //Inputs:
      //  A_dm: m-by-n general or symmetric matrix.
      //  _alpha: double scalar.
   void MatrixTimesVector(TSdvector *x_dv, TSdmatrix *A_dm, const TSdvector *b_dv, const double _alpha, const double _beta, const char tn);
      //Output: x_dv->v = _alpha*A_dm'*b_dv + _beta*x_dv  for tn=='T'; x_dv = _alpha*A_dm*b_dv + _beta*x_dv  for tn=='N'
      //  where x_dv->v is ncols-by-1 or nrows-by-1 and needs not be initialized if _beta is set to 0.0.
      //Inputs:
      //  A_dm->M: nrows-by-ncols;
      //  b_dv->v: nrows-by-1 or ncols-by-1;
      //  _alpha: double scalar;
      //  _beta:  double scalar;
      //  tn: if =='t' or 'T', transpose of A_dm is used; otherwise, A_dm itself (no transpose) is used.
   void TrimatrixTimesVector(TSdvector *x_dv, TSdmatrix *A_dm, TSdvector *b_dv, const char tn, const char un);
      //Output: x_dv = A_dm'*b_dv  for tn=='T'; x_dv = A_dm*b_dv  for tn=='N' where x_dv->v is _n-by-1.
      //  If x_dv = b_dv (which gives the fastest return, so try to use this option), x_dv will be relaced by A*b or A'*b.
      //Inputs:
      //  A_dm->M: _n-by-_n triangular matrix.
      //  b_dv->v: _n-by-1 vector.
      //  tn: if =='T' or 't', transpose of A_dm is used; otherwise, A_dm itself (no transpose) is used.
      //  un: if =='U' or 'u', A_dm is unit triangular; otherwise, A_dm is non-unit triangular (i.e., a regular triangular matrix).
   void VectorTimesMatrix(TSdvector *x_dv, const TSdvector *b_dv, TSdmatrix *A_dm, const double _alpha, const double _beta, const char tn);
      //Output: x_dv->v = _alpha*b_dv*A_dm + _beta*x_dv  for tn=='n'; x_dv = _alpha*b_dv*A_dm' + _beta*x_dv  for tn=='t'
      //  where x_dv->v is 1-by-ncols or 1-by-nrows and needs not be initialized if _beta is set to 0.0.
      //Inputs:
      //  A_dm->M: nrows-by-ncols;
      //  b_dv->v: 1-by-nrows or 1-by-ncols;
      //  _alpha: double scalar;
      //  _beta:  double scalar;
      //  tn: if =='T' or 't', transpose of A_dm is used; otherwise (=='N' or 'n'), A_dm itself (no transpose) is used.
   void ScalarTimesMatrix(TSdmatrix *x_dm, const double _alpha, TSdmatrix *a_dm, const double _beta);
      //$$$$$ If a_dm or x_dm (when _beta!=0) is only upper or lower symmetric, it will be always converted to a general (and symmetric) matrix.  $$$$$$
      //Output: x_dm = alpha*a_dm + beta*x_dm where x_dm is m-by-n.
      //  When beta=0.0 and x_dm->M = a_dm->M, x_dm->M will be replaced by new values.
      //Inputs:
      //  a_dm: m-by-n.
      //  _alpha: a double scalar.
      //  _beta: a double scalar.
   void MatrixTimesMatrix(TSdmatrix *C_dm, TSdmatrix *A_dm, TSdmatrix *B_dm, const double _alpha, const double _beta, const char tn1, const char tn2);
      //?????????? NOT fully tested for all possibilities.
      //Output is C and all other arguments are inputs.
      //Computes C = alpah*op(A)*op(B) + beta*C where op() is either transpose or not, depending on 't' or 'n',
      //  op(A) is m-by-k,
      //  op(B) is k-by-n,
      //  C is m-by-n,
      //  alpha is a double scalar,
      //  beta is a double scalar.
      //  tn1: if == 'T' or 't', the transpose of A is used; otherwise (== 'N' or 'n"), A itself (no transpose) is used.
      //  tn2: if == 'T' or 't', the transpose of B is used; otherwise (== 'N' or 'n"), B itself (no transpose) is used.
   void SolveTriSysVector(TSdvector *x_dv, const TSdmatrix *T_dm, TSdvector *b_dv, const char tn, const char un);
      //Output --- computes x_dv = inv(T_dm)*b_dv by solving a triangular system of equation T_dm * x_dv = b_dv.
      //  x_dv(_n-by-1) = inv(T_dm)*b_v.  If x_dv->v = b_dv->v, x_dv->v will be replaced by new values.




//   #define ScalarTimesVector(x_v, a, b_v, _n)       cblas_daxpy(_n, a, b_v, 1, x_v, 1)
//      //Output: x_v = a * b_v + x_v where double *x_v (_n-by-1) must be initialized.
//      //Inputs:  a -- double scalar;  b_v -- pointer (_n-by-1) to double.
//   #define VectorDotVector(a_v, b_v, _n)       cblas_ddot(_n, a_v, 1, b_v, 1)
//      //Output:  x=a_v'*b_v: double scalar.
//      //Inputs:  a_v, b_v: pointer (_n-by-1) to double.

   void SymmetricMatrixTimesVector(double *x_v, const double a, const double *A_m, const double *a_v, const double b,  const int _n, const char ul);
      //Output: x_v = a*A_m*a_v + b*X_m where x_v (_n-by-1) must be allocated (but needs not be initialized).
      //Inputs:
      //  A_m: _n-by-_n symmetric matrix;
      //  a_v: _n-by-1;
      //  a, b: scalars;
      //  ul: if =='u' or 'U', upper triangular elements in A_m are filled; if =='l' or 'L', lower triangular elements in A_m are filled.
   void SolveTriangularSystemVector(double *x_v, const double *A_m, const double *b_v, const int _n, const char ul, const char tn, const char un);
      //Outputs:
      //  x_v(_n-by-1) = inv(A_m)*b_v.  If x_v=b_v, b_v will be overwritten by x_v.
      //-------
      //Inputs:
      //  A_m: _n-by-_n upper or lower triangular matrix;
      //  b_v: _n-by-1 vector.
      //  ul: if =='u' or 'U', A_m is upper triangular; if =='l' or 'L', A_m is lower triangular.
      //  tn: if =='t' or 'T', A_m' (transpose), instead of A_m, will be used; if =='n', A_m itself (no transpose) will be used.
      //  un: if =='u' or 'U', A_m is a unit upper triangular (i.e., the diagonal being 1);
      //      if =='n' or 'N', A_m is a non-unit upper triangular.
      //
      // Computes x_v = inv(A_m)*b_v by solving a triangular system of equation A_m * x_v = b_v.
      // Note I: Intel MLK cblas_dtrsv() does not test for singularity or near-singulariy of the system.
      //   Such tests must be performed before calling this BLAS routine.
      // Note II: if x_v=b_v, b_v will be overwritten by x_v.






   //------------------------------------------------------
   // MKL Vector Mathematical Library with default using my own routines.
   //------------------------------------------------------
   void VectorDotDivByVector(TSdvector *x_dv, const TSdvector *a_dv, const TSdvector *b_dv);
      //????????? NOT tested yet.  06/13/03.





   //------------------------------------------------------
   // Matrix routines (my own).
   //------------------------------------------------------
   void VectorPlusVector(TSdvector *x_dv, const TSdvector *a_dv, const TSdvector *b_dv);
      //Output: x_dv = a_dv + b_dv where x_dv is _n-by-1.
      //Inputs:
      //  a_dv: _n-by-1 double vector.
      //  b_dv: _n-by-1 double vector.
   void VectorMinusVector(TSdvector *x_dv, const TSdvector *a_dv, const TSdvector *b_dv);
      //Output: x_dv = a_dv - b_dv where x_dv is _n-by-1.
      //  If x_dv = a_dv, x_dv will be replaced by x_dv - b_dv.
      //  If x_dv = b_dv, x_dv will be replaced by a_dv - x_dv.
      //Inputs:
      //  a_dv: _n-by-1 double vector.
      //  b_dv: _n-by-1 double vector.
   void VectorPlusVectorUpdate(TSdvector *x_dv, const TSdvector *b_dv);
      //Output: x_dv = b_dv +  x_dv where x_dv is _n-by-1.
      //Inputs:
      //  b_dv: _n-by-1 double vector.
   void VectorDotTimesVector(TSdvector *x_dv, const TSdvector *a_dv, TSdvector *b_dv, const double _alpha, const double _beta);
      //Output:
      //  x_dv is _n-by-1.
      //  x_dv = _alpha * a_dv .* b_dv + _beta * x_dv if x_dv != b_dv.
      //  x_dv = _alpha * a_dv .* x_dv + _beta * x_dv if x_dv = b_dv.
      //Inputs:
      //  a_dv: _n-by-1 double vector.
      //  b_dv: _n-by-1 double vector.
      //  _alpha: double scalar.
      //  _beta: a double scalar.
   void ElementwiseInverseofMatrix(TSdmatrix *Y_dm, TSdmatrix *X_dm);
   void ElementwiseInverseofVector(TSdvector *y_dv, TSdvector *x_dv);
   void SwapColsofMatrix(TSdmatrix *X_dm, int j1, int j2);
      //??????? NOT tested yet.
   void SwapColsofMatrices(TSdmatrix *X1_dm, int j1, TSdmatrix *X2_dm, int j2);
   void SwapPositionsofMatrix(TSdmatrix *X_dm, int j1, int j2);
//   void SwapMatricesofCell(TSdcell *X_dc, int c1, int c2);
//      //??????? NOT tested yet.
   void PermuteColsofMatrix(TSdmatrix *A_dm, const TSivector *indx_iv);
      //???????? NOT tested yet.  20 Oct 03.
   void PermuteRowsofMatrix(TSdmatrix *A_dm, const TSivector *indx_iv);
      //???????? Default option is NOT tested yet.  20 Oct 03.
   void PermuteMatrix(TSdmatrix *A_dm, const TSivector *indx_iv);
   void PermuteMatricesofCell(TSdcell *A_dm, const TSivector *indx_iv);
   void ScalarTimesColofMatrix(TSdvector *y_dv, double _alpha, TSdmatrix *x_dm, int _j);
          //????????? Default option, in the #else, has NOT been tested yet!
   void ScalarTimesColofMatrix2ColofMatrix(TSdmatrix *y_dm, int jy, double _alpha, TSdmatrix *x_dm, int jx);
   void ScalarTimesColofMatrixPlusVector2ColofMatrix(TSdmatrix *Y_dm, int jy, double _alpha, TSdmatrix *X_dm, int jx, double _beta, TSdvector *x_dv);
//   void ColofMatrixDotTimesVector(TSdvector *y_dv, TSdmatrix *X_dm, int jx, TSdvector *x_dv, double _alpha, double _beta);
   void MatrixDotDivideVector_row(TSdmatrix *Y_dm, TSdmatrix *X_dm, TSdvector *x_dv, double _alpha, double _beta);
   void RowofMatrixDotDivideVector(TSdvector *y_dv, TSdmatrix *X_dm, int ix, TSdvector *x_dv, double _alpha, double _beta);
      //??????? NOT tested yet, 01/02/04.
   void ColofMatrixDotTimesVector(TSdvector *y_dv, TSdmatrix *X_dm, int jx, TSdvector *x_dv, double _alpha, double _beta);
   void ColofMatrixDotTimesColofMatrix(TSdvector *y_dv, TSdmatrix *X1_dm, int jx1, TSdmatrix *X2_dm, int jx2, double _alpha, double _beta);
   void ColofMatrixDotTimesColofMatrix2ColofMatrix(TSdmatrix *Y_dm, int jy, TSdmatrix *X1_dm, int jx1, TSdmatrix *X2_dm, int jx2, double _alpha, double _beta);
   void MatrixPlusMatrixUpdate(TSdmatrix *X_dm, TSdmatrix *A_dm);
      //Output: X = X + A where X_dm is an m-by-n general matrix.
      //Inputs:
      //  A_dm: m-by-n general matrix.
      //  B_dm: m-by-n general matrix.
   void MatrixPlusMatrix(TSdmatrix *X_dm, TSdmatrix *A_dm, TSdmatrix *B_dm);
      //Output: X = A + B where X_dm is an m-by-n general matrix.
      //Inputs:
      //  A_dm: m-by-n general matrix.
      //  B_dm: m-by-n general matrix.
   void MatrixMinusMatrix(TSdmatrix *X_dm, TSdmatrix *A_dm, TSdmatrix *B_dm);
      //Output: X = A - B where X_dm is an m-by-n general matrix.
      //Inputs:
      //  A_dm: m-by-n general matrix.
      //  B_dm: m-by-n general matrix.
   void Matrix2PlusMinusMatrix(TSdmatrix *X_dm, TSdmatrix *A_dm, TSdmatrix *B_dm, TSdmatrix *C_dm, const double _alpha, const double _beta, const double _gamma);
      //????? Not yet exhaust all possibilities of alpha, beta, and gamma to get most efficiency.  Add more as required.  10 February 2003.
      //Output: X = alpha*A + beta*B + gamma*C where X_dm is an m-by-n general matrix.
      //Inputs:
      //  A_dm: m-by-n general matrix.
      //  B_dm: m-by-n general matrix.
      //  C_dm: m-by-n general matrix.
      //  _alpha: a double scalar for A_dm.
      //  _beta: a double scalar for B_dm.
      //  _gamma: a double scalar for C_dm.
   void MatrixPlusConstantDiagUpdate(TSdmatrix *X_dm, const double _alpha);
      //Output: X = X + diag([_alpha, ..., _alpha]) where X is an n-by-n square real matrix.
   void MatrixDotTimesMatrix(TSdmatrix *X_dm, TSdmatrix *A_dm, TSdmatrix *B_dm, const double _alpha, const double _beta);
      //$$$$$ If A_dm or B_dm or X_dm (when _beta!=0) is only upper or lower symmetric, it will be always converted to a general (and symmetric) matrix.  $$$$$$
      //Output:
      //  X_dm is m-by-n.
      //  X_dm = _alpha * A_dm .* B_dm + _beta * X_dm if X_dm != B_dm.
      //  X_dm = _alpha * A_dm .* X_dm + _beta * X_dm if X_dm = B_dm.
   void CopyVector0(TSdvector *x1_dv, const TSdvector *x2_dv);
   void CopyMatrix0(TSdmatrix *x1_dm, TSdmatrix *x2_dm);
   void CopyCellvec0(TSdcellvec *x1_dcv, TSdcellvec *x2_dcv);
   void CopyCell0(TSdcell *x1_dc, TSdcell *x2_dc);
   void CopySubmatrix0(TSdmatrix *x1_dm, TSdmatrix *x2_dm, const int br, const int bc, const int nrs, const int ncs);
   void CopySubmatrix(TSdmatrix *x1_dm, const int br1, const int bc1, TSdmatrix *x2_dm, const int br2, const int bc2, const int nrs, const int ncs);
   void CopySubvector(TSdvector *x1_dv, const int ptrloc1, const TSdvector *x2_dv, const int ptrloc2, const int nels);
   void CopySubmatrix2vector(TSdvector *x1_dv, const int ptrloc1, TSdmatrix *x2_dm, const int br, const int bc, const int nels);
   void CopySubmatrix2vector_sub(TSdvector *x1_dv, const int ptrloc1, TSdmatrix *x2_dm, const int br, const int bc, const int nrs, const int ncs);
   void CopySubmatrix2vector_int(TSivector *x1_iv, const int ptrloc1, TSimatrix *x2_im, const int br, const int bc, const int nels);
   void CopySubmatrix2vector_row(TSdvector *x1_dv, const int ptrloc1, TSdmatrix *x2_dm, const int br, const int bc, const int nels);
   void CopySubvector2matrix(TSdmatrix *x1_dm, const int br, const int bc, const TSdvector *x2_dv, const int ptrloc2, const int nels);
   void CopySubvector2matrix_sub(TSdmatrix *x1_dm, const int br, const int bc, const int nrs, const int ncs, TSdvector *x2_dv, const int ptrloc2);
   void CopySubvector2matrix_unr(TSdmatrix *x1_dm, const int br, const int bc, const TSdvector *x2_dv, const int ptrloc2, const int nels);
   void TransposeSquare(TSdmatrix *B_dm, TSdmatrix *A_dm);
             //???????? Some options are NOT test yet.  2/27/03.  ???????????
   void TransposeRegular(TSdmatrix *B_dm, const TSdmatrix *A_dm);
   void SUtoGE(TSdmatrix *x_dm);
      //Output: x_dm (nrows<=ncols) becomes a general matrix in addition to being upper symmetric.
      //Input: x_dm (nrows<=ncols) is upper symmetric.
   void SLtoGE(TSdmatrix *x_dm);
      //Output: x_dm (nrows>=ncols) becomes a general matrix in addition to being lower symmetric.
      //Input: x_dm (nrows>=ncols) is lower symmetric.
   double SumVector(TSdvector *x_dv);
   double MinVector(TSdvector *x_dv);
   int MaxVector_int(TSivector *x_iv);
   //+
   void diagdv(TSdvector *x_dv, TSdmatrix *x_dm);
   double tracefabs(TSdmatrix *x_dm);
   double tracelogfabs(TSdmatrix *x_dm);
   double tracelog(TSdmatrix *x_dm);





   //=== Self-written routines calling the above two routines.
   void ergodicp(TSdvector *p_dv, TSdmatrix *P_dm);
   //double *fn_ergodp2(const double *cp_m, const int _n);
   double *alloc_ergodp2(const double *cp_m, const int _n);
#endif
