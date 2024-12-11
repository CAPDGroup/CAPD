/*
**  CXSC is a C++ library for eXtended Scientific Computing (V 2.5.4)
**
**  Copyright (C) 1990-2000 Institut fuer Angewandte Mathematik,
**                          Universitaet Karlsruhe, Germany
**            (C) 2000-2014 Wiss. Rechnen/Softwaretechnologie
**                          Universitaet Wuppertal, Germany   
**
**  This library is free software; you can redistribute it and/or
**  modify it under the terms of the GNU Library General Public
**  License as published by the Free Software Foundation; either
**  version 2 of the License, or (at your option) any later version.
**
**  This library is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
**  Library General Public License for more details.
**
**  You should have received a copy of the GNU Library General Public
**  License along with this library; if not, write to the Free
**  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/* CVS $Id: cxsc_blas.hpp,v 1.11 2014/01/30 17:23:44 cxsc Exp $ */

/*
**  FastPLSS: A library of (parallel) verified linear (interval) system 
**  solvers using C-XSC (V 0.2)
**  
**  Author: Michael Zimmer
**
**  This software is based on:
**    - Module LinSys of the C-XSC-Toolbox
**      Authors: Rolf Hammer, Matthias Hocks, Dietmar Ratz
**    - Self-verifying solver for a dense system of linear equations
**      Authors: Carlos Holbig, Walter Kraemer, Paulo Sergio Morandi Junior,
**               Bernardo Frederes Kramer Alcalde, 
*/

 
#ifndef _CXSC_BLAS_HEADER_INCLUDED
#define _CXSC_BLAS_HEADER_INCLUDED

namespace cxsc {

enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102};
enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113};
enum CBLAS_UPLO {CblasUpper=121, CblasLower=122};
enum CBLAS_DIAG {CblasNonUnit=131, CblasUnit=132};
enum CBLAS_SIDE {CblasLeft=141, CblasRight=142};

//Declaration of BLAS-routines
extern "C" {
   double cblas_ddot(const int N, const double *X, const int incX,
                  const double *Y, const int incY);

   void cblas_zdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu);

   void cblas_daxpy(const int N, const double alpha, const double *X,
                 const int incX, double *Y, const int incY);

   void cblas_zaxpy(const int N, const double alpha, const double *X,
                 const int incX, double *Y, const int incY);

   void cblas_dgemv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 const double *X, const int incX, const double beta,
                 double *Y, const int incY);

   void cblas_zgemv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *X, const int incX, const void *beta,
                 void *Y, const int incY);

   void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc);

   void cblas_zgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const void *alpha, const void *A,
                 const int lda, const void *B, const int ldb,
                 const void *beta, void *C, const int ldc);

}

inline void blasdot(const rvector& x, const rvector& y, real& res);

inline void blasdot(const rvector& x, const rvector& y, interval& res);

inline void blasdot(const rvector& x, const cvector& y, complex& res);

inline void blasdot(const cvector& x, const rvector& y, complex& res);

inline void blasdot(const rvector& x, const cvector& y, cinterval& res);

inline void blasdot(const cvector& x, const rvector& y, cinterval& res);

inline void blasdot(const cvector& x, const cvector& y, complex& res);

inline void blasdot(const cvector& x, const cvector& y, cinterval& res);

inline void bsort(const ivector &x, const ivector &y, rvector& x_inf, rvector& y_inf, rvector &x_sup, rvector &y_sup, int n, int lb1, int lb2);

//Sorts inf and sup of an interval vector and a real vector into two real vector for the computation of the infimum 
//and two real vectors for the computation of the supremum. Rounding mode must be set to upwards!
inline void bsort(const ivector &x, const rvector &y, rvector& x_inf, rvector& y_inf, rvector &x_sup, rvector &y_sup, int n, int lb1, int lb2);

//Sorts inf and sup of an interval vector and a real vector into two real vector for the computation of the infimum 
//and two real vectors for the computation of the supremum. Rounding mode must be set to upwards!
inline void bsort(const rvector &y, const ivector &x, rvector& x_inf, rvector& y_inf, rvector &x_sup, rvector &y_sup, int n, int lb1, int lb2);

inline void blasdot(const ivector& x, const ivector& y, interval& res);

inline void blasdot(const ivector& x, const rvector& y, interval& res);

inline void blasdot(const rvector& x, const ivector& y, interval& res);

inline void blasdot(const rvector& x, const civector& y, cinterval& res);

inline void blasdot(const civector& x, const rvector& y, cinterval& res);

inline void blasdot(const cvector& x, const ivector& y, cinterval& res);

inline void blasdot(const ivector& x, const cvector& y, cinterval& res);

inline void blasdot(const civector& x, const civector& y, cinterval& res);

inline void blasdot(const civector& x, const cvector& y, cinterval& res);

inline void blasdot(const cvector& x, const civector& y, cinterval& res);

inline void blasdot(const civector& x, const ivector& y, cinterval& res);

inline void blasdot(const ivector& x, const civector& y, cinterval& res);

/***************************************************************************/

inline void blasmvmul(const rmatrix& A, const rvector& x, rvector& r);

inline void blasmvmul(const rmatrix& A, const rvector& x, ivector& r);

inline void blasmvmul(const rmatrix& A, const ivector& x, ivector& r);

inline void blasmvmul(const imatrix& A, const rvector& x, ivector& r);

inline void blasmvmul(const imatrix& A, const ivector& x, ivector& r);

inline void blasmvmul(const cmatrix& A, const ivector& x, civector& r);

inline void blasmvmul(const imatrix& A, const cvector& x, civector& r);

inline void blasmvmul(const rmatrix& A, const civector& x, civector& r);

inline void blasmvmul(const cimatrix& A, const rvector& x, civector& r);

inline void blasmvmul(const cmatrix& A, const civector& x, civector& r);

inline void blasmvmul(const cimatrix& A, const cvector& x, civector& r);

inline void blasmvmul(const imatrix& A, const civector& x, civector& r);

inline void blasmvmul(const cimatrix& A, const ivector& x, civector& r);

inline void blasmvmul(const cimatrix& A, const civector& x, civector& r);

inline void blasmvmul(const cmatrix& A, const cvector& x, cvector& r);

inline void blasmvmul(const cmatrix& A, const cvector& x, civector& r);

inline void blasmvmul(const rmatrix& A, const cvector& x, cvector& r);

inline void blasmvmul(const cmatrix& A, const rvector& x, cvector& r);

inline void blasmvmul(const rmatrix& A, const cvector& x, civector& r);

inline void blasmvmul(const cmatrix& A, const rvector& x, civector& r);

/***************************************************************************/

inline void blasmatmul(const rmatrix &A, const rmatrix &B, imatrix &C);

inline void blasmatmul(const rmatrix &A, const rmatrix &B, rmatrix &C);

inline void blasmatmul(const imatrix &A, const imatrix &B, imatrix &C);

inline void blasmatmul(const rmatrix &A, const imatrix &B, imatrix &C);

inline void blasmatmul(const imatrix &A, const rmatrix &B, imatrix &C);

inline void blasmatmul(const cmatrix &A, const cmatrix &B, cmatrix &C);

inline void blasmatmul(const cmatrix &A, const cmatrix &B, cimatrix &C);

inline void blasmatmul(const cmatrix &A, const rmatrix &B, cmatrix &C);

inline void blasmatmul(const rmatrix &A, const cmatrix &B, cmatrix &C);

inline void blasmatmul(const cmatrix &A, const rmatrix &B, cimatrix &C);

inline void blasmatmul(const rmatrix &A, const cmatrix &B, cimatrix &C);

inline void blasmatmul(const cmatrix &A, const imatrix &B, cimatrix &C);

inline void blasmatmul(const imatrix &A, const cmatrix &B, cimatrix &C);

inline void blasmatmul(const rmatrix &A, const cimatrix &B, cimatrix &C);

inline void blasmatmul(const cimatrix &A, const rmatrix &B, cimatrix &C);

inline void blasmatmul(const imatrix &A, const cimatrix &B, cimatrix &C);

inline void blasmatmul(const cimatrix &A, const imatrix &B, cimatrix &C);

inline void blasmatmul(const cimatrix &A, const cmatrix &B, cimatrix &C);

inline void blasmatmul(const cmatrix &A, const cimatrix &B, cimatrix &C);

inline void blasmatmul(const cimatrix &A, const cimatrix &B, cimatrix &C);

}

#endif
