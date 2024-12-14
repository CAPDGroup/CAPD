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

/* CVS $Id: matinv_aprx.hpp,v 1.14 2014/01/30 17:23:39 cxsc Exp $ */

#ifndef _CXSC_MATINVAPRX_HEADER
#define _CXSC_MATINVAPRX_HEADER

#include <rmatrix.hpp>
#include <cmatrix.hpp>
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef Unchanged
  #define dgetri_ dgetri
  #define dgetrf_ dgetrf
  #define zgetri_ ZGETRI
  #define zgetrf_ ZGETRF
#elif defined(Add_)
  //Nothing to do
#elif defined(Uppercase)
  #define dgetri_ DGETRI
  #define dgetrf_ DGETRF
  #define zgetri_ ZGETRI
  #define zgetrf_ ZGETRF
#endif

//Declaration of LAPACK-routines
extern "C" {
  void dgetri_(int *n, double* A, int *lda, int* ipiv, double* work, int *lwork, int *info);
  void dgetrf_(int *m, int *n, double* A, int *lda, int* ipiv, int *info);
  void sgetri_(int *n, float* A, int *lda, int* ipiv, float* work, int *lwork, int *info);
  void sgetrf_(int *m, int *n, float* A, int *lda, int* ipiv, int *info);
  void zgetri_(int *n, double* A, int *lda, int* ipiv, double* work, int *lwork, int *info);
  void zgetrf_(int *m, int *n, double* A, int *lda, int* ipiv, int *info);
} 


namespace cxsc {

//! Computes an approximate inverse of a real matrix
/*
    This function computes an approximate inverse of a real matrix using
    LAPACK-Routines. First the given Matrix A of type rmatrix must be 
    converted into a double array which can be used by the LAPACK routines.
    After the computation the result ie copied back into an rmatrix.
    
    \param A The matrix to be inverted
    \param R Will be overwritten by the approximate inverse
    \param error Errorcode, 1 if inversion failed, else 0
*/
inline void MatInvAprx( rmatrix &A, rmatrix &R, int& error ) 
{
#ifdef CXSC_USE_LAPACK

   int n = RowLen(A);
   int nb = 256;
   int ld = n;
   int info;
   int lwork = n*nb;
   
   double *DA = new double[n*n];
   double *work = new double[lwork];
   int *ipiv  = new int[n];
   
   // Copy real Matrix into double array 
   for(int i=1 ; i<=n ; i++) {
      for(int j=1 ; j<=n ; j++) {
         DA[(j-1)*n+(i-1)] = _double(A[Lb(A,1)+i-1][Lb(A,2)+j-1]);
      }
   }        
   
   error = 0;
   
   //LU-decomp.
   dgetrf_(&n, &n, DA, &ld, ipiv, &info);
   
   if(info != 0) {
     error = 1;
     delete[] DA;
     delete[] ipiv;
     delete[] work;
     return;
   }

   //Calculation of Inverse
   dgetri_(&n, DA, &ld, ipiv, work, &lwork, &info);
   
   //Copy computed inverse into rmatrix R
   for(int i=1 ; i<=n ; i++) {
      for(int j=1 ; j<=n ; j++) {
         R[Lb(R,1)+j-1][Lb(R,2)+i-1] = DA[(i-1)*n+(j-1)];
      }
   }        

   if(info != 0) {
     error = 1;
   }

   //delete double arrays
   delete[] DA;
   delete[] ipiv;
   delete[] work;
   
   return;  


#else

  const real
    Tiny = 1E-200;          // A divisor less than 'Tiny' is handled like zero

  real    Accu;                      
  real    Max, Temp;                       // Help variables
  int     p1 = Lb(A,1), q1 = Lb(A,2);      // Lower bounds of 'A'.
  int     pn = Ub(A,1), qm = Ub(A,2);      // Upper bounds of 'A'.
  int     n = pn-p1+1;                     // Length of the rows of 'A'
  int     m = qm-q1+1;                     // Length of the columns of 'A'
  int     i, j, k, l, kk;                  // For loops

  rmatrix LU(n,n);
  LU = A;

  error = 0;                           // Initialization

  if (n != m)                              // Error: 'A' is not square
    {error = 1; return; }


  if (n == 2) {                            // Special case: (2,2)-matrix
    Accu = 0.0;                            // Compute the determinant of 'A'
    Accu += LU[1][1]*LU[2][2] - LU[2][1]*LU[1][2];

    if (abs(Accu) < Tiny)             // Pivot element < 'Tiny', inversion
    { error = 2; return; }            // failed, matrix 'A' assumed to be
                                      // singular

    R[p1][q1] =  LU[2][2] / Accu;  R[p1][qm] = -LU[1][2] / Accu;
    R[pn][q1] = -LU[2][1] / Accu;  R[pn][qm] =  LU[1][1] / Accu;
    return;
  }

  // Start usual case: Dimension of 'A' > 2
  //---------------------------------------
  rvector v(n), x(n);          // Help vectors
  int*    p = new int [n+1];   // Dynamic array of type integer
                               // Note: p[0] not used !

  for (i = 1; i <= n; i++)  p[i] = i; // Initializing permutation vector

  // Start LU factorization
  //-----------------------
  i = 1;
  while ( (error == 0) && (i <= n) ) {
    // Compute the numerators of those elements of the i-th column
    // of L which are not updated so far and store them in 'v'.
    //------------------------------------------------------------
    #ifdef _OPENMP
    #pragma omp parallel for firstprivate(Accu) private(k,j) default(shared)
    #endif
    for (k = i; k <= n; k++) {   
      Accu = LU[k][i];            
      for (j = 1; j < i; j++)
        Accu -= LU[k][j]*LU[j][i];
      v[k] = Accu;
    }

    // Look for the column pivot
    //--------------------------
    j = i;  Max = abs(v[i]);
    for (k = i+1; k <= n; k++)
      if ( (Temp = abs(v[k])) > Max ) { j = k; Max = Temp; }

    // Swap rows of 'A' and 'v', store the permutation in 'p'
    //-------------------------------------------------------
    if (j != i) {
      x = LU[i];    LU[i] = LU[j];  LU[j] = x;
      k = p[i];     p[i] = p[j];  p[j] = k;
      Temp = v[i];  v[i] = v[j];  v[j] = Temp;
    }

    if (Max < Tiny)                   // Pivot element < 'Tiny', inversion
      { error = 2; return; }         // failed, matrix 'A' assumed to be
                                      // singular
    Temp = v[i];
    LU[i][i] = Temp;    // Update the diagonal element of U

    // Update U's row and L's column elements
    //---------------------------------------
    #ifdef _OPENMP
    #pragma omp parallel for firstprivate(Accu) private(k,j) default(shared)
    #endif
    for (k = i+1; k <= n; k++) { 
      Accu = LU[i][k];            
      for (j = 1; j < i; j++) 
        Accu -= LU[i][j] * LU[j][k];
      LU[i][k] = Accu;
      LU[k][i] = v[k] / Temp;
    }
    i++;
  } // while

  // Now 'A' is overwritten with the subdiagonal elements of L in its
  // lower left triangle and with the elements of U in its diagonal and
  // its upper right triangle. The elements of the inverse matrix are
  // computed column by column using forward/backward substitution.
  //-------------------------------------------------------------------
  #ifdef _OPENMP
  #pragma omp parallel for firstprivate(Accu,x) private(i,j,k,l,kk) default(shared)
  #endif
  for (k = 1; k <= n; k++) {
    // Forward substitution: L*x = P*e_k, where e_k is the k-th unit
    // vector. Note: If P*e_k has m leading zeros, this results in
    // x_i = 0 for 1,..,l-1 and x_l = 1. Thus, forward substitution
    // starts at index l+1.
    //--------------------------------------------------------------
    l = 1;
    while (p[l] != k) { x[l] = 0.0; l++; }
    x[l] = 1.0;
    for (i = l+1; i <= n; i++) {   
      Accu = 0.0;                  
      for (j = l; j < i; j++) 
         Accu += LU[i][j] * x[j];
      x[i] = -Accu;
    }

    // Backward substitution: U * x = x, where the right-hand side is
    // the result of the forward substitution. It will be overwritten
    // by the solution of the system, the k-th column of the inverse
    // matrix.
    //---------------------------------------------------------------
    kk = q1+k-1;                 // Index offset for the k-th column of R
    for (i = n; i >= 1; i--) {  
      Accu = x[i];              
      for (j = i+1; j <= n; j++) 
        Accu -= LU[i][j]*x[j];
      x[i] = Accu / LU[i][i];
      R[p1+i-1][kk] = x[i];      // Remember index offset !
    }
  } // for (k = 1,...

  delete [] p;   // Free dynamically allocated memory

#endif //ifdef CXSC_USE_LAPACK
}



//! Computes an approximate inverse of a complex matrix
/*
    This function computes an approximate inverse of a complex matrix using
    LAPACK-Routines. First the given Matrix A of type cmatrix must be 
    converted into a double array which can be used by the LAPACK routines.
    After the computation the result is copied back into an rmatrix.
    
    \param A The matrix to be inverted
    \param R Will be overwritten by the approximate inverse
    \param error Errorcode, 1 if inversion failed, else 0
*/
inline void MatInvAprx( cmatrix& A, cmatrix &R, int& error ) {
#ifdef CXSC_USE_LAPACK

   int n = RowLen(A);
   int nb = 256;
   int ld = n;
   int info;
   int lwork = n*nb;
   
   //Arrays for LAPACK. Real an imag. part of a complex number are stored in 
   //succeding elements of the double array.
   double *DA = new double[2*n*n];
   double *work = new double[2*lwork];
   int *ipiv  = new int[2*n];
   
   //Copy A into double array
   for(int i=1 ; i<=n ; i++) {
      for(int j=1 ; j<=n ; j++) {
         int ind1 = 2*((j-1)*n+(i-1));
         int ind2 = ind1 + 1;              
         DA[ind1] = _double(Re(A[Lb(A,1)+i-1][Lb(A,2)+j-1]));
         DA[ind2] = _double(Im(A[Lb(A,1)+i-1][Lb(A,2)+j-1]));         
      }
   }        
   
   error = 0;
   
   //LU-decomp.
   zgetrf_(&n, &n, DA, &ld, ipiv, &info);
   
   if(info != 0) {
     error = 1;
     delete[] DA;
     delete[] ipiv;
     delete[] work;
     return;
   }
   
   //Inversion
   zgetri_(&n, DA, &ld, ipiv, work, &lwork, &info);
   
   //Copy result into rmatrix R
   for(int i=1 ; i<=n ; i++) {
      for(int j=1 ; j<=n ; j++) {
         int ind1 = 2*((i-1)*n+(j-1));
         int ind2 = ind1 + 1;             
         SetRe(R[Lb(R,1)+j-1][Lb(R,2)+i-1], DA[ind1]);
         SetIm(R[Lb(R,1)+j-1][Lb(R,2)+i-1], DA[ind2]);
      }
   }        

   if(info != 0) {
     error = 1;
   }

   //Delete data for LAPACK-routines
   delete[] DA;
   delete[] ipiv;
   delete[] work;
   
   return;
 
#else

  const real
    Tiny = 1E-200;          // A divisor less than 'Tiny' is handled like zero

  complex    Accu;                      
  real    Max, Temp;                       // Help variables
  complex cTemp;
  int     p1 = Lb(A,1), q1 = Lb(A,2);      // Lower bounds of 'A'.
  int     pn = Ub(A,1), qm = Ub(A,2);      // Upper bounds of 'A'.
  int     n = pn-p1+1;                     // Length of the rows of 'A'
  int     m = qm-q1+1;                     // Length of the columns of 'A'
  int     i, j, k, l, kk;                  // For loops

  cmatrix LU(n,n);
  LU = A;

  error = 0;                           // Initialization

  if (n != m)                              // Error: 'A' is not square
    {error = 1; return; }


  if (n == 2) {                            // Special case: (2,2)-matrix
    Accu = 0.0;                            // Compute the determinant of 'A'
    Accu += LU[1][1]*LU[2][2] - LU[2][1]*LU[1][2];

    if (abs(Accu) < Tiny)             // Pivot element < 'Tiny', inversion
    { error = 2; return; }            // failed, matrix 'A' assumed to be
                                      // singular

    R[p1][q1] =  LU[2][2] / Accu;  R[p1][qm] = -LU[1][2] / Accu;
    R[pn][q1] = -LU[2][1] / Accu;  R[pn][qm] =  LU[1][1] / Accu;
    return;
  }

  // Start usual case: Dimension of 'A' > 2
  //---------------------------------------
  cvector v(n), x(n);          // Help vectors
  int*    p = new int [n+1];   // Dynamic array of type integer
                               // Note: p[0] not used !

  for (i = 1; i <= n; i++)  p[i] = i; // Initializing permutation vector

  // Start LU factorization
  //-----------------------
  i = 1;
  while ( (error == 0) && (i <= n) ) {
    // Compute the numerators of those elements of the i-th column
    // of L which are not updated so far and store them in 'v'.
    //------------------------------------------------------------
    #ifdef _OPENMP
    #pragma omp parallel for firstprivate(Accu) private(k,j) default(shared)
    #endif
    for (k = i; k <= n; k++) {   
      Accu = LU[k][i];            
      for (j = 1; j < i; j++)
        Accu -= complex(Re(LU[k][j])*Re(LU[j][i])-Im(LU[k][j])*Im(LU[j][i]), 
                        Im(LU[k][j])*Re(LU[j][i])+Re(LU[k][j])*Im(LU[j][i]));
      v[k] = Accu;
    }

    // Look for the column pivot
    //--------------------------
    j = i;  Max = abs(v[i]);
    for (k = i+1; k <= n; k++)
      if ( (Temp = abs(v[k])) > Max ) { j = k; Max = Temp; }

    // Swap rows of 'A' and 'v', store the permutation in 'p'
    //-------------------------------------------------------
    if (j != i) {
      x = LU[i];    LU[i] = LU[j];  LU[j] = x;
      k = p[i];     p[i] = p[j];  p[j] = k;
      cTemp = v[i];  v[i] = v[j];  v[j] = cTemp;
    }

    if (Max < Tiny)                   // Pivot element < 'Tiny', inversion
      { error = 2; return; }         // failed, matrix 'A' assumed to be
                                      // singular
    cTemp = v[i];
    LU[i][i] = cTemp;    // Update the diagonal element of U

    // Update U's row and L's column elements
    //---------------------------------------
    #ifdef _OPENMP
    #pragma omp parallel for firstprivate(Accu) private(k,j) default(shared)
    #endif
    for (k = i+1; k <= n; k++) { 
      Accu = LU[i][k];            
      for (j = 1; j < i; j++) 
        Accu -= complex(Re(LU[i][j])*Re(LU[j][k])-Im(LU[i][j])*Im(LU[j][k]),
                        Im(LU[i][j])*Re(LU[j][k])+Re(LU[i][j])*Im(LU[j][k]));
      LU[i][k] = Accu;
      LU[k][i] = v[k] / cTemp;
    }
    i++;
  } // while

  // Now 'A' is overwritten with the subdiagonal elements of L in its
  // lower left triangle and with the elements of U in its diagonal and
  // its upper right triangle. The elements of the inverse matrix are
  // computed column by column using forward/backward substitution.
  //-------------------------------------------------------------------
  #ifdef _OPENMP
  #pragma omp parallel for firstprivate(Accu,x) private(i,j,k,l,kk) default(shared)
  #endif
  for (k = 1; k <= n; k++) {
    // Forward substitution: L*x = P*e_k, where e_k is the k-th unit
    // vector. Note: If P*e_k has m leading zeros, this results in
    // x_i = 0 for 1,..,l-1 and x_l = 1. Thus, forward substitution
    // starts at index l+1.
    //--------------------------------------------------------------
    l = 1;
    while (p[l] != k) { x[l] = 0.0; l++; }
    x[l] = 1.0;
    for (i = l+1; i <= n; i++) {   
      Accu = 0.0;                  
      for (j = l; j < i; j++) 
         Accu += complex(Re(LU[i][j])*Re(x[j])-Im(LU[i][j])*Im(x[j]),
                         Im(LU[i][j])*Re(x[j])+Re(LU[i][j])*Im(x[j]));
      x[i] = -Accu;
    }

    // Backward substitution: U * x = x, where the right-hand side is
    // the result of the forward substitution. It will be overwritten
    // by the solution of the system, the k-th column of the inverse
    // matrix.
    //---------------------------------------------------------------
    kk = q1+k-1;                 // Index offset for the k-th column of R
    for (i = n; i >= 1; i--) {  
      Accu = x[i];              
      for (j = i+1; j <= n; j++) 
        Accu -=  complex(Re(LU[i][j])*Re(x[j])-Im(LU[i][j])*Im(x[j]),
                         Im(LU[i][j])*Re(x[j])+Re(LU[i][j])*Im(x[j]));
      x[i] = Accu / LU[i][i];
      R[p1+i-1][kk] = x[i];      // Remember index offset !
    }
  } // for (k = 1,...

  delete [] p;   // Free dynamically allocated memory

#endif //ifdef CXSC_USE_LAPACK

}


} //namespace cxsc

#endif

