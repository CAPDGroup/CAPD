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

/* CVS $Id: fastlss.hpp,v 1.23 2014/01/30 17:23:39 cxsc Exp $ */

/*
**  FastLSS: A fast verified linear (interval) system solver
**
**  Author: Michael Zimmer
**
**  Computes a verified solution of linear system by bounding the residual of a 
**  computed approximate solution using the Krawczyk-Operator. 
**
**  Theoretical background:
**  S.M. Rump: Verification Methods for Dense and Sparse Systems of Equations.
**  In: Topics in Validated Numerics, Studies in Computational Mathematics, pages
**  63-136. Elsevier, Amsterdam, 1994.
**
**  Software is partially based on:
**    - Module LinSys of the C-XSC-Toolbox
**      Authors: Rolf Hammer, Matthias Hocks, Dietmar Ratz
**    - Self-verifying solver for a dense system of linear equations
**      Authors: Carlos Holbig, Walter Kraemer, Paulo Sergio Morandi Junior,
**               Bernardo Frederes Kramer Alcalde
*/

 
#if !defined(_CXSC_FASTLSS_HPP) || (defined(_CXSC_FASTLSS_HPP) && (defined(_CXSC_LSS_INCLUDE) || defined(_CXSC_ILSS_INCLUDE) || defined(_CXSC_CLSS_INCLUDE) || defined(_CXSC_CILSS_INCLUDE)))
#define _CXSC_FASTLSS_HPP

#include <rmatrix.hpp>
#include <imatrix.hpp>
#include <cmatrix.hpp>
#include <cimatrix.hpp> 
//Computation of approximate inverse
#include <matinv_aprx.hpp>   
#include <vector>


#if !defined(_CXSC_LSS_HPP) && !defined(_CXSC_ILSS_HPP) && !defined(_CXSC_CLSS_HPP) && !defined(_CXSC_CILSS_HPP)
#define _CXSC_LSS_HPP
#define _CXSC_ILSS_HPP
#define _CXSC_CLSS_HPP
#define _CXSC_CILSS_HPP
#define _CXSC_LSS_INCLUDE
#define _CXSC_ILSS_INCLUDE
#define _CXSC_CLSS_INCLUDE
#define _CXSC_CILSS_INCLUDE
#define _CXSC_LSS_UNDEFINED
#endif

namespace cxsc {

#ifndef _CXSC_LSS_CONSTANTS_DEFINED
#define _CXSC_LSS_CONSTANTS_DEFINED
//Control constants
static const int
  LSS_ONLY_PART_ONE = 0, //Use only first part of algorithm
  LSS_ONLY_PART_TWO = 1, //Use only second part (using R1+R2 as approximate inverse)
                         //second part can solve systems with higher condition numbers
                         //but is much slower
  LSS_BOTH_PARTS    = 2; //Try part one first, if it fails try part two

/*!
 * This struct can be used to configure the verified linear system solver. If you want to change a default setting,
 * create an instance of this struct, change the desired settings, and use the instance of the 
 * struct as the last parameter of the solver function call.
 */
struct lssconfig {
  int   K;              /*! Dot product precision for residual computations. K==0 uses the long accumulator (slow, max. precision),
                         *  K==1 uses double precision and a splitting algorithm by Rump/Dekker, K>=2 uses simulated K-fold double 
                         *  precision. The default is K=2.
                         */
  bool  msg;            //! Status message output during solver run? The default is false
  int   lssparts;       /*! Solver stages to use. Possible values: LSS_ONLY_PART_ONE (use part one of algorithm with normal floating point
                         *  approximate inverse), LSS_ONLY_PART_TWO (use part two of algorithm with inverse R1+R2), LSS_BOTH_PARTS (try 
			 *  part one first, if it fails, try part two). Part one is much faster and work for condition numbers up to about 
			 *  10^15, Part two is slow but can solve systems with condition numbers of up to 10^30. The default is LSS_BOTH_PARTS
			 */
  int   threads;        /*! Number of threads to use for OpenMP (only has an effect if OpenMP ist activated during compilation).
                         *  Using a value <=0 means use all available threads (set for example by OMP_NUM_THREADS environment variable).
			 *  The default is -1.
			 */
  int   maxIterResCorr; /*! Maximum number of iterations during residual correction (not available for K=1).
                         *  Higher settings might increase result quality, but also increase time costs. The default is10.
			 */
  int   maxIterVer;     /*! Maximum number of iterations during the verification step.
                         *  For systems with condition numbers approaching 10^15 for part one or 10^30 for part two
                         *  higher than default settings might be necessary to find a verified result. The default is 5.
			 */
  bool  refinement;     /*! Perform an iterative refinement? This can sometimes improve the computed enclosure slightly with little cost.
                         *  The default is false.
			 */
  int   maxIterRef;     /*! Maximum number of iterations during the refinement step. Only has an effect if refinement is set to true. Maximum
                         *  number of iterations is only used if the stopping criterion using epsRef is not reached. The default
                         *  is 5.
			 */
  real  epsVer;         //! Epsilon for the verification step. This is used for the Epsilon inflation during verification. The default is 0.1.
  real  epsRef;         /*! Epsilon for the stopping criterion during the refinement step. The iteration stops, is the difference between to 
                         *  successive iterates is lower than epsRef. The default is 1e-5.
			 */
  bool  matrixMode;     /*! Activates matrix mode, which uses only operators for all computations. This is faster for multiple right hand sides,
                         *  especially when using K=1, but might lose some accuracy. Use matrix mode and K=1 if you want to compute a verified
			 *  inverse of a matrix! The default is false.
			 */

  //! Default constructor setting all configuration variables to their default settings
  lssconfig() : K(2),msg(false),lssparts(LSS_BOTH_PARTS),threads(-1),maxIterResCorr(10), 
                maxIterVer(5), refinement(false), maxIterRef(5), epsVer(0.1), epsRef(1e-5), matrixMode(false) {}
};
#endif

//! Translates the error codes of the solver into corresponding error messages
std::string LinSolveErrMsg(int);

#ifdef _CXSC_LSS_UNDEFINED
//Legacy versions
//! Entry function of real linear sytems solver
static inline void  lss(const cxsc::rmatrix&, cxsc::rmatrix&, cxsc::imatrix&, int&, int, bool=false, int=LSS_BOTH_PARTS, int=-1);
//! Entry function of real linear sytems solver
static inline void  lss(const cxsc::rmatrix&, cxsc::rvector&, cxsc::ivector&, int&, int, bool=false, int=LSS_BOTH_PARTS, int=-1);
//! Entry function of real linear sytems solver
static inline void  lss(const cxsc::rmatrix&, cxsc::imatrix&, cxsc::imatrix&, int&, int, bool=false, int=LSS_BOTH_PARTS, int=-1);
//! Entry function of real linear sytems solver
static inline void  lss(const cxsc::rmatrix&, cxsc::ivector&, cxsc::ivector&, int&, int, bool=false, int=LSS_BOTH_PARTS, int=-1);

//New versions using struct
//! Entry function of real linear sytems solver
static inline void  lss(const cxsc::rmatrix&, cxsc::rmatrix&, cxsc::imatrix&, int&, struct lssconfig=lssconfig());
//! Entry function of real linear sytems solver
static inline void  lss(const cxsc::rmatrix&, cxsc::rvector&, cxsc::ivector&, int&, struct lssconfig=lssconfig());
//! Entry function of real linear sytems solver
static inline void  lss(const cxsc::rmatrix&, cxsc::imatrix&, cxsc::imatrix&, int&, struct lssconfig=lssconfig());
//! Entry function of real linear sytems solver
static inline void  lss(const cxsc::rmatrix&, cxsc::ivector&, cxsc::ivector&, int&, struct lssconfig=lssconfig());

//New versions using struct and computing an inner enclosure
//! Entry function of real linear sytems solver
static inline void  lss(const cxsc::rmatrix&, cxsc::rmatrix&, cxsc::imatrix&, cxsc::imatrix&, int&, struct lssconfig=lssconfig());
//! Entry function of real linear sytems solver
static inline void  lss(const cxsc::rmatrix&, cxsc::rvector&, cxsc::ivector&, cxsc::ivector&,int&, struct lssconfig=lssconfig());
//! Entry function of real linear sytems solver
static inline void  lss(const cxsc::rmatrix&, cxsc::imatrix&, cxsc::imatrix&, cxsc::imatrix&,int&, struct lssconfig=lssconfig());
//! Entry function of real linear sytems solver
static inline void  lss(const cxsc::rmatrix&, cxsc::ivector&, cxsc::ivector&, cxsc::ivector&,int&, struct lssconfig=lssconfig());
#endif

#ifdef _CXSC_LSS_UNDEFINED
//Legacy versions
//! Entry function of interval linear sytems solver
static inline void  ilss(const cxsc::imatrix&, cxsc::imatrix&, cxsc::imatrix&, int&, int, bool=false, int=LSS_BOTH_PARTS, int=-1);
//! Entry function of interval linear sytems solver
static inline void  ilss(const cxsc::imatrix&, cxsc::ivector&, cxsc::ivector&, int&, int, bool=false, int=LSS_BOTH_PARTS, int=-1);

//New versions using struct
//! Entry function of interval linear sytems solver
static inline void  ilss(const cxsc::imatrix&, cxsc::imatrix&, cxsc::imatrix&, int&, struct lssconfig=lssconfig());
//! Entry function of interval linear sytems solver
static inline void  ilss(const cxsc::imatrix&, cxsc::ivector&, cxsc::ivector&, int&, struct lssconfig=lssconfig());

//New versions using struct and computing an inner enclosure
//! Entry function of interval linear sytems solver
static inline void  ilss(const cxsc::imatrix&, cxsc::imatrix&, cxsc::imatrix&, cxsc::imatrix&, int&, struct lssconfig=lssconfig());
//! Entry function of interval linear sytems solver
static inline void  ilss(const cxsc::imatrix&, cxsc::ivector&, cxsc::ivector&, cxsc::ivector&, int&, struct lssconfig=lssconfig());
#endif

#ifdef _CXSC_LSS_UNDEFINED
//Legacy versions
//! Entry function of complex linear sytems solver
static inline void  clss(const cxsc::cmatrix&, cxsc::cmatrix&, cxsc::cimatrix&, int&, int, bool=false, int=LSS_BOTH_PARTS, int=-1);
//! Entry function of complex linear sytems solver
static inline void  clss(const cxsc::cmatrix&, cxsc::cvector&, cxsc::civector&, int&, int, bool=false, int=LSS_BOTH_PARTS, int=-1);
//! Entry function of complex linear sytems solver
static inline void  clss(const cxsc::cmatrix&, cxsc::cimatrix&, cxsc::cimatrix&, int&, int, bool=false, int=LSS_BOTH_PARTS, int=-1);
//! Entry function of complex linear sytems solver
static inline void  clss(const cxsc::cmatrix&, cxsc::civector&, cxsc::civector&, int&, int, bool=false, int=LSS_BOTH_PARTS, int=-1);

//New versions using struct
//! Entry function of complex linear sytems solver
static inline void  clss(const cxsc::cmatrix&, cxsc::cmatrix&, cxsc::cimatrix&, int&, struct lssconfig=lssconfig());
//! Entry function of complex linear sytems solver
static inline void  clss(const cxsc::cmatrix&, cxsc::cvector&, cxsc::civector&, int&, struct lssconfig=lssconfig());
//! Entry function of complex linear sytems solver
static inline void  clss(const cxsc::cmatrix&, cxsc::cimatrix&, cxsc::cimatrix&, int&, struct lssconfig=lssconfig());
//! Entry function of complex linear sytems solver
static inline void  clss(const cxsc::cmatrix&, cxsc::civector&, cxsc::civector&, int&, struct lssconfig=lssconfig());

//New versions using struct and computing an inner enclosure
//! Entry function of complex linear sytems solver
static inline void  clss(const cxsc::cmatrix&, cxsc::cmatrix&, cxsc::cimatrix&, cxsc::cimatrix&, int&, struct lssconfig=lssconfig());
//! Entry function of complex linear sytems solver
static inline void  clss(const cxsc::cmatrix&, cxsc::cvector&, cxsc::civector&, cxsc::civector&,int&, struct lssconfig=lssconfig());
//! Entry function of complex linear sytems solver
static inline void  clss(const cxsc::cmatrix&, cxsc::cimatrix&, cxsc::cimatrix&, cxsc::cimatrix&,int&, struct lssconfig=lssconfig());
//! Entry function of complex linear sytems solver
static inline void  clss(const cxsc::cmatrix&, cxsc::civector&, cxsc::civector&, cxsc::civector&,int&, struct lssconfig=lssconfig());
#endif

#ifdef _CXSC_LSS_UNDEFINED
//Legacy versions
//! Entry function of complex interval linear sytems solver
static inline void  cilss(const cxsc::cimatrix&, cxsc::cimatrix&, cxsc::cimatrix&, int&, int, bool=false, int=LSS_BOTH_PARTS, int=-1);
//! Entry function of complex interval linear sytems solver
static inline void  cilss(const cxsc::cimatrix&, cxsc::civector&, cxsc::civector&, int&, int, bool=false, int=LSS_BOTH_PARTS, int=-1);

//New versions using struct
//! Entry function of complex interval linear sytems solver
static inline void  cilss(const cxsc::cimatrix&, cxsc::cimatrix&, cxsc::cimatrix&, int&, struct lssconfig=lssconfig());
//! Entry function of complex interval linear sytems solver
static inline void  cilss(const cxsc::cimatrix&, cxsc::civector&, cxsc::civector&, int&, struct lssconfig=lssconfig());

//New versions using struct and computing an inner enclosure
//! Entry function of complex interval linear sytems solver
static inline void  cilss(const cxsc::cimatrix&, cxsc::cimatrix&, cxsc::cimatrix&, cxsc::cimatrix&, int&, struct lssconfig=lssconfig());
//! Entry function of complex interval linear sytems solver
static inline void  cilss(const cxsc::cimatrix&, cxsc::civector&, cxsc::civector&, cxsc::civector&,int&, struct lssconfig=lssconfig());
#endif

#ifndef _CXSC_IMINUSRA_DEFINED
#define _CXSC_IMINUSRA_DEFINED

#ifdef CXSC_USE_BLAS
//Computes [C]=[I-RA] using BLAS routines
static inline void IminusRA(const rmatrix& A, const rmatrix& B, imatrix& C) {
  int rnd = getround();

   int lbA_i = Lb(A,1); int lbA_j = Lb(A,2);
   int lbB_i = Lb(B,1); int lbB_j = Lb(B,2); 
   int ubA_i = Ub(A,1); int ubA_j = Ub(A,2);
   int lbC_i = Lb(C,1); int lbC_j = Lb(C,2); 
   int ubC_j = Ub(C,2);     
   
   int m = ubA_i - lbA_i + 1;
   int n = ubA_j - lbA_j + 1;
   int o = ubC_j - lbC_j + 1;
   
   double *DA = new double[m*n];
   double *DB = new double[n*o];
   double *DC = new double[m*o];
   
   //Copy A and B into double array   
   #ifdef _OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for(int i=1 ; i<=m ; i++) {
      for(int j=1 ; j<=n ; j++) {
         DA[(i-1)*n+(j-1)] = _double(-A[lbA_i+i-1][lbA_j+j-1]);
      }
   }        

   #ifdef _OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for(int i=1 ; i<=n ; i++) {
      for(int j=1 ; j<=o ; j++) {
         DB[(i-1)*o+(j-1)] = _double(B[lbB_i+i-1][lbB_j+j-1]);         
      }
   }        


   int LDA=n,LDB=o,LDC=m;
   double alpha = 1.0, beta = 0.0;
    
   
   //Compute Infimum
   setround(-1);
   cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, o, n, alpha, DA, 
          LDA, DB, LDB, beta, DC, LDC);
   
   //Copy infimum into C
   #ifdef _OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for(int i=1 ; i<=o ; i++) {
      for(int j=1 ; j<=m ; j++) {
           C[lbC_i+j-1][lbC_j+i-1] = DC[(i-1)*m+(j-1)];
      }
      C[lbC_i+i-1][lbC_j+i-1] = DC[(i-1)*m+(i-1)] + 1.0;
   }        
   
   //Compute supremum
   setround(1);
   cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, o, n, alpha, DA, 
          LDA, DB, LDB, beta, DC, LDC);


   //copy supremum into C
   #ifdef _OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for(int i=1 ; i<=o ; i++) {
      for(int j=1 ; j<=m ; j++) {
         SetSup(C[lbC_i+j-1][lbC_j+i-1], DC[(i-1)*m+(j-1)]);
      }
      SetSup(C[lbC_i+i-1][lbC_j+i-1], DC[(i-1)*m+(i-1)] + 1.0);
   }        
  
   delete[] DA;
   delete[] DB;
   delete[] DC; 

  setround(rnd);
}

//Computes [C]=[I-R[A]] using BLAS routines
static inline void IminusRA(const rmatrix& R, const imatrix& A, imatrix& C) {
  int rnd = getround();

  C = -R*A;

  for(int i=Lb(C,1) ; i<=Ub(C,1) ; i++)
    C[i][i] += 1.0;

  setround(rnd);
}

//Computes [C]=[I-RA] using BLAS routines
static inline void IminusRA(const cmatrix& A, const cmatrix& B, cimatrix& C) {
   int rnd = getround();

   int lbA_i = Lb(A,1); int lbA_j = Lb(A,2);
   int lbB_i = Lb(B,1); int lbB_j = Lb(B,2); 
   int ubA_i = Ub(A,1); int ubA_j = Ub(A,2);
   int lbC_i = Lb(C,1); int lbC_j = Lb(C,2); 
   int ubC_j = Ub(C,2);     
   
   int m = ubA_i - lbA_i + 1;
   int n = ubA_j - lbA_j + 1;
   int o = ubC_j - lbC_j + 1;
   
   double *DAR  = new double[m*n]; //A.re
   double *DAI  = new double[m*n]; //A.Im
   double *DAIm = new double[m*n]; //-A.Im
   double *DBR = new double[n*o];  //B.re
   double *DBI = new double[n*o];  //B.Im
   double *DCR = new double[m*o];  //C.re
   double *DCI = new double[m*o];  //C.Im
   
   //Copy A and B into corresponding double arrays
   #ifdef _OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for(int i=1 ; i<=m ; i++) {
      for(int j=1 ; j<=n ; j++) {
         int ind = ((i-1)*n+(j-1));
         DAR[ind] = _double(-Re(A[lbA_i+i-1][lbA_j+j-1]));
         DAI[ind] = _double(-Im(A[lbA_i+i-1][lbA_j+j-1]));
         DAIm[ind] =  -DAI[ind];         
      }
   }        

   #ifdef _OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for(int i=1 ; i<=n ; i++) {
      for(int j=1 ; j<=o ; j++) {
         int ind = ((i-1)*o+(j-1));
         DBR[ind] = _double(Re(B[lbB_i+i-1][lbB_j+j-1]));         
         DBI[ind] = _double(Im(B[lbB_i+i-1][lbB_j+j-1]));                  
      }
   }        

   int LDA=n,LDB=o,LDC=m;
   double alpha = 1.0, beta = 0.0;
    
   
   //Compute lower bound of result    
   setround(-1);
   cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, o, n, alpha, DAR, 
          LDA, DBR, LDB, beta, DCR, LDC);
   beta = 1.0;
   cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, o, n, alpha, DAIm, 
          LDA, DBI, LDB, beta, DCR, LDC);
   beta = 0.0;
   cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, o, n, alpha, DAR, 
          LDA, DBI, LDB, beta, DCI, LDC);
   beta = 1.0;
   cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, o, n, alpha, DAI, 
          LDA, DBR, LDB, beta, DCI, LDC);

   //Copy lower bound into C  
   #ifdef _OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for(int i=1 ; i<=o ; i++) {
      for(int j=1 ; j<=m ; j++) {
         int ind = ((i-1)*m+(j-1));
         complex c = complex(DCR[ind],DCI[ind]);
         if(i==j) c += 1.0;
         C[lbC_i+j-1][lbC_j+i-1] = c;
      }
   }        
 
   //Compute upper bound of result
   setround(1);
   beta = 0.0;
   cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, o, n, alpha, DAR, 
          LDA, DBR, LDB, beta, DCR, LDC);
   beta = 1.0;
   cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, o, n, alpha, DAIm, 
          LDA, DBI, LDB, beta, DCR, LDC);
   beta = 0.0;
   cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, o, n, alpha, DAR, 
          LDA, DBI, LDB, beta, DCI, LDC);
   beta = 1.0;
   cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, o, n, alpha, DAI, 
          LDA, DBR, LDB, beta, DCI, LDC);

   //Copy upper bound into C
   #ifdef _OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for(int i=1 ; i<=o ; i++) {
      for(int j=1 ; j<=m ; j++) {
         int ind = ((i-1)*m+(j-1));
         complex c = complex(DCR[ind],DCI[ind]);
         if(i==j) c += 1.0;         
         SetSup(C[lbC_i+j-1][lbC_j+i-1], c);
      }
   }       

   setround(rnd);
   
   delete[] DAR;
   delete[] DAI;
   delete[] DAIm;
   delete[] DBR;
   delete[] DBI;   
   delete[] DCR;    
   delete[] DCI; 
}

//Computes [C]=[I-R[A]] using BLAS routines
static inline void IminusRA(const cmatrix& R, const cimatrix& A, cimatrix& C) {
  int rnd = getround();

  C = -R*A;

  for(int i=Lb(C,1) ; i<=Ub(C,1) ; i++)
    C[i][i] += 1.0;

  setround(rnd);
}

#endif
#endif

} //namespace cxsc

#ifndef _CXSC_FAST_LSS_INLINE
#define _CXSC_FAST_LSS_INLINE
#include <fastlss.inl>
#endif

namespace cxsc {

#if defined(_CXSC_LSS_HPP) && defined(_CXSC_LSS_INCLUDE)
#undef _CXSC_LSS_INCLUDE
/*!
 *  Computes a verified solution of the real linear point system AX=B by bounding the residual of a 
 *  computed approximate solution using the Krawczyk-Operator. 
 * 
 *  This is a legacy version of the starting function. You should use the new versions which can be configured
 *  by a struct (giving a lot more options than this versions).
 * 
 *  For best performance, compile with max. compiler optimizations, activate OpenMP, and set the preprocessor variables
 *  CXSC_USE_BLAS and CXSC_USE_LAPACK (compiler options -DCXSC_USE_BLAS -DCXSC_USE_LAPACK for most compilers) and link
 *  to optimize BLAS and LAPACK libraries.
 *
 *   \param A      System matrix
 *   \param B      Right-hand side
 *   \param X      If successful this is an enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1
 *   \param Err    Error code (0=no error occured). Meaning of error codes can be parsed by calling LinSolveErrMsg(Err)
 *   \param K      Precision to use for residual calculations
 *                 If K==0 the C-XSC accumulator will be used (=maximum precision, but slower), K==1 uses normal double precision
 *                 and a splitting algorithm by Rump/Dekker to compute an improved residual, K>2 uses simulated K-fold double precision
 *                 using the DotK algorithm 
 *   \param np     Number of thread to be used for OpenMP.
 *   \param msg    If true, status messages will be put out to cout during 
 *                 computation
 *   \param lsscfg Determines, if only part 1 (LSS_ONLY_PART_ONE), only part 2 (LSS_ONLY_PART_TWO) or both parts
 *                 of the solver are to be executed. Part one is faster and works for condition numbers up to about 10^15. Part two
 *                 is much slower but can solve systems with higher condition numbers (up to about 10^30)
*/
static inline void lss( const rmatrix& A, rmatrix& B, imatrix& X, int& Err, int K, bool msg, int lsscfg, int np) {
  imatrix dummy;
  struct lssconfig cfg;
  cfg.K = K;
  cfg.msg = msg;
  cfg.lssparts = lsscfg;
  cfg.threads = np;
  SolverStart<rmatrix, rmatrix, imatrix, rmatrix, imatrix, rmatrix, rvector, dotprecision, idotprecision, ivector, interval>(A, B, X, Err, cfg, false, dummy);
}

/*!
 *  Computes a verified solution of the real linear point system Ax=b by bounding the residual of a 
 *  computed approximate solution using the Krawczyk-Operator. 
 * 
 *  This is a legacy version of the starting function. You should use the new versions which can be configured
 *  by a struct (giving a lot more options than this versions).
 * 
 *  For best performance, compile with max. compiler optimizations, activate OpenMP, and set the preprocessor variables
 *  CXSC_USE_BLAS and CXSC_USE_LAPACK (compiler options -DCXSC_USE_BLAS -DCXSC_USE_LAPACK for most compilers) and link
 *  to optimize BLAS and LAPACK libraries.
 *
 *   \param A      System matrix
 *   \param b      Right-hand side
 *   \param x      If successful this is an enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1
 *   \param Err    Error code (0=no error occured). Meaning of error codes can be parsed by calling LinSolveErrMsg(Err)
 *   \param K      Precision to use for residual calculations
 *                 If K==0 the C-XSC accumulator will be used (=maximum precision, but slower), K==1 uses normal double precision
 *                 and a splitting algorithm by Rump/Dekker to compute an improved residual, K>2 uses simulated K-fold double precision
 *                 using the DotK algorithm 
 *   \param np     Number of thread to be used for OpenMP.
 *   \param msg    If true, status messages will be put out to cout during 
 *                 computation
 *   \param lsscfg Determines, if only part 1 (LSS_ONLY_PART_ONE), only part 2 (LSS_ONLY_PART_TWO) or both parts
 *                 of the solver are to be executed. Part one is faster and works for condition numbers up to about 10^15. Part two
 *                 is much slower but can solve systems with higher condition numbers (up to about 10^30)
*/
static inline void lss( const rmatrix& A, rvector& b, ivector& x, int& Err, int K, bool msg, int lsscfg, int np) {
  rmatrix B(VecLen(b),1);
  imatrix X(VecLen(x),1);
  B[Col(1)] = b;
  lss(A,B,X,Err,K,msg,lsscfg,np);
  x = X[Col(1)];
}

/*!
 *  Computes a verified solution of the real linear point system AX=B by bounding the residual of a 
 *  computed approximate solution using the Krawczyk-Operator. 
 * 
 *  This is a legacy version of the starting function. You should use the new versions which can be configured
 *  by a struct (giving a lot more options than this versions).
 * 
 *  For best performance, compile with max. compiler optimizations, activate OpenMP, and set the preprocessor variables
 *  CXSC_USE_BLAS and CXSC_USE_LAPACK (compiler options -DCXSC_USE_BLAS -DCXSC_USE_LAPACK for most compilers) and link
 *  to optimize BLAS and LAPACK libraries.
 *
 *   \param A      System matrix
 *   \param B      Right-hand side
 *   \param X      If successful this is an enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1
 *   \param Err    Error code (0=no error occured). Meaning of error codes can be parsed by calling LinSolveErrMsg(Err)
 *   \param K      Precision to use for residual calculations
 *                 If K==0 the C-XSC accumulator will be used (=maximum precision, but slower), K==1 uses normal double precision
 *                 and a splitting algorithm by Rump/Dekker to compute an improved residual, K>2 uses simulated K-fold double precision
 *                 using the DotK algorithm 
 *   \param np     Number of thread to be used for OpenMP.
 *   \param msg    If true, status messages will be put out to cout during 
 *                 computation
 *   \param lsscfg Determines, if only part 1 (LSS_ONLY_PART_ONE), only part 2 (LSS_ONLY_PART_TWO) or both parts
 *                 of the solver are to be executed. Part one is faster and works for condition numbers up to about 10^15. Part two
 *                 is much slower but can solve systems with higher condition numbers (up to about 10^30)
*/
static inline void lss( const rmatrix& A, imatrix& B, imatrix& X, int& Err, int K, bool msg, int lsscfg, int np) {
  imatrix dummy;
  struct lssconfig cfg;
  cfg.K = K;
  cfg.msg = msg;
  cfg.lssparts = lsscfg;
  cfg.threads = np;
  SolverStart<rmatrix, imatrix, imatrix, rmatrix, imatrix, rmatrix, rvector, dotprecision, idotprecision, ivector, interval>(A, B, X, Err, cfg, false, dummy);
}

/*!
 *  Computes a verified solution of the real linear point system AX=B by bounding the residual of a 
 *  computed approximate solution using the Krawczyk-Operator. 
 * 
 *  This is a legacy version of the starting function. You should use the new versions which can be configured
 *  by a struct (giving a lot more options than this versions).
 * 
 *  For best performance, compile with max. compiler optimizations, activate OpenMP, and set the preprocessor variables
 *  CXSC_USE_BLAS and CXSC_USE_LAPACK (compiler options -DCXSC_USE_BLAS -DCXSC_USE_LAPACK for most compilers) and link
 *  to optimize BLAS and LAPACK libraries.
 *
 *   \param A      System matrix
 *   \param b      Right-hand side
 *   \param x      If successful this is an enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1
 *   \param Err    Error code (0=no error occured). Meaning of error codes can be parsed by calling LinSolveErrMsg(Err)
 *   \param K      Precision to use for residual calculations
 *                 If K==0 the C-XSC accumulator will be used (=maximum precision, but slower), K==1 uses normal double precision
 *                 and a splitting algorithm by Rump/Dekker to compute an improved residual, K>2 uses simulated K-fold double precision
 *                 using the DotK algorithm 
 *   \param np     Number of thread to be used for OpenMP.
 *   \param msg    If true, status messages will be put out to cout during 
 *                 computation
 *   \param lsscfg Determines, if only part 1 (LSS_ONLY_PART_ONE), only part 2 (LSS_ONLY_PART_TWO) or both parts
 *                 of the solver are to be executed. Part one is faster and works for condition numbers up to about 10^15. Part two
 *                 is much slower but can solve systems with higher condition numbers (up to about 10^30)
*/
static inline void lss( const rmatrix& A, ivector& b, ivector& x, int& Err, int K, bool msg, int lsscfg, int np) {
  imatrix B(VecLen(b),1);
  imatrix X(VecLen(x),1);
  B[Col(1)] = b;
  lss(A,B,X,Err,K,msg,lsscfg,np);
  x = X[Col(1)];
}


/*!
 *  Computes a verified solution of the real linear point system AX=B by bounding the residual of a 
 *  computed approximate solution using the Krawczyk-Operator. 
 * 
 *  New version of the starting function. Configuration options can be set by passing an instance of the 
 *  struct lssconfig.
 * 
 *  For best performance, compile with max. compiler optimizations, activate OpenMP, and set the preprocessor variables
 *  CXSC_USE_BLAS and CXSC_USE_LAPACK (compiler options -DCXSC_USE_BLAS -DCXSC_USE_LAPACK for most compilers) and link
 *  to optimize BLAS and LAPACK libraries.
 *
 *   \param A      System matrix
 *   \param B      Right-hand side
 *   \param X      If successful this is an enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1
 *   \param Err    Error code (0=no error occured). Meaning of error codes can be parsed by calling LinSolveErrMsg(Err)
 *   \param cfg    An instance of the struct lssconfig containing configuration information for the solver. This parameter
 *                 is optional, default settings are applied if it is not passed.
*/
static inline void lss( const rmatrix& A, rmatrix& B, imatrix& X, int& Err, struct lssconfig cfg) {
  imatrix ienc;
  SolverStart<rmatrix, rmatrix, imatrix, rmatrix, imatrix, rmatrix, rvector, dotprecision, idotprecision, ivector, interval>(A, B, X, Err, cfg, false, ienc);
}

/*!
 *  Computes a verified solution of the real linear point system Ax=b by bounding the residual of a 
 *  computed approximate solution using the Krawczyk-Operator. 
 * 
 *  New version of the starting function. Configuration options can be set by passing an instance of the 
 *  struct lssconfig.
 * 
 *  For best performance, compile with max. compiler optimizations, activate OpenMP, and set the preprocessor variables
 *  CXSC_USE_BLAS and CXSC_USE_LAPACK (compiler options -DCXSC_USE_BLAS -DCXSC_USE_LAPACK for most compilers) and link
 *  to optimize BLAS and LAPACK libraries.
 *
 *   \param A      System matrix
 *   \param b      Right-hand side
 *   \param x      If successful this is an enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1
 *   \param Err    Error code (0=no error occured). Meaning of error codes can be parsed by calling LinSolveErrMsg(Err)
 *   \param cfg    An instance of the struct lssconfig containing configuration information for the solver. This parameter
 *                 is optional, default settings are applied if it is not passed.
*/
static inline void lss( const rmatrix& A, rvector& b, ivector& x, int& Err, struct lssconfig cfg) {
  rmatrix B(VecLen(b),1);
  imatrix X(VecLen(x),1);
  B[Col(1)] = b;
  lss(A,B,X,Err,cfg);
  x = X[Col(1)];
}

/*!
 *  Computes a verified solution of the real linear point system AX=B by bounding the residual of a 
 *  computed approximate solution using the Krawczyk-Operator. 
 * 
 *  New version of the starting function. Configuration options can be set by passing an instance of the 
 *  struct lssconfig.
 * 
 *  For best performance, compile with max. compiler optimizations, activate OpenMP, and set the preprocessor variables
 *  CXSC_USE_BLAS and CXSC_USE_LAPACK (compiler options -DCXSC_USE_BLAS -DCXSC_USE_LAPACK for most compilers) and link
 *  to optimize BLAS and LAPACK libraries.
 *
 *   \param A      System matrix
 *   \param B      Right-hand side
 *   \param X      If successful this is an enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1
 *   \param Err    Error code (0=no error occured). Meaning of error codes can be parsed by calling LinSolveErrMsg(Err)
 *   \param cfg    An instance of the struct lssconfig containing configuration information for the solver. This parameter
 *                 is optional, default settings are applied if it is not passed.
*/
static inline void lss( const rmatrix& A, imatrix& B, imatrix& X, int& Err, struct lssconfig cfg) {
  imatrix ienc;
  SolverStart<rmatrix, imatrix, imatrix, rmatrix, imatrix, rmatrix, rvector, dotprecision, idotprecision, ivector, interval>(A, B, X, Err, cfg, false, ienc);
}

/*!
 *  Computes a verified solution of the real linear point system Ax=b by bounding the residual of a 
 *  computed approximate solution using the Krawczyk-Operator. 
 * 
 *  New version of the starting function. Configuration options can be set by passing an instance of the 
 *  struct lssconfig.
 * 
 *  For best performance, compile with max. compiler optimizations, activate OpenMP, and set the preprocessor variables
 *  CXSC_USE_BLAS and CXSC_USE_LAPACK (compiler options -DCXSC_USE_BLAS -DCXSC_USE_LAPACK for most compilers) and link
 *  to optimize BLAS and LAPACK libraries.
 *
 *   \param A      System matrix
 *   \param b      Right-hand side
 *   \param x      If successful this is an enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1
 *   \param Err    Error code (0=no error occured). Meaning of error codes can be parsed by calling LinSolveErrMsg(Err)
 *   \param cfg    An instance of the struct lssconfig containing configuration information for the solver. This parameter
 *                 is optional, default settings are applied if it is not passed.
*/

static inline void lss( const rmatrix& A, ivector& b, ivector& x, int& Err, struct lssconfig cfg) {
  imatrix B(VecLen(b),1);
  imatrix X(VecLen(x),1);
  B[Col(1)] = b;
  lss(A,B,X,Err,cfg);
  x = X[Col(1)];
}


/*!
 *  Computes a verified solution of the real linear point system AX=B by bounding the residual of a 
 *  computed approximate solution using the Krawczyk-Operator. 
 * 
 *  New version of the starting function. Configuration options can be set by passing an instance of the 
 *  struct lssconfig.
 * 
 *  For best performance, compile with max. compiler optimizations, activate OpenMP, and set the preprocessor variables
 *  CXSC_USE_BLAS and CXSC_USE_LAPACK (compiler options -DCXSC_USE_BLAS -DCXSC_USE_LAPACK for most compilers) and link
 *  to optimize BLAS and LAPACK libraries.
 *
 *   \param A      System matrix
 *   \param B      Right-hand side
 *   \param X      If successful this is an enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1
 *   \param Y      If successful this is an inner enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1. Elements for which no inner enclosure could be found are set to SignalingNaN.
 *   \param Err    Error code (0=no error occured). Meaning of error codes can be parsed by calling LinSolveErrMsg(Err)
 *   \param cfg    An instance of the struct lssconfig containing configuration information for the solver. This parameter
 *                 is optional, default settings are applied if it is not passed.
*/
static inline void lss( const rmatrix& A, rmatrix& B, imatrix& X, imatrix& Y, int& Err, struct lssconfig cfg) {
  SolverStart<rmatrix, rmatrix, imatrix, rmatrix, imatrix, rmatrix, rvector, dotprecision, idotprecision, ivector, interval>(A, B, X, Err, cfg, true, Y);
}

/*!
 *  Computes a verified solution of the real linear point system Ax=b by bounding the residual of a 
 *  computed approximate solution using the Krawczyk-Operator. 
 * 
 *  New version of the starting function. Configuration options can be set by passing an instance of the 
 *  struct lssconfig.
 * 
 *  For best performance, compile with max. compiler optimizations, activate OpenMP, and set the preprocessor variables
 *  CXSC_USE_BLAS and CXSC_USE_LAPACK (compiler options -DCXSC_USE_BLAS -DCXSC_USE_LAPACK for most compilers) and link
 *  to optimize BLAS and LAPACK libraries.
 *
 *   \param A      System matrix
 *   \param b      Right-hand side
 *   \param x      If successful this is an enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1
 *   \param y      If successful this is an inner enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1. Elements for which no inner enclosure could be found are set to SignalingNaN.
 *   \param Err    Error code (0=no error occured). Meaning of error codes can be parsed by calling LinSolveErrMsg(Err)
 *   \param cfg    An instance of the struct lssconfig containing configuration information for the solver. This parameter
 *                 is optional, default settings are applied if it is not passed.
*/
static inline void lss( const rmatrix& A, rvector& b, ivector& x, ivector& y, int& Err, struct lssconfig cfg) {
  rmatrix B(VecLen(b),1);
  imatrix X(VecLen(x),1);
  imatrix Y(VecLen(y),1);
  B[Col(1)] = b;
  lss(A,B,X,Y,Err,cfg);
  x = X[Col(1)];
  y = Y[Col(1)];
}

/*!
 *  Computes a verified solution of the real linear point system AX=B by bounding the residual of a 
 *  computed approximate solution using the Krawczyk-Operator. 
 * 
 *  New version of the starting function. Configuration options can be set by passing an instance of the 
 *  struct lssconfig.
 * 
 *  For best performance, compile with max. compiler optimizations, activate OpenMP, and set the preprocessor variables
 *  CXSC_USE_BLAS and CXSC_USE_LAPACK (compiler options -DCXSC_USE_BLAS -DCXSC_USE_LAPACK for most compilers) and link
 *  to optimize BLAS and LAPACK libraries.
 *
 *   \param A      System matrix
 *   \param B      Right-hand side
 *   \param X      If successful this is an enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1
 *   \param Y      If successful this is an inner enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1. Elements for which no inner enclosure could be found are set to SignalingNaN.
 *   \param Err    Error code (0=no error occured). Meaning of error codes can be parsed by calling LinSolveErrMsg(Err)
 *   \param cfg    An instance of the struct lssconfig containing configuration information for the solver. This parameter
 *                 is optional, default settings are applied if it is not passed.
*/
static inline void lss( const rmatrix& A, imatrix& B, imatrix& X, imatrix& Y, int& Err, struct lssconfig cfg) {
  SolverStart<rmatrix, imatrix, imatrix, rmatrix, imatrix, rmatrix, rvector, dotprecision, idotprecision, ivector, interval>(A, B, X, Err, cfg, true, Y);
}

/*!
 *  Computes a verified solution of the real linear point system Ax=b by bounding the residual of a 
 *  computed approximate solution using the Krawczyk-Operator. 
 * 
 *  New version of the starting function. Configuration options can be set by passing an instance of the 
 *  struct lssconfig.
 * 
 *  For best performance, compile with max. compiler optimizations, activate OpenMP, and set the preprocessor variables
 *  CXSC_USE_BLAS and CXSC_USE_LAPACK (compiler options -DCXSC_USE_BLAS -DCXSC_USE_LAPACK for most compilers) and link
 *  to optimize BLAS and LAPACK libraries.
 *
 *   \param A      System matrix
 *   \param b      Right-hand side
 *   \param x      If successful this is an enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1
 *   \param y      If successful this is an inner enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1. Elements for which no inner enclosure could be found are set to SignalingNaN.
 *   \param Err    Error code (0=no error occured). Meaning of error codes can be parsed by calling LinSolveErrMsg(Err)
 *   \param cfg    An instance of the struct lssconfig containing configuration information for the solver. This parameter
 *                 is optional, default settings are applied if it is not passed.
*/
static inline void lss( const rmatrix& A, ivector& b, ivector& x, ivector& y, int& Err, struct lssconfig cfg) {
  imatrix B(VecLen(b),1);
  imatrix X(VecLen(x),1);
  imatrix Y(VecLen(y),1);
  B[Col(1)] = b;
  lss(A,B,X,Y,Err,cfg);
  x = X[Col(1)];
  y = Y[Col(1)];
}

#endif

//----------------------------------------------------------------------------                      

#if defined(_CXSC_ILSS_HPP) && defined(_CXSC_ILSS_INCLUDE)
#undef _CXSC_ILSS_INCLUDE
/*!
 *  Computes a verified solution of the real linear interval system AX=B by bounding the residual of a 
 *  computed approximate solution using the Krawczyk-Operator. 
 * 
 *  This is a legacy version of the starting function. You should use the new versions which can be configured
 *  by a struct (giving a lot more options than this versions).
 * 
 *  For best performance, compile with max. compiler optimizations, activate OpenMP, and set the preprocessor variables
 *  CXSC_USE_BLAS and CXSC_USE_LAPACK (compiler options -DCXSC_USE_BLAS -DCXSC_USE_LAPACK for most compilers) and link
 *  to optimize BLAS and LAPACK libraries.
 *
 *   \param A      System matrix
 *   \param B      Right-hand side
 *   \param X      If successful this is an enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1
 *   \param Err    Error code (0=no error occured). Meaning of error codes can be parsed by calling LinSolveErrMsg(Err)
 *   \param K      Precision to use for residual calculations
 *                 If K==0 the C-XSC accumulator will be used (=maximum precision, but slower), K==1 uses normal double precision
 *                 and a splitting algorithm by Rump/Dekker to compute an improved residual, K>2 uses simulated K-fold double precision
 *                 using the DotK algorithm 
 *   \param np     Number of thread to be used for OpenMP.
 *   \param msg    If true, status messages will be put out to cout during 
 *                 computation
 *   \param lsscfg Determines, if only part 1 (LSS_ONLY_PART_ONE), only part 2 (LSS_ONLY_PART_TWO) or both parts
 *                 of the solver are to be executed. Part one is faster and works for condition numbers up to about 10^15. Part two
 *                 is much slower but can solve systems with higher condition numbers (up to about 10^30)
*/
static inline void ilss( const imatrix& A, imatrix& B, imatrix& X, int& Err, int K, bool msg, int lsscfg, int np) {
  imatrix dummy;
  struct lssconfig cfg;
  cfg.K = K;
  cfg.msg = msg;
  cfg.lssparts = lsscfg;
  cfg.threads = np;  
  SolverStart<imatrix, imatrix, imatrix, rmatrix, imatrix, rmatrix, rvector, dotprecision, idotprecision, ivector, interval>(A, B, X, Err, cfg, false, dummy);
}


/*!
 *  Computes a verified solution of the real linear interval system Ax=b by bounding the residual of a 
 *  computed approximate solution using the Krawczyk-Operator. 
 * 
 *  This is a legacy version of the starting function. You should use the new versions which can be configured
 *  by a struct (giving a lot more options than this versions).
 * 
 *  For best performance, compile with max. compiler optimizations, activate OpenMP, and set the preprocessor variables
 *  CXSC_USE_BLAS and CXSC_USE_LAPACK (compiler options -DCXSC_USE_BLAS -DCXSC_USE_LAPACK for most compilers) and link
 *  to optimize BLAS and LAPACK libraries.
 *
 *   \param A      System matrix
 *   \param b      Right-hand side
 *   \param x      If successful this is an enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1
 *   \param Err    Error code (0=no error occured). Meaning of error codes can be parsed by calling LinSolveErrMsg(Err)
 *   \param K      Precision to use for residual calculations
 *                 If K==0 the C-XSC accumulator will be used (=maximum precision, but slower), K==1 uses normal double precision
 *                 and a splitting algorithm by Rump/Dekker to compute an improved residual, K>2 uses simulated K-fold double precision
 *                 using the DotK algorithm 
 *   \param np     Number of thread to be used for OpenMP.
 *   \param msg    If true, status messages will be put out to cout during 
 *                 computation
 *   \param lsscfg Determines, if only part 1 (LSS_ONLY_PART_ONE), only part 2 (LSS_ONLY_PART_TWO) or both parts
 *                 of the solver are to be executed. Part one is faster and works for condition numbers up to about 10^15. Part two
 *                 is much slower but can solve systems with higher condition numbers (up to about 10^30)
*/
static inline void ilss( const imatrix& A, ivector& b, ivector& x, int& Err, int K, bool msg, int lsscfg, int np) {
  imatrix B(VecLen(b),1);
  imatrix X(VecLen(x),1);
  B[Col(1)] = b;
  ilss(A,B,X,Err,K,msg,lsscfg,np);
  x = X[Col(1)];
}

/*!
 *  Computes a verified solution of the real linear interval system AX=B by bounding the residual of a 
 *  computed approximate solution using the Krawczyk-Operator. 
 * 
 *  New version of the starting function. Configuration options can be set by passing an instance of the 
 *  struct lssconfig.
 * 
 *  For best performance, compile with max. compiler optimizations, activate OpenMP, and set the preprocessor variables
 *  CXSC_USE_BLAS and CXSC_USE_LAPACK (compiler options -DCXSC_USE_BLAS -DCXSC_USE_LAPACK for most compilers) and link
 *  to optimize BLAS and LAPACK libraries.
 *
 *
 *   \param A      System matrix
 *   \param B      Right-hand side
 *   \param X      If successful this is an enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1
 *   \param Err    Error code (0=no error occured). Meaning of error codes can be parsed by calling LinSolveErrMsg(Err)
 *   \param cfg    An instance of the struct lssconfig containing configuration information for the solver. This parameter
 *                 is optional, default settings are applied if it is not passed.*/
static inline void ilss( const imatrix& A, imatrix& B, imatrix& X, int& Err, struct lssconfig cfg) {
  imatrix ienc;  
  SolverStart<imatrix, imatrix, imatrix, rmatrix, imatrix, rmatrix, rvector, dotprecision, idotprecision, ivector, interval>(A, B, X, Err, cfg, false, ienc);
}

/*!
 *  Computes a verified solution of the real linear interval system AX=B by bounding the residual of a 
 *  computed approximate solution using the Krawczyk-Operator. 
 * 
 *  New version of the starting function. Configuration options can be set by passing an instance of the 
 *  struct lssconfig.
 * 
 *  For best performance, compile with max. compiler optimizations, activate OpenMP, and set the preprocessor variables
 *  CXSC_USE_BLAS and CXSC_USE_LAPACK (compiler options -DCXSC_USE_BLAS -DCXSC_USE_LAPACK for most compilers) and link
 *  to optimize BLAS and LAPACK libraries.
 *
 *   \param A      System matrix
 *   \param b      Right-hand side
 *   \param x      If successful this is an enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1
 *   \param Err    Error code (0=no error occured). Meaning of error codes can be parsed by calling LinSolveErrMsg(Err)
 *   \param cfg    An instance of the struct lssconfig containing configuration information for the solver. This parameter
 *                 is optional, default settings are applied if it is not passed.
*/
static inline void ilss( const imatrix& A, ivector& b, ivector& x, int& Err, struct lssconfig cfg) {
  imatrix B(VecLen(b),1);
  imatrix X(VecLen(x),1);
  B[Col(1)] = b;
  ilss(A,B,X,Err,cfg);
  x = X[Col(1)];
}

/*!
 *  Computes a verified solution of the real linear interval system AX=B by bounding the residual of a 
 *  computed approximate solution using the Krawczyk-Operator. 
 * 
 *  New version of the starting function. Configuration options can be set by passing an instance of the 
 *  struct lssconfig.
 * 
 *  For best performance, compile with max. compiler optimizations, activate OpenMP, and set the preprocessor variables
 *  CXSC_USE_BLAS and CXSC_USE_LAPACK (compiler options -DCXSC_USE_BLAS -DCXSC_USE_LAPACK for most compilers) and link
 *  to optimize BLAS and LAPACK libraries.
 *
 *   \param A      System matrix
 *   \param B      Right-hand side
 *   \param X      If successful this is an enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1
 *   \param Y      If successful this is an inner enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1. Elements for which no inner enclosure could be found are set to SignalingNaN.
 *   \param Err    Error code (0=no error occured). Meaning of error codes can be parsed by calling LinSolveErrMsg(Err)
 *   \param cfg    An instance of the struct lssconfig containing configuration information for the solver. This parameter
 *                 is optional, default settings are applied if it is not passed.
*/
static inline void ilss( const imatrix& A, imatrix& B, imatrix& X, imatrix& Y, int& Err, struct lssconfig cfg) {
  SolverStart<imatrix, imatrix, imatrix, rmatrix, imatrix, rmatrix, rvector, dotprecision, idotprecision, ivector, interval>(A, B, X, Err, cfg, true, Y);
}

/*!
 *  Computes a verified solution of the real linear interval system Ax=b by bounding the residual of a 
 *  computed approximate solution using the Krawczyk-Operator. 
 * 
 *  New version of the starting function. Configuration options can be set by passing an instance of the 
 *  struct lssconfig.
 * 
 *  For best performance, compile with max. compiler optimizations, activate OpenMP, and set the preprocessor variables
 *  CXSC_USE_BLAS and CXSC_USE_LAPACK (compiler options -DCXSC_USE_BLAS -DCXSC_USE_LAPACK for most compilers) and link
 *  to optimize BLAS and LAPACK libraries.
 *
 *   \param A      System matrix
 *   \param b      Right-hand side
 *   \param x      If successful this is an enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1
 *   \param y      If successful this is an inner enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1. Elements for which no inner enclosure could be found are set to SignalingNaN.
 *   \param Err    Error code (0=no error occured). Meaning of error codes can be parsed by calling LinSolveErrMsg(Err)
 *   \param cfg    An instance of the struct lssconfig containing configuration information for the solver. This parameter
 *                 is optional, default settings are applied if it is not passed.
*/
static inline void ilss( const imatrix& A, ivector& b, ivector& x, ivector& y, int& Err, struct lssconfig cfg) {
  imatrix B(VecLen(b),1);
  imatrix X(VecLen(x),1);
  imatrix Y(VecLen(y),1);
  B[Col(1)] = b;
  ilss(A,B,X,Y,Err,cfg);
  x = X[Col(1)];
  y = Y[Col(1)];
}
#endif

//----------------------------------------------------------------------------                      

#if defined(_CXSC_CLSS_HPP) && defined(_CXSC_CLSS_INCLUDE)
#undef _CXSC_CLSS_INCLUDE
/*!
 *  Computes a verified solution of the complex linear point system AX=B by bounding the residual of a 
 *  computed approximate solution using the Krawczyk-Operator. 
 * 
 *  This is a legacy version of the starting function. You should use the new versions which can be configured
 *  by a struct (giving a lot more options than this versions).
 * 
 *  For best performance, compile with max. compiler optimizations, activate OpenMP, and set the preprocessor variables
 *  CXSC_USE_BLAS and CXSC_USE_LAPACK (compiler options -DCXSC_USE_BLAS -DCXSC_USE_LAPACK for most compilers) and link
 *  to optimize BLAS and LAPACK libraries.
 *
 *   \param A      System matrix
 *   \param B      Right-hand side
 *   \param X      If successful this is an enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1
 *   \param Err    Error code (0=no error occured). Meaning of error codes can be parsed by calling LinSolveErrMsg(Err)
 *   \param K      Precision to use for residual calculations
 *                 If K==0 the C-XSC accumulator will be used (=maximum precision, but slower), K==1 uses normal double precision
 *                 and a splitting algorithm by Rump/Dekker to compute an improved residual, K>2 uses simulated K-fold double precision
 *                 using the DotK algorithm 
 *   \param np     Number of thread to be used for OpenMP.
 *   \param msg    If true, status messages will be put out to cout during 
 *                 computation
 *   \param lsscfg Determines, if only part 1 (LSS_ONLY_PART_ONE), only part 2 (LSS_ONLY_PART_TWO) or both parts
 *                 of the solver are to be executed. Part one is faster and works for condition numbers up to about 10^15. Part two
 *                 is much slower but can solve systems with higher condition numbers (up to about 10^30)
*/
static inline void clss( const cmatrix& A, cmatrix& B, cimatrix& X, int& Err, int K, bool msg, int lsscfg, int np) {
  cimatrix dummy;
  struct lssconfig cfg;
  cfg.K = K;
  cfg.msg = msg;
  cfg.lssparts = lsscfg;
  cfg.threads = np;    
  SolverStart<cmatrix, cmatrix, cimatrix, cmatrix, cimatrix, cmatrix, cvector, cdotprecision, cidotprecision, civector, cinterval>(A, B, X, Err, cfg, false, dummy);
}

/*!
 *  Computes a verified solution of the complex linear point system Ax=b by bounding the residual of a 
 *  computed approximate solution using the Krawczyk-Operator. 
 * 
 *  This is a legacy version of the starting function. You should use the new versions which can be configured
 *  by a struct (giving a lot more options than this versions).
 * 
 *  For best performance, compile with max. compiler optimizations, activate OpenMP, and set the preprocessor variables
 *  CXSC_USE_BLAS and CXSC_USE_LAPACK (compiler options -DCXSC_USE_BLAS -DCXSC_USE_LAPACK for most compilers) and link
 *  to optimize BLAS and LAPACK libraries.
 *
 *   \param A      System matrix
 *   \param b      Right-hand side
 *   \param x      If successful this is an enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1
 *   \param Err    Error code (0=no error occured). Meaning of error codes can be parsed by calling LinSolveErrMsg(Err)
 *   \param K      Precision to use for residual calculations
 *                 If K==0 the C-XSC accumulator will be used (=maximum precision, but slower), K==1 uses normal double precision
 *                 and a splitting algorithm by Rump/Dekker to compute an improved residual, K>2 uses simulated K-fold double precision
 *                 using the DotK algorithm 
 *   \param np     Number of thread to be used for OpenMP.
 *   \param msg    If true, status messages will be put out to cout during 
 *                 computation
 *   \param lsscfg Determines, if only part 1 (LSS_ONLY_PART_ONE), only part 2 (LSS_ONLY_PART_TWO) or both parts
 *                 of the solver are to be executed. Part one is faster and works for condition numbers up to about 10^15. Part two
 *                 is much slower but can solve systems with higher condition numbers (up to about 10^30)
*/
static inline void clss( const cmatrix& A, cvector& b, civector& x, int& Err, int K, bool msg, int lsscfg, int np) {
  cmatrix B(VecLen(b),1);
  cimatrix X(VecLen(x),1);
  B[Col(1)] = b;
  clss(A,B,X,Err,K,msg,lsscfg,np);
  x = X[Col(1)];
}

/*!
 *  Computes a verified solution of the complex linear point system AX=B by bounding the residual of a 
 *  computed approximate solution using the Krawczyk-Operator. 
 * 
 *  This is a legacy version of the starting function. You should use the new versions which can be configured
 *  by a struct (giving a lot more options than this versions).
 * 
 *  For best performance, compile with max. compiler optimizations, activate OpenMP, and set the preprocessor variables
 *  CXSC_USE_BLAS and CXSC_USE_LAPACK (compiler options -DCXSC_USE_BLAS -DCXSC_USE_LAPACK for most compilers) and link
 *  to optimize BLAS and LAPACK libraries.
 *
 *   \param A      System matrix
 *   \param B      Right-hand side
 *   \param X      If successful this is an enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1
 *   \param Err    Error code (0=no error occured). Meaning of error codes can be parsed by calling LinSolveErrMsg(Err)
 *   \param K      Precision to use for residual calculations
 *                 If K==0 the C-XSC accumulator will be used (=maximum precision, but slower), K==1 uses normal double precision
 *                 and a splitting algorithm by Rump/Dekker to compute an improved residual, K>2 uses simulated K-fold double precision
 *                 using the DotK algorithm 
 *   \param np     Number of thread to be used for OpenMP.
 *   \param msg    If true, status messages will be put out to cout during 
 *                 computation
 *   \param lsscfg Determines, if only part 1 (LSS_ONLY_PART_ONE), only part 2 (LSS_ONLY_PART_TWO) or both parts
 *                 of the solver are to be executed. Part one is faster and works for condition numbers up to about 10^15. Part two
 *                 is much slower but can solve systems with higher condition numbers (up to about 10^30)
*/
static inline void clss( const cmatrix& A, cimatrix& B, cimatrix& X, int& Err, int K, bool msg, int lsscfg, int np) {
  cimatrix dummy;
  struct lssconfig cfg;
  cfg.K = K;
  cfg.msg = msg;
  cfg.lssparts = lsscfg;
  cfg.threads = np;    
  SolverStart<cmatrix, cimatrix, cimatrix, cmatrix, cimatrix, cmatrix, cvector, cdotprecision, cidotprecision, civector, cinterval>(A, B, X, Err, cfg, false, dummy);
}

/*!
 *  Computes a verified solution of the complex linear point system Ax=b by bounding the residual of a 
 *  computed approximate solution using the Krawczyk-Operator. 
 * 
 *  This is a legacy version of the starting function. You should use the new versions which can be configured
 *  by a struct (giving a lot more options than this versions).
 * 
 *  For best performance, compile with max. compiler optimizations, activate OpenMP, and set the preprocessor variables
 *  CXSC_USE_BLAS and CXSC_USE_LAPACK (compiler options -DCXSC_USE_BLAS -DCXSC_USE_LAPACK for most compilers) and link
 *  to optimize BLAS and LAPACK libraries.
 *
 *   \param A      System matrix
 *   \param b      Right-hand side
 *   \param x      If successful this is an enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1
 *   \param Err    Error code (0=no error occured). Meaning of error codes can be parsed by calling LinSolveErrMsg(Err)
 *   \param K      Precision to use for residual calculations
 *                 If K==0 the C-XSC accumulator will be used (=maximum precision, but slower), K==1 uses normal double precision
 *                 and a splitting algorithm by Rump/Dekker to compute an improved residual, K>2 uses simulated K-fold double precision
 *                 using the DotK algorithm 
 *   \param np     Number of thread to be used for OpenMP.
 *   \param msg    If true, status messages will be put out to cout during 
 *                 computation
 *   \param lsscfg Determines, if only part 1 (LSS_ONLY_PART_ONE), only part 2 (LSS_ONLY_PART_TWO) or both parts
 *                 of the solver are to be executed. Part one is faster and works for condition numbers up to about 10^15. Part two
 *                 is much slower but can solve systems with higher condition numbers (up to about 10^30)
*/
static inline void clss( const cmatrix& A, civector& b, civector& x, int& Err, int K, bool msg, int lsscfg, int np) {
  cimatrix B(VecLen(b),1);
  cimatrix X(VecLen(x),1);
  B[Col(1)] = b;
  clss(A,B,X,Err,K,msg,lsscfg,np);
  x = X[Col(1)];
}


/*!
 *  Computes a verified solution of the complex linear point system AX=B by bounding the residual of a 
 *  computed approximate solution using the Krawczyk-Operator. 
 * 
 *  New version of the starting function. Configuration options can be set by passing an instance of the 
 *  struct lssconfig.
 * 
 *  For best performance, compile with max. compiler optimizations, activate OpenMP, and set the preprocessor variables
 *  CXSC_USE_BLAS and CXSC_USE_LAPACK (compiler options -DCXSC_USE_BLAS -DCXSC_USE_LAPACK for most compilers) and link
 *  to optimize BLAS and LAPACK libraries.
 *
 *   \param A      System matrix
 *   \param b      Right-hand side
 *   \param x      If successful this is an enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1
 *   \param Err    Error code (0=no error occured). Meaning of error codes can be parsed by calling LinSolveErrMsg(Err)
 *   \param cfg    An instance of the struct lssconfig containing configuration information for the solver. This parameter
 *                 is optional, default settings are applied if it is not passed.
*/
static inline void clss( const cmatrix& A, cmatrix& B, cimatrix& X, int& Err, struct lssconfig cfg) {
  cimatrix ienc;
  SolverStart<cmatrix, cmatrix, cimatrix, cmatrix, cimatrix, cmatrix, cvector, cdotprecision, cidotprecision, civector, cinterval>(A, B, X, Err, cfg, false, ienc);
}

/*!
 *  Computes a verified solution of the complex linear point system Ax=b by bounding the residual of a 
 *  computed approximate solution using the Krawczyk-Operator. 
 * 
 *  New version of the starting function. Configuration options can be set by passing an instance of the 
 *  struct lssconfig.
 * 
 *  For best performance, compile with max. compiler optimizations, activate OpenMP, and set the preprocessor variables
 *  CXSC_USE_BLAS and CXSC_USE_LAPACK (compiler options -DCXSC_USE_BLAS -DCXSC_USE_LAPACK for most compilers) and link
 *  to optimize BLAS and LAPACK libraries.
 *
 *   \param A      System matrix
 *   \param b      Right-hand side
 *   \param x      If successful this is an enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1
 *   \param Err    Error code (0=no error occured). Meaning of error codes can be parsed by calling LinSolveErrMsg(Err)
 *   \param cfg    An instance of the struct lssconfig containing configuration information for the solver. This parameter
 *                 is optional, default settings are applied if it is not passed.
*/
static inline void clss( const cmatrix& A, cvector& b, civector& x, int& Err, struct lssconfig cfg) {
  cmatrix B(VecLen(b),1);
  cimatrix X(VecLen(x),1);
  B[Col(1)] = b;
  clss(A,B,X,Err,cfg);
  x = X[Col(1)];
}

/*!
 *  Computes a verified solution of the complex linear point system AX=B by bounding the residual of a 
 *  computed approximate solution using the Krawczyk-Operator. 
 * 
 *  New version of the starting function. Configuration options can be set by passing an instance of the 
 *  struct lssconfig.
 * 
 *  For best performance, compile with max. compiler optimizations, activate OpenMP, and set the preprocessor variables
 *  CXSC_USE_BLAS and CXSC_USE_LAPACK (compiler options -DCXSC_USE_BLAS -DCXSC_USE_LAPACK for most compilers) and link
 *  to optimize BLAS and LAPACK libraries.
 *
 *   \param A      System matrix
 *   \param b      Right-hand side
 *   \param x      If successful this is an enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1
 *   \param Err    Error code (0=no error occured). Meaning of error codes can be parsed by calling LinSolveErrMsg(Err)
 *   \param cfg    An instance of the struct lssconfig containing configuration information for the solver. This parameter
 *                 is optional, default settings are applied if it is not passed.
*/
static inline void clss( const cmatrix& A, cimatrix& B, cimatrix& X, int& Err, struct lssconfig cfg) {
  cimatrix ienc;
  SolverStart<cmatrix, cimatrix, cimatrix, cmatrix, cimatrix, cmatrix, cvector, cdotprecision, cidotprecision, civector, cinterval>(A, B, X, Err, cfg, false, ienc);
}

/*!
 *  Computes a verified solution of the complex linear point system Ax=b by bounding the residual of a 
 *  computed approximate solution using the Krawczyk-Operator. 
 * 
 *  New version of the starting function. Configuration options can be set by passing an instance of the 
 *  struct lssconfig.
 * 
 *  For best performance, compile with max. compiler optimizations, activate OpenMP, and set the preprocessor variables
 *  CXSC_USE_BLAS and CXSC_USE_LAPACK (compiler options -DCXSC_USE_BLAS -DCXSC_USE_LAPACK for most compilers) and link
 *  to optimize BLAS and LAPACK libraries.
 *
 *   \param A      System matrix
 *   \param b      Right-hand side
 *   \param x      If successful this is an enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1
 *   \param Err    Error code (0=no error occured). Meaning of error codes can be parsed by calling LinSolveErrMsg(Err)
 *   \param cfg    An instance of the struct lssconfig containing configuration information for the solver. This parameter
 *                 is optional, default settings are applied if it is not passed.
*/
static inline void clss( const cmatrix& A, civector& b, civector& x, int& Err, struct lssconfig cfg) {
  cimatrix B(VecLen(b),1);
  cimatrix X(VecLen(x),1);
  B[Col(1)] = b;
  clss(A,B,X,Err,cfg);
  x = X[Col(1)];
}


/*!
 *  Computes a verified solution of the complex linear point system AX=B by bounding the residual of a 
 *  computed approximate solution using the Krawczyk-Operator. 
 * 
 *  New version of the starting function. Configuration options can be set by passing an instance of the 
 *  struct lssconfig.
 * 
 *  For best performance, compile with max. compiler optimizations, activate OpenMP, and set the preprocessor variables
 *  CXSC_USE_BLAS and CXSC_USE_LAPACK (compiler options -DCXSC_USE_BLAS -DCXSC_USE_LAPACK for most compilers) and link
 *  to optimize BLAS and LAPACK libraries.
 *
 *   \param A      System matrix
 *   \param B      Right-hand side
 *   \param X      If successful this is an enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1
 *   \param Y      If successful this is an inner enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1. Elements for which no inner enclosure could be found are set to SignalingNaN.* 
 *   \param Err    Error code (0=no error occured). Meaning of error codes can be parsed by calling LinSolveErrMsg(Err)
 *   \param cfg    An instance of the struct lssconfig containing configuration information for the solver. This parameter
 *                 is optional, default settings are applied if it is not passed.
*/
static inline void clss( const cmatrix& A, cmatrix& B, cimatrix& X, cimatrix& Y, int& Err, struct lssconfig cfg) {
  SolverStart<cmatrix, cmatrix, cimatrix, cmatrix, cimatrix, cmatrix, cvector, cdotprecision, cidotprecision, civector, cinterval>(A, B, X, Err, cfg, true, Y);
}

/*!
 *  Computes a verified solution of the complex linear point system Ax=b by bounding the residual of a 
 *  computed approximate solution using the Krawczyk-Operator. 
 * 
 *  New version of the starting function. Configuration options can be set by passing an instance of the 
 *  struct lssconfig.
 * 
 *  For best performance, compile with max. compiler optimizations, activate OpenMP, and set the preprocessor variables
 *  CXSC_USE_BLAS and CXSC_USE_LAPACK (compiler options -DCXSC_USE_BLAS -DCXSC_USE_LAPACK for most compilers) and link
 *  to optimize BLAS and LAPACK libraries.
 *
 *   \param A      System matrix
 *   \param b      Right-hand side
 *   \param x      If successful this is an enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1
 *   \param y      If successful this is an inner enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1. Elements for which no inner enclosure could be found are set to SignalingNaN. 
 *   \param Err    Error code (0=no error occured). Meaning of error codes can be parsed by calling LinSolveErrMsg(Err)
 *   \param cfg    An instance of the struct lssconfig containing configuration information for the solver. This parameter
 *                 is optional, default settings are applied if it is not passed.
*/
static inline void clss( const cmatrix& A, cvector& b, civector& x, civector& y, int& Err, struct lssconfig cfg) {
  cmatrix B(VecLen(b),1);
  cimatrix X(VecLen(x),1);
  cimatrix Y(VecLen(y),1);  
  B[Col(1)] = b;
  clss(A,B,X,Y,Err,cfg);
  x = X[Col(1)];
  y = Y[Col(1)];
}

/*!
 *  Computes a verified solution of the complex linear point system AX=B by bounding the residual of a 
 *  computed approximate solution using the Krawczyk-Operator. 
 * 
 *  New version of the starting function. Configuration options can be set by passing an instance of the 
 *  struct lssconfig.
 * 
 *  For best performance, compile with max. compiler optimizations, activate OpenMP, and set the preprocessor variables
 *  CXSC_USE_BLAS and CXSC_USE_LAPACK (compiler options -DCXSC_USE_BLAS -DCXSC_USE_LAPACK for most compilers) and link
 *  to optimize BLAS and LAPACK libraries.
 *
 *   \param A      System matrix
 *   \param B      Right-hand side
 *   \param X      If successful this is an enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1
 *   \param Y      If successful this is an inner enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1. Elements for which no inner enclosure could be found are set to SignalingNaN. 
 *   \param Err    Error code (0=no error occured). Meaning of error codes can be parsed by calling LinSolveErrMsg(Err)
 *   \param cfg    An instance of the struct lssconfig containing configuration information for the solver. This parameter
 *                 is optional, default settings are applied if it is not passed.
*/
static inline void clss( const cmatrix& A, cimatrix& B, cimatrix& X, cimatrix& Y, int& Err, struct lssconfig cfg) {
  SolverStart<cmatrix, cimatrix, cimatrix, cmatrix, cimatrix, cmatrix, cvector, cdotprecision, cidotprecision, civector, cinterval>(A, B, X, Err, cfg, true, Y);
}

/*!
 *  Computes a verified solution of the complex linear point system Ax=b by bounding the residual of a 
 *  computed approximate solution using the Krawczyk-Operator. 
 * 
 *  New version of the starting function. Configuration options can be set by passing an instance of the 
 *  struct lssconfig.
 * 
 *  For best performance, compile with max. compiler optimizations, activate OpenMP, and set the preprocessor variables
 *  CXSC_USE_BLAS and CXSC_USE_LAPACK (compiler options -DCXSC_USE_BLAS -DCXSC_USE_LAPACK for most compilers) and link
 *  to optimize BLAS and LAPACK libraries.
 *
 *   \param A      System matrix
 *   \param b      Right-hand side
 *   \param x      If successful this is an enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1
 *   \param y      If successful this is an inner enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1. Elements for which no inner enclosure could be found are set to SignalingNaN. 
 *   \param Err    Error code (0=no error occured). Meaning of error codes can be parsed by calling LinSolveErrMsg(Err)
 *   \param cfg    An instance of the struct lssconfig containing configuration information for the solver. This parameter
 *                 is optional, default settings are applied if it is not passed.
*/
static inline void clss( const cmatrix& A, civector& b, civector& x, civector& y, int& Err, struct lssconfig cfg) {
  cimatrix B(VecLen(b),1);
  cimatrix X(VecLen(x),1);
  cimatrix Y(VecLen(y),1);  
  B[Col(1)] = b;
  clss(A,B,X,Y,Err,cfg);
  x = X[Col(1)];
  y = Y[Col(1)];
}
#endif

//----------------------------------------------------------------------------                      

#if defined(_CXSC_CILSS_HPP) && defined(_CXSC_CILSS_INCLUDE)
#undef _CXSC_CILSS_INCLUDE
/*!
 *  Computes a verified solution of the complex linear interval system AX=B by bounding the residual of a 
 *  computed approximate solution using the Krawczyk-Operator. 
 * 
 *  This is a legacy version of the starting function. You should use the new versions which can be configured
 *  by a struct (giving a lot more options than this versions).
 * 
 *  For best performance, compile with max. compiler optimizations, activate OpenMP, and set the preprocessor variables
 *  CXSC_USE_BLAS and CXSC_USE_LAPACK (compiler options -DCXSC_USE_BLAS -DCXSC_USE_LAPACK for most compilers) and link
 *  to optimize BLAS and LAPACK libraries.
 *
 *   \param A      System matrix
 *   \param B      Right-hand side
 *   \param X      If successful this is an enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1
 *   \param Err    Error code (0=no error occured). Meaning of error codes can be parsed by calling LinSolveErrMsg(Err)
 *   \param K      Precision to use for residual calculations
 *                 If K==0 the C-XSC accumulator will be used (=maximum precision, but slower), K==1 uses normal double precision
 *                 and a splitting algorithm by Rump/Dekker to compute an improved residual, K>2 uses simulated K-fold double precision
 *                 using the DotK algorithm 
 *   \param np     Number of thread to be used for OpenMP.
 *   \param msg    If true, status messages will be put out to cout during 
 *                 computation
 *   \param lsscfg Determines, if only part 1 (LSS_ONLY_PART_ONE), only part 2 (LSS_ONLY_PART_TWO) or both parts
 *                 of the solver are to be executed. Part one is faster and works for condition numbers up to about 10^15. Part two
 *                 is much slower but can solve systems with higher condition numbers (up to about 10^30)
*/
static inline void cilss( const cimatrix& A, cimatrix& B, cimatrix& X, int& Err, int K, bool msg, int lsscfg, int np) {
  cimatrix dummy;
  struct lssconfig cfg;
  cfg.K = K;
  cfg.msg = msg;
  cfg.lssparts = lsscfg;
  cfg.threads = np;      
  SolverStart<cimatrix, cimatrix, cimatrix, cmatrix, cimatrix, cmatrix, cvector, cdotprecision, cidotprecision, civector, cinterval>(A, B, X, Err, cfg, false, dummy);
}

/*!
 *  Computes a verified solution of the complex linear point system Ax=b by bounding the residual of a 
 *  computed approximate solution using the Krawczyk-Operator. 
 * 
 *  This is a legacy version of the starting function. You should use the new versions which can be configured
 *  by a struct (giving a lot more options than this versions).
 * 
 *  For best performance, compile with max. compiler optimizations, activate OpenMP, and set the preprocessor variables
 *  CXSC_USE_BLAS and CXSC_USE_LAPACK (compiler options -DCXSC_USE_BLAS -DCXSC_USE_LAPACK for most compilers) and link
 *  to optimize BLAS and LAPACK libraries.
 *
 *   \param A      System matrix
 *   \param b      Right-hand side
 *   \param x      If successful this is an enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1
 *   \param Err    Error code (0=no error occured). Meaning of error codes can be parsed by calling LinSolveErrMsg(Err)
 *   \param K      Precision to use for residual calculations
 *                 If K==0 the C-XSC accumulator will be used (=maximum precision, but slower), K==1 uses normal double precision
 *                 and a splitting algorithm by Rump/Dekker to compute an improved residual, K>2 uses simulated K-fold double precision
 *                 using the DotK algorithm 
 *   \param np     Number of thread to be used for OpenMP.
 *   \param msg    If true, status messages will be put out to cout during 
 *                 computation
 *   \param lsscfg Determines, if only part 1 (LSS_ONLY_PART_ONE), only part 2 (LSS_ONLY_PART_TWO) or both parts
 *                 of the solver are to be executed. Part one is faster and works for condition numbers up to about 10^15. Part two
 *                 is much slower but can solve systems with higher condition numbers (up to about 10^30)
*/
static inline void cilss( const cimatrix& A, civector& b, civector& x, int& Err, int K, bool msg, int lsscfg, int np) {
  cimatrix B(VecLen(b),1);
  cimatrix X(VecLen(x),1);
  B[Col(1)] = b;
  cilss(A,B,X,Err,K,msg,lsscfg,np);
  x = X[Col(1)];
}

/*!
 *  Computes a verified solution of the complex linear point system AX=B by bounding the residual of a 
 *  computed approximate solution using the Krawczyk-Operator. 
 * 
 *  New version of the starting function. Configuration options can be set by passing an instance of the 
 *  struct lssconfig.
 * 
 *  For best performance, compile with max. compiler optimizations, activate OpenMP, and set the preprocessor variables
 *  CXSC_USE_BLAS and CXSC_USE_LAPACK (compiler options -DCXSC_USE_BLAS -DCXSC_USE_LAPACK for most compilers) and link
 *  to optimize BLAS and LAPACK libraries.
 *
 *   \param A      System matrix
 *   \param B      Right-hand side
 *   \param X      If successful this is an enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1
 *   \param Err    Error code (0=no error occured). Meaning of error codes can be parsed by calling LinSolveErrMsg(Err)
 *   \param cfg    An instance of the struct lssconfig containing configuration information for the solver. This parameter
 *                 is optional, default settings are applied if it is not passed.
*/
static inline void cilss( const cimatrix& A, cimatrix& B, cimatrix& X, int& Err, struct lssconfig cfg) {
  cimatrix ienc;
  SolverStart<cimatrix, cimatrix, cimatrix, cmatrix, cimatrix, cmatrix, cvector, cdotprecision, cidotprecision, civector, cinterval>(A, B, X, Err, cfg, false, ienc);
}

/*!
 *  Computes a verified solution of the complex linear point system Ax=b by bounding the residual of a 
 *  computed approximate solution using the Krawczyk-Operator. 
 * 
 *  New version of the starting function. Configuration options can be set by passing an instance of the 
 *  struct lssconfig.
 * 
 *  For best performance, compile with max. compiler optimizations, activate OpenMP, and set the preprocessor variables
 *  CXSC_USE_BLAS and CXSC_USE_LAPACK (compiler options -DCXSC_USE_BLAS -DCXSC_USE_LAPACK for most compilers) and link
 *  to optimize BLAS and LAPACK libraries.
 *
 *   \param A      System matrix
 *   \param b      Right-hand side
 *   \param x      If successful this is an enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1
 *   \param Err    Error code (0=no error occured). Meaning of error codes can be parsed by calling LinSolveErrMsg(Err)
 *   \param cfg    An instance of the struct lssconfig containing configuration information for the solver. This parameter
 *                 is optional, default settings are applied if it is not passed.
*/
static inline void cilss( const cimatrix& A, civector& b, civector& x, int& Err, struct lssconfig cfg) {
  cimatrix B(VecLen(b),1);
  cimatrix X(VecLen(x),1);
  B[Col(1)] = b;
  cilss(A,B,X,Err,cfg);
  x = X[Col(1)];
}

/*!
 *  Computes a verified solution of the complex linear point system AX=B by bounding the residual of a 
 *  computed approximate solution using the Krawczyk-Operator. 
 * 
 *  New version of the starting function. Configuration options can be set by passing an instance of the 
 *  struct lssconfig.
 * 
 *  For best performance, compile with max. compiler optimizations, activate OpenMP, and set the preprocessor variables
 *  CXSC_USE_BLAS and CXSC_USE_LAPACK (compiler options -DCXSC_USE_BLAS -DCXSC_USE_LAPACK for most compilers) and link
 *  to optimize BLAS and LAPACK libraries.
 *
 *   \param A      System matrix
 *   \param B      Right-hand side
 *   \param X      If successful this is an enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1
 *   \param Y      If successful this is an inner enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1. Elements for which no inner enclosure could be found are set to SignalingNaN. 
 *   \param Err    Error code (0=no error occured). Meaning of error codes can be parsed by calling LinSolveErrMsg(Err)
 *   \param cfg    An instance of the struct lssconfig containing configuration information for the solver. This parameter
 *                 is optional, default settings are applied if it is not passed.
*/
static inline void cilss( const cimatrix& A, cimatrix& B, cimatrix& X, cimatrix& Y, int& Err, struct lssconfig cfg) {
  SolverStart<cimatrix, cimatrix, cimatrix, cmatrix, cimatrix, cmatrix, cvector, cdotprecision, cidotprecision, civector, cinterval>(A, B, X, Err, cfg, true, Y);
}

/*!
 *  Computes a verified solution of the complex linear point system Ax=b by bounding the residual of a 
 *  computed approximate solution using the Krawczyk-Operator. 
 * 
 *  New version of the starting function. Configuration options can be set by passing an instance of the 
 *  struct lssconfig.
 * 
 *  For best performance, compile with max. compiler optimizations, activate OpenMP, and set the preprocessor variables
 *  CXSC_USE_BLAS and CXSC_USE_LAPACK (compiler options -DCXSC_USE_BLAS -DCXSC_USE_LAPACK for most compilers) and link
 *  to optimize BLAS and LAPACK libraries.
 *
 *   \param A      System matrix
 *   \param b      Right-hand side
 *   \param x      If successful this is an enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1
 *   \param y      If successful this is an inner enclosure of the unique solution, resized to standard index range 
 *                 with lower index bound 1. Elements for which no inner enclosure could be found are set to SignalingNaN. 
 *   \param Err    Error code (0=no error occured). Meaning of error codes can be parsed by calling LinSolveErrMsg(Err)
 *   \param cfg    An instance of the struct lssconfig containing configuration information for the solver. This parameter
 *                 is optional, default settings are applied if it is not passed.
*/
static inline void cilss( const cimatrix& A, civector& b, civector& x, civector& y, int& Err, struct lssconfig cfg) {
  cimatrix B(VecLen(b),1);
  cimatrix X(VecLen(x),1);
  cimatrix Y(VecLen(y),1);  
  B[Col(1)] = b;
  cilss(A,B,X,Y,Err,cfg);
  x = X[Col(1)];
  y = Y[Col(1)];
}
#endif

} //namespace cxsc

#endif 


