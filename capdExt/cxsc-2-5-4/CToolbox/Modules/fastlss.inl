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

/* CVS $Id: fastlss.inl,v 1.21 2014/01/30 17:23:39 cxsc Exp $ */

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


#include <string.h>        
#include <mvi_util.hpp>    

#include <rmatrix.hpp>
#include <imatrix.hpp>
#include <cmatrix.hpp>
#include <cimatrix.hpp> 
#include <idot.hpp>
#include <cidot.hpp>
#include <fstream>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif


namespace cxsc {


//! Error codes for verified linear solver
const int
  NoError      = 0,   // No error occurred
  NotSquare    = 1,   // System to be solved is not square
  DimensionErr = 2,   // Dimensions of A and b are not compatible
  InvFailed    = 3,   // System is probably singular
  VerivFailed  = 4;   // Verification failed, system is probably
                      // ill-conditioned


//----------------------------------------------------------------------------                      

/*!
  This function translates the error codes used in the verified linear system 
  solver into a corresponding textmessage, which can be used as feedback for the
  user.
  \param Err The error code
  \return A textmessage corresponding to the error code
*/
inline string LinSolveErrMsg ( int Err ) {
  string ret;

  if (Err != NoError) {
    switch (Err) {
      case NotSquare:
        ret = "System to be solved is not square"; break;
      case DimensionErr:
        ret = "Dimensions of A and b are not compatible"; break;
      case InvFailed:
        ret = "System is probably singular"; break;
      case VerivFailed:
        ret = "Verification failed, system is probably ill-conditioned";
        break;
      default:
        ret = "Error-Code not defined";
    }
  }

  return ret;
} // LinSolveErrMsg


//----------------------------------------------------------------------------                      

// Computes the transpose of a real matrix
/*
  \param A Matrix that is to be transposed
  \return The transposed of A
*/
inline rmatrix transpherm ( const rmatrix& A )
{                                            
  rmatrix  res(Lb(A,2),Ub(A,2),Lb(A,1),Ub(A,1));

  #ifdef _OPENMP
  #pragma omp parallel for default(shared)
  #endif
  for (int i = Lb(A,1); i <= Ub(A,1); i++) 
    for (int j = Lb(A,2); j <= Ub(A,2); j++) 
      res[j][i] = A[i][j];
      
  return res;
}

//----------------------------------------------------------------------------                      

// Computes the transpose of an interval matrix
/*
  \param A Matrix that is to be transposed
  \return The transposed of A
*/
inline imatrix transpherm ( const imatrix& A )
{                                            
  imatrix  res(Lb(A,2),Ub(A,2),Lb(A,1),Ub(A,1));

  #ifdef _OPENMP
  #pragma omp parallel for default(shared)
  #endif
  for (int i = Lb(A,1); i <= Ub(A,1); i++) 
    for (int j = Lb(A,2); j <= Ub(A,2); j++) 
      res[j][i] = A[i][j];
      
  return res;
}

//----------------------------------------------------------------------------                      

// Computes the hermitian of a complex matrix
/*
  \param A Matrix for which the hermitian is to be computed
  \return The hermitian of A
*/
inline cmatrix transpherm ( const cmatrix& A )
{                                            
  cmatrix  res(Lb(A,2),Ub(A,2),Lb(A,1),Ub(A,1));

  #ifdef _OPENMP
  #pragma omp parallel for default(shared)
  #endif
  for (int i = Lb(A,1); i <= Ub(A,1); i++) 
    for (int j = Lb(A,2); j <= Ub(A,2); j++) {
      res[j][i] = A[i][j];
      SetIm(res[j][i], -Im(res[j][i]));
    }
  
  return res;
}

//----------------------------------------------------------------------------                      

// Computes the hermitian of a complex interval matrix
/*
  \param A Matrix for which the hermitian is to be computed
  \return The hermitian of A
*/
inline cimatrix transpherm ( const cimatrix& A )
{                                            
  cimatrix  res(Lb(A,2),Ub(A,2),Lb(A,1),Ub(A,1));

  #ifdef _OPENMP
  #pragma omp parallel for default(shared)
  #endif
  for (int i = Lb(A,1); i <= Ub(A,1); i++) 
    for (int j = Lb(A,2); j <= Ub(A,2); j++) {
      res[j][i] = A[i][j];
      SetIm(res[j][i], -Im(res[j][i]));
    }
  
  return res;
}

//----------------------------------------------------------------------------                      

// Tries to guess zero values in successive iterates
/* The vectors x and y are successive approximations for the solution of a
   linear system of equations computed by iterative refinement. If a compo-
   nent of y is diminished by more than 'Factor', it is a good candidate for
   a zero entry. Thus, it is set to zero.
   \param x An iterate
   \param y Iterate succeding x
*/
template<typename Tvec>
inline void CheckForZeros ( Tvec& x, Tvec& y )
{
  const real Factor = 1E+5;

  for (int i = Lb(y); i <= Ub(y); i++)
    if ( abs(y[i])*Factor < abs(x[i]) ) y[i] = 0.0;
}

//----------------------------------------------------------------------------

/*
 * RelErr computes componentwise the maximum relative error of A w.r.t B. 
 * if A[i] and B[i] do not have the same sign or if B[i] = 0, then 
 * rel. error = 0 for this component.
 * @param A The new value of an iteration 
 * @param B The old value of an iteration
 * @return Relative error
 */
template <typename Tvec>
inline real RelErr(const Tvec &A, const Tvec &B)
{
  real r;
  const real zerotest = 1e6;
 
  for(int i=Lb(A);i<=Ub(A);i++) {
    if( zerotest*abs(A[i]) < abs(B[i]) || B[i] == 0 ) 
      r = 0.0;
    else 
      r = abs( (A[i]-B[i])/B[i] );
  }
  return r;
} 

//----------------------------------------------------------------------------                      


// Checks if an iterate is sufficiently accurate
/*
   The vectors x and y are successive iterates. The function returns true if
   the relative error of all components x_i and y_i is <= 10^(-12), i.e. y_i
   has about 12 correct decimals. If x_i or y_i vanishes, the relative error
   is implicitly set to zero.
   \param x An iterate
   \param y Iterate succeding x
   \return true, if y hast about 12 correct decimals, else false
*/
template<typename Tvec>
inline int Accurate ( const Tvec& x, const Tvec& y )
{
  const real Delta = 1E-12;   // Relative error bound
  int        i, n, ok;

  i = Lb(y); n = Ub(y);
  do {
    // Relative error > Delta ?
    ok = (abs(y[i] - x[i]) <= Delta * abs(y[i]) );
    i++;
  } while (ok && (i <= n));

  return ok;
} // Accurate

//----------------------------------------------------------------------------                      

//Blow function for class civector (epsilon inflation)
inline civector Blow(const civector &x, const real &eps) {
  civector y(Lb(x), Ub(x));
  
  #ifdef _OPENMP
  #pragma omp parallel for default(shared)
  #endif
  for(int i=Lb(x) ; i<=Ub(x) ; i++)
    y[i] = Blow(x[i], eps);
  
  return y;
}

//----------------------------------------------------------------------------                      

// Computes componentwise midpoint of a real matrix
/*
   Computes componentwise the midpoint of a matrix
   Real case: Just returns the matrix itself
   \param A A matrix
   \return Midpoint matrix
*/
inline rmatrix& midpoint(rmatrix &A) {
  return A;
}

//----------------------------------------------------------------------------                      

// Computes componentwise midpoint of an interval matrix
/*
   Computes componentwise the midpoint of a matrix
   Interval case: Computes the midpoint and returns it in a
   new rmatrix
   \param A A matrix
   \return Midpoint matrix
*/
inline rmatrix& midpoint(imatrix &A) {
  int m = Ub(A,ROW) - Lb(A,ROW) + 1;
  int n = Ub(A,COL) - Lb(A,COL) + 1;

  rmatrix* Am = new rmatrix(m,n);

  #ifdef _OPENMP
  #pragma omp parallel for default(shared)
  #endif
  for(int i=1 ; i<=m ; i++) {
    for(int j=1 ; j<=n ; j++) {
      (*Am)[i][j] = (Sup(A[i][j]) + Inf(A[i][j])) / 2.0;
    }
  }

  return *Am;
}

//----------------------------------------------------------------------------                      

// Computes componentwise midpoint of a complex matrix
/*
   Computes componentwise the midpoint of a matrix
   Complex case: Just returns the matrix itself
   \param A A matrix
   \return Midpoint matrix
*/
inline cmatrix& midpoint(cmatrix &A) {
  return A;
}

//----------------------------------------------------------------------------                      

// Computes componentwise midpoint of a complex interval matrix
/*
   Computes componentwise the midpoint of a matrix
   Complex interval case: Computes the midpoint and returns it in a
   new rmatrix
   \param A A matrix
   \return Midpoint matrix
*/
inline cmatrix& midpoint(cimatrix &A) {
  int m = Ub(A,ROW) - Lb(A,ROW) + 1;
  int n = Ub(A,COL) - Lb(A,COL) + 1;

  cmatrix* Am = new cmatrix(m,n);

  #ifdef _OPENMP
  #pragma omp parallel for default(shared)
  #endif
  for(int i=1 ; i<=m ; i++) {
    for(int j=1 ; j<=n ; j++) {
      (*Am)[i][j] = complex( (SupRe(A[i][j]) + InfRe(A[i][j])) / 2.0,
                             (SupIm(A[i][j]) + InfIm(A[i][j])) / 2.0);
    }
  }

  return *Am;
}

//----------------------------------------------------------------------------                      

// Determines if all elements of an interval vector are all 0
/*
  Determines if all elements of an interval vector are equal to 0
  \param zz Interval vector whose elements are to be checked
  \return true, if all elements of zz are zero
*/
inline bool isZero(const ivector &zz) {
  return (Inf(zz) == 0 && Sup(zz) == 0);
}

//----------------------------------------------------------------------------                      

// Determines if all elements of a complex interval vector are all 0
/*
  Determines if all elements of a complex interval vector are equal to 0
  \param zz Interval vector whose elements are to be checked
  \return true, if all elements of zz are zero
*/
inline bool isZero(const civector &zz) {
  return (InfRe(zz) == 0 && InfIm(zz) == 0 && SupRe(zz) == 0 && SupIm(zz) == 0);
}

//----------------------------------------------------------------------------                      

// Determines if an interval is empty
/*
  Determines if an interval x is empty, i.e. Inf(x) > Sup(x)
  \param zz Interval to check
  \return true if Inf(x) > Sup(x)
*/
inline bool Empty(const interval &i) {
  return Inf(i) > Sup(i);
}

//----------------------------------------------------------------------------                      

// Determines if a complex interval is empty
/*
  Determines if a complex interval x is empty, i.e. Inf > Sup for the real or imaginary part
  \param zz Interval to check
  \return true if interval is empty
*/
inline bool Empty(const cinterval &i) {
  return InfRe(i) > SupRe(i) || InfIm(i) > SupIm(i);
}

//----------------------------------------------------------------------------           

// Compute the max distance between two interval vectors
/*'maxDist()' computes the maximal distance between two 
   real/complex interval vectors.
   For interval x, y; dist(x,y)= Max(Abs(x.inf-y.inf), Abs(x.sup-y.sup))
   \param xx First interval
   \param yy Second interval
   \return Maximal distance between xx and yy
*/
template<class TVec1, class TVec2>
static real maxDist(const TVec1& xx, const TVec2& yy) {
  real dist = 0.0, tmp;
  int n = VecLen(xx);

  for(int i=1 ; i<=n ; i++) {
    tmp = abs(Inf(xx[i]) - Inf(yy[i]));
    if(tmp > dist) dist = tmp;
    tmp = abs(Sup(xx[i]) - Sup(yy[i]));
    if(tmp > dist) dist = tmp;
  }

  return dist;
}

// Compute the max distance between two interval matrices
/*'maxDist()' computes the maximal distance between two 
   real/complex interval matrices.
   For interval x, y; dist(x,y)= Max(Abs(x.inf-y.inf), Abs(x.sup-y.sup))
   \param xx First interval
   \param yy Second interval
   \return Maximal distance between xx and yy
*/
template<class TMat1, class TMat2>
static real maxDistMat(const TMat1& xx, const TMat2& yy) {
  real dist = 0.0, tmp;
  int n = RowLen(xx);

  for(int i=1 ; i<=n ; i++) {
    tmp = maxDist(xx[Col(i)], yy[Col(i)]);
    if(tmp > dist) dist = tmp;
  }

  return dist;
}

//----------------------------------------------------------------------------

//Fast Residual computation based on splitting algorithm by Rump/Dekker
template<typename TB>
inline TB FastResidual(const TB& B, const rmatrix& A, const rmatrix& X) {
  static real factor = 68719476737.0;
  rmatrix T = factor * A;
  rmatrix T2 = T - A;
  T = T - T2;                          
  T2 = A - T;                            

  rmatrix Y = factor * X;
  rmatrix Y2 = Y - X;
  Y = Y - Y2;                          
  Y2 = X - Y ;                            

  return (B-TB(T)*Y)-(TB(T)*Y2+TB(T2)*X);
  //return (B-T*Y)-(T*Y2+T2*X);  
}

//Fast Residual computation based on splitting algorithm by Rump/Dekker
template<typename TB>
inline TB FastResidual(const TB& B, const cmatrix& A, const cmatrix& X) {
  static real factor = 68719476737.0;
  cmatrix T = factor * A;
  cmatrix T2 = T - A;
  T = T - T2;                          
  T2 = A - T;                            

  cmatrix MX = -X;
  cmatrix Y = factor * MX;
  cmatrix Y2 = Y - MX;
  Y = Y - Y2;                          
  Y2 = MX - Y ;                            

  
  return (TB(T)*Y+B)+(TB(T)*Y2+TB(T2)*MX);
}

//Fast Residual for intervals not necessary in general
template<typename TB>
inline imatrix FastResidual(const TB& B, const imatrix& A, const rmatrix& X) {
  return B-A*X;
}

//Fast Residual for intervals not necessary in general
template<typename TB>
inline cimatrix FastResidual(const TB& B, const cimatrix& A, const cmatrix& X) {
  return B-A*X;
}


//----------------------------------------------------------------------------

// Performs the verification step of the solver
/* 
   This function 'VerificationStep()' performs the iteration
   [y] = Blow([x],Eps), [x] = [z] + [C]*[y] for k = 1,2,... until the new
   iterate [x] lies in the interior of [y] or until the maximum number of
   iterations is exceeded. The flag 'IsVerified' is set if an inclusion in
   the interior could be established.
   \param xx Overwritten with solution (if verification is possible)
   \param zz First iterate
   \param C  Matrix for iteration 
   \param IsVerified true, if a verified solution could be found, else false
   \param s The current right hand side
   \param cfg Configuration options of the solver
*/
template <typename TSysMat, typename TRhs, typename TSolution, typename TInverse, typename TC,
          typename TMidMat, typename TMidVec, typename TDot, typename TIDot, typename TVerivVec, typename TInterval>
inline void VerificationStep ( TSolution& xx, TSolution& zz, TC& C, int& IsVerified, int s, struct lssconfig cfg ) {
  TIDot dot(0.0);
 
  if(cfg.K<1) dot.set_k(0); else dot.set_k(1);
  

  int     p;
  TVerivVec yy(Ub(xx,ROW)-Lb(xx,ROW)+1);

  xx[Col(s)] = zz[Col(s)]; p = 0;                      // Initialize:  [x] = [z]
  do {
    if(cfg.msg) std::cout << " " << p+1 << ". verification step" << std::endl;        
    yy = Blow(xx[Col(s)],cfg.epsVer);             // Epsilon inflation
    int n = VecLen(yy);
    //xx = zz + C*yy;                    // New iterate: [x] = [z]+[C]*[y]
    #ifdef _OPENMP
    #pragma omp parallel for firstprivate(dot) default(shared)
    #endif
    for(int i=1 ; i<=n ; i++) {
       dot = zz[i][s];
       accumulate(dot,C[Row(i)], yy);
       xx[i][s] = rnd(dot);
    }    
    IsVerified = in((TVerivVec)xx[Col(s)],yy);            // Inclusion in the interior?
    p++;
  } while (!IsVerified && (p < cfg.maxIterVer));
}

//----------------------------------------------------------------------------                      

// Computes a verified solution of a square linear (interval) system of equations A*x=b
/* 
    An approximate inverse 'R' of 'A' is computed by calling the
    function 'MatInv_Aprx()'. Then an approximate solution 'x' is computed
    applying a conventional residual iteration. For the final verification, 
    an interval residual iteration is performed. An enclosure of the
    unique solution is returned in the interval vector 'xx'. 
    
    If this first step fails, the solver will try to compute a verified
    solution by using an approximate inverse of double length. This second
    step takes considerably longer than the first step, because all computations
    must me performed using high precision scalar products.
    
    The user can choose if only one of the two parts or both should be used.

    The solver can use BLAS, LAPACK and OpenMP if activated. Higher precision
    dot products are used for the residual computation.

    \param A Matrix of the system (must be square!)
    \param b Right-hand side of the system
    \param xx Enclosure of the unique solution, resized to standard index range 
              with lower index bound 1
    \param Err Error code
    \param K Precision to use for residual computation
    \param cfg A struct containing configuration options for the solver
    \param inner true, if an inner enclosure should be computed
    \param ienc A vector containing the inner enclosure on exit
*/
template <typename TSysMat, typename TRhs, typename TSolution, typename TInverse, typename TC,
          typename TMidMat, typename TMidVec, typename TDot, typename TIDot, typename TVerivVec, typename TInterval>
inline void LinSolveMain ( TSysMat &A, TRhs &b, TSolution& xx, int& Err, struct lssconfig cfg, bool inner, TSolution& ienc) {
  int      n = Ub(A,ROW) - Lb(A,ROW) + 1;        // Length of the columns of 'A'
  int      m = Ub(A,COL) - Lb(A,COL) + 1;        // Length of the rows of 'A'
  int      o = Ub(b,ROW) - Lb(b,ROW) + 1;        // Length of the columns of 'A'
  int      rhs = Ub(b,COL) - Lb(b,COL) + 1;      // Right hand sides
  int      IsVerified, i, j=0, k;
  int      bRowLow, bColLow, ARowLow, AColLow;   // Lower index bounds od 'A' and 'b'
  bool     C_computed = false;                   // Has matrix [C] been computed?
  TInverse R(n,n);                               // To store the inverse of 'A' or mid(A)

  if (n != m) { 
    // Error: 'A' is not square
    Err = NotSquare;
    return; 
  }

  if (n != o) { 
    // Error: Dimensions of 'A' and 'b' are not compatible
    Err = DimensionErr;
    return; 
  }     
  
  if(n <= 0) {
    Err = DimensionErr;
    return;
  }

  if(cfg.lssparts < 0  ||  cfg.lssparts > 2)
    cfg.lssparts = LSS_BOTH_PARTS;

#ifdef _OPENMP
    if(cfg.threads > 0) {
       omp_set_num_threads(cfg.threads);
    } else {
       cfg.threads = omp_get_max_threads(); 
       omp_set_num_threads(cfg.threads);
    }
    if(cfg.msg) std::cout << "Using " << cfg.threads << " thread(s) for computations" << std::endl;    
#endif

  //Compute mid of A and b (handled as references, point matrices do 
  //not need this computation)
  TMidMat& Am = midpoint(A);
  TMidMat& bm = midpoint(b);
  
  //Inversion of A  
  if(cfg.msg) std::cout << "Computing approximate inverse of A" << std::endl;
  MatInvAprx(Am,R,Err);  
  if (Err != NoError) { 
    // Error: Inversion failed
    Err = InvFailed; 
    return; 
  }


  // Start algorithm
  //----------------
  TMidVec       x(n), y(n), d(n);       // Allocate dynamic arrays
  TSolution     yy(n,rhs), zz(n,rhs);
  TC            C(n,n);
  TDot          dot;              //DotK-Objects for accurate dot product
  TIDot         idot;

  dot.set_k(cfg.K);
  idot.set_k(cfg.K);
  opdotprec = cfg.K;

  // Normalize index range of 'A', 'R' and 'b' to standard range 1..n
  // and resize 'xx' to length 'n'. Save lower bounds of 'b' and 'A'
  // to restore them before leaving 'LinSolveMain()'.
  //-----------------------------------------------------------------
  Resize(xx,n,rhs);
  ARowLow = Lb(A,ROW); SetLb(A,ROW,1);
  AColLow = Lb(A,COL); SetLb(A,COL,1);
  bRowLow = Lb(b,ROW); SetLb(b,ROW,1);
  bColLow = Lb(b,COL); SetLb(b,COL,1);
  SetLb(R,ROW,1); SetLb(R,COL,1); 
  
  //Start part one of the algorithm
  if(cfg.lssparts != LSS_ONLY_PART_TWO) {
    
    //Matrix Mode, especially useful for multiple right hand sides
    if(cfg.matrixMode) {
      opdotprec = cfg.K;
      TMidMat X(n,rhs);
      if(cfg.msg) std::cout << "Computing approximate solution x" << std::endl;     
      X = R*bm;
      if(cfg.K == 1) 
	X += R*FastResidual(bm,Am,X);
      else
	X += R * (bm - Am * X);
      if(cfg.msg) std::cout << "Computing [z] = R*(b - A*x)" << std::endl;      
      if(cfg.K == 1)
	zz = R*FastResidual(TSolution(b),A,X);
      else
	zz = R * (b - A * TSolution(X));
      if(cfg.msg) std::cout << "Computing [C] = I-R*A" << std::endl;
      // Compute [C] = ##(I-R*A)
      #ifdef CXSC_USE_BLAS
        IminusRA(R,A,C);
      #else
        idot = 0.0;
        idot.set_k(1);
        #ifdef _OPENMP
        #pragma omp parallel for firstprivate(idot) private(i,j) default(shared)
        #endif
        for(i=1 ; i<=n ; i++) {
          for(j=1 ; j<=n ; j++) {
            if(i==j) idot=1.0; else idot=0.0;
            accumulate(idot,-R[i],A[Col(j)]);
            C[i][j] = rnd(idot);
          } 
        }
        idot.set_k(cfg.K);
      #endif    
      int p = 0;
      if(cfg.msg) std::cout << "Starting verification" << std::endl;        
      xx = zz;
      opdotprec = 1;
      do {
        if(cfg.msg) std::cout << " " << p+1 << ". verification step" << std::endl;
	for(i=1 ; i<=rhs ; i++)
          yy[Col(i)] = Blow(xx[Col(i)],cfg.epsVer);             // Epsilon inflation
        xx = zz + C*yy;                    // New iterate: [x] = [z]+[C]*[y]
        IsVerified = in((TVerivVec)xx[Col(1)],(TVerivVec)yy[Col(1)]);            // Inclusion in the interior?
	i = 2;
	while(IsVerified && i<=rhs) {
          IsVerified = in((TVerivVec)xx[Col(i)],(TVerivVec)yy[Col(i)]);          // Inclusion in the interior?
	  i++;
	}
        p++;
      } while (!IsVerified && (p < cfg.maxIterVer));

      if(IsVerified) {
	if(inner) {
	  std::cerr <<  "Warning: Computations of inner enclosure not possible in Matrix Mode!" << std::endl;
	  ienc = SignalingNaN;
	}
	
        if(cfg.refinement) {
	  if(cfg.msg) std::cout << "Starting iterative refinement" << std::endl;
          p = 0; 
          
          do {	   
            if(cfg.msg) std::cout << " " << p+1 << ". refinement step " << std::endl;
            yy = xx;
            xx = zz + C * xx;
            xx &= yy;
            p++;
          } while (maxDistMat(xx, yy) > cfg.epsRef && (p < cfg.maxIterRef)); 
	}
	xx += X;	
      }	
      opdotprec = cfg.K;
    }
    //End Matrix Mode version

    //Loop over all right hand sides if matrix mode is not used
    for(int s=1 ; s<=rhs && Err!=VerivFailed && !cfg.matrixMode ; s++) {

	if(cfg.msg) {
           std::cout << std::endl << "Computing result for right hand side " << s << std::endl;
	   std::cout << "Residual iteration..." << std::endl;
        }

        // Real residual iteration
	x = R*bm[Col(s)]; k = 0; 
	if(cfg.K == 1)
	  x += TMidVec(R*FastResidual(TMidMat(bm[Col(s)]),Am,TMidMat(x)));
	else {
	  do {              
	      k++;
	      if(cfg.msg) std::cout << " "  << k << ". iteration step" << std::endl;
	      y = x;
	      // Compute: d = #*(b-A*y)
	      #ifdef _OPENMP
	      #pragma omp parallel for firstprivate(dot) default(shared)
	      #endif
	      for (i = 1; i <= n; i++) {   
	        dot = bm[i][s];
	        accumulate_approx(dot,-Am[i],y);
	        d[i] = rnd(dot);
	      } 
	    
	      // Compute: x = #*(y+R*d)
	      #ifdef _OPENMP
	      #pragma omp parallel for firstprivate(dot) default(shared)
	      #endif
	      for (i = 1; i <= n; i++) {   
	        dot = y[i];
	        accumulate_approx(dot,R[i],d);
	        x[i] = rnd(dot);
	      }
	    
	      CheckForZeros(y,x);
	  } while (!Accurate(y,x) && (k < cfg.maxIterResCorr));
	} 

	// [z] = R*(b - A*x).
	//-------------------------------------------------------	
	if(cfg.msg) std::cout << "Computing [z] = R*(b - A*x)" << std::endl;
	if(cfg.K == 1) {
	  zz[Col(s)] = TVerivVec(R * FastResidual(TSolution(b[Col(s)]), A, TMidMat(x)));
	} else {
	  // Compute: d = #*(b-A*x)
          #ifdef _OPENMP
          #pragma omp parallel for firstprivate(idot) default(shared)
          #endif
	  for (i = 1; i <= n; i++) {              
	    idot = b[i][s];
	    accumulate(idot,-A[i],x);
	    yy[i][s] = rnd(idot);
	    //TInterval tmp = rnd(idot);
            //d[i] = Inf(tmp) + 0.5 * (Sup(tmp)-Inf(tmp));
	  }

// 	  // Compute: [y] = ##(b-A*x-d)
//           #ifdef _OPENMP
//           #pragma omp parallel for firstprivate(idot) default(shared)
//           #endif
// 	  for (i = 1; i <= n; i++) {              
// 	    idot = b[i][s];
// 	    accumulate(idot,-A[i],x);
// 	    idot += -d[i];
// 	    yy[i][s] = rnd(idot);
// 	  }

	  // Now b-A*x lies in d+[y]. Thus R*(b-A*x) lies in [z] = ##(R*d+R*[y])
	  //--------------------------------------------------------------------
          #ifdef _OPENMP
          #pragma omp parallel for firstprivate(idot) default(shared)
          #endif
	  for (i = 1; i <= n; i++) {
	    idot = 0.0;
	    //accumulate(idot,R[i],d);
	    accumulate(idot,R[i],yy[Col(s)]);
	    zz[i][s] = rnd(idot);
	  }
	}

	// If R*(b-A*x) = 0 then x is an exact solution else try to find a
	// verified enclosure [x] for the absolute error of x.
	//----------------------------------------------------------------
	
	if(cfg.msg) std::cout << "Verifying..." << std::endl;
	if (isZero(zz[Col(s)])) {
	  xx[Col(s)] = x;
          IsVerified = 1;
	} else {
          if(!C_computed) {
	    if(cfg.msg) std::cout << " Computing [C] = I-R*A" << std::endl;
	    // Compute [C] = ##(I-R*A)
            #ifdef CXSC_USE_BLAS
              opdotprec = 1;
              IminusRA(R,A,C);
              opdotprec = cfg.K;
            #else
              idot = 0.0;
              idot.set_k(1);
              #ifdef _OPENMP
              #pragma omp parallel for firstprivate(idot) private(i,j) default(shared)
              #endif
              for(i=1 ; i<=n ; i++) {
                for(j=1 ; j<=n ; j++) {
                 if(i==j) idot=1.0; else idot=0.0;
                 accumulate(idot,-R[i],A[Col(j)]);
                 C[i][j] = rnd(idot);
                } 
              }
              idot.set_k(cfg.K);
            #endif
            C_computed = true;
          }

          // Attempt to compute [x]
	  VerificationStep<TSysMat, TRhs, TSolution, TInverse, TC, TMidMat, TMidVec, TDot, TIDot, TVerivVec, TInterval>(xx,zz,C,IsVerified,s,cfg);   

	  if ( !IsVerified )
		Err = VerivFailed;
	  else {
	        //Iterative refinement
                if(cfg.refinement) {
		  if(cfg.msg) std::cout << "Starting iterative refinement" << std::endl;
                  int        p=0; 
                  TVerivVec yy(n); 
  
                  do {	   
                    if(cfg.msg) std::cout << " " << p+1 << ". refinement step " << std::endl;
                    yy = xx[Col(s)];
                    xx[Col(s)] = zz[Col(s)] + C * TVerivVec(xx[Col(s)]);
                    xx[Col(s)] &= yy;
                    p++;
                  } while (maxDist(xx[Col(s)], yy) > cfg.epsRef && (p < cfg.maxIterRef)); 
		}
		
		//Inne enclosure
		if(inner) {
	          if(cfg.msg) std::cout << "Computing inner enclosure" << std::endl;
		  //Compute inner enclosure
		  TVerivVec D(n);
		  TVerivVec zzi(n), zze(n);
		  idot.set_dotprec(0);
		  opdotprec = 1;

		  D = C * TVerivVec(xx[Col(s)]);

		  #ifdef _OPENMP
        	  #pragma omp parallel for firstprivate(idot) private(i) default(shared)
		  #endif
		  for(i=1 ; i<=n ; i++) {
		    idot = b[i][s];
		    accumulate(idot, -A[i], x);
		    UncheckedSetInf(zzi[i], rnd(Inf(idot),RND_UP));
		    UncheckedSetSup(zzi[i], rnd(Sup(idot),RND_DOWN));
		  }

		  #ifdef _OPENMP
        	  #pragma omp parallel for firstprivate(idot) private(i) default(shared)
		  #endif
		  for(i=1 ; i<=n ; i++) {
		    idot = 0.0;
		    accumulate(idot, R[i], zzi);
		    UncheckedSetInf(zze[i], rnd(Inf(idot), RND_UP));
		    UncheckedSetSup(zze[i], rnd(Sup(idot), RND_DOWN));		    
		  }

                  idot.set_dotprec(cfg.K);			
		  opdotprec = cfg.K;

		  
		  for(i=1 ; i<=n ; i++) {
		    UncheckedSetInf(ienc[i][s], addu(addu(x[i],Inf(zze[i])), Sup(D[i])));
		    UncheckedSetSup(ienc[i][s], addd(addd(x[i],Sup(zze[i])), Inf(D[i])));
		    if(Empty(ienc[i][s]))
		      ienc[i][s] = TInterval(SignalingNaN);
		  }
                }
		
		xx[Col(s)] = x + xx[Col(s)];                          // The exact solution lies x+[x]

          }
	}

    }
  
    if(IsVerified) {
      // Restore lower index bounds of 'A' and 'b'
      //------------------------------------------
      SetLb(b,ROW,bRowLow); SetLb(b,COL,bColLow);
      SetLb(A,ROW,ARowLow); SetLb(A,COL,AColLow);
      if(cfg.msg) std::cout << "Verified solution found" << std::endl;
      if((TMidMat*)&A != &Am) delete &Am;
      if((TMidMat*)&b != &bm) delete &bm;
      return; //done
    } else if(cfg.lssparts == LSS_ONLY_PART_ONE) {
      // Restore lower index bounds of 'A' and 'b'
      //------------------------------------------
      SetLb(b,ROW,bRowLow); SetLb(b,COL,bColLow);
      SetLb(A,ROW,ARowLow); SetLb(A,COL,AColLow);
      if((TMidMat*)&A != &Am) delete &Am;
      if((TMidMat*)&b != &bm) delete &bm;
      return; //done
    }
    
  } //end lsscfg!=LSS_ONLY_PART_TWO
  
//******************************************************************************

  //If verification not successful or part two selected:
  //Restart using inverse of double length (Rumps Device)
  if(cfg.lssparts == LSS_BOTH_PARTS) { 
     if(cfg.msg) std::cout << std::endl << "Part 1 not successful, starting Part 2" << std::endl;
  } else {
     if(cfg.msg) std::cout << "Starting Part 2" << std::endl; 
  }
     
  Err = NoError;
  C_computed = false;

  if(cfg.msg) std::cout << "Computing R2 = R*A" << std::endl;
  //R2=R*A
  TInverse R2(n,n); 
#ifdef _OPENMP
#pragma omp parallel for firstprivate(dot,j) default(shared)
#endif
  for(i=1 ; i<=n ; i++) {
     for(j=1 ; j<=n ; j++) {
       dot = 0.0;
       accumulate_approx(dot,R[Row(i)], Am[Col(j)]);
       R2[i][j] = rnd(dot);
     }
  } 

  if(cfg.msg) std::cout << "Computing approximate inverse of R2" << std::endl;
  //Invert R2    
  MatInvAprx(R2,R2,Err);
  if (Err != NoError)                   // Error: Inversion failed
    { Err = InvFailed; return; }
  
  if(cfg.msg) std::cout << "Computing R2 = R2*R - D" << std::endl;
  //R2:=#*(R2*R - D)
  TInverse R1temp(n,n);
  TInverse R2temp(n,n);
#ifdef _OPENMP
#pragma omp parallel for firstprivate(dot,j) default(shared)
#endif
  for(i=1 ; i<=n ; i++) {
    for(j=1 ; j<=n ; j++) {
      dot = 0.0;
      accumulate_approx(dot,R2[i],R[Col(j)]);
      R1temp[i][j] = rnd(dot);
      dot += -R1temp[i][j];
      R2temp[i][j] = rnd(dot);
    }
  }
      
  R2 = R2temp;
  R = R1temp;

  //delete temp matrices  
  R1temp = TInverse();
  R2temp = TInverse();
  
  TMidVec x1(n), x1temp(n);
  TSolution yytemp(n,rhs); 
  real p(0.0);
  

  //loop over all right hand sides	
  for(int s=1 ; s<=rhs && Err!=VerivFailed ; s++) {
     if(cfg.msg) {
       std::cout << std::endl << "Computing result for right hand side " << s << std::endl;
       std::cout << "Preparing residual iteration" << std::endl;
     }

     //floating point defect iteration
     p = 0.0;	

     // x1 := #*( R1*bm + R2*bm );
#ifdef _OPENMP
#pragma omp parallel for firstprivate(dot) default(shared)
#endif
     for (i=1 ; i<=n ; i++) {  
        dot = 0.0;
        accumulate_approx(dot,R[Row(i)],bm[Col(s)]);
        accumulate_approx(dot,R2[Row(i)],bm[Col(s)]);
        x1[i] = rnd(dot);
     }

     // x = #*( R1*bm + R2*bm - x1):
#ifdef _OPENMP
#pragma omp parallel for firstprivate(dot) default(shared)
#endif
     for (i=1 ; i<=n ; i++) {  
        dot = -x1[i];
        accumulate_approx(dot,R[Row(i)],bm[Col(s)]);
        accumulate_approx(dot,R2[Row(i)],bm[Col(s)]);
        x[i] = rnd(dot);
     }

     k = 0;
     if(cfg.msg) std::cout << "Residual iteration..." << std::endl;  
     // iteration x = x + (R+R2)*(b-Ax), x = x1 + x)
     do {
        k++;
        if(cfg.msg) std::cout << " "  << k << ". iteration step" << std::endl;

        // d = #*(bm - AM*x1 - AM*x)
#ifdef _OPENMP
#pragma omp parallel for firstprivate(dot) default(shared)
#endif
        for (i=1 ; i<=n ; i++) {  
           dot = bm[i][s];
           accumulate_approx(dot,-Am[Row(i)],x1);
           accumulate_approx(dot,-Am[Row(i)],x);
           d[i] = rnd(dot);
        }

        // y0 := #*(x + R*d + R2*d);
#ifdef _OPENMP
#pragma omp parallel for firstprivate(dot) default(shared)
#endif
        for (i=1 ; i<=n ; i++) {  
           dot = x[i];
           accumulate_approx(dot,R[Row(i)],d);
           accumulate_approx(dot,R2[Row(i)],d);
           x1temp[i]= rnd(dot);
        }
        d = x1temp;

        p = RelErr(x1+d,x1+x);
        d = x1 + d;

        // x0 := #*(x1 + x - d)
#ifdef _OPENMP
#pragma omp parallel for firstprivate(dot) default(shared)
#endif
        for (i=1 ; i<=n ; i++) {
           dot = x1[i];
           dot += x[i];
           dot += -d[i];
           x[i] = rnd(dot); 
        }
        x1 = d;

     } while ( p>=1e-6 && k<cfg.maxIterResCorr );  
     

     // compute enclosure y+Y1 of the residuum b-A*x1 of the approximation x1 and 
     // initialize Y1:= Z:= (R1+R2)*(b-A*x1), C:= I-(R1+R2)*A   
     if(cfg.msg) std::cout << "Computing enclosure of the residuum" << std::endl;

     
     // y := #*(b-A*x1)  ====> y:=MID(##(b - A*x1))
#ifdef _OPENMP
#pragma omp parallel for firstprivate(idot) default(shared)
#endif
     for (i=1 ; i<=n ; i++) {
        idot = b[i][s];
        accumulate(idot,-A[Row(i)],x1);
        TInterval tmp = rnd(idot);
        d[i] = Inf(tmp) + (Sup(tmp)-Inf(tmp)) * 0.5;
     } 

     // Y1 := ##(b-A*x1-y)           
      idot.set_k(0);
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
     for (i=1 ; i<=n ; i++) {
#ifdef _OPENMP
          TIDot id; //private(idot) leads to wrong results!?
	  id.set_k(0);
          id = b[i][s];
          accumulate(id,-A[Row(i)],x1);
          id += -d[i];
          yy[i][s] = rnd(id);
#endif
#ifndef _OPENMP
          idot = b[i][s];
          accumulate(idot,-A[Row(i)],x1);
          idot += -d[i];
          yy[i][s] = rnd(idot);
#endif
     }
     idot.set_k(cfg.K);
     
     // Y1 := ##(R*y + R2*y + R*Y1 + R2*Y1 );   
#ifdef _OPENMP
#pragma omp parallel for firstprivate(dot,idot) default(shared)
#endif
     for (i=1 ; i<=n ; i++) {  
        dot = 0.0;
        accumulate(dot,R[Row(i)],d);
        accumulate(dot,R2[Row(i)],d);
        idot = TInterval(rnd(dot,RND_DOWN), rnd(dot,RND_UP));
        accumulate(idot,R[Row(i)],yy[Col(s)]);
        accumulate(idot,R2[Row(i)],yy[Col(s)]);
        yytemp[i][s] = rnd(idot);
     }

     yy[Col(s)] = yytemp[Col(s)];
     zz[Col(s)] = yy[Col(s)];  
     
     dot.set_k(cfg.K);
     idot.set_k(cfg.K);

     if(cfg.msg) std::cout << "Verifying..." << std::endl;

     if (isZero(zz[Col(s)])) {
        xx[Col(s)] = x1; // exact solution! (however, not necessarily unique!)
        IsVerified = 1;
     } else {
        if(!C_computed) {
           if(cfg.msg) std::cout << " Computing [C] = I - R*A - R2*A" << std::endl;
           // C:= ## (ID(A) - R*A - R2*A);
#ifdef _OPENMP
#pragma omp parallel for firstprivate(j) default(shared) 
#endif
           for (i=1 ; i<=n ; i++) {
#ifdef _OPENMP
              TIDot id; //private(idot) leads to wrong results!?
	      id.set_k(cfg.K);
              for (j=1 ; j<=n ; j++) {
                 id = ( i == j) ? 1.0 : 0.0;
                 accumulate(id,-R[i],A[Col(j)]);
                 accumulate(id,-R2[i],A[Col(j)]);
                 C[i][j] = rnd(id);
              }
#endif
#ifndef _OPENMP
              for (j=1 ; j<=n ; j++) {
                 idot = ( i == j) ? 1.0 : 0.0;
                 accumulate(idot,-R[i],A[Col(j)]);
                 accumulate(idot,-R2[i],A[Col(j)]);
                 C[i][j] = rnd(idot);
              }
#endif
           }
           
           C_computed = true;
        }

        // interval iteration until inclusion is obtained (or max. iteration count)
        VerificationStep<TSysMat, TRhs, TSolution, TInverse, TC, TMidMat, TMidVec, TDot, TIDot, TVerivVec, TInterval>(yy,zz,C,IsVerified,s,cfg);   // Attempt to compute [x]
        if ( !IsVerified )
           Err = VerivFailed;
        else {

	   //Refinement step
           if(cfg.refinement) {
	     if(cfg.msg) std::cout << "Starting iterative refinement" << std::endl;
             int        p=0; 
             TVerivVec yy(n); 
  
             do {	   
               if(cfg.msg) std::cout << " " << p+1 << ". refinement step " << std::endl;
               yy = xx[Col(s)];
               xx[Col(s)] = zz[Col(s)] + C * TVerivVec(xx[Col(s)]);
               xx[Col(s)] &= yy;
               p++;
             } while (maxDist(xx[Col(s)], yy) > cfg.epsRef && (p < cfg.maxIterRef)); 
	   }
	  
	   //Inner enclosure
	   if(inner) {
	     
	     if(cfg.msg) std::cout << "Computing inner enclosure" << std::endl;

	     //Compute inner enclosure
	     TVerivVec D(n);
	     TVerivVec zzi(n), zze(n);
             idot.set_dotprec(0);
	     opdotprec = 0;

	     D = C * TVerivVec(xx[Col(s)]);

	     #ifdef _OPENMP
             #pragma omp parallel for firstprivate(idot) private(i) default(shared)
	     #endif
	     for(i=1 ; i<=n ; i++) {
		idot = b[i][s];
		accumulate(idot, -A[i], x1);
		UncheckedSetInf(zzi[i], rnd(Inf(idot),RND_UP));
		UncheckedSetSup(zzi[i], rnd(Sup(idot),RND_DOWN));
	     }

             #ifdef _OPENMP
             #pragma omp parallel for firstprivate(idot) private(i) default(shared)
             #endif
             for(i=1 ; i<=n ; i++) {
                idot = 0.0;
                accumulate(idot, R[i], zzi);
                accumulate(idot, R2[i], zzi);		    
                UncheckedSetInf(zze[i], rnd(Inf(idot), RND_UP));
                UncheckedSetSup(zze[i], rnd(Sup(idot), RND_DOWN));		    
             }

             idot.set_dotprec(cfg.K);			
             opdotprec = cfg.K;

             for(i=1 ; i<=n ; i++) {
                UncheckedSetInf(ienc[i][s], addu(addu(x1[i],Inf(zze[i])), Sup(D[i])));
                UncheckedSetSup(ienc[i][s], addd(addd(x1[i],Sup(zze[i])), Inf(D[i])));
                if(Empty(ienc[i][s]))
                  ienc[i][s] = TInterval(SignalingNaN);
             }
	  }

          xx[Col(s)] = x1 + yy[Col(s)];                          // The exact solution lies x+[x]     

      }
    }
  }

  

  if(cfg.msg) {
    if(IsVerified) 
      std::cout << "Verified solution found!" << std::endl;
  }

  if((TMidMat*)&A != &Am) delete &Am;
  if((TMidMat*)&b != &bm) delete &bm;

  // Restore lower index bounds of 'A' and 'b'
  //------------------------------------------
  SetLb(A,ROW,ARowLow); SetLb(A,COL,AColLow);
  SetLb(b,ROW,bRowLow); SetLb(b,COL,bColLow); 
} //LinSolveMain

//----------------------------------------------------------------------------                      

// Entry function for solver
/*
  This function starts the solver. If A ist not square it will transform the System into an equivalent 
  square system. A part of the solution of this new System then is the solution of the original 
  under- or overdetermined System (=the smallest vector minimizing the square of the norm of Ax-b)

    \param A Matrix of the system
    \param b Right-hand side of the system
    \param xx Enclosure of the unique solution, resized to standard index range 
              with lower index bound 1
    \param Err Error code
    \param cfg struct with configuration options
    \param inner Shall an inner enclosure be computed?
    \param ienc Vector containing inner enclosure on exit
*/
template <typename TSysMat, typename TRhs, typename TSolution, typename TInverse, typename TC,
          typename TMidMat, typename TMidVec, typename TDot, typename TIDot, typename TVerivVec, 
          typename TInterval>
inline void SolverStart ( const TSysMat& A, TRhs& b, TSolution& xx, int& Err, struct lssconfig cfg, bool inner, TSolution& ienc) {
  int m   = Ub(A,ROW) - Lb(A,ROW) + 1;
  int n   = Ub(A,COL) - Lb(A,COL) + 1;
  int rhs = Ub(b,COL) - Lb(b,COL) + 1;
  int dim = m+n;
  int olddotprec = opdotprec;
  
  if(m == n) {

    //Square
    LinSolveMain<TSysMat, TRhs, TSolution, TInverse, TC, TMidMat, TMidVec, TDot, TIDot, TVerivVec, TInterval>(const_cast<TSysMat&>(A), b, xx, Err, cfg, inner, ienc);

  } else if (m > n) {

    //Overdetermined
    if(cfg.msg) std::cout << "Overdetermined system, generating equivalent square system" << std::endl; 

    TSysMat BIG_A(1,dim,1,dim);
    TRhs BIG_b(1,dim,1,rhs);
    TSolution  BIG_xx(1,dim,1,rhs), BIG_ienc;

    if(inner) Resize(BIG_ienc, 1, dim, 1, rhs);

    BIG_A( 1,m, 1,n ) = A;
    BIG_A( 1,m, n+1,n+m ) = -Id(rmatrix(m,m));
    BIG_A( m+1,m+n, 1,n ) = 0.0;
    BIG_A( m+1,m+n, n+1,n+m ) = transpherm(A);
    
    BIG_b( 1,m,1,rhs ) = b;
    BIG_b( m+1,m+n,1,rhs ) = 0.0;
    
    LinSolveMain<TSysMat, TRhs, TSolution, TInverse, TC, TMidMat, TMidVec, TDot, TIDot, TVerivVec, TInterval>( BIG_A, BIG_b, BIG_xx, Err, cfg, inner, BIG_ienc); 
 
    xx = BIG_xx( 1,n,1,rhs );    

    if(inner) ienc = BIG_ienc( 1,n,1,rhs );

  } else if (m < n) {

    //Underdetermined
    if(cfg.msg) std::cout << "Underdetermined system, generating equivalent square system" << std::endl; 

    TSysMat BIG_A(1,dim,1,dim);
    TRhs BIG_b(1,dim,1,rhs);
    TSolution  BIG_xx(1,dim,1,rhs), BIG_ienc;

    if(inner) Resize(BIG_ienc, 1, dim, 1, rhs);
    
    BIG_A( 1,n, m+1,m+n ) = -Id(rmatrix(n,n));
    BIG_A( 1,n, 1,m ) = transpherm(A);
    BIG_A( n+1,n+m, 1,m ) = 0.0;
    BIG_A( n+1,n+m, m+1,m+n ) = A;
    
    BIG_b( 1,n,1,rhs ) = 0.0;
    BIG_b( n+1,n+m,1,rhs ) = b;

    LinSolveMain<TSysMat, TRhs, TSolution, TInverse, TC, TMidMat, TMidVec, TDot, TIDot, TVerivVec, TInterval>( BIG_A, BIG_b, BIG_xx, Err, cfg, inner, BIG_ienc);    

    xx = BIG_xx( m+1,m+n,1,rhs );    
    
    SetLb(xx,ROW,1);

    if(inner) ienc = BIG_ienc( m+1,m+n,1,rhs );
  }
  
  opdotprec = olddotprec;
}


} //namespace cxsc



