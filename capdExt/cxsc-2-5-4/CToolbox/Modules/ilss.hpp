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

/* CVS $Id: ilss.hpp,v 1.15 2014/01/30 17:23:39 cxsc Exp $ */

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

#ifndef _CXSC_ILSS_HPP
#define _CXSC_ILSS_HPP

#ifndef _CXSC_LSS_CONSTANTS_DEFINED
#define _CXSC_LSS_CONSTANTS_DEFINED

#include <real.hpp>

namespace cxsc {
  
//Control constants
static const int
  LSS_ONLY_PART_ONE = 0, //Use only first part of algorithm
  LSS_ONLY_PART_TWO = 1, //Use only second part (using R1+R2 as approximate inverse)
                         //second part can solve systems with higher condition numbers
                         //but is much slower
  LSS_BOTH_PARTS    = 2; //Try part one first, if it fails try part two
  
/**
 * This struct can be used to configure the solver. If you want to change a default setting,
 * create an instance of this struct, change the desired settings, and use the instance of the 
 * strut as the last parameter of the solver function call (see below).
 */
struct lssconfig {
  int   K;              //Dot product precision for residual computations
  bool  msg;            //Status message output during solver run?
  int   lssparts;       //Solver stages to use (see above)
  int   threads;        //Number of threads to use for OpenMP (only has an effect if OpenMP ist activated during compilation).
                        //Using a value <=0 means use all available threads (set for example by OMP_NUM_THREADS environment variable)
  int   maxIterResCorr; //Maximum number of iterations during residual correction (not available for K=1).
                        //Higher settings might increase result quality, but also increases time costs
  int   maxIterVer;     //Maximum number of iterations during the verification step.
                        //For systems with condition numbers approaching 10^15 for part one or 10^30 for part two
                        //higher than default settings might be necessary to find a verified result
  bool  refinement;     //Perform an iterative refinement? This can sometimes improve the computed enclosure slightly with little cost
  int   maxIterRef;     //Maximum number of iterations during the refinement step
  real  epsVer;         //Epsilon for the verification step (used for the Epsilon inflation during verification)
  real  epsRef;         //Epsilon for the refinement step (stopping criterion for refinement step)
  bool  matrixMode;     //Use only operators for computations. Faster for multiple right hand sides,
                        //but might lose some accuracy. Use matrix mode and K=1 if you want to compute a verified inverse of a matrix!

  lssconfig() : K(2),msg(false),lssparts(LSS_BOTH_PARTS),threads(-1),maxIterResCorr(10), 
                maxIterVer(5), refinement(false), maxIterRef(5), epsVer(0.1), epsRef(1e-5), matrixMode(false) {}
};

}
#endif

#include <imatrix.hpp>

/**
 * All versions of the function ilss compute the verified solution of the real interval linear system AX=B. The variable
 * Err returns an error code (0 meaning no error occured). The according error message is given by calling LinSolveErrMsg(Err).
 * 
 * For best performance, compile with max. compiler optimizations, activate OpenMP, and set the preprocessor variables
 * CXSC_USE_BLAS and CXSC_USE_LAPACK (compiler options -DCXSC_USE_BLAS -DCXSC_USE_LAPACK for most compilers) and link
 * to optimize BLAS and LAPACK libraries.
 */

namespace cxsc {

// Legacy versions. Additional parameters have same meaning as in struct described above  
//! Entry function of interval linear sytems solver
static inline void  ilss(const cxsc::imatrix& A, cxsc::imatrix& B, cxsc::imatrix& X, int& Err, int K, bool msg=false, int lssparts=LSS_BOTH_PARTS, int threads=-1);
//! Entry function of interval linear sytems solver
static inline void  ilss(const cxsc::imatrix& A, cxsc::ivector& B, cxsc::ivector& X, int& Err, int K, bool msg=false, int lssparts=LSS_BOTH_PARTS, int threads=-1);

//New versions using struct. Pass an instance of the struct as last parameter to use different than default settings.
//! Entry function of interval linear sytems solver
static inline void  ilss(const cxsc::imatrix& A, cxsc::imatrix& B, cxsc::imatrix& X, int& Err, struct lssconfig=lssconfig());
//! Entry function of interval linear sytems solver
static inline void  ilss(const cxsc::imatrix& A, cxsc::ivector& B, cxsc::ivector& X, int& Err, struct lssconfig=lssconfig());

//Additional parameter Y can be passed to function to compute an inner enclosure. For elements of the solution for which no
//inner enclousre could be computed, the according element of Y will have the value SignalingNaN.
//! Entry function of interval linear sytems solver
static inline void  ilss(const cxsc::imatrix& A, cxsc::imatrix& B, cxsc::imatrix& X, cxsc::imatrix& Y, int& Err, struct lssconfig=lssconfig());
//! Entry function of interval linear sytems solver
static inline void  ilss(const cxsc::imatrix& A, cxsc::ivector& B, cxsc::ivector& X, cxsc::ivector& Y, int& Err, struct lssconfig=lssconfig());
} //namespace cxsc

#define _CXSC_ILSS_INCLUDE
#include <fastlss.hpp> 

#endif
 
