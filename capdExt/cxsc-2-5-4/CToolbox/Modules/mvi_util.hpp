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

/* CVS $Id: mvi_util.hpp,v 1.14 2014/01/30 17:49:27 cxsc Exp $ */

//============================================================================
//
//                              Program/Module
//                                   from
//                 C++ TOOLBOX FOR VERIFIED COMPUTING I
//                         Basic Numerical Problems
//
//      Copyright (c) 1995   Rolf Hammer, Matthias Hocks, Dietmar Ratz
//
// For details on theory, algorithms, and programs, see the book
//
//  R. Hammer, M. Hocks, U. Kulisch, D. Ratz:  C++ Toolbox for
//  Verified Computing I - Basic Numerical Problems. Springer-Verlag,
//  Heidelberg, New York, 1995.
//
//============================================================================
//----------------------------------------------------------------------------
// File: mvi_util (header)
// Purpose: Utilities of type 'ivector' and 'imatrix'.
// Global functions:
//    in()        : contained-in-the-interior relation for two interval
//                  vectors and for an integer in an interval vector
//    Blow()      : epsilon inflation for an interval vector
//    Disjoint()  : test for disjointedness of two interval vectors
//    Id()        : identity interval matrix
//    Zero()      : check for a zero interval matrix
//    mid()       : midpoint vector and matrix
//    MaxRelDiam(): to get the maximum of the relative diameters of the
//                  components of an interval vector.
//    UlpAcc()    : to check whether all components of an interval vector
//                  have width less than a certain number of ulps
//                  (ulp = unit in the last place of the mantissa).
//    VecLen()    : length of an interval vector
//    RowLen()    : length of the rows of an interval matrix
//    ColLen()    : length of the columns of an interval matrix
//    DoubleSize(): for doubling the size of an interval vector or matrix
//----------------------------------------------------------------------------
#ifndef __MVI_UTIL_HPP
#define __MVI_UTIL_HPP

#include <i_util.hpp>      // Interval utility functions
#include <rmatrix.hpp>     // Real matrix/vector arithmetic
#include <imatrix.hpp>     // Interval matrix/vector arithmetic
#include <mv_util.hpp>     // Utility functions for real matrices, real
 
//using namespace cxsc;
//using namespace std;
                          // vectors, and integer vectors

//extern int      in         ( const ivector&, const ivector& );
//extern int      in         ( int,      ivector& );
//extern ivector  Blow       ( const ivector&, real );
//extern int      Disjoint   ( ivector&, ivector& );
//extern imatrix  Id         ( imatrix& );
//extern int      Zero       ( ivector& );
//extern rmatrix  mid        ( imatrix& );
//extern rvector  mid        ( ivector& );
//extern real     MaxRelDiam ( const ivector& );
//extern real     MaxRelDiam ( const ivector_slice& );
//extern real     MaxRelDiam ( const imatrix_subv& );
//extern int      UlpAcc     ( ivector&, int );
//extern int      VecLen     ( ivector& );
//extern int      RowLen     ( imatrix& );
//extern int      ColLen     ( imatrix& );
//extern void     DoubleSize ( ivector& );
//extern void     DoubleSize ( imatrix& );
#endif




