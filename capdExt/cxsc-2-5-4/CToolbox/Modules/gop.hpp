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

/* CVS $Id: gop.hpp,v 1.16 2014/01/30 17:49:26 cxsc Exp $ */

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
// File: gop (header)
// Purpose: Computing enclosures for all global minimizers and for the global
//    minimum value of a twice continuously differentiable multi-dimensional,
//    scalar valued function, assuming that the global minimum is a stationa-
//    ry point. If it is a boundary point of the search area with gradient of
//    the function being different from zero, the method fails in its form
//    presented here.
// Global functions:
//    AllGOp()      : computes enclosures for all zeros
//    AllGOpErrMsg(): delivers an error message text
//----------------------------------------------------------------------------
#ifndef __GOP_HPP
#define __GOP_HPP

#include <intvector.hpp>    // Integer vector type
#include <interval.hpp>     // Interval arithmetic
#include <hess_ari.hpp>     // Hessian differentiation arithmetic

// Additional header files to get access to interval standard functions and
// interval utility functions which often are needed for implementing user
// defined functions of type 'HessType'
//--------------------------------------------------------------------------
#include <i_util.hpp>       // Interval utility functions
#include <imath.hpp>        // Interval mathematical functions

using namespace cxsc;
using namespace std;

const int MaxCountGOp = 10000; // Maximum count of result components

extern char* AllGOpErrMsg ( int );
extern void  AllGOp ( HTscalar_FctPtr, ivector, real, imatrix&, intvector&,
                      int&, interval&, int&, int = MaxCountGOp );
#endif




