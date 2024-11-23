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

/* CVS $Id: nlfzero.hpp,v 1.14 2014/01/30 17:49:27 cxsc Exp $ */

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
// File: nlfzero (header)
// Purpose: Computing enclosures for all zeros of a continuously
//    differentiable one-dimensional, real-valued function.
// Global functions:
//    AllZeros()      : computes enclosures for all zeros
//    AllZerosErrMsg(): delivers an error message text
//----------------------------------------------------------------------------
#ifndef __NLFZERO_HPP
#define __NLFZERO_HPP

#include <intvector.hpp>     // Integer vector type
#include <ivector.hpp>       // Interval vector arithmetic
#include <mvi_util.hpp>      // Interval matrix/vector utility functions
#include <ddf_ari.hpp>       // Differentiation arithmetic

using namespace cxsc;
using namespace std;

const int MaxCount = 10000;  // Maximum count of result components

extern char* AllZerosErrMsg ( int );
extern void  AllZeros ( ddf_FctPtr, interval, real, ivector&, intvector&,
                        int&, int&, int = MaxCount );
#endif




