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

/* CVS $Id: gop1.hpp,v 1.15 2014/01/30 17:49:26 cxsc Exp $ */

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
// File: gop1 (header)
// Purpose: Computing enclosures for all global minimizers and for the global
//    minimum value of a twice continuously differentiable one-dimensional,
//    scalar valued function.
// Global functions:
//    AllGOp1()      : computes enclosures for all global optimizers
//    AllGOp1ErrMsg(): delivers an error message text
//----------------------------------------------------------------------------
#ifndef __GOP1_HPP
#define __GOP1_HPP

#include <interval.hpp>      // Interval arithmetic
#include <intvector.hpp>     // Integer vector type
#include <ivector.hpp>       // Interval vector arithmetic
#include <ddf_ari.hpp>       // Differentiation arithmetic

using namespace cxsc;
using namespace std;

const int MaxCountGOp1 = 10000;  // Maximum count of result components

extern char* AllGOp1ErrMsg ( int );
extern void  AllGOp1 ( ddf_FctPtr, interval, real, ivector&, intvector&,
                       int&, interval&, int&, int = MaxCountGOp1 );
#endif




