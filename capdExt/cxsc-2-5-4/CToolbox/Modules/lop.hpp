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

/* CVS $Id: lop.hpp,v 1.15 2014/01/30 17:49:27 cxsc Exp $ */

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
// File: lop (header)
// Purpose: Determine enclosures for the optimal value 'z_opt', for the set
//    of optimal basic index sets 'V_opt', and for the set of solution
//    vectors 'X_opt' for a linear programming problem P = (A,b,c) given in
//    the standard form:
//               ( z = c^t * x = max! )
//         (LP)  (     A * x = b      )
//               (      x >= 0        )
//    with an initial optimal basic index set.
// Global functions:
//    LinOpt()      : determines the enclosures of z_opt, V_opt, and X_opt
//                    for LP P = (A,b,c)
//    LinOptErrMsg(): delivers an error message text
//----------------------------------------------------------------------------
#ifndef __LOP_HPP
#define __LOP_HPP

#include <intmatrix.hpp>    // Integer matrix/vector types
#include <rmatrix.hpp>      // Real matrix/vector arithmetic
#include <imatrix.hpp>      // Interval matrix/vector arithmetic
#include <mv_util.hpp>      // Real matrix/vector utility functions
#include <mvi_util.hpp>     // Interval matrix/vector utility functions

using namespace cxsc;
using namespace std;

extern char* LinOptErrMsg ( int );
extern void  LinOpt ( rmatrix, rvector, rvector, intvector,
                      interval&, intmatrix&, imatrix&, int&, int& );
#endif




