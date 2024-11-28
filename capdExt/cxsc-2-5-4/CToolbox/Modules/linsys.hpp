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

/* CVS $Id: linsys.hpp,v 1.14 2014/01/30 17:49:26 cxsc Exp $ */

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
// File: linsys (header)
// Purpose: Computation of a verified solution of a square linear system of
//    equations A*x = b with full real matrix A and real right-hand side b.
// Global functions:
//    LinSolve()      : to get a verified enclosure of the solution (two
//                      versions)
//    LinSolveErrMsg(): to get an error message text
//----------------------------------------------------------------------------
#ifndef __LINSYS_HPP
#define __LINSYS_HPP

#include <rmatrix.hpp>     // Real matrix/vector arithmetic
#include <ivector.hpp>     // Interval vector arithmetic

using namespace cxsc;
using namespace std;

extern char* LinSolveErrMsg ( int );
extern void  LinSolve ( const rmatrix&, const rvector&, ivector&, int& );
extern void  LinSolve ( const rmatrix&, const rvector&, ivector&, real&, int& );
#endif




