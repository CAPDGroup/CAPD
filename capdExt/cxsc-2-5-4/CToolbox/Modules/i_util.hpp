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

/* CVS $Id: i_util.hpp,v 1.14 2014/01/30 17:49:26 cxsc Exp $ */

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
// File: i_util (header)
// Purpose: Utilities of type 'interval'
// Global functions:
//    in()      : contained-in relation for a real
//    in()      : contained-in-the-interior relation for two intervals
//    rnd()     : to round a dotprecision argument to an interval
//    Blow()    : epsilon inflation
//    Disjoint(): test for disjointedness of two intervals
//    AbsMin()  : smallest absolute value of an interval
//    AbsMax()  : greatest absolute value of an interval
//    RelDiam() : relative diameter of an interval
//    UlpAcc()  : to check whether the width of an interval is less than a
//                certain number of ulps (ulp = unit in the last place of
//                the mantissa)
//    Power()   : exponentiation by an integer for intervals
//----------------------------------------------------------------------------
#ifndef __I_UTIL_HPP
#define __I_UTIL_HPP

#include <interval.hpp>     // Interval arithmetic

// obsolete file

//using namespace cxsc;
//using namespace std;

//extern int      in         ( real, const interval& );
//extern int      in         ( const interval&, const interval& );
//extern void     rnd        ( dotprecision&, interval& );
//extern interval Blow       ( interval, real );
//extern int      Disjoint   ( interval, interval );
//extern real     AbsMin     ( interval );
//extern real     AbsMax     ( interval );
//extern real     RelDiam    ( interval );
//extern int      UlpAcc     ( interval, int );
//extern interval Power      ( interval, int );
//extern interval Pi         ( );
#endif




