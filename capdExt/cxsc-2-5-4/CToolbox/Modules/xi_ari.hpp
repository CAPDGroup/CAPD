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

/* CVS $Id: xi_ari.hpp,v 1.15 2014/01/30 17:49:27 cxsc Exp $ */

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
// File: xi_ari (header)
// Purpose: Definition of an extended interval arithmetic which allows the
//    division by an interval containing zero.
// Global type:
//    KindType         : component type of extended intervals
// Global function:
//    EmptyIntval()    : empty set as irregular interval
// Class xinterval:
//    operators %, -, &: operators of extended interval arithmetic
//    operator =       : assignment operator
// Updates:
//    04.03.1996 : Modification of extended interval operations and kind type
//----------------------------------------------------------------------------
#ifndef __XI_ARI_HPP
#define __XI_ARI_HPP

#include <interval.hpp>     // Interval arithmetic
#include <ivector.hpp>      // Interval vector arithmetic

using namespace cxsc;
using namespace std;

typedef enum { Finite, PlusInfty, MinusInfty, Double, Empty } KindType;

extern interval EmptyIntval ( );   // Irregular (empty) interval

class xinterval {                  // Extended intervals according
  private:                         // to the above definition
    KindType  kind;                //-----------------------------
    real      inf, sup;

  public:
    xinterval ( );
    xinterval ( const KindType&, const real&, const real& );
    xinterval ( const xinterval& );

    xinterval& operator= ( const xinterval& );

    friend xinterval operator- ( const real&, const xinterval& );
    friend xinterval operator% ( const interval&, const interval& );
    friend ivector   operator& ( const interval&, const xinterval& );
};

xinterval operator- ( const real&, const xinterval& );
xinterval operator% ( const interval&, const interval& );
ivector   operator& ( const interval&, const xinterval& );
#endif
