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

/* CVS $Id: ddf_ari.hpp,v 1.16 2014/01/30 17:49:26 cxsc Exp $ */

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
// File: ddf_ari (header)
// Purpose: Definition of an interval differentiation arithmetic which allows
//    function evaluation with automatic differentiation up to second order.
// Type:
//    ddf_FctPtr              : pointer for a function of type 'DerivType'
// Class DerivType:
//    DerivType()             : constructors
//    operators +, -, *, /    : operators of diff. arithmetic
//    operator =              : assignment operator
//    DerivConst()
//    DerivVar()              : to define derivative constants/variables
//    fValue()
//    dfValue()
//    ddfValue()              : to get function and derivative values
//    sqr(), sqrt(), power(),
//    exp(), sin(), cos(), ...: elementary functions of diff. arithmetic
//    fEval()                 : to compute function value only
//    dfEval()                : to compute function and first derivative
//                              value
//    ddfEval()               : to compute function, first, and second
//                              derivative value
//----------------------------------------------------------------------------
#ifndef __DDF_ARI_HPP
#define __DDF_ARI_HPP

#include <interval.hpp>     // Interval arithmetic

using namespace cxsc;
using namespace std;

class DerivType;

typedef DerivType (*ddf_FctPtr)(const DerivType&);

class DerivType {
  private:
    interval f, df, ddf;

  public:
    DerivType ( );
    DerivType ( const interval&, const interval&, const interval& );
    DerivType ( const DerivType& );

    DerivType& operator= ( const DerivType& );

    friend DerivType DerivConst ( const real& );
    friend DerivType DerivConst ( const interval& );
    friend DerivType DerivVar   ( const real& );
    friend DerivType DerivVar   ( const interval& );

    friend inline const interval fValue   ( const DerivType& );  
    friend inline const interval dfValue  ( const DerivType& );  
    friend inline const interval ddfValue ( const DerivType& );  

    friend DerivType operator+ ( const DerivType& );
    friend DerivType operator- ( const DerivType& );

    friend DerivType operator+ ( const DerivType&, const DerivType& );
    friend DerivType operator- ( const DerivType&, const DerivType& );
    friend DerivType operator* ( const DerivType&, const DerivType& );
    friend DerivType operator/ ( const DerivType&, const DerivType& );

    friend DerivType operator+ ( const DerivType&, const interval& );
    friend DerivType operator- ( const DerivType&, const interval& );
    friend DerivType operator/ ( const DerivType&, const interval& );
    friend DerivType operator* ( const DerivType&, const interval& );

    friend DerivType operator+ ( const interval&, const DerivType& );
    friend DerivType operator- ( const interval&, const DerivType& );
    friend DerivType operator* ( const interval&, const DerivType& );
    friend DerivType operator/ ( const interval&, const DerivType& );

    friend DerivType operator+ ( const DerivType&, const real& );
    friend DerivType operator- ( const DerivType&, const real& );
    friend DerivType operator* ( const DerivType&, const real& );
    friend DerivType operator/ ( const DerivType&, const real& );

    friend DerivType operator+ ( const real&, const DerivType& );
    friend DerivType operator- ( const real&, const DerivType& );
    friend DerivType operator* ( const real&, const DerivType& );
    friend DerivType operator/ ( const real&, const DerivType& );

    friend DerivType sqr   ( const DerivType& );
    friend DerivType power ( const DerivType&, int );
    friend DerivType sqrt  ( const DerivType& );
    friend DerivType exp   ( const DerivType& );
    friend DerivType ln    ( const DerivType& );

    friend DerivType sin    ( const DerivType& );
    friend DerivType cos    ( const DerivType& );
    friend DerivType tan    ( const DerivType& );
    friend DerivType cot    ( const DerivType& );
    friend DerivType asin   ( const DerivType& );
    friend DerivType acos   ( const DerivType& );
    friend DerivType atan   ( const DerivType& );
    friend DerivType acot   ( const DerivType& );

    friend DerivType sinh   ( const DerivType& );
    friend DerivType cosh   ( const DerivType& );
    friend DerivType tanh   ( const DerivType& );
    friend DerivType coth   ( const DerivType& );
    friend DerivType asinh  ( const DerivType& );
    friend DerivType acosh  ( const DerivType& );
    friend DerivType atanh  ( const DerivType& );
    friend DerivType acoth  ( const DerivType& );

    friend void fEval  ( ddf_FctPtr f, interval, interval&);
    friend void dfEval ( ddf_FctPtr f, interval, interval&, interval& );
    friend void ddfEval( ddf_FctPtr f, interval,  interval&, interval&,
                                                             interval& );
};

DerivType sqr    ( const DerivType& );
DerivType power  ( const DerivType&, int );
DerivType sqrt   ( const DerivType& );
DerivType exp    ( const DerivType& );
DerivType ln     ( const DerivType& );
DerivType sin    ( const DerivType& );
DerivType cos    ( const DerivType& );
DerivType tan    ( const DerivType& );
DerivType cot    ( const DerivType& );
DerivType asin   ( const DerivType& );
DerivType acos   ( const DerivType& );
DerivType atan   ( const DerivType& );
DerivType acot   ( const DerivType& );
DerivType sinh   ( const DerivType& );
DerivType cosh   ( const DerivType& );
DerivType tanh   ( const DerivType& );
DerivType coth   ( const DerivType& );
DerivType asinh  ( const DerivType& );
DerivType acosh  ( const DerivType& );
DerivType atanh  ( const DerivType& );
DerivType acoth  ( const DerivType& );
void fEval  ( ddf_FctPtr f, interval, interval&);
void dfEval ( ddf_FctPtr f, interval, interval&, interval& );
void ddfEval( ddf_FctPtr f, interval,  interval&, interval&, interval& );

DerivType DerivConst ( const real& );
DerivType DerivConst ( const interval& );
DerivType DerivVar   ( const real& );
DerivType DerivVar   ( const interval& );

#include "ddf_ari.inl"
#endif




