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

/* CVS $Id: xi_ari.cpp,v 1.16 2014/01/30 18:13:37 cxsc Exp $ */

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
// File: xi_ari (implementation)
// Purpose: Definition of an extended interval arithmetic which allows the
//    division by an interval containing zero.
// Method: Overloading of operators for arithmetic and lattice operations
//    of data type 'xinterval'.
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
#include <i_util.hpp>     // Interval utility functions
#include <xi_ari.hpp>

using namespace cxsc;
using namespace std;

//----------------------------------------------------------------------------
// An extended interval 'x', represented by the type 'xinterval', is defined
// according to the following rules (a <= b):
//
// x = [a, b]              :  x.kind = Finite,     x.inf = a, x.sup = b
// x = [a, +oo]            :  x.kind = PlusInfty,  x.inf = a, x.sup undefined
// x = [-oo, a]            :  x.kind = MinusInfty, x.sup = a, x.inf undefined
// x = [-oo, a] v [b, +oo] :  x.kind = Double,     x.inf = b, x.sup = a
// x = [-oo, +oo]          :  x.kind = Double,     x.inf = a, x.sup = a
// x = [/]                 :  x.kind = Empty,      x.inf  and x.sup undefined
//
// In this definition, 'v' stands for the set union and 'oo' for infinity.
//----------------------------------------------------------------------------
interval EmptyIntval ( )                         // Irregular (empty) interval
{                                                //---------------------------
  interval x;
  x = _unchecked_interval(999999999.0,-999999999.0);
  return x;
}

//----------------------------------------------------------------------------
// Constructors and assignment operator
//----------------------------------------------------------------------------
xinterval::xinterval ( )
{
  kind = Finite;
  inf  = 0.0;
  sup  = 0.0;
}

xinterval::xinterval ( const KindType& k, const real& i, const real& s )
{
  kind = k;
  inf  = i;
  sup  = s;
}

xinterval::xinterval ( const xinterval& a )
{
  kind = a.kind;
  inf  = a.inf;
  sup  = a.sup;
}

xinterval& xinterval::operator= ( const xinterval& a )
{
  kind = a.kind;
  inf  = a.inf;
  sup  = a.sup;
  return *this;
}

//----------------------------------------------------------------------------
// Extended interval division 'A / B' where 0 in 'B' is allowed.
//----------------------------------------------------------------------------
xinterval operator% ( const interval& A, const interval& B )
{
  interval  c;
  xinterval Q;

  if ( in(0.0, B) ) {
    if ( in(0.0, A) ) {
      Q.kind = Double;                    // Q = [-oo,+oo] = [-oo,0] v [0,+oo]
      Q.sup  = 0.0;                       //----------------------------------
      Q.inf  = 0.0;
    }
    else if ( B == 0.0 ) {                                          // Q = [/]
      Q.kind = PlusInfty;                                           //--------
      Q.inf  = divd(Sup(A),Inf(B));
    }
    else if ( (Sup(A) < 0.0) && (Sup(B) == 0.0) ) {         // Q = [Q.inf,+oo]
      Q.kind = PlusInfty;                                   //----------------
      Q.inf  = divd(Sup(A),Inf(B));
    }
    else if ( (Sup(A) < 0.0) && (Inf(B) < 0.0) && (Sup(B) > 0.0) ) {
      Q.kind = Double;                        // Q = [-oo,Q.sup] v [Q.inf,+oo]
      Q.sup  = divu(Sup(A),Sup(B));           //------------------------------
      Q.inf  = divd(Sup(A),Inf(B));
    }
    else if ( (Sup(A) < 0.0) && (Inf(B) == 0.0) ) {         // Q = [-oo,Q.sup]
      Q.kind = MinusInfty;                                  //----------------
      Q.sup  = divu(Sup(A),Sup(B));
    }
    else if ( (Inf(A) > 0.0) && (Sup(B) == 0.0) ) {         // Q = [-oo,Q.sup]
      Q.kind = MinusInfty;                                  //----------------
      Q.sup  = divu(Inf(A),Inf(B));
    }
    else if ( (Inf(A) > 0.0) && (Inf(B) < 0.0) && (Sup(B) > 0.0) ) {
      Q.kind = Double;                        // Q = [-oo,Q.sup] v [Q.inf,+oo]
      Q.sup  = divu(Inf(A),Inf(B));           //------------------------------
      Q.inf  = divd(Inf(A),Sup(B));
    }
    else { // if ( (Inf(A) > 0.0) && (Inf(B) == 0.0) )
      Q.kind = PlusInfty;                                   // Q = [Q.inf,+oo]
      Q.inf  = divd(Inf(A),Sup(B));                         //----------------
    }
  } // in(0.0,B)
  else {  // !in(0.0,B)
    c = A / B;                                            // Q = [C.inf,C.sup]
    Q.kind = Finite;                                      //------------------
    Q.inf  = Inf(c);
    Q.sup  = Sup(c);
  }

  return Q;
} // operator%

//----------------------------------------------------------------------------
// Subtraction of an extended interval 'B' from a real value 'a'.
//----------------------------------------------------------------------------
xinterval operator- ( const real& a, const xinterval& B )
{
  xinterval D;

  switch (B.kind) {
    case Finite     : D.kind = Finite;                    // D = [D.inf,D.sup]
                      D.inf  = subd(a,B.sup);             //------------------
                      D.sup  = subu(a,B.inf);
                      break;
    case PlusInfty  : D.kind = MinusInfty;                    // D = [inf,+oo]
                      D.sup  = subu(a,B.inf);                 //--------------
                      break;
    case MinusInfty : D.kind = PlusInfty;                     // D = [-oo,sup]
                      D.inf  = subd(a,B.sup);                 //--------------
                      break;
    case Double     : D.kind = Double;        // D = [-oo,D.sup] v [D.inf,+oo]
                      D.inf  = subd(a,B.sup); //------------------------------
                      D.sup  = subu(a,B.inf);
                      if (D.inf < D.sup) D.inf = D.sup;
                      break;
    case Empty      : D.kind = Empty;                               // D = [/]
                      D.inf  = subd(a,B.sup);                       //--------
                      break;
  } // switch
  return D;
}

//----------------------------------------------------------------------------
// Intersection of an interval 'X' and an extended interval 'Y'. The result
// is given as a pair (vector) of intervals, where one or both of them can
// be empty intervals.
//----------------------------------------------------------------------------
ivector operator& ( const interval& X, const xinterval& Y )
{
  interval H;
  ivector  IS(2);

  IS[1] = EmptyIntval();
  IS[2] = EmptyIntval();

  switch (Y.kind) {
    case Finite     : // [X.inf,X.sup] & [Y.inf,Y.sup]
                      //------------------------------
                      H = interval(Y.inf,Y.sup);
                      if ( !Disjoint(X,H) ) IS[1] = X & H;
                      break;
    case PlusInfty  : // [X.inf,X.sup] & [Y.inf,+oo]
                      //----------------------------
                      if (Sup(X) >= Y.inf) {
                        if (Inf(X) > Y.inf)
                          IS[1] = X;
                        else
                          IS[1] = interval(Y.inf,Sup(X));
		      }  
                      break;
    case MinusInfty : // [X.inf,X.sup] & [-oo,Y.sup]
                      //----------------------------
                      if (Y.sup >= Inf(X)) {
                        if (Sup(X)<Y.sup)
                          IS[1] = X;
                        else
                          IS[1] = interval(Inf(X),Y.sup);
		      }	  
                      break;
    case Double     : if ( (Inf(X) <= Y.sup) && (Y.inf <= Sup(X)) ) {
                        IS[1] = interval(Inf(X),Y.sup);    // X & [-oo,Y.sup]
                        IS[2] = interval(Y.inf,Sup(X));    // X & [Y.inf,+oo]
                      }
                      else if (Y.inf <= Sup(X)) // [X.inf,X.sup] & [Y.inf,+oo]
                        if (Inf(X) >= Y.inf)    //----------------------------
                          IS[1] = X;
                        else
                          IS[1] = interval(Y.inf,Sup(X));
                      else if (Inf(X) <= Y.sup) { // [X.inf,X.sup] & [-oo,Y.sup]
                        if (Sup(X) <= Y.sup)      //----------------------------
                          IS[1] = X;
                        else
                          IS[1] = interval(Inf(X),Y.sup);
		      }
                      break;
    case Empty      : break;                           // [X.inf,X.sup] ** [/]
  } // switch                                          //---------------------

  return IS;
} // operator&




