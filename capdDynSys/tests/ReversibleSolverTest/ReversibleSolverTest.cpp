//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file ReversibleSolverTest.cpp
///
/// @author Daniel Wilczak
//
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) CAPD group
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#include "doTest.h"

using namespace capd;

/**
 * The Michelson system is reversible with respect to involution (x,y,z,t) -> (-x,y,-z,-t)
 * We check rigorous C^0-C^n solvers for preserving this property.
 */

BOOST_AUTO_TEST_SUITE(SolverSuite)

BOOST_AUTO_TEST_CASE(C0Test)
{
  IMap f="par:c;var:x,y,z;fun:y,z,1-y-0.5*x*x;";
  IOdeSolver solver(f,10);
  ITimeMap tm(solver);
  IVector x(0.,1.563,0.);

  doC0test<C0Rect2Set>(tm,x);
  doC0test<C0Pped2Set>(tm,x);
  doC0test<C0TripletonSet>(tm,x);
  doC0test<C0HORect2Set>(tm,x);
  doC0test<C0HOTripletonSet>(tm,x);
}

BOOST_AUTO_TEST_CASE(C1Test)
{
  IMap f="par:c;var:x,y,z;fun:y,z,1-y-0.5*x*x;";
  IOdeSolver solver(f,10);
  ITimeMap tm(solver);
  IVector x(0.,1.563,0.);

  doC1test<C1Rect2Set>(tm,x);
  doC1test<C1Pped2Set>(tm,x);
  doC1test<C1HORect2Set>(tm,x);
  doC1test<C1HOPped2Set>(tm,x);
}

BOOST_AUTO_TEST_CASE(C2Test)
{
  IMap f("par:c;var:x,y,z;fun:y,z,1-y-0.5*x*x;",2);
  IC2OdeSolver solver(f,10);
  IC2TimeMap tm(solver);
  IVector x(0.,1.563,0.);

  doC0test<C0Rect2Set>(tm,x);
  doC0test<C0Pped2Set>(tm,x);
  doC0test<C0TripletonSet>(tm,x);

  doC1test<C1Rect2Set>(tm,x);
  doC1test<C1Pped2Set>(tm,x);

  doC2test<C2Rect2Set>(tm,x);
}

BOOST_AUTO_TEST_CASE(CnTest)
{
  IMap f("par:c;var:x,y,z;fun:y,z,1-y-0.5*x*x;",3);
  ICnOdeSolver solver(f,10);
  ICnTimeMap tm(solver);
  IVector x(0.,1.563,0.);

  doC0test<C0Rect2Set>(tm,x);
  doC0test<C0Pped2Set>(tm,x);
  doC0test<C0TripletonSet>(tm,x);

  doC1test<C1Rect2Set>(tm,x);
  doC1test<C1Pped2Set>(tm,x);

  doC2test<C2Rect2Set>(tm,x);

  doCntest<CnRect2Set>(tm,x);
  doCntest<CnMultiMatrixRect2Set>(tm,x);
}

BOOST_AUTO_TEST_SUITE_END()
