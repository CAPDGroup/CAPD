//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file IFixStep.cpp
///
/// @author Daniel Wilczak
//
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) CAPD group
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#include "checkGetStep.h"
#include "capd/dynsys/FadOdeSolver.hpp"
#include "capd/dynsys/FadMap.h"

BOOST_AUTO_TEST_SUITE(IFixStepSuite)

template<class S>
S getInstance(IVector x){
  return S(x);
}

template<>
CnRect2Set getInstance(IVector x){
  return CnRect2Set(x,3);
}

template<>
CnMultiMatrixRect2Set getInstance(IVector x){
  return CnMultiMatrixRect2Set(x,3);
}

template<class S, class TM>
void dotest(TM& tm){
  IVector x(3);
  x[0] = 1.;
  x[1] = 2.;
  x[2] = 1.;

  S s = getInstance<S>(x);
  const double fixedTimeStep = 3./256;
  tm.getSolver().setStep(fixedTimeStep);

  checkFixedTime(tm,s,fixedTimeStep);

  tm(s.getCurrentTime(),s);
  BOOST_REQUIRE(tm.getSolver().getStep()==0.0);

  s.setCurrentTime(0.0);
  checkFixedTime(tm,s,fixedTimeStep);

  BOOST_REQUIRE(tm.getSolver().getStep()!=0.0);
  tm.getSolver().setStep(1.0);
  BOOST_REQUIRE(tm.getSolver().getStep()!=1.0);
}

BOOST_AUTO_TEST_CASE(C0C1_IFixStepTest)
{
  IMap f("var:x,y,z;fun:-y,x,z-x;");
  IOdeSolver solver(f,10);
  ITimeMap tm(solver);
  BOOST_REQUIRE(tm.getSolver().getStep()==0.0);

  dotest<C0Rect2Set>(tm);
  dotest<C0Pped2Set>(tm);
  dotest<C0TripletonSet>(tm);
  dotest<C0HORect2Set>(tm);
  dotest<C0HOTripletonSet>(tm);

  dotest<C1Rect2Set>(tm);
  dotest<C1Pped2Set>(tm);
  dotest<C1HORect2Set>(tm);
  dotest<C1HOPped2Set>(tm);
}

BOOST_AUTO_TEST_CASE(C2_IFixStepTest)
{
  IMap f("var:x,y,z;fun:-y,x,z-x;");
  IC2OdeSolver solver(f,10);
  IC2TimeMap tm(solver);
  BOOST_REQUIRE(tm.getSolver().getStep()==0.0);

  dotest<C0Rect2Set>(tm);
  dotest<C0Pped2Set>(tm);
  dotest<C0TripletonSet>(tm);

  dotest<C1Rect2Set>(tm);
  dotest<C1Pped2Set>(tm);

  dotest<C2Rect2Set>(tm);
}

BOOST_AUTO_TEST_CASE(Cn_IFixStepTest)
{
  IMap f("var:x,y,z;fun:-y,x,z-x;",3);
  ICnOdeSolver solver(f,10);
  ICnTimeMap tm(solver);
  BOOST_REQUIRE(tm.getSolver().getStep()==0.0);

  dotest<C0Rect2Set>(tm);
  dotest<C0Pped2Set>(tm);
  dotest<C0TripletonSet>(tm);

  dotest<C1Rect2Set>(tm);
  dotest<C1Pped2Set>(tm);

  dotest<C2Rect2Set>(tm);

  dotest<CnRect2Set>(tm);
  dotest<CnMultiMatrixRect2Set>(tm);
}

BOOST_AUTO_TEST_CASE(Fad_IFixStepTest)
{
  typedef capd::dynsys::LorenzFadMap<interval,0> TestMap;
  typedef capd::dynsys::FadOdeSolver<TestMap> IFadSolver;
  TestMap f(10.,28,8./3.);
  IFadSolver solver(f,10);
  capd::poincare::TimeMap<IFadSolver> tm(solver);
  BOOST_REQUIRE(tm.getSolver().getStep()==0.0);

  dotest<C0Rect2Set>(tm);
  dotest<C0Pped2Set>(tm);
  dotest<C0TripletonSet>(tm);

  dotest<C1Rect2Set>(tm);
  dotest<C1Pped2Set>(tm);
}

BOOST_AUTO_TEST_SUITE_END()
