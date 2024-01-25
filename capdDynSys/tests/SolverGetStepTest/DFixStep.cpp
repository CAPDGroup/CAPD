//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file DFixStep.cpp
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
#include "capd/dynsys/BasicFadOdeSolver.hpp"
#include "capd/dynsys/FadMap.h"

BOOST_AUTO_TEST_SUITE(DFixStepSuite)

template<class TM>
void dotest(TM& tm){
  DVector x(3);
  x[0] = 1.;
  x[1] = 2.;
  x[2] = 1.;
  const double fixedTimeStep = 0.3;
  tm.getSolver().setStep(fixedTimeStep);

  BOOST_REQUIRE(tm.getSolver().getStep()==0.0);
  checkFixedTime(tm,x,fixedTimeStep);
  checkFixedTime(tm,x,fixedTimeStep);

  BOOST_REQUIRE(tm.getSolver().getStep()!=0.0);
  tm.getSolver().setStep(1.0);
  BOOST_REQUIRE(tm.getSolver().getStep()!=1.0);
}

BOOST_AUTO_TEST_CASE(C0C1_DFixStepTest)
{
  DMap f("var:x,y,z;fun:-y,x,z-x;");
  DOdeSolver solver(f,10);
  DTimeMap tm(solver);
  dotest(tm);
}


BOOST_AUTO_TEST_CASE(C2_DFixStepTest)
{
  DMap f("var:x,y,z;fun:-y,x,z-x;");
  DC2OdeSolver solver(f,10);
  DC2TimeMap tm(solver);
  dotest(tm);
}

BOOST_AUTO_TEST_CASE(Cn_DFixStepTest)
{
  DMap f("var:x,y,z;fun:-y,x,z-x;");
  DCnOdeSolver solver(f,10);
  DCnTimeMap tm(solver);
  dotest(tm);
}

BOOST_AUTO_TEST_CASE(Fad_DFixStepTest)
{
  typedef capd::dynsys::LorenzFadMap<double,0> TestMap;
  typedef capd::dynsys::BasicFadOdeSolver<TestMap> DFadSolver;
  TestMap f(10.,28,8./3.);
  DFadSolver solver(f,10);
  capd::poincare::TimeMap<DFadSolver> tm(solver);
  dotest(tm);
}

BOOST_AUTO_TEST_SUITE_END()
