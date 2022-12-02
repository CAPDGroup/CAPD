//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file AddMaskTest.cpp
///
/// @author Daniel Wilczak
//
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) CAPD group
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "capd/capdlib.h"
#include "capd/map/Map.hpp"
#include "capd/diffAlgebra/Jet.hpp"
#include "capd/vectalg/Matrix.hpp"
using namespace capd;

typedef capd::vectalg::Matrix<bool,0,0> BMatrix;
typedef capd::diffAlgebra::Jet<BMatrix,0> BJet;

void doTest(DMap& f, const BJet& mask, const DJet& expected);
DJet getPolyExpected();
DJet getIdExpected();
DJet getSwapExpected();
extern const char formula[];

BOOST_AUTO_TEST_SUITE(AddMaskSuite)

// directional derivatives -----------------------------------------------------------

void addDirectionalDerivativesTest(DMap& f, DJet expected){
  BJet m(2,3,3);

  // dx^3, dy^2, dz
  Multiindex m1[] = {{3,0,0},{0,2,0},{0,0,1}};
  f.setMask(m1,m1+1);
  m(0,0,0,0) = m(0,0,0) = m(0,0) = true;
  m(1,0,0,0) = m(1,0,0) = m(1,0) = true;
  m(0) = m(1) = true;
  doTest(f,m,expected);

  f.addMultiindexToMask(m1[1]);
  m(0,1,1) = m(0,1) = true;
  m(1,1,1) = m(1,1) = true;
  doTest(f,m,expected);

  f.addMultiindexToMask(m1[2]);
  m(0,2) = true;
  m(1,2) = true;
  doTest(f,m,expected);
}

BOOST_AUTO_TEST_CASE(xaddDirPolyTest)
{
  DMap f(formula,3);
  addDirectionalDerivativesTest(f,getPolyExpected());
}

BOOST_AUTO_TEST_CASE(xaddDirIdTest)
{
  DMap f("var:x,y,z;fun:x,y;",3);
  addDirectionalDerivativesTest(f,getIdExpected());
}

BOOST_AUTO_TEST_CASE(xaddDirSwapTest)
{
  DMap f("var:x,y,z;fun:y,x;",3);
  addDirectionalDerivativesTest(f,getSwapExpected());
}

// mixed derivatives -----------------------------------------------------------

void addMixedDerivativesTest(DMap& f, DJet expected){
  BJet m(2,3,3);

  // dx^2dy, dxdydz, dxdz^2
  Multiindex m1[] = {{2,1,0},{1,1,1},{1,0,2}};
  f.setMask(m1,m1+1);
  m(0,0,0,1) = m(0,0,0) = m(0,0,1) = m(0,0) = m(0,1) = true;
  m(1,0,0,1) = m(1,0,0) = m(1,0,1) = m(1,0) = m(1,1) = true;
  m(0) = m(1) = true;
  doTest(f,m,expected);

  f.addMultiindexToMask(m1[1]);
  m(0,0,1,2) = m(0,0,2) = m(0,1,2) = m(0,2) = true;
  m(1,0,1,2) = m(1,0,2) = m(1,1,2) = m(1,2) = true;
  doTest(f,m,expected);

  f.addMultiindexToMask(m1[2]);
  m(0,0,2,2) = m(0,2,2) = true;
  m(1,0,2,2) = m(1,2,2) = true;
  doTest(f,m,expected);
}

BOOST_AUTO_TEST_CASE(xaddMixPolyTest)
{
  DMap f(formula,3);
  addMixedDerivativesTest(f,getPolyExpected());
}

BOOST_AUTO_TEST_CASE(xaddMixIdTest)
{
  DMap f("var:x,y,z;fun:x,y;",3);
  addMixedDerivativesTest(f,getIdExpected());
}

BOOST_AUTO_TEST_CASE(xaddMixSwapTest)
{
  DMap f("var:x,y,z;fun:y,x;",3);
  addMixedDerivativesTest(f,getSwapExpected());
}

BOOST_AUTO_TEST_SUITE_END()
