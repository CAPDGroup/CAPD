//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file SetMaskTest.cpp
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

BOOST_AUTO_TEST_SUITE(SetMaskSuite)

// directional derivatives -----------------------------------------------------------

void directionalDerivativesTest(DMap& f, DJet expected){
  BJet m(2,3,3);

  // dx^3, dy^2, dz
  Multiindex m1[] = {{3,0,0},{0,2,0},{0,0,1}};
  f.setMask(m1,m1+3);
  m(0,0,0,0) = m(0,0,0) = m(0,0) = true;
  m(1,0,0,0) = m(1,0,0) = m(1,0) = true;
  m(0,1,1) = m(0,1) = true;
  m(1,1,1) = m(1,1) = true;
  m(0,2) = true;
  m(1,2) = true;
  m(0) = m(1) = true;

  doTest(f,m,expected);
}

BOOST_AUTO_TEST_CASE(xdirPolyTest)
{
  DMap f(formula,3);
  directionalDerivativesTest(f,getPolyExpected());
}

BOOST_AUTO_TEST_CASE(xdirIdTest)
{
  DMap f("var:x,y,z;fun:x,y;",3);
  directionalDerivativesTest(f,getIdExpected());
}

BOOST_AUTO_TEST_CASE(xdirSwapTest)
{
  DMap f("var:x,y,z;fun:y,x;",3);
  directionalDerivativesTest(f,getSwapExpected());
}

// mixed derivatives -----------------------------------------------------------

void mixedDerivativesTest(DMap& f, DJet expected){
  BJet m(2,3,3);

  // dx^2dy, dxdydz, dxdz^2
  Multiindex m1[] = {{2,1,0},{1,1,1},{1,0,2}};
  f.setMask(m1,m1+3);
  m(0,0,0,1) = m(0,0,0) = m(0,0,1) = m(0,0) = m(0,1) = true;
  m(1,0,0,1) = m(1,0,0) = m(1,0,1) = m(1,0) = m(1,1) = true;
  m(0,0,1,2) = m(0,0,2) = m(0,1,2) = m(0,2) = true;
  m(1,0,1,2) = m(1,0,2) = m(1,1,2) = m(1,2) = true;
  m(0,0,2,2) = m(0,2,2) = true;
  m(1,0,2,2) = m(1,2,2) = true;
  m(0) = m(1) = true;

  doTest(f,m,expected);
}

BOOST_AUTO_TEST_CASE(xmixPolyTest)
{
  DMap f(formula,3);
  mixedDerivativesTest(f,getPolyExpected());
}

BOOST_AUTO_TEST_CASE(xmixIdTest)
{
  DMap f("var:x,y,z;fun:x,y;",3);
  mixedDerivativesTest(f,getIdExpected());
}

BOOST_AUTO_TEST_CASE(xmixSwapTest)
{
  DMap f("var:x,y,z;fun:y,x;",3);
  mixedDerivativesTest(f,getSwapExpected());
}

// overwrite mask -----------------------------------------------------------

void overwriteMaskTest(DMap& f, DJet expected){
  BJet m(2,3,3);

  Multiindex tmp[] = {{2,1,0},{1,0,2},{1,1,0}};
  f.setMask(tmp,tmp+3);

  // dx^2dy, dxdydz, dxdz^2
  Multiindex m1[] = {{2,1,0},{1,1,1},{1,0,2}};
  f.setMask(m1,m1+3);
  m(0,0,0,1) = m(0,0,0) = m(0,0,1) = m(0,0) = m(0,1) = true;
  m(1,0,0,1) = m(1,0,0) = m(1,0,1) = m(1,0) = m(1,1) = true;
  m(0,0,1,2) = m(0,0,2) = m(0,1,2) = m(0,2) = true;
  m(1,0,1,2) = m(1,0,2) = m(1,1,2) = m(1,2) = true;
  m(0,0,2,2) = m(0,2,2) = true;
  m(1,0,2,2) = m(1,2,2) = true;
  m(0) = m(1) = true;

  doTest(f,m,expected);
}

BOOST_AUTO_TEST_CASE(xoverPolyTest)
{
  DMap f(formula,3);
  overwriteMaskTest(f,getPolyExpected());
}

BOOST_AUTO_TEST_CASE(xoverIdTest)
{
  DMap f("var:x,y,z;fun:x,y;",3);
  overwriteMaskTest(f,getIdExpected());
}

BOOST_AUTO_TEST_CASE(xoverSwapTest)
{
  DMap f("var:x,y,z;fun:y,x;",3);
  overwriteMaskTest(f,getSwapExpected());
}


BOOST_AUTO_TEST_SUITE_END()
