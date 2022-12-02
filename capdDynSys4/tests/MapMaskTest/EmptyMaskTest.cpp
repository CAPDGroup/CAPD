//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file EmptyMaskTest.cpp
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

using namespace capd;

DJet getPolyExpected();
DJet getIdExpected();
DJet getSwapExpected();
extern const DVector x;
extern const char formula[];

BOOST_AUTO_TEST_SUITE(EmptyMaskSuite)

// C^0 tests -----------------------------------------------------------

void c0Test(DMap& f, DVector x, DJet expected){
  BOOST_CHECK_EQUAL(f(x),DVector(expected));
  BOOST_CHECK_EQUAL(f(10,x),DVector(expected));
  f.setMask((Multiindex*)0,(Multiindex*)0);
  BOOST_CHECK_EQUAL(f(x),DVector(expected));
  BOOST_CHECK_EQUAL(f(20,x),DVector(expected));

  f.setOrder(10);
  BOOST_CHECK_EQUAL(f(x),DVector(expected));
  BOOST_CHECK_EQUAL(f(30,x),DVector(expected));

  f.resetMask();
  BOOST_CHECK_EQUAL(f(x),DVector(expected));
  BOOST_CHECK_EQUAL(f(30,x),DVector(expected));
}

BOOST_AUTO_TEST_CASE(xc0PolyTest)
{
  DMap f(formula,3);
  c0Test(f,x,getPolyExpected());
}

BOOST_AUTO_TEST_CASE(xc0IdTest)
{
  DMap f("var:x,y,z;fun:x,y;",3);
  c0Test(f,x,getIdExpected());
}

BOOST_AUTO_TEST_CASE(xc0SwapTest)
{
  DMap f("var:x,y,z;fun:y,x;",3);
  c0Test(f,x,getSwapExpected());
}

// C^1 tests -----------------------------------------------------------

void c1Test(DMap& f, DVector x, DJet expected){
  DMatrix df(2,3), zero(2,3);
  DVector y(2);

  BOOST_CHECK_EQUAL(f.derivative(x),DMatrix(expected));
  BOOST_CHECK_EQUAL(f.derivative(10,x),DMatrix(expected));
  BOOST_CHECK_EQUAL(f[x],DMatrix(expected));
  y=f(x,df);
  BOOST_CHECK_EQUAL(y,DVector(expected));
  BOOST_CHECK_EQUAL(df,DMatrix(expected));
  y=f(20,x,df);
  BOOST_CHECK_EQUAL(y,DVector(expected));
  BOOST_CHECK_EQUAL(df,DMatrix(expected));

  f.setMask((Multiindex*)0,(Multiindex*)0);
  BOOST_CHECK_EQUAL(f.derivative(x),zero);
  BOOST_CHECK_EQUAL(f.derivative(10,x),zero);
  BOOST_CHECK_EQUAL(f[x],zero);
  y=f(x,df);
  BOOST_CHECK_EQUAL(y,DVector(expected));
  BOOST_CHECK_EQUAL(df,zero);
  y=f(20,x,df);
  BOOST_CHECK_EQUAL(y,DVector(expected));
  BOOST_CHECK_EQUAL(df,zero);

  f.setOrder(10);
  BOOST_CHECK_EQUAL(f.derivative(x),zero);
  BOOST_CHECK_EQUAL(f.derivative(10,x),zero);
  BOOST_CHECK_EQUAL(f[x],zero);
  y=f(x,df);
  BOOST_CHECK_EQUAL(y,DVector(expected));
  BOOST_CHECK_EQUAL(df,zero);
  y=f(20,x,df);
  BOOST_CHECK_EQUAL(y,DVector(expected));
  BOOST_CHECK_EQUAL(df,zero);

  f.resetMask();
  BOOST_CHECK_EQUAL(f.derivative(x),DMatrix(expected));
  BOOST_CHECK_EQUAL(f.derivative(10,x),DMatrix(expected));
  BOOST_CHECK_EQUAL(f[x],DMatrix(expected));
  y=f(x,df);
  BOOST_CHECK_EQUAL(y,DVector(expected));
  BOOST_CHECK_EQUAL(df,DMatrix(expected));
  y=f(20,x,df);
  BOOST_CHECK_EQUAL(y,DVector(expected));
  BOOST_CHECK_EQUAL(df,DMatrix(expected));

  f.setMask((Multiindex*)0,(Multiindex*)0);
  BOOST_CHECK_EQUAL(f.derivative(x),zero);
  BOOST_CHECK_EQUAL(f.derivative(10,x),zero);
  BOOST_CHECK_EQUAL(f[x],zero);
  y=f(x,df);
  BOOST_CHECK_EQUAL(y,DVector(expected));
  BOOST_CHECK_EQUAL(df,zero);
  y=f(20,x,df);
  BOOST_CHECK_EQUAL(y,DVector(expected));
  BOOST_CHECK_EQUAL(df,zero);
}

BOOST_AUTO_TEST_CASE(xc1PolyTest)
{
  DMap f(formula,3);
  c1Test(f,x,getPolyExpected());
}

BOOST_AUTO_TEST_CASE(xc1IdTest)
{
  DMap f("var:x,y,z;fun:x,y;",3);
  c1Test(f,x,getIdExpected());
}

BOOST_AUTO_TEST_CASE(xc1SwapTest)
{
  DMap f("var:x,y,z;fun:y,x;",3);
  c1Test(f,x,getSwapExpected());
}


// C^2 tests -----------------------------------------------------------

void c2Test(DMap& f, DVector x, DJet expected){
  DMatrix df(2,3), zero(2,3);
  DHessian d2f(2,3), z(2,3);
  DVector y(2);

  y=f(x,df,d2f);
  BOOST_CHECK_EQUAL(y,DVector(expected));
  BOOST_CHECK_EQUAL(df,DMatrix(expected));
  BOOST_CHECK_EQUAL(d2f,DHessian(expected));
  y=f(20,x,df,d2f);
  BOOST_CHECK_EQUAL(y,DVector(expected));
  BOOST_CHECK_EQUAL(df,DMatrix(expected));
  BOOST_CHECK_EQUAL(d2f,DHessian(expected));

  f.setMask((Multiindex*)0,(Multiindex*)0);
  y=f(x,df,d2f);
  BOOST_CHECK_EQUAL(y,DVector(expected));
  BOOST_CHECK_EQUAL(df,zero);
  BOOST_CHECK_EQUAL(d2f,z);
  y=f(20,x,df,d2f);
  BOOST_CHECK_EQUAL(y,DVector(expected));
  BOOST_CHECK_EQUAL(df,zero);
  BOOST_CHECK_EQUAL(d2f,z);

  f.setOrder(10);
  y=f(x,df,d2f);
  BOOST_CHECK_EQUAL(y,DVector(expected));
  BOOST_CHECK_EQUAL(df,zero);
  BOOST_CHECK_EQUAL(d2f,z);
  y=f(20,x,df,d2f);
  BOOST_CHECK_EQUAL(y,DVector(expected));
  BOOST_CHECK_EQUAL(df,zero);
  BOOST_CHECK_EQUAL(d2f,z);

  f.resetMask();
  y=f(x,df,d2f);
  BOOST_CHECK_EQUAL(y,DVector(expected));
  BOOST_CHECK_EQUAL(df,DMatrix(expected));
  BOOST_CHECK_EQUAL(d2f,DHessian(expected));
  y=f(20,x,df,d2f);
  BOOST_CHECK_EQUAL(y,DVector(expected));
  BOOST_CHECK_EQUAL(df,DMatrix(expected));
  BOOST_CHECK_EQUAL(d2f,DHessian(expected));

  f.setMask((Multiindex*)0,(Multiindex*)0);
  y=f(x,df,d2f);
  BOOST_CHECK_EQUAL(y,DVector(expected));
  BOOST_CHECK_EQUAL(df,zero);
  BOOST_CHECK_EQUAL(d2f,z);
  y=f(20,x,df,d2f);
  BOOST_CHECK_EQUAL(y,DVector(expected));
  BOOST_CHECK_EQUAL(df,zero);
  BOOST_CHECK_EQUAL(d2f,z);
}

BOOST_AUTO_TEST_CASE(xc2PolyTest)
{
  DMap f(formula,3);
  c2Test(f,x,getPolyExpected());
}

BOOST_AUTO_TEST_CASE(xc2IdTest)
{
  DMap f("var:x,y,z;fun:x,y;",3);
  c2Test(f,x,getIdExpected());
}

BOOST_AUTO_TEST_CASE(xc2SwapTest)
{
  DMap f("var:x,y,z;fun:y,x;",3);
  c2Test(f,x,getSwapExpected());
}

// C^n tests -----------------------------------------------------------

void cnTest(DMap& f, DVector x, DJet expected){
  DJet i(3,3,3);
  i() = x;
  i(0,0) = i(1,1) = i(2,2) = 1;
  DJet df(2,3,3), zero(2,3,3);
  zero() = expected();
  DVector y(2);

  df=f(i);
  BOOST_CHECK_EQUAL(df,expected);
  df=f(10,i);
  BOOST_CHECK_EQUAL(df,expected);
  y=f(x,df);
  BOOST_CHECK_EQUAL(y,DVector(expected));
  BOOST_CHECK_EQUAL(df,expected);
  y=f(20,x,df);
  BOOST_CHECK_EQUAL(y,DVector(expected));
  BOOST_CHECK_EQUAL(df,expected);

  f.setMask((Multiindex*)0,(Multiindex*)0);
  df=f(i);
  BOOST_CHECK_EQUAL(df,zero);
  df=f(10,i);
  BOOST_CHECK_EQUAL(df,zero);
  y=f(x,df);
  BOOST_CHECK_EQUAL(y,DVector(expected));
  BOOST_CHECK_EQUAL(df,zero);
  y=f(20,x,df);
  BOOST_CHECK_EQUAL(y,DVector(expected));
  BOOST_CHECK_EQUAL(df,zero);

  f.setOrder(10);
  df=f(i);
  BOOST_CHECK_EQUAL(df,zero);
  df=f(10,i);
  BOOST_CHECK_EQUAL(df,zero);
  y=f(x,df);
  BOOST_CHECK_EQUAL(y,DVector(expected));
  BOOST_CHECK_EQUAL(df,zero);
  y=f(20,x,df);
  BOOST_CHECK_EQUAL(y,DVector(expected));
  BOOST_CHECK_EQUAL(df,zero);

  f.resetMask();
  df=f(i);
  BOOST_CHECK_EQUAL(df,expected);
  df=f(10,i);
  BOOST_CHECK_EQUAL(df,expected);
  y=f(x,df);
  BOOST_CHECK_EQUAL(y,DVector(expected));
  BOOST_CHECK_EQUAL(df,expected);
  y=f(20,x,df);
  BOOST_CHECK_EQUAL(y,DVector(expected));
  BOOST_CHECK_EQUAL(df,expected);

  f.setMask((Multiindex*)0,(Multiindex*)0);
  df=f(i);
  BOOST_CHECK_EQUAL(df,zero);
  df=f(10,i);
  BOOST_CHECK_EQUAL(df,zero);
  y=f(x,df);
  BOOST_CHECK_EQUAL(y,DVector(expected));
  BOOST_CHECK_EQUAL(df,zero);
  y=f(20,x,df);
  BOOST_CHECK_EQUAL(y,DVector(expected));
  BOOST_CHECK_EQUAL(df,zero);
}

BOOST_AUTO_TEST_CASE(xcnPolyTest)
{
  DMap f(formula,3);
  cnTest(f,x,getPolyExpected());
}

BOOST_AUTO_TEST_CASE(xcnIdTest)
{
  DMap f("var:x,y,z;fun:x,y;",3);
  cnTest(f,x,getIdExpected());
}

BOOST_AUTO_TEST_CASE(xcnSwapTest)
{
  DMap f("var:x,y,z;fun:y,x;",3);
  cnTest(f,x,getSwapExpected());
}

BOOST_AUTO_TEST_SUITE_END()
