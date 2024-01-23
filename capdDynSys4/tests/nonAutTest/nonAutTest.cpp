/// @addtogroup nonAutTest
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file nonAutTest.cpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2009 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.


/**
 * This test compares solutions to a nonautonomous system
 * x' = f(t,x)
 * to solutions to the autonomous system
 * y' = F(y)
 * with y = (t,x) and F(t,x) = (1,f(t,x))
 * 
 * Rigorous enclosures to both systems should have nonempty intersection.
 */ 


#define BOOST_TEST_MODULE NonAutTest
#define BOOST_TEST_DYN_LINK

#include <iostream>
#include <stdexcept>
#include "capd/capdlib.h"
#include <boost/test/unit_test.hpp>

using namespace capd;
using namespace std;

// -------------------------------------------------------------------------

void checkVectorIntersection(const IVector& u1, const IVector& u2, interval time)
{
  interval temp;
  BOOST_CHECK_EQUAL(intersection(u1[0],u2[0],temp),true);
  BOOST_CHECK_EQUAL(intersection(u1[1],u2[1],temp),true);
  BOOST_CHECK_EQUAL(intersection(u1[2],time,temp),true);
}

void checkMatrixIntersection(const IMatrix& m1, const IMatrix& m2)
{
  interval temp;
  BOOST_CHECK_EQUAL(intersection(m1(1,1),m2(1,1),temp),true);
  BOOST_CHECK_EQUAL(intersection(m1(1,2),m2(1,2),temp),true);
  BOOST_CHECK_EQUAL(intersection(m1(2,1),m2(2,1),temp),true);
  BOOST_CHECK_EQUAL(intersection(m1(2,2),m2(2,2),temp),true);
}

template<class C2CoeffType>
void checkHessianIntersection(const C2CoeffType& c1, const C2CoeffType& c2)
{
  interval temp;
  for(int i=0;i<2;++i)
    for(int j=0;j<2;++j)
      for(int k=j;k<2;++k)
        BOOST_CHECK_EQUAL(intersection(c1(i,j,k),c2(i,j,k),temp),true);
}

template<class CnCoeffType>
void checkThirdOrderIntersection(const CnCoeffType& c1, const CnCoeffType& c2)
{
  interval temp;
  for(int i=0;i<2;++i)
    for(int j=0;j<2;++j)
      for(int k=j;k<2;++k)
        for(int s=k;s<2;++s)
          BOOST_CHECK_EQUAL(intersection(c1(i,j,k,s),c2(i,j,k,s),temp),true);
}

void checkIntersection(C0Rect2Set&s, C0Rect2Set& nas)
{
  checkVectorIntersection(IVector(s),IVector(nas),nas.getCurrentTime());
}

void checkIntersection(C1Rect2Set&s, C1Rect2Set& nas)
{
  checkVectorIntersection(IVector(s),IVector(nas),nas.getCurrentTime());
  checkMatrixIntersection(IMatrix(s),IMatrix(nas));
}

void checkIntersection(C2Rect2Set&s, C2Rect2Set& nas)
{
  checkVectorIntersection(IVector(s),IVector(nas),nas.getCurrentTime());
  checkMatrixIntersection(IMatrix(s),IMatrix(nas));
  checkHessianIntersection(IHessian(s),IHessian(nas));
}

void checkIntersection(CnRect2Set&s, CnRect2Set& nas)
{
  checkVectorIntersection(IVector(s),IVector(nas),nas.getCurrentTime());
  checkMatrixIntersection(IMatrix(s),IMatrix(nas));
  checkHessianIntersection(s.currentSet(),nas.currentSet());
  checkThirdOrderIntersection(s.currentSet(),nas.currentSet());
}

// -------------------------------------------------------------------------

template<class Solver,class Set1, class Set2>
void Test(Solver& T, Solver& naT, Set1& s, Set2& nas, int order, int numberOfSteps, string description)
{
  naT.setOrder(order);
  T.setOrder(order);

  cout << "\n------------------------------------------------------------------------------------------------------------\n\n";
  cout << "Computing " << description << " trajectory using Solver and NonAutSolver for the forced pendulum.\nInitial point: " << IVector(s) << endl;
  for(int i=0;i<numberOfSteps;++i)
  {
    s.move(T);
    nas.move(naT);
    checkIntersection(s,nas);
  }
}

BOOST_AUTO_TEST_SUITE(TestSuite)

BOOST_AUTO_TEST_CASE(xNonAutTest)
{
  cout.precision(12);
  int order=6;
  interval step(0.125);

  string nonAutPendulumFormula = "time:t;par:beta;var:x,dx;fun:dx,cos(t)-beta*x-sin(x);";
  string pendulumFormula = "par:beta;var:x,dx,t;fun:dx,cos(t)-beta*x-sin(x),1;";
  IMap naVF(nonAutPendulumFormula,3);
  IMap VF(pendulumFormula,3);
  naVF.setParameter("beta",0.1);
  VF.setParameter("beta",0.1);

  ICnOdeSolver T(VF,order);
  ICnOdeSolver naT(naVF,order);
  T.setStep(step);
  naT.setStep(step);

  interval coeff[] = {interval(1),interval(2),interval(0)};
  IVector v(3,coeff);
  IVector nav(2,coeff);

// C0 Test
  C0Rect2Set s(v);
  C0Rect2Set nas(nav,0.);
  Test(T,naT,s,nas,10,100,"C0");

// C1Test
  C1Rect2Set c1s(v);
  C1Rect2Set c1nas(nav,0.);
  Test(T,naT,c1s,c1nas,10,100,"C1");

// C2Test
  C2Rect2Set c2s(v);
  C2Rect2Set c2nas(nav,0.);
  Test(T,naT,c2s,c2nas,10,100,"C2");

// C3Test
  CnRect2Set c3s(v,3);
  CnRect2Set c3nas(nav,3,0.);
  Test(T,naT,c3s,c3nas,10,100,"C3");
}

BOOST_AUTO_TEST_SUITE_END()
