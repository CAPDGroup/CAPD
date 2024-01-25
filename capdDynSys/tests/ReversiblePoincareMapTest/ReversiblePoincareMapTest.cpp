//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file ReversiblePoincareMapTest.cpp
///
/// @author Daniel Wilczak
//
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) CAPD group
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#define BOOST_TEST_MODULE ReversiblePoincareMapTest
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include "capd/capdlib.h"
#include "capd/dynset/C2DoubletonSet.hpp"
#include "capd/dynsys/C2OdeSolver.hpp"
using namespace capd;

/**
 * The Michelson system is reversible with respect to involution (x,y,z,t) -> (-x,y,-z,-t)
 * We check rigorous C^0-C^n solvers for preserving this property.
 */

BOOST_AUTO_TEST_SUITE(ReversiblePoincareMapTestSuite)

template<class S, class PM>
void doC0test(PM& pm, IVector x){
  S s(x);

  IVector z = pm(s);
  // apply reversing symmetry
  z[2] = -z[2];
  z[0] = 0;

  s = S(z);
  z = pm(s);

  // apply reversing symmetry
  z[2] = -z[2];
  z[0] = 0;

  // we should obtain identity map
  BOOST_REQUIRE_MESSAGE(subset(x,z),"x=" << x << ", z=" << z);
}

BOOST_AUTO_TEST_CASE(C0Test)
{
  IMap f="par:c;var:x,y,z;fun:y,z,1-y-0.5*x*x;";
  IOdeSolver solver(f,10);
  INonlinearSection section("var:x,y,z;fun:x;");
  IPoincareMap pm(solver,section);
  IVector x(0.,1.563,0.);

  doC0test<C0Rect2Set>(pm,x);
  doC0test<C0Pped2Set>(pm,x);
  doC0test<C0TripletonSet>(pm,x);
  doC0test<C0HORect2Set>(pm,x);
  doC0test<C0HOTripletonSet>(pm,x);
}


template<class S, class PM>
void doC1test(PM& pm, IVector x){
  S s(x);

  // matrix of reversing symmetry
  IMatrix R(3,3);
  R[0][0] = 0;
  R[2][2] = -1;
  R[1][1] = 1;

  IMatrix R2 = R*R; // Id on y,z, 0 on x

  IMatrix D(3,3);
  IVector z = pm(s,D);
  D = pm.computeDP(z,D);

  // apply reversing symmetry
  D = R*D*R2;
  z = R*z;
  s = S(typename S::C0BaseSet(z),typename S::C1BaseSet(D));

  // compute second iterate
  z = pm(s,D);
  D = pm.computeDP(z,D);

  // apply reversing symmetry
  D = R*D*R2;
  z = R*z;

  // we should obtain identity map
  BOOST_REQUIRE_MESSAGE(subset(x,z),"x=" << x << ", z=" << z);
  BOOST_REQUIRE_MESSAGE(subset(R2,D),"D=" << D);
}


BOOST_AUTO_TEST_CASE(C1Test)
{
  IMap f="par:c;var:x,y,z;fun:y,z,1-y-0.5*x*x;";
  IOdeSolver solver(f,10);
  INonlinearSection section("var:x,y,z;fun:x;");
  IPoincareMap pm(solver,section);
  IVector x(0.,1.563,0.);

  doC1test<C1Rect2Set>(pm,x);
  doC1test<C1Pped2Set>(pm,x);
  doC1test<C1HORect2Set>(pm,x);
  doC1test<C1HOPped2Set>(pm,x);
}

template<class S, class PM>
void doC2test(PM& pm, IVector x){
  S s(x);

  // matrix of reversing symmetry
  IMatrix R(3,3);
  R[0][0] = 0;
  R[2][2] = -1;
  R[1][1] = 1;

  IMatrix R2 = R*R; // Id on y,z, 0 on x

  IMatrix monodromyMatrix(3,3), D(3,3);
  IHessian hessianFlowH(3,3), H(3,3);
  IVector z = pm(s,monodromyMatrix,hessianFlowH);
  pm.computeDP(z,monodromyMatrix,hessianFlowH,D,H);

  // apply reversing symmetry
  z = R*z;
  D = R*D*R2;
  H = R*H*R2;
  s = S(typename S::C0BaseSet(z),typename S::C1BaseSet(D),H);

  // compute second iterate
  z = pm(s,monodromyMatrix,hessianFlowH);
  pm.computeDP(z,monodromyMatrix,hessianFlowH,D,H);

  // apply reversing symmetry
  D.column(0).clear();
  z = R*z;
  D = R*D*R2;
  H = R*H*R2;

  // we should obtain identity map
  BOOST_REQUIRE_MESSAGE(subset(x,z),"x=" << x << ", z=" << z);
  BOOST_REQUIRE_MESSAGE(subset(R2,D),"D=" << D);
  BOOST_REQUIRE(subset(IHessian(3,3),H));
}

BOOST_AUTO_TEST_CASE(C2Test)
{
  IMap f="par:c;var:x,y,z;fun:y,z,1-y-0.5*x*x;";
  IC2OdeSolver solver(f,10);
  INonlinearSection section("var:x,y,z;fun:x;");
  IC2PoincareMap pm(solver,section);
  IVector x(0.,1.563,0.);

  doC0test<C0Rect2Set>(pm,x);
  doC0test<C0Pped2Set>(pm,x);
  doC0test<C0TripletonSet>(pm,x);

  doC1test<C1Rect2Set>(pm,x);
  doC1test<C1Pped2Set>(pm,x);

  doC2test<C2Rect2Set>(pm,x);
}



template<class S, class PM>
void doCntest(PM& pm, IVector x){
  const int degree = pm.getSolver().getVectorField().degree();
  S s = S(x,degree);

  // matrix of reversing symmetry
  IMatrix R(3,3);
  R[0][0] = 0;
  R[2][2] = -1;
  R[1][1] = 1;

  IMatrix R2 = R*R; // Id on y,z, 0 on x

  // integrate over time 1
  IJet D(3,degree);
  IVector z = pm(s,D);
  D = pm.computeDP(D);

  // apply reversing symmetry
  D = R*D;
  // clear x-dependency
  for(int i=1;i<3;++i){
    D(i,0) = interval(0.);
    for(int j=0;j<3;++j){
      D(i,0,j) = interval(0.0);
      for(int k=j;k<3;++k)
        D(i,0,j,k) = interval(0.0);
    }
  }
  s = S(D);

  // integrate over time 1
  z = pm(s,D);
  D = pm.computeDP(D);

  // apply reversing symmetry
  D = R*D;
  // clear x-dependency
  for(int i=1;i<3;++i){
    D(i,0) = interval(0.);
    for(int j=0;j<3;++j){
      D(i,0,j) = interval(0.0);
      for(int k=j;k<3;++k)
        D(i,0,j,k) = interval(0.0);
    }
  }

  IJet J(x.dimension(),degree);
  J() = x;
  J.setMatrix(IMatrix::Identity(x.dimension()));
  J(0,0) = 0;
  // we should obtain identity map
  BOOST_REQUIRE_MESSAGE(subset(J,D),"J=" << J.toString() << ", D=" << D.toString());
}

// BOOST_AUTO_TEST_CASE(CnTest)
// {
//   IMap f("par:c;var:x,y,z;fun:y,z,1-y-0.5*x*x;",3);
//   ICnOdeSolver solver(f,10);
//   INonlinearSection section("var:x,y,z;fun:x;");
//   ICnPoincareMap pm(solver,section);
//   IVector x(0.,1.563,0.);

//   doC0test<C0Rect2Set>(pm,x);
//   doC0test<C0Pped2Set>(pm,x);
//   doC0test<C0TripletonSet>(pm,x);

//   doC1test<C1Rect2Set>(pm,x);
//   doC1test<C1Pped2Set>(pm,x);

//   doC2test<C2Rect2Set>(pm,x);

//   doCntest<CnRect2Set>(pm,x);
//   doCntest<CnMultiMatrixRect2Set>(pm,x);
// }

BOOST_AUTO_TEST_SUITE_END()
