#ifndef __CAPD_REVERSIVBLE_TEST_H__
#define __CAPD_REVERSIVBLE_TEST_H__

#include "capd/capdlib.h"
using namespace capd;

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>


template<class S, class TM>
void doC0test(TM& tm, IVector x){
  S s(x);

  // integrate over time 1
  IVector z = tm(1.,s);

  // apply reversing symmetry
  z[2] = -z[2];
  z[0] = -z[0];

  s = S(z);
  // integrate over time 1
  z = tm(1.,s);

  // apply reversing symmetry
  z[2] = -z[2];
  z[0] = -z[0];

  // we should obtain identity map
  BOOST_REQUIRE_MESSAGE(subset(x,z),"x=" << x << ", z=" << z);
}

template<class S, class TM>
void doC1test(TM& tm, IVector x){
  S s(x);

  // integrate over time 1
  IVector z = tm(1.,s);
  IMatrix D = IMatrix(s);

  // matrix of reversing symmetry
  IMatrix R(3,3);
  R[0][0] = -1;
  R[2][2] = -1;
  R[1][1] = 1;

  // apply reversing symmetry
  z = R*z;
  D = R*D;

  s = S(typename S::C0BaseSet(z),typename S::C1BaseSet(D));
  // integrate over time 1
  z = tm(1.,s);
  D = IMatrix(s);

  // apply reversing symmetry
  z = R*z;
  D = R*D;

  // we should obtain identity map
  BOOST_REQUIRE_MESSAGE(subset(x,z),"x=" << x << ", z=" << z);
  BOOST_REQUIRE_MESSAGE(subset(IMatrix::Identity(3),D),"D=" << D);
}

template<class S, class TM>
void doC2test(TM& tm, IVector x){
  S s = S(x);

  // matrix of reversing symmetry
  IMatrix R(3,3);
  R[0][0] = -1;
  R[2][2] = -1;
  R[1][1] = 1;

  // integrate over time 1
  IVector z = tm(1.,s);

  IMatrix D = IMatrix(s);
  IHessian H = IHessian(s);

  // apply reversing symmetry
  z = R*z;
  D = R*D;
  H = R*H;

  s = S(typename S::C0BaseSet(z),typename S::C1BaseSet(D),H);
  // integrate over time 1
  z = tm(1.,s);
  D = IMatrix(s);
  H = IHessian(s);

  // apply reversing symmetry
  z = R*z;
  D = R*D;
  H = R*H;

  // we should obtain identity map
  BOOST_REQUIRE_MESSAGE(subset(x,z),"x=" << x << ", z=" << z);
  BOOST_REQUIRE_MESSAGE(subset(IMatrix::Identity(3),D),"D=" << D);
  BOOST_REQUIRE(subset(IHessian(3,3),H));
}

template<class S, class TM>
void doCntest(TM& tm, IVector x){
  const int degree = tm.getSolver().getVectorField().degree();
  S s = S(x,degree);

  // matrix of reversing symmetry
  IMatrix R(3,3);
  R[0][0] = R[2][2] = -1;
  R[1][1] = 1;

  // integrate over time 1
  IVector z = tm(1.,s);
  IJet D = s.currentSet();

  // apply reversing symmetry
  D = R*D;
  s = S(D);

  // integrate over time 1
  z = tm(1.,s);
  D = s.currentSet();

  // apply reversing symmetry
  D = R*D;

  IJet J(x.dimension(),degree);
  J() = x;
  J.setMatrix(IMatrix::Identity(x.dimension()));
  // we should obtain identity map
  BOOST_REQUIRE_MESSAGE(subset(J,D),"J=" << J.toString() << ", D=" << D.toString());
}

#endif
