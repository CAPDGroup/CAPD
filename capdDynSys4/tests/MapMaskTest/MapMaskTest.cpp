#define BOOST_TEST_MODULE MapMaskTest
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include "capd/capdlib.h"
using namespace capd;
typedef capd::vectalg::Matrix<bool,0,0> BMatrix;
typedef capd::diffAlgebra::Jet<BMatrix,0> BJet;

extern const char formula[]="time:t;par:p;var:x,y,z;fun:\
    1+2*x+3*y-4*z\
    +5*x^2+2*y^2+3*z^2-x*y-2*x*z+5*y*z\
    +x^2*(2*x+y-2*z)+y^2*(3*y-z-x)+z^2*(z+2*x+3*y)+x*y*z,\
    2+x-3*y-2*z\
    +2*x^2+3*y^2+4*z^2-x*y-7*x*z+2*y*z\
    +x^2*(x+2*y-3*z)+y^2*(y-4*z-2*x)+z^2*(2*z+3*x-2*y)-x*y*z;";
extern const DVector x(1.,2.,-3.);

void check(const BJet& mask, const DJet& expected, const DJet& out){
  BJet::const_iterator m = mask.begin();
  DJet::const_iterator o = out.begin();
  for(DJet::const_iterator e=expected.begin();e!=expected.end();++m,++e,++o){
    if(*m){
      BOOST_CHECK_EQUAL(*o,*e);
    } else {
      BOOST_CHECK_EQUAL(*o,0.);
    }
  }
}

void doTest(DMap& f, const BJet& mask, const DJet& expected){
  DJet df(2,3,3);
  DJet i(3,3,3);
  i() = x;
  i(0,0) = i(1,1) = i(2,2) = 1;

  f(x,df);
  check(mask,expected,df);
  f(10,x,df);
  check(mask,expected,df);
  df = f(i);
  check(mask,expected,df);
  df = f(10,i);
  check(mask,expected,df);
}

DJet getPolyExpected(){
  const double out[] = {116,46,64,-39,19,-6,-16,22,-16,2,2,1,-2,-1,1,2,3,-1,3,1,65,78,41,10,18,-2,-33,19,-3,-15,1,2,-3,-2,-1,3,1,-4,-2,2};
  DJet expected(2,3,3);
  std::copy(out,out+40,expected.begin());
  return expected;
}

DJet getIdExpected(){
  DJet expected(2,3,3);
  expected(0) = 1.;
  expected(1) = 2.;
  expected(0,0) = 1;
  expected(1,1) = 1;
  return expected;
}

DJet getSwapExpected(){
  DJet expected(2,3,3);
  expected(0) = 2.;
  expected(1) = 1.;
  expected(0,1) = 1;
  expected(1,0) = 1;
  return expected;
}
