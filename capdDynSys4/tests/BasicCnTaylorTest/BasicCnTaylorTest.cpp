/// @addtogroup test
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file BasicCnTaylorTest.cpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2013 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#define BOOST_TEST_MODULE BasicCnTaylorTest
#define BOOST_TEST_DYN_LINK

#include <iostream>
#include <iterator>
#include <boost/test/unit_test.hpp>
#include "capd/capdlib.h"
using namespace capd;
using namespace std;

// this is exported from Mathematica - exact solution
constexpr double value[] = {-2.5316314956714900, 1.8949517065413926, 3.4903429574618414};
constexpr double derivative[] = {
      -0.15841556424007908, -2.3701983992616307, 0,
      0.31607674543771402, -1.5834012593583489, 0,
      0.87258573936546036, 2.6177572180963811, 3.4903429574618414
  };
constexpr double hessian[] = {
  -0.39500163391683646, -0.31588814967710271, 0,  0.23866062304348090, 0, 0, -0.29684058718953776,
  -0.35605782801895652, 0, -1.5611475689403926, 0, 0, 0.54536608710341272, 0.65443930452409527,
  0.87258573936546036, 1.4179518264688731, 2.6177572180963811, 0};

constexpr double c3[] = {
  -0.044428646935448214, 0.21287091444760703, 0, 0.18297940039198052,
  0, 0, 0.59806587347974383, 0, 0, 0, -0.042864924531560857,
  -0.42472812668192761, 0, -0.16782147313347984, 0, 0,
  -0.20695387479197153, 0, 0, 0, 0.11816265220573942,
  0.40902456532755954, 0.54536608710341272, 0.35448795661721827,
  0.65443930452409527, 0, 0.57263439145858330, 1.4179518264688731, 0, 0};

constexpr double c4[] = {0.019385967054582082, 0.062625614738937477, 0, 0.24455695887102799,
0, 0, 0.088902478271876140, 0, 0, 0, 0.16416017749989503, 0, 0, 0, 0,
-0.035083622888334415, -0.040703994913533721, 0,
0.026791566449043028, 0, 0, 0.038378056699432450, 0, 0, 0,
0.13634464261114032, 0, 0, 0, 0, 0.041470546206822008,
0.088621989154304567, 0.11816265220573942, 0.22155497288576143,
0.40902456532755954, 0, 0.14315859786464583, 0.35448795661721827, 0,
0, 0.19599093755278890, 0.57263439145858330, 0, 0, 0};


constexpr double c5[] = {0.0041563948134421252, 0.034803930351018676, 0,
0.023390554801000430, 0, 0, 0.034964278720380194, 0, 0, 0,
0.00097020584133128179, 0, 0, 0, 0, -0.0089063093630541244, 0, 0, 0,
0, 0, -0.0029503766309958565, 0.015834312503170186, 0,
0.029359870162769891, 0, 0, 0.090516545165687523, 0, 0, 0,
0.035051577929882012, 0, 0, 0, 0, 0.061538593597794712, 0, 0, 0, 0,
0, 0.0079816599206280724, 0.031102909655116508, 0.041470546206822008,
0.048003577458581637, 0.088621989154304567, 0, 0.089474123665403651,
0.22155497288576143, 0, 0, 0.048997734388197225, 0.14315859786464583,
0, 0, 0, 0.058030360205847516, 0.19599093755278890, 0, 0, 0, 0};


constexpr double c6[] =  {0.0019094844343481466, 0.0025171660982330666, 0,
-0.0020379782822336464, 0, 0, -0.0064975747688835331, 0, 0, 0,
-0.020675969261991292, 0, 0, 0, 0, -0.0089294578491597582, 0, 0, 0,
0, 0, -0.014559460217323789, 0, 0, 0, 0, 0, 0,
0.00093020114033524036, 0.0050657028363276322, 0,
0.018658924707740599, 0, 0, 0.013028421936501985, 0, 0, 0,
0.025023343729418203, 0, 0, 0, 0, 0.0048900781120480696, 0, 0, 0, 0,
0, 0.0068025838171620187, 0, 0, 0, 0, 0, 0, 0.0020605085886437535,
0.0059862449404710543, 0.0079816599206280724, 0.016847409396521441,
0.031102909655116508, 0, 0.019386060127504125, 0.048003577458581637,
0, 0, 0.030623583992623271, 0.089474123665403651, 0, 0, 0,
0.014507590051461879, 0.048997734388197225, 0, 0, 0, 0,
0.015420084090430480, 0.058030360205847509, 0, 0, 0, 0, 0};


template<class iterator1, class iterator2>
void checkSmall(iterator1 b, iterator1 e, iterator2 i, double errorTolerance)
{
  for(;b!=e;++b,++i) 
  {
    BOOST_CHECK_SMALL(*b-*i,errorTolerance);
  }
}

// compare numerical solution to exact solution
void checkSmall(const DVector& x, const DMatrix& der, const DHessian& h, double errorTolerance)
{
  checkSmall(x.begin(),x.end(),value,errorTolerance);
  checkSmall(der.begin(),der.end(),derivative,errorTolerance);

  // we have different indexing in Hessian and CnCoeff!
  const double* p = hessian;
  typedef DVector::size_type size_type;
  for(size_type i=0;i<h.dimension();++i)
    for(size_type j=0;j<h.dimension();++j)
      for(size_type c=j;c<h.dimension();++c,++p)
      {
        BOOST_CHECK_SMALL(*p-h(i,j,c),errorTolerance);
      }
}

// compare numerical solution to exact solution
void checkSmall(const DJet& c, double errorTolerance)
{
  DVector x = c();
  checkSmall(x.begin(),x.end(),value,errorTolerance);
  DMatrix der = (DMatrix)c;
  checkSmall(der.begin(),der.end(),derivative,errorTolerance);

  // no breaks in this switch!
  switch(c.degree())
  {
    case 6:
      checkSmall(c.begin(0,6),c.end(0,6),c6,errorTolerance);
      checkSmall(c.begin(1,6),c.end(1,6),c6+28,errorTolerance);
      checkSmall(c.begin(2,6),c.end(2,6),c6+56,errorTolerance);
    case 5:
      checkSmall(c.begin(0,5),c.end(0,5),c5,errorTolerance);
      checkSmall(c.begin(1,5),c.end(1,5),c5+21,errorTolerance);
      checkSmall(c.begin(2,5),c.end(2,5),c5+42,errorTolerance);
    case 4:
      checkSmall(c.begin(0,4),c.end(0,4),c4,errorTolerance);
      checkSmall(c.begin(1,4),c.end(1,4),c4+15,errorTolerance);
      checkSmall(c.begin(2,4),c.end(2,4),c4+30,errorTolerance);
    case 3:
      checkSmall(c.begin(0,3),c.end(0,3),c3,errorTolerance);
      checkSmall(c.begin(1,3),c.end(1,3),c3+10,errorTolerance);
      checkSmall(c.begin(2,3),c.end(2,3),c3+20,errorTolerance);
    case 2:
      checkSmall(c.begin(0,2),c.end(0,2),hessian,errorTolerance);
      checkSmall(c.begin(1,2),c.end(1,2),hessian+6,errorTolerance);
      checkSmall(c.begin(2,2),c.end(2,2),hessian+12,errorTolerance);
  }
}

void doTest(std::string vfFormula)
{
  int order=15;         // of Taylor method
  const unsigned degree = 6; // max derivative in test
  const unsigned dim = 3;
  const double errorTolerance = 1e-13;

  double T=0.125;
  double step(T/32);    // fixed time step
  int numberOfSteps = T/step;

  DMap vf(vfFormula,degree);

  DC2Taylor c2solver(vf,order);
  DCnTaylor cnsolver(vf,order);
  c2solver.setStep(step);
  cnsolver.setStep(step);

  DVector u1(1,3,1);
  DVector u2 = u1;
  DMatrix D1 = DMatrix::Identity(dim);
  DMatrix D2 = D1;
  DHessian h1(u1.dimension());
  DHessian h2(u1.dimension());

  DCnTaylor::JetType c(u1.dimension(),degree);
  c() = u1;
  c.setMatrix(D1);

  double t=1;

  for(int i=0;i<numberOfSteps;++i)
    u1 = c2solver(t,u1,D1,h1,D1,h1);
  checkSmall(u1,D1,h1,errorTolerance);

  t=1;
  for(int i=0;i<numberOfSteps;++i)
    u2 = cnsolver(t,u2,D2,h2,D2,h2);
  checkSmall(u2,D2,h2,errorTolerance);

  t=1;
  for(int i=0;i<numberOfSteps;++i)
    cnsolver(t,c);
  checkSmall(c,errorTolerance);
}

BOOST_AUTO_TEST_SUITE(TestSuite)


BOOST_AUTO_TEST_CASE(x1VarMulTime)
{
  // 1. very basic test for nodes: NODE_MUL, UNARY_MINUS, NODE_ADD
  doTest("var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z;");
  // 2. time added
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z;");
  // 3. TIME_PLUS_VAR, VAR_MINUS_TIME
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),t+x*(x*x+y*y)-t,(x*x+y*y)*z;");
  // 4. TIME_MINUS_VAR, VAR_MINUS_TIME, TIME_PLUS_VAR
  doTest("time:t;var:x,y,z;fun:t-y*(x*x+y*y)-t,t+x*(x*x+y*y)-t,(x*x+y*y)*z;");
  // 5. TIME_PLUS_CONST, TIME_MINUS_CONST, MUL_FUNTIME_BY_FUNTIME, MUL_TIME_BY_FUNTIME, FUNTIME_MINUS_CONST
  doTest("time:t;var:x,y,z;fun:(t+1)*(t-1)-y*(x*x+y*y)-(t*t-1),x*(x*x+y*y),(x*x+y*y)*z;");
  // 6. TIME_PLUS_CONST, TIME_MINUS_CONST, MUL_FUNTIME_BY_FUNTIME, MUL_TIME_BY_FUNTIME, FUNTIME_MINUS_CONST, FUNTIME_MINUS_VAR, VAR_MINUS_FUNTIME
  doTest("time:t;var:x,y,z;fun:(t+1)*(t-1)-y*(x*x+y*y)-(t*t-1),x*(x*x+y*y),(x*x+y*y)*z;");
  // 7. ADDITIONALY FUNTIME_MINUS_FUNTIME
  doTest("time:t;var:x,y,z;fun:((t+1)*(t-1)-(t*t-1))-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z;");
  // 8. NODE_MUL_CONST_BY_CONST
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),(2*0.5)*x*(x*x+y*y),(x*x+y*y)*z;");
  // 9. NODE_MUL_CONST_BY_VAR
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),2*x*(x*x+y*y)*0.5,(x*x+y*y)*z;");
}

BOOST_AUTO_TEST_CASE(x2Div)
{
  // 10. DIV
  doTest("time:t;var:x,y,z;fun:-y*(1/y)*y*(x*x+y*y),x*(x*x+y*y)*((x*x+y*y+1)/(x*x+y*y+10))*((x*x+y*y+10)/(x*x+y*y+1)),(x*x+y*y)*z;");
  // 11. DIV_VAR_BY_CONST
  doTest("time:t;var:x,y,z;fun:-y*(1/y)*y*(x*x+y*y),4*x*(x*x+y*y)/2/2,(x*x+y*y)*z;");
  // 12. DIV_VAR_BY_FUNTIME
  doTest("time:t;var:x,y,z;fun:-sin(2*t)*y*(x*x+y*y)/(2*sin(t)*cos(t)),x*(x*x+y*y),(t*t-0.25)*(x*x+y*y)*z/(t-0.5)/(t+0.5);");
  // 13. DIV_VAR_BY_TIME
  doTest("time:t;var:x,y,z;fun:-y*(1/y)*y*(x*x+y*y)+(x+y)/t-x/t-y/t,t*x*(x*x+y*y)/t,(x*x+y*y)*z;");
  // 14. DIV_TIME_BY_CONST
  doTest("time:t;var:x,y,z;fun:-(t/3)*y*(x*x+y*y)*3/t,x*(x*x+y*y),(x*x+y*y)*z;");
  // 15. DIV_FUNTIME_BY_CONST
  doTest("time:t;var:x,y,z;fun:-(sin(2*t)/2)*y*(x*x+y*y)/sin(t)/cos(t),x*(x*x+y*y),(x*x+y*y)*z;");
  // 16. DIV_FUNTIME_BY_TIME
  doTest("time:t;var:x,y,z;fun:-(sin(t)/t)*y*(x*x+y*y)*t/sin(t),x*(x*x+y*y),(x*x+y*y)*z;");
  // 17. DIV_FUNTIME_BY_FUNTIME
  doTest("time:t;var:x,y,z;fun:-exp(2*sin(t))/exp(sin(t))/exp(sin(t))*y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z;");
  // 18. DIV_CONST_BY_CONST
  doTest("time:t;var:x,y,z;fun:-(7/5)*y*(x*x+y*y)*(5/7),x*(x*x+y*y),(x*x+y*y)*z;");
}

BOOST_AUTO_TEST_CASE(x2ExpLog)
{
  // 21. EXP
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+exp(x+y)-exp(x)*exp(y);");
  // 22. EXP_CONST
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+exp(6)-exp(3)*exp(3);");
  // 23. EXP_TIME
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+sqrt(exp(t)*exp(t))-exp(t);");
  // 24. EXP_FUNTIME
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+exp(t*t-t)-exp(t^2)/exp(t);");
  // 25. LOG
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+log((x*x+y*y)*(x*x+z*z))-log(x*x+y*y)-log(x*x+z*z);");
  // 26. LOG_CONST
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+log(10)-log(2)-log(5);");
  // 27. LOG_TIME
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),sqr(log(t+1))*(x*x+y*y)*z/log(t+1)/log(t+1);");
  // 28. LOG_FUNTIME
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+log((t+1)^2)-2*log(t+1);");
}

BOOST_AUTO_TEST_CASE(x3SqrSqrt)
{
  // 31. SQR 
  doTest("var:x,y,z;fun:-y*(x^2+y^2),x*(x^2+y^2),(x^2+y^2)*z;");
  // 32. SQR_CONST
  doTest("var:x,y,z;fun:3^2-y*(x^2+y^2)-9,x*(x^2+y^2),(x^2+y^2)*z;");
  // 33. SQR_TIME
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(t^2-t*t)+(x*x+y*y)*z;");
  // 34. SQR_FUNTIME
  doTest("time:t;var:x,y,z;fun:(t+1)^2-y*(x^2+y^2)-t^2-2*t-1,x*(x^2+y^2),(x^2+y^2)*z;");
  // 35. SQRT
  doTest("var:x,y,z;fun:sqrt((x*x+2*x+y*y*y*y+1)*(x*x+2*x+y*y*y*y+1))-(x+1)*(x+1)-(y*y)*(y*y)-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z;");
  // 36. SQRT_CONST
  doTest("var:x,y,z;fun:3-y*(x^2+y^2)-sqrt(9),x*(x^2+y^2),(x^2+y^2)*z;");
  // 37. SQRT_TIME
  doTest("time:t;var:x,y,z;fun:(sqrt(t)+0.5)/(sqrt(t)-0.5)-(t+sqrt(t)+0.25)/(t-0.25)-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z;");
  // 38. SQRT_FUNTIME
  doTest("time:t;var:x,y,z;fun:t+1-sqrt(t^2+2*t+1)-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z;");
}

BOOST_AUTO_TEST_CASE(x4Pow)
{
  // 41. POW
  doTest("var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+(sqr(1/sqr(5+x*x+y*y)))^0.25-1/(5+x*x+y*y);");
  // 42. POW_CONST
  doTest("var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+(256^0.125)-2;");
  // 43. POW_TIME
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+sqr(sqr(t^.25))-t;");
  // 44. POW_FUNTIME
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),t+1-(t*t*t*t+4*t*t*t+6*t*t+4*t+1)^0.25+x*(x*x+y*y),(x*x+y*y)*z;");
}

BOOST_AUTO_TEST_CASE(x5SinCos)
{
  // 51. SIN/COS
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+sin(x+y)-sin(x)*cos(y)-sin(y)*cos(x);");
  // 52. SIN/COS_CONST
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+sin(3+4)-sin(3)*cos(4)-sin(4)*cos(3);");
  // 53. SIN/COS_TIME
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+sin(t)^2+cos(t)^2-1;");
  // 54. SIN/COS_FUNTIME
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+sin(2*t)-2*sin(t)*cos(t);");
}

BOOST_AUTO_TEST_CASE(x5Atan)
{
  // 55. ATAN
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+sin(atan(x*y+1))/cos(atan(x*y+1))-x*y-1;");
  // 56. ATAN_CONST
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+2*sin(atan(1))-sqrt(2);");
  // 57. ATAN_TIME
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+sin(atan(t))/cos(atan(t))-t;");
  // 58. ATAN_FUNTIME
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+atan(t)-2*atan(t/(1+sqrt(1+t*t)));");
}

BOOST_AUTO_TEST_CASE(x6ASinACos)
{
  // 61. ASIN
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+asin(1/(x^2+2*y^2+2))-2*atan((1/(x^2+2*y^2+2))/(1+sqrt(1-1/(x^2+2*y^2+2)^2)));");
  // 62. ASIN_CONST
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+asin(1./sqrt(2.))-atan(1);");
  // 63. ASIN_TIME
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+asin(t-1)-2*atan((t-1)/(1+sqrt(1-(t-1)^2)));");
  // 64. ASIN_FUNTIME
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+atan(t)-asin(t/sqrt(1+t^2));");
  // 65. ACOS
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+acos(1/(x^2+2*y^2+2))-2*atan(sqrt(1-1/(x^2+2*y^2+2)^2)/(1+1/(x^2+2*y^2+2)));");
  // 66. ACOS_CONST
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+acos(0.5)-atan(sqrt(3));");
  // 67. ACOS_TIME
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+acos(t-1)-2*atan(sqrt(1-(t-1)^2)/t);");
  // 68. ACOS_FUNTIME
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+atan(t)-acos(1/sqrt(1+t^2));");
}

BOOST_AUTO_TEST_CASE(x7CubeQuartic)
{
  // 71. CUBE
  doTest("var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+1/(sqr(5+x*x+y*y))^3-(5+x*x+y*y)^(-6);");
  // 72. CUBE_CONST
  doTest("var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+3^3+2^3-35;");
  // 73. CUBE_TIME
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+t^3-t*t*t;");
  // 74. CUBE_FUNTIME
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),(t+1)^3-(t^3+3*t*t+3*t+1)+x*(x*x+y*y),(x*x+y*y)*z;");

  // 75. QUARTIC
  doTest("var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+1/((1+x*x+y*y))^4-sqr(sqr(1/(1+x*x+y*y)));");
  // 76. QUARTIC_CONST
  doTest("var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+3^4+2^4-97;");
  // 77. QUARTIC_TIME
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+t^4-sqr(sqr(t));");
  // 78. QUARTIC_FUNTIME
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),(t+1)^4-(t^4+4*t*t*t+6*t*t+4*t+1)+x*(x*x+y*y),(x*x+y*y)*z;");
}

BOOST_AUTO_TEST_SUITE_END()

/// @}
