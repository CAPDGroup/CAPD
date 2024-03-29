#define BOOST_TEST_MODULE C1MaskTest
#define BOOST_TEST_DYN_LINK
#include "compare.h"

template<class TM>
  void doC0Test(TM& tm){
  DVector u(1,3,1);
  double t0 = 1;
  u = tm(t0+0.125,u,t0);
  check(u);
}

template<class TM>
void doC1Test(TM& tm, const BMatrix& mask){
  DVector u(1,3,1);
  DMatrix D = DMatrix::Identity(3);
  double t0=1;
  u = tm(t0+0.125,u,D,t0);
  check(u);
  check(mask,D);
}

template<class TM>
void doC1Test(std::string vf){
  DMap f(vf,6);
  typename TM::Solver solver(f,15);
  TM tm(solver);
  doC0Test(tm);

  BMatrix m(3,3);
  solver.setMask((Multiindex*)0,(Multiindex*)0);
  doC1Test(tm,m);

  Multiindex m1[] = {{1,0,0},{0,1,0},{0,0,1}};
  for(int j=0;j<3;++j){
    solver.addMultiindexToMask(m1[j]);
    m.column(j) = {true,true,true};
    doC1Test(tm,m);
  }
}

void doTest(std::string vf){
  doC1Test<DTimeMap>(vf);
  doC1Test<DC2TimeMap>(vf);
  doC1Test<DCnTimeMap>(vf);
}

BOOST_AUTO_TEST_SUITE(C1TestSuite)

BOOST_AUTO_TEST_CASE(xC1Test)
{
  // 0. very basic test for nodes: NODE_MUL, UNARY_MINUS, NODE_ADD
  doTest("var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z;");
  // 1. Add/subtract const
  doTest("var:x,y,z;fun:5-(y*(x*x+y*y)+5),x*(x*x+y*y),(x*x+y*y)*z;");
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

  // ------------------------- DIV ----------------------------------------------
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

   // ------------------------- EXP-LOG ----------------------------------------------
  // 21. EXP
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+exp(x+y)-exp(x)*exp(y);");
  // 22. EXP_CONST
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+exp(6)-exp(3)*exp(3);");
  // 23. EXP_TIME
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+sqrt(exp(t)*exp(t))-exp(t);");
  // 24. EXP_FUNTIME
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+exp(t*t-t)-exp(t^2)/exp(t);");
  // 25. LOG
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),log((x*x+y*y)*(x*x+z*z))+(x*x+y*y)*z-log(x*x+y*y)-log(x*x+z*z);");
  // 26. LOG_CONST
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),log(10)+(x*x+y*y)*z-log(2)-log(5);");
  // 27. LOG_TIME
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),sqr(log(t+1))*(x*x+y*y)*z/log(t+1)/log(t+1);");
  // 28. LOG_FUNTIME
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+log((t+1)^2)-2*log(t+1);");

  // ------------------------- SQR ----------------------------------------------
  // 31. SQR
  doTest("var:x,y,z;fun:-y*(x^2+y^2),x*(x^2+y^2),(x^2+y^2)*z;");
  // 32. SQR_CONST
  doTest("var:x,y,z;fun:sqr(3)-y*(x^2+y^2)-9,x*(x^2+y^2),(x^2+y^2)*z;");
  // 33. SQR_TIME
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(t^2-t*t)+(x*x+y*y)*z;");
  // 34. SQR_FUNTIME
  doTest("time:t;var:x,y,z;fun:(t+1)^2-y*(x^2+y^2)-t^2-2*t-1,x*(x^2+y^2),(x^2+y^2)*z;");

  // ------------------------- SQRT ----------------------------------------------
  // 35. SQRT
  doTest("var:x,y,z;fun:sqrt((x*x+2*x+y*y*y*y+1)*(x*x+2*x+y*y*y*y+1))-(x+1)*(x+1)-(y*y)*(y*y)-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z;");
  // 36. SQRT_CONST
  doTest("var:x,y,z;fun:3-y*(x^2+y^2)-sqrt(9),x*(x^2+y^2),(x^2+y^2)*z;");
  // 37. SQRT_TIME
  doTest("time:t;var:x,y,z;fun:(sqrt(t)+0.5)/(sqrt(t)-0.5)-(t+sqrt(t)+0.25)/(t-0.25)-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z;");
  // 38. SQRT_FUNTIME
  doTest("time:t;var:x,y,z;fun:t+1-sqrt(t^2+2*t+1)-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z;");

  // ------------------------- POW ----------------------------------------------
  // 41. POW
  doTest("var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+(sqr(1/sqr(5+x*x+y*y)))^0.25-1/(5+x*x+y*y);");
  // 42. POW_CONST
  doTest("var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+(256^0.125)-2;");
  // 43. POW_TIME
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+sqr(sqr(t^.25))-t;");
  // 44. POW_FUNTIME
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),t+1-(t*t*t*t+4*t*t*t+6*t*t+4*t+1)^0.25+x*(x*x+y*y),(x*x+y*y)*z;");

  // ------------------------- SIN-COS ----------------------------------------------
  // 51. SIN/COS
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+sin(x+y)-sin(x)*cos(y)-sin(y)*cos(x);");
  // 52. SIN/COS_CONST
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+sin(3+4)-sin(3)*cos(4)-sin(4)*cos(3);");
  // 53. SIN/COS_TIME
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+sin(t)^2+cos(t)^2-1;");
  // 54. SIN/COS_FUNTIME
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+sin(2*t)-2*sin(t)*cos(t);");

  // ------------------------- ATAN ----------------------------------------------
  // 61. ATAN
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+sin(atan(x*y+1))/cos(atan(x*y+1))-x*y-1;");
  // 62. ATAN_CONST
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+2*sin(atan(1))-sqrt(2);");
  // 63. ATAN_TIME
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+sin(atan(t))/cos(atan(t))-t;");
  // 64. ATAN_FUNTIME
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+atan(t)-2*atan(t/(1+sqrt(1+t*t)));");

  // ------------------------- ASIN ----------------------------------------------

  // 65. ASIN
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+asin(1/(x^2+2*y^2+2))-2*atan((1/(x^2+2*y^2+2))/(1+sqrt(1-1/(x^2+2*y^2+2)^2)));");
  // 66. ASIN_CONST
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+asin(1./sqrt(2.))-atan(1);");
  // 67. ASIN_TIME
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+asin(t-1)-2*atan((t-1)/(1+sqrt(1-(t-1)^2)));");
  // 68. ASIN_FUNTIME
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+atan(t)-asin(t/sqrt(1+t^2));");

  // ------------------------- ACOS ----------------------------------------------

  // 69. ACOS
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+acos(1/(x^2+2*y^2+2))-2*atan(sqrt(1-1/(x^2+2*y^2+2)^2)/(1+1/(x^2+2*y^2+2)));");
  // 70. ACOS_CONST
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+acos(0.5)-atan(sqrt(3));");
  // 71. ACOS_TIME
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+acos(t-1)-2*atan(sqrt(1-(t-1)^2)/t);");
  // 72. ACOS_FUNTIME
  doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+atan(t)-acos(1/sqrt(1+t^2));");
}

BOOST_AUTO_TEST_SUITE_END()
