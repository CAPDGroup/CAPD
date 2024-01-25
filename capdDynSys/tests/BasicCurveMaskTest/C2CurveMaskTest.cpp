#define BOOST_TEST_MODULE C2MaskTest
#define BOOST_TEST_DYN_LINK
#include "compare.h"

template<class TM>
void doC2Test(TM& tm, const BMatrix& dMask, const BHessian& hMask){
  DVector u(1,3,1);
  DMatrix D = DMatrix::Identity(3);
  DHessian H(3,3);
  double t0=1;
  typename TM::SolutionCurve solution(t0);
  double finalTime = 0.125;
  u = tm(t0+finalTime,u,D,H,solution);
  check(dMask,D);
  check(hMask,H);
  double t[] = {0.0,finalTime/3,finalTime/2,finalTime*0.75,finalTime};
  for(int i=0;i<5;++i){
    check(dMask,solution.derivative(t0+t[i]));
    check(hMask,solution.hessian(t0+t[i]));
  }
}

template<class TM>
void doC2Test(std::string vf){
  DMap f(vf,6);
  typename TM::Solver solver(f,20);
  TM tm(solver);

  BMatrix dm(3,3);
  BHessian hm(3,3);
  solver.setMask((Multiindex*)0,(Multiindex*)0);
  doC2Test(tm,dm,hm);

  Multiindex m1[] = {{1,0,0},{0,1,0}};
  for(int j=0;j<2;++j){
    solver.addMultiindexToMask(m1[j]);
    dm.column(j) = {true,true,true};
    doC2Test(tm,dm,hm);
  }

  solver.addMultiindexToMask({0,0,2});
  for(int i=0;i<3;++i){
    dm(i+1,3) = true;
    hm(i,2,2) = true;
  }
  doC2Test(tm,dm,hm);

  solver.addMultiindexToMask({1,1,0});
  for(int i=0;i<3;++i){
    hm(i,0,1) = true;
  }
  doC2Test(tm,dm,hm);

  solver.addMultiindexToMask({2,0,0});
  for(int i=0;i<3;++i){
    hm(i,0,0) = true;
  }
  doC2Test(tm,dm,hm);

  solver.addMultiindexToMask({0,1,1});
  for(int i=0;i<3;++i){
    hm(i,1,2) = true;
  }
  doC2Test(tm,dm,hm);
}


void doTest(std::string vf){
  doC2Test<DC2TimeMap>(vf);
  doC2Test<DCnTimeMap>(vf);
}

BOOST_AUTO_TEST_SUITE(C2TestSuite)

BOOST_AUTO_TEST_CASE(xC2Test)
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
    doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+log((x*x+y*y)*(x*x+z*z))-log(x*x+y*y)-log(x*x+z*z);");
    // 26. LOG_CONST
    doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+log(10)-log(2)-log(5);");
    // 27. LOG_TIME
    doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),sqr(log(t))*(x*x+y*y)*z/log(t)/log(t);");
    // 28. LOG_FUNTIME
    doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+log((t+1)^2)-2*log(t+1);");

    // ------------------------- SQR ----------------------------------------------
    // 31. SQR
    doTest("var:x,y,z;fun:-y*(x^2+y^2),x*(x^2+y^2),(x^2+y^2)*z;");
    // 32. SQR_CONST
    doTest("var:x,y,z;fun:3^2-y*(x^2+y^2)-9,x*(x^2+y^2),(x^2+y^2)*z;");
    // 33. SQR_TIME
    doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(t^2-t*t)+(x*x+y*y)*z;");
    // 34. SQR_FUNTIME
    doTest("time:t;var:x,y,z;fun:(t+1)^2-y*(x^2+y^2)-t^2-2*t-1,x*(x^2+y^2),(x^2+y^2)*z;");

    // ------------------------- SQRT ----------------------------------------------
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
    // 55. ATAN
    doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+sin(atan(x*y+1))/cos(atan(x*y+1))-x*y-1;");
    // 56. ATAN_CONST
    doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+2*sin(atan(1))-sqrt(2);");
    // 57. ATAN_TIME
    doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+sin(atan(t))/cos(atan(t))-t;");
    // 58. ATAN_FUNTIME
    doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+atan(t)-2*atan(t/(1+sqrt(1+t*t)));");

    // ------------------------- ASIN ----------------------------------------------
    // 61. ASIN
    doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+asin(1/(x^2+2*y^2+2))-2*atan((1/(x^2+2*y^2+2))/(1+sqrt(1-1/(x^2+2*y^2+2)^2)));");
    // 62. ASIN_CONST
    doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+asin(1./sqrt(2.))-atan(1);");
    // 63. ASIN_TIME
    doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+asin(t-1)-2*atan((t-1)/(1+sqrt(1-(t-1)^2)));");
    // 64. ASIN_FUNTIME
    doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+atan(t)-asin(t/sqrt(1+t^2));");

    // ------------------------- ACOS ----------------------------------------------

    // 65. ACOS
    doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+acos(1/(x^2+2*y^2+2))-2*atan(sqrt(1-1/(x^2+2*y^2+2)^2)/(1+1/(x^2+2*y^2+2)));");
    // 66. ACOS_CONST
    doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+acos(0.5)-atan(sqrt(3));");
    // 67. ACOS_TIME
    doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+acos(t-1)-2*atan(sqrt(1-(t-1)^2)/t);");
    // 68. ACOS_FUNTIME
    doTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+atan(t)-acos(1/sqrt(1+t^2));");

    // ------------------------- CUBE ----------------------------------------------
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
