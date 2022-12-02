#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE CnMaskTest
#include "compare.h"

// compare numerical solution to exact solution
void check(const BJet& m, const DJet& c)
{
  check(DVector(c));
  check(BMatrix(m),DMatrix(c));

  check(m.begin(0,2),c.begin(0,2),c.end(0,2),hessian);
  check(m.begin(0,2),c.begin(1,2),c.end(1,2),hessian+6);
  check(m.begin(0,2),c.begin(2,2),c.end(2,2),hessian+12);

  check(m.begin(0,3),c.begin(0,3),c.end(0,3),c3);
  check(m.begin(0,3),c.begin(1,3),c.end(1,3),c3+10);
  check(m.begin(0,3),c.begin(2,3),c.end(2,3),c3+20);

  check(m.begin(0,4),c.begin(0,4),c.end(0,4),c4);
  check(m.begin(0,4),c.begin(1,4),c.end(1,4),c4+15);
  check(m.begin(0,4),c.begin(2,4),c.end(2,4),c4+30);

  check(m.begin(0,5),c.begin(0,5),c.end(0,5),c5);
  check(m.begin(0,5),c.begin(1,5),c.end(1,5),c5+21);
  check(m.begin(0,5),c.begin(2,5),c.end(2,5),c5+42);

  check(m.begin(0,6),c.begin(0,6),c.end(0,6),c6);
  check(m.begin(0,6),c.begin(1,6),c.end(1,6),c6+28);
  check(m.begin(0,6),c.begin(2,6),c.end(2,6),c6+56);
}


void doCnTest(DCnTimeMap& tm, const BJet& mask){
  DVector u(1,3,1);
  DJet D(3,3,6);
  double t0 = 1;
  double t = t0+0.125;
  u = tm(t,u,D,t0);
  check(mask,D);
}

void doCnTest(std::string vf){
  DMap f(vf,8);
  DCnOdeSolver solver(f,15);
  solver.setStep(0.125/32);
  DCnTimeMap tm(solver);

  BJet m(3,3,6);
  solver.setMask((Multiindex*)0,(Multiindex*)0);
  doCnTest(tm,m);

  Multiindex m1[] = {{1,0,0},{0,1,0},{0,6,0},{4,2,0},{1,4,1},{4,0,2}};
  for(int j=0;j<6;++j){
    solver.addMultiindexToMask(m1[j]);
    for(int a=0;a<=m1[j][0];++a)
      for(int b=0;b<=m1[j][1];++b)
        for(int c=0;c<=m1[j][2];++c)
        {
          Multiindex p({a,b,c});
          for(int i=0;i<3;++i)
            m(i,p) = true;
        }
    doCnTest(tm,m);
  }
}


BOOST_AUTO_TEST_SUITE(CnTestSuite)

BOOST_AUTO_TEST_CASE(xCnTest)
{
  // 1. very basic test for nodes: NODE_MUL, UNARY_MINUS, NODE_ADD
  doCnTest("var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z;");
  // 2. time added
  doCnTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z;");
  // 3. TIME_PLUS_VAR, VAR_MINUS_TIME
  doCnTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),t+x*(x*x+y*y)-t,(x*x+y*y)*z;");
  return;
  // 4. TIME_MINUS_VAR, VAR_MINUS_TIME, TIME_PLUS_VAR
  doCnTest("time:t;var:x,y,z;fun:t-y*(x*x+y*y)-t,t+x*(x*x+y*y)-t,(x*x+y*y)*z;");
  // 5. TIME_PLUS_CONST, TIME_MINUS_CONST, MUL_FUNTIME_BY_FUNTIME, MUL_TIME_BY_FUNTIME, FUNTIME_MINUS_CONST
  doCnTest("time:t;var:x,y,z;fun:(t+1)*(t-1)-y*(x*x+y*y)-(t*t-1),x*(x*x+y*y),(x*x+y*y)*z;");
  // 6. TIME_PLUS_CONST, TIME_MINUS_CONST, MUL_FUNTIME_BY_FUNTIME, MUL_TIME_BY_FUNTIME, FUNTIME_MINUS_CONST, FUNTIME_MINUS_VAR, VAR_MINUS_FUNTIME
  doCnTest("time:t;var:x,y,z;fun:(t+1)*(t-1)-y*(x*x+y*y)-(t*t-1),x*(x*x+y*y),(x*x+y*y)*z;");
  // 7. ADDITIONALY FUNTIME_MINUS_FUNTIME
  doCnTest("time:t;var:x,y,z;fun:((t+1)*(t-1)-(t*t-1))-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z;");
  // 8. NODE_MUL_CONST_BY_CONST
  doCnTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),(2*0.5)*x*(x*x+y*y),(x*x+y*y)*z;");
  // 9. NODE_MUL_CONST_BY_VAR
  doCnTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),2*x*(x*x+y*y)*0.5,(x*x+y*y)*z;");
/*
  // ------------------------- DIV ----------------------------------------------
  // 10. DIV
  doCnTest("time:t;var:x,y,z;fun:-y*(1/y)*y*(x*x+y*y),x*(x*x+y*y)*((x*x+y*y+1)/(x*x+y*y+10))*((x*x+y*y+10)/(x*x+y*y+1)),(x*x+y*y)*z;");
  // 11. DIV_VAR_BY_CONST
  doCnTest("time:t;var:x,y,z;fun:-y*(1/y)*y*(x*x+y*y),4*x*(x*x+y*y)/2/2,(x*x+y*y)*z;");
  // 12. DIV_VAR_BY_FUNTIME
  doCnTest("time:t;var:x,y,z;fun:-sin(2*t)*y*(x*x+y*y)/(2*sin(t)*cos(t)),x*(x*x+y*y),(t*t-0.25)*(x*x+y*y)*z/(t-0.5)/(t+0.5);");
  // 13. DIV_VAR_BY_TIME
  doCnTest("time:t;var:x,y,z;fun:-y*(1/y)*y*(x*x+y*y)+(x+y)/t-x/t-y/t,t*x*(x*x+y*y)/t,(x*x+y*y)*z;");
  // 14. DIV_TIME_BY_CONST
  doCnTest("time:t;var:x,y,z;fun:-(t/3)*y*(x*x+y*y)*3/t,x*(x*x+y*y),(x*x+y*y)*z;");
  // 15. DIV_FUNTIME_BY_CONST
  doCnTest("time:t;var:x,y,z;fun:-(sin(2*t)/2)*y*(x*x+y*y)/sin(t)/cos(t),x*(x*x+y*y),(x*x+y*y)*z;");
  // 16. DIV_FUNTIME_BY_TIME
  doCnTest("time:t;var:x,y,z;fun:-(sin(t)/t)*y*(x*x+y*y)*t/sin(t),x*(x*x+y*y),(x*x+y*y)*z;");
  // 17. DIV_FUNTIME_BY_FUNTIME
  doCnTest("time:t;var:x,y,z;fun:-exp(2*sin(t))/exp(sin(t))/exp(sin(t))*y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z;");
  // 18. DIV_CONST_BY_CONST
  doCnTest("time:t;var:x,y,z;fun:-(7/5)*y*(x*x+y*y)*(5/7),x*(x*x+y*y),(x*x+y*y)*z;");

    // ------------------------- EXP-LOG ----------------------------------------------
  // 21. EXP
  doCnTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+exp(x+y)-exp(x)*exp(y);");
  // 22. EXP_CONST
  doCnTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+exp(6)-exp(3)*exp(3);");
  // 23. EXP_TIME
  doCnTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+sqrt(exp(t)*exp(t))-exp(t);");
  // 24. EXP_FUNTIME
  doCnTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+exp(t*t-t)-exp(t^2)/exp(t);");
  // 25. LOG
  doCnTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+log((x*x+y*y)*(x*x+z*z))-log(x*x+y*y)-log(x*x+z*z);");
  // 26. LOG_CONST
  doCnTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+log(10)-log(2)-log(5);");
  // 27. LOG_TIME
  doCnTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),sqr(log(t))*(x*x+y*y)*z/log(t)/log(t);");
  // 28. LOG_FUNTIME
  doCnTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+log((t+1)^2)-2*log(t+1);");
*/
  // ------------------------- SQR ----------------------------------------------
  // 31. SQR
  doCnTest("var:x,y,z;fun:-y*(x^2+y^2),x*(x^2+y^2),(x^2+y^2)*z;");
  // 32. SQR_CONST
  doCnTest("var:x,y,z;fun:3^2-y*(x^2+y^2)-9,x*(x^2+y^2),(x^2+y^2)*z;");
  // 33. SQR_TIME
  doCnTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(t^2-t*t)+(x*x+y*y)*z;");
  // 34. SQR_FUNTIME
  doCnTest("time:t;var:x,y,z;fun:(t+1)^2-y*(x^2+y^2)-t^2-2*t-1,x*(x^2+y^2),(x^2+y^2)*z;");

  // ------------------------- SQRT ----------------------------------------------
  // 35. SQRT
  doCnTest("var:x,y,z;fun:sqrt((x*x+2*x+y*y*y*y+1)*(x*x+2*x+y*y*y*y+1))-(x+1)*(x+1)-(y*y)*(y*y)-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z;");
  // 36. SQRT_CONST
  doCnTest("var:x,y,z;fun:3-y*(x^2+y^2)-sqrt(9),x*(x^2+y^2),(x^2+y^2)*z;");
  // 37. SQRT_TIME
  doCnTest("time:t;var:x,y,z;fun:(sqrt(t)+0.5)/(sqrt(t)-0.5)-(t+sqrt(t)+0.25)/(t-0.25)-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z;");
  // 38. SQRT_FUNTIME
  doCnTest("time:t;var:x,y,z;fun:t+1-sqrt(t^2+2*t+1)-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z;");

  // ------------------------- POW ----------------------------------------------
  // 41. POW
  doCnTest("var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+(sqr(1/sqr(5+x*x+y*y)))^0.25-1/(5+x*x+y*y);");
  // 42. POW_CONST
  doCnTest("var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+(256^0.125)-2;");
  // 43. POW_TIME
  doCnTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+sqr(sqr(t^.25))-t;");
  // 44. POW_FUNTIME
  doCnTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),t+1-(t*t*t*t+4*t*t*t+6*t*t+4*t+1)^0.25+x*(x*x+y*y),(x*x+y*y)*z;");
/*
  // ------------------------- SIN-COS ----------------------------------------------
  // 51. SIN/COS
  doCnTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+sin(x+y)-sin(x)*cos(y)-sin(y)*cos(x);");
  // 52. SIN/COS_CONST
  doCnTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+sin(3+4)-sin(3)*cos(4)-sin(4)*cos(3);");
  // 53. SIN/COS_TIME
  doCnTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+sin(t)^2+cos(t)^2-1;");
  // 54. SIN/COS_FUNTIME
  doCnTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+sin(2*t)-2*sin(t)*cos(t);");

  // ------------------------- ATAN ----------------------------------------------
  // 61. ATAN
  doCnTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+sin(atan(x*y+1))/cos(atan(x*y+1))-x*y-1;");
  // 62. ATAN_CONST
  doCnTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+2*sin(atan(1))-sqrt(2);");
  // 63. ATAN_TIME
  doCnTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+sin(atan(t))/cos(atan(t))-t;");
  // 64. ATAN_FUNTIME
  doCnTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+atan(t)-2*atan(t/(1+sqrt(1+t*t)));");

  // ------------------------- ASIN ----------------------------------------------

  // 65. ASIN
  doCnTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+asin(1/(x^2+2*y^2+2))-2*atan((1/(x^2+2*y^2+2))/(1+sqrt(1-1/(x^2+2*y^2+2)^2)));");
  // 66. ASIN_CONST
  doCnTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+asin(1./sqrt(2.))-atan(1);");
  // 67. ASIN_TIME
  doCnTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+asin(t-1)-2*atan((t-1)/(1+sqrt(1-(t-1)^2)));");
  // 68. ASIN_FUNTIME
  doCnTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+atan(t)-asin(t/sqrt(1+t^2));");

  // ------------------------- ACOS ----------------------------------------------

  // 69. ACOS
  doCnTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+acos(1/(x^2+2*y^2+2))-2*atan(sqrt(1-1/(x^2+2*y^2+2)^2)/(1+1/(x^2+2*y^2+2)));");
  // 70. ACOS_CONST
  doCnTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+acos(0.5)-atan(sqrt(3));");
  // 71. ACOS_TIME
  doCnTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+acos(t-1)-2*atan(sqrt(1-(t-1)^2)/t);");
  // 72. ACOS_FUNTIME
  doCnTest("time:t;var:x,y,z;fun:-y*(x*x+y*y),x*(x*x+y*y),(x*x+y*y)*z+atan(t)-acos(1/sqrt(1+t^2));");
*/
}

BOOST_AUTO_TEST_SUITE_END()
