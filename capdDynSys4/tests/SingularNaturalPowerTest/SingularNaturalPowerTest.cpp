#define BOOST_TEST_MODULE SingularNaturalPowTest
#define BOOST_TEST_DYN_LINK

#include <iostream>
#include <boost/test/unit_test.hpp>
#include "capd/capdlib.h"
using namespace capd;

constexpr long double tol = 1e-13;

template<class E, class I>
void check(I b, I e, E r){
  for(;b!=e;++b,++r)
    BOOST_CHECK_SMALL(*b-*r,tol);
}

struct TestData{
  TestData(std::string vf, int order, LDVector x, int degree=1) 
    : f(vf,degree), order(order) , degree(degree)
  {
    alloc(x.dimension());
  }

  void setOrder(int order){
    LDVector x = c[0];
    free();
    this->order = order;
    alloc(x.dimension());
    init(x);
  }
  
  void init(LDVector x){
    c[0] = x;
    c1[0] = x;
    j[0]() = x;    
  }
  
  void alloc(int dim){
    f.setOrder(order);
    f.setCurrentTime(0.0);
    c = new LDVector[order+1];
    c1 = new LDVector[order+1];
    m = new LDMatrix[order+1];
    j = new LDJet[order+1];
    for(int i=0;i<=order;++i){
      c[i] = LDVector(dim);
      c1[i] = LDVector(dim);
      m[i] = LDMatrix::Identity(dim);
      j[i] = LDJet(dim,dim,degree);
    }
    m[0] = LDMatrix::Identity(dim);
    f.computeODECoefficients(c,order);
    f.computeODECoefficients(c1,m,order);
    j[0].setMatrix(LDMatrix::Identity(dim));    
  }
  
  void free(){
    delete[] c; 
    delete[] c1;
    delete[] m;
    delete[] j;    
  }
  
  ~TestData() { 
    free();
  }
  
  LDMap f;
  LDVector* c;
  LDVector* c1;
  LDMatrix* m;
  LDJet* j;
  int order, degree;
};

void makeC0Test(std::string f1, std::string f2, int order, LDVector u){
  TestData td1(f1,order,u);
  TestData td2(f2,order,u);
  for(int i=0;i<=order;++i){
    std::cout << "i:" << i << " " << td1.c[i] << " " << td2.c[i] << std::endl;
    check(td1.c[i].begin(),td1.c[i].end(),td2.c[i].begin());
  }
  std::cout << std::endl;

  // change order
  order += 5;
  td1.setOrder(order);
  td2.setOrder(order);
  for(int i=0;i<=order;++i){
    std::cout << "i:" << i << " " << td1.c[i] << " " << td2.c[i] << std::endl;
    check(td1.c[i].begin(),td1.c[i].end(),td2.c[i].begin());
  }
  std::cout << std::endl;
}

void makeC1Test(std::string f1, std::string f2, int order, LDVector u){
  TestData td1(f1,order,u);
  TestData td2(f2,order,u);
  for(int i=0;i<=order;++i){
    std::cout << "i:" << i << " " << td1.c[i] << " " << td2.c[i] << std::endl;
    std::cout << "i:" << i << " " << td1.m[i] << " " << td2.m[i] << std::endl;
    check(td1.c1[i].begin(),td1.c1[i].end(),td2.c1[i].begin());
    check(td1.m[i].begin(),td1.m[i].end(),td2.m[i].begin());
  }
  std::cout << std::endl;

  // change order
  order += 5;
  td1.setOrder(order);
  td2.setOrder(order);
  for(int i=0;i<=order;++i){
    std::cout << "i:" << i << " " << td1.c[i] << " " << td2.c[i] << std::endl;
    std::cout << "i:" << i << " " << td1.m[i] << " " << td2.m[i] << std::endl;
    check(td1.c1[i].begin(),td1.c1[i].end(),td2.c1[i].begin());
    check(td1.m[i].begin(),td1.m[i].end(),td2.m[i].begin());
  }
  std::cout << std::endl;
}

void makeCnTest(std::string f1, std::string f2, int order, LDVector u){
  int degree = 10;
  int dim = u.dimension();
  TestData td1(f1,order,u,degree);
  TestData td2(f2,order,u,degree);

  td1.f.computeODECoefficients(td1.j,degree,order/2);
  td2.f.computeODECoefficients(td2.j,degree,order/2);
  
  for(int i=0;i<=order;++i)
  {
    check(td1.j[i].begin(),td1.j[i].end(),td2.j[i].begin());
  }

  td1.f.computeODECoefficients(td1.j,degree,order);
  td2.f.computeODECoefficients(td2.j,degree,order);
  
  for(int i=0;i<=order;++i)
  {
    check(td1.j[i].begin(),td1.j[i].end(),td2.j[i].begin());
  }
  std::cout << std::endl;

  // change order
  order += 5;
  td1.setOrder(order);
  td2.setOrder(order);

  td1.f.computeODECoefficients(td1.j,degree,order/2);
  td2.f.computeODECoefficients(td2.j,degree,order/2);
  
  for(int i=0;i<=order;++i)
  {
    check(td1.j[i].begin(),td1.j[i].end(),td2.j[i].begin());
  }

  td1.f.computeODECoefficients(td1.j,degree,order);
  td2.f.computeODECoefficients(td2.j,degree,order);
  
  for(int i=0;i<=order;++i)
  {
    check(td1.j[i].begin(),td1.j[i].end(),td2.j[i].begin());
  }
  std::cout << std::endl;
}
// --------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(TestSuite)

BOOST_AUTO_TEST_CASE(xC0TimeTest)
{
  std::cout << "xC0TimeTest" << std::endl;
  // here we test Power of t^c
  makeC0Test("var:x,y;time:t;fun:sin(y)*t^11,cos(x)*t^16;",
             "var:x,y;time:t;fun:sin(y)*t^2*t^2*t^2*t^2*t^3,cos(x)*t^2*t^2*t^2*t^2*t^2*t^3*t^3;",
             20, {1.L,1.L}
             );
}

BOOST_AUTO_TEST_CASE(xC0FunTimeTest)
{
  std::cout << "xC0FunTimeTest" << std::endl;
  // here we test Power of f(t)^c, only
  makeC0Test("var:x,y;time:t;fun:sin(y)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t),cos(x)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t);",
             "var:x,y;time:t;fun:sin(y)*(((sin(t))^2)^2)^2*(sin(t))^2*sin(t),cos(x)*((((sin(t))^2)^2)^2)^2;",
             20, {1.L,1.L}
             );
  makeC0Test("var:x,y;time:t;fun:sin(y)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t),cos(x)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t);",
             "var:x,y;time:t;fun:sin(y)*(sin(t))^11,cos(x)*((((sin(t))^2)^2)^2)^2;",
             20, {1.L,1.L}
             );
  makeC0Test("var:x,y;time:t;fun:sin(y)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t),cos(x)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t);",
             "var:x,y;time:t;fun:sin(y)*(sin(t))^11,cos(x)*(sin(t))^16;",
             20, {1.L,1.L}
             );
}

BOOST_AUTO_TEST_CASE(xC0Test)
{
  std::cout << "xC0Test" << std::endl;
  // here we test Power of t^c, only
  makeC0Test("var:x,y;time:t;fun:(1+y*(1+y*(1+y*(1+y*(1+y*(1+y*(1+y)))))))+sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t),cos(x)-sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t);",
             "var:x,y;time:t;fun:1+y^7+y^6+y^5+y^4+y^3+y^2+y+(sin(t))^11,cos(x)-(sin(t))^16;",
             20, {0.L,0.L}
             );
}

BOOST_AUTO_TEST_CASE(xC1Test)
{
  std::cout << "xC1Test" << std::endl;
  // here we test Power of t^c, only
  makeC1Test("var:x,y;time:t;fun:(1+y*(1+y*(1+y*(1+y*(1+y*(1+y*(1+y)))))))+sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t),cos(x)-sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t);",
             "var:x,y;time:t;fun:1+y^7+y^6+y^5+y^4+y^3+y^2+y+(sin(t))^11,cos(x)-(sin(t))^16;",
             20, {0.L,0.L}
             );
}

BOOST_AUTO_TEST_CASE(xCnTest)
{
  std::cout << "xCnTest" << std::endl;
  // here we test Power of t^c, only
  makeCnTest("var:x,y;time:t;fun:(1+y*(1+y*(1+y*(1+y*(1+y*(1+y*(1+y)))))))+sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t),cos(x)-sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t)*sin(t);",
             "var:x,y;time:t;fun:1+y^7+y^6+y^5+y^4+y^3+y^2+y+(sin(t))^11,cos(x)-(sin(t))^16;",
             10, {0.L,0.L}
             );
}

BOOST_AUTO_TEST_SUITE_END()
