#define BOOST_TEST_MODULE DPoincareMapMaskTest
#define BOOST_TEST_DYN_LINK

#include <iostream>
#include <boost/test/unit_test.hpp>
#include "capd/capdlib.h"
#include "capd/normalForms/planarMaps.hpp"

#include <chrono>
std::chrono::time_point<std::chrono::system_clock> tic(){
  auto start = std::chrono::system_clock::now();
  return start;
}

void tac(std::chrono::time_point<std::chrono::system_clock>& start){
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::cout << "elapsed time: " << elapsed_seconds.count() << " seconds\n";
}


#define CAPD_USER_NAMESPACE capd2
#define CAPD_DEFAULT_DIMENSION 2
#include "capd/fdcapdlib.h"
#undef CAPD_USER_NAMESPACE
#undef CAPD_DEFAULT_DIMENSION

using namespace std;
using namespace capd;
template class capd::map::Map<DMatrix>;
template class capd::dynsys::BasicCnOdeSolver<DMap>;
template class capd::poincare::BasicPoincareMap<DCnOdeSolver>;

// global data
const double tol = 1e-12;
const DVector x({0.,1.55,0.05});
const DMatrix R({{-1,0,0},{0,1,0},{0,0,-1}});
const DMatrix I({{0,0,0},{0,1,0},{0,0,1}});
DVector x1(3), x2(3);
DMatrix D1(3,3), D2(3,3), D(3,3);
DHessian H1(3,3), H2(3,3), H(3,3);
DJet J1(3,3,3), J2(3,3,3), J(3,3,3);

template<class E, class O>
void check(O b, O e, E r){
  for(;b!=e;++b,++r)
    BOOST_CHECK_SMALL(*b-*r,tol);
}

template<class I>
void check(I b, I e){
  for(;b!=e;++b)
    BOOST_CHECK_EQUAL(*b,0.);
}

// --------------------------------------------------------------

struct C1Test{
  template<class PM>
  static void computePM(PM& pm){
    x1 = pm(R*x,D1);
    D1 = pm.computeDP(x1,D1);
    x2 = pm(R*x1,D2);
    D2 = pm.computeDP(x2,D2);
    D = D2*R*D1*R;
  }

  static void checkNoMask(){
    check(x2.begin(),x2.end(),x.begin());
    for(int i=1;i<3;++i)
      check(D.column(i).begin(),D.column(i).end(),I.column(i).begin());
  }

  static void checkMask(){
    checkNoMask();
    check(D.begin(),D.end(),I.begin());
    check(D.column(0).begin(),D.column(0).end());
  }

  template<class Solver>
  static void setMask(Solver& solver){
    Multiindex i[] = {{0,1,0},{0,0,1}};
    solver.setMask(i,i+2);
  }
};

// --------------------------------------------------------------

struct C2Test{
  template<class PM>
  static void computePM(PM& pm){
    DMatrix m(3,3);
    DHessian h(3,3);
    x1 = pm(R*x,m,h);
    pm.computeDP(x1,m,h,D1,H1);
    x2 = pm(R*x1,m,h);
    pm.computeDP(x2,m,h,D2,H2);
    D = D2*R*D1*R;
    H = (H2*R)*(D1*R) + (D2*R)*(H1*R);
  }

  template<class Solver>
  static void setMask(Solver& solver){
    Multiindex i[] = {{0,2,0},{0,1,1},{0,0,2}};
    solver.setMask(i,i+3);
  }

  static void checkNoMask(){
    C1Test::checkNoMask();
    check(H(1,1).begin(),H(1,1).end(),I.column(0).begin());
    check(H(1,2).begin(),H(1,2).end(),I.column(0).begin());
    check(H(2,2).begin(),H(2,2).end(),I.column(0).begin());
  }

  static void checkMask(){
    checkNoMask();
    check(D.begin(),D.end(),I.begin());
    check(H(0,0).begin(),H(0,0).end());
    check(H(0,1).begin(),H(0,1).end());
    check(H(0,2).begin(),H(0,2).end());
  }
};

// --------------------------------------------------------------

struct CnTest{
  template<class PM>
  static void computePM(PM& pm){
    J1.clear();
    J2.clear();

    J1() = x;
    J1.setMatrix(DMatrix::Identity(3));
    x1 = R*pm(J1);
    J1 = R*pm.computeDP(J1);

    J2() = x1;
    J2.setMatrix(DMatrix::Identity(3));
    x2 = R*pm(J2);
    J2 = R*pm.computeDP(J2);
    substitutionPowerSeries(J2,J1,J,false);
    D = DMatrix(J);
    H = DHessian(J);
  }

  template<class Solver>
  static void setMask(Solver& solver){
    Multiindex i[] = {{0,3,0},{0,2,1},{0,1,2},{0,0,3}};
    solver.setMask(i,i+4);
  }

  static void checkNoMask(){
    C2Test::checkNoMask();
    Multiindex m[] = {{0,3,0},{0,2,1},{0,1,2},{0,0,3}};
    for(auto i : m)
      check(J(i).begin(),J(i).end(),I.column(0).begin());
  }

  static void checkMask(){
    checkNoMask();
    C2Test::checkMask();
    Multiindex m[] = {{3,0,0},{2,1,0},{2,0,1},{1,2,0},{1,1,1},{1,0,2}};
    for(auto i : m)
      check(J(i).begin(),J(i).end());
  }
};

// --------------------------------------------------------------

template<class PM,class Test>
void doTest(){
  DMap f("var:x,y,z;fun:y,z,1-0.5*x^2-y;",5);
  typename PM::Solver solver(f,15);
  solver.setStep(1e-2);
  DCoordinateSection section(3,0);
  PM pm(solver,section);

  Test::computePM(pm);
  Test::checkNoMask();

  Test::setMask(solver);
  Test::computePM(pm);
  Test::checkMask();
}

// --------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(TestSuite)

BOOST_AUTO_TEST_CASE(xC1Test)
{
  // test based on reversibility of the system
  doTest<DPoincareMap,  C1Test>();
  doTest<DC2PoincareMap,C1Test>();
  doTest<DCnPoincareMap,C1Test>();
}

BOOST_AUTO_TEST_CASE(xC2Test)
{
  // test based on reversibility of the system
  doTest<DC2PoincareMap,C2Test>();
  doTest<DCnPoincareMap,C2Test>();
}

BOOST_AUTO_TEST_CASE(xCnTest)
{
  // test based on reversibility of the system
  doTest<DCnPoincareMap,CnTest>();
}

void normalFormTest(LDCnPoincareMap& pm){
  // numerical approximation of an elliptic periodic point
  LDVector x({0.,0.24464743045793616028,0.});
  const int degree = 5;
  const int dim = 3;
  LDJet c(dim,dim,degree);
  c() = x;
  c.setMatrix(LDMatrix::Identity(dim));
  auto t = tic();
  pm(c);
  pm(c);
  LDJet DP = pm.computeDP(c);
  tac(t);

  capd2::DJet d(2,2,degree);
  Multiindex mi(2), mi2(dim);
  mi2[0]=0;
  for(int r=1;r<=degree;++r)
  {
     for(int j=0;j<=r;++j)
     {
        mi2[1]=mi[0]=j;
        mi2[2]=mi[1]=r-j;
        d(0,mi) = DP(1,mi2);
        d(1,mi) = DP(2,mi2);
     }
  }

  // Here we compute normal form coefficients and check if the imaginary part is close to zero.
  capd::vectalg::Vector<std::complex<double>,0> nForm = capd::normalForms::computePlanarEllipticNormalForm(d);
  std::cout << nForm[1] << ", " << nForm[2] << std::endl;
  BOOST_CHECK_SMALL(nForm[1].imag(),1e-14);
  BOOST_CHECK_SMALL(nForm[2].imag(),2e-12);
}

BOOST_AUTO_TEST_CASE(xNormalFormTest)
{
  // Birkhoff Normal Form at an elliptic periodic point is computed using complex numbers.
  // Form the theory it is known that the coefficients are reals.
  // We check if the imaginary part is almost zero.
  LDMap f("par:c;var:x,y,z;fun:y,z,c^2-0.5*x^2-y;",5);
  f.setParameter("c",0.125);
  LDCnOdeSolver solver(f,10);
  LDCoordinateSection section(3,0);
  LDCnPoincareMap pm(solver,section);

  normalFormTest(pm);
  Multiindex m[] = {{0,5,0},{0,4,1},{0,3,2},{0,2,3},{0,1,4},{0,0,5}};
  solver.setMask(m,m+6);
  normalFormTest(pm);
}

BOOST_AUTO_TEST_SUITE_END()
