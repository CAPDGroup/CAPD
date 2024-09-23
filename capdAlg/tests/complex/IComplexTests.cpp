
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#define BOOST_TEST_MODULE IComplexTest
#define BOOST_TEST_DYN_LINK
//#include <boost/test/test_case_template.hpp>
#include <boost/test/unit_test.hpp>

//#include "capd/rounding/DoubleRounding.h"
//#include "capd/intervals/DoubleInterval.h"
//#include "capd/intervals/IntervalError.h"
#include "capd/intervals/lib.h"
#include "capd/fields/Complex.h"
using namespace std;

using IComplex = capd::fields::Complex<capd::interval>;
using interval = capd::interval;

BOOST_AUTO_TEST_SUITE(IComplexTest)

BOOST_AUTO_TEST_CASE(placeholder) {}

BOOST_AUTO_TEST_CASE(constructorsTest)
{
  IComplex a;
  BOOST_CHECK_EQUAL(a.real(), interval(0.0,0.0));
  BOOST_CHECK_EQUAL(a.imag(), interval(0.0,0.0));
  IComplex b(1.0);
  BOOST_CHECK_EQUAL(b.real(), interval(1.0));
  BOOST_CHECK_EQUAL(b.imag(), interval(0.0));
  IComplex  c(2.0, 3.0);
  BOOST_CHECK_EQUAL(c.real(), interval(2.0,2.0));
  BOOST_CHECK_EQUAL(c.imag(), interval(3.0,3.0));
  IComplex  d(c);
  BOOST_CHECK_EQUAL(d.real(), interval(2.0,2.0));
  BOOST_CHECK_EQUAL(d.imag(), interval(3.0,3.0));
  IComplex  e({2.5, 3.0},{3.0, 5.0});
  BOOST_CHECK_EQUAL(e.real(), interval(2.5,3.0));
  BOOST_CHECK_EQUAL(e.imag(), interval(3.0,5.0));
}
 
BOOST_AUTO_TEST_CASE( operatorsTest)
{
  interval a{1.0, 4.0},
           b{2.0, 5.0},
		   c{2.0, 3.0},
		   d{1.0, 7.0};

  IComplex  x(a, b);
  IComplex  y(c, d);

  {
	IComplex z = x+y;
	BOOST_CHECK_EQUAL(z.real(), a + c);
	BOOST_CHECK_EQUAL(z.imag(), b + d);
  }

  {
	IComplex z = x-y;
	BOOST_CHECK_EQUAL(z.real(), a - c);
	BOOST_CHECK_EQUAL(z.imag(), b - d);
  }

  {
	IComplex z = x*y;
	BOOST_CHECK_EQUAL(z.real(), a*c-b*d);
	BOOST_CHECK_EQUAL(z.imag(), a*d + b*c);
  }
  {
	IComplex z = x/y;
//	BOOST_CHECK_EQUAL(z.real(), a*c-b*d);
//	BOOST_CHECK_EQUAL(z.imag(), a*d + b*c);
  }

}

BOOST_AUTO_TEST_CASE( operatorsTest2)
{
  interval a{2.0, 4.0},
	  b{0.0, 0.0},
	  c{2.0, 2.0},
	  d{0.0, 0.0};

  IComplex  x(a, b);
  IComplex  y(c, d);

  {
	IComplex z = x+y;
	BOOST_CHECK_EQUAL(z.real(), a + c);
	BOOST_CHECK_EQUAL(z.imag(), b + d);
  }

  {
	IComplex z = x-y;
	BOOST_CHECK_EQUAL(z.real(), a - c);
	BOOST_CHECK_EQUAL(z.imag(), b - d);
  }

  {
	IComplex z = x*y;
	BOOST_CHECK_EQUAL(z.real(), a*c-b*d);
	BOOST_CHECK_EQUAL(z.imag(), a*d + b*c);
  }
  {
	IComplex z = x/y;
    BOOST_CHECK(z.real().contains(interval{1.0,2.0}) );
    BOOST_CHECK(z.imag().contains(interval{0.0}));
  }

}
BOOST_AUTO_TEST_CASE( functionsTest)
{
}


  void test2(
    double intervalSize = 0.1,
    double stepSize =  0.1,
    interval r = interval(1.1, 1.6),
    double eps = 1.e-14
  ) {
    if(stepSize<=0) stepSize = intervalSize;
    bool testPassed = true;
    int count = 0;
    for (interval x = -interval::pi() + eps + interval(0., intervalSize); x < interval::pi() + 0.3; x += stepSize) {

      IComplex z(r * interval(cos(x)), r * interval(sin(x)));
      auto a = arg(z);

      BOOST_CHECK((a.contains(x)) or ((a + 2 * interval::pi()).contains(x)) or ((a - 2 * interval::pi()).contains(x)));
    }
  }

  BOOST_AUTO_TEST_CASE(Atan2Test)
  {
    test2(0.1);
    test2(0.1, 0.1, interval(0.000000001, 1.0));
    test2(0.1, 0.1, interval(0.1, 6.0));
    test2(0.1, 0.1, interval(1.0, 100000.0));
    test2(0.2);
  }

BOOST_AUTO_TEST_SUITE_END()