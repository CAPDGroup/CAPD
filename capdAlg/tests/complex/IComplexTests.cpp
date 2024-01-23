
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

BOOST_AUTO_TEST_SUITE_END()