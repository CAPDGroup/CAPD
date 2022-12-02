#define BOOST_TEST_MODULE DComplexTest
#define BOOST_TEST_DYN_LINK
//#include <boost/test/test_case_template.hpp>
#include <boost/test/unit_test.hpp>

#include "capd/fields/Complex.h"

using DComplex = capd::fields::Complex<double>;

BOOST_AUTO_TEST_SUITE(DComplexTest)

//BOOST_AUTO_TEST_CASE(placeholder) {}

BOOST_AUTO_TEST_CASE(constructorsTest)
{
  
	DComplex a;
	BOOST_CHECK_EQUAL(a.real(), 0.0);
	BOOST_CHECK_EQUAL(a.imag(), 0.0);
	DComplex b(1.0);
	BOOST_CHECK_EQUAL(b.real(), 1.0);
	BOOST_CHECK_EQUAL(b.imag(), 0.0);
	DComplex  c(2.0, 3.0);
	BOOST_CHECK_EQUAL(c.real(), 2.0);
	BOOST_CHECK_EQUAL(c.imag(), 3.0);
	DComplex  d(c);
	BOOST_CHECK_EQUAL(d.real(), 2.0);
	BOOST_CHECK_EQUAL(d.imag(), 3.0);
}

  BOOST_AUTO_TEST_CASE( operatorsTest)
  {
	double a{1.0},
		b{2.0},
		c{3.0},
		d{-4.0};

	DComplex  x(a, b);
	DComplex  y(c, d);

	{
	  DComplex z = x+y;
	  BOOST_CHECK_EQUAL(z.real(), a + c);
	  BOOST_CHECK_EQUAL(z.imag(), b + d);
	}

	{
	  DComplex z = x-y;
	  BOOST_CHECK_EQUAL(z.real(), a - c);
	  BOOST_CHECK_EQUAL(z.imag(), b - d);
	}

	{
	  DComplex z = x*y;
	  BOOST_CHECK_EQUAL(z.real(), a*c-b*d);
	  BOOST_CHECK_EQUAL(z.imag(), a*d + b*c);
	}
	{
	  DComplex z = x/y;
//	BOOST_CHECK_EQUAL(z.real(), a*c-b*d);
//	BOOST_CHECK_EQUAL(z.imag(), a*d + b*c);
	}

  }

  BOOST_AUTO_TEST_CASE( operatorsTest2)
  {
	double a{4.0},
		b{0.0},
		c{2.0},
		d{0.0};

	DComplex  x(a, b);
	DComplex  y(c, d);

	{
	  DComplex z = x+y;
	  BOOST_CHECK_EQUAL(z.real(), a + c);
	  BOOST_CHECK_EQUAL(z.imag(), b + d);
	}

	{
	  DComplex z = x-y;
	  BOOST_CHECK_EQUAL(z.real(), a - c);
	  BOOST_CHECK_EQUAL(z.imag(), b - d);
	}

	{
	  DComplex z = x*y;
	  BOOST_CHECK_EQUAL(z.real(), a*c-b*d);
	  BOOST_CHECK_EQUAL(z.imag(), a*d + b*c);
	}
	{
	  DComplex z = x/y;
	  BOOST_CHECK_EQUAL(z.real() ,2.0 );
	  BOOST_CHECK_EQUAL(z.imag(), 0.0);
	}

  }
  BOOST_AUTO_TEST_CASE( functionsTest)
  {
  }



BOOST_AUTO_TEST_SUITE_END() 