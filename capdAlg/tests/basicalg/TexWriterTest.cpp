/////////////////////////////////////////////////////////////////////////////
/// @file TexWriterTest.cpp
///
/// @author Tomasz Kapela
///
/// @date 2015-11-04
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2015 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.edu.pl/ for details.



#define BOOST_TEST_MODULE TexWriterTest
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <sstream>

#include <capd/auxil/Logger.h>
#include <capd/basicalg/TexWriter.h>

using namespace capd;
using namespace capd::basicalg;

BOOST_AUTO_TEST_SUITE(TexWriterDouble)

BOOST_AUTO_TEST_CASE(numberOfDigits)
{

  std::ostringstream sout;
  capd::TexWriter out(sout);

  BOOST_CHECK_EQUAL(3,out.numberOfDigits(112.1432453254));
  BOOST_CHECK_EQUAL(31,out.numberOfDigits(1.1432453254e30));
  BOOST_CHECK_EQUAL(256,out.numberOfDigits(1.1432453254e255));
  BOOST_CHECK_EQUAL(1,out.numberOfDigits(0.123));
  BOOST_CHECK_EQUAL(1,out.numberOfDigits(0.000123));
  BOOST_CHECK_EQUAL(1,out.numberOfDigits(1.23e-42));
  BOOST_CHECK_EQUAL(3,out.numberOfDigits(-112.1432453254));
  BOOST_CHECK_EQUAL(31,out.numberOfDigits(-1.1432453254e30));
  BOOST_CHECK_EQUAL(256,out.numberOfDigits(-1.1432453254e255));
  BOOST_CHECK_EQUAL(1,out.numberOfDigits(-0.123));
  BOOST_CHECK_EQUAL(1,out.numberOfDigits(-0.000123));
  BOOST_CHECK_EQUAL(1,out.numberOfDigits(-1.23e-42));
  BOOST_CHECK_EQUAL(1,out.numberOfDigits(0.0));

}

BOOST_AUTO_TEST_CASE(totalNumberOfDigits)
{
  {
  int precision = 10;
  int numDecDig = 0, length = 0, exponent=-1;
  double left = 1.0, right=1.0;
  capd::TexWriter::computeNumberOfDigitsForInterval(left, right, precision, TexWriter::FloatSci, numDecDig, length, exponent );
  BOOST_CHECK_EQUAL(1, numDecDig);
  BOOST_CHECK_EQUAL(11, length);
  BOOST_CHECK_EQUAL(0, exponent);
  BOOST_CHECK_EQUAL(1.0, left);
  BOOST_CHECK_EQUAL(1.0, right);
  }
  {
  int precision = 10;
  int numDecDig = 0, length = 0, exponent=-1;
  double left = 1.0, right=0.01;
  capd::TexWriter::computeNumberOfDigitsForInterval(left, right, precision, TexWriter::FloatSci, numDecDig, length, exponent );
  BOOST_CHECK_EQUAL(1, numDecDig);
  BOOST_CHECK_EQUAL(11, length);
  BOOST_CHECK_EQUAL(0, exponent);
  BOOST_CHECK_EQUAL(1.0, left);
  BOOST_CHECK_EQUAL(0.01, right);
  }
  {
  int precision = 10;
  int numDecDig = 0, length = 0, exponent=-1;
  double left = 0.01, right=1.0;
  capd::TexWriter::computeNumberOfDigitsForInterval(left, right, precision, TexWriter::FloatSci, numDecDig, length, exponent );
  BOOST_CHECK_EQUAL(1, numDecDig);
  BOOST_CHECK_EQUAL(11, length);
  BOOST_CHECK_EQUAL(0, exponent);
  BOOST_CHECK_EQUAL(0.01, left);
  BOOST_CHECK_EQUAL(1.0, right);
  }
  {
  int precision = 10;
  int numDecDig = 0, length = 0, exponent=-1;
  double left = 10., right=10.0;
  capd::TexWriter::computeNumberOfDigitsForInterval(left, right, precision, TexWriter::FloatSci, numDecDig, length, exponent );
  BOOST_CHECK_EQUAL(1, numDecDig);
  BOOST_CHECK_EQUAL(11, length);
  BOOST_CHECK_EQUAL(1, exponent);
  BOOST_CHECK_CLOSE(1.0, left, 1./precision);
  BOOST_CHECK_CLOSE(1.0, right, 1./precision);
  }
  {
  int precision = 10;
  int numDecDig = 0, length = 0, exponent=-1;
  double left = 10.0, right=1000.0;
  capd::TexWriter::computeNumberOfDigitsForInterval(left, right, precision, TexWriter::FloatSci, numDecDig, length, exponent );
  BOOST_CHECK_EQUAL(1, numDecDig);
  BOOST_CHECK_EQUAL(11, length);
  BOOST_CHECK_EQUAL(3, exponent);
  BOOST_CHECK_CLOSE( 0.01, left, 1e-13 );
  BOOST_CHECK_CLOSE(1.0, right,1e-13);
  }
  {
  int precision = 10;
  int numDecDig = 0, length = 0, exponent=-1;
  double left = 0.01, right=1.0;
  capd::TexWriter::computeNumberOfDigitsForInterval(left, right, precision, TexWriter::FloatSci, numDecDig, length, exponent );
  BOOST_CHECK_EQUAL(1, numDecDig);
  BOOST_CHECK_EQUAL(11, length);
  BOOST_CHECK_EQUAL(0, exponent);
  BOOST_CHECK_CLOSE( 0.01, left, 1e-13 );
//  BOOST_CHECK_EQUAL(0.01, left);
  BOOST_CHECK_CLOSE(1.0, right, 1e-13);
  }
  {
  int precision = 10;
  int numDecDig = 0, length = 0, exponent=-1;
  double left = -1.0, right=-1.0;
  capd::TexWriter::computeNumberOfDigitsForInterval(left, right, precision, TexWriter::FloatSci, numDecDig, length, exponent );
  BOOST_CHECK_EQUAL(1, numDecDig);
  BOOST_CHECK_EQUAL(11, length);
  BOOST_CHECK_EQUAL(0, exponent);
  BOOST_CHECK_EQUAL(-1.0, left);
  BOOST_CHECK_EQUAL(-1.0, right);
  }
  {
  int precision = 10;
  int numDecDig = 0, length = 0, exponent=-1;
  double left = -1.0, right=1.0;
  capd::TexWriter::computeNumberOfDigitsForInterval(left, right, precision, TexWriter::FloatSci, numDecDig, length, exponent );
  BOOST_CHECK_EQUAL(1, numDecDig);
  BOOST_CHECK_EQUAL(11, length);
  BOOST_CHECK_EQUAL(0, exponent);
  BOOST_CHECK_EQUAL(-1.0, left);
  BOOST_CHECK_EQUAL(1.0, right);
  }
  {
  int precision = 10;
  int numDecDig = 0, length = 0, exponent=-1;
  double left = 1.0e123, right=1.0e122;
  capd::TexWriter::computeNumberOfDigitsForInterval(left, right, precision, TexWriter::FloatSci, numDecDig, length, exponent );
  BOOST_CHECK_EQUAL(1, numDecDig);
  BOOST_CHECK_EQUAL(11, length);
  BOOST_CHECK_EQUAL(123, exponent);
  BOOST_CHECK_CLOSE( 1.0, left, 1e-13);
  BOOST_CHECK_CLOSE( 0.1, right, 1e-13 );
  }
  {
  int precision = 10;
  int numDecDig = 0, length = 0, exponent=-1;
  double left = 1.0e-124, right=1.0e-123;
  capd::TexWriter::computeNumberOfDigitsForInterval(left, right, precision, TexWriter::FloatSci, numDecDig, length, exponent );
  BOOST_CHECK_EQUAL(1, numDecDig);
  BOOST_CHECK_EQUAL(11, length);
  BOOST_CHECK_EQUAL(-123, exponent);
  BOOST_CHECK_CLOSE( 0.1, left, 1e-12);
  BOOST_CHECK_CLOSE( 1., right, 1e-12 );
  }
}


//  {
//
//    Vector<int, 0> v1 = {1, 2, 3};
//
//    BOOST_CHECK_EQUAL(3, v1.dimension());
//    BOOST_CHECK_EQUAL(1, v1[0]);
//    BOOST_CHECK_EQUAL(2, v1[1]);
//    BOOST_CHECK_EQUAL(3, v1[2]);
//  }
//
//  {
//    Vector<int, 0> v1({1, 2, 3});
//
//    BOOST_CHECK_EQUAL(3, v1.dimension());
//    BOOST_CHECK_EQUAL(1, v1[0]);
//    BOOST_CHECK_EQUAL(2, v1[1]);
//    BOOST_CHECK_EQUAL(3, v1[2]);
//  }
//
//  {
//    typedef Vector<int, 2> V2;
//    BOOST_CHECK_THROW(V2 v1({1, 2, 3}) , std::range_error);
//  }




BOOST_AUTO_TEST_SUITE_END()

//INIT_CAPD_CXX_LOGGER;

//INIT_CAPD_CXX_LOGGER;
