/////////////////////////////////////////////////////////////////////////////
/// @file StaticArrayInitTest.cpp
///
/// @author Mateusz Juda <mateusz.juda@gmail.com>
///
/// @date 2014-06-12
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD::RedHom Group.
//
// This file constitutes a part of the CAPD::RedHom library,
// distributed under the terms of the GNU General Public License.
// Consult  http://redhom.ii.edu.pl/ for details.

#include <capd/vectalg/Matrix.hpp>
#include <capd/vectalg/Vector.hpp>

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
//#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/assign.hpp>
#include <boost/foreach.hpp>
 
#include <map>
#include <vector>
#include <algorithm>
#include <utility>
#include <sstream>

namespace bl = boost::lambda;
namespace b = boost;

using namespace boost::assign;
using namespace boost::lambda;

using namespace capd;
using namespace capd::vectalg;

BOOST_AUTO_TEST_SUITE(StaticArrayInitTestSuite)

BOOST_AUTO_TEST_CASE(vectorInit)
{
  {
    int data[] = {1, 2, 3};
    Vector<int, 0> v1(data);

    BOOST_CHECK_EQUAL(3u, v1.dimension());
    BOOST_CHECK_EQUAL(1, v1[0]);
    BOOST_CHECK_EQUAL(2, v1[1]);
    BOOST_CHECK_EQUAL(3, v1[2]);
  }

  {
    Vector<int, 0> v1((int[]){1, 2, 3});

    BOOST_CHECK_EQUAL(3u, v1.dimension());
    BOOST_CHECK_EQUAL(1, v1[0]);
    BOOST_CHECK_EQUAL(2, v1[1]);
    BOOST_CHECK_EQUAL(3, v1[2]);
  }
  {
    Vector<int, 3> v1((int[]){1, 2, 3});

    BOOST_CHECK_EQUAL(3u, v1.dimension());
    BOOST_CHECK_EQUAL(1, v1[0]);
    BOOST_CHECK_EQUAL(2, v1[1]);
    BOOST_CHECK_EQUAL(3, v1[2]);
  }

  {
    Vector<int, 2> v1((int[]){1, 2, 3});

    BOOST_CHECK_EQUAL(2u, v1.dimension());
    BOOST_CHECK_EQUAL(1, v1[0]);
    BOOST_CHECK_EQUAL(2, v1[1]);
  }

}


void foo(std::string){}

BOOST_AUTO_TEST_CASE(matrixInit)
{
  {
    int data[][2] = { {1, 2}, {4, 5} };
    Matrix<int, 0, 0> v1(data);

    BOOST_CHECK_EQUAL(2u, v1.numberOfRows());
    BOOST_CHECK_EQUAL(2u, v1.numberOfColumns());
    BOOST_CHECK_EQUAL(1, v1[0][0]);
    BOOST_CHECK_EQUAL(2, v1[0][1]);
    BOOST_CHECK_EQUAL(4, v1[1][0]);
    BOOST_CHECK_EQUAL(5, v1[1][1]);
  }

  {
    int data[2][2] = { {1}, {4} };
    Matrix<int, 0, 0> v1(data);

    BOOST_CHECK_EQUAL(2u, v1.numberOfRows());
    BOOST_CHECK_EQUAL(2u, v1.numberOfColumns());
    BOOST_CHECK_EQUAL(1, v1[0][0]);
    BOOST_CHECK_EQUAL(0, v1[0][1]);
    BOOST_CHECK_EQUAL(4, v1[1][0]);
    BOOST_CHECK_EQUAL(0, v1[1][1]);
  }

}


BOOST_AUTO_TEST_SUITE_END()
