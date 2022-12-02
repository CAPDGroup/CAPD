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

//
//#include <boost/test/test_case_template.hpp>
#include <boost/test/unit_test.hpp>

// namespace bl = boost::lambda;
// namespace b = boost;

using namespace capd;
using namespace capd::vectalg;

BOOST_AUTO_TEST_SUITE(C11StaticArrayInitTestSuite)

BOOST_AUTO_TEST_CASE(placeholder) {}

BOOST_AUTO_TEST_CASE(vectorInit)
{
   {
      Vector<int, 0> v1 = {1, 2, 3};

      BOOST_CHECK_EQUAL(3u, v1.dimension());
      BOOST_CHECK_EQUAL(1, v1[0]);
      BOOST_CHECK_EQUAL(2, v1[1]);
      BOOST_CHECK_EQUAL(3, v1[2]);
   }

   {
      Vector<int, 0> v1({1, 2, 3});

      BOOST_CHECK_EQUAL(3u, v1.dimension());
      BOOST_CHECK_EQUAL(1, v1[0]);
      BOOST_CHECK_EQUAL(2, v1[1]);
      BOOST_CHECK_EQUAL(3, v1[2]);
   }
   {
      Vector<int, 3> v1({1, 2, 3});

      BOOST_CHECK_EQUAL(3u, v1.dimension());
      BOOST_CHECK_EQUAL(1, v1[0]);
      BOOST_CHECK_EQUAL(2, v1[1]);
      BOOST_CHECK_EQUAL(3, v1[2]);
   }

   {
      typedef Vector<int, 2> V2;
      BOOST_CHECK_THROW(V2 v1({1, 2, 3}), std::range_error);
   }
}

BOOST_AUTO_TEST_CASE(matrixInit)
{
   {
      //    int data[][2] = { {1, 2}, {4, 5} };
      Matrix<int, 0, 0> v1({{1, 2}, {4, 5}});

      BOOST_CHECK_EQUAL(2u, v1.numberOfRows());
      BOOST_CHECK_EQUAL(2u, v1.numberOfColumns());
      BOOST_CHECK_EQUAL(1, v1[0][0]);
      BOOST_CHECK_EQUAL(2, v1[0][1]);
      BOOST_CHECK_EQUAL(4, v1[1][0]);
      BOOST_CHECK_EQUAL(5, v1[1][1]);
   }
   {
      //    int data[][2] = { {1, 2}, {4, 5} };
      Matrix<double, 0, 0> v1({{1.0, 2.0}, {4.0, 5.0}});

      BOOST_CHECK_EQUAL(2u, v1.numberOfRows());
      BOOST_CHECK_EQUAL(2u, v1.numberOfColumns());
      BOOST_CHECK_EQUAL(1., v1[0][0]);
      BOOST_CHECK_EQUAL(2., v1[0][1]);
      BOOST_CHECK_EQUAL(4., v1[1][0]);
      BOOST_CHECK_EQUAL(5., v1[1][1]);
   }
   {
      Matrix<int, 0, 0> v1({{1}, {4}});
      BOOST_CHECK_EQUAL(2u, v1.numberOfRows());
      BOOST_CHECK_EQUAL(1u, v1.numberOfColumns());
      BOOST_CHECK_EQUAL(1, v1[0][0]);
      BOOST_CHECK_EQUAL(4, v1[1][0]);
   }
}
BOOST_AUTO_TEST_CASE(matrixInitThrow)
{
   {
      typedef Matrix<int, 0, 0> MyZMatrix;
      BOOST_CHECK_THROW(MyZMatrix A1({{1, 2}, {3}}), std::range_error);
      BOOST_CHECK_THROW(MyZMatrix A2({{1}, {3, 4}}), std::runtime_error);
   	  BOOST_CHECK_THROW(MyZMatrix A2({{1}, {3, 4}}), std::range_error);
	}
}

// typedef Matrix<int,0,0> MyMatrix;
// int * entry = 0;

// MyMatrix foo(){
//  MyMatrix B({{1,3,4},{3,4,5}});
//  entry = &B[0][0];
//  return B;
//}
//
BOOST_AUTO_TEST_CASE(matrixMove)
{
   {
      typedef Matrix<int, 0, 0> MyMatrix;
      MyMatrix A = {{1, 2, 4}, {3, 4, 5}};
      MyMatrix A2(std::move(A));
      BOOST_CHECK_EQUAL(2u, A2.numberOfRows());
      BOOST_CHECK_EQUAL(3u, A2.numberOfColumns());
      BOOST_CHECK_EQUAL(1, A2(1, 1));
      BOOST_CHECK_EQUAL(2, A2(1, 2));
      BOOST_CHECK_EQUAL(4, A2(1, 3));
      BOOST_CHECK_EQUAL(3, A2(2, 1));
      BOOST_CHECK_EQUAL(4, A2(2, 2));
      BOOST_CHECK_EQUAL(5, A2(2, 3));
      BOOST_CHECK_EQUAL(0, (long)&A[0][0]);
      //    BOOST_CHECK_EQUAL(0, A.numberOfRows());
      //    BOOST_CHECK_EQUAL(0, A.numberOfColumns());
      //    std::cout << A.numberOfColumns() << "\n";
   }
   {
      typedef Matrix<int, 0, 0> MyMatrix;
      MyMatrix A = {{1, 2, 4}, {3, 4, 5}};
      MyMatrix B = {{-1, -2}, {-3, -4}};
      //    std::cout << "B c = " <<  B.numberOfColumns() << " r= " << B.numberOfRows() << "\n";
      //    std::cout << A << "\n" << B << "\n=========\n";
      B = std::move(A);
      BOOST_CHECK_EQUAL(2u, B.numberOfRows());
      BOOST_CHECK_EQUAL(3u, B.numberOfColumns());
      BOOST_CHECK_EQUAL(1, B(1, 1));
      BOOST_CHECK_EQUAL(2, B(1, 2));
      BOOST_CHECK_EQUAL(4, B(1, 3));
      BOOST_CHECK_EQUAL(3, B(2, 1));
      BOOST_CHECK_EQUAL(4, B(2, 2));
      BOOST_CHECK_EQUAL(5, B(2, 3));

      BOOST_CHECK_EQUAL(2u, A.numberOfRows());
      BOOST_CHECK_EQUAL(2u, A.numberOfColumns());
      BOOST_CHECK_EQUAL(-1, A(1, 1));
      BOOST_CHECK_EQUAL(-2, A(1, 2));
      BOOST_CHECK_EQUAL(-3, A(2, 1));
      BOOST_CHECK_EQUAL(-4, A(2, 2));

      //  std::cout << A << "\n" << B << "\n";
      //    std::cout << "A c = " <<  A.numberOfColumns() << " r= " << A.numberOfRows() << "\n";
      //    std::cout << "B c = " <<  B.numberOfColumns() << " r= " << B.numberOfRows() << "\n";
   }
}

BOOST_AUTO_TEST_SUITE_END()
