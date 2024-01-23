/////////////////////////////////////////////////////////////////////////////
/// @file smithTest.cpp
///
/// @author Mateusz Juda <mateusz.juda@gmail.com>
///
/// @date 2014-06-03
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD::RedHom Group.
//
// This file constitutes a part of the CAPD::RedHom library,
// distributed under the terms of the GNU General Public License.
// Consult  http://redhom.ii.edu.pl/ for details.



#define BOOST_TEST_DYN_LINK

#include "SmithFixtures.h"
#include "LinearEquationFixtures.h"
#include "QuotientBaseMatrixFixtures.h"

#include <capd/matrixAlgorithms/CAPDSmithForm.h>
#include <capd/matrixAlgorithms/intMatrixAlgorithms.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/assign.hpp>

#include <map>
#include <vector>
#include <algorithm>
#include <utility>
#include <capd/matrixAlgorithms/PARIIntMatrixAlgorithms.h>
#include <capd/matrixAlgorithms/PARIInterface.h>

namespace bl = boost::lambda;
namespace b = boost;

using namespace boost::assign;
using namespace boost::lambda;

using namespace capd::test;

template<class matrix>
matrix capdSmithForm(const matrix& A)
{
  matrix B(A);
  int m=B.numberOfRows();
  int n=B.numberOfColumns();
  matrix Q(m,m);
  matrix Qinv(m,m);
  matrix R(n,n);
  matrix Rinv(n,n);
  int s,t;
  capd::matrixAlgorithms::smithForm(B,Q,Qinv,R,Rinv,s,t);

  return B;
}


BOOST_AUTO_TEST_SUITE(smithTestSuite)

BOOST_AUTO_TEST_CASE_TEMPLATE(capdSmithForm1, Fixture, SmithFixtures)
{
  Fixture fixture;

  for (size_t i = 0; i < fixture.matrices_matrix.size(); ++i) {
    typename Fixture::MatrixType B = capdSmithForm(fixture.matrices_matrix[i]);
    BOOST_CHECK_EQUAL(fixture.matrices_smithMatrix[i], B);
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(smithFormMatrices, Fixture, SmithFixtures)
{
  typedef typename Fixture::MatrixType MatrixType;
  Fixture fixture;

  for (size_t caseIdx = 0; caseIdx < fixture.matrices_matrix.size(); ++caseIdx) {
    const MatrixType& matrixA = fixture.matrices_matrix[caseIdx];
    MatrixType matrixB = matrixA;
    const int m=matrixA.numberOfRows();
    const int n=matrixA.numberOfColumns();
    MatrixType matrixU(m, m), matrixV(n, n), matrixUInv(m, m), matrixVInv(n, n);
    int s, t;

    capd::matrixAlgorithms::smithForm(matrixB, matrixU, matrixUInv,
				      matrixV, matrixVInv, s, t);

    BOOST_CHECK_EQUAL(matrixB, matrixUInv * matrixA * matrixV);
    BOOST_CHECK_EQUAL(matrixU.Identity(matrixU.numberOfRows()), matrixU * matrixUInv);
    BOOST_CHECK_EQUAL(matrixU.Identity(matrixU.numberOfRows()), matrixUInv * matrixU);
    BOOST_CHECK_EQUAL(matrixV.Identity(matrixV.numberOfRows()), matrixV * matrixVInv);
    BOOST_CHECK_EQUAL(matrixV.Identity(matrixV.numberOfRows()), matrixVInv * matrixV);
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(CAPDSmithFormTest, Fixture, SmithFixtures)
{
  typedef typename Fixture::MatrixType MatrixType;
  Fixture fixture;

  for (size_t caseIdx = 0; caseIdx < fixture.matrices_matrix.size(); ++caseIdx) {
    const MatrixType& matrixA = fixture.matrices_matrix[caseIdx];
    MatrixType matrixB = matrixA;
    const int m=matrixA.numberOfRows();
    const int n=matrixA.numberOfColumns();

    capd::matrixAlgorithms::CAPDSmithForm<MatrixType> smithForm(matrixB, true, true, true, true);

    BOOST_CHECK_EQUAL(MatrixType::Identity(m), smithForm.getQ() * smithForm.getQinv());
    BOOST_CHECK_EQUAL(MatrixType::Identity(m), smithForm.getQinv() * smithForm.getQ());
    BOOST_CHECK_EQUAL(MatrixType::Identity(n), smithForm.getR() * smithForm.getRinv());
    BOOST_CHECK_EQUAL(MatrixType::Identity(n), smithForm.getRinv() * smithForm.getR());
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(capdSolveLinearEquation, Fixture, linearEquation::LinearEquationFixtures)
{
  Fixture fixture;

  for (size_t i = 0; i < fixture.data_le_matrix.size(); ++i) {
    typename Fixture::VectorType result;
    typename Fixture::VectorType& rightSide = fixture.data_le_vector[i];
    typename Fixture::MatrixType& matrix = fixture.data_le_matrix[i];
    capd::matrixAlgorithms::solveLinearEquation(matrix, rightSide, result);

    BOOST_CHECK_EQUAL(fixture.data_le_solution[i], result);
    BOOST_CHECK_EQUAL(rightSide, matrix * result);
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(capdQuotientBaseMatrix, Fixture, quotientBaseMatrix::QuotientBaseMatrixFixtures)
{
  Fixture fixture;

  for (size_t i = 0; i < fixture.data_A_U.size(); ++i) {
    typename Fixture::MatrixType result_A_U;
    typename Fixture::VectorType result_A_orders;

    capd::matrixAlgorithms::quotientBaseMatrix(fixture.data_A_W[i], fixture.data_A_V[i], result_A_U, result_A_orders);

    BOOST_CHECK_EQUAL(fixture.data_A_U[i], result_A_U);
    BOOST_CHECK_EQUAL(fixture.data_A_orders[i], result_A_orders);
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(PARIQuotientBaseMatrix, Fixture, quotientBaseMatrix::QuotientBaseMatrixFixtures)
{
  if (! capd::matrixAlgorithms::PARIInterface::enabled()) {
    return;
  }

  Fixture fixture;
  typedef capd::matrixAlgorithms::PARIIntMatrixAlgorithms IntMatrixAlgorithms;

  for (size_t i = 0; i < fixture.data_A_U.size(); ++i) {
    typename Fixture::MatrixType result_A_U;
    typename Fixture::VectorType result_A_orders;

    typedef typename capd::matrixAlgorithms::DefaultIntMatrixAlgorithms IntMatrixAlgorithms;
    typedef typename IntMatrixAlgorithms::template QuotientBaseMatrix<typename Fixture::MatrixType>::type QuotientBaseMatrix;
    IntMatrixAlgorithms intMatrixAlgorithms;
    QuotientBaseMatrix quotientBaseMatrix = intMatrixAlgorithms.quotientBaseMatrix(fixture.data_A_W[i], fixture.data_A_V[i]);
    quotientBaseMatrix();
    result_A_U = quotientBaseMatrix.pseudoBasis();;
    quotientBaseMatrix.getOrders(result_A_orders);

    BOOST_CHECK_EQUAL(fixture.data_A_U[i], result_A_U);
    BOOST_CHECK_EQUAL(fixture.data_A_orders[i], result_A_orders);
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(cpadInvert, Fixture, SmithFixtures)
{
  typedef typename Fixture::MatrixType MatrixType;
  Fixture fixture;

  for (size_t caseIdx = 0; caseIdx < fixture.matrices_matrix.size(); ++caseIdx) {
    const MatrixType& matrixA = fixture.matrices_matrix[caseIdx];
    MatrixType matrixB = matrixA;
    const int m=matrixA.numberOfRows();
    const int n=matrixA.numberOfColumns();
    MatrixType matrixU(m, m), matrixV(n, n), matrixUInv(m, m), matrixVInv(n, n);
    int s, t;

    capd::matrixAlgorithms::smithForm(matrixB, matrixU, matrixUInv,
				      matrixV, matrixVInv, s, t);


    MatrixType tmpU(m, m);
    capd::matrixAlgorithms::invert(matrixU, tmpU);
    BOOST_CHECK_EQUAL(matrixUInv, tmpU);

    capd::matrixAlgorithms::invert(matrixUInv, tmpU);
    BOOST_CHECK_EQUAL(matrixU, tmpU);

    MatrixType tmpV(m, m);
    capd::matrixAlgorithms::invert(matrixV, tmpV);
    BOOST_CHECK_EQUAL(matrixVInv, tmpV);

    capd::matrixAlgorithms::invert(matrixVInv, tmpV);
    BOOST_CHECK_EQUAL(matrixV, tmpV);

  }
}


BOOST_AUTO_TEST_SUITE_END()
