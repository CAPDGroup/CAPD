/////////////////////////////////////////////////////////////////////////////
/// @file PARIInterfaceTest.cpp
///
/// @author Mateusz Juda <mateusz.juda@gmail.com>
///
/// @date 2014-06-13
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD::RedHom Group.
//
// This file constitutes a part of the CAPD::RedHom library,
// distributed under the terms of the GNU General Public License.
// Consult  http://redhom.ii.edu.pl/ for details.


#include <capd/matrixAlgorithms/PARISmithForm.h>
#include <capd/auxil/Logger.h>
#include <capd/matrixAlgorithms/intMatrixAlgorithms.hpp>
#include <capd/matrixAlgorithms/PARIInterface.h>
#include <capd/multiPrec/MpInt.h>

#include "SmithFixtures.h"

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
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
#include <stdexcept>

namespace bl = boost::lambda;
namespace b = boost;

using namespace boost::assign;
using namespace boost::lambda;

using namespace capd;
using namespace capd::matrixAlgorithms;
using namespace capd::test;


BOOST_AUTO_TEST_SUITE(PARIInterfaceTestSuite)

typedef capd::vectalg::Vector<long long, 0> ZVector;
typedef capd::vectalg::Matrix<long long, 0, 0> ZMatrix;

typedef capd::vectalg::Vector<capd::multiPrec::MpInt, 0> MpVector;
typedef capd::vectalg::Matrix<capd::multiPrec::MpInt, 0, 0> MpMatrix;


BOOST_AUTO_TEST_CASE_TEMPLATE(smithFormWithoutMatrices, Fixture, SmithFixtures)
{
  Fixture fixture;

  PARIInterface& pariInterface = PARIInterface::instance();

  for (size_t caseIdx = 0; caseIdx < fixture.matrices_matrix.size(); ++caseIdx) {
    const ZMatrix& matrixCAPD = fixture.matrices_matrix[caseIdx];
    const ZVector divisors = pariInterface.smithForm(matrixCAPD);
    ZMatrix smithMatrix(matrixCAPD.numberOfRows(), matrixCAPD.numberOfColumns());

    BOOST_REQUIRE(divisors.dimension() <= smithMatrix.numberOfRows());
    BOOST_REQUIRE(divisors.dimension() <= smithMatrix.numberOfColumns());

    for (size_t i = 0; i < divisors.dimension(); ++i) {
      smithMatrix[i][i] = divisors[divisors.dimension() - i - 1];
    }

    BOOST_CHECK_EQUAL(fixture.matrices_smithMatrix[caseIdx], smithMatrix);
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(smithFormWithTwoMatrices, Fixture, SmithFixtures)
{
  Fixture fixture;
  PARIInterface& pariInterface = PARIInterface::instance();

  for (size_t caseIdx = 0; caseIdx < fixture.matrices_matrix.size(); ++caseIdx) {
    const ZMatrix& matrixCAPD = fixture.matrices_matrix[caseIdx];
    const int m = matrixCAPD.numberOfRows();
    const int n = matrixCAPD.numberOfColumns();
    ZMatrix matrixU(m, m), matrixV(n, n);

    try {
      const ZVector divisors = pariInterface.smithForm(matrixCAPD, matrixU, matrixV);

      for (size_t i = 0, end = matrixV.numberOfColumns();
	   i < end / 2; ++i) {
	capd::matrixAlgorithms::columnExchange(matrixV, end - i, i + 1);
      }

      for (size_t i = 0, end = matrixU.numberOfRows();
	   i < end / 2; ++i) {
	capd::matrixAlgorithms::rowExchange(matrixU, end - i, i + 1);
      }

      const ZMatrix diagonal = matrixU * matrixCAPD * matrixV;
      BOOST_CHECK_EQUAL(fixture.matrices_smithMatrix[caseIdx], diagonal);
    } catch (std::overflow_error& ex) {
      if (!fixture.bigInt)
	throw ex;
    }

  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(pariSmithFormMatrices, Fixture, SmithFixtures)
{
  typedef typename Fixture::MatrixType MatrixType;
  Fixture fixture;
  PARIInterface& pariInterface = PARIInterface::instance();

  for (size_t caseIdx = 0; caseIdx < fixture.matrices_matrix.size(); ++caseIdx) {
    const MatrixType& matrixA = fixture.matrices_matrix[caseIdx];
    const int m=matrixA.numberOfRows();
    const int n=matrixA.numberOfColumns();
    MatrixType matrixU(m, m), matrixV(n, n), matrixUInv(m, m), matrixVInv(n, n);

    try {
      pariInterface.smithForm(matrixA, matrixU, matrixUInv,
			      matrixV, matrixVInv);

      BOOST_CHECK_EQUAL(matrixU.Identity(matrixU.numberOfRows()), matrixU * matrixUInv);
      BOOST_CHECK_EQUAL(matrixU.Identity(matrixU.numberOfRows()), matrixUInv * matrixU);
      BOOST_CHECK_EQUAL(matrixV.Identity(matrixV.numberOfRows()), matrixV * matrixVInv);
      BOOST_CHECK_EQUAL(matrixV.Identity(matrixV.numberOfRows()), matrixVInv * matrixV);
    } catch (std::overflow_error& ex) {
      if (!fixture.bigInt)
	throw ex;
    }
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(PARISmithFormTest, Fixture, SmithFixtures)
{
  typedef typename Fixture::MatrixType MatrixType;
  Fixture fixture;

  for (size_t caseIdx = 0; caseIdx < fixture.matrices_matrix.size(); ++caseIdx) {
    const MatrixType& matrixA = fixture.matrices_matrix[caseIdx];
    MatrixType matrixB = matrixA;
    const int m=matrixA.numberOfRows();
    const int n=matrixA.numberOfColumns();

    capd::matrixAlgorithms::PARISmithForm<typename MatrixType::ScalarType> smithForm(matrixB, true, true, true, true);

    BOOST_CHECK_EQUAL(MatrixType::Identity(m), smithForm.getQ() * smithForm.getQinv());
    BOOST_CHECK_EQUAL(MatrixType::Identity(m), smithForm.getQinv() * smithForm.getQ());
    BOOST_CHECK_EQUAL(MatrixType::Identity(n), smithForm.getR() * smithForm.getRinv());
    BOOST_CHECK_EQUAL(MatrixType::Identity(n), smithForm.getRinv() * smithForm.getR());
  }
}

#ifdef HAVE_MPCAPD_ALG

MpMatrix promote(const ZMatrix& orgData)
{
  const size_t rows = orgData.numberOfRows();
  const size_t cols = orgData.numberOfColumns();

  MpMatrix promoted(rows, cols);

  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < cols; ++j) {
      ZMatrix::ScalarType scalar = orgData[i][j];
      promoted[i][j] = (long) scalar; // for now it is enough, there is no support for long long.
      //      mpz_import(promoted[i][j].get_mpz_t(), 1, -1, sizeof(scalar), 0, 0, &scalar);
    }
  }

  return promoted;
}

BOOST_AUTO_TEST_CASE_TEMPLATE(smithFormWithMatricesMp, Fixture, SmithFixtures)
{
  Fixture fixture;
  PARIInterface& pariInterface = PARIInterface::instance();

  for (size_t caseIdx = 0; caseIdx < fixture.matrices_matrix.size(); ++caseIdx) {
    const ZMatrix& orgData = fixture.matrices_matrix[caseIdx];
    const MpMatrix matrixCAPD = promote(orgData);
    const int m=matrixCAPD.numberOfRows();
    const int n=matrixCAPD.numberOfColumns();
    MpMatrix matrixU(m, m), matrixV(n, n);

    const MpVector divisors = pariInterface.smithForm(matrixCAPD, matrixU, matrixV);

    for (size_t i = 0, end = matrixV.numberOfColumns();
	 i < end / 2; ++i) {
      capd::matrixAlgorithms::columnExchange(matrixV, end - i, i + 1);
    }

    for (size_t i = 0, end = matrixU.numberOfRows();
	 i < end / 2; ++i) {
      capd::matrixAlgorithms::rowExchange(matrixU, end - i, i + 1);
    }

    const MpMatrix diagonal = matrixU * matrixCAPD * matrixV;
    BOOST_CHECK_EQUAL(promote(fixture.matrices_smithMatrix[caseIdx]), diagonal);
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(PARISmithFormMpTest, Fixture, SmithFixtures)
{
  typedef typename Fixture::MatrixType MatrixType;
  Fixture fixture;

  for (size_t caseIdx = 0; caseIdx < fixture.matrices_matrix.size(); ++caseIdx) {
    const MpMatrix matrixA = promote(fixture.matrices_matrix[caseIdx]);
    MpMatrix matrixB = matrixA;
    const int m=matrixA.numberOfRows();
    const int n=matrixA.numberOfColumns();

    capd::matrixAlgorithms::PARISmithForm<typename MpMatrix::ScalarType> smithForm(matrixB, true, true, true, true);

    BOOST_CHECK_EQUAL(MpMatrix::Identity(m), smithForm.getQ() * smithForm.getQinv());
    BOOST_CHECK_EQUAL(MpMatrix::Identity(m), smithForm.getQinv() * smithForm.getQ());
    BOOST_CHECK_EQUAL(MpMatrix::Identity(n), smithForm.getR() * smithForm.getRinv());
    BOOST_CHECK_EQUAL(MpMatrix::Identity(n), smithForm.getRinv() * smithForm.getR());
  }
}

#endif // HAVE_MPCAPD_ALG

BOOST_AUTO_TEST_SUITE_END()
