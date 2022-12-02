/////////////////////////////////////////////////////////////////////////////
/// @file PARIConvert.h
///
/// @author Mateusz Juda <mateusz.juda@gmail.com>
///
/// @date 2014-06-17
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.edu.pl/ for details.

#ifndef CAPD_FILE_PARICONVERT_H
#define CAPD_FILE_PARICONVERT_H

#include <capd/multiPrec/MpInt.h>

namespace
{
  template<typename T>
  struct PARIConvert
  {
    static GEN to();//const T& v);
    static T from();//const GEN& v);
  };

#define INSTANCE_CONVERT_POD(type) \
  template<> \
  struct PARIConvert<type> \
  {\
    static GEN to(const type& v)\
    {\
      return stoi(v);\
    }\
    static type from(const GEN& v)\
    {\
      if (is_bigint(v)) {\
	throw std::overflow_error("PARI requires big int");\
      }\
      return itos(v);\
    }\
  };

    INSTANCE_CONVERT_POD(short);
    INSTANCE_CONVERT_POD(int);
    INSTANCE_CONVERT_POD(long);
    typedef long long llong;
    INSTANCE_CONVERT_POD(llong);

  template<typename Scalar>
  struct PARIConvert<capd::vectalg::Matrix<Scalar, 0, 0> >
  {
    static GEN to(const capd::vectalg::Matrix<Scalar, 0, 0>& capdMatrix)
    {
      const size_t rows = capdMatrix.numberOfRows();
      const size_t cols = capdMatrix.numberOfColumns();
      GEN pariMatrix = zeromatcopy(rows, cols);

      for (size_t j = 0; j < cols; ++j) { // a pari matrix is column oriented
	for (size_t i = 0; i < rows; ++i) {
	  gcoeff(pariMatrix, i + 1, j + 1) = PARIConvert<Scalar>::to(capdMatrix[i][j]);
	}
      }

      return pariMatrix;
    }

    static capd::vectalg::Matrix<Scalar, 0, 0> from(const GEN& pariMatrix)
    {
      const size_t cols = lg(pariMatrix) - 1;
      const size_t rows = (cols > 0 ? lg(gel(pariMatrix, 1)) - 1 : 0);
      capd::vectalg::Matrix<Scalar, 0, 0> capdMatrix(rows, cols);

      for (size_t i = 0; i < rows; ++i) {
	for (size_t j = 0; j < cols; ++j) {
	  GEN coeff = gcoeff(pariMatrix, i + 1, j + 1);
	  capdMatrix[i][j] = PARIConvert<Scalar>::from(coeff);
	}
      }

      return capdMatrix;
    }
  };


  template<typename Scalar>
  struct PARIConvert<capd::vectalg::Vector<Scalar, 0> >
  {
    static GEN to(const capd::vectalg::Vector<Scalar, 0>& capdVector)
    {
      const size_t rows = capdVector.dimension();
      GEN pariVector = zerovec(rows);

      for (size_t i = 0; i < rows; ++i) {
	gel(pariVector, i + 1) = PARIConvert<Scalar>::to(capdVector[i]);
      }

      return pariVector;
    }

    static capd::vectalg::Vector<Scalar, 0> from(const GEN& pariVector)
    {
      const size_t rows = lg(pariVector) - 1;
      capd::vectalg::Vector<Scalar, 0> capdVector(rows);

      for (size_t i = 0; i < rows; ++i) {
	capdVector[i] = PARIConvert<Scalar>::from(gel(pariVector, i + 1));
      }

      return capdVector;
    }
  };
}

#endif // CAPD_FILE_PARICONVERT_H
