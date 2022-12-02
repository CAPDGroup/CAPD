/////////////////////////////////////////////////////////////////////////////
///
/// @file DoubletonData.h
///
/// @author Daniel Wilczak
/// Created on: Apr 15, 2014
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSET_DOUBLETONDATA_H_
#define _CAPD_DYNSET_DOUBLETONDATA_H_

#include "capd/basicalg/factrial.h"
#include "capd/matrixAlgorithms/floatMatrixAlgorithms.h"

namespace capd{
namespace dynset{
/// @addtogroup dynset 
/// @{
// @addtogroup capd
/// @{

/**
 * This class is a data structure used in implementation of all types of Doubleton and Tripleton sets.
 * It stores temporary objects that do not need to be allocated in each call to move function.
 */

template <class MatrixT>
struct DoubletonData {

  typedef MatrixT MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename MatrixType::size_type size_type;

  explicit DoubletonData(size_type dimension)
    : x(dimension), deltaX(dimension), y(dimension), deltaY(dimension), rem(dimension), enc(dimension),
      jacPhi(dimension, dimension), deltaC(dimension,dimension), B(dimension,dimension)
  {}

  VectorType x, deltaX, y, deltaY, rem, enc;
  MatrixType jacPhi, deltaC, B;
};

template<class MatrixT>
struct TripletonData : public DoubletonData<MatrixT>{
  typedef MatrixT MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename MatrixType::size_type size_type;

  explicit TripletonData(size_type dimension)
    : DoubletonData<MatrixType>(dimension),
      qr(dimension), Q(dimension,dimension)
  {}

  VectorType qr;
  MatrixType Q;
};

template<class MatrixT>
struct C1DoubletonData : public DoubletonData<MatrixT>{
  typedef MatrixT MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename MatrixType::size_type size_type;

  explicit C1DoubletonData(size_type dimension)
    : DoubletonData<MatrixType>(dimension),
      jacRem(dimension,dimension), jacEnc(dimension,dimension)
  {}

  MatrixType jacRem, jacEnc;
};

/// @}
}} // namespace capd::dynset

#endif // _CAPD_DYNSET_DOUBLETONSET_H_

