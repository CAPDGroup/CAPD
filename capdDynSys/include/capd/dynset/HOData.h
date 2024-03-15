/////////////////////////////////////////////////////////////////////////////
///
/// @file HOData.h
///
/// @author Daniel Wilczak
/// Created on: Apr 15, 2014
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSET_HODATA_H_
#define _CAPD_DYNSET_HODATA_H_

#include "capd/basicalg/factrial.h"
#include "capd/matrixAlgorithms/floatMatrixAlgorithms.h"
#include "capd/dynset/DoubletonData.h"

namespace capd{
namespace dynset{
/// @addtogroup dynset 
/// @{
// @addtogroup capd
/// @{

/**
 * This class is a data structure used in implementation of all types of HO sets.
 * It stores temporary objects that do not need to be allocated in each call to move function.
 * It also helps in implementation of C^1 algorithms.
 */

template <class BaseData>
struct HOData : public BaseData {

  typedef typename BaseData::MatrixType MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename MatrixType::size_type size_type;

  explicit HOData(size_type dimension)
    : BaseData(dimension),
      psiPlus(dimension), psiMinus(dimension), pBound(dimension), cBound(dimension),
      JPlus(dimension,dimension), JMinus(dimension,dimension),
      midJMinusInverse(dimension,dimension), T(dimension,dimension)
  {}

  void computeC0HORemainder(size_type p, size_type q){
    double sign = q%2 ? -1.0 : 1.0;
    if(p+q>60){
      throw std::runtime_error("Hermite-Obreshkov remainder: the order of the method exceeds limit 60! ");
    }

    const ScalarType r = ScalarType(sign)/binomial(p+q,q);
    const size_type dimension = this->x.dimension();
    for(size_type i=0;i<dimension;++i)
      this->rem[i] = this->psiPlus[i] - this->psiMinus[i] + r*this->rem[i];
  }

  void computeC1HORemainder(size_type p, size_type q, MatrixType& jacRem){
    double sign = q%2 ? -1.0 : 1.0;
    if(p+q>60){
      throw std::runtime_error("Hermite-Obreshkov remainder: the order of the method exceeds limit 60! ");
    }
    jacRem *= ScalarType(sign)/binomial(p+q,q);
  }

  void computeC0HOCoefficients(){
    mid(this->JMinus,this->midJMinusInverse);
    mid(capd::matrixAlgorithms::gaussInverseMatrix(this->midJMinusInverse),this->midJMinusInverse);
    matrixByMatrix(this->midJMinusInverse,this->JPlus,this->jacPhi);
    matrixByMatrix(this->midJMinusInverse,this->JMinus,this->T);
    for(size_type i=1;i<=this->T.numberOfRows();++i)
      this->T(i,i) -= TypeTraits<typename ScalarType::BoundType>::one();
  }

  VectorType psiPlus, psiMinus, pBound, cBound;
  MatrixType JPlus, JMinus, midJMinusInverse, T;
};

/// @}
}} // namespace capd::dynset

#endif // _CAPD_DYNSET_HOBASESET_H_

