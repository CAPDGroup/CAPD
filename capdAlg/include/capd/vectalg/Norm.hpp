/// @addtogroup vectalg
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file Norm.hpp
///
/// @author The CAPD Group
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2013 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#include <sstream>
#include "capd/vectalg/Norm.h"
#include "capd/basicalg/minmax.h"
#include "capd/matrixAlgorithms/floatMatrixAlgorithms.hpp"

#ifndef _CAPD_VECTALG_NORM_HPP_ 
#define _CAPD_VECTALG_NORM_HPP_ 

namespace capd{
namespace vectalg{

/**
 *  Computes norm of the vector 
 *
 *  \f$ |x| = \sqrt{ x_1^2 + x_2^2 + \dots + x_n^2} \f$
 */
template<typename VectorType, typename MatrixType>
typename EuclNorm<VectorType, MatrixType>::ScalarType
EuclNorm<VectorType, MatrixType>::operator()(const VectorType & x) const{
  
  ScalarType s(0.);
  for(typename VectorType::const_iterator i = x.begin(); i!=x.end();++i)
    s+=power(*i, 2);
  return sqrt(nonnegativePart(s));
}

/**
 *  Computes norm of the matrix
 *
 *  \f$ |A| = \sqrt{ \rho(A^T A) } \f$
 *  where $\rho$ is spectral radius of a matrix (the biggest modulus of an eigenvalue)
 */
template<typename VectorType, typename MatrixType>
typename EuclNorm<VectorType, MatrixType>::ScalarType
EuclNorm<VectorType, MatrixType>::operator()(const MatrixType &A) const {
  return sqrt(
          nonnegativePart(
           capd::matrixAlgorithms::spectralRadiusOfSymMatrix(
             Transpose(A) * A, 
             capd::TypeTraits<ScalarType>::epsilon()*1.0e-4
         )));
}

template<typename VectorType, typename MatrixType>
std::string EuclNorm<VectorType, MatrixType>::show(void) const {
  return "Eucl norm ";
}

template<typename VectorType, typename MatrixType>
Norm<VectorType, MatrixType> * EuclNorm<VectorType, MatrixType>::clone(void) const {
  return (Norm<VectorType, MatrixType> *) new EuclNorm<VectorType, MatrixType>();
}

// --------------------------------------------------------------------------- //

template<typename VectorType, typename MatrixType>
typename MaxNorm<VectorType, MatrixType>::ScalarType
MaxNorm<VectorType, MatrixType>::operator()(const VectorType &x) const {
  ScalarType s(0.);
  for(typename VectorType::const_iterator i=x.begin();i!=x.end();++i) {
    s = capd::max(s, capd::abs(*i));
  }
  return s;
}

template<typename VectorType, typename MatrixType>
typename MaxNorm<VectorType, MatrixType>::ScalarType
MaxNorm<VectorType, MatrixType>::operator()(const MatrixType &A) const {
  ScalarType maximum(0.);
  for(size_type i=0;i<A.numberOfRows();i++) {
    // Computes sum of absolute values of entries in given row
    ScalarType rowSum(0.);
    for(size_type j=0; j<A.numberOfColumns(); ++j){
      rowSum += capd::abs(A[i][j]);
    }
    maximum = capd::max(maximum, rowSum);  
  }
  return maximum;
}

template<typename VectorType, typename MatrixType>
std::string MaxNorm<VectorType, MatrixType>::show(void) const {
  return "Max norm ";
}

template<typename VectorType, typename MatrixType>
Norm<VectorType, MatrixType> * MaxNorm<VectorType, MatrixType>::clone(void) const {
  return (Norm<VectorType, MatrixType> *) new MaxNorm<VectorType, MatrixType>();
}

// --------------------------------------------------------------------------- //

template<typename VectorType, typename MatrixType>
typename SumNorm<VectorType, MatrixType>::ScalarType
SumNorm<VectorType, MatrixType>::operator()(const VectorType &x) const {
  ScalarType s(0);
  for(typename VectorType::const_iterator i = x.begin();i!=x.end();++i)
    s += capd::abs(*i);
  return s;
}

template<typename VectorType, typename MatrixType>
typename SumNorm<VectorType, MatrixType>::ScalarType
SumNorm<VectorType, MatrixType>::operator()(const MatrixType &A) const {
  VectorType x(A.numberOfRows());
  MaxNorm<VectorType, MatrixType> maxNorm;
  SumNorm<VectorType, MatrixType> sumNorm;
  for(size_type i=0;i<A.numberOfRows();++i)
    x[i] = sumNorm(A.column(i));
  return maxNorm(x);
}

template<typename VectorType, typename MatrixType>
std::string SumNorm<VectorType, MatrixType>::show(void) const {
  return "Sum norm ";
}

template<typename VectorType, typename MatrixType>
Norm<VectorType, MatrixType> * SumNorm<VectorType, MatrixType>::clone(void) const {
  return (Norm<VectorType, MatrixType> *) new SumNorm<VectorType, MatrixType>();
}

// --------------------------------------------------------------------------- //

template<typename VectorType, typename MatrixType>
typename EuclLNorm<VectorType, MatrixType>::ScalarType
EuclLNorm<VectorType, MatrixType>::operator()(const VectorType &x) const {
  ScalarType s(0.);
  for(typename VectorType::const_iterator i = x.begin();i!=x.end();++i)
    s += power(*i, 2);
  return sqrt(nonnegativePart(s));
}

template<typename VectorType, typename MatrixType>
typename EuclLNorm<VectorType, MatrixType>::ScalarType
EuclLNorm<VectorType, MatrixType>::operator()(const MatrixType &A) const {
  return capd::matrixAlgorithms::maxEigenValueOfSymMatrix(
          (Transpose(A) + A)/static_cast<ScalarType>(2),
          capd::TypeTraits<ScalarType>::epsilon()*1.0e-6
         );
}

template<typename VectorType, typename MatrixType>
std::string EuclLNorm<VectorType, MatrixType>::show(void) const {
 return "Eucl log norm ";
}

template<typename VectorType, typename MatrixType>
Norm<VectorType, MatrixType> * EuclLNorm<VectorType, MatrixType>::clone(void) const {
  return (Norm<VectorType, MatrixType> *) new EuclLNorm<VectorType, MatrixType>();
}

// --------------------------------------------------------------------------- //

template<typename VectorType, typename MatrixType>
typename MaxLNorm<VectorType, MatrixType>::ScalarType
MaxLNorm<VectorType, MatrixType>::operator()(const VectorType &x) const {
  ScalarType s(0.);
  for(typename VectorType::const_iterator i=x.begin();i!=x.end();++i) {
    s = capd::max(s, capd::abs(*i));
  }
  return s;
}

template<typename VectorType, typename MatrixType>
typename MaxLNorm<VectorType, MatrixType>::ScalarType
MaxLNorm<VectorType, MatrixType>::operator()(const MatrixType & A) const {
  ScalarType res(0.0);
  for(size_type i=0; i<A.numberOfRows(); i++) {
    ScalarType s=A[i][i];
    
    for(size_type j=0; j<A.numberOfColumns(); j++)
      if(j != i)
        s += capd::abs(A[i][j]);
    
    if(i==0)
      res = s;
    else {
      res = capd::max(res, s);
    }
  }
  return(res);
}

template<typename VectorType, typename MatrixType>
std::string MaxLNorm<VectorType, MatrixType>::show(void) const {
  return "Max log norm ";
}

template<typename VectorType, typename MatrixType>
Norm<VectorType, MatrixType> * MaxLNorm<VectorType, MatrixType>::clone(void) const {
  return (Norm<VectorType, MatrixType> *) new MaxLNorm<VectorType, MatrixType>();
}

// --------------------------------------------------------------------------- //

template<typename VectorType, typename MatrixType>
typename SumLNorm<VectorType, MatrixType>::ScalarType
SumLNorm<VectorType, MatrixType>::operator()(const VectorType &x) const {
  ScalarType s(0.0);
  for(typename VectorType::const_iterator i=x.begin();i!=x.end();++i)
    s += capd::abs(*i);
  return s;
}

template<typename VectorType, typename MatrixType>
typename SumLNorm<VectorType, MatrixType>::ScalarType
SumLNorm<VectorType, MatrixType>::operator()(const MatrixType &A) const {
  ScalarType res(0.0);
  for(size_type j=0; j<A.numberOfColumns(); j++) {
    ScalarType s=A[j][j];
    
    for(size_type i=0;i<A.numberOfRows();i++)
      if(i != j)
        s += capd::abs(A[i][j]);
    
    s=right(s);
    if(j==0)
      res=s;
    else {
      res = capd::max(res, s);
    }
  }
  return(res);
}

template<typename VectorType, typename MatrixType>
std::string SumLNorm<VectorType, MatrixType>::show(void) const {
  
 return "Sum log norm ";
}

template<typename VectorType, typename MatrixType>
Norm<VectorType, MatrixType> * SumLNorm<VectorType, MatrixType>::clone(void) const {
  
  return (Norm<VectorType, MatrixType> *) new SumLNorm<VectorType, MatrixType>();
}

}} // namespace capd::vectalg

#endif // _CAPD_VECTALG_NORM_HPP_ 

/// @}
