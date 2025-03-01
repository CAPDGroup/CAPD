/////////////////////////////////////////////////////////////////////////////
/// @file QRPolicy.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef CAPD_DYNSET_QRPOLICY_H
#define CAPD_DYNSET_QRPOLICY_H

#include "capd/vectalg/iobject.hpp"
#include "capd/matrixAlgorithms/floatMatrixAlgorithms.hpp"
#include "capd/dynset/DefaultPolicy.h"
#include "capd/diffAlgebra/Hessian.h"
#include "capd/vectalg/Matrix_Interval.hpp"

namespace capd{
namespace dynset{
/// @addtogroup dynset
/// @{

template<class BasePolicy = DefaultPolicy>
class InverseQRPolicy : public BasePolicy
{
public:
  template<class VectorT, class MatrixT>
  void computeBinvB(MatrixT& B, MatrixT& invB, const VectorT& v) const
  {
    try{
      mid(B,B);
      invB = capd::matrixAlgorithms::krawczykInverse(B);
    }catch(...)
    {
      BasePolicy::computeBinvB(B,invB,v);
    }
  }
  virtual std::string toString() const {
    return "InverseQRPolicy " + BasePolicy::toString();
  }
};

template<class BasePolicy = DefaultPolicy>
class FullQRWithPivoting : public BasePolicy
{
public:
  template<class VectorT, class MatrixT>
  void computeBinvB(MatrixT& B, MatrixT& invB, const VectorT& v) const
  {
    try{
      mid(B,B);
      capd::matrixAlgorithms::orthonormalize(B,v);
      invB = Transpose(B);
    }catch(...)
    {
      BasePolicy::computeBinvB(B,invB,v);
    }
  }

  template<class MatrixT>
  void computeBinvB(MatrixT& B, MatrixT& invB,
                    const capd::diffAlgebra::Hessian<typename MatrixT::ScalarType,MatrixT::RowVectorType::csDim,MatrixT::RowVectorType::csDim>& H
                    ) const
  {
    typename MatrixT::RowVectorType d(B.numberOfRows());
    typedef typename MatrixT::size_type size_type;
    try{
      mid(B,B);
      for(size_type i=0;i<B.numberOfRows();++i)
      {
        for(size_type j=0;j<B.numberOfRows();++j)
          for(size_type c=j;c<B.numberOfRows();++c)
            d[i] = intervalHull(H(i,j,c),d[i]);
      }
      capd::matrixAlgorithms::orthonormalize(B,d);
      invB = Transpose(B);
    }catch(...)
    {
      BasePolicy::computeBinvB(B,invB,d);
    }
  }
  virtual std::string toString() const {
    return "FullQRWithPivoting " + BasePolicy::toString();
  }
};

template<int N, class BasePolicy = DefaultPolicy>
class PartialQRWithPivoting : public BasePolicy
{
public:
  template<class MatrixT>
  void reduceMatrix(MatrixT& B) const
  {
    mid(B,B);
    typedef typename MatrixT::size_type size_type;
    for(size_type i=N+1;i<=B.numberOfColumns();++i)
      for(size_type j=1;j<i;++j)
        B(j,i)= typename MatrixT::ScalarType(0.);
  }

  template<class VectorT, class MatrixT>
  void orthonormalize(MatrixT& B, const VectorT& v)
  {
    reduceMatrix(B);
    capd::matrixAlgorithms::orthonormalize(B,v);
  }

  template<class MatrixT>
  void orthonormalize(MatrixT& B)
  {
    reduceMatrix(B);
    capd::matrixAlgorithms::orthonormalize(B);
  }

  template<class VectorT, class MatrixT>
  void computeBinvB(MatrixT& B, MatrixT& invB, const VectorT& v) const
  {
    try{
      reduceMatrix(B);
      capd::matrixAlgorithms::orthonormalize(B,v);
      invB = Transpose(B);
      invB = capd::matrixAlgorithms::krawczykInverse(B);
      // TODO: Maybe it should be replaced by
      // capd::matrixAlgorithms::krawczykCorrection(A,invA);
    }catch(...)
    {
      BasePolicy::computeBinvB(B,invB,v);
    }
  }
  template<class MatrixT>
  void computeBinvB(MatrixT& B, MatrixT& invB,
                    const capd::diffAlgebra::Hessian<typename MatrixT::ScalarType,MatrixT::RowVectorType::csDim,MatrixT::RowVectorType::csDim>& H
                    ) const
  {
    typename MatrixT::RowVectorType d(B.numberOfRows());
    typedef typename MatrixT::size_type size_type;
    try{
      reduceMatrix(B);
      for(size_type i=0;i<B.numberOfRows();++i)
      {
        for(size_type j=0;j<B.numberOfRows();++j)
          for(size_type c=j;c<B.numberOfRows();++c)
            d[i] = intervalHull(H(i,j,c),d[i]);
      }
      capd::matrixAlgorithms::orthonormalize(B,d);
      invB = Transpose(B);
    }catch(...)
    {
      BasePolicy::computeBinvB(B,invB,d);
    }
  }
  virtual std::string toString() const {
    return "PartialQRWithPivoting " + BasePolicy::toString();
  }
};

/**
 *  Vectors are orthogonalized only if they are close to be parallel
 */
template<class BasePolicy = DefaultPolicy>
class SelectiveQRWithPivoting : public BasePolicy
{
    double level;

public:
    SelectiveQRWithPivoting(double parallelTolerance = 1.0e-15)
            : level(parallelTolerance){}
  ///  if vectors are closer to be parallel more than given tolerance level
  ///  one of them will be orthogonalized
  void setParallelTolerance(const double & level) { this->level = level; }
  double getParallelTolerance() const { return this->level; }

  template<class VectorType, class MatrixType>
  void orthonormalize(MatrixType & B, const VectorType & v) const
  {
    typedef typename VectorType::template rebind<int>::other IntVectorType;
    typedef typename MatrixType::size_type size_type;

    capd::vectalg::EuclNorm<VectorType, MatrixType> norm;

    mid(B,B);

    const size_type dim = v.dimension();
    VectorType sizes(dim);
    IntVectorType p(dim);
    MatrixType Q(dim, dim), R(dim, dim);
    capd::matrixAlgorithms::QRdecomposeWithPivoting(B, v, Q, R, sizes, p);

    for(size_type i=0; i<dim; i++){
      if(sizes[i].rightBound()/norm(B.column(i))<= level){
        B.column(i) = Q.column(i);
      }
      else{
        B.column(i).normalize();
      }
    }
  }

  template<class MatrixT>
  void orthonormalize(MatrixT& /*B*/)
  {
           throw std::runtime_error("SelectiveQRWithPivoting::orthonormalize - NOT IMPLEMENTED METHOD!");
  }

  template<class VectorT, class MatrixT>
  void computeBinvB(MatrixT& B, MatrixT& invB, const VectorT& v) const
  {
    mid(B,B);
    try{
      this->orthonormalize(B, v);
    }catch(...) {
      std::cout << "we failed to orthonormalize!";
      BasePolicy::computeBinvB(B,invB,v);
    }
    try{
      mid(B,B);
      invB = capd::matrixAlgorithms::krawczykInverse(B);
    }catch(...) {
      std::cout << "we failed to invert!\n\n";
      BasePolicy::computeBinvB(B,invB,v);
    }
  }
  virtual std::string toString() const {
    return "SelectiveQRWithPivoting " + BasePolicy::toString();
  }
};
/// @}
}} // namespace capd::dynset



#endif // CAPD_DYNSET_QRPOLICY_H




