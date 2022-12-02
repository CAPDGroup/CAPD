/////////////////////////////////////////////////////////////////////////////
/// @file CenteredTripletonSet.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_GEOMSET_CENTEREDTRIPLETONSET_H_
#define _CAPD_GEOMSET_CENTEREDTRIPLETONSET_H_

#include "capd/geomset/CenteredDoubletonSet.h"

namespace capd{
namespace geomset{
/// @addtogroup geomset
/// @{

/**
 * Class \b CenteredTripletonSet \b represents a subset of R^n in the following form
 *
 * x + C*r0 + intersection(B*r , Q*q )
 *
 * where
 *  x is a point vector
 *  C,B,Q are point matrices, where B and Q are invertible and Q is close to orthogonal
 *  r0,r are interval vectors
 *
 *  Moreover it stores rigorous inverse matrices of B and Q
 *
 */

template<typename MatrixT>
class CenteredTripletonSet : public capd::geomset::CenteredDoubletonSet<MatrixT>
{
public:
  typedef MatrixT MatrixType;
  typedef typename MatrixT::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename MatrixType::size_type size_type;
  typedef capd::geomset::CenteredDoubletonSet<MatrixT> BaseSet;

  CenteredTripletonSet(const VectorType& x);
  CenteredTripletonSet(const VectorType& x, const VectorType& r0);
  CenteredTripletonSet(const VectorType& x, const MatrixType& C, const VectorType& r0);
  CenteredTripletonSet(const VectorType& x, const MatrixType& C, const VectorType& r0, const VectorType& r);
  CenteredTripletonSet(const VectorType& x, const MatrixType& C,
      const VectorType& r0, const MatrixType& B,
      const VectorType& r
    );

  /// This method computes value of an affine map f at interval vector represented by this set.
  /// This set is represented as tripleton X=x+C*r0+intersection(B*r,Q*q). Then f(X) = M*(X-x0) can be computed as
  /// M*(x-x0) + (M*C)*r0 + intersection( M*B)*r, (M*Q)*q )
  virtual VectorType affineTransformation(const MatrixType& M, const VectorType& x0) const;

  /// This method computes value of an affine functional f at interval vector represented by this set.
  /// This set is represented as tripleton X=x+C*r0+intersection(B*r,Q*q). Then f(X) = grad*(X-x0) can be computed as
  /// grad*(x-x0) + (grad*C)*r0 + intersection( grad*B)*r, (grad*Q)*q )
  virtual ScalarType evalAffineFunctional(const VectorType& gradient, const VectorType& x0) const;

  operator VectorType() const{
    return this->m_x + this->m_C*this->m_r0 + intersection(this->m_B*this->m_r, this->m_Q*this->m_q);
  }

  std::string name() const { return "CenteredTripletonSet"; }

  using BaseSet::get_r0;
  using BaseSet::getElement_r0;
  using BaseSet::get_x;
  using BaseSet::getElement_x;
  using BaseSet::get_B;
  using BaseSet::get_invB;
  using BaseSet::getElement_B;
  using BaseSet::getRow_B;
  using BaseSet::getColumn_B;
  using BaseSet::get_C;
  using BaseSet::getElement_C;
  using BaseSet::getRow_C;
  using BaseSet::getColumn_C;

  // representation of set
  VectorType m_q;
  MatrixType m_Q, m_invQ;

protected:
  CenteredTripletonSet(
                  const VectorType& x,
                  const MatrixType& C, const VectorType& r0,
                  const MatrixType& B, const VectorType& r,
                  const MatrixType& Q, const VectorType& q
                  );
  void copyQ();
};

/// @}
}} //namespace capd::geomset

#endif // _CAPD_GEOMSET_CENTEREDTRIPLETONSET_H_
