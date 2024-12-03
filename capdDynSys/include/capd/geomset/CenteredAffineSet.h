
/////////////////////////////////////////////////////////////////////////////
/// @file CenteredAffineSet.h
///
/// @author kapela
/// Created on: Oct 21, 2009
// ///////////////////////////////////////////////////////////////////////////

// Copyright (C) 2009 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_GEOMSET_CENTEREDAFFINESET_H_
#define _CAPD_GEOMSET_CENTEREDAFFINESET_H_

#include "capd/geomset/AffineSet.h"
namespace capd {
namespace geomset {
/// @addtogroup geomset
/// @{

/**
 * Affine set representation of the form x + B * r which assures that r contains zero.
 *
 * CenteredAffineSet represents set of the form
 *
 *   x + B * r
 *
 * where
 * - the vector x is a center,
 * - the matrix B is a coordinate system
 * - the vector r is a product of intervals and represents the set in a given coordinate system.
 * Constructors assures that r contains zero.
 */
template<typename MatrixT>
class CenteredAffineSet : public capd::geomset::AffineSet<MatrixT> {

public:
  typedef capd::geomset::AffineSet<MatrixT> BaseSet;
  typedef MatrixT MatrixType;
  typedef  typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ColumnVectorType ColumnVectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename MatrixType::size_type size_type;

  /// x:=0 r:=0 B:=Id
  explicit CenteredAffineSet(size_type dim) : BaseSet(dim) {
  }
  /// x:=mid(v)  r:=[-radius(v),radius(v)]  B:=Id
  CenteredAffineSet(const VectorType& v) : BaseSet(v){
     // splitting already done in AffineSet
  }
  /// We do not split x
  /// x:= x, r:=0, B:=Id
  CenteredAffineSet(const VectorType& x, bool) : BaseSet(x, true){
  }
  CenteredAffineSet(const VectorType& x, const VectorType& r) : BaseSet(x,r){
    if (!subset(VectorType(x.dimension()), m_r)) { // if m_r does not contain zero
      m_x += m_r;                                  //   we center the set
      split(m_x, m_r);
    }
  }

  CenteredAffineSet(const VectorType& x,
                    const MatrixType& B, const VectorType& r)
    :  BaseSet(x,B,r) {
    if (!subset(VectorType(x.dimension()), m_r)) {   // if m_r does not contain zero
      VectorType centerR = m_r;                      //   we center the set
      split(centerR, m_r);
      // m_x will be not a single point, but we must assure m_r contains zero
      m_x += m_B * centerR;
    }
  }

  /// This method computes value of functor f at interval vector represented by this set.
  /// This set is represented as X=x+B*r, where r contains zero. Then f(X) can be evaluated as
  /// f(x) + (Df(X)*B)*r
  template<class Functional>
  ScalarType evalAt(const Functional& f, const VectorType& gradient) const{
    ScalarType r = f(this->m_x);
    for(size_type i=0;i<this->m_x.dimension();++i)
      r += (gradient*this->m_B.column(i))*this->m_r[i];
    return r;
  }

  virtual std::string name() const { return "CenteredAffineSet"; }

  using BaseSet::m_x;
  using BaseSet::m_B;
  using BaseSet::m_invB;
  using BaseSet::m_r;

};

/// @}
}} //capd::geomset

#endif // _CAPD_GEOMSET_CENTEREDAFFINESET_H_

