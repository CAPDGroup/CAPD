
/////////////////////////////////////////////////////////////////////////////
/// @file DoubletonSet.h
///
/// @author kapela
/// Created on: Oct 21, 2009
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2009 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.


#ifndef _CAPD_GEOMSET_DOUBLETON_H_
#define _CAPD_GEOMSET_DOUBLETON_H_

#include "capd/geomset/AffineSet.h"

namespace capd{
namespace geomset{
/// @addtogroup geomset
/// @{


/**
 * Doubleton representation of the form x0 + C*r0 + B*r .
 *
 * We define sets representation of the form
 *
 * x + C*r0 + B * r
 *
 * where
 *  - the vector x is a center,
 *  - the matrices C and B are coordinate systems,
 *  - the vector r0 is a product of intervals and represents the set in a given coordinate system,
 *  - the vector r is used mainly to store computational errors
 *
 *  \ingroup geomset
 */
template <typename MatrixT>
class DoubletonSet : public AffineSet<MatrixT>{

public:
  typedef AffineSet<MatrixT> BaseSet;
  typedef typename BaseSet::MatrixType MatrixType;
  typedef typename BaseSet::VectorType VectorType;
  typedef typename BaseSet::ColumnVectorType ColumnVectorType;
  typedef typename BaseSet::ScalarType ScalarType;
  typedef typename MatrixType::size_type size_type;

  explicit DoubletonSet(size_type);
  explicit DoubletonSet(const VectorType& x);
  DoubletonSet(const VectorType& r0, bool);
  DoubletonSet(const VectorType& x, const VectorType& r0);
  DoubletonSet(const VectorType& x, const MatrixType& C, const VectorType& r0);
  DoubletonSet(const VectorType& x,
               const MatrixType& C, const VectorType& r0,
               const VectorType& r
  );
  DoubletonSet(const VectorType& x,
               const MatrixType& C, const VectorType& r0,
               const MatrixType& B, const VectorType& r
  );

  operator VectorType() const;
  DoubletonSet& operator=(const VectorType & v);

  /// This method computes value of an affine map f at interval vector represented by this set.
  /// This set is represented as tripleton X=x+C*r0+intersection(B*r,Q*q). Then f(X) = M*(X-x0) can be computed as
  /// M*(x-x0) + (M*C)*r0 + intersection( M*B)*r, (M*Q)*q )
  virtual VectorType affineTransformation(const MatrixType& M, const VectorType& x0) const;

  /// This method computes value of an affine functional f at the vector represented by this set.
  /// This set is represented as doubleton X=x+C*r0+B*r. Then f(X) = grad*(X-x0) can be computed as
  /// grad*(x-x0) + (grad*C)*r0 + (grad*B)*r
  virtual ScalarType evalAffineFunctional(const VectorType& gradient, const VectorType& x0) const;

  std::string toString() const;
  virtual std::string name() const { return "DoubletonSet"; }

  using BaseSet::m_x;
  using BaseSet::m_B;
  using BaseSet::m_invB;
  using BaseSet::m_r;

  /// r0 is a interval set in given coordinate system that represents 'Lipschitz' part
protected:
  VectorType m_r0;
public:
  const VectorType & get_r0 () const {return m_r0;}                            ///< returns r0
  const ScalarType & getElement_r0 (size_type i) const {return m_r0[i];}             ///< returns i-th coordinate of r0
  inline void set_r0 (const VectorType & r) {m_r0 = r;}                        ///< sets r0
  inline void setElement_r0 (size_type i, const ScalarType & s) {m_r0[i] = s;}       ///< sets i-th coordinate of r0

  /// C is a coordinate system of the 'Lipschitz' part
protected:
  MatrixType m_C;
public:
  const MatrixType & get_C () const {return m_C;}                              ///< returns matrix C
  const ScalarType & getElement_C (size_type i, size_type j) const {return m_C(i+1,j+1);}  ///< returns element C[i][j]
  VectorType getRow_C (size_type i) const {return m_C.row(i);}                       ///< returns i-th row of C
  ColumnVectorType getColumn_C (size_type j) const {return m_C.column(j);}           ///< returns j-th column of C
  void set_C (const MatrixType & B) {m_C = B;}                                 ///< sets matrix C:= B
  void setElement_C (size_type i, size_type j, const ScalarType & s) {m_C(i+1, j+1) = s;}  ///< sets C[i][j] := s
  template <typename VectorT>
  void setRow_C (size_type i, const VectorT & v) {m_C.row(i) = v;}                   ///< sets i-th row of C to v
  template <typename VectorT>
  void setColumn_C (size_type j, const VectorT & v) {m_C.column(j) = v;}             ///< sets j-th columne of C to v

};

template<typename MatrixT>
inline DoubletonSet<MatrixT>::operator typename DoubletonSet<MatrixT>::VectorType(void) const {
  return m_x + m_C*m_r0 + m_B*m_r;
}

/// @}
}} // namespace capd::geomset
#endif /* DOUBLETON_H_ */


