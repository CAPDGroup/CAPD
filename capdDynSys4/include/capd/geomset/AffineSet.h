/////////////////////////////////////////////////////////////////////////////
// Package: CAPD

/////////////////////////////////////////////////////////////////////////////
/// @file AffineSet.h
///
/// @author kapela
/// @date Oct 21, 2009
// ///////////////////////////////////////////////////////////////////////////

// Copyright (C) 2009 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_GEOMSET_AFFINESET_H_
#define _CAPD_GEOMSET_AFFINESET_H_

#include "capd/matrixAlgorithms/floatMatrixAlgorithms.h"


namespace capd {
/// Definitions of sets as geometrical objects which can represent different shapes.
namespace geomset {

/// @addtogroup geomset
/// @{


/**
 * Affine set representanion of the form x0 + B*r .
 *
 * We define sets representation of the form
 *
 * x + B * r
 *
 * where
 *   - the vector x is a center,
 *   - the matrix B is a coordinate system
 *   - the vector r is a product of intervals and represents the set in a given coordinate system.
 */
template<typename MatrixT>
class AffineSet {

public:
  typedef MatrixT MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ColumnVectorType ColumnVectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename MatrixType::size_type size_type;
  typedef typename MatrixType::template rebind<typename ScalarType::BoundType>::other RealMatrix;
  typedef typename VectorType::template rebind<typename ScalarType::BoundType>::other RealVector;

  explicit AffineSet(size_type);                                                     // x:=0,      B:=Id  r:=0,
  explicit AffineSet(const VectorType& v);                                     // x:=mid(v)  B:=Id  r:=[-radius(v), radius(v)]
  AffineSet(const VectorType& x, bool);                                        // x:=x,      B:=Id  r:=0
  AffineSet(const VectorType& x, const VectorType& r);                         // x:=x       B:=Id  r:=r
  AffineSet(const VectorType& x, const MatrixType& B, const VectorType& r);    // x:=x       B:=B   r:=r
  virtual ~AffineSet() {}
  /// returns dimension of a set
  size_type dimension() const { return m_x.dimension(); }
  /// returns interval vector that contains set.
  operator VectorType() const;
  /// returns set detailed information
  virtual std::string toString() const;
  /// returns set's name
  virtual std::string name() const { return "AffineSet"; }
  /// returns set image after affine transformation
  virtual VectorType affineTransformation(const MatrixType&, const VectorType&) const;

  /// This method computes value of an affine functional f at the vector represented by this set.
  /// This set is represented as doubleton X=x+B*r. Then f(X) = grad*(X-x0) can be computed as
  /// grad*(x-x0) + (grad*B)*r
  virtual ScalarType evalAffineFunctional(const VectorType& gradient, const VectorType& x0) const;

  protected:
  /// x is a center of the set
  VectorType m_x;
public:
  const VectorType & get_x () const {return m_x;}
  const ScalarType & getElement_x (int i) const {return m_x[i];}
  void set_x (const VectorType & x) {m_x = x;}
  void setElement_x (size_type i, const ScalarType & s) {m_x[i] = s;}

//ADD_VECTOR_MEMBER(x, protected, public, public);


protected:
  /// r is a interval set in given coordinate system
  VectorType m_r;
public:
  const VectorType & get_r () const {return m_r;}
  const ScalarType & getElement_r (size_type i) const {return m_r[i];}
  inline void set_r (const VectorType & r) {m_r = r;}
  inline void setElement_r (size_type i, const ScalarType & s) {m_r[i] = s;}


protected:
  /// B is a coordinate system
  MatrixType m_B;
  MatrixType m_invB;
public:
  const MatrixType & get_B () const {return m_B;}
  const MatrixType & get_invB () const {return m_invB;}
  void setToIdentity(){
    m_B.setToIdentity();
    m_invB.setToIdentity();
    m_r.clear();
  }
  const ScalarType & getElement_B (size_type i, size_type j) const {return m_B(i+1,j+1);}
  VectorType getRow_B (size_type i) const {return m_B.row(i);}
  ColumnVectorType getColumn_B (size_type j) const {return m_B.column(j);}
  void set_B (const MatrixType & B) { m_B = B;  }
  void set_invB (const MatrixType & B) { m_invB = B;  }
  void setElement_B (size_type i, size_type j, const ScalarType & s) {m_B(i+1, j+1) = s;}
  template <typename VectorT>
  void setRow_B (size_type i, const VectorT & v) {m_B.row(i) = v;}
  template <typename VectorT>
  void setColumn_B (size_type j, const VectorT & v) {m_B.column(j) = v;}

};
/// @}

template<typename MatrixT>
inline AffineSet<MatrixT>::operator typename AffineSet<MatrixT>::VectorType(void) const {
  return m_x + m_B*m_r;
}

}} //end of namespace capd::geomset
#endif

