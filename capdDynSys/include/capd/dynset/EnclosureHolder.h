/////////////////////////////////////////////////////////////////////////////
/// @file EnclosureHolder.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#include "capd/diffAlgebra/Hessian.h"
#include "capd/dynset/AbstractSet.h"

#ifndef _CAPD_DYNSET_ENCLOSUREHOLDER_H_
#define _CAPD_DYNSET_ENCLOSUREHOLDER_H_

namespace capd{
namespace dynset{
/// @addtogroup dynset
/// @{

/**
 * These classes are used as base classes for all types of C^0-C^n sets.
 * They store last computed enclosure for the flow and provide methods for manipulating them.
 * They store as well the current set and/or derivatives represented in the canonical coordinates.
 */


template<class VectorT>
class C0EnclosureHolder : public AbstractSet<VectorT>{
public:
  C0EnclosureHolder(const VectorT& x, const VectorT& enc)
    : m_currentSet(x), m_lastEnclosure(enc)
  {}
  virtual ~C0EnclosureHolder() {}

  const VectorT& getLastEnclosure() const{
    return m_lastEnclosure;
  }
  /// returns an enclosure of the set in the canonical coordinates
  virtual operator VectorT() const{
    return m_currentSet;
  }
  void setCurrentSet(const VectorT& x) {
    m_currentSet = x;
  }

protected:
  void setLastEnclosure(const VectorT& enc){
    m_lastEnclosure = enc;
  }
  VectorT m_currentSet;
  VectorT m_lastEnclosure;
};

template<class MatrixT>
class C1EnclosureHolder : public C0EnclosureHolder<typename MatrixT::RowVectorType>{
public:
  typedef typename MatrixT::RowVectorType VectorType;
  C1EnclosureHolder(const VectorType& x, const VectorType& enc, const MatrixT& M, const MatrixT& mEnc)
    : C0EnclosureHolder<VectorType>(x,enc),
      m_currentMatrix(M), m_lastMatrixEnclosure(mEnc)
  {}

  const MatrixT& getLastMatrixEnclosure() const{
    return m_lastMatrixEnclosure;
  }
  /// returns an enclosure of derivative in the canonical coordinates
  virtual operator MatrixT() const{
    return m_currentMatrix;
  }

protected:
  void setLastMatrixEnclosure(const MatrixT& M){
    m_lastMatrixEnclosure = M;
  }

  MatrixT m_currentMatrix;
  MatrixT m_lastMatrixEnclosure;
};

template<class MatrixT>
class C2EnclosureHolder : public C1EnclosureHolder<MatrixT>{
public:
  typedef typename MatrixT::ScalarType ScalarType;
  typedef typename MatrixT::RowVectorType VectorType;
  typedef capd::diffAlgebra::Hessian<ScalarType,VectorType::csDim,VectorType::csDim> HessianType;

  C2EnclosureHolder(const VectorType& x, const VectorType& enc,
                    const MatrixT& M, const MatrixT& mEnc,
                    const HessianType& h, const HessianType& hEnc
  ) : C1EnclosureHolder<MatrixT>(x,enc,M,mEnc),
      m_currentHessian(h), m_lastHessianEnclosure(hEnc)
  {}

  inline const HessianType& getLastHessianEnclosure() const{
    return m_lastHessianEnclosure;
  }
  /// returns an enclosure of second order derivative in the canonical coordinates
  inline virtual operator HessianType() const{
    return m_currentHessian;
  }

protected:
  void setLastHessianEnclosure(const HessianType& h){
    m_lastHessianEnclosure = h;
  }

  HessianType m_currentHessian;
  HessianType m_lastHessianEnclosure;
};

/// @}
}} //namespace capd::dynset

#endif // _CAPD_DYNSET_ENCLOSUREHOLDER_H_
