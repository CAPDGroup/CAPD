/////////////////////////////////////////////////////////////////////////////
/// @file CnSet.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2020 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

/* Author: Daniel Wilczak 2006-2020 */

#ifndef _CAPD_DYNSET_CNSET_H_
#define _CAPD_DYNSET_CNSET_H_

#include "capd/diffAlgebra/TimeRange.h"
#include "capd/diffAlgebra/Jet.h"
#include "capd/vectalg/Multiindex.h"
#include "capd/vectalg/iobject.h" // for maxDiam
#include "capd/dynset/SetTraits.h"
#include "capd/dynset/AbstractSet.h"

namespace capd{
namespace dynset{
using capd::vectalg::__size_type;

/// @addtogroup dynset
/// @{

/**
 * Common interface of all sets that store Cn information (set position and derivatives to order n)
 */
template<typename MatrixT, __size_type DEGREE=0>
class CnSet : public capd::diffAlgebra::TimeRange<typename MatrixT::ScalarType>,
              public capd::dynset::AbstractSet<typename MatrixT::RowVectorType>
{
public:
  typedef MatrixT MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename MatrixType::RefColumnVectorType RefVectorType;
  typedef typename MatrixType::size_type size_type;

  typedef capd::vectalg::Multipointer Multipointer;
  typedef capd::vectalg::Multiindex Multiindex;
  typedef capd::diffAlgebra::Jet<MatrixType,DEGREE> JetType;
  typedef CnSet SetType;

  CnSet(const VectorType& x, size_type degree, ScalarType t)
    : capd::diffAlgebra::TimeRange<typename MatrixT::ScalarType>(t),
      m_currentSet(x.dimension(),degree), m_lastEnclosure(x.dimension(),degree)
  {
    m_currentSet() = x;
    m_lastEnclosure() = x;
  }

  CnSet(const JetType& x, ScalarType t)
    : capd::diffAlgebra::TimeRange<typename MatrixT::ScalarType>(t),
      m_currentSet(x), m_lastEnclosure(x)
  {}

  /// returns maximal order of partial derivative stored in the jet
  size_type degree() const     { return this->m_currentSet.degree(); }
  /// returns number of variables in the jet
  size_type dimension() const  {    return this->m_currentSet.dimension(); }
  /// returns actual maximal diameter of the current set
  ScalarType maxDiam() const { return capd::vectalg::maxDiam(this->m_currentSet()); }

  /// returns actual set as a read-only jet.
  const JetType& currentSet() const { return this->m_currentSet; }
  /// returns actual set as a jet
  JetType& currentSet() { return this->m_currentSet; }

  /// returns an enclosure for trajectories for last performed step of integration. On initialization is equal to the current set.
  RefVectorType getLastEnclosure() const    { return this->m_lastEnclosure(); }
  /// returns an enclosure for first order variational equations for last performed step of integration. On initialization is equal to the current derivatives set.
  MatrixType getLastMatrixEnclosure() const { return (MatrixType)this->m_lastEnclosure; }
  /// returns jet of enclosures computed in the last performed step of integration.
  virtual const JetType& getLastJetEnclosure() const { return this->m_lastEnclosure; }

  virtual void setLastJetEnclosure(const JetType& enc) { this->m_lastEnclosure = enc; }

  /// returns actual value of function represented by jet
  operator VectorType() const { return (VectorType) this->m_currentSet; }
  /// returns actual derivative of function represented by jet
  operator MatrixType() const { return (MatrixType) this->m_currentSet; }

  /// returns a Taylor coefficients corresponding to multipointer, (1/mp!)d^{mp}f_i
  const ScalarType& operator()(size_type i,const Multipointer& mp) const { return this->m_currentSet(i,mp); }

  /// returns a Taylor coefficients corresponding to multiindex, (1/mp!)d^{mp}f_i
  const ScalarType& operator()(size_type i,const Multiindex& mp) const { return this->m_currentSet(i,mp); }

  /// returns vector of Taylor coefficients corresponding to multipointer, i.e. result[i] = (1/mp!)d^{mp}f
  RefVectorType operator()(const Multipointer& mp) const {  return this->m_currentSet(mp); }

  /// returns Taylor coefficient corresponding to partial derivative d^2f_i/dx_jdx_c
  const ScalarType& operator()(size_type i, size_type j, size_type c) const{
    return this->m_currentSet(i,j,c);
  }

  /// returns df_j/dx_c
  const ScalarType& operator()(size_type j, size_type c) const{
    return this->m_currentSet(j,c);
  }
  /// returns actual value of f_i
  const ScalarType& operator()(size_type i) const{
    return this->m_currentSet(i);
  }

  virtual ~CnSet(){}
protected:
  JetType m_currentSet;
  JetType m_lastEnclosure;
};

template<class MatrixT,__size_type DEGREE>
struct SetTraits< CnSet<MatrixT,DEGREE> >{
	const static bool isC0Set=false;
	const static bool isC1Set=false;
	const static bool isC2Set=false;
	const static bool isCnSet=true;
};

/// @}
}} // the end of the namespace capd::dynset

#endif // _CAPD_DYNSET_CNSET_H_
