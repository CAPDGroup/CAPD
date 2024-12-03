/////////////////////////////////////////////////////////////////////////////
/// @file C1DoubletonSet.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSET_C1DOUBLETONSET_H_
#define _CAPD_DYNSET_C1DOUBLETONSET_H_

#include "capd/dynset/C1Set.h"
#include "capd/geomset/CenteredDoubletonSet.h"
#include "capd/geomset/MatrixDoubletonSet.h"
#include "capd/dynset/DoubletonData.h"
#include "capd/dynset/C0DoubletonSet.h"

namespace capd{
namespace dynset{
/// @addtogroup dynset
/// @{

/////////////////////////////////////////////////////////////////////
///
///  The C1 set is represented as doubleton: x + C*r0 + B*r;
///
///  Internal representation :
///        C*r0 - basic 'Lipschitz part'
///        B*r  - remaining errors
///
///  The class represents derivatives in a doubleton form as described in the paper C^1-Lohner algorithm by Piotr Zgliczynski (FoCM 2001).
///////////////////////////////////////////////////////////////////////

template<typename MatrixT, typename Policies>
class C1DoubletonSet : public Policies, public C1Set<MatrixT>,
                       public capd::geomset::CenteredDoubletonSet<MatrixT>,
                       public capd::geomset::MatrixDoubletonSet<MatrixT>,
                       public C1DoubletonData<MatrixT>
{
public:
  typedef MatrixT MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename MatrixType::size_type size_type;
  typedef C1Set<MatrixT> SetType;
  typedef typename C1Set<MatrixT>::DynSysType DynSysType;
  typedef capd::geomset::CenteredDoubletonSet<MatrixT> C0BaseSet;
  typedef capd::geomset::MatrixDoubletonSet<MatrixT> C1BaseSet;
  typedef C1DoubletonData<MatrixT> Data;
  typedef Policies Policy;

  C1DoubletonSet(const VectorType& x, ScalarType t = TypeTraits<ScalarType>::zero());
  C1DoubletonSet(const VectorType& x, const VectorType& r0, ScalarType t = TypeTraits<ScalarType>::zero());
  C1DoubletonSet(const VectorType& x, const MatrixType& C, const VectorType& r0, ScalarType t = TypeTraits<ScalarType>::zero());
  C1DoubletonSet(const VectorType& x, const MatrixType& C, const VectorType& r0, const VectorType& r, ScalarType t = TypeTraits<ScalarType>::zero());
  C1DoubletonSet(const VectorType& x, const MatrixType& C, const VectorType& r0, const MatrixType& B, const VectorType& r, ScalarType t = TypeTraits<ScalarType>::zero());
  C1DoubletonSet(const C0BaseSet & c0part, const C1BaseSet& c1part, ScalarType t = TypeTraits<ScalarType>::zero());

  /// computes image of the set after one step/iterate of the dynamical system
  void move(DynSysType & dynsys);
  /// computes image of the set after one step/iterate of the dynamical system and stores it in result
  void move(DynSysType & dynsys, C1DoubletonSet& result) const;

  /// This method computes value of functor f at interval vector represented by this set.
  template<class Functional>
  ScalarType evalAt(const Functional& f) const {
    ScalarType r;
    VectorType x0 = midVector(this->m_currentSet);
    VectorType gradient = f.gradient(this->m_currentSet);
    if(!intersection(f(x0) + C0BaseSet::evalAffineFunctional(gradient,x0), f(this->m_currentSet), r)){
      throw std::logic_error("C1DoubletonSet::evalAt - empty intersection!. Report this error to CAPD developers!");
    }
    return r;
  }

  /// this computes next representation of the set given computed one-step enclosure of the form y + jacPhi*deltaX + rem
//  void moveC0Part(const VectorType& y, const VectorType& deltaX, const MatrixType& jacPhi, const VectorType& rem, C1DoubletonSet& result) const;

  /// this computes next representation of the set given computed one-step enclosure of the form y + jacPhi*deltaX + rem
  static void move(const C1BaseSet& set, C1BaseSet& result, MatrixType& bound, Data& data);

  using SetType::operator VectorType;
  using SetType::operator MatrixType;

  std::string show() const;
  virtual std::string name() const { return "C1DoubletonSet"; }

  using C0BaseSet::get_r0;
  using C0BaseSet::getElement_r0;
  using C0BaseSet::get_x;
  using C0BaseSet::getElement_x;
  using C0BaseSet::get_B;
  using C0BaseSet::get_invB;
  using C0BaseSet::getElement_B;
  using C0BaseSet::getRow_B;
  using C0BaseSet::getColumn_B;
  using C0BaseSet::get_C;
  using C0BaseSet::getElement_C;
  using C0BaseSet::getRow_C;
  using C0BaseSet::getColumn_C;

  using C1BaseSet::m_D;
  using C1BaseSet::m_R;
  using C1BaseSet::m_Bjac;
  using C1BaseSet::m_R0;
  using C1BaseSet::m_Cjac;
  using C1BaseSet::get_D;
  using C1BaseSet::getElement_D;
  using C1BaseSet::getRow_D;
  using C1BaseSet::getColumn_D;
  using C1BaseSet::get_R;
  using C1BaseSet::getElement_R;
  using C1BaseSet::getRow_R;
  using C1BaseSet::getColumn_R;
  using C1BaseSet::get_invBjac;
  using C1BaseSet::get_Bjac;
  using C1BaseSet::getElement_Bjac;
  using C1BaseSet::getRow_Bjac;
  using C1BaseSet::getColumn_Bjac;
  using C1BaseSet::get_R0;
  using C1BaseSet::getElement_R0;
  using C1BaseSet::getRow_R0;
  using C1BaseSet::getColumn_R0;
  using C1BaseSet::get_Cjac;
  using C1BaseSet::getElement_Cjac;
  using C1BaseSet::getRow_Cjac;
  using C1BaseSet::getColumn_Cjac;

protected:
  using C0BaseSet::m_x;
  using C0BaseSet::m_r;
  using C0BaseSet::m_r0;
  using C0BaseSet::m_B;
  using C0BaseSet::m_C;

  using C1BaseSet::m_invBjac;
  // why does it not compile????
/*
  using C1BaseSet::m_D;
  using C1BaseSet::m_R;
  using C1BaseSet::m_R0;
  using C1BaseSet::m_Bjac;
  using C1BaseSet::m_Cjac;
*/
};

/// @}

}} // namespace capd::dynset

#endif // _CAPD_DYNSET_C1DOUBLETONSET_H_
