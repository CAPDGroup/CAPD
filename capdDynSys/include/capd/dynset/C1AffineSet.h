/////////////////////////////////////////////////////////////////////////////
/// @file C1AffineSet.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSET_C1AFFINESET_H_
#define _CAPD_DYNSET_C1AFFINESET_H_

#include "capd/dynset/C1Set.h"
#include "capd/geomset/CenteredAffineSet.h"
#include "capd/geomset/MatrixAffineSet.h"

namespace capd{
namespace dynset{
/// @addtogroup dynset
/// @{

/////////////////////////////////////////////////////////////////////
///
///  The C1 set is represented as doubleton: x + B*r;
///
///  The class represents derivatives in a doubleton form as described in the paper C^1-Lohner algorithm by Piotr Zgliczy�ski (FoCM 2001).
///////////////////////////////////////////////////////////////////////

template<typename MatrixT, typename Policies>
class C1AffineSet : public Policies, public C1Set<MatrixT>,
                       protected capd::geomset::CenteredAffineSet<MatrixT>,
                       protected capd::geomset::MatrixAffineSet<MatrixT>
{
public:
  typedef MatrixT MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename MatrixType::size_type size_type;
  typedef C1Set<MatrixT> SetType;
  typedef typename C1Set<MatrixT>::DynSysType DynSysType;
  typedef capd::geomset::CenteredAffineSet<MatrixT> C0BaseSet;
  typedef capd::geomset::MatrixAffineSet<MatrixT> C1BaseSet;

  C1AffineSet(const VectorType& x, ScalarType t = TypeTraits<ScalarType>::zero());
  C1AffineSet(const VectorType& x, const VectorType& r, ScalarType t = TypeTraits<ScalarType>::zero());
  C1AffineSet(const VectorType& x, const MatrixType& B, const VectorType& r, ScalarType t = TypeTraits<ScalarType>::zero());
  C1AffineSet(const C0BaseSet & c0part, const C1BaseSet& c1part, ScalarType t = TypeTraits<ScalarType>::zero());

  static C1AffineSet create(const VectorType& x, ScalarType t = TypeTraits<ScalarType>::zero());
  static C1AffineSet create(const VectorType& x, const VectorType& r, ScalarType t = TypeTraits<ScalarType>::zero());
  static C1AffineSet create(const VectorType& x, const MatrixType& B, const VectorType& r, ScalarType t = TypeTraits<ScalarType>::zero());
  static C1AffineSet create(const C0BaseSet & c0part, const C1BaseSet& c1part, ScalarType t = TypeTraits<ScalarType>::zero());

  /// computes image of the set after one step/iterate of the dynamical system
  void move(DynSysType & dynsys);
  /// computes image of the set after one step/iterate of the dynamical system and stores it in result
  void move(DynSysType & dynsys, C1AffineSet& result) const;

  using SetType::operator VectorType;
  using SetType::operator MatrixType;

  std::string show() const;
  virtual std::string name() const { return "C1AffineSet"; }

  /// This method computes value of functor f at interval vector represented by this set.
  template<class Functional>
  ScalarType evalAt(const Functional& f) const {
    ScalarType r;
    VectorType gradient = f.gradient((VectorType)(*this));
    if(!intersection(C0BaseSet::evalAt(f,gradient), f(this->m_currentSet), r)){
      throw std::logic_error("C1AffineSet::evalAt - empty intersection!. Report this error to CAPD developers!");
    }
    return r;
  }
  using C0BaseSet::get_x;
  using C0BaseSet::getElement_x;
  using C0BaseSet::get_B;
  using C0BaseSet::get_invB;
  using C0BaseSet::getElement_B;
  using C0BaseSet::getRow_B;
  using C0BaseSet::getColumn_B;

  using C1BaseSet::m_D;
  using C1BaseSet::m_R;
  using C1BaseSet::m_Bjac;
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

protected:
  using C0BaseSet::m_x;
  using C0BaseSet::m_r;
  using C0BaseSet::m_B;

  using C1BaseSet::m_invBjac;
  // why does it not compile????
/*
  using C1BaseSet::m_D;
  using C1BaseSet::m_R;
  using C1BaseSet::m_Bjac;
s*/
};

/// @}

}} // namespace capd::dynset

#endif // _CAPD_DYNSET_C1AFFINESET_H_
