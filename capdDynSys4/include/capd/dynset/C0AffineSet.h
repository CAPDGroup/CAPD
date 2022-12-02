/////////////////////////////////////////////////////////////////////////////
/// @file C0AffineSet.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSET_C0AFFINESET_H_
#define _CAPD_DYNSET_C0AFFINESET_H_

#include "capd/dynset/C0Set.h"
#include "capd/geomset/CenteredAffineSet.h"

namespace capd{
namespace dynset{
/// @addtogroup dynset
/// @{

/////////////////////////////////////////////////////////////////////
///
///  The set is represented as: x + B*r;
///  and is moved by the following method.
///
///  internal representation :
///        B*r  - depending on the QRPolicy
///
///////////////////////////////////////////////////////////////////////
template<typename MatrixT, typename Policies>
class C0AffineSet : public Policies, public C0Set<MatrixT>, public capd::geomset::CenteredAffineSet<MatrixT> {
public:
  typedef MatrixT MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename MatrixType::size_type size_type;
  typedef capd::vectalg::Norm<VectorType,MatrixType> NormType;
  typedef capd::geomset::CenteredAffineSet<MatrixT> BaseSet;
  typedef typename C0Set<MatrixT>::SetType SetType;
  typedef typename C0Set<MatrixT>::DynSysType DynSysType;

  C0AffineSet(const VectorType& x, ScalarType t = TypeTraits<ScalarType>::zero());
  C0AffineSet(const VectorType& x, const VectorType& r, ScalarType t = TypeTraits<ScalarType>::zero());
  C0AffineSet(const VectorType& x, const MatrixType& B, const VectorType& r, ScalarType t = TypeTraits<ScalarType>::zero());

  /// computes image of the set after one step/iterate of the dynamical system
  void move(DynSysType & dynsys);
  /// computes image of the set after one step/iterate of the dynamical system and stores it in result
  void move(DynSysType & dynsys, C0AffineSet& result) const;
  using SetType::operator VectorType;

  /// This method computes value of functor f at interval vector represented by this set.
  template<class Funcitonal>
  ScalarType evalAt(const Funcitonal& f) const {
    ScalarType r;
    VectorType gradient = f.gradient((VectorType)(*this));
    if(!intersection(BaseSet::evalAt(f,gradient), f(this->m_currentSet), r)){
      throw std::logic_error("C0AffineSet::evalAt - empty intersection!. Report this error to CAPD developers!");
    }
    return r;
  }

  std::string show() const;
  std::string name() const { return "C0AffineSet"; }

  using BaseSet::get_x;
  using BaseSet::getElement_x;
  using BaseSet::get_B;
  using BaseSet::get_invB;
  using BaseSet::getElement_B;
  using BaseSet::getRow_B;
  using BaseSet::getColumn_B;

protected:
  using BaseSet::m_x;
  using BaseSet::m_r;
  using BaseSet::m_B;
};

/// @}

}} // namespace capd::dynset

#endif // _CAPD_DYNSET_C0AFFINESET_H_
