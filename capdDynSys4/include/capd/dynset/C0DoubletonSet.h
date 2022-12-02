/////////////////////////////////////////////////////////////////////////////
/// @file C0DoubletonSet.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSET_C0DOUBLETONSET_H_
#define _CAPD_DYNSET_C0DOUBLETONSET_H_

#include "capd/dynset/C0Set.h"
#include "capd/geomset/CenteredDoubletonSet.h"
#include "capd/dynset/DoubletonData.h"

namespace capd{
namespace dynset{
/// @addtogroup dynset
/// @{

/////////////////////////////////////////////////////////////////////
///
///  The set is represented as doubleton: x + C*r0 + B*r;
///  and is moved by the following method.
///
///  internal representation :
///        C*r0 - basic 'Lipschitz part'
///        B*r  - depending on the QRPolicy
///
///////////////////////////////////////////////////////////////////////
template<typename MatrixT, typename Policies>
class C0DoubletonSet : public Policies,
                       public C0Set<MatrixT>,
                       public capd::geomset::CenteredDoubletonSet<MatrixT>,
                       public DoubletonData<MatrixT> {
public:
  typedef MatrixT MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename MatrixType::size_type size_type;
  typedef capd::vectalg::Norm<VectorType,MatrixType> NormType;
  typedef capd::geomset::CenteredDoubletonSet<MatrixT> BaseSet;
  typedef typename C0Set<MatrixT>::SetType SetType;
  typedef typename C0Set<MatrixT>::DynSysType DynSysType;
  typedef DoubletonData<MatrixT> Data;
  typedef Policies Policy;

  C0DoubletonSet(const BaseSet& set, ScalarType t);
  C0DoubletonSet(const VectorType& x, ScalarType t = TypeTraits<ScalarType>::zero());
  C0DoubletonSet(const VectorType& x, const VectorType& r0, ScalarType t = TypeTraits<ScalarType>::zero());
  C0DoubletonSet(const VectorType& x, const MatrixType& C, const VectorType& r0, ScalarType t = TypeTraits<ScalarType>::zero());
  C0DoubletonSet(const VectorType& x, const MatrixType& C, const VectorType& r0, const VectorType& r, ScalarType t = TypeTraits<ScalarType>::zero());
  C0DoubletonSet(const VectorType& x, const MatrixType& C, const VectorType& r0, const MatrixType& B, const VectorType& r, ScalarType t = TypeTraits<ScalarType>::zero());

  /// computes image of the set after one step/iterate of the dynamical system
  void move(DynSysType & dynsys);
  /// computes image of the set after one step/iterate of the dynamical system and stores it in result
  void move(DynSysType & dynsys, C0DoubletonSet& result) const;

  /// this computes next representation of the set given computed one-step enclosure of the form y + jacPhi*deltaX + rem
  static void move(const BaseSet& set, BaseSet& result, VectorType& bound, Data& data);

  /// This method computes value of functor f at interval vector represented by this set.
  template<class Functional>
  ScalarType evalAt(const Functional& f) const {
    ScalarType r = f(this->m_currentSet);
    VectorType gradient = f.gradient(this->m_currentSet);
    if(subset(this->m_x,this->m_currentSet)){
      intersection( r, f(this->m_x) + this->evalAffineFunctional(gradient,this->m_x), r);
    }else{
      VectorType x0 = midVector(this->m_currentSet);
      intersection( r, f(x0) + this->evalAffineFunctional(gradient,x0), r);
    }

    return r;
  }


  using SetType::operator VectorType;

  std::string show() const;
  std::string name() const { return "C0DoubletonSet"; }

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
  using BaseSet::affineTransformation;
protected:
  using BaseSet::m_x;
  using BaseSet::m_r;
  using BaseSet::m_r0;
  using BaseSet::m_B;
  using BaseSet::m_C;
};

/// @}

}} // namespace capd::dynset

#endif // _CAPD_DYNSET_C0DOUBLETONSET_H_
