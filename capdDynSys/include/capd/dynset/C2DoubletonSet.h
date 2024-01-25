/////////////////////////////////////////////////////////////////////////////
/// @file C2DoubletonSet.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details. 

#ifndef _CAPD_DYNSET_C2DOUBLETONSET_H_
#define _CAPD_DYNSET_C2DOUBLETONSET_H_

#include "capd/dynset/C2Set.h"
#include "capd/dynsys/C2DynSys.h"
#include "capd/geomset/CenteredDoubletonSet.h"
#include "capd/geomset/MatrixDoubletonSet.h"
#include "capd/dynset/DoubletonData.h"


namespace capd{
namespace dynset{
/// @addtogroup dynset
/// @{

/**
 * C2 set in doubleton form.
 *
 *  C^2-Lohner algorithm.
 */

template<typename MatrixT, class Policies>
class C2DoubletonSet: public Policies, public C2Set<MatrixT>,
                      public capd::geomset::CenteredDoubletonSet<MatrixT>,
                      public capd::geomset::MatrixDoubletonSet<MatrixT>,
                      public C1DoubletonData<MatrixT>
{
public:
  typedef MatrixT MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename MatrixType::size_type size_type;
  typedef C2Set<MatrixT> SetType;
  typedef typename SetType::HessianType HessianType;
  typedef capd::dynsys::C2DynSys<MatrixType> DynSysType;
  typedef capd::geomset::CenteredDoubletonSet<MatrixT> C0BaseSet;
  typedef capd::geomset::MatrixDoubletonSet<MatrixT> C1BaseSet;
  typedef C1DoubletonData<MatrixT> Data;

  C2DoubletonSet(const VectorType& x, ScalarType t = TypeTraits<ScalarType>::zero());
  C2DoubletonSet(const VectorType& x, const VectorType& r0, ScalarType t = TypeTraits<ScalarType>::zero());
  C2DoubletonSet(const VectorType& x, const MatrixType& C, const VectorType& r0, ScalarType t = TypeTraits<ScalarType>::zero());
  C2DoubletonSet(const VectorType& x, const MatrixType& C, const VectorType& r0, const VectorType& r, ScalarType t = TypeTraits<ScalarType>::zero());
  C2DoubletonSet(const VectorType& x, const MatrixType& C, const VectorType& r0, const MatrixType& B, const VectorType& r, ScalarType t = TypeTraits<ScalarType>::zero());
  C2DoubletonSet(const C0BaseSet & c0part, const C1BaseSet& c1part, ScalarType t = TypeTraits<ScalarType>::zero());
  C2DoubletonSet(const C0BaseSet & c0part, const C1BaseSet& c1part, const HessianType& H, ScalarType t = TypeTraits<ScalarType>::zero());

  using SetType::operator VectorType;
  using SetType::operator MatrixType;
  using SetType::operator HessianType;

  void move(DynSysType& c2dynsys);
  void move(DynSysType& c2dynsys, C2DoubletonSet& result);

  /// This method computes value of functor f at interval vector represented by this set.
  template<class Functional>
  ScalarType evalAt(const Functional& f) const {
    ScalarType r;
    VectorType x0 = midVector(this->m_currentSet);
    VectorType gradient = f.gradient(this->m_currentSet);
    if(!intersection(f(x0) + C0BaseSet::evalAffineFunctional(gradient,x0), f(this->m_currentSet), r)){
      throw std::logic_error("C2DoubletonSet::evalAt - empty intersection!. Report this error to CAPD developers!");
    }
    return r;
  }

  virtual std::string name() const { return "C2DoubletonSet"; }
  std::string show(void) const;

protected:

// C^2 part represented as alpha + Bhess*HR + Chess*HR0

  MatrixType m_Chess, m_Bhess, m_invBhess;
  HessianType m_alpha, m_HR, m_HR0;
};
/// @}

// -----------------------------------------------------------------------------

template<typename MatrixType, class Policies>
inline void C2DoubletonSet<MatrixType,Policies>::move(DynSysType& c2dynsys)
{
  move(c2dynsys,*this);
}

}} // the end of the namespace capd::dynset

#endif // _CAPD_DYNSET_C2DOUBLETONSET_H_
