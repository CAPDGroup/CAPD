/////////////////////////////////////////////////////////////////////////////
/// @file SectionDerivativesEnclosure.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2020 by the CAPD Group.
//
// Distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_POINCARE_SECTION_DERIVATIVES_ENCLOSURE_H_
#define _CAPD_POINCARE_SECTION_DERIVATIVES_ENCLOSURE_H_

#include "capd/dynset/C0Set.h"
#include "capd/dynset/C1Set.h"
#include "capd/dynset/C2Set.h"
#include "capd/dynset/CnSet.h"

namespace capd{
/// This namespace contains classes to compute Poincare Maps and Time Maps
namespace poincare{
/// @addtogroup poincare 
/// @{

template<class MatrixT>
struct SectionDerivativesEnclosure{
  typedef MatrixT MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename VectorType::size_type size_type;             ///< integral type used to index containers (vectors, matrices, etc)
  typedef capd::dynset::AbstractSet<VectorType> Set;   ///< type of abstract base class for all sets
  typedef capd::diffAlgebra::Hessian<ScalarType,VectorType::csDim,VectorType::csDim> HessianType;
  typedef capd::diffAlgebra::Jet<MatrixT,0> JetType;

  typedef capd::dynset::C0Set<MatrixType> C0Set;   ///< type of abstract base class for all C0 sets
  typedef capd::dynset::C1Set<MatrixType> C1Set;   ///< type of abstract base class for all C1 sets
  typedef capd::dynset::C2Set<MatrixType> C2Set;   ///< type of abstract base class for all C2 sets
  typedef capd::dynset::CnSet<MatrixType,0> CnSet;   ///< type of abstract base class for all C2 sets

  // These function are used in crossSection only
  // the sets of different types should be computed in a different way
  // the C^0 sets have not defined the MatrixType operators

  template<class Pointer, class Value>
  void checkAndAssignPointer(Pointer& p, const Value& v){
    if(p!=nullptr) *p = v;
  }
  void saveTime(capd::diffAlgebra::TimeRange<ScalarType>& t){ *returnTime = t.getCurrentTime(); }
  void saveEnclosure(C0Set& set){ saveTime(set); }
  void saveEnclosure(C1Set& set){
    this->saveTime(set);
    checkAndAssignPointer(derivativeOfFlow,MatrixType(set));
  }
  void saveEnclosure(C2Set& set){
    this->saveTime(set);
    checkAndAssignPointer(derivativeOfFlow,MatrixType(set));
    checkAndAssignPointer(hessianOfFlow,HessianType(set));
  }
  void saveEnclosure(CnSet& set){
    this->saveTime(set);
    checkAndAssignPointer(jetOfFlow,set.currentSet());
  }

  void updateTime(capd::diffAlgebra::TimeRange<ScalarType>& t){
    *returnTime = intervalHull(*returnTime,t.getCurrentTime());
  }

  void updateEnclosure(C0Set&, C0Set& next){ updateTime(next); }
  void updateEnclosure(C1Set& prev, C1Set& next);
  void updateEnclosure(C2Set& prev, C2Set& next);
  void updateEnclosure(CnSet& prev, CnSet& next);

  template<class SetT, class SolverT>
  void computeOneStepSectionEnclosure(SetT& set, SolverT& solver, VectorType& bound, ScalarType oneStepReturnTime)
  {
    *(this->returnTime) = set.getCurrentTime() + oneStepReturnTime;
    if(this->derivativeOfFlow!=nullptr)
      *(this->derivativeOfFlow) = solver.getCurve().derivative(oneStepReturnTime);
    if(this->hessianOfFlow!=nullptr)
      *(this->hessianOfFlow) = solver.getCurve().hessian(oneStepReturnTime);
    if(this->jetOfFlow!=nullptr){
      solver.getCurve().eval(oneStepReturnTime,*jetOfFlow);
      (*jetOfFlow)() = (typename JetType::VectorType)bound;
    }
  }

  void init(ScalarType* rt, MatrixType* der, HessianType* hess, JetType* jet){
    this->returnTime = rt;
    this->derivativeOfFlow = der;
    this->hessianOfFlow = hess;
    this->jetOfFlow =jet;
  }

  ScalarType* returnTime;
  MatrixType* derivativeOfFlow; ///< stores derivative of the flow (if needed)
  HessianType* hessianOfFlow;   ///< stores hessian of the flow (if needed)
  JetType* jetOfFlow;           ///< stores derivative of the flow (for C^n sets only)
};

template <class MatrixT>
void SectionDerivativesEnclosure<MatrixT>::updateEnclosure(C1Set& prev, C1Set& next)
{
  if(derivativeOfFlow!=nullptr)
  {
    MatrixType oneStepBound = next.getLastMatrixEnclosure()*MatrixType(prev);
    intervalHull(*derivativeOfFlow,oneStepBound,*derivativeOfFlow);
  }
  this->updateTime(next);
}

template <class MatrixT>
void SectionDerivativesEnclosure<MatrixT>::updateEnclosure(C2Set& prev, C2Set& next)
{
  if(this->derivativeOfFlow!=NULL)
  {
    intervalHull(*derivativeOfFlow,next.getLastMatrixEnclosure()*MatrixType(prev),*derivativeOfFlow);

    // computation of the hessian of Poincare map
    if(hessianOfFlow!=NULL)
      intervalHull(*hessianOfFlow,next.getLastMatrixEnclosure()*HessianType(prev) + next.getLastHessianEnclosure()*MatrixType(prev),*hessianOfFlow);
  }
  this->updateTime(next);
}


/*__________________________________________________________________________*/

template <class MatrixT>
void SectionDerivativesEnclosure<MatrixT>::updateEnclosure(CnSet& prev, CnSet& next)
{
  if(jetOfFlow!=nullptr){
    const size_type dim = prev.currentSet().dimension();
    JetType enc(dim, prev.degree());

    substitutionPowerSeries(next.getLastJetEnclosure(), prev.currentSet(),enc, false);

    for (size_type p = 0; p < dim; ++p)
    {
      typename JetType::iterator
          b = jetOfFlow->begin(p, 1),
          e = jetOfFlow->end(p, prev.degree());
      typename JetType::iterator i = enc.begin(p, 1);
      while (b != e)
      {
        *b = intervalHull(*b, *i);
        ++b;
        ++i;
      }
    }
  }
  this->updateTime(next);
}

/// @}
}} // namespace capd::poincare
#endif
