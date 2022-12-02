/////////////////////////////////////////////////////////////////////////////
/// @file PdeSectionDerivativesEnclosure.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2020 by the CAPD Group.
//
// Distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_PDES_PDE_SECTION_DERIVATIVES_ENCLOSURE_H_
#define _CAPD_PDES_PDE_SECTION_DERIVATIVES_ENCLOSURE_H_


namespace capd{
namespace pdes{
/// @addtogroup pdes
/// @{

template<bool isC1>
struct ComputeOneStepSectionEnclosure{
  template<class SetT, class SolverT, class V, class S, class Section>
  void computeOneStepSectionEnclosure(Section& section, SetT& set, SolverT& solver, V& bound, S oneStepReturnTime)
  {
    *(section.returnTime) = set.getCurrentTime() + oneStepReturnTime;
    bound = solver.getCurve()(oneStepReturnTime);
    if(section.derivativeOfFlow!=nullptr){
      *(section.derivativeOfFlow) = solver.getCurve().derivative(oneStepReturnTime);
      section.encDyx = set.encDyx;
      section.encDy = set.encDy;
    }
  }

  template<class Section,class SetT>
  void saveEnclosure(Section& section, SetT& set){
    section.saveTime(set);
    if(section.derivativeOfFlow!=nullptr){
      *(section.derivativeOfFlow) = typename SetT::MatrixType(set);
      section.encDyx = set.Dyx;
      section.encDy = set.Dy;
    }
  }

  template<class SetT, class Section>
  void updateEnclosure(Section& section, SetT& prev, SetT& next)
  {
    if(section.derivativeOfFlow!=nullptr)
    {
      typedef typename SetT::MatrixType MatrixType;
      intervalHull(*(section.derivativeOfFlow),next.encDxx,*(section.derivativeOfFlow));
      section.encDyx = intervalHull(section.encDyx,next.encDyx);
      intervalHull(section.encDy,next.encDy,section.encDy);
    }
    section.updateTime(next);
  }
};

template<>
struct ComputeOneStepSectionEnclosure<false>{
  template<class SetT, class SolverT, class V, class S, class Section>
  void computeOneStepSectionEnclosure(Section& section, SetT& set, SolverT& solver, V& bound, S oneStepReturnTime)
  {
    *(section.returnTime) = set.getCurrentTime() + oneStepReturnTime;
  }
  template<class Section,class SetT>
  void saveEnclosure(Section& section, SetT& set){
    section.saveTime(set);
  }
  template<class SetT, class Section>
  void updateEnclosure(Section& section, SetT& prev, SetT& next)
  {
    section.updateTime(next);
  }
};

template<class VectorT, class MatrixT>
struct PdeSectionDerivativesEnclosure{
  typedef MatrixT MatrixType;
  typedef VectorT VectorType;
  typedef typename MatrixType::RowVectorType FiniteVectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename VectorType::size_type size_type;             ///< integral type used to index containers (vectors, matrices, etc)
  typedef capd::dynset::AbstractSet<FiniteVectorType> Set;   ///< type of abstract base class for all sets

  typedef capd::diffAlgebra::Hessian<ScalarType,FiniteVectorType::csDim,FiniteVectorType::csDim> HessianType;
  typedef capd::diffAlgebra::Jet<MatrixT,0> JetType;

  typedef capd::dynset::C0Set<MatrixType> C0Set;   ///< type of abstract base class for all C0 sets
  typedef capd::dynset::C1Set<MatrixType> C1Set;   ///< type of abstract base class for all C1 sets

  void saveTime(capd::diffAlgebra::TimeRange<ScalarType>& t){ *returnTime = t.getCurrentTime(); }
  void updateTime(capd::diffAlgebra::TimeRange<ScalarType>& t){
    *returnTime = intervalHull(*returnTime,t.getCurrentTime());
  }

  template<class SetT>
  void saveEnclosure(SetT& set){
    ComputeOneStepSectionEnclosure<capd::dynset::SetTraits<typename SetT::SetType>::isC1Set>().saveEnclosure(*this,set);
  }

  template<class SetT>
  void updateEnclosure(SetT& prev, SetT& next)
  {
    ComputeOneStepSectionEnclosure<capd::dynset::SetTraits<typename SetT::SetType>::isC1Set>().updateEnclosure(*this,prev,next);
  }

  template<class SetT, class SolverT>
  void computeOneStepSectionEnclosure(SetT& set, SolverT& solver, VectorType& bound, ScalarType oneStepReturnTime){
    ComputeOneStepSectionEnclosure<capd::dynset::SetTraits<typename SetT::SetType>::isC1Set>().computeOneStepSectionEnclosure(*this,set,solver,bound,oneStepReturnTime);
  }

  // for compatilibilty with ODEs
  void init(ScalarType* rt, MatrixType* der, HessianType* H, JetType* jet){
    this->returnTime = rt;
    this->derivativeOfFlow = der;
    if(H!=nullptr or jet!=nullptr)
      throw std::logic_error("PdeSectionDerivativesEnclosure::init error");
  }

  ScalarType* returnTime;
  MatrixType* derivativeOfFlow; ///< stores derivative of the flow (if needed)
  FiniteVectorType encDy; ///< stores norm of block Vxy and Vyy of derivativeOfFlow
  ScalarType encDyx; ///< stores norm of block Vyx of derivativeOfFlow
};

/// @}
}} // namespace capd::pdes
#endif
