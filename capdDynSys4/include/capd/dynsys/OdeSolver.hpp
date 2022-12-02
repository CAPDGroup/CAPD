/// @addtogroup dynsys
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file OdeSolver.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_ODESOLVER_HPP_
#define _CAPD_DYNSYS_ODESOLVER_HPP_

#include <sstream>
#include <string>
#include <stdexcept>

#include "capd/basicalg/factrial.h"

#include "capd/dynsys/OdeSolver.h"
#include "capd/dynsys/SolverException.h"
#include "capd/dynsys/DynSys.hpp"
#include "capd/dynsys/BasicOdeSolver.hpp"
#include "capd/dynsys/FirstOrderEnclosure.hpp"
#include "capd/dynsys/HighOrderEnclosure.h"
#include "capd/dynsys/approveRemainder.h"

namespace capd{
namespace dynsys{

//###########################################################//

template<typename MapType, typename StepControlPolicy,typename EnclosurePolicy, typename CurveType>
OdeSolver<MapType, StepControlPolicy,EnclosurePolicy,CurveType>::OdeSolver(VectorFieldType& vectorField, size_type a_order, const StepControlPolicy& stepControl)
  : BaseTaylor(vectorField,a_order,stepControl), implicitCurve(1.0,1.0,vectorField.dimension(),a_order,vectorField.degree())
{}

//###########################################################//

template <typename MapType, typename StepControlPolicy,typename EnclosurePolicy, typename CurveT>
void OdeSolver<MapType, StepControlPolicy, EnclosurePolicy,CurveT>::setOrder(size_type order)
{
  BaseTaylor::setOrder(order);
  this->implicitCurve.setOrder(order);
}

//###########################################################//

template<typename MapType, typename StepControlPolicy,typename EnclosurePolicy, typename CurveType>
typename OdeSolver<MapType,StepControlPolicy,EnclosurePolicy,CurveType>::VectorType
OdeSolver<MapType,StepControlPolicy,EnclosurePolicy,CurveType>::Phi(const ScalarType& t, const VectorType& v)
{
  VectorType result(v.dimension());
  this->setCurrentTime(t);
  this->computeCoefficientsAtCenter(v,this->getOrder());
  BaseTaylor::sumTaylorSeries(result,this->getCoefficientsAtCenter(),this->getOrder());
  return result;
}

//###########################################################//

template<typename MapType, typename StepControlPolicy,typename EnclosurePolicy, typename CurveType>
typename OdeSolver<MapType, StepControlPolicy,EnclosurePolicy,CurveType>::MatrixType
OdeSolver<MapType, StepControlPolicy,EnclosurePolicy,CurveType>::JacPhi(const ScalarType& t, const VectorType &iv)
{
  BaseTaylor::computeCoefficients(t,iv,this->getOrder());

  // the summation of the Taylor series
  MatrixType* matrixCoeff = this->getMatrixCoefficients();
  MatrixType result = matrixCoeff[this->getOrder()];
  for(int r=this->getOrder()-1;r>=0;--r)
    capd::vectalg::multiplyAssignObjectScalarAddObject(result,this->m_step,matrixCoeff[r]);

  return result;
}

//###########################################################//

template<typename MapType, typename StepControlPolicy,typename EnclosurePolicy, typename CurveType>
typename OdeSolver<MapType, StepControlPolicy,EnclosurePolicy,CurveType>::VectorType
OdeSolver<MapType, StepControlPolicy,EnclosurePolicy,CurveType>::Remainder(const ScalarType& t, const VectorType &x, VectorType& o_enc)
{
  capd::poincare::SaveStepControl<OdeSolver> saveStepControl(*this);
  this->turnOffStepControl();
  VectorType rem;
  EnclosurePolicy::computeEnclosureAndRemainder(*this,t,x,o_enc,rem);
  return rem;
}

//###########################################################//

template<typename MapType, typename StepControlPolicy,typename EnclosurePolicy, typename CurveType>
typename OdeSolver<MapType, StepControlPolicy,EnclosurePolicy,CurveType>::VectorType
OdeSolver<MapType, StepControlPolicy,EnclosurePolicy,CurveType>::enclosure(const ScalarType& t, const VectorType &x)
///< the function finds an enclosure for \varphi([0,step],x)
{
  return EnclosurePolicy::enclosure(*this,t,x);
}

//###########################################################//

template<typename MapType, typename StepControlPolicy,typename EnclosurePolicy, typename CurveType>
void OdeSolver<MapType,StepControlPolicy,EnclosurePolicy,CurveType>::encloseC0Map(
      const ScalarType& t,
      const VectorType& x,
      const VectorType& xx,
      VectorType& o_phi,
      VectorType& o_rem,
      VectorType& o_enc,
      MatrixType& o_jacPhi
  )
{
  this->computeTaylorCoefficients(t,x,xx);
  capd::dynsys::computeAndApproveRemainder(*this,t,xx,o_rem,o_enc);
  this->sumTaylorSeries(o_phi,o_jacPhi,this->getCoefficientsAtCenter(),this->getMatrixCoefficients(),this->getOrder());
}

//###########################################################//

template<typename MapType, typename StepControlPolicy,typename EnclosurePolicy, typename CurveType>
void OdeSolver<MapType,StepControlPolicy,EnclosurePolicy,CurveType>::encloseC1Map(
      const ScalarType& t,
      const VectorType& x,
      const VectorType& xx,
      VectorType& o_phi,
      VectorType& o_rem,
      VectorType& o_enc,
      MatrixType& o_jacPhi,
      MatrixType& o_jacRem,
      MatrixType& o_jacEnc
  )
{
  capd::diffAlgebra::C1TimeJet<MatrixType> rem(&o_rem,&o_jacRem);
  capd::diffAlgebra::C1TimeJet<MatrixType> enc(&o_enc,&o_jacEnc);

  this->computeTaylorCoefficients(t,x,xx);
  capd::dynsys::computeAndApproveRemainder(*this,t,xx,rem,enc);
  this->sumTaylorSeries(o_phi,o_jacPhi,this->getCoefficientsAtCenter(),this->getMatrixCoefficients(),this->getOrder());
}

//###########################################################//

template<typename MapType, typename StepControlPolicy,typename EnclosurePolicy, typename CurveType>
typename OdeSolver<MapType, StepControlPolicy,EnclosurePolicy,CurveType>::ScalarType
OdeSolver<MapType, StepControlPolicy,EnclosurePolicy,CurveType>::getCoeffNorm(size_type r, size_type degree) const
{
  typename TypeTraits<ScalarType>::Real result = 0;
  for(size_type i=0;i<this->dimension();++i){
    result = capd::max(result,rightBound(abs(this->remainderCoefficient(i,r))));
    result = capd::max(result,rightBound(abs(this->coefficient(i,r))));
  }
  if(degree)
  {
    for(size_type i=0;i<this->dimension();++i){
      for(size_type j=0;j<this->dimension();++j){
        result = capd::max(result,rightBound(abs(this->remainderCoefficient(i,j,r))));
        result = capd::max(result,rightBound(abs(this->coefficient(i,j,r))));
      }
    }
  }
  return ScalarType(result);
}

}} //namespace capd::dynsys

#endif // _CAPD_DYNSYS_ODESOLVER_HPP_

/// @}
