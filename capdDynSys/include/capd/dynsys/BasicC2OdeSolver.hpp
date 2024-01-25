

/////////////////////////////////////////////////////////////////////////////
/// @file BasicC2OdeSolver.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_BASICC2ODESOLVER_HPP_
#define _CAPD_DYNSYS_BASICC2ODESOLVER_HPP_

#include <string>
#include <stdexcept>

#include "capd/dynsys/BasicC2OdeSolver.h"
#include "capd/dynsys/BasicOdeSolver.hpp"
#include "capd/diffAlgebra/BasicC2Curve.hpp"
#include "capd/diffAlgebra/C2Curve.hpp"

namespace capd{
namespace dynsys{
/// @addtogroup dynsys
/// @{

//###########################################################//

template <typename MapType, typename StepControlType,typename CurveT>
BasicC2OdeSolver<MapType,StepControlType,CurveT>::BasicC2OdeSolver(
      VectorFieldType& vectorField,
      size_type order,
      const StepControlType& stepControl
  ) : BaseTaylor(vectorField,order,stepControl)
{
  if(this->getVectorField().degree()<2)
    this->getVectorField().setDegree(2);
}


//###########################################################//

template <typename MapType, typename StepControlType,typename CurveT>
typename BasicC2OdeSolver<MapType,StepControlType,CurveT>::VectorType
BasicC2OdeSolver<MapType,StepControlType,CurveT>::operator()(VectorType v,MatrixType& der, HessianType& hessian)
{
  VectorType* coeff = this->getCoefficientsAtCenter();
  MatrixType* matrixCoeff = this->getMatrixCoefficients();
  HessianType* hessianCoeff = this->getHessianCoefficients();
  coeff[0] = v;
  matrixCoeff[0].setToIdentity();
  hessianCoeff[0].clear();

  this->m_vField->computeODECoefficients(coeff,matrixCoeff,hessianCoeff,this->getOrder());
  this->computeTimeStep(v);

  this->sumTaylorSeries(v,der,hessian,coeff,matrixCoeff,hessianCoeff,this->getOrder());
  return v;
}

//###########################################################//

template <typename MapType, typename StepControlType,typename CurveT>
typename BasicC2OdeSolver<MapType,StepControlType,CurveT>::VectorType
BasicC2OdeSolver<MapType,StepControlType,CurveT>::operator()(
      VectorType x, const MatrixType& D, const HessianType& H,
      MatrixType& out_der, HessianType& out_hessian
    )
{
  VectorType* coeff = this->getCoefficientsAtCenter();
  MatrixType* matrixCoeff = this->getMatrixCoefficients();
  HessianType* hessianCoeff = this->getHessianCoefficients();
  coeff[0] = x;
  matrixCoeff[0] = D;
  hessianCoeff[0] = H;

  this->m_vField->computeODECoefficients(coeff,matrixCoeff,hessianCoeff,this->getOrder());
  this->computeTimeStep(x);

  this->sumTaylorSeries(x,out_der,out_hessian,coeff,matrixCoeff,hessianCoeff,this->getOrder());
  return x;
}

//###########################################################//

template <typename MapType, typename StepControlType,typename CurveT>
void BasicC2OdeSolver<MapType,StepControlType,CurveT>::sumTaylorSeries(
      VectorType& v,
      MatrixType& der,
      HessianType& hessian,
      VectorType* coeff,
      MatrixType* matrixCoeff,
      HessianType* hessianCoeff,
      size_type order
  )
{
  capd::vectalg::evalPolynomial(coeff,this->m_step,v,order);
  if(this->getMask()==0){
    capd::vectalg::evalPolynomial(matrixCoeff,this->m_step,der,order);
    capd::vectalg::evalPolynomial(hessianCoeff,this->m_step,hessian,order);
  }
  else {
    der = matrixCoeff[order];
    hessian = hessianCoeff[order];
    for(size_type j=0;j<this->dimension();++j){
      if(this->getMask(j)){
        for(int r=order-1;r>=0;--r){
          for(size_type i=1;i<=this->dimension();++i){
            der(i,j+1) = der(i,j+1)*(this->m_step) + matrixCoeff[r](i,j+1);
          }
        }
        for(size_type c = j;c<this->dimension();++c){
          if(this->getMask(j,c)){
            for(int r=order-1;r>=0;--r){
              for(size_type i=0;i<this->dimension();++i){
                hessian(i,j,c) = hessian(i,j,c)*(this->m_step) + hessianCoeff[r](i,j,c);
              }
            }
          }
        }
      }
    }
  }
}
/// @}
}} // namespace capd::dynsys

#endif // _CAPD_DYNSYS_BASICC2ODESOLVER_HPP_


