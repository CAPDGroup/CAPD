/////////////////////////////////////////////////////////////////////////////
/// @file C2Curve.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DIFFALGEBRA_C2CURVE_HPP_
#define _CAPD_DIFFALGEBRA_C2CURVE_HPP_

#include <stdexcept>
#include "capd/vectalg/algebraicOperations.hpp"
#include "capd/diffAlgebra/C2Curve.h"
#include "capd/diffAlgebra/Curve.hpp"

namespace capd{
namespace diffAlgebra{
/// @addtogroup diffAlgebra
/// @{

template<class BaseCurveT, bool isInterval>
typename C2Curve<BaseCurveT,isInterval>::HessianType C2Curve<BaseCurveT,isInterval>::hessian(const ScalarType& h) const
{
  if((h>=this->m_left) and (h<=this->m_right))
  {
    int p = this->getOrder();
    HessianType result = this->getHessianCoefficients(p);
    --p;
    for(;p>=0;--p)
      capd::vectalg::multiplyAssignObjectScalarAddObject(result,h,this->getHessianCoefficients(p));
    return result;
  }
  throw this->domainErrorMessage("C2Curve::hessian(h)",h,this->m_left,this->m_right);
}

template<class BaseCurveT>
typename C2Curve<BaseCurveT,true>::HessianType C2Curve<BaseCurveT,true>::oneStepHessian(const ScalarType& h) const
{
  if((h>=this->m_left) and (h<=this->m_right))
  {
    int p = this->getOrder();
    HessianType result = this->getHessianCoefficients(p);
    --p;
    for(;p>=0;--p)
      capd::vectalg::multiplyAssignObjectScalarAddObject(result,h,this->getHessianCoefficients(p));

    ScalarType c = power(h,this->getOrder()+1);
    typename HessianType::iterator i1 = result.begin(), i2 = result.end();
    typename HessianType::const_iterator j1 = this->getHessianRemainderCoefficients(this->getOrder()+1).begin();
    while(i1!=i2)
    {
      *i1 += c*(*j1);
      ++i1;
      ++j1;
    }

    return result;
  }
  throw this->domainErrorMessage("C2Curve::oneStepHessian(h)",h,this->m_left,this->m_right);
}

///@}
}} // namespace capd::diffAlgebra

#endif // _CAPD_DIFFALGEBRA_C2CURVE_H_
