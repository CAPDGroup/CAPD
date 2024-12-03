/////////////////////////////////////////////////////////////////////////////
/// @file SwapReorganization.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_SWAPREORGANIZATION_H_
#define _CAPD_DYNSYS_SWAPREORGANIZATION_H_

#include "capd/dynset/reorganization/FactorPolicy.h"
#include "capd/matrixAlgorithms/floatMatrixAlgorithms.hpp"
#include "capd/dynset/DefaultPolicy.h"

namespace capd{
namespace dynset{
/// @addtogroup dynset
/// @{

/**
 *   Reorganization is performed if r is bigger than r0 but in coordinate system of r.
 *
 *   Works for doubleton sets, represented as: x + C*r0 + B*r.
 *
 *   Reorganization takes place if size of r  is greater then size of B^{-1}*C*r0 times given factor.
 *
 */
template <typename BasePolicy = DefaultPolicy>
class SwapReorganization: public BasePolicy, public FactorPolicy
{
public:

  template<class SetType>
  void reorganizeIfNeeded(SetType & result) const
  {
    typedef typename SetType::VectorType Vector;
    typedef typename SetType::MatrixType Matrix;
    if(this->isReorganizationEnabled())
    {
      Vector r0 = (result.get_invB()*result.get_C()) * result.get_r0();
      if( (capd::vectalg::size(result.get_r()) > this->getFactor() * capd::vectalg::size(r0)) )
      {
        r0 = result.get_r();
        result.set_r(result.get_r0());
        result.set_r0(r0);

        Matrix  M = result.get_C();
        result.set_C(result.get_B());
        result.set_B(M);
        // this method sets B = invB = Id and r=0
        //result.setToIdentity();
      }
    }
  }

  std::string name() const{
      return "doubleton reorganization (when size(r) > factor * B^-1 * C *size(r0))";
  }
};

/// @}
}} // end of capd::dynset

#endif /* _CAPD_SWAPREORGANIZATION_H_ */
