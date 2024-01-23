/////////////////////////////////////////////////////////////////////////////
/// @file InvBByCFactorReorganization.h
///
/// @author kapela
/// Created on: Oct 25, 2009
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2009 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_C1FACTORREORGANIZATION_H_
#define _CAPD_DYNSYS_C1FACTORREORGANIZATION_H_

#include "capd/dynset/reorganization/FactorPolicy.h"
#include "capd/vectalg/Vector_Interval.hpp"

namespace capd{
namespace dynset{
/// @addtogroup dynset 
/// @{
// @addtogroup capd
/// @{

/**
 *   Factor based reorganization for C1 sets.
 *
 *   Works for C1 doubleton sets.
 *   We assume that C1 part is represented as: X + Cjac*R0 + Bjac*R.
 *
 *   Reorganization of C1 part takes place if size of R  is greater then size of R0 times given C1factor.
 *
 *   Reorganization of C0 part depends on C0ReorganizationT base class.
 *
 *   Previously it was built-in into C1Rect2.
 */

template <typename BasePolicy = DefaultPolicy>
class InvBByCFactorReorganization : public FactorPolicy<BasePolicy>{
public:

  template<class Matrix, class Vector>
  bool isReorganizationNeeded(const Matrix& invB, const Matrix& C, const Vector& r, const Vector& r0) const {
     return (this->isReorganizationEnabled() &&
         (capd::vectalg::maxDiam(r) > this->getC0Factor() * capd::vectalg::maxDiam((invB*C)*r0))
         );
  }

  template<class SetType>
  bool isReorganizationNeeded(const SetType & result) const {
     return isReorganizationNeeded(result.get_invB(),result.get_C(),result.get_r(),result.get_r0());
  }

  /// makes reorganization if needed.
  template<class SetType>
  bool reorganizeIfNeeded(SetType & result) const {
    if(this->isReorganizationNeeded(result)){
      this->reorganize(result);
      return true;
    }
    return false;
  }

  template<class Matrix, class Vector>
  bool reorganizeIfNeeded(Matrix& B, Matrix& invB, Vector& r, Matrix& C, Vector& r0) const
  {
    if(this->isReorganizationNeeded(invB,C,r,r0)){
      this->reorganize(B,invB,r,C,r0);
      return true;
    }
    return false;
  }

  template<class SetType>
  bool isC1ReorganizationNeeded(const SetType & result) const {
   return isC1ReorganizationNeeded(result.get_invBjac(), result.get_Cjac(),result.get_R(),result.get_R0());
  }

  template<class Matrix>
  bool isC1ReorganizationNeeded(const Matrix& invB, const Matrix&C, const Matrix& R, const Matrix& R0) const {
   return (this->isReorganizationEnabled() &&
         (capd::vectalg::maxDiam(R) > this->getC1Factor() * capd::vectalg::maxDiam((invB*C)*R0)));
  }

  template<class Matrix>
  bool reorganizeC1IfNeeded(Matrix& B, Matrix& invB, Matrix& R, Matrix& C, Matrix& R0) const{
    if(this->isC1ReorganizationNeeded(invB,C,R,R0)){
      this->reorganizeC1(B,invB,R,C,R0);
      return true;
    }
    return false;
  }

  /// makes reorganization if needed.
  /// return true if reorganization was performed
  template<class SetType>
  bool reorganizeC1IfNeeded(SetType & result) const{
    if(isC1ReorganizationNeeded(result)){
      this->reorganizeC1(result);
      return true;
    }
    return false;
  }

  std::string name() const{
    return "doubleton reorganization (when maxDiam(r) > factor * maxDiam((invB*C)*r0))";
  }
};

}} // capd::dynset
#endif // _CAPD_DYNSYS_C1FACTORREORGANIZATION_H_
