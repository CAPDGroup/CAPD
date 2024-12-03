
/////////////////////////////////////////////////////////////////////////////
/// @file FactorReorganization.h
///
/// @author The CAPD Group
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSET_FACTORREORGANIZATION_H_
#define _CAPD_DYNSET_FACTORREORGANIZATION_H_

#include "capd/matrixAlgorithms/floatMatrixAlgorithms.hpp"
#include "capd/dynset/reorganization/FactorPolicy.h"
#include "capd/diffAlgebra/Hessian.h"
#include "capd/vectalg/Vector_Interval.hpp"

namespace capd{
namespace dynset{

/**
 *   Factor based reorganization.
 *   It reorganizes a doubleton set when size(r) > factor * size(r0)
 */

template <typename BasePolicy = DefaultPolicy>
class FactorReorganization : public FactorPolicy<BasePolicy>{
public:

  template<class Vector>
  bool isReorganizationNeeded(const Vector& r, const Vector& r0) const {
     return (this->isReorganizationEnabled() &&
         (capd::vectalg::maxDiam(r) > this->getC0Factor() * capd::vectalg::maxDiam(r0))
         );
  }

  template<class SetType>
  bool isReorganizationNeeded(const SetType & result) const {
     return isReorganizationNeeded(result.get_r(),result.get_r0());
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

  template<class Matrix1, class Matrix2, class Vector>
  bool reorganizeIfNeeded(Matrix1& B, Matrix2& invB, Vector& r, Matrix1& C, Vector& r0) const
  {
    if(this->isReorganizationNeeded(r,r0)){
      this->reorganize(B,invB,r,C,r0);
      return true;
    }
    return false;
  }

  template<class SetType>
  bool isC1ReorganizationNeeded(const SetType & result) const {
   return isC1ReorganizationNeeded(result.get_R(),result.get_R0());
  }

  template<class Matrix>
  bool isC1ReorganizationNeeded(const Matrix& R, const Matrix& R0) const {
   return (this->isReorganizationEnabled() &&
         (capd::vectalg::maxDiam(R) > this->getC1Factor() * capd::vectalg::maxDiam(R0)));
  }

  template<class Matrix>
  bool reorganizeC1IfNeeded(Matrix& B, Matrix& invB, Matrix& R, Matrix& C, Matrix& R0) const{
    if(this->isC1ReorganizationNeeded(R,R0)){
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

  template<class HessianT>
  bool isC2ReorganizationNeeded(const HessianT& R, const HessianT& R0) const {
   return (this->isReorganizationEnabled() &&
         (capd::vectalg::maxDiam(R) > this->getC2Factor() * capd::vectalg::maxDiam(R0)));
  }

  template<class MatrixT, class HessianT>
  bool reorganizeC2IfNeeded(MatrixT& B, MatrixT& invB, HessianT& R, MatrixT& C, HessianT& R0) const{
    if(isC2ReorganizationNeeded(R,R0)){
      this->reorganizeC2(B,invB,R,C,R0);
      return true;
    }
    return false;

  }

  std::string name() const{
    return "doubleton reorganization (when maxDiam(r) > factor * maxDiam(r0))";
  }
};

}} // capd::dynset

#endif /* FACTORREORGANIZE_H_ */
