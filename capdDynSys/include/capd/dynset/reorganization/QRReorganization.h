/////////////////////////////////////////////////////////////////////////////
/// @file QRReorganization.h
///
/// @author kapela
/// Created on: Oct 23, 2009
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2009 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_QRREORGANIZATION_H_
#define _CAPD_QRREORGANIZATION_H_

#include "capd/dynset/reorganization/FactorPolicy.h"
#include "capd/dynset/DefaultPolicy.h"

namespace capd{
namespace dynset{
// @addtogroup capd
/// @{

/**
 *   During reorganization we orthogonalize B
 */
template <typename BasePolicy = DefaultPolicy>
class QRReorganization : public FactorPolicy<BasePolicy>{

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

  template<class Matrix, class Vector>
  void reorganize(Matrix& B, Matrix& invB, Vector& r, Matrix& /*C*/, Vector& /*r0*/) const
  {
    Matrix Q = B;
    capd::matrixAlgorithms::orthonormalize(Q, r);
    invB = Transpose(Q);
    r = (invB*B) * r;
    B = Q;
  }

  template<class SetType>
  void reorganize(SetType& result) const
  {
    typename SetType::MatrixType Q = result.get_B();
    capd::matrixAlgorithms::orthonormalize(Q, result.get_r());
    result.set_invB(Transpose(Q));
    result.set_r( (result.get_invB()*result.get_B()) * result.get_r() );
    result.set_B(Q);
  }

  template<class Matrix, class Vector>
  bool reorganizeIfNeeded(Matrix& B, Matrix& invB, Vector& r, Matrix& C, Vector& r0) const
  {
    if(this->isReorganizationNeeded(r,r0)){
      this->reorganize(B,invB,r,C,r0);
      return true;
    }
    return false;
  }

  /// makes reorganization if needed.
  /// return true if reorganization was performed
  template<class SetType>
  void reorganizeIfNeeded(SetType & result) const{
    if(!isReorganizationNeeded(result)) return;
    reorganize(result);
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
  void reorganizeC1(Matrix& B, Matrix& invB, Matrix& R, Matrix& /*C*/, Matrix& /*rR0*/) const
  {
    Matrix Q = B;
    capd::matrixAlgorithms::orthonormalize(Q, R);
    invB = Transpose(Q);
    R = (invB*B) * R;
    B = Q;
  }

  template<class SetType>
  void reorganizeC1(SetType& result) const
  {
    typename SetType::MatrixType Q = result.get_Bjac();
    capd::matrixAlgorithms::orthonormalize(Q, result.get_R());
    result.set_invBjac(Transpose(Q));
    result.set_rR( (result.get_invBjac()*result.get_Bjac()) * result.get_R() );
    result.set_Bjac(Q);
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

  std::string name() const{
    return "QR reorganization (it replaces B with Q where B=QR)";
  }
};

/// @}
}} // capd::dynset
#endif // _CAPD_QRREORGANIZATION_H_
