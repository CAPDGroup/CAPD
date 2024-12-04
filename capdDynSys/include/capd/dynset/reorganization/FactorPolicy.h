/////////////////////////////////////////////////////////////////////////////
/// @file FactorPolicy.h
///
/// @author kapela
/// Created on: Oct 21, 2009
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2009 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_FACTORPOLICY_H_
#define _CAPD_FACTORPOLICY_H_

#include <sstream>
#include "capd/dynset/DefaultPolicy.h"
#include "capd/vectalg/Matrix.h"

namespace capd{
namespace dynset{

template <typename BasePolicy = DefaultPolicy>
class FactorPolicy : public BasePolicy{
public:
  FactorPolicy()
    : m_c0Factor(20.0),
      m_c1Factor(20.0),
      m_c2Factor(20.0),
      m_onoff(true)
  {}

  template<class SetType>
  void reorganize(SetType & result) const
  {
    result.set_r0(result.get_r() + (result.get_invB() * result.get_C()) * result.get_r0());
    result.set_C(result.get_B());
    // this method sets B = invB = Id and r=0
    result.setToIdentity();
  }

  template<class Matrix1, class Matrix2, class Vector>
  void reorganize(Matrix1& B, Matrix2& invB, Vector& r, Matrix1& C, Vector& r0) const
  {
    r0 = r + capd::vectalg::matrixByMatrix<Matrix2>(invB,C) * r0;
    C = B;
    B.setToIdentity();
    invB.setToIdentity();
    r.clear();
  }

  template<class SetType>
  void reorganizeC1(SetType & result) const{
    result.set_R0(result.get_R() + (result.get_invBjac()*result.get_Cjac()) * result.get_R0());
    result.set_Cjac(result.get_Bjac());
    result.setToIdentity();
  }

  template<class Matrix>
  void reorganizeC1(Matrix& B, Matrix& invB, Matrix& R, Matrix& C, Matrix& R0) const{
    R0 =  R + (invB*C) * R0;
    C = B;
    B.setToIdentity();
    invB.setToIdentity();
    R.clear();
  }

  template<class Matrix, class HessianT>
  void reorganizeC2(Matrix& B, Matrix& invB, HessianT& R, Matrix& C, HessianT& R0) const{
    R0 =  R + (invB*C) * R0;
    C = B;
    B.setToIdentity();
    invB.setToIdentity();
    R.clear();
  }

  /// sets c0,c1,c2 factors to new value
  void setFactor(double A_factor){
    m_c0Factor = A_factor;
    m_c1Factor = A_factor;
    m_c2Factor = A_factor;
  }

  /// sets new value of c0Factor;
  void setC0Factor(double A_factor){
    m_c0Factor = A_factor;
  }

  /// sets new value of c1Factor;
  void setC1Factor(double A_factor){
    m_c1Factor = A_factor;
  }

  /// sets new value of c2Factor;
  void setC2Factor(double A_factor){
    m_c2Factor = A_factor;
  }

  /// returns current value of c0Factor;
  double getC0Factor() const{
    return m_c0Factor;
  }

  /// returns current value of c1Factor;
  double getC1Factor() const{
    return m_c1Factor;
  }

  /// returns current value of c2Factor;
  double getC2Factor() const{
    return m_c2Factor;
  }

  /// sets the flag which controls possibility of reorganization to new value
  void onoffReorganization(bool flag){
    m_onoff = flag;
  }

  /// returns current value of the flag which controls possibility of reorganization
  bool isReorganizationEnabled() const{
    return m_onoff;
  }

  std::string toString() const {
     std::ostringstream str;
     str << " m_c0Factor = " << m_c0Factor;
     str << ", m_c1Factor = " << m_c1Factor;
     str << ", m_c2Factor = " << m_c2Factor;
     return str.str();
  }

protected:
  double m_c0Factor;
  double m_c1Factor;
  double m_c2Factor;
  bool m_onoff;
};

}}
#endif // _CAPD_FACTORPOLICY_H_
