/////////////////////////////////////////////////////////////////////////////
/// @file NoReorganization.h
///
/// @author The CAPD Group
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details. 

#ifndef _CAPD_DYNSET_NOREORGANIZATION_H_
#define _CAPD_DYNSET_NOREORGANIZATION_H_

namespace capd{
namespace dynset{

///
///  Reorganization that does nothing.
///
class NoReorganization  {

public:
  template<class SetType>
  void reorganize(SetType& /*result*/) const{}

  template<class SetType>
  bool isReorganizationNeeded(const SetType& /*result*/) const{
    return false;
  }

  template<class SetType>
  bool reorganizeIfNeeded(SetType& /*result*/) const
  { return false; }

  template<class Matrix1, class Matrix2, class Vector>
  bool reorganizeIfNeeded(Matrix1& /*B*/, Matrix2& /*invB*/, Vector& /*r*/, Matrix1& /*C*/, Vector& /*r0*/) const
  { return false; }

  template<class Matrix>
  bool reorganizeC1IfNeeded(Matrix& /*B*/, Matrix& /*invB*/, Matrix& /*r*/, Matrix& /*C*/, Matrix& /*r0*/) const
  { return false; }

  template<class Matrix, class Hessian>
  bool reorganizeC2IfNeeded(Matrix& /*B*/, Matrix& /*invB*/, Hessian& /*r*/, Matrix& /*C*/, Hessian& /*r0*/) const
  { return false; }

  void disableReorganization(){}

  std::string name() const{
      return "no reorganization";
  }

  std::string toString() const {
       return "";
  }
};

}}

#endif /* _CAPD_DYNSET_NOREORGANIZATION_H_ */
