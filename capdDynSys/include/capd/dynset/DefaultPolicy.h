/////////////////////////////////////////////////////////////////////////////
/// @file DefaultPolicy.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSET_DEFAULT_POLICY_H_
#define _CAPD_DYNSET_DEFAULT_POLICY_H_

#include "capd/dynset/reorganization/NoReorganization.h"

namespace capd{
namespace dynset{
/// @addtogroup dynset 
/// @{
// @addtogroup capd
/// @{

class IdQRPolicy
{
public:
  template<class VectorT, class MatrixT>
  void computeBinvB(MatrixT& B, MatrixT& invB, const VectorT& /*v*/) const
  {
    B.setToIdentity();
    invB.setToIdentity();
  }
  std::string toString() const {
    return "";
  }
};

class DefaultPolicy : public IdQRPolicy, public NoReorganization{
  public:
  std::string toString() const {
    return "";
  }
};


/// @}
}} // namespace capd::dynset

#endif // _CAPD_DYNSET_DEFAULT_POLICY_H_

