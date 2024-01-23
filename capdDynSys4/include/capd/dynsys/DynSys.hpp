

/////////////////////////////////////////////////////////////////////////////
/// @file DynSys.hpp
///
/// @author The CAPD Group
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.wsb-nlu.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_DYNSYS_HPP_
#define _CAPD_DYNSYS_DYNSYS_HPP_

#include "capd/dynsys/DynSys.h"
#include "capd/vectalg/Norm.hpp"

namespace capd{
namespace dynsys{
/// @addtogroup dynsys
/// @{

/// this function returns Lipschitz contants for maps. It should be overriden in classes implementing flows.
template<typename MatrixType>
typename DynSys<MatrixType>::ScalarType
DynSys<MatrixType>::Lipschitz(const ScalarType& t, const VectorType &iv, NormType &n)
{
  return n(JacPhi(t,iv)).rightBound();
}
/// @}
}} // namespace capd::dynsys

#endif // _CAPD_DYNSYS_DYNSYS_HPP_


