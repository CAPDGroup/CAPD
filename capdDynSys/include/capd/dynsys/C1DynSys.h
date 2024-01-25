

/////////////////////////////////////////////////////////////////////////////
/// @file C1DynSys.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2015 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_C1DYNSYS_H_
#define _CAPD_DYNSYS_C1DYNSYS_H_

#include <string>
#include "capd/dynsys/DynSys.h"

namespace capd{
namespace dynsys{
/// @addtogroup dynsys
/// @{

template<typename MatrixType>
class C1DynSys : public capd::dynsys::DynSys<MatrixType>{
public:
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;

  virtual void encloseC1Map(
      const ScalarType& t,
      const VectorType& x,
      const VectorType& xx,
      VectorType& o_phi,
      VectorType& o_rem,
      VectorType& o_enc,
      MatrixType& o_jacPhi,
      MatrixType& o_jacRem,
      MatrixType& o_jacEnc
  ) = 0;
};
/// @}
}} // namespace capd::dynsys

#endif // _CAPD_DYNSYS_C1DYNSYS_H_

