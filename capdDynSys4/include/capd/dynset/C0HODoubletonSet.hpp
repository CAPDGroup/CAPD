/////////////////////////////////////////////////////////////////////////////
///
/// @file C0HODoubletonSet.hpp
///
/// @author Daniel Wilczak
///
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2013 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSET_C0HODOUBLETONSET_HPP_
#define _CAPD_DYNSET_C0HODOUBLETONSET_HPP_

#include "capd/dynset/C0HODoubletonSet.h"
#include "capd/dynset/C0DoubletonSet.hpp"
#include "capd/geomset/CenteredDoubletonSet.hpp"
#include "capd/basicalg/factrial.h"

namespace capd{
namespace dynset{
// @addtogroup capd
/// @{


template <class MatrixT, class Policies>
C0HODoubletonSet<MatrixT,Policies>::C0HODoubletonSet(BaseSet & set)
  : BaseSet(set), predictor(set)
{}

template <class MatrixT, class Policies>
C0HODoubletonSet<MatrixT,Policies>::C0HODoubletonSet(const VectorType& x, ScalarType t)
  : BaseSet(x,t), predictor(x,t)
{}

template <class MatrixT, class Policies>
C0HODoubletonSet<MatrixT,Policies>::C0HODoubletonSet(const VectorType& x, const VectorType& r0, ScalarType t)
  : BaseSet(x,r0,t), predictor(x,r0,t)
{}

template <class MatrixT, class Policies>
C0HODoubletonSet<MatrixT,Policies>::C0HODoubletonSet(
      const VectorType& x, const MatrixType& C, const VectorType& r0, ScalarType t
    ) : BaseSet(x,C,r0,t), predictor(x,C,r0,t)
{}

template <class MatrixT, class Policies>
C0HODoubletonSet<MatrixT,Policies>::C0HODoubletonSet(
      const VectorType& x, const MatrixType& C, const VectorType& r0, const VectorType& r, ScalarType t
    ) : BaseSet(x,C,r0,r,t), predictor(x,C,r0,r,t)
{}

template <class MatrixT, class Policies>
C0HODoubletonSet<MatrixT,Policies>::C0HODoubletonSet(const VectorType& x, const MatrixType& C,
    const VectorType& r0, const MatrixType& B,
    const VectorType& r,
    ScalarType t
    ) : BaseSet(x,C,r0,B,r,t), predictor(x,C,r0,B,r,t)
{}


/// @}
}} // namespace capd::dynset

#endif // _CAPD_DYNSET_C0HODOUBLETONSET_HPP_

