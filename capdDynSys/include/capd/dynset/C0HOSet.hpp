/////////////////////////////////////////////////////////////////////////////
///
/// @file C0HOSet.h
///
/// @author Daniel Wilczak
/// Created on: Dec 28, 2011
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2011 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSET_C0HOSET_HPP_
#define _CAPD_DYNSET_C0HOSET_HPP_

#include "capd/dynset/C0HOSet.h"

namespace capd{
namespace dynset{
// @addtogroup capd
/// @{

template <class BaseSetT>
C0HOSet<BaseSetT>::C0HOSet(const BaseSet & set)
  : SetType((VectorType)set,set.getLastEnclosure(),set.getCurrentTime()),
    Data(set.dimension()),
    predictor(set), corrector(set)
{}

template <class BaseSetT>
C0HOSet<BaseSetT>::C0HOSet(const VectorType& x, ScalarType t)
  : SetType(x,VectorType(x.dimension()),t),
    Data(x.dimension()),
    predictor(x), corrector(x)
{}

template <class BaseSetT>
C0HOSet<BaseSetT>::C0HOSet(const VectorType& x, const VectorType& r0, ScalarType t)
  : SetType(x+r0,VectorType(x.dimension()),t),
    Data(x.dimension()),
    predictor(x,r0), corrector(x,r0)
{}

template <class BaseSetT>
C0HOSet<BaseSetT>::C0HOSet(const VectorType& x, const MatrixType& C, const VectorType& r0, ScalarType t)
  : SetType(x+C*r0,VectorType(x.dimension()),t),
    Data(x.dimension()),
    predictor(x,C,r0), corrector(x,C,r0)
{}

template <class BaseSetT>
C0HOSet<BaseSetT>::C0HOSet(const VectorType& x, const MatrixType& C, const VectorType& r0, const VectorType& r, ScalarType t)
  : SetType(x+C*r0+r,VectorType(x.dimension()),t),
    Data(x.dimension()),
    predictor(x,C,r0,r), corrector(x,C,r0,r)
{}

template <class BaseSetT>
C0HOSet<BaseSetT>::C0HOSet(const VectorType& x, const MatrixType& C,
    const VectorType& r0, const MatrixType& B,
    const VectorType& r,
    ScalarType t
    )
  : SetType(x+C*r0+B*r,VectorType(x.dimension()),t),
    Data(x.dimension()),
    predictor(x,C,r0,B,r), corrector(x,C,r0,B,r)
{}

/// @}
}} // namespace capd::dynset

#endif // _CAPD_DYNSET_C0HOSET_HPP_

