/////////////////////////////////////////////////////////////////////////////
///
/// @file C1HOSet.hpp
///
/// @author Daniel Wilczak
/// Created on: Apr 15, 2014
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSET_C1HOSET_HPP_
#define _CAPD_DYNSET_C1HOSET_HPP_

#include "capd/dynset/C1HOSet.h"

namespace capd{
namespace dynset{
// @addtogroup capd
/// @{
/*
template <class BaseSetT>
C1HOSet<BaseSetT>::C1HOSet(const BaseSet & set)
  : SetType((VectorType)set,set.getLastEnclosure(),set.getCurrentTime()),
    HOBaseSet<MatrixType>(set.dimension()),
    predictor(set), corrector(set)
{}
*/
template <class BaseSetT>
C1HOSet<BaseSetT>::C1HOSet(const C0BaseSet& c0part, const C1BaseSet& c1part, ScalarType t)
  : SetType(
      VectorType(c0part),
      VectorType(c0part.dimension()),
      MatrixType(c1part),
      MatrixType(c0part.dimension(),c0part.dimension()),
      t),
    C1BaseSet(c1part),
    Data(c0part.dimension()),
    predictor(c0part), corrector(c0part)
{}

template <class BaseSetT>
C1HOSet<BaseSetT>::C1HOSet(const VectorType& x, ScalarType t)
  : SetType(
      x,
      VectorType(x.dimension()),
      MatrixType::Identity(x.dimension()),
      MatrixType(x.dimension(),x.dimension()),
      t),
    C1BaseSet(x.dimension()),
    Data(x.dimension()),
    predictor(x), corrector(x)
{}

template <class BaseSetT>
C1HOSet<BaseSetT>::C1HOSet(const VectorType& x, const VectorType& r0, ScalarType t)
  : SetType(
      x+r0,
      VectorType(x.dimension()),
      MatrixType::Identity(x.dimension()),
      MatrixType(x.dimension(),x.dimension()),
      t),
    C1BaseSet(x.dimension()),
    Data(x.dimension()),
    predictor(x,r0), corrector(x,r0)
{}

template <class BaseSetT>
C1HOSet<BaseSetT>::C1HOSet(const VectorType& x, const MatrixType& C, const VectorType& r0, ScalarType t)
  : SetType(
      x+C*r0,
      VectorType(x.dimension()),
      MatrixType::Identity(x.dimension()),
      MatrixType(x.dimension(),x.dimension()),
      t),
    C1BaseSet(x.dimension()),
    Data(x.dimension()),
    predictor(x,C,r0), corrector(x,C,r0)
{}

template <class BaseSetT>
C1HOSet<BaseSetT>::C1HOSet(const VectorType& x, const MatrixType& C, const VectorType& r0, const VectorType& r, ScalarType t)
  : SetType(
      x+C*r0+r,
      VectorType(x.dimension()),
      MatrixType::Identity(x.dimension()),
      MatrixType(x.dimension(),x.dimension()),
      t),
    C1BaseSet(x.dimension()),
    Data(x.dimension()),
    predictor(x,C,r0,r), corrector(x,C,r0,r)
{}

template <class BaseSetT>
C1HOSet<BaseSetT>::C1HOSet(const VectorType& x, const MatrixType& C,
    const VectorType& r0, const MatrixType& B,
    const VectorType& r,
    ScalarType t
    )
  : SetType(
      x+C*r0+B*r,
      VectorType(x.dimension()),
      MatrixType::Identity(x.dimension()),
      MatrixType(x.dimension(),x.dimension()),
      t),
    C1BaseSet(x.dimension()),
    Data(x.dimension()),
    predictor(x,C,r0,B,r), corrector(x,C,r0,B,r)
{}

/// @}
}} // namespace capd::dynset

#endif // _CAPD_DYNSET_C1HOSET_HPP_

