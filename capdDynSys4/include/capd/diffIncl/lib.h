/////////////////////////////////////////////////////////////////////////////
//
/// @file diffIncl/lib.h
///
/// @author Tomasz Kapela   @date 2010-01-23
//
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) Tomasz Kapela 2010
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DIFFINCL_LIB_H_
#define _CAPD_DIFFINCL_LIB_H_
#include "capd/map/lib.h"
#include "capd/poincare/lib.h"
#include "capd/diffIncl/MultiMap.h"
#include "capd/diffIncl/DiffInclusion.h"
#include "capd/diffIncl/DiffInclusionCW.h"
#include "capd/diffIncl/DiffInclusionLN.h"
#include "capd/diffIncl/InclRect2Set.h"

namespace capd{
typedef capd::diffIncl::MultiMap<capd::IMap> IMultiMap;
typedef capd::diffIncl::DiffInclusion<IMultiMap> DiffInclSolver;
typedef capd::diffIncl::DiffInclusionCW<IMultiMap> CWDiffInclSolver;
typedef capd::diffIncl::DiffInclusionLN<IMultiMap> LNDiffInclSolver;
typedef capd::diffIncl::InclRect2Set<capd::IMatrix> InclRect2Set;
typedef capd::poincare::TimeMap<CWDiffInclSolver> CWTimeMap;
typedef capd::poincare::TimeMap<LNDiffInclSolver> LNTimeMap;
} // end of namespace capd

#endif // _CAPD_DIFFINCL_LIB_H_
