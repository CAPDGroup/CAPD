//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file map/lib.h
///
/// @author Tomasz Kapela   @date 2010-01-23
//
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) Tomasz Kapela 2010
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_MAP_LIB_H_
#define _CAPD_MAP_LIB_H_
#include "capd/vectalg/lib.h"
#include "capd/map/Function.h"
#include "capd/map/Map.h"

namespace capd{

typedef capd::map::Function<IVector> IFunction;
typedef capd::map::Map<IMatrix> IMap;

typedef capd::map::Function<DVector> DFunction;
typedef capd::map::Map<DMatrix> DMap;

typedef capd::map::Function<LDVector> LDFunction;
typedef capd::map::Map<LDMatrix> LDMap;

} // end of namespace capd

#endif // _CAPD_MAP_LIB_H_
