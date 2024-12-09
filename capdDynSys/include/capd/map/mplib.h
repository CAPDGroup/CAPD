//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file map/mplib.h
///
/// @author Tomasz Kapela   @date 2010-01-22
//
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) Tomasz Kapela 2010
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_MAP_MPLIB_H_
#define _CAPD_MAP_MPLIB_H_

#include "capd/vectalg/mplib.h"
#include "capd/map/Function.h"
#include "capd/map/Map.h"

#ifdef __HAVE_MPFR__

namespace capd{
  typedef capd::map::Function<MpIVector> MpIFunction;
  typedef capd::map::Map<MpIMatrix> MpIMap;

  typedef capd::map::Function<MpVector> MpFunction;
  typedef capd::map::Map<MpMatrix> MpMap;

} // end of namespace capd

#endif //__HAVE_MPFR__

#endif // _CAPD_MAP_MPLIB_H_
