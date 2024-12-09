//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file poincare/lib.h
///
/// @author Tomasz Kapela   @date 2010-01-23
//
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) Tomasz Kapela 2010
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef CAPD_POINCARE_LIB_H
#define CAPD_POINCARE_LIB_H

#include "capd/dynsys/lib.h"
#include "capd/dynset/lib.h"

#include "capd/poincare/PoincareMap.h"
#include "capd/poincare/TimeMap.h"
#include "capd/poincare/AffineSection.h"
#include "capd/poincare/CoordinateSection.h"
#include "capd/poincare/NonlinearSection.h"

#define CAPD_USER_NAMESPACE capd
#include "capd/poincare/typedefs.h"
#undef CAPD_USER_NAMESPACE

#endif // CAPD_POINCARE_LIB_H
