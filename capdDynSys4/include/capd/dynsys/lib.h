//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file dynsys/lib.h
///
/// @author Tomasz Kapela   @date 2010-01-23
//
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) Tomasz Kapela 2010
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_LIB_H_
#define _CAPD_DYNSYS_LIB_H_

#include "capd/basicalg/factrial.h"
#include "capd/dynsys/SolverException.h"
#include "capd/dynsys/OdeSolver.h"
#include "capd/dynsys/C2OdeSolver.h"
#include "capd/dynsys/CnOdeSolver.h"
#include "capd/vectalg/lib.h"
#include "capd/map/lib.h"

#define CAPD_USER_NAMESPACE capd
#include "capd/dynsys/typedefs.h"
#undef CAPD_USER_NAMESPACE

#endif // _CAPD_DYNSYS_LIB_H_
