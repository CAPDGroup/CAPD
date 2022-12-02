//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file dynset/lib.h
///
/// @author Tomasz Kapela   @date 2010-01-23
//
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) Tomasz Kapela 2010
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSET_LIB_H_
#define _CAPD_DYNSET_LIB_H_

#include "capd/vectalg/lib.h"
#include "capd/dynset/reorganization/InvBByCFactorReorganization.h"
#include "capd/dynset/reorganization/FactorReorganization.h"
#include "capd/dynset/reorganization/QRReorganization.h"
#include "capd/dynset/QRPolicy.h"

#include "capd/dynset/C0BallSet.h"
#include "capd/dynset/C0FlowballSet.h"
#include "capd/dynset/C0AffineSet.h"
#include "capd/dynset/C0DoubletonSet.h"
#include "capd/dynset/C0TripletonSet.h"
#include "capd/dynset/C0HOSet.h"

#include "capd/dynset/C1AffineSet.h"
#include "capd/dynset/C1DoubletonSet.h"
#include "capd/dynset/C11Rect2Set.h"
#include "capd/dynset/C1HOSet.h"

#include "capd/dynset/C2DoubletonSet.h"
#include "capd/dynset/CnDoubletonSet.h"
#include "capd/dynset/CnRect2Set.h"

#define CAPD_USER_NAMESPACE capd
#include "capd/dynset/typedefs.h"
#undef CAPD_USER_NAMESPACE



#endif // _CAPD_DYNSET_LIB_H_
