//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file dynset/mplib.h
///
/// @author Tomasz Kapela   @date 2010-01-22
//
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) Tomasz Kapela 2010
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSET_MPLIB_H_
#define _CAPD_DYNSET_MPLIB_H_

#include "capd/vectalg/mplib.h"

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

#ifdef __HAVE_MPFR__

namespace capd{

/// @addtogroup mpcapdlib
/// @{

typedef capd::dynset::FactorReorganization< capd::dynset::InverseQRPolicy< > > C0Pped2Policies;
typedef capd::dynset::FactorReorganization< capd::dynset::FullQRWithPivoting<> > C0Rect2Policies;
typedef capd::dynset::FactorReorganization<> C0Intv2Policies;
typedef capd::dynset::InverseQRPolicy<> C0PpedPolicies;
typedef capd::dynset::FullQRWithPivoting<>  C0RectPolicies;

typedef capd::dynset::InverseQRPolicy<> C1PpedPolicies;
typedef capd::dynset::FullQRWithPivoting<>  C1RectPolicies;
typedef capd::dynset::QRReorganization<capd::dynset::InverseQRPolicy<> > C1Pped2Policies;
typedef capd::dynset::FactorReorganization<capd::dynset::FullQRWithPivoting<> > C1Rect2Policies;

typedef capd::dynset::QRReorganization<capd::dynset::InverseQRPolicy<> > C2Pped2Policies;
typedef capd::dynset::FactorReorganization<capd::dynset::FullQRWithPivoting<> > C2Rect2Policies;

typedef capd::dynset::C0Set<MpIMatrix> MpC0Set;

typedef capd::dynset::C0DoubletonSet<MpIMatrix,C0Intv2Policies> MpC0Intv2Set;
typedef capd::dynset::C0DoubletonSet<MpIMatrix,C0Pped2Policies> MpC0Pped2Set;
typedef capd::dynset::C0DoubletonSet<MpIMatrix,C0Rect2Policies> MpC0Rect2Set;
typedef capd::dynset::C0TripletonSet<MpIMatrix,C0Rect2Policies> MpC0TripletonSet;
typedef capd::dynset::C0AffineSet<MpIMatrix,C0PpedPolicies> MpC0PpedSet;
typedef capd::dynset::C0AffineSet<MpIMatrix,C0RectPolicies> MpC0RectSet;
typedef capd::dynset::C0HOSet<MpC0Rect2Set> MpC0HORect2Set;
typedef capd::dynset::C0HOSet<MpC0TripletonSet> MpC0HOTripletonSet;

 typedef capd::dynset::C1Set<MpIMatrix> MpC1Set;
 typedef capd::dynset::C1AffineSet<MpIMatrix,C1RectPolicies> MpC1RectSet;
 typedef capd::dynset::C1AffineSet<MpIMatrix,C1PpedPolicies> MpC1PpedSet;
 typedef capd::dynset::C1DoubletonSet<MpIMatrix,C1Rect2Policies> MpC1Rect2Set;
 typedef capd::dynset::C1DoubletonSet<MpIMatrix,C1Pped2Policies> MpC1Pped2Set;
 typedef capd::dynset::C11Rect2Set<MpIMatrix> MpC11Rect2Set;
 typedef capd::dynset::C1HOSet<MpC1Rect2Set> MpC1HORect2Set;
 typedef capd::dynset::C1HOSet<MpC1Pped2Set> MpC1HOPped2Set;

 typedef capd::dynset::C2Set<MpIMatrix> MpC2Set;
 typedef capd::dynset::C2DoubletonSet<MpIMatrix,C1Rect2Policies> MpC2Rect2Set;

 typedef capd::dynset::CnSet<MpIMatrix> MpCnSet;
 typedef capd::dynset::CnRect2Set<MpIMatrix,C2Rect2Policies> MpCnRect2Set;
 typedef capd::dynset::CnRect2Set<MpIMatrix,C2Rect2Policies> MpCnMultiMatrixRect2Set;

 /// @}
} // end of namespace capd

#endif //__HAVE_MPFR__

#endif // _CAPD_DYNSET_MPLIB_H_
