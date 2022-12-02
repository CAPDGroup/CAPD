//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file poincare/mplib.h
///
/// @author Daniel Wilczak
//
/////////////////////////////////////////////////////////////////////////////
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_POINCARE_MPLIB_H_
#define _CAPD_POINCARE_MPLIB_H_

#include "capd/vectalg/mplib.h"
#include "capd/map/mplib.h"
#include "capd/dynsys/mplib.h"
#include "capd/dynset/mplib.h"

#include "capd/poincare/PoincareMap.h"
#include "capd/poincare/TimeMap.h"
#include "capd/poincare/AffineSection.h"
#include "capd/poincare/CoordinateSection.h"
#include "capd/poincare/NonlinearSection.h"

#ifdef __HAVE_MPFR__

namespace capd{

typedef capd::poincare::AffineSection<capd::MpMatrix> MpAffineSection;
typedef capd::poincare::AffineSection<capd::MpIMatrix> MpIAffineSection;

typedef capd::poincare::CoordinateSection<capd::MpMatrix> MpCoordinateSection;
typedef capd::poincare::CoordinateSection<capd::MpIMatrix> MpICoordinateSection;

typedef capd::poincare::NonlinearSection<capd::MpMatrix> MpNonlinearSection;
typedef capd::poincare::NonlinearSection<capd::MpIMatrix> MpINonlinearSection;

typedef capd::poincare::BasicPoincareMap<capd::MpTaylor> MpPoincareMap;
typedef capd::poincare::BasicPoincareMap<capd::MpC2Taylor> MpC2PoincareMap;
typedef capd::poincare::BasicPoincareMap<capd::MpCnTaylor> MpCnPoincareMap;

typedef capd::poincare::PoincareMap<capd::MpITaylor> MpIPoincareMap;
typedef capd::poincare::PoincareMap<capd::MpIC2Taylor> MpIC2PoincareMap;
typedef capd::poincare::PoincareMap<capd::MpICnTaylor> MpICnPoincareMap;


typedef capd::poincare::TimeMap<capd::MpTaylor> MpTimeMap;
typedef capd::poincare::TimeMap<capd::MpC2Taylor> MpC2TimeMap;
typedef capd::poincare::TimeMap<capd::MpCnTaylor> MpCnTimeMap;

typedef capd::poincare::TimeMap<capd::MpITaylor> MpITimeMap;
typedef capd::poincare::TimeMap<capd::MpIC2Taylor> MpIC2TimeMap;
typedef capd::poincare::TimeMap<capd::MpICnTaylor> MpICnTimeMap;

} // end of namespace capd

#endif //__HAVE_MPFR__

#endif // _CAPD_POINCARE_MPLIB_H_
