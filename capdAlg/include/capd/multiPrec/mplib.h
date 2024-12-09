//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file multiPrec/mplib.h
///
/// @author Tomasz Kapela   @date 2010-01-22
//
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) Tomasz Kapela 2010
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef CAPD_MULTIPREC_MPLIB_H
#define CAPD_MULTIPREC_MPLIB_H

#ifdef __HAVE_MPFR__
#include "capd/multiPrec/MpReal.h"
#include "capd/multiPrec/MpInt.h"

namespace capd{
  typedef ::capd::multiPrec::MpReal MpFloat;
  typedef ::capd::multiPrec::MpInt MpInt;
} // end of namespace capd

#endif //__HAVE_MPFR__

#endif // CAPD_MULTIPREC_MPLIB_H
