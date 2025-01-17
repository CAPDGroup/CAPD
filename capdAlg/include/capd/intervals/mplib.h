/////////////////////////////////////////////////////////////////////////////
//
/// @file intervals/mplib.h
///
/// @author Tomasz Kapela   @date 2010-01-22
//
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) Tomasz Kapela 2010
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef CAPD_INTERVAL_MPLIB_H
#define CAPD_INTERVAL_MPLIB_H

#include "capd/multiPrec/mplib.h"
#include "capd/intervals/MpInterval.h"

#ifdef __HAVE_MPFR__

namespace capd{
 typedef capd::intervals::MpInterval MpInterval;
} // end of namespace capd

#endif //__HAVE_MPFR__

#endif // CAPD_INTERVAL_MPLIB_H
