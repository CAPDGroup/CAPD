
/////////////////////////////////////////////////////////////////////////////
/// @file C0HOSet.cpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2013 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#include "capd/vectalg/mplib.h"
#include "capd/dynset/mplib.h"
#include "capd/dynset/C1HOSet.hpp"

#ifdef __HAVE_MPFR__
template class capd::dynset::C1HOSet< capd::MpC1Rect2Set >;
template class capd::dynset::C1HOSet< capd::MpC1Pped2Set >;
#endif

