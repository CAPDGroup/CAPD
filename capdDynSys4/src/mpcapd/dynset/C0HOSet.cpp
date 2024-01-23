
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
#include "capd/dynset/C0HOSet.hpp"

#ifdef __HAVE_MPFR__
template class capd::dynset::C0HOSet< capd::MpC0Rect2Set >;
template class capd::dynset::C0HOSet< capd::MpC0TripletonSet >;
#endif

