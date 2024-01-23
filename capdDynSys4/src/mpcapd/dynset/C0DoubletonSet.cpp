/////////////////////////////////////////////////////////////////////////////
/// @file C0DoubletonSet.cpp
///
/// @author kapela
/// Created on: Oct 21, 2009
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2009-2013 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.
#include "capd/vectalg/mplib.h"
#include "capd/dynset/mplib.h"
#include "capd/dynset/C0DoubletonSet.hpp"


#ifdef __HAVE_MPFR__
template class capd::dynset::C0DoubletonSet< capd::MpIMatrix, capd::C0Intv2Policies >;
template class capd::dynset::C0DoubletonSet< capd::MpIMatrix, capd::C0Pped2Policies >;
template class capd::dynset::C0DoubletonSet< capd::MpIMatrix, capd::C0Rect2Policies >;
#endif
