/////////////////////////////////////////////////////////////////////////////
/// @file CnDoubletonSet.cpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2006 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifdef __HAVE_MPFR__

#include "capd/vectalg/mplib.h"
#include "capd/dynset/mplib.h"
#include "capd/dynset/CnDoubletonSet.hpp"

using namespace capd::dynset;

template class capd::dynset::CnDoubletonSet< capd::MpIMatrix, capd::C2Pped2Policies >;
template class capd::dynset::CnDoubletonSet< capd::MpIMatrix, capd::C2Rect2Policies >;

#endif


