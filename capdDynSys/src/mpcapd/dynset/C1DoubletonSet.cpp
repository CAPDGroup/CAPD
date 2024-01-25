// @addtogroup capd


/////////////////////////////////////////////////////////////////////////////
/// @file mpcapd/dynset/C1DoubletonSet.cpp
///
/// @author kapela
/// Created on: Oct 24, 2009
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2009 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.wsb-nlu.edu.pl/ for details.

#include "capd/vectalg/mplib.h"
#include "capd/dynset/mplib.h"
#include "capd/dynset/C1DoubletonSet.hpp"


#ifdef __HAVE_MPFR__
template class capd::dynset::C1DoubletonSet< capd::MpIMatrix , capd::C1Pped2Policies >;
template class capd::dynset::C1DoubletonSet< capd::MpIMatrix , capd::C1Rect2Policies >;
#endif



