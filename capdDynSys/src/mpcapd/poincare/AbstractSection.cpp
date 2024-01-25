
/////////////////////////////////////////////////////////////////////////////
/// @file AbstractSection.cpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2013 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifdef __HAVE_MPFR__

#include "capd/vectalg/mplib.h"
#include "capd/poincare/AbstractSection.h"
#include "capd/poincare/AbstractSection.hpp"
#include "capd/poincare/AffineSection.h"
#include "capd/poincare/CoordinateSection.h"
#include "capd/poincare/NonlinearSection.h"

using namespace capd;

template class poincare::AbstractSection<MpMatrix>;
template class poincare::AbstractSection<MpIMatrix>;

template class poincare::AffineSection<MpMatrix>;
template class poincare::AffineSection<MpIMatrix>;

template class poincare::CoordinateSection<MpMatrix>;
template class poincare::CoordinateSection<MpIMatrix>;

template class poincare::NonlinearSection<MpMatrix>;
template class poincare::NonlinearSection<MpIMatrix>;

#endif
