
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

#include "capd/vectalg/lib.h"
#include "capd/poincare/AbstractSection.hpp"
#include "capd/poincare/AffineSection.h"
#include "capd/poincare/CoordinateSection.h"
#include "capd/poincare/NonlinearSection.h"

using namespace capd;

template class poincare::AbstractSection<LDMatrix>;
template class poincare::AbstractSection<DMatrix>;
template class poincare::AbstractSection<IMatrix>;

template class poincare::AffineSection<DMatrix>;
template class poincare::AffineSection<LDMatrix>;
template class poincare::AffineSection<IMatrix>;

template class poincare::CoordinateSection<DMatrix>;
template class poincare::CoordinateSection<LDMatrix>;
template class poincare::CoordinateSection<IMatrix>;

template class poincare::NonlinearSection<DMatrix>;
template class poincare::NonlinearSection<LDMatrix>;
template class poincare::NonlinearSection<IMatrix>;
