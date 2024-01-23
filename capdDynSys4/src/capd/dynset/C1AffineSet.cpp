
/////////////////////////////////////////////////////////////////////////////
/// @file C1AffineSet.cpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2013 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#include "capd/vectalg/lib.h"
#include "capd/dynset/lib.h"
#include "capd/dynset/C1AffineSet.hpp"

using namespace capd::dynset;

template class capd::dynset::C1AffineSet< capd::IMatrix , capd::C1PpedPolicies >;
template class capd::dynset::C1AffineSet< capd::IMatrix , capd::C1RectPolicies >;
