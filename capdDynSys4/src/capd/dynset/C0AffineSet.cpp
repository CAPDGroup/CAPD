
/////////////////////////////////////////////////////////////////////////////
/// @file C0AffineSet.cpp
///
/// @author kapela
/// Created on: Oct 21, 2009
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2009 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.
#include "capd/vectalg/lib.h"
#include "capd/dynset/lib.h"
#include "capd/dynset/C0AffineSet.hpp"

template class capd::dynset::C0AffineSet< capd::IMatrix, capd::C0PpedPolicies >;
template class capd::dynset::C0AffineSet< capd::IMatrix, capd::C0RectPolicies >;
