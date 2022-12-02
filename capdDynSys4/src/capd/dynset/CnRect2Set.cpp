
/////////////////////////////////////////////////////////////////////////////
/// @file CnRect2Set.cpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2006 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details. 

//#include "capd/capdlib.h"
#include "capd/vectalg/lib.h"
#include "capd/dynset/lib.h"
#include "capd/dynset/CnRect2Set.hpp"

template class capd::dynset::CnRect2Set< capd::IMatrix, capd::C2Pped2Policies >;
template class capd::dynset::CnRect2Set< capd::IMatrix, capd::C2Rect2Policies >;
