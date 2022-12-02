
/////////////////////////////////////////////////////////////////////////////
/// @file Function.cpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.wsb-nlu.edu.pl/ for details. 

#include "capd/vectalg/lib.h"
#include "capd/map/Function.hpp"
template class capd::map::BasicFunction<double>;
template class capd::map::BasicFunction<long double>;
template class capd::map::BasicFunction<capd::DInterval>;
template class capd::map::Function<capd::DVector>;
template class capd::map::Function<capd::LDVector>;
template class capd::map::Function<capd::IVector>;

