/////////////////////////////////////////////////////////////////////////////
/// @file mpcapd/diffAlgebra/C2TimeJet.cpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.wsb-nlu.edu.pl/ for details. 

#include "capd/intervals/mplib.h"
#include "capd/vectalg/mplib.h"
#include "capd/vectalg/Vector.hpp"
#include "capd/vectalg/Matrix.hpp"
#include "capd/diffAlgebra/C2TimeJet.h"

template class capd::diffAlgebra::C2TimeJet<capd::MpMatrix>;
template class capd::diffAlgebra::C2TimeJet<capd::MpIMatrix>;
