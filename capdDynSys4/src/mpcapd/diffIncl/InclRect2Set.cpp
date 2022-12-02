/////////////////////////////////////////////////////////////////////////////
/// @file mpcapd/diffIncl/InclRect2Set.cpp
///
/// @author Tomasz Kapela
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details. 

#include "capd/vectalg/mplib.h"
#include "capd/vectalg/Matrix.hpp"
#include "capd/vectalg/Norm.hpp"
#include "capd/diffIncl/InclRect2Set.hpp"

#ifdef __HAVE_MPFR__  
template class capd::diffIncl::InclRect2Set<capd::MpIMatrix>;
#endif

