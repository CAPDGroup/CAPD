/////////////////////////////////////////////////////////////////////////////
/// @file mpcapd/diffIncl/MultiMap.cpp
///
/// @author Tomasz Kapela
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details. 

#include "capd/map/mplib.h"
#include "capd/diffIncl/MultiMap.h"
#include "capd/map/Map.hpp"
#include "capd/vectalg/Matrix.hpp"

#ifdef __HAVE_MPFR__
  template class capd::diffIncl::MultiMap<capd::MpIMap>;
#endif

