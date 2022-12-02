/////////////////////////////////////////////////////////////////////////////
/// @file Map.cpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#include "capd/vectalg/mplib.h"
#include "capd/map/Map.hpp"


#ifdef __HAVE_MPFR__
namespace capd{
namespace map{
  template class Map<MpMatrix>;
  template class Map<MpIMatrix>;
}}
#endif

