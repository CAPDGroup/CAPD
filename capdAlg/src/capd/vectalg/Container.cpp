/////////////////////////////////////////////////////////////////////////////
/// @file Container.cpp
///
/// @author The CAPD Group
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2013 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#include "capd/vectalg/Dimension.h"
#include "capd/vectalg/Container.hpp"
#include "capd/intervals/lib.h"

namespace capd{
  namespace vectalg{
/// @addtogroup vectalg 
/// @{
template class capd::vectalg::Container<int,CAPD_DEFAULT_DIMENSION>;
template class capd::vectalg::Container<double,CAPD_DEFAULT_DIMENSION>;
template class capd::vectalg::Container<long double,CAPD_DEFAULT_DIMENSION>;
template class capd::vectalg::Container<capd::DInterval,CAPD_DEFAULT_DIMENSION>;

  }
}

