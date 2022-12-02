/////////////////////////////////////////////////////////////////////////////
/// @file mpcapd/dynsys/OdeSolver.cpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#include "capd/vectalg/mplib.h"
#include "capd/map/mplib.h"
#include "capd/dynsys/BasicOdeSolver.hpp"
#include "capd/dynsys/OdeSolver.hpp"
#include "capd/vectalg/Matrix.hpp"

#ifdef __HAVE_MPFR__

  template class capd::dynsys::OdeSolver<capd::MpIMap, capd::dynsys::ILastTermsStepControl>;
  template class capd::dynsys::OdeSolver<capd::MpIMap, capd::dynsys::IEncFoundStepControl>;

  template class capd::dynsys::BasicOdeSolver<capd::MpIMap, capd::dynsys::ILastTermsStepControl>;
  template class capd::dynsys::BasicOdeSolver<capd::MpIMap, capd::dynsys::IEncFoundStepControl>;
  template class capd::dynsys::BasicOdeSolver<capd::MpIMap, capd::dynsys::DLastTermsStepControl>;

  template class capd::dynsys::BasicOdeSolver<capd::MpMap, capd::dynsys::DLastTermsStepControl>;


#endif
