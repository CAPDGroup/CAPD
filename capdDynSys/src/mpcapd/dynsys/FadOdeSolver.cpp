/////////////////////////////////////////////////////////////////////////////
/// @file mpcapd/dynsys/FadOdeSolver.cpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#include "capd/vectalg/mplib.h"
#include "capd/map/mplib.h"
#include "capd/dynsys/BasicFadOdeSolver.hpp"
#include "capd/dynsys/FadOdeSolver.hpp"
#include "capd/dynsys/FadMap.h"

using namespace capd;

#ifdef __HAVE_MPFR__
  typedef dynsys::LorenzFadMap<MpFloat,0> MpLorenz;
  typedef dynsys::LorenzFadMap<MpInterval,0> IMpLorenz;

  template class dynsys::FadOdeSolver<IMpLorenz>;
//  template class dynsys::BasicFadTaylor<IMpLorenz, dynsys::IMpLastTermsStepControl<capd::MpInterval> >;
  template class dynsys::BasicFadOdeSolver<IMpLorenz, dynsys::ILastTermsStepControl>;
  template class dynsys::BasicFadOdeSolver<IMpLorenz, dynsys::IEncFoundStepControl>;
  template class dynsys::BasicFadOdeSolver<MpLorenz, dynsys::DLastTermsStepControl>;
//  template class dynsys::BasicFadTaylor<MpLorenz, dynsys::MpLastTermsStepControl<capd::MpFloat> >;
#endif
