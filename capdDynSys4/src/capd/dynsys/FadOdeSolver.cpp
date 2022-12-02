/////////////////////////////////////////////////////////////////////////////
/// @file FadOdeSolver.cpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.


#include "capd/vectalg/lib.h"
#include "capd/map/lib.h"
#include "capd/dynsys/FadMap.h"
#include "capd/dynsys/BasicFadOdeSolver.hpp"
#include "capd/dynsys/FadOdeSolver.hpp"

using namespace capd;

template class dynsys::FadOdeSolver<dynsys::LorenzFadMap<capd::DInterval,0>, dynsys::IEncFoundStepControl >;
template class dynsys::FadOdeSolver<dynsys::LorenzFadMap<capd::DInterval,0>,dynsys::ILastTermsStepControl >;

template class capd::dynsys::BasicFadOdeSolver<dynsys::LorenzFadMap<capd::DInterval,0>,dynsys::IEncFoundStepControl>;
template class capd::dynsys::BasicFadOdeSolver<dynsys::LorenzFadMap<capd::DInterval,0>,dynsys::ILastTermsStepControl>;

template class capd::dynsys::BasicFadOdeSolver<dynsys::LorenzFadMap<double,0> >;
template class capd::dynsys::BasicFadOdeSolver<dynsys::LorenzFadMap<long double,0> >;
