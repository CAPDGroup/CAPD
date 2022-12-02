
/////////////////////////////////////////////////////////////////////////////
/// @file PoincareMap.cpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2013 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifdef __HAVE_MPFR__

#include "capd/intervals/mplib.h"
#include "capd/vectalg/mplib.h"
#include "capd/map/mplib.h"
#include "capd/dynsys/mplib.h"
#include "capd/poincare/PoincareMap.hpp"
#include "capd/dynsys/FadMap.h"
#include "capd/dynsys/BasicFadOdeSolver.hpp"
#include "capd/dynsys/FadOdeSolver.hpp"

using namespace capd;

template class poincare::BasicPoincareMap<MpTaylor>;
template class poincare::BasicPoincareMap<MpC2Taylor>;
template class poincare::BasicPoincareMap<MpCnTaylor>;
template class poincare::BasicPoincareMap< dynsys::BasicFadOdeSolver<dynsys::LorenzFadMap<MpFloat,0>,dynsys::DLastTermsStepControl > >;

template class poincare::PoincareMap<MpITaylor>;
template class poincare::PoincareMap<MpIC2Taylor>;
template class poincare::PoincareMap<MpICnTaylor>;

template class poincare::PoincareMap< dynsys::FadOdeSolver<dynsys::LorenzFadMap<MpInterval,0>,dynsys::ILastTermsStepControl > >;

#endif
