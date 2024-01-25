
/////////////////////////////////////////////////////////////////////////////
/// @file TimeMap.cpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2013 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifdef __HAVE_MPFR__

#include "capd/map/mplib.h"
#include "capd/dynsys/mplib.h"
#include "capd/poincare/TimeMap.hpp"
#include "capd/dynsys/FadMap.h"
#include "capd/dynsys/BasicFadOdeSolver.hpp"
#include "capd/dynsys/FadOdeSolver.hpp"
#include "capd/vectalg/Matrix.hpp"

using namespace capd;

template class poincare::TimeMap<MpITaylor>;
template class poincare::TimeMap<MpIC2Taylor>;
template class poincare::TimeMap<MpICnTaylor>;

template class poincare::TimeMap<MpTaylor>;
template class poincare::TimeMap<MpC2Taylor>;
template class poincare::TimeMap<MpCnTaylor>;

template class poincare::TimeMap< dynsys::FadOdeSolver<dynsys::LorenzFadMap<MpInterval,0>,dynsys::ILastTermsStepControl > >;

#endif
