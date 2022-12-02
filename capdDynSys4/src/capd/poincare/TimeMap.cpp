
/////////////////////////////////////////////////////////////////////////////
/// @file TimeMap.cpp
///
/// @author Tomasz Kapela
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2013 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#include "capd/map/lib.h"
#include "capd/dynsys/lib.h"
#include "capd/poincare/TimeMap.hpp"
#include "capd/dynsys/FadMap.h"
#include "capd/dynsys/BasicFadOdeSolver.hpp"
#include "capd/dynsys/FadOdeSolver.hpp"

using namespace capd;

template class poincare::TimeMap<IOdeSolver>;
template class poincare::TimeMap<IC2OdeSolver>;
template class poincare::TimeMap<ICnOdeSolver>;

template class poincare::TimeMap<DOdeSolver>;
template class poincare::TimeMap<DC2OdeSolver>;
template class poincare::TimeMap<DCnOdeSolver>;

template class poincare::TimeMap<LDOdeSolver>;
template class poincare::TimeMap<LDC2OdeSolver>;
template class poincare::TimeMap<LDCnOdeSolver>;

template class poincare::TimeMap< dynsys::BasicFadOdeSolver<dynsys::LorenzFadMap<long double,0>,dynsys::DLastTermsStepControl > >;
template class poincare::TimeMap< dynsys::FadOdeSolver<dynsys::LorenzFadMap<DInterval,0>,dynsys::ILastTermsStepControl > >;
