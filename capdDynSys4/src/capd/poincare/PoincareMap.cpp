
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


#include "capd/intervals/lib.h"
#include "capd/vectalg/lib.h"
#include "capd/map/lib.h"
#include "capd/dynsys/lib.h"
#include "capd/poincare/PoincareMap.hpp"
#include "capd/dynsys/FadMap.h"
#include "capd/dynsys/BasicFadOdeSolver.hpp"
#include "capd/dynsys/FadOdeSolver.hpp"

using namespace capd;

template class poincare::BasicPoincareMap<DOdeSolver>;
template class poincare::BasicPoincareMap<DC2OdeSolver>;
template class poincare::BasicPoincareMap<DCnOdeSolver>;
template class poincare::BasicPoincareMap< dynsys::BasicFadOdeSolver<dynsys::LorenzFadMap<double,0>,dynsys::DLastTermsStepControl > >;

template class poincare::BasicPoincareMap<LDOdeSolver>;
template class poincare::BasicPoincareMap<LDC2OdeSolver>;
template class poincare::BasicPoincareMap<LDCnOdeSolver>;
template class poincare::BasicPoincareMap< dynsys::BasicFadOdeSolver<dynsys::LorenzFadMap<long double,0>,dynsys::DLastTermsStepControl > >;

template class poincare::PoincareMap<IOdeSolver>;
template class poincare::PoincareMap<IC2OdeSolver>;
template class poincare::PoincareMap<ICnOdeSolver>;

template class poincare::PoincareMap< dynsys::FadOdeSolver<dynsys::LorenzFadMap<DInterval,0>,dynsys::ILastTermsStepControl > >;
