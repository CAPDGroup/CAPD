/////////////////////////////////////////////////////////////////////////////
/// @file CnOdeSolver.cpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details. 

#include "capd/map/lib.h"

#include "capd/dynsys/BasicCnOdeSolver.hpp"
#include "capd/dynsys/CnOdeSolver.hpp"
#include "capd/diffAlgebra/Jet.hpp"
#include "capd/diffAlgebra/CnTimeJet.h"
#include "capd/diffAlgebra/CnCurve.hpp"

template class capd::dynsys::BasicCnOdeSolver<capd::DMap>;
template class capd::dynsys::BasicCnOdeSolver<capd::LDMap>;
template class capd::dynsys::BasicCnOdeSolver<capd::IMap,capd::dynsys::ILastTermsStepControl>;
template class capd::dynsys::CnOdeSolver<capd::IMap,capd::dynsys::ILastTermsStepControl>;

template class capd::dynsys::BasicCnOdeSolver<capd::IMap,capd::dynsys::IEncFoundStepControl>;
template class capd::dynsys::CnOdeSolver<capd::IMap,capd::dynsys::IEncFoundStepControl>;
