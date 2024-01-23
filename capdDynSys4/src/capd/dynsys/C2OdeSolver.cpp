/////////////////////////////////////////////////////////////////////////////
/// @file C2OdeSolver.cpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details. 

#include "capd/vectalg/lib.h"
#include "capd/map/lib.h"
#include "capd/dynsys/OdeSolver.hpp"
#include "capd/dynsys/BasicC2OdeSolver.hpp"
#include "capd/dynsys/C2OdeSolver.hpp"

template class capd::dynsys::C2OdeSolver<capd::IMap>;

template class capd::dynsys::BasicC2OdeSolver<capd::IMap,capd::dynsys::IEncFoundStepControl>;
template class capd::dynsys::BasicC2OdeSolver<capd::DMap>;
template class capd::dynsys::BasicC2OdeSolver<capd::LDMap>;

template class capd::dynsys::BasicOdeSolver<capd::DMap,capd::dynsys::DLastTermsStepControl,capd::diffAlgebra::C2Curve< capd::diffAlgebra::BasicC2Curve<capd::DMatrix> > >;
template class capd::dynsys::BasicOdeSolver<capd::LDMap,capd::dynsys::DLastTermsStepControl,capd::diffAlgebra::C2Curve< capd::diffAlgebra::BasicC2Curve<capd::LDMatrix> > >;
template class capd::dynsys::BasicOdeSolver<capd::IMap,capd::dynsys::ILastTermsStepControl,capd::diffAlgebra::C2Curve< capd::diffAlgebra::BasicC2Curve<capd::IMatrix> > >;
template class capd::dynsys::BasicOdeSolver<capd::IMap,capd::dynsys::IEncFoundStepControl,capd::diffAlgebra::C2Curve< capd::diffAlgebra::BasicC2Curve<capd::IMatrix> > >;
