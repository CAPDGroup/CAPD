/////////////////////////////////////////////////////////////////////////////
/// @file mpcapd/dynsys/C2OdeSolver.cpp
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
#include "capd/dynsys/OdeSolver.hpp"
#include "capd/dynsys/BasicC2OdeSolver.hpp"
#include "capd/dynsys/C2OdeSolver.hpp"


#ifdef __HAVE_MPFR__

template class capd::dynsys::C2OdeSolver<capd::MpIMap>;

template class capd::dynsys::BasicC2OdeSolver<capd::MpIMap,capd::dynsys::IEncFoundStepControl>;
template class capd::dynsys::BasicC2OdeSolver<capd::MpIMap>;
template class capd::dynsys::BasicC2OdeSolver<capd::MpMap>;


template class capd::dynsys::BasicOdeSolver<capd::MpMap,capd::dynsys::DLastTermsStepControl,capd::diffAlgebra::C2Curve< capd::diffAlgebra::BasicC2Curve<capd::MpMatrix> > >;
template class capd::dynsys::BasicOdeSolver<capd::MpIMap,capd::dynsys::ILastTermsStepControl,capd::diffAlgebra::C2Curve< capd::diffAlgebra::BasicC2Curve<capd::MpIMatrix> > >;
template class capd::dynsys::BasicOdeSolver<capd::MpIMap,capd::dynsys::IEncFoundStepControl,capd::diffAlgebra::C2Curve< capd::diffAlgebra::BasicC2Curve<capd::MpIMatrix> > >;

#endif

