/////////////////////////////////////////////////////////////////////////////
/// @file mpcapd/diffAlgebra/SolutionCurve.cpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2013 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.


#include "capd/vectalg/mplib.h"
#include "capd/diffAlgebra/BasicCurve.hpp"
#include "capd/diffAlgebra/FadCurve.hpp"
#include "capd/diffAlgebra/Curve.hpp"

#include "capd/diffAlgebra/BasicC2Curve.hpp"
#include "capd/diffAlgebra/C2Curve.hpp"

#include "capd/diffAlgebra/BasicCnCurve.hpp"
#include "capd/diffAlgebra/CnCurve.hpp"

#include "capd/diffAlgebra/SolutionCurve.h"

namespace capd{
namespace diffAlgebra{
/// @addtogroup diffAlgebra 
/// @{
template class BaseSolutionCurve< Curve < BasicCurve <MpIMatrix> > >;
template class BaseSolutionCurve< Curve < FadCurve <MpIMatrix> > >;
template class BaseSolutionCurve< Curve < BasicC2Curve <MpIMatrix> > >;
template class BaseSolutionCurve< C2Curve < BasicC2Curve <MpIMatrix> > >;

template class BaseSolutionCurve< Curve < BasicCurve <MpMatrix> > >;
template class BaseSolutionCurve< Curve < FadCurve <MpMatrix> > >;
template class BaseSolutionCurve< Curve < BasicC2Curve <MpMatrix> > >;
template class BaseSolutionCurve< C2Curve < BasicC2Curve <MpMatrix> > >;

template class SolutionCurve< Curve < BasicCurve <MpIMatrix> > >;
template class SolutionCurve< Curve < FadCurve <MpIMatrix> > >;
template class SolutionCurve< Curve < BasicC2Curve <MpIMatrix> > >;
template class SolutionCurve< C2Curve < BasicC2Curve <MpIMatrix> > >;

template class SolutionCurve< Curve < BasicCurve <MpMatrix> > >;
template class SolutionCurve< Curve < FadCurve <MpMatrix> > >;
template class SolutionCurve< Curve < BasicC2Curve <MpMatrix> > >;
template class SolutionCurve< C2Curve < BasicC2Curve <MpMatrix> > >;

}}
