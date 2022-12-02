/////////////////////////////////////////////////////////////////////////////
/// @file SolutionCurve.cpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2013 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.


#include "capd/vectalg/lib.h"
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
template class BaseSolutionCurve< Curve < BasicCurve <IMatrix> > >;
template class BaseSolutionCurve< Curve < FadCurve <IMatrix> > >;
template class BaseSolutionCurve< Curve < BasicC2Curve <IMatrix> > >;
template class BaseSolutionCurve< C2Curve < BasicC2Curve <IMatrix> > >;

template class BaseSolutionCurve< Curve < BasicCurve <DMatrix> > >;
template class BaseSolutionCurve< Curve < FadCurve <DMatrix> > >;
template class BaseSolutionCurve< Curve < BasicC2Curve <DMatrix> > >;
template class BaseSolutionCurve< C2Curve < BasicC2Curve <DMatrix> > >;

template class BaseSolutionCurve< Curve < BasicCurve <LDMatrix> > >;
template class BaseSolutionCurve< Curve < FadCurve <LDMatrix> > >;
template class BaseSolutionCurve< Curve < BasicC2Curve <LDMatrix> > >;
template class BaseSolutionCurve< C2Curve < BasicC2Curve <LDMatrix> > >;

template class SolutionCurve< Curve < BasicCurve <IMatrix> > >;
template class SolutionCurve< Curve < FadCurve <IMatrix> > >;
template class SolutionCurve< Curve < BasicC2Curve <IMatrix> > >;
template class SolutionCurve< C2Curve < BasicC2Curve <IMatrix> > >;

template class SolutionCurve< Curve < BasicCurve <DMatrix> > >;
template class SolutionCurve< Curve < FadCurve <DMatrix> > >;
template class SolutionCurve< Curve < BasicC2Curve <DMatrix> > >;
template class SolutionCurve< C2Curve < BasicC2Curve <DMatrix> > >;

template class SolutionCurve< Curve < BasicCurve <LDMatrix> > >;
template class SolutionCurve< Curve < FadCurve <LDMatrix> > >;
template class SolutionCurve< Curve < BasicC2Curve <LDMatrix> > >;
template class SolutionCurve< C2Curve < BasicC2Curve <LDMatrix> > >;

}}
