/////////////////////////////////////////////////////////////////////////////
/// @file Curve.cpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
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

using namespace capd;

template class diffAlgebra::BasicCurve<DMatrix>;
template class diffAlgebra::BasicCurve<LDMatrix>;
template class diffAlgebra::BasicCurve<IMatrix>;

template class diffAlgebra::FadCurve<DMatrix>;
template class diffAlgebra::FadCurve<LDMatrix>;
template class diffAlgebra::FadCurve<IMatrix>;

template class diffAlgebra::BasicC2Curve<DMatrix>;
template class diffAlgebra::BasicC2Curve<LDMatrix>;
template class diffAlgebra::BasicC2Curve<IMatrix>;

template class diffAlgebra::BasicCnCurve<DMatrix>;
template class diffAlgebra::BasicCnCurve<LDMatrix>;
template class diffAlgebra::BasicCnCurve<IMatrix>;



template class diffAlgebra::Curve< diffAlgebra::BasicCurve<DMatrix> >;
template class diffAlgebra::Curve< diffAlgebra::BasicCurve<LDMatrix> >;
template class diffAlgebra::Curve< diffAlgebra::BasicCurve<IMatrix> >;

template class diffAlgebra::Curve< diffAlgebra::FadCurve<DMatrix> >;
template class diffAlgebra::Curve< diffAlgebra::FadCurve<LDMatrix> >;
template class diffAlgebra::Curve< diffAlgebra::FadCurve<IMatrix> >;

template class diffAlgebra::Curve< diffAlgebra::BasicC2Curve<DMatrix> >;
template class diffAlgebra::Curve< diffAlgebra::BasicC2Curve<LDMatrix> >;
template class diffAlgebra::Curve< diffAlgebra::BasicC2Curve<IMatrix> >;

template class diffAlgebra::C2Curve< diffAlgebra::BasicC2Curve<DMatrix> >;
template class diffAlgebra::C2Curve< diffAlgebra::BasicC2Curve<LDMatrix> >;
template class diffAlgebra::C2Curve< diffAlgebra::BasicC2Curve<IMatrix> >;

template class diffAlgebra::Curve< diffAlgebra::BasicCnCurve<DMatrix> >;
template class diffAlgebra::Curve< diffAlgebra::BasicCnCurve<LDMatrix> >;
template class diffAlgebra::Curve< diffAlgebra::BasicCnCurve<IMatrix> >;

template class diffAlgebra::CnCurve< diffAlgebra::BasicCnCurve<DMatrix> >;
template class diffAlgebra::CnCurve< diffAlgebra::BasicCnCurve<LDMatrix> >;
template class diffAlgebra::CnCurve< diffAlgebra::BasicCnCurve<IMatrix> >;

