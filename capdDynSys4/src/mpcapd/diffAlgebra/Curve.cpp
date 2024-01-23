/////////////////////////////////////////////////////////////////////////////
/// @file mpcapd/diffAlgebra/Curve.cpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
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

using namespace capd;

#ifdef __HAVE_MPFR__

  template class diffAlgebra::BasicCurve<MpMatrix>;
  template class diffAlgebra::BasicCurve<MpIMatrix>;

  template class diffAlgebra::FadCurve<MpMatrix>;
  template class diffAlgebra::FadCurve<MpIMatrix>;

  template class diffAlgebra::BasicC2Curve<MpMatrix>;
  template class diffAlgebra::BasicC2Curve<MpIMatrix>;

  template class diffAlgebra::BasicCnCurve<MpMatrix>;
  template class diffAlgebra::BasicCnCurve<MpIMatrix>;


  template class diffAlgebra::Curve< diffAlgebra::BasicCurve<MpMatrix> >;
  template class diffAlgebra::Curve< diffAlgebra::BasicCurve<MpIMatrix> >;

  template class diffAlgebra::Curve< diffAlgebra::FadCurve<MpMatrix> >;
  template class diffAlgebra::Curve< diffAlgebra::FadCurve<MpIMatrix> >;

  template class diffAlgebra::Curve< diffAlgebra::BasicC2Curve<MpMatrix> >;
  template class diffAlgebra::Curve< diffAlgebra::BasicC2Curve<MpIMatrix> >;

  template class diffAlgebra::C2Curve< diffAlgebra::BasicC2Curve<MpMatrix> >;
  template class diffAlgebra::C2Curve< diffAlgebra::BasicC2Curve<MpIMatrix> >;

  template class diffAlgebra::Curve< diffAlgebra::BasicCnCurve<MpMatrix> >;
  template class diffAlgebra::Curve< diffAlgebra::BasicCnCurve<MpIMatrix> >;

  template class diffAlgebra::CnCurve< diffAlgebra::BasicCnCurve<MpMatrix> >;
  template class diffAlgebra::CnCurve< diffAlgebra::BasicCnCurve<MpIMatrix> >;

  #endif
