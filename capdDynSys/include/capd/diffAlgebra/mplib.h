//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file diffAlgebra/mplib.h
///
/// @author Daniel Wilczak
//
/////////////////////////////////////////////////////////////////////////////


// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DIFFALGEBRA_MPLIB_H_
#define _CAPD_DIFFALGEBRA_MPLIB_H_

#include "capd/vectalg/mplib.h"
#include "capd/diffAlgebra/C0TimeJet.h"
#include "capd/diffAlgebra/C1TimeJet.h"
#include "capd/diffAlgebra/C2TimeJet.h"
#include "capd/diffAlgebra/CnTimeJet.h"

#ifdef __HAVE_MPFR__

namespace capd{

typedef capd::diffAlgebra::C0TimeJet<capd::MpIVector> MpIC0TimeJet;
typedef capd::diffAlgebra::C0TimeJet<capd::MpVector> MpC0TimeJet;

typedef capd::diffAlgebra::C1TimeJet<capd::MpIMatrix> MpIC1TimeJet;
typedef capd::diffAlgebra::C1TimeJet<capd::MpMatrix> MpC1TimeJet;

typedef capd::diffAlgebra::C2TimeJet<capd::MpIMatrix> MpIC2TimeJet;
typedef capd::diffAlgebra::C2TimeJet<capd::MpMatrix> MpC2TimeJet;

typedef capd::diffAlgebra::CnTimeJet<capd::MpIMatrix, 0> MpICnTimeJet;
typedef capd::diffAlgebra::CnTimeJet<capd::MpMatrix, 0> MpCnTimeJet;

typedef capd::diffAlgebra::Hessian<capd::MpFloat,0,0> MpHessian;
typedef capd::diffAlgebra::Hessian<capd::MpInterval,0,0> MpIHessian;

typedef capd::diffAlgebra::Jet<capd::MpMatrix,0> MpJet;
typedef capd::diffAlgebra::Jet<capd::MpIMatrix,0> MpIJet;

} // end of namespace capd

#endif //__HAVE_MPFR__

#endif // _CAPD_DIFFALGEBRA_MPLIB_H_
