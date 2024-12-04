//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file diffAlgebra/lib.h
///
/// @author Daniel Wilczak
//
/////////////////////////////////////////////////////////////////////////////


// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DIFFALGEBRA_LIB_H_
#define _CAPD_DIFFALGEBRA_LIB_H_

#include "capd/vectalg/lib.h"
#include "capd/diffAlgebra/C0TimeJet.h"
#include "capd/diffAlgebra/C1TimeJet.h"
#include "capd/diffAlgebra/C2TimeJet.h"
#include "capd/diffAlgebra/CnTimeJet.h"

namespace capd{

typedef capd::diffAlgebra::C0TimeJet<capd::IVector> IC0TimeJet;
typedef capd::diffAlgebra::C0TimeJet<capd::DVector> DC0TimeJet;
typedef capd::diffAlgebra::C0TimeJet<capd::LDVector> LDC0TimeJet;

typedef capd::diffAlgebra::C1TimeJet<capd::IMatrix> IC1TimeJet;
typedef capd::diffAlgebra::C1TimeJet<capd::DMatrix> DC1TimeJet;
typedef capd::diffAlgebra::C1TimeJet<capd::LDMatrix> LDC1TimeJet;

typedef capd::diffAlgebra::C2TimeJet<capd::IMatrix> IC2TimeJet;
typedef capd::diffAlgebra::C2TimeJet<capd::DMatrix> DC2TimeJet;
typedef capd::diffAlgebra::C2TimeJet<capd::LDMatrix> LDC2TimeJet;

typedef capd::diffAlgebra::CnTimeJet<capd::IMatrix, 0> ICnTimeJet;
typedef capd::diffAlgebra::CnTimeJet<capd::DMatrix, 0> DCnTimeJet;
typedef capd::diffAlgebra::CnTimeJet<capd::LDMatrix, 0> LDCnTimeJet;

typedef capd::diffAlgebra::Hessian<double,0,0> DHessian;
typedef capd::diffAlgebra::Hessian<long double,0,0> LDHessian;
typedef capd::diffAlgebra::Hessian<interval,0,0> IHessian;

typedef capd::diffAlgebra::Jet<capd::DMatrix,0> DJet;
typedef capd::diffAlgebra::Jet<capd::LDMatrix,0> LDJet;
typedef capd::diffAlgebra::Jet<capd::IMatrix,0> IJet;

} // end of namespace capd

#endif // _CAPD_DIFFALGEBRA_LIB_H_
