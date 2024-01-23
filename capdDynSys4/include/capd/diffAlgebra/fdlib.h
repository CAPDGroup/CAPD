//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file diffAlgebra/fdlib.h
///
/// @author Daniel Wilczak
//
/////////////////////////////////////////////////////////////////////////////


// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#include "capd/diffAlgebra/C0TimeJet.h"
#include "capd/diffAlgebra/C1TimeJet.h"
#include "capd/diffAlgebra/C2TimeJet.h"
#include "capd/diffAlgebra/CnTimeJet.h"

namespace CAPD_USER_NAMESPACE{

typedef capd::diffAlgebra::C0TimeJet<CAPD_USER_NAMESPACE::IVector> IC0TimeJet;
typedef capd::diffAlgebra::C0TimeJet<CAPD_USER_NAMESPACE::DVector> DC0TimeJet;
typedef capd::diffAlgebra::C0TimeJet<CAPD_USER_NAMESPACE::LDVector> LDC0TimeJet;

typedef capd::diffAlgebra::C1TimeJet<CAPD_USER_NAMESPACE::IMatrix> IC1TimeJet;
typedef capd::diffAlgebra::C1TimeJet<CAPD_USER_NAMESPACE::DMatrix> DC1TimeJet;
typedef capd::diffAlgebra::C1TimeJet<CAPD_USER_NAMESPACE::LDMatrix> LDC1TimeJet;

typedef capd::diffAlgebra::C2TimeJet<CAPD_USER_NAMESPACE::IMatrix> IC2TimeJet;
typedef capd::diffAlgebra::C2TimeJet<CAPD_USER_NAMESPACE::DMatrix> DC2TimeJet;
typedef capd::diffAlgebra::C2TimeJet<CAPD_USER_NAMESPACE::LDMatrix> LDC2TimeJet;

typedef capd::diffAlgebra::CnTimeJet<CAPD_USER_NAMESPACE::IMatrix, 0> ICnTimeJet;
typedef capd::diffAlgebra::CnTimeJet<CAPD_USER_NAMESPACE::DMatrix, 0> DCnTimeJet;
typedef capd::diffAlgebra::CnTimeJet<CAPD_USER_NAMESPACE::LDMatrix, 0> LDCnTimeJet;

typedef capd::diffAlgebra::Hessian<CAPD_USER_NAMESPACE::DVector::ScalarType,CAPD_USER_NAMESPACE::DVector::csDim,CAPD_USER_NAMESPACE::DMatrix::ColumnVectorType::csDim> DHessian;
typedef capd::diffAlgebra::Hessian<CAPD_USER_NAMESPACE::LDVector::ScalarType,CAPD_USER_NAMESPACE::LDVector::csDim,CAPD_USER_NAMESPACE::LDMatrix::ColumnVectorType::csDim> LDHessian;
typedef capd::diffAlgebra::Hessian<CAPD_USER_NAMESPACE::IVector::ScalarType,CAPD_USER_NAMESPACE::IVector::csDim,CAPD_USER_NAMESPACE::IMatrix::ColumnVectorType::csDim> IHessian;

typedef capd::diffAlgebra::Jet<CAPD_USER_NAMESPACE::DMatrix,0> DJet;
typedef capd::diffAlgebra::Jet<CAPD_USER_NAMESPACE::LDMatrix,0> LDJet;
typedef capd::diffAlgebra::Jet<CAPD_USER_NAMESPACE::IMatrix,0> IJet;

} // end of namespace CAPD_USER_NAMESPACE

