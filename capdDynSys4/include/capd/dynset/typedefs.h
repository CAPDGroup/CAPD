
#include "QRPolicy.h"

//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file dynset/typedefs.h
///
/// @author Daniel Wilczak   @date 2013-01-09
//
/////////////////////////////////////////////////////////////////////////////

// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

namespace CAPD_USER_NAMESPACE{

typedef capd::DInterval (*v_form)(CAPD_USER_NAMESPACE::IVector &);
/// @addtogroup capdlib
/// @{
typedef capd::dynset::FactorReorganization< capd::dynset::InverseQRPolicy<> > C0Pped2Policies;
typedef capd::dynset::FactorReorganization< capd::dynset::FullQRWithPivoting<> > C0Rect2Policies;
typedef capd::dynset::FactorReorganization<> C0Intv2Policies;
typedef capd::dynset::InverseQRPolicy<> C0PpedPolicies;
typedef capd::dynset::FullQRWithPivoting<>  C0RectPolicies;
typedef capd::dynset::SelectiveQRWithPivoting<> SelectiveQRPolicy;

typedef capd::dynset::InverseQRPolicy<> C1PpedPolicies;
typedef capd::dynset::FullQRWithPivoting<>  C1RectPolicies;
typedef capd::dynset::QRReorganization<capd::dynset::InverseQRPolicy<> > C1Pped2Policies;
typedef capd::dynset::FactorReorganization<capd::dynset::FullQRWithPivoting<> > C1Rect2Policies;

typedef capd::dynset::QRReorganization<capd::dynset::InverseQRPolicy<> > C2Pped2Policies;
typedef capd::dynset::FactorReorganization<capd::dynset::FullQRWithPivoting<> > C2Rect2Policies;

typedef capd::dynset::C0Set<CAPD_USER_NAMESPACE::IMatrix> C0Set;
typedef capd::dynset::C0BallSet<CAPD_USER_NAMESPACE::IMatrix> C0BallSet;
typedef capd::dynset::C0FlowballSet<CAPD_USER_NAMESPACE::IMatrix> C0FlowballSet;
typedef capd::dynset::C0AffineSet<CAPD_USER_NAMESPACE::IMatrix,C0PpedPolicies> C0PpedSet;
typedef capd::dynset::C0AffineSet<CAPD_USER_NAMESPACE::IMatrix,C0RectPolicies> C0RectSet;
typedef capd::dynset::C0DoubletonSet<CAPD_USER_NAMESPACE::IMatrix,C0Intv2Policies> C0Intv2Set;
typedef capd::dynset::C0DoubletonSet<CAPD_USER_NAMESPACE::IMatrix,C0Pped2Policies> C0Pped2Set;
typedef capd::dynset::C0DoubletonSet<CAPD_USER_NAMESPACE::IMatrix,C0Rect2Policies> C0Rect2Set;
typedef capd::dynset::C0TripletonSet<CAPD_USER_NAMESPACE::IMatrix,C0Rect2Policies> C0TripletonSet;
typedef capd::dynset::C0HOSet<CAPD_USER_NAMESPACE::C0Rect2Set> C0HORect2Set;
typedef capd::dynset::C0HOSet<CAPD_USER_NAMESPACE::C0TripletonSet> C0HOTripletonSet;

typedef capd::dynset::C1Set<CAPD_USER_NAMESPACE::IMatrix> C1Set;
typedef capd::dynset::C1AffineSet<CAPD_USER_NAMESPACE::IMatrix,C1RectPolicies> C1RectSet;
typedef capd::dynset::C1AffineSet<CAPD_USER_NAMESPACE::IMatrix,C1PpedPolicies> C1PpedSet;
typedef capd::dynset::C1DoubletonSet<CAPD_USER_NAMESPACE::IMatrix,C1Rect2Policies> C1Rect2Set;
typedef capd::dynset::C1DoubletonSet<CAPD_USER_NAMESPACE::IMatrix,C1Pped2Policies> C1Pped2Set;
typedef capd::dynset::C11Rect2Set<CAPD_USER_NAMESPACE::IMatrix> C11Rect2Set;
typedef capd::dynset::C1HOSet<CAPD_USER_NAMESPACE::C1Rect2Set> C1HORect2Set;
typedef capd::dynset::C1HOSet<CAPD_USER_NAMESPACE::C1Pped2Set> C1HOPped2Set;

typedef capd::dynset::C2Set<CAPD_USER_NAMESPACE::IMatrix> C2Set;
typedef capd::dynset::C2DoubletonSet<CAPD_USER_NAMESPACE::IMatrix,C2Rect2Policies> C2Rect2Set;
typedef capd::dynset::C2DoubletonSet<CAPD_USER_NAMESPACE::IMatrix,C2Pped2Policies> C2Pped2Set;

typedef capd::dynset::CnSet<CAPD_USER_NAMESPACE::IMatrix> CnSet;
typedef capd::dynset::CnRect2Set<CAPD_USER_NAMESPACE::IMatrix,C2Rect2Policies> CnRect2Set;
typedef capd::dynset::CnDoubletonSet<CAPD_USER_NAMESPACE::IMatrix,C2Rect2Policies> CnMultiMatrixRect2Set;

/// @}
} // end of namespace
