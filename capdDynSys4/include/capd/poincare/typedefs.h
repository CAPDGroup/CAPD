//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file poincare/typedefs.h
///
/// @author Daniel Wilczak   @date 2013-01-09
//
/////////////////////////////////////////////////////////////////////////////

// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

namespace CAPD_USER_NAMESPACE{

// classes for nonrigorous computations

typedef capd::poincare::AbstractSection<CAPD_USER_NAMESPACE::DMatrix> DAbstractSection;
typedef capd::poincare::AbstractSection<CAPD_USER_NAMESPACE::LDMatrix> LDAbstractSection;
typedef capd::poincare::AbstractSection<CAPD_USER_NAMESPACE::IMatrix> IAbstractSection;

typedef capd::poincare::AffineSection<CAPD_USER_NAMESPACE::DMatrix> DAffineSection;
typedef capd::poincare::AffineSection<CAPD_USER_NAMESPACE::LDMatrix> LDAffineSection;
typedef capd::poincare::AffineSection<CAPD_USER_NAMESPACE::IMatrix> IAffineSection;

typedef capd::poincare::CoordinateSection<CAPD_USER_NAMESPACE::DMatrix> DCoordinateSection;
typedef capd::poincare::CoordinateSection<CAPD_USER_NAMESPACE::LDMatrix> LDCoordinateSection;
typedef capd::poincare::CoordinateSection<CAPD_USER_NAMESPACE::IMatrix> ICoordinateSection;

typedef capd::poincare::NonlinearSection<CAPD_USER_NAMESPACE::DMatrix> DNonlinearSection;
typedef capd::poincare::NonlinearSection<CAPD_USER_NAMESPACE::LDMatrix> LDNonlinearSection;
typedef capd::poincare::NonlinearSection<CAPD_USER_NAMESPACE::IMatrix> INonlinearSection;

typedef capd::poincare::BasicPoincareMap<CAPD_USER_NAMESPACE::DOdeSolver> DPoincareMap;
typedef capd::poincare::BasicPoincareMap<CAPD_USER_NAMESPACE::DC2OdeSolver> DC2PoincareMap;
typedef capd::poincare::BasicPoincareMap<CAPD_USER_NAMESPACE::DCnOdeSolver> DCnPoincareMap;

typedef capd::poincare::BasicPoincareMap<CAPD_USER_NAMESPACE::LDOdeSolver> LDPoincareMap;
typedef capd::poincare::BasicPoincareMap<CAPD_USER_NAMESPACE::LDC2OdeSolver> LDC2PoincareMap;
typedef capd::poincare::BasicPoincareMap<CAPD_USER_NAMESPACE::LDCnOdeSolver> LDCnPoincareMap;

typedef capd::poincare::PoincareMap<CAPD_USER_NAMESPACE::IOdeSolver> IPoincareMap;
typedef capd::poincare::PoincareMap<CAPD_USER_NAMESPACE::IC2OdeSolver> IC2PoincareMap;
typedef capd::poincare::PoincareMap<CAPD_USER_NAMESPACE::ICnOdeSolver> ICnPoincareMap;


typedef capd::poincare::TimeMap<CAPD_USER_NAMESPACE::DOdeSolver> DTimeMap;
typedef capd::poincare::TimeMap<CAPD_USER_NAMESPACE::DC2OdeSolver> DC2TimeMap;
typedef capd::poincare::TimeMap<CAPD_USER_NAMESPACE::DCnOdeSolver> DCnTimeMap;

typedef capd::poincare::TimeMap<CAPD_USER_NAMESPACE::LDOdeSolver> LDTimeMap;
typedef capd::poincare::TimeMap<CAPD_USER_NAMESPACE::LDC2OdeSolver> LDC2TimeMap;
typedef capd::poincare::TimeMap<CAPD_USER_NAMESPACE::LDCnOdeSolver> LDCnTimeMap;

typedef capd::poincare::TimeMap<CAPD_USER_NAMESPACE::IOdeSolver> ITimeMap;
typedef capd::poincare::TimeMap<CAPD_USER_NAMESPACE::IC2OdeSolver> IC2TimeMap;
typedef capd::poincare::TimeMap<CAPD_USER_NAMESPACE::ICnOdeSolver> ICnTimeMap;

} // end of namespace capd
