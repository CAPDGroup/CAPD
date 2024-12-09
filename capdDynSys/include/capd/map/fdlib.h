//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file map/fdlib.h
///
/// @author Daniel Wilczak   @date 2013-01-09
//
/////////////////////////////////////////////////////////////////////////////

// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#include "capd/map/Function.hpp"
#include "capd/map/Map.hpp"

namespace CAPD_USER_NAMESPACE{

typedef capd::map::Function<CAPD_USER_NAMESPACE::IVector> IFunction;
typedef capd::map::Map<CAPD_USER_NAMESPACE::IMatrix> IMap;

typedef capd::map::Function<CAPD_USER_NAMESPACE::DVector> DFunction;
typedef capd::map::Map<CAPD_USER_NAMESPACE::DMatrix> DMap;

typedef capd::map::Function<CAPD_USER_NAMESPACE::LDVector> LDFunction;
typedef capd::map::Map<CAPD_USER_NAMESPACE::LDMatrix> LDMap;

}
