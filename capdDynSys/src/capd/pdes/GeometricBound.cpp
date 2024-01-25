/// @addtogroup pdes
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file GeometricBound.cpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2008-2016 by the CAPD Group.
//
// Distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.


#include "capd/pdes/GeometricBound.hpp"
#include "capd/vectalg/Vector.hpp"

namespace capd {
namespace pdes {

template class capd::pdes::GeometricBound<capd::interval>;

template GeometricBound<interval> operator+<interval>(const GeometricBound<interval>& x, const GeometricBound<interval>& y);
template GeometricBound<interval> operator-<interval>(const GeometricBound<interval>& x, const GeometricBound<interval>& y);
template GeometricBound<interval> operator*<interval>(const interval& s, const GeometricBound<interval>& x);
template GeometricBound<interval> operator*<interval>(const GeometricBound<interval>& x, const interval& s);
template GeometricBound<interval> operator*<interval>(const capd::vectalg::Matrix<interval,0,0>& A, const GeometricBound<interval>& v);
template std::ostream& operator<< <interval>(std::ostream& out, const GeometricBound<interval>& x);
template void split<interval>(const GeometricBound<interval>& X, GeometricBound<interval>& x, GeometricBound<interval>& dx);
template GeometricBound<interval> intersection<interval>(const GeometricBound<interval>& x, const GeometricBound<interval>& y);

} // end of namespace pdes
} // end of namespace capd

/// @}


