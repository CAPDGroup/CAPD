/////////////////////////////////////////////////////////////////////////////
/// @file minmax_interval.h
///
/// @author Mateusz Juda <mateusz.juda@gmail.com>
///
/// @date 2014-07-23
/////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library (capdAux),
// distributed under the terms of the GNU General Public License.
// Consult http://capd.ii.uj.edu.pl and  http://redhom.ii.edu.pl/ for details.
/////////////////////////////////////////////////////////////////////////////

#ifndef CAPD_FILE_CAPDAUX_AUXIL_MINMAX_INTERVAL_H
#define CAPD_FILE_CAPDAUX_AUXIL_MINMAX_INTERVAL_H

#include <capd/basicalg/minmax.h>

namespace capd{
namespace intervals{
/// @addtogroup intervals 
/// @{
template < typename T_Bound, typename T_Rnd >
class Interval;
}  //namespace intervals


} // namespace capd

//#ifdef  __USE_CXSC__

namespace capd{
namespace cxsc{
class Interval;
} // end of namespace cxsc


} // namespace capd
//#endif // __USE_FILIB__


#ifdef  __USE_FILIB__

#include <interval/interval.hpp>

namespace capd{
namespace filib{

template <typename T, ::filib::rounding_strategy R, ::filib::interval_mode M>
class Interval;

} // end of namespace filib

} // namespace capd
#endif // __USE_FILIB__

#endif // CAPD_FILE_CAPDAUX_AUXIL_MINMAX_INTERVAL_H
