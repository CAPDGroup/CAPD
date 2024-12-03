

/////////////////////////////////////////////////////////////////////////////
/// @file minmax.h
///
/// @author The CAPD Group
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.edu.pl/ for details.

/* min, max and abs definitions */

#ifndef _CAPD_BASICALG_MINMAX_H_
#define _CAPD_BASICALG_MINMAX_H_

#include "capd/basicalg/TypeTraits.h"

#undef max
#undef min

namespace capd {
/// @addtogroup capd
/// @{

//
// The following lines was prone to errors
//
//template<typename T>
//inline T min(const T& x, const T& y) {
//  return (x<y ? x : y);
//}
//
//template<typename T>
//inline T max(const T& x, const T& y) {
//  return (x<y ? y : x);
//}


  template<typename T>
  inline T min(const T &x, const T &y) {
	return ::capd::TypeTraits<T>::min(x, y);
  }

  template<typename T>
  inline T max(const T &x, const T &y) {
	return ::capd::TypeTraits<T>::max(x, y);
  }

  template<typename T>
  inline T abs(const T &x) {
	return ::capd::TypeTraits<T>::abs(x);
  }
//
//  inline long double abs(long double x) {
//	return (x < 0.) ? -x : x;
//  }
//
//  inline double abs(double x) {
//	return (x < 0.) ? -x : x;
//  }
//
//  inline int abs(int x) {
//	return (x < 0.) ? -x : x;
//  }

/// @}
}
#endif // _CAPD_BASICALG_MINMAX_H_

