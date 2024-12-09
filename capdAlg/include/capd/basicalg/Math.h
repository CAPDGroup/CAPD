//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file Math.h
///
/// @author Daniel Wilczak   @date 2024-11-22
//
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) CAPD
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_CAPD_MATH_H_
#define _CAPD_CAPD_MATH_H_

#include <cmath>
#include "capd/basicalg/doubleFun.h"
namespace capd {

/**
 *  Defines requested interface for all numeric types. 
 *
 *  Specializations will be proivided for double, float, long double, Complex, MpFloat and related interval types
 *
 */
template <typename T>
class Math {};

template<class T>
class PrimitiveFloatMath{
public:
  typedef T Real;
  static constexpr inline Real _sqr(Real x) noexcept{	return x*x;  }
  static constexpr inline Real _log(Real x) noexcept{	return std::log(x);  }
  static constexpr inline Real _pow(Real x, int c) noexcept{	return std::pow(x,c);  }  
  static constexpr inline Real _sqrt(Real x) noexcept{	return std::sqrt(x);  }  
  static constexpr inline Real _exp(Real x) noexcept{	return std::exp(x);  }  
  static constexpr inline Real _sin(Real x) noexcept{	return std::sin(x);  }  
  static constexpr inline Real _cos(Real x) noexcept{	return std::cos(x);  }  
  static constexpr inline Real _atan(Real x) noexcept{	return std::atan(x);  }  
  static constexpr inline Real _asin(Real x) noexcept{	return std::asin(x);  }  
  static constexpr inline Real _acos(Real x) noexcept{	return std::acos(x);  }  
};

template<>
class Math<float> : public PrimitiveFloatMath<float>{};

template<>
class Math<double> : public PrimitiveFloatMath<double>{};

template<>
class Math<long double> : public PrimitiveFloatMath<long double>{};

} // end of namespace capd

#endif // _CAPD_CAPD_MATH_H_
