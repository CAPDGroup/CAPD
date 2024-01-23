//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file TypeTraits.h
///
/// @author Tomasz Kapela   @date 2010-03-08
//
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) Tomasz Kapela 2010
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_CAPD_TYPETRAITS_H_
#define _CAPD_CAPD_TYPETRAITS_H_

//#include "capd/interval/Interval.h"
//#include "capd/vectalg/Z2.h"
#include <limits>
#include <algorithm>
#include <cmath>

namespace capd {

/**
 *  Defines type traits such as their values zero, one etc.
 *
 *  CAPD types should define specialization in their header files.
 *
 *  Known specialization are in
 *    - capd/interval/Interval.h
 *    - capd/filib/Interval.h
 *    - capd/interval/IComplex.h
 *    - capd/multiPrec/MpReal.h
 *    - capd/rings/Z2.h
 *    - capd/rings/Zp.h
 *
 */
template <typename T>
class TypeTraits {
public:
  typedef T Real;

  /// returns object set to zero
  static constexpr inline T zero() noexcept{
	return static_cast<T>(0.0);
  }

  /// returns object set to one
  static constexpr inline T one() noexcept{
	return static_cast<T>(1.0);
  }
  static constexpr T max(T a, T b) noexcept;
  static constexpr T min(T a, T b) noexcept;
/* Expected interface
  /// number of decimal digits
  static inline int numberOfDigits() noexcept;

  /// Machine epsilon (the difference between 1 and the least value greater than 1 that is representable).
  static inline T epsilon() noexcept;

  template <typename S>
  static inline T convert(const S & obj){
    return static_cast<T>(obj);
  }
  /// this flag is true for all interval types
  static const bool isInterval = false;

  static constexpr T max(T a, T b) noexcept;
  static constexpr T min(T a, T b) noexcept;
  static constexpr T abs(T a) noexcept;
  static bool isSingular(T a);
  static constexpr bool isInf(T a) noexcept;
  static constexpr bool isNaN(T a) noexcept;
   */

};


/// for given type returns object that represents zero
template <typename T>
inline T zero(){
  return TypeTraits<T>::zero();
}

/// for given type returns object that represents one (identity)
template <typename T>
inline T one(){
  return TypeTraits<T>::one();
}


template <typename T>
class IntegralTypeTraits {
 public:
  typedef T Real;

  /// returns object set to zero
  static inline T zero() {
	return static_cast<T>(0);
  }
  /// returns object set to one
  static inline T one() {
	return static_cast<T>(1);
  }
  /// number of decimal digits
  static inline int numberOfDigits() throw() {
	return std::numeric_limits<T>::digits10;
  }
  /// Machine epsilon (the difference between 1 and the least value greater than 1 that is representable).
  static inline T epsilon() throw() {
	return std::numeric_limits<T>::epsilon();
  }

  template<typename S>
  static inline T convert(const S &obj) {
	return static_cast<T>(obj);
  }
  /// this flag is true for all interval types
  static const bool isInterval = false;

  static constexpr T max(T a, T b) noexcept{
	return std::max(a, b);
  }
  static constexpr T min(T a, T b) noexcept{
	return std::min(a, b);
  }
  static constexpr T abs(T a) noexcept{
	return std::abs(a);
  }

  static bool isSingular(T a){
	return a == TypeTraits<T>::zero();
  }
};

template <typename T>
struct FloatingTypeTraits : public IntegralTypeTraits<T>{

  typedef T Real;

  /// returns object set to zero
  static inline T zero() {
	return static_cast<T>(0.0);
  }
  /// returns object set to one
  static inline T one() {
	return static_cast<T>(1.0);
  }
  /// number of decimal digits
  static inline int numberOfDigits() throw() {
	return std::numeric_limits<T>::digits10;
  }
  /// Machine epsilon (the difference between 1 and the least value greater than 1 that is representable).
  static inline T epsilon() throw() {
	return std::numeric_limits<T>::epsilon();
  }

  using IntegralTypeTraits<T>::convert;
  using IntegralTypeTraits<T>::isInterval;

  static constexpr bool isInf(T a) noexcept{
	return std::isinf(a);
  }
  static constexpr bool isNaN(T a) noexcept{
	return std::isnan(a);
  }
};

/**
 * Traits of type int
 */
template<>
struct TypeTraits<int> : public IntegralTypeTraits<int>{
  using Real = int;
  using T = int;
  static constexpr T max(T a, T b) noexcept{
	return std::max(a, b);
  }
  static constexpr T min(T a, T b) noexcept{
	return std::min(a, b);
  }
};

/**
 * Traits of type short
 */
template<>
struct TypeTraits<short> : public IntegralTypeTraits<short>{
//  using IntegralTypeTraits<short>::Real;
//  using IntegralTypeTraits<short>::convert;
};

/**
 * Traits of type long
 */
template<>
struct TypeTraits<long> : public IntegralTypeTraits<long>{
//  using IntegralTypeTraits<long>::Real;
//  using IntegralTypeTraits<long>::convert;
};

/**
 * Traits of type long long
 */
template<>
struct TypeTraits<long long> : public IntegralTypeTraits<long long>{
//  using IntegralTypeTraits<long long>::Real;
//  using IntegralTypeTraits<long long>::convert;
};

/**
 * Traits of type double
 */
template<>
struct TypeTraits<double> : public FloatingTypeTraits<double>{
  //using FloatingTypeTraits<double>::Real;
//  using FloatingTypeTraits<double>::convert;
};

/**
 * Traits of type float
 */
template<>
struct TypeTraits<float> : public FloatingTypeTraits<float>{
//  using FloatingTypeTraits<float>::Real;
//  using FloatingTypeTraits<float>::convert;
};

/**
 * Traits of type long double
 */
template<>
struct TypeTraits<long double> : public FloatingTypeTraits<long double>{
//  using FloatingTypeTraits<long double>::Real;
//  using FloatingTypeTraits<long double>::convert;
};


template <typename T>
class TypeTraits<T*> {
public:

  /// returns object set to zero
  static inline T* zero(){
    return static_cast<T*>(0);
  }
};

} // end of namespace capd

#endif // _CAPD_CAPD_TYPETRAITS_H_
