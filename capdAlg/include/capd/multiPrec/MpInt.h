//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file MpInt.h
///
///  Defines C++ wrapper for multiple precision integers
///  from gmpxx library
///
/// @author Tomasz Kapela   @date 2010-03-08
//
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) Tomasz Kapela 2010
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

// Protects against compilations in systems without mpfr and gmp package
#ifdef __HAVE_MPFR__

#ifndef _CAPD_MULTIPREC_MPINT_H_
#define _CAPD_MULTIPREC_MPINT_H_

#include <cstddef>
#include <gmpxx.h>
#include <gmp.h>
#include <stdexcept>

namespace capd{
namespace multiPrec{

typedef mpz_class MpInt;


}} // end of namespace capd::multiPrec

inline capd::multiPrec::MpInt power(const capd::multiPrec::MpInt & base, int exp){
  capd::multiPrec::MpInt result;
  mpz_pow_ui(result.get_mpz_t(), base.get_mpz_t(), exp);
  return result;
}

inline capd::multiPrec::MpInt nonnegativePart(const capd::multiPrec::MpInt & x) {
  if(x>=0)
    return x;
  else
    throw "non negative part is empty!";
}

#include "capd/basicalg/minmax.h"
#include "capd/basicalg/doubleFun.h"
#include "capd/basicalg/TypeTraits.h"
namespace capd{

template<>
struct TypeTraits<capd::multiPrec::MpInt> {
  typedef capd::multiPrec::MpInt Real;

  /// returns object set to zero
  static inline Real zero(){
    return Real(0);
  }
  /// returns object set to one
  static inline Real one(){
    return Real(1);
  }
  /// number of decimal digits
  static inline int numberOfDigits(){
    return 0; // constant = log(2)/log(10);
  }
  /// Machine epsilon (the difference between 1 and the least value greater than 1 that is representable).
  static inline Real epsilon() {
      return 1;
  }
  template <typename S>
  static inline Real convert(const S & obj){
    return static_cast<Real>(obj);
  }
  /// this flag is true for all interval types
  static const bool isInterval = false;
};

template<>
inline capd::multiPrec::MpInt abs(const capd::multiPrec::MpInt & x) {
  return ::abs(x);
}

template<>
inline long convertToLong(const capd::multiPrec::MpInt& x)
{
  if (mpz_fits_slong_p(x.get_mpz_t()) == 0) {
    throw std::logic_error("MpInt value does not fit into long");
  }

  return mpz_get_si(x.get_mpz_t());
}

}
#endif // _CAPD_MULTIPREC_MPINT_H_

#endif   // __HAVE_MPFR__
