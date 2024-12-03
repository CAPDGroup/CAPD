//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file MpReal_Fun.hpp
///
/// @author Tomasz Kapela   @date 2010-03-15
//
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) Tomasz Kapela 2010
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_MULTIPREC_MPREAL_FUN_HPP_
#define _CAPD_MULTIPREC_MPREAL_FUN_HPP_

namespace capd {
namespace multiPrec {
//----------------------------------------------------------------------
//
//   Functions
//
//----------------------------------------------------------------------

inline MpReal MpReal::pow(const unsigned long int e, MpReal::RoundingMode rnd) const {
  MpReal res;
  mpfr_pow_ui(res.mpfr_rep, mpfr_rep, e, rnd);
  return res;
}

inline MpReal MpReal::pow(const long int e, MpReal::RoundingMode rnd) const {
  MpReal res;
  mpfr_pow_si(res.mpfr_rep, mpfr_rep, e, rnd);
  return res;
}

inline MpReal MpReal::pow(const MpReal& e, MpReal::RoundingMode rnd) const {
  MpReal res;
  mpfr_pow(res.mpfr_rep, mpfr_rep, e.mpfr_rep, rnd);
  return res;
}


//inline MpReal factorial(const unsigned long int n, MpReal::RoundingMode rnd) {
//  MpReal res;
//  mpfr_fac_ui(res.mpfr_rep, n, rnd);
//  return res;
//}

//--------------------------------------------------------------
//
// Mathematical and miscellaneous functions
//
//--------------------------------------------------------------
// NaN of infinity handled by MPFR
// => no test on the value (sign...) of the operand
//void MpReal::random(PrecisionType prec) // member to avoid conflict with GMP random
//{
//  if(nbref->decr() <= 0) // the memory can be reused
//  {
//    nbref->refvalue() = 1;
//    mpfr_set_prec(mpfr_rep, prec);
//  } else // the previous value must be preserved
//  {
//    mpfr_init2(mpfr_rep, prec);
//    nbref = new RefCounter(1);
//    inexact = new InexactFlag();
//  }
//  mpfr_random(mpfr_rep);
//  inexact->refvalue() = EXACT_FLAG;
//}


inline MpReal round(const MpReal& r) {
  MpReal res;
  mpfr_round(res.mpfr_rep, r.mpfr_rep);
  return res;
}

inline MpReal floor(const MpReal& r) {
  MpReal res;
  mpfr_floor(res.mpfr_rep, r.mpfr_rep);
  return res;
}

inline MpReal trunc(const MpReal& r) {
  MpReal res;
  mpfr_trunc(res.mpfr_rep, r.mpfr_rep);
  return res;
}

inline MpReal ceil(const MpReal& r) {
  MpReal res;
  mpfr_ceil(res.mpfr_rep, r.mpfr_rep);
  return res;
}

inline MpReal frac(const MpReal& r) {
  MpReal res;
  mpfr_frac(res.mpfr_rep, r.mpfr_rep, MpReal::getDefaultRndMode());
  return res;
}

inline const MpReal & nonnegativePart(const MpReal & x) {
  if(x >= 0)
    return x;
  else
    throw std::runtime_error(" MpReal : nonnegative part is empty");
}

inline MpReal right(const MpReal& x) {
  return x;
}

inline MpReal left(const MpReal& x) {
  return x;
}

inline MpReal mid(const MpReal& x) {
  return x;
}
// end new functions
//
//
}
} // end of namespace capd::multiPrec

#endif // _CAPD_MULTIPREC_MPREAL_FUN_HPP_
