//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file MpReal_FriendInternal.hpp
///
/// @author Tomasz Kapela   @date 2010-03-15
//
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) Tomasz Kapela 2010
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.



inline friend MpReal agm(const MpReal& r1, const MpReal& r2, MpReal::RoundingMode rnd= MpReal::getDefaultRndMode()) {
  MpReal res;
  mpfr_agm(res.mpfr_rep, r1.mpfr_rep, r2.mpfr_rep, rnd);
  return res;
}

inline friend MpReal sqr(const MpReal& r, MpReal::RoundingMode rnd= MpReal::getDefaultRndMode()) {
  MpReal res;
  // NaN handled by MPFR in case r < 0
  mpfr_sqr(res.mpfr_rep, r.mpfr_rep, rnd);
  return res;
}

inline friend MpReal sqrt(const MpReal& r, MpReal::RoundingMode rnd= MpReal::getDefaultRndMode()) {
  MpReal res;
  // NaN handled by MPFR in case r < 0
  mpfr_sqrt(res.mpfr_rep, r.mpfr_rep, rnd);
  return res;
}

inline friend MpReal exp(const MpReal& r, MpReal::RoundingMode rnd= MpReal::getDefaultRndMode()) {
  MpReal res;
  mpfr_exp(res.mpfr_rep, r.mpfr_rep, rnd);
  return res;
}

inline friend MpReal expm1(const MpReal& r, MpReal::RoundingMode rnd= MpReal::getDefaultRndMode()) {
  MpReal res;
  mpfr_expm1(res.mpfr_rep, r.mpfr_rep, rnd);
  return res;
}

inline friend MpReal exp2(const MpReal& r, MpReal::RoundingMode rnd= MpReal::getDefaultRndMode()) {
  MpReal res;
  mpfr_exp2(res.mpfr_rep, r.mpfr_rep, rnd);
  return res;
}

inline friend MpReal log(const MpReal& r, MpReal::RoundingMode rnd= MpReal::getDefaultRndMode()) {
  MpReal res;
  mpfr_log(res.mpfr_rep, r.mpfr_rep, rnd);
  return res;
}

inline friend MpReal log2(const MpReal& r, MpReal::RoundingMode rnd= MpReal::getDefaultRndMode()) {
  MpReal res;
  mpfr_log2(res.mpfr_rep, r.mpfr_rep, rnd);
  return res;
}

inline friend MpReal log10(const MpReal& r, MpReal::RoundingMode rnd= MpReal::getDefaultRndMode()) {
  MpReal res;
  mpfr_log10(res.mpfr_rep, r.mpfr_rep, rnd);
  return res;
}

inline friend MpReal log1p(const MpReal& r, MpReal::RoundingMode rnd= MpReal::getDefaultRndMode()) {
  MpReal res;
  mpfr_log1p(res.mpfr_rep, r.mpfr_rep, rnd);
  return res;
}

inline friend MpReal sin(const MpReal& r, MpReal::RoundingMode rnd= MpReal::getDefaultRndMode()) {
  MpReal res;
  mpfr_sin(res.mpfr_rep, r.mpfr_rep, rnd);
  return res;
}

inline friend MpReal cos(const MpReal& r, MpReal::RoundingMode rnd= MpReal::getDefaultRndMode()) {
  MpReal res;
  mpfr_cos(res.mpfr_rep, r.mpfr_rep, rnd);
  return res;
}

inline friend void sin_cos(MpReal& res_sin, MpReal& res_cos, const MpReal& r,
    MpReal::RoundingMode rnd= MpReal::getDefaultRndMode()) {
  mpfr_sin_cos(res_sin.mpfr_rep, res_cos.mpfr_rep, r.mpfr_rep, rnd);
}

inline friend MpReal tan(const MpReal& r, MpReal::RoundingMode rnd= MpReal::getDefaultRndMode()) {
  MpReal res;
  mpfr_tan(res.mpfr_rep, r.mpfr_rep, rnd);
  return res;
}

inline friend MpReal acos(const MpReal& r, MpReal::RoundingMode rnd= MpReal::getDefaultRndMode()) {
  MpReal res;
  mpfr_acos(res.mpfr_rep, r.mpfr_rep, rnd);
  return res;
}

inline friend MpReal asin(const MpReal& r, MpReal::RoundingMode rnd= MpReal::getDefaultRndMode()) {
  MpReal res;
  mpfr_asin(res.mpfr_rep, r.mpfr_rep, rnd);
  return res;
}

inline friend MpReal atan(const MpReal& r, MpReal::RoundingMode rnd= MpReal::getDefaultRndMode()) {
  MpReal res;
  mpfr_atan(res.mpfr_rep, r.mpfr_rep, rnd);
  return res;
}

inline friend MpReal cosh(const MpReal& r, MpReal::RoundingMode rnd= MpReal::getDefaultRndMode()) {
  MpReal res;
  mpfr_cosh(res.mpfr_rep, r.mpfr_rep, rnd);
  return res;
}

inline friend MpReal sinh(const MpReal& r, MpReal::RoundingMode rnd= MpReal::getDefaultRndMode()) {
  MpReal res;
  mpfr_sinh(res.mpfr_rep, r.mpfr_rep, rnd);
  return res;
}

inline friend MpReal tanh(const MpReal& r, MpReal::RoundingMode rnd= MpReal::getDefaultRndMode()) {
  MpReal res;
  mpfr_tanh(res.mpfr_rep, r.mpfr_rep, rnd);
  return res;
}

inline friend MpReal asinh(const MpReal& r, MpReal::RoundingMode rnd= MpReal::getDefaultRndMode()) {
  MpReal res;
  mpfr_asinh(res.mpfr_rep, r.mpfr_rep, rnd);
  return res;
}

inline friend MpReal acosh(const MpReal& r, MpReal::RoundingMode rnd= MpReal::getDefaultRndMode()) {
  MpReal res;
  mpfr_acosh(res.mpfr_rep, r.mpfr_rep, rnd);
  return res;
}

inline friend MpReal atanh(const MpReal& r, MpReal::RoundingMode rnd= MpReal::getDefaultRndMode()) {
  MpReal res;
  mpfr_atanh(res.mpfr_rep, r.mpfr_rep, rnd);
  return res;
}

inline friend MpReal pow(const MpReal & x, const MpReal &e) {
  return x.pow(e, MpReal::getDefaultRndMode());
}
inline friend MpReal pow(const MpReal & x, const long int e) {
  return x.pow(e, MpReal::getDefaultRndMode());
}
inline friend MpReal power(const MpReal & x, const MpReal &e) {
  return x.pow(e, MpReal::getDefaultRndMode());
}
inline friend MpReal power(const MpReal & x, const long int e) {
  return x.pow(e, MpReal::getDefaultRndMode());
}


inline friend MpReal cbrt(const MpReal& r, MpReal::RoundingMode rnd= MpReal::getDefaultRndMode()) {
  MpReal res;
  mpfr_cbrt(res.mpfr_rep, r.mpfr_rep, rnd);
  return res;
}


inline friend MpReal gamma(const MpReal& r, MpReal::RoundingMode rnd= MpReal::getDefaultRndMode()) {
  MpReal res;
  mpfr_gamma(res.mpfr_rep, r.mpfr_rep, rnd);
  return res;
}

inline friend MpReal erf(const MpReal& r, MpReal::RoundingMode rnd= MpReal::getDefaultRndMode()) {
  MpReal res;
  mpfr_erf(res.mpfr_rep, r.mpfr_rep, rnd);
  return res;
}

inline friend MpReal hypot(const MpReal& r1, const MpReal& r2,
    MpReal::RoundingMode rnd= MpReal::getDefaultRndMode()) {
  MpReal res;
  mpfr_hypot(res.mpfr_rep, r1.mpfr_rep, r2.mpfr_rep, rnd);
  return res;
}

inline friend MpReal zeta(const MpReal& r, MpReal::RoundingMode rnd= MpReal::getDefaultRndMode()) {
  MpReal res;
  mpfr_zeta(res.mpfr_rep, r.mpfr_rep, rnd);
  return res;
}

inline friend long int toInt (const MpReal& r, RoundingMode rnd = MpReal::RoundToZero){
  return mpfr_get_si(r.mpfr_rep, rnd);
}

inline friend long int toLongInt (const MpReal& r, RoundingMode rnd = MpReal::RoundToZero){
  return mpfr_get_si(r.mpfr_rep, rnd);
}

inline friend unsigned long int toUInt (const MpReal& r, RoundingMode rnd = MpReal::RoundToZero){
  return mpfr_get_ui(r.mpfr_rep, rnd);
}

inline friend double toDouble (const MpReal& r, RoundingMode rnd = MpReal::getDefaultRndMode()){
  return mpfr_get_d(r.mpfr_rep, rnd);
}

inline friend long double toLongDouble (const MpReal& r, RoundingMode rnd = MpReal::getDefaultRndMode()){
  return mpfr_get_ld(r.mpfr_rep, rnd);
}

inline friend MpReal relDiff(const MpReal& r1, const MpReal& r2,
    MpReal::RoundingMode rnd = MpReal::getDefaultRndMode()) {
  MpReal res;
  mpfr_reldiff(res.mpfr_rep, r1.mpfr_rep, r2.mpfr_rep, rnd);
  return res;
}

inline friend MpReal nextAbove(const MpReal & r) {
  MpReal res(r);
  mpfr_nextabove(res.mpfr_rep);
  return res;
}

inline friend MpReal nextBelow(const MpReal& r) {
  MpReal res(r);
  mpfr_nextbelow(res.mpfr_rep);
  return res;
}

inline friend MpReal nextToward(const MpReal& r, const MpReal& dir) {
  MpReal res(r);
  mpfr_nexttoward(res.mpfr_rep, dir.mpfr_rep);
  return res;
}
