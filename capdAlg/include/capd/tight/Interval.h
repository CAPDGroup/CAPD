/////////////////////////////////////////////////////////////////////////////
//
/// @file tight/Interval.h
///
/// CAPD interface for the TIGHT interval library
/// (correctly rounded interval arithmetic, https://gitlab.com/fsichetti/tight).
///
/// It is an adapter that has the same interface as capd::filib::Interval,
/// so that the whole CAPD can be compiled with TIGHT intervals
/// (cmake -DCAPD_INTERVAL_TYPE=TIGHT -DCAPD_TIGHT_DIR=<path to tight>).
//
/////////////////////////////////////////////////////////////////////////////

// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef CAPD_TIGHT_INTERVAL_H
#define CAPD_TIGHT_INTERVAL_H

#include <cmath>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <sstream>
#include <string>
#include <typeinfo>
#include <iostream>
// elementary functions of TIGHT are declared only if TIGHT_CORE_MATH is defined
// (CAPD builds TIGHT with CORE-MATH)
#ifndef TIGHT_CORE_MATH
#define TIGHT_CORE_MATH
#endif
#include <tight/Interval.hpp>
#include "capd/intervals/IntervalError.h"
#include "capd/intervals/intra/interval_atan2.h"
#include "capd/basicalg/minmax.h"
#include "capd/basicalg/TypeTraits.h"
#include "capd/basicalg/Math.h"

namespace capd{
namespace tight{

/**
 *  CAPD interface for the TIGHT interval (double endpoints, correct rounding).
 *
 *  The TIGHT library requires the FPU to be in the rounding-up mode.
 *  This is done automatically when a program including this header starts
 *  (see capd::tight::Interval::setup()). Any code that changes the rounding mode
 *  (e.g. capd::rounding::DoubleRounding::roundNearest()) must restore round-up
 *  before TIGHT intervals are used again.
 *
 *  As in capd::filib::Interval, the TIGHT interval is not a base class and is
 *  accessible only by getBaseInterval(). Functions throw IntervalError for
 *  arguments outside of the domain (TIGHT itself would return NaN intervals).
 */
class Interval{
public:
  typedef ::tight::Interval BaseInterval;
  typedef double BoundType;
  typedef capd::intervals::IntervalError<BoundType> IntervalError;

protected:
  BaseInterval m_interval;

  /// converts decimal representation ([+-]digits[.digits][e[+-]digits]) to a rigorous enclosure
  static BaseInterval fromDecimal(const std::string & s){
    std::string mant, digits;
    long exp10 = 0;
    size_t i = 0;
    while(i < s.size() && isspace((unsigned char)s[i])) ++i;
    bool neg = false;
    if(i < s.size() && (s[i]=='+' || s[i]=='-')){ neg = s[i]=='-'; ++i; }
    long fracDigits = 0;
    bool seenDot = false;
    for(; i < s.size(); ++i){
      char c = s[i];
      if(c >= '0' && c <= '9'){ digits += c; if(seenDot) ++fracDigits; }
      else if(c == '.' && !seenDot) seenDot = true;
      else break;
    }
    if(i < s.size() && (s[i]=='e' || s[i]=='E')){
      exp10 = std::strtol(s.c_str() + i + 1, nullptr, 10);
    } else {
      while(i < s.size() && isspace((unsigned char)s[i])) ++i;
      if(i < s.size()) throw IntervalError("Interval(string): cannot parse number", 0, 0);
    }
    if(digits.empty())
      throw IntervalError("Interval(string): cannot parse number", 0, 0);
    exp10 -= fracDigits;
    std::string p = (neg ? "-" : "") + digits, q = "1";
    if(exp10 >= 0) p += std::string((size_t)exp10, '0');
    else           q += std::string((size_t)(-exp10), '0');
    return BaseInterval::fromRational(p, q);
  }

public:

  inline Interval() : m_interval(0.0) {}

  inline Interval(const Interval & A_iv) : m_interval(A_iv.m_interval) {}

  Interval(const BoundType & A_scalar) : m_interval(A_scalar) {}

  Interval(const BoundType & A_left, const BoundType & A_right) : m_interval(A_left, A_right) {}

  Interval(const BaseInterval & interval) : m_interval(interval) {}

  /// rigorous enclosure [left, right] of two decimal numbers
  Interval(const char left[], const char right[]) : m_interval(
      lowerOf(fromDecimal(left)), upperOf(fromDecimal(right))) {}

  Interval(const std::string & left, const std::string & right) : m_interval(
      lowerOf(fromDecimal(left)), upperOf(fromDecimal(right))) {}

  /// returns reference to the base (TIGHT) interval
  BaseInterval & getBaseInterval(){ return m_interval; }
  const BaseInterval & getBaseInterval() const { return m_interval; }

  /// sets up the floating point unit (round up). Called automatically at program start.
  static void setup(){ BaseInterval::init(); }
  /// restores the default rounding mode (to the nearest)
  static void reset(){ BaseInterval::deinit(); }

  BoundType leftBound() const { return lower(m_interval); }
  BoundType rightBound() const { return upper(m_interval); }
  friend inline BoundType leftBound(const Interval & x){ return lower(x.m_interval); }
  friend inline BoundType rightBound(const Interval & x){ return upper(x.m_interval); }

  BoundType inf() const { return lower(m_interval); }
  BoundType sup() const { return upper(m_interval); }
  friend inline BoundType inf(const Interval & x){ return lower(x.m_interval); }
  friend inline BoundType sup(const Interval & x){ return upper(x.m_interval); }

  /// @deprecated
  void setLeftBound(const BoundType & A_left){ m_interval = BaseInterval(A_left, upper(m_interval)); }
  /// @deprecated
  void setRightBound(const BoundType & A_right){ m_interval = BaseInterval(lower(m_interval), A_right); }

  Interval left() const { return Interval(leftBound()); }
  Interval right() const { return Interval(rightBound()); }
  friend inline Interval left(const Interval & x){ return Interval(x.leftBound()); }
  friend inline Interval right(const Interval & x){ return Interval(x.rightBound()); }

  template <typename T_Scalar>
  bool contains(const T_Scalar & A_X) const { return leftBound() <= (BoundType)A_X && (BoundType)A_X <= rightBound(); }
  bool contains(const Interval & A_iv) const { return A_iv.subset(*this); }

  template <typename T_Scalar>
  bool containsInInterior(const T_Scalar & A_X) const { return (leftBound() < (BoundType)A_X) && ((BoundType)A_X < rightBound()); }
  bool containsInInterior(const Interval & A_iv) const { return (leftBound() < A_iv.leftBound()) && (A_iv.rightBound() < rightBound()); }

  bool subset(const Interval & A_iv) const { return A_iv.leftBound() <= leftBound() && rightBound() <= A_iv.rightBound(); }
  bool subsetInterior(const Interval & A_iv) const { return A_iv.containsInInterior(*this); }
  friend bool subset(const Interval & A_iv1, const Interval & A_iv2){ return A_iv1.subset(A_iv2); }
  friend bool subsetInterior(const Interval & A_iv1, const Interval & A_iv2){ return A_iv1.subsetInterior(A_iv2); }

  /// returns the middle point of the interval
  Interval mid() const { return Interval(midPoint()); }
  Interval abs() const { return Interval(::tight::abs(m_interval)); }

  /// Splits interval into the form  mid + remainder, where mid is the middle point
  void split(Interval & A_rMid, Interval & A_rRemainder) const {
    if(std::isfinite(leftBound()) && std::isfinite(rightBound())){
      BoundType m = midPoint();
      A_rRemainder.m_interval = m_interval - m;
      A_rMid.m_interval = BaseInterval(m);
    } else {
      A_rMid = Interval(0.0);
      A_rRemainder = Interval(-std::numeric_limits<BoundType>::infinity(), std::numeric_limits<BoundType>::infinity());
    }
  }
  void split(BoundType & A_rMid, Interval & A_rRemainder) const {
    if(std::isfinite(leftBound()) && std::isfinite(rightBound())){
      A_rMid = midPoint();
      A_rRemainder.m_interval = m_interval - A_rMid;
    } else {
      A_rMid = 0.0;
      A_rRemainder = Interval(-std::numeric_limits<BoundType>::infinity(), std::numeric_limits<BoundType>::infinity());
    }
  }
  void split(Interval & r){ split(*this, r); }

  static Interval pi(){ return Interval(::tight::constants::pi); }
  static Interval euler(){ return Interval(exp(BaseInterval(1.0))); }

  Interval & operator = (const Interval & A_iv){ m_interval = A_iv.m_interval; return *this; }
  Interval & operator = (const BoundType & A_x){ m_interval = BaseInterval(A_x); return *this; }
  Interval & operator += (const Interval & A_iv){ m_interval += A_iv.m_interval; return *this; }
  Interval & operator -= (const Interval & A_iv){ m_interval -= A_iv.m_interval; return *this; }
  Interval & operator *= (const Interval & A_iv){ m_interval *= A_iv.m_interval; return *this; }
  Interval & operator /= (const Interval & A_iv){
    if(A_iv.contains(0.0))
      throw IntervalError("Division by 0 in operator/=(Interval, Interval)", A_iv.inf(), A_iv.sup());
    m_interval /= A_iv.m_interval;
    return *this;
  }

  friend bool operator == (const Interval & a, const Interval & b){ return a.leftBound()==b.leftBound() && a.rightBound()==b.rightBound(); }
  friend bool operator <= (const Interval & a, const Interval & b){ return a.rightBound() <= b.leftBound(); }
  friend bool operator >= (const Interval & a, const Interval & b){ return a.leftBound() >= b.rightBound(); }
  friend bool operator <  (const Interval & a, const Interval & b){ return a.rightBound() < b.leftBound(); }
  friend bool operator >  (const Interval & a, const Interval & b){ return a.leftBound() > b.rightBound(); }
  friend bool operator != (const Interval & a, const Interval & b){ return a.leftBound()!=b.leftBound() || a.rightBound()!=b.rightBound(); }

  friend inline bool operator == (const Interval & a, const BoundType & b){ return a.leftBound()==b && a.rightBound()==b; }
  friend inline bool operator == (const BoundType & a, const Interval & b){ return a==b.leftBound() && a==b.rightBound(); }
  friend inline bool operator != (const Interval & a, const BoundType & b){ return a.leftBound()!=b || a.rightBound()!=b; }
  friend inline bool operator != (const BoundType & a, const Interval & b){ return a!=b.leftBound() || a!=b.rightBound(); }
  friend inline bool operator >  (const Interval & a, const BoundType & b){ return a.leftBound() > b; }
  friend inline bool operator >  (const BoundType & a, const Interval & b){ return a > b.rightBound(); }
  friend inline bool operator >= (const Interval & a, const BoundType & b){ return a.leftBound() >= b; }
  friend inline bool operator >= (const BoundType & a, const Interval & b){ return a >= b.rightBound(); }
  friend inline bool operator <  (const Interval & a, const BoundType & b){ return a.rightBound() < b; }
  friend inline bool operator <  (const BoundType & a, const Interval & b){ return a < b.leftBound(); }
  friend inline bool operator <= (const Interval & a, const BoundType & b){ return a.rightBound() <= b; }
  friend inline bool operator <= (const BoundType & a, const Interval & b){ return a <= b.leftBound(); }

  //==== arithmetic ======================================================================

  friend inline Interval operator - (const Interval & a){ return Interval(-a.m_interval); }
  friend inline Interval operator + (const Interval & a, const Interval & b){ return Interval(a.m_interval + b.m_interval); }
  friend inline Interval operator - (const Interval & a, const Interval & b){ return Interval(a.m_interval - b.m_interval); }
  friend inline Interval operator * (const Interval & a, const Interval & b){ return Interval(a.m_interval * b.m_interval); }
  friend inline Interval operator / (const Interval & a, const Interval & b){
    if(!b.contains(0.0))
      return Interval(a.m_interval / b.m_interval);
    throw IntervalError("Division by 0 in operator/(Interval, Interval)", b.inf(), b.sup());
  }

  /// this template prevents implicit conversion of the exponent type
  template<class E>
  friend inline Interval operator ^ (const Interval & x, E){
    std::ostringstream out;
    out << " power(interval x, E n): exponent type typeid(E).name()=" << typeid(E).name() << " is not supported. \n";
    throw IntervalError(out.str().c_str(), x.leftBound(), x.rightBound());
  }

  friend inline Interval operator ^ (const Interval & x, unsigned int n){
    if(n <= (unsigned)INT_MAX)
      return Interval(::tight::pown(x.m_interval, (int)n));
    int m = (int)(n >> 1);
    Interval h = x ^ (unsigned int)m;
    return (n & 1) ? x * sqr(h) : sqr(h);
  }

  friend inline Interval operator ^ (const Interval & x, int n){
    if(n >= 0 || !isSingular(x))
      return Interval(::tight::pown(x.m_interval, n));
    throw IntervalError(" power(interval x, negativeInt): Interval x must be nonzero. \n", x.leftBound(), x.rightBound());
  }

  friend inline Interval operator ^ (const Interval & x, BoundType b){
    if(x.leftBound() > 0.0)
      return Interval(::tight::pow(x.m_interval, BaseInterval(b)));
    throw IntervalError(" power(interval x, double): Interval x must be positive. \n", x.leftBound(), x.rightBound());
  }

  friend inline Interval operator ^ (const Interval & a, const Interval & b){
    if(a.leftBound() > 0)
      return Interval(::tight::pow(a.m_interval, b.m_interval));
    throw IntervalError(" power(interval A, interval B): Interval A must be positive. \n", a.leftBound(), a.rightBound());
  }

  //==== arithmetic between interval and BoundType =======================================

  friend inline Interval operator + (const Interval & a, const BoundType & x){ return Interval(a.m_interval + x); }
  friend inline Interval operator + (const BoundType & x, const Interval & a){ return Interval(a.m_interval + x); }
  friend inline Interval operator - (const Interval & a, const BoundType & x){ return Interval(a.m_interval - x); }
  friend inline Interval operator - (const BoundType & x, const Interval & a){ return Interval(x - a.m_interval); }
  friend inline Interval operator * (const Interval & a, const BoundType & x){ return Interval(a.m_interval * x); }
  friend inline Interval operator * (const BoundType & x, const Interval & a){ return Interval(a.m_interval * x); }
  friend inline Interval operator / (const Interval & a, const BoundType & x){
    if(x != 0.0)
      return Interval(a.m_interval / x);
    throw IntervalError("Division by 0 in operator/(Interval, BoundType)", x, x);
  }
  friend inline Interval operator / (const BoundType & x, const Interval & a){
    if(!a.contains(0.0))
      return Interval(x / a.m_interval);
    throw IntervalError("Division by 0 in operator/(BoundType, Interval)", a.inf(), a.sup());
  }

  //==== input / output ==================================================================

  friend std::ostream & operator << (std::ostream & s, const Interval & A_iv){
    s << "[" << A_iv.leftBound() << ", " << A_iv.rightBound() << "]";
    return s;
  }

  /// reads an interval written as [a, b], a single number a or "[a]"; decimal numbers are rounded outwards
  friend std::istream & operator >> (std::istream & s, Interval & A_iv){
    std::string tok;
    char c;
    s >> std::ws;
    if(s.peek() != '['){
      s >> tok;
      A_iv = Interval(tok, tok);
      return s;
    }
    s.get(c);
    std::string a, b;
    s >> std::ws;
    while(s.get(c) && c != ',' && c != ']') a += c;
    if(c == ',' ){
      while(s.get(c) && c != ']') b += c;
    } else b = a;
    try {
      A_iv = Interval(a, b);
    } catch(IntervalError &) {
      s.setstate(std::ios::failbit);
    }
    return s;
  }

  friend std::ostream & bitWrite(std::ostream & out, const Interval & iv){
    return out << bitsOf(iv.leftBound()) << ' ' << bitsOf(iv.rightBound());
  }
  friend std::istream & bitRead(std::istream & in, Interval & iv){
    std::string a, b;
    in >> a >> b;
    iv = Interval(fromBits(a), fromBits(b));
    return in;
  }
  friend std::ostream & hexWrite(std::ostream & out, const Interval & iv){
    char buf[64];
    std::snprintf(buf, sizeof buf, "%a %a", iv.leftBound(), iv.rightBound());
    return out << buf;
  }
  friend std::istream & hexRead(std::istream & in, Interval & iv){
    std::string a, b;
    in >> a >> b;
    iv = Interval(std::strtod(a.c_str(), nullptr), std::strtod(b.c_str(), nullptr));
    return in;
  }
  friend std::ostream & binWrite(std::ostream & out, const Interval & iv){
    BoundType t[2] = { iv.leftBound(), iv.rightBound() };
    return out.write((const char *)t, sizeof t);
  }
  friend std::istream & binRead(std::istream & in, Interval & iv){
    BoundType t[2];
    in.read((char *)t, sizeof t);
    if(in) iv = Interval(t[0], t[1]);
    return in;
  }

  //==== functions on intervals ==========================================================

  friend inline Interval diam(const Interval & ix){
    return Interval(BaseInterval(ix.rightBound()) - BaseInterval(ix.leftBound()));
  }
  friend inline BoundType width(const Interval & ix){
    return ix.rightBound() - ix.leftBound();
  }
  friend inline Interval mid(const Interval & A_iv){ return A_iv.mid(); }

  /// Intersection of two intervals; returns false if it is empty
  friend bool intersection(const Interval & a, const Interval & b, Interval & r){
    BoundType l = (a.leftBound() > b.leftBound()) ? a.leftBound() : b.leftBound();
    BoundType u = (a.rightBound() < b.rightBound()) ? a.rightBound() : b.rightBound();
    if(l <= u){
      r = Interval(l, u);
      return true;
    }
    return false;
  }

  /// returns an interval containing ix and iy
  friend Interval intervalHull(const Interval & ix, const Interval & iy){
    BoundType l = (ix.leftBound() < iy.leftBound()) ? ix.leftBound() : iy.leftBound(),
              u = (ix.rightBound() > iy.rightBound()) ? ix.rightBound() : iy.rightBound();
    return Interval(l, u);
  }

  //  On output:  iv \subset Mid + [-diam , diam]
  friend inline void split(Interval & A_iv, Interval & A_rMid, BoundType & A_diam){
    BoundType m = A_iv.midPoint();
    Interval t = A_iv - m;
    A_rMid = Interval(m);
    A_diam = (-t.leftBound() > t.rightBound()) ? -t.leftBound() : t.rightBound();
  }

  friend inline void split(Interval & A_rIv, BoundType & A_diam){
    BoundType m = A_rIv.midPoint();
    Interval t = A_rIv - m;
    A_rIv = Interval(m);
    A_diam = (-t.leftBound() > t.rightBound()) ? -t.leftBound() : t.rightBound();
  }

  friend inline bool isSingular(const Interval & A_x){
    return (A_x.leftBound() <= 0) && (A_x.rightBound() >= 0);
  }

  template<class E>
  friend inline Interval power(const Interval & x, E n){ return x ^ n; }

  friend inline Interval sqrt(const Interval & x){
    if(x.leftBound() < 0)
      throw IntervalError(" sqrt(x): Interval x must be nonnegative. \n", x.leftBound(), x.leftBound());
    return Interval(sqrt(x.m_interval));
  }

  friend inline Interval sin(const Interval & x){ return Interval(sin(x.m_interval)); }
  friend inline Interval cos(const Interval & x){ return Interval(cos(x.m_interval)); }
  friend inline Interval tan(const Interval & x){ return Interval(tan(x.m_interval)); }

  /// cot x = 1/tan x. Throws if x contains a zero of tan (a pole of cot) or a pole of tan.
  friend inline Interval cot(const Interval & x){
    BaseInterval t = tan(x.m_interval);
    if(std::isnan(lower(t)) || std::isnan(upper(t)) || (lower(t) <= 0.0 && upper(t) >= 0.0))
      throw IntervalError(" cot(A): Interval A contains a singularity of cot or tan. \n", x.leftBound(), x.rightBound());
    return Interval(BaseInterval(1.0) / t);
  }

  friend inline Interval atan(const Interval & x){ return Interval(atan(x.m_interval)); }
  friend inline Interval atan2(const Interval & x, const Interval & y){
    return capd::intervals::intra::atan2(x, y);
  }

  friend inline Interval asin(const Interval & x){
    if((x.leftBound() < -1) || (x.rightBound() > 1))
      throw IntervalError(" asin(A): Interval A must satisfy -1 <= A <= 1. \n", x.leftBound(), x.rightBound());
    return Interval(asin(x.m_interval));
  }

  friend inline Interval acos(const Interval & x){
    if((x.leftBound() < -1) || (x.rightBound() > 1))
      throw IntervalError(" acos(A): Interval A must satisfy -1 <= A <= 1. \n", x.leftBound(), x.rightBound());
    return Interval(acos(x.m_interval));
  }

  friend inline Interval sinh(const Interval & x){ return Interval(sinh(x.m_interval)); }
  friend inline Interval cosh(const Interval & x){ return Interval(cosh(x.m_interval)); }
  friend inline Interval tanh(const Interval & x){ return Interval(tanh(x.m_interval)); }

  /// coth x = 1/tanh x. Throws if x contains 0.
  friend inline Interval coth(const Interval & x){
    if(x.leftBound() <= 0.0 && x.rightBound() >= 0.0)
      throw IntervalError(" coth(A): Interval A must not contain 0. \n", x.leftBound(), x.rightBound());
    return Interval(BaseInterval(1.0) / tanh(x.m_interval));
  }

  friend inline Interval exp(const Interval & x){ return Interval(exp(x.m_interval)); }

  friend inline Interval log(const Interval & x){
    if(x.leftBound() <= 0.0)
      throw IntervalError(" log(A): Interval A must satisfy 0 < A. \n", x.leftBound(), x.rightBound());
    return Interval(log(x.m_interval));
  }

  friend inline Interval sqr(const Interval & x){
    if(capd::abs(x.leftBound()) > 1.3407807929942596e154 || capd::abs(x.rightBound()) > 1.3407807929942596e154)
      throw IntervalError(" sqr(A): sqr overflow. \n", x.leftBound(), x.rightBound());
    return Interval(sqr(x.m_interval));
  }

  /// returns true if any of the endpoints is infinite
  friend inline bool isinf(const Interval & r){ return std::isinf(r.leftBound()) || std::isinf(r.rightBound()); }
  /// returns true if any of the endpoints is a NaN
  friend inline bool isnan(const Interval & r){ return std::isnan(r.leftBound()) || std::isnan(r.rightBound()); }

  /// returns interval containing copysign(a, b) for each a in x and b in y
  friend inline Interval copysign(const Interval & x, const Interval & y){
    if(y.contains(0.0)){
      BoundType m = capd::max(-x.leftBound(), x.rightBound());
      return Interval(-m, m);
    }
    return y < 0 ? -x : x;
  }

  /// returns the nonnegative part of the interval
  friend inline Interval nonnegativePart(const Interval & iv){
    if(iv.rightBound() < 0.0)
      throw IntervalError(" Nonnegative part is empty! ", iv.leftBound(), iv.rightBound());
    return Interval(capd::max(iv.leftBound(), BoundType(0.0)), iv.rightBound());
  }

  /// Ball with center iv and radius r
  friend inline Interval ball(const Interval & iv, const Interval & r){
    return Interval(iv.leftBound() - r.rightBound(), iv.rightBound() + r.rightBound());
  }
  friend inline Interval ball(const Interval & iv, const BoundType & r){
    return Interval(iv.leftBound() - r, iv.rightBound() + r);
  }

  /// solves inclusion a+[0,t]*p \subset c for t
  friend Interval solveAffineInclusion(const Interval & a, const Interval & p, const Interval & c){
    Interval t;
    if(!(a.subset(c))){
      throw IntervalError("Cannot solve affine inclusion", a.leftBound(), a.rightBound());
    } else if(p >= 0){
      t = ((c.right() - a.right()) / p.right()).left();
    } else if(p <= 0){
      t = ((c.left() - a.left()) / p.left()).left();
    } else {
      BoundType t1 = ((c.right() - a.right()) / p.right()).leftBound();
      BoundType t2 = ((c.left() - a.left()) / p.left()).leftBound();
      t = Interval((t1 < t2 ? t1 : t2));
    }
    return t;
  }

  /// solves inclusion a+[0,t]*p \subset c for t
  friend Interval solveAffineInclusion(const Interval & a, const Interval & p, const Interval & c, int & dir){
    Interval t;
    if(!(a.subset(c))){
      throw IntervalError("Cannot solve affine inclusion", a.leftBound(), a.rightBound());
    } else if(p >= 0){
      t = ((c.right() - a.right()) / p.right()).left();
      dir = 1;
    } else if(p <= 0){
      t = ((c.left() - a.left()) / p.left()).left();
      dir = -1;
    } else {
      BoundType t1 = ((c.right() - a.right()) / p.right()).leftBound();
      BoundType t2 = ((c.left() - a.left()) / p.left()).leftBound();
      if(t1 < t2){
        dir = 1;
        t = Interval(t1);
      } else {
        dir = -1;
        t = Interval(t2);
      }
    }
    return t;
  }

private:
  static BoundType lowerOf(const BaseInterval & i){ return lower(i); }
  static BoundType upperOf(const BaseInterval & i){ return upper(i); }
  BoundType midPoint() const { return m_interval.middle(); }

  static std::string bitsOf(double d){
    unsigned long long u;
    static_assert(sizeof u == sizeof d, "unexpected double size");
    std::memcpy(&u, &d, sizeof u);
    std::string s(64, '0');
    for(int i = 0; i < 64; ++i) if((u >> (63 - i)) & 1ULL) s[i] = '1';
    return s;
  }
  static double fromBits(const std::string & s){
    unsigned long long u = 0;
    for(size_t i = 0; i < s.size() && i < 64; ++i) u = (u << 1) | (s[i] == '1');
    double d;
    std::memcpy(&d, &u, sizeof d);
    return d;
  }
};

/// the TIGHT library needs the FPU in the rounding-up mode: set it when the program starts
/// (inline variable, so it is done once for the whole program)
inline const bool intervalSetupDone = (Interval::setup(), true);

}} // namespace capd::tight

namespace capd{

/**
 * Specialization of TypeTraits for TIGHT intervals
 */
template <>
class TypeTraits< capd::tight::Interval > {
public:
  typedef double Real;
  typedef ::capd::tight::Interval IntervalType;

  static inline const IntervalType & zero(){
    static const IntervalType z(0.0);
    return z;
  }
  static inline const IntervalType & one(){
    static const IntervalType o(1.0);
    return o;
  }
  static inline int numberOfDigits(){
    return TypeTraits<double>::numberOfDigits();
  }
  /// Machine epsilon
  static inline double epsilon() noexcept {
    return TypeTraits<double>::epsilon();
  }

  template <typename S>
  static inline IntervalType convert(const S & obj){
    return static_cast<IntervalType>(TypeTraits<double>::convert(obj));
  }
  static inline IntervalType convert(const double & obj){
    return static_cast<IntervalType>(obj);
  }

  static const bool isInterval = true;

  /// an absolute value
  static inline IntervalType abs(const IntervalType & A_inter){
    return IntervalType(A_inter.abs());
  }

  /// maximum
  static inline IntervalType max(const IntervalType & A_iv1, const IntervalType & A_iv2){
    return IntervalType(
      (A_iv1.leftBound() > A_iv2.leftBound() ? A_iv1.leftBound() : A_iv2.leftBound()),
      (A_iv1.rightBound() > A_iv2.rightBound() ? A_iv1.rightBound() : A_iv2.rightBound())
    );
  }

  /// minimum
  static inline IntervalType min(const IntervalType & A_iv1, const IntervalType & A_iv2){
    return IntervalType(
      (A_iv1.leftBound() < A_iv2.leftBound() ? A_iv1.leftBound() : A_iv2.leftBound()),
      (A_iv1.rightBound() < A_iv2.rightBound() ? A_iv1.rightBound() : A_iv2.rightBound())
    );
  }

  static inline bool isSingular(const IntervalType & x){
    return (x.leftBound() <= 0.0 && x.rightBound() >= 0.0);
  }

  static inline int _int(const IntervalType & z) noexcept { return TypeTraits<double>::_int(z.leftBound()); }
  static inline double _double(const IntervalType & z) noexcept { return TypeTraits<double>::_double(z.leftBound()); }
};

template <>
class Math< tight::Interval > {
public:
  typedef tight::Interval Interval;
  static inline Interval _sqr(const Interval & z){ return sqr(z); }
  static inline Interval _log(const Interval & x){ return log(x); }
  static inline Interval _pow(const Interval & x, int c){ return power(x, c); }
  static inline Interval _sqrt(const Interval & z){ return sqrt(z); }
  static inline Interval _exp(const Interval & x){ return exp(x); }
  static inline Interval _sin(const Interval & x){ return sin(x); }
  static inline Interval _cos(const Interval & x){ return cos(x); }
  static inline Interval _atan(const Interval & x){ return atan(x); }
  static inline Interval _asin(const Interval & x){ return asin(x); }
  static inline Interval _acos(const Interval & x){ return acos(x); }
};

} // namespace capd

#endif // CAPD_TIGHT_INTERVAL_H
