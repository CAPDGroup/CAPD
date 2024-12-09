/////////////////////////////////////////////////////////////////////////////
//
/// @file cxsc/Interval.h
///
/// @author Tomasz Kapela   @date 2010-06-08
//
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) Tomasz Kapela 2010
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_CXSC_INTERVAL_H_
#define _CAPD_CXSC_INTERVAL_H_
#include <interval.hpp>
#include <except.hpp>
#undef ZERO
#include "capd/intervals/IntervalError.h"
#include "capd/basicalg/minmax.h"
#include "capd/basicalg/TypeTraits.h"
#include "capd/intervals/IntervalTraits.h"



namespace capd{
namespace cxsc{

class Interval;

inline Interval diam(const Interval & ix);
}
/// absolute values of all elements of a given interval
template<>
cxsc::Interval abs (const cxsc::Interval & A_inter);

/// maximum
template<>
inline cxsc::Interval max(const cxsc::Interval& A_iv1, const cxsc::Interval& A_iv2);

///minimum
template<>
inline cxsc::Interval min (const cxsc::Interval& A_iv1, const cxsc::Interval& A_iv2);

/// fast interval library
namespace cxsc{

////////////////////////////////////////////////////////////////////////////////////////
//  capd::cxsc::Interval
///
///
/// CAPD interface for interval library cxsc
///
/// It works as an adapter for the other CAPD routines.
/// It has exactly the same interface as native CAPD interval
/// and it just calls the corresponding functions from the cxsc.
///
///
/// To prevent misusing and enable easy interchange of a different interval implementation
/// the cxsc interval T is not a base class of the Interval class
/// and can be accesed only by getBaseInterval() method.
///
//////////////////////////////////////////////////////////////////////////////////////////

class Interval{
public:
  typedef ::cxsc::interval BaseInterval;
  typedef double BoundType;
//  typedef ::cxsc::fp_traits<T, R> RoundingPolicy;
  typedef capd::intervals::IntervalError<BoundType> IntervalError;
  
protected:
  typedef ::cxsc::real Real;
  BaseInterval m_interval;
  Interval(const Real & A_scalar) : m_interval(A_scalar){
  }
  Interval(const Real & left, const Real & right) : m_interval(left, right){
  }
public:

  inline Interval(){}

  /// copying constructor
  inline Interval( const Interval & A_iv ) : m_interval(A_iv.m_interval){
  }

  /// constructor from any class that can be coverted to BoundType
  // template < typename T_Scalar >
  // Interval( const T_Scalar & A_scalar );
  Interval(const BoundType & A_scalar) : m_interval(A_scalar){
  }
  /// constructor from any class that can be coverted to BoundType
  //template < typename T_Scalar1, typename T_Scalar2 >
  //Interval( const T_Scalar1 & A_left, const T_Scalar2 & A_right );
  Interval( const BoundType & A_left, const BoundType & A_right ) : m_interval(A_left, A_right){
  }

  Interval( const BaseInterval & interval) : m_interval(interval){
  }

  Interval(const char left[], const char right[]){
    Real l, r;
    std::string(left) >> l;
    std::string(right) >> r;
    m_interval = BaseInterval(l,r);
  }

  Interval(const std::string & left, const std::string & right) {
    Real l, r;
    left >> l;
    right >> r;
    m_interval = BaseInterval(l,r);

  }

  /**
   *  returns reference to base (cxsc) interval
   */
  BaseInterval & getBaseInterval(){
    return m_interval;
  }

  /**
   *  returns const reference to base (cxsc) interval
   */
  const BaseInterval & getBaseInterval() const {
    return m_interval;
  }
 BoundType inf() const  ///<  returns the left end of the interval
      {  return ::cxsc::_double(::cxsc::Inf(m_interval)); }
  BoundType sup() const ///<  returns the right end of the interval
      { return ::cxsc::_double(::cxsc::Sup(m_interval)); }
  
  BoundType leftBound() const  ///<  returns the left end of the interval
      {  return inf(); }
  BoundType rightBound() const ///<  returns the right end of the interval
      { return sup(); }

  friend inline BoundType leftBound(const Interval & x)  ///<  returns the left end of the interval
      {  return x.inf(); }
  friend inline BoundType rightBound(const Interval & x)  ///<  returns the right end of the interval
      { return x.sup(); }

  

  friend inline BoundType inf(const Interval & x)  ///<  returns the left end of the interval
      {  return x.inf(); }
  friend inline BoundType sup(const Interval & x)  ///<  returns the right end of the interval
      { return x.sup(); }

  /// @deprecated
  void setLeftBound(const BoundType & A_left)
      { ::cxsc::SetInf(m_interval, A_left);}
  /// @deprecated
  void setRightBound(const BoundType & A_right)
      { ::cxsc::SetSup(m_interval, A_right);}

  Interval left() const      ///< returns interval containing left end
      { return Interval(::cxsc::Inf(m_interval)); }
  Interval right() const     ///< returns interval containing right end
      { return Interval(::cxsc::Sup(m_interval)); }

  friend inline Interval left(const Interval & x)      ///< returns interval containing left end
      { return Interval(::cxsc::Inf(x.m_interval)); }
  friend inline Interval right(const Interval & x)     ///< returns interval containing right end
      { return Interval(::cxsc::Sup(x.m_interval)); }

  template <typename T_Scalar>
  bool contains( const T_Scalar & A_X ) const     ///< checks if interval contains given point X
      { return ( inf() <= A_X ) and (A_X <= sup());}
  bool contains( const Interval & A_iv ) const  ///< checks if interval contains given interval iv
  {  
 //   std::cout << " L : " << ( inf() - A_iv.inf()) << " R : " << (sup() - A_iv.sup())<< "   " <<  ::cxsc::in( A_iv.m_interval, m_interval) << std::endl;  
    return (A_iv.inf() >= inf() ) and (A_iv.sup() <= sup());
//     std::cout << ((*this) - A_iv) << "   " <<  ::cxsc::in( A_iv.m_interval, m_interval) << std::endl;  
//    return ::cxsc::in( A_iv.m_interval, m_interval);
  }
  template <typename T_Scalar>
  bool containsInInterior( const T_Scalar & A_X ) const     ///< checks if interval contains in interior given point X
      { return (leftBound() < (BoundType)A_X) && ((BoundType)A_X < rightBound()); }

  bool containsInInterior( const Interval & A_iv ) const    ///< checks if interval contains in interior given interval iv
      { return ((leftBound() < A_iv.leftBound()) && (A_iv.rightBound() < rightBound())); }

  bool subset( const Interval & A_iv ) const          ///< checks if interval is subset of iv
      { return  (A_iv.inf() <= inf() ) and ( sup() <= A_iv.sup()); }
  bool subsetInterior( const Interval & A_iv ) const  ///< checks if interval is subset of interior of iv
      { return A_iv.containsInInterior(*this); }
  friend bool subset(const Interval & A_iv1, const Interval & A_iv2)
      { return A_iv1.subset(A_iv2); }
  friend bool subsetInterior(const Interval & A_iv1, const Interval & A_iv2)
      { return A_iv1.subsetInterior(A_iv2); }

  BoundType midPoint() const           ///< returns middle point of interval
      { return ::cxsc::_double(::cxsc::mid(m_interval)); }
  Interval mid() const           ///< returns middle point of interval
      { return ::cxsc::mid(m_interval); }
  Interval abs() const           ///<Returns the interval of absolute values of this interval, i.e.
      { return ::cxsc::abs(m_interval); }
  /// Splits interval into the form  mid + remainder, where mid - is middle point
  void split( Interval & A_rMid, Interval & A_rRemainder ) const {
    if(std::isfinite(m_interval.inf()) && std::isfinite(m_interval.sup())){
      BoundType m = midPoint();
      A_rRemainder.m_interval = m_interval - m;
      A_rMid.m_interval = m;
    } else {
      A_rMid = Interval(0.0, 0.0);
      A_rRemainder = Interval(-std::numeric_limits<BoundType>::infinity(), std::numeric_limits<BoundType>::infinity() );
    }
  }
  void split( BoundType & A_rMid, Interval & A_rRemainder ) const {
    if(std::isfinite(m_interval.inf()) && std::isfinite(m_interval.sup())){
      A_rMid = midPoint();
      A_rRemainder.m_interval = m_interval - A_rMid;
    } else {
      A_rMid = 0.0;
      A_rRemainder = Interval(-std::numeric_limits<BoundType>::infinity(), std::numeric_limits<BoundType>::infinity() );
    }
  }
  void split(Interval &r)
      { split(*this, r); }

  // "Constants" (but they depend on the bound type and the precision)
  static Interval pi();     ///< returns pi constant
  static Interval euler();  ///< returns euler constant

  Interval & operator = ( const Interval & A_iv ){
    m_interval = A_iv.m_interval;
    return *this;
  }
  Interval & operator = ( const BoundType & A_x){
    m_interval = A_x;
    return *this;
  }
  Interval & operator += ( const Interval & A_iv ){
    m_interval += A_iv.m_interval;
    return *this;
  }
  Interval & operator -= ( const Interval & A_iv ){
    m_interval -= A_iv.m_interval;
    return *this;
  }
  Interval & operator *= ( const Interval & A_iv ){
    m_interval *= A_iv.m_interval;
    return *this;
  }
  Interval & operator /= ( const Interval & A_iv ){
    m_interval /= A_iv.m_interval;
    return *this;
  }

  friend bool operator ==  ( const Interval & A_iv1, const Interval & A_iv2 ) {
    return ((A_iv1.leftBound()==A_iv2.leftBound()) && (A_iv1.rightBound()==A_iv2.rightBound()));
  }
  friend bool operator <=  ( const Interval & A_iv1, const Interval & A_iv2 ) {
    return (A_iv1.rightBound() <= A_iv2.leftBound());
  }
  friend bool operator >=  ( const Interval & A_iv1, const Interval & A_iv2 ) {
    return (A_iv1.leftBound()>= A_iv2.rightBound() );
  }
  friend bool operator <   ( const Interval & A_iv1, const Interval & A_iv2 ) {
    return (A_iv1.rightBound() < A_iv2.leftBound());
  }
  friend bool operator >   ( const Interval & A_iv1, const Interval & A_iv2 ) {
    return (A_iv1.leftBound()> A_iv2.rightBound());
  }
  friend bool operator !=  ( const Interval & A_iv1, const Interval & A_iv2 ) {
    return ((A_iv1.leftBound()!= A_iv2.leftBound()) || (A_iv1.rightBound() != A_iv2.rightBound()));
  }


  ///  operator == (interval, scalar)
  friend inline bool operator == (const Interval & A_iVal1, const BoundType & A_Val2) {
    return ((A_iVal1.leftBound() == (A_Val2)) && (A_iVal1.rightBound() == (A_Val2)));
  }

  ///  operator == (scalar, interval)
  friend inline bool operator == (const BoundType & A_Val1, const Interval & A_iVal2) {
     return ((A_Val1 == A_iVal2.leftBound()) && (A_Val1 == A_iVal2.rightBound()));
  }

  ///  operator !=  (interval, scalar)
  friend inline bool operator != (const Interval & A_iVal1, const BoundType & A_Val2) {
    return ((A_iVal1.leftBound() != A_Val2) || (A_iVal1.rightBound() != A_Val2));
  }

  ///  operator !=   (scalar, interval)
  friend inline bool operator != (const BoundType & A_Val1, const Interval & A_iVal2) {
     return ((A_Val1 != A_iVal2.leftBound()) || (A_Val1 != A_iVal2.rightBound()));
  }

  ///  operator >  (interval, scalar)
  friend inline bool operator > (const Interval & A_iVal1,  const BoundType & A_Val2) {
     return (A_iVal1.leftBound() > A_Val2);
  }

  ///  operator >   (scalar, interval)
  friend inline bool operator > (const BoundType & A_Val1, const Interval & A_iVal2) {
     return (A_Val1 > A_iVal2.rightBound());
  }

  ///  operator >= (interval, scalar)
  friend inline bool operator >= (const Interval & A_iVal1, const BoundType & A_Val2) {
     return (A_iVal1.leftBound() >= A_Val2);
  }

  ///  operator >=  (scalar, interval)
  friend inline bool operator >= (const BoundType & A_Val1, const Interval & A_iVal2) {
     return (A_Val1 >= A_iVal2.rightBound());
  }

  ///  operator  < (interval, scalar)
  friend inline bool operator < (const Interval & A_iVal1, const BoundType & A_Val2) {
     return (A_iVal1.rightBound() < A_Val2);
  }

  ///  operator  <  (scalar, interval)
  friend inline bool operator < (const BoundType & A_Val1, const Interval & A_iVal2) {
     return (A_Val1 < A_iVal2.leftBound());
  }

  /// operator <=  (interval, scalar)
  friend inline bool operator <= (const Interval & A_iVal1, const BoundType & A_Val2) {
     return (A_iVal1.rightBound() <= A_Val2);
  }

  /// operator <=  (scalar, interval)
  friend inline bool operator <= (const BoundType & A_Val1, const Interval & A_iVal2) {
     return (A_Val1 <= A_iVal2.leftBound());
  }

//==== declarations in intervalFriend.h  definitions in intervalOp.hpp ==================

  friend Interval operator - (const Interval& A_iv){
    return Interval(-A_iv.m_interval);
  }
  friend Interval operator + (const Interval& A_iv1, const Interval& A_iv2) {
    return Interval(A_iv1.m_interval + A_iv2.m_interval);
  }
  friend Interval operator - (const Interval& A_iv1, const Interval& A_iv2) {
    return Interval(A_iv1.m_interval - A_iv2.m_interval);
  }
  friend Interval operator * (const Interval& A_iv1, const Interval& A_iv2) {
    return Interval(A_iv1.m_interval * A_iv2.m_interval);
  }
  friend Interval operator / (const Interval& A_iv1, const Interval& A_iv2) {
   try{
      return Interval(A_iv1.m_interval / A_iv2.m_interval);
    } catch(::cxsc::DIV_BY_ZERO & e){
      throw typename Interval::IntervalError("Division by 0 in operator/(Interval, Interval)", A_iv2.inf(),  A_iv2.sup());
    }

  }
  friend Interval operator ^ (const Interval& A_iv1, int i) {
    return Interval(power(A_iv1.m_interval, i));
  }


  //////////////////////////////////////////////////////////////////////////
  //  ARITHMETIC OPERATORS between interval and BoundType
  //////////////////////////////////////////////////////////////////////////


  /// operator +  (interval, scalar)
  friend inline Interval operator+ (const Interval & A_iVal, const BoundType &A_x) {
    return Interval(A_iVal.m_interval + A_x);
  }

  /// operator + (scalar, interval)
  friend inline Interval operator+ (const BoundType & A_x, const Interval & A_iVal) {
    return Interval(A_x + A_iVal.m_interval);
  }

  /// operator -  (interval, scalar)
  friend inline Interval operator- (const Interval & A_iVal, const BoundType &A_x) {
    return Interval(A_iVal.m_interval - A_x);
  }

  /// operator - (scalar, interval)
  friend inline Interval operator- (const BoundType & A_x, const Interval & A_iVal) {
    return Interval(A_x - A_iVal.m_interval);
  }

  /// operator * (interval, scalar)
  friend inline Interval operator* (const Interval & A_iVal, const BoundType &A_x) {
    return Interval(A_iVal.m_interval * A_x);
  }

  /// operator * (scalar, interval)
  friend inline Interval operator* (const BoundType & A_x, const Interval & A_iVal) {
    return Interval(A_x * A_iVal.m_interval);
  }

  /// operator / (scalar, interval)
  friend inline Interval operator/ (const Interval & A_iVal, const BoundType &A_x) {
    if(A_x != capd::TypeTraits<BoundType>::zero()){
      return Interval(A_iVal.m_interval / A_x);
    } else {
      throw typename Interval::IntervalError("Division by 0 in operator/(Interval, BoundType)", A_x, A_x);
    }
  }

  /// operator / (interval, scalar)
  friend inline Interval operator/ (const BoundType & A_x, const Interval & A_iVal) {
    if(! A_iVal.contains(capd::TypeTraits<BoundType>::zero())){
      return Interval(A_x / A_iVal.m_interval);
    } else {
      throw typename Interval::IntervalError("Division by 0 in operator/(BoundType, Interval)", A_iVal.inf(),  A_iVal.sup());
    }
  }

  friend std::ostream & operator << (std::ostream& s, const Interval & A_iv){
	//Interval::precision(s.precision());
    s << "[" << A_iv.inf() << ", " << A_iv.sup() <<"]";
    return s;
  }

  friend std::istream & operator >> (std::istream& s, Interval & A_iv){
    s >> A_iv.m_interval;
    return s;
  }
/* TODO: 
  friend std::ostream & bitWrite(std::ostream & out, const Interval & iv){
	  return iv.m_interval.bitImage(out);
  }
  friend std::istream & bitRead(std::istream & in, Interval & iv){
	  iv.m_interval = iv.m_interval.readBitImage(in);
    return in;
  }
  friend std::ostream & hexWrite(std::ostream & out, const Interval & iv){
	  return iv.m_interval.hexImage(out);
  }
  friend std::istream & hexRead(std::istream & in, Interval & iv){
	  iv.m_interval = iv.m_interval.readHexImage(in);
    return in;
  }
  friend std::ostream & binWrite(std::ostream & out, const Interval & iv){
	  return out.write((const char *)&iv.m_interval, sizeof(Interval::BaseInterval));
  }
  friend std::istream & binRead(std::istream & in, Interval & iv){
	  return in.read((char *)&iv.m_interval, sizeof(Interval::BaseInterval));
  }
*/

  
  

friend Interval diam(const Interval & ix);
friend BoundType width(const Interval & ix);
  
// {
//    return ix.m_interval.diam();
//  }

friend inline Interval mid(const Interval& A_iv){
  return  A_iv.mid();
}

  ///  Intersection of two intervals
  friend bool  intersection( const Interval & A_iv1,
                             const Interval & A_iv2,
                             Interval & A_rIntersection){
     BoundType left = (A_iv1.leftBound() > A_iv2.leftBound())? A_iv1.leftBound() : A_iv2.leftBound();
     BoundType right = (A_iv1.rightBound() < A_iv2.rightBound())? A_iv1.rightBound(): A_iv2.rightBound();

     if(left <= right) // is intersection nonempty ?
     {
     A_rIntersection = Interval(left, right);
     return true;
     }
     else
       return false;
   } // intersection

  /// returns an interval containing ix and iy
  friend Interval intervalHull(const Interval & ix,
                                       const Interval & iy){
    BoundType left = (ix.leftBound() < iy.leftBound())? ix.leftBound() : iy.leftBound(),
             right = (ix.rightBound() >iy.rightBound())? ix.rightBound() : iy.rightBound();
     return Interval(left, right);
  }

//  On output:  iv \subset Mid + [-diam , diam]

  friend inline void split( Interval & A_iv,
            Interval & A_rMid,
            BoundType & A_diam) {
    /// TODO : correct to not switch rounding
//     T_Rnd::roundUp();
//     T_Bound m = (A_iv.rightBound() + A_iv.leftBound()) /2,
//             l = m - A_iv.leftBound(),
//             r = A_iv.rightBound() - m;
//     A_rMid = Interval<T_Bound, T_Rnd>(m);
//     A_diam = ( l > r) ? l : r;

     BoundType m = A_iv.midPoint();
     Interval t = A_iv - m;
     A_rMid = Interval(m);
     A_diam = ( t.leftBound() > t.rightBound()) ? t.leftBound() : t.rightBound();
  }


  friend inline void split(Interval& A_rIv, BoundType & A_diam) {
    //split(A_rIv, A_rIv, A_diam);
    BoundType m = A_rIv.midPoint();
        Interval t = A_rIv - m;
        A_rIv = Interval(m);
        A_diam = ( t.leftBound() > t.rightBound()) ? t.leftBound() : t.rightBound();
  }


  friend inline bool isSingular(const Interval& A_x) {
    return ((A_x.leftBound()<=0) && (A_x.rightBound()>=0));
  }


// returns x^n

  friend inline Interval power (const Interval & x, int n){
     return Interval(power(x.m_interval, n));

  }

// returns a^b

  friend inline Interval power (const Interval & a,
                         const Interval & b){
    if( a.leftBound() >= 0 ){
       return Interval(pow(a.m_interval, b.m_interval));
    } else
          throw typename Interval::IntervalError(" power(A, B): Interval A must be nonnegative. \n", a.leftBound(), a.rightBound());

  }


// square root of x

  friend inline Interval sqrt (const Interval &x){

    if( x.leftBound() < 0 )
      throw typename Interval::IntervalError(" sqrt(x): Interval x must be nonnegative. \n", x.leftBound(), x.leftBound());

    return Interval(sqrt(x.m_interval));
  }


// sin x

  friend inline Interval  sin (const Interval& x){
    return Interval(sin(x.m_interval));
  }

// cos x

  friend inline Interval cos (const Interval& x){
    return Interval(cos(x.m_interval));
  }

// tan x

friend inline Interval tan (const Interval& x){
  return Interval(tan(x.m_interval));
}

// cot x

friend inline Interval cot (const Interval& x){
  return Interval(cot(x.m_interval));
}

// atan x

friend inline Interval atan (const Interval& x){
  return Interval(atan(x.m_interval));
}

// asin x

friend inline Interval asin (const Interval& x){
  // x must satisfy -1 <= x <= 1
  if ( ( x.leftBound() < -1 ) || ( x.rightBound() > 1) )
    throw typename Interval::IntervalError(" asin(A): Interval A must satisfy -1 <= A <= 1. \n", x.leftBound(), x.rightBound());

  return Interval(asin(x.m_interval));
}

// acos x

friend inline Interval acos (const Interval& x){
  // x must satisfy -1 <= x <= 1
  if ( ( x.leftBound() < -1 ) || ( x.rightBound() > 1) )
    throw typename Interval::IntervalError(" acos(A): Interval A must satisfy -1 <= A <= 1. \n", x.leftBound(), x.rightBound());

  return Interval(acos(x.m_interval));
}

// sinh x

friend inline Interval sinh (const Interval& x){
  return Interval(sinh(x.m_interval));
}

// cosh x

friend inline Interval cosh (const Interval& x){
  return Interval(cosh(x.m_interval));
}

// tanh x

friend inline Interval tanh (const Interval& x){
  return Interval(tanh(x.m_interval));
}

// coth x

friend inline Interval coth (const Interval& x){
  return Interval(coth(x.m_interval));
}

// exp x

friend inline Interval exp (const Interval & x){
  return Interval(exp(x.m_interval));
}

// natural logarithm of x

friend inline Interval log (const Interval& x){
  if( x.leftBound()<=0.0 )
    throw typename Interval::IntervalError(" log(A): Interval A must satisfy 0 < A. \n", x.leftBound(), x.rightBound());

  return Interval(::cxsc::ln(x.m_interval));
}


// square of x

friend inline Interval sqr (const Interval &x){
    if( capd::abs(x.leftBound()) > 1.3407807929942596e154 || capd::abs(x.rightBound()) > 1.3407807929942596e154 )
    throw typename Interval::IntervalError(" sqr(A): sqr overflow. \n", x.leftBound(), x.rightBound());

  return Interval(sqr(x.m_interval));
}


/// returns nonnegative part of interval
/// @remark if nonnegative part is empty throws exception

friend inline Interval nonnegativePart(const Interval &iv){
  if(iv.rightBound() < 0.0)
    throw typename Interval::IntervalError(" Nonnegative part is empty! ", iv.leftBound(), iv.rightBound());
  return Interval(capd::max(iv.leftBound(), BoundType(0.0)), iv.rightBound());
}

/// Ball with center iv and radius r

friend inline Interval ball(const Interval &iv,
                            const Interval & r)
{
    return Interval(iv.leftBound() - r.rightBound(),
                    iv.rightBound() + r.rightBound());
}

/// Ball with center iv and radius r

friend inline Interval ball(const Interval &iv,
                            const BoundType &r)
{
    return Interval(iv.leftBound() - r, iv.rightBound() + r);
}


/// solves inclusion a+[0,t]*p\subset c for t
friend Interval solveAffineInclusion(const Interval & a,
                              const Interval & p,
                              const Interval & c){
  Interval t;
  if ( !(a.subset(c)) ) {

    throw typename Interval::IntervalError("Cannot solve affine inclusion", a.leftBound(), a.rightBound());

  } else if(p >= 0) {
      t = ((c.right() - a.right()) / p.right()).left();
  } else if(p<=0) {
      t = ((c.left() - a.left()) / p.left()).left();
  } else {
      typename Interval::BoundType t1 = ((c.right() - a.right()) / p.right()).leftBound();
      typename Interval::BoundType t2 = ((c.left() - a.left()) / p.left()).leftBound();
      t = Interval( ( t1<t2 ? t1 : t2 ) );
  }
  return t;
} // solveAffineInclusion(Interval&, Interval&, Interval&)

/// solves inclusion a+[0,t]*p\subset c for t
friend Interval solveAffineInclusion(const Interval & a,
                              const Interval & p,
                              const Interval & c,
                              int & dir) {
  Interval t;
  if ( !(a.subset(c)) ) {
    throw typename Interval::IntervalError("Cannot solve affine inclusion", a.leftBound(), a.rightBound());

  } else if(p >= 0) {

    t = ((c.right() - a.right()) / p.right()).left();
    dir = 1;

  } else if(p<=0) {

    t = ((c.left() - a.left()) / p.left()).left();
    dir = -1;

  } else {

    typename Interval::BoundType t1 = ((c.right() - a.right()) / p.right()).leftBound();
    typename Interval::BoundType t2 = ((c.left() - a.left()) / p.left()).leftBound();
    if( t1 < t2 ){
      dir=1;
      t=Interval(t1);
    } else {
      dir=-1;
      t=Interval(t2);
    }
  }
  return t;
}



};  // end of class Interval
 


inline std::ostream & bitWrite(std::ostream & out, const Interval &iv){
	 out << "[";
     capd::intervals::IntervalIOTraits<double>::bitWrite(out, iv.leftBound());
     out << ",";
     capd::intervals::IntervalIOTraits<double>::bitWrite(out, iv.rightBound());
	 out << "]";
	 return out;
}

inline std::istream & bitRead(std::istream & inp, Interval &iv){
	double ll, rr;
	inp >> std::ws;
	if('['==inp.get()){
		inp >> std::ws;
		capd::intervals::IntervalIOTraits<double>::bitRead(inp, ll);
		// read white spaces
		inp >> std::ws;
		if(inp.get()==','){
			inp >> std::ws;
			capd::intervals::IntervalIOTraits<double>::bitRead(inp, rr);

			// read white spaces
			inp >> std::ws;
			if(inp.get()==']'){
				iv = Interval(ll, rr);
        return inp;
			}
		}
	} 
	// TODO : set fail bit in inp
  // iv.m_left = iv.m_right = 0.0;
	return inp;
}

inline std::ostream & hexWrite(std::ostream & out, const Interval &iv){
	out << "[";
	capd::intervals::IntervalIOTraits<double>::hexWrite(out, iv.leftBound());
	out << ",";
	capd::intervals::IntervalIOTraits<double>::hexWrite(out, iv.rightBound());
	out << "]";
	return out;
}

inline std::istream & hexRead(std::istream & inp, Interval &iv){
	double ll, rr;
	inp >> std::ws;
	if('['==inp.get()){
		inp >> std::ws;
		capd::intervals::IntervalIOTraits<double>::hexRead(inp, ll);
		// read white spaces
		inp >> std::ws;
		if(inp.get()==',')
		{
			inp >> std::ws;
			capd::intervals::IntervalIOTraits<double>::hexRead(inp, rr);

			// read white spaces
			inp >> std::ws;
			if(inp.get()==']')
			{
				//checkInterval(" bitRead ", ll, rr);
        iv = Interval(ll, rr);
				return inp;
			}
		}
	}
	// TODO : set fail bit in inp
  // iv.m_left = iv.m_right = 0.0;
	return inp;
}

inline std::ostream & binWrite(std::ostream & out, const Interval &iv){
	return out.write((const char *)&iv, sizeof(Interval));
}

inline std::istream & binRead(std::istream & in, Interval &iv){
	return in.read((char *)&iv, sizeof(Interval));
}


inline Interval diam(const Interval & ix){
    return ::cxsc::diam(ix.m_interval);
}

inline Interval::BoundType width(const Interval & ix){
    return ::cxsc::_double(::cxsc::diam(ix.m_interval));
}

}} // namespace capd::intervals

namespace capd{

/**
 * Specialization for intervals
 */
template<>
class TypeTraits < capd::cxsc::Interval > {
public:
  typedef  double Real;
  typedef capd::cxsc::Interval IntervalType;
//  static inline const ::capd::cxsc::Interval &zero(){
//    return S_zero;
//  }
  static inline const ::capd::cxsc::Interval zero(){
    return ::capd::cxsc::Interval(0.0);
  }
//  static inline const ::capd::cxsc::Interval & one(){
//    return S_one;
  static inline const ::capd::cxsc::Interval one(){
    return ::capd::cxsc::Interval(1.0);
  }
  static inline int numberOfDigits(){
    return TypeTraits<Real>::numberOfDigits();
  }
  /// Machine epsilon (the difference between 1 and the least value greater than 1 that is representable).
  static inline Real epsilon() throw(){
    return TypeTraits<Real>::epsilon();
  }

  template <typename S>
  static inline IntervalType convert(const S & obj){
    return static_cast<IntervalType>(TypeTraits<T>::convert(obj));
  }
  static inline IntervalType convert(const T & obj){
    return static_cast<IntervalType>(obj);
  }
  
  static const bool isInterval = true;

private:
  static const  ::capd::cxsc::Interval S_zero ;// = ::capd::cxsc::Interval<T,R>(T(0.0));
  static const  ::capd::cxsc::Interval S_one ;
};

//const ::capd::cxsc::Interval TypeTraits< ::capd::cxsc::Interval >::S_zero = ::capd::cxsc::Interval(0.0);
//const ::capd::cxsc::Interval TypeTraits< ::capd::cxsc::Interval >::S_one = ::capd::cxsc::Interval(1.0);


/// an absolute value
template<>
inline cxsc::Interval abs (const cxsc::Interval & A_inter){
   return cxsc::Interval(A_inter.abs());
} // abs


///maximum
template<>
inline cxsc::Interval max (      
  const cxsc::Interval & A_iv1,
  const cxsc::Interval & A_iv2
){
   return cxsc::Interval(
                  (A_iv1.leftBound()>A_iv2.leftBound() ? A_iv1.leftBound() : A_iv2.leftBound()),
                  (A_iv1.rightBound()>A_iv2.rightBound() ? A_iv1.rightBound() : A_iv2.rightBound())
          );
}

///minimum
template<>
inline cxsc::Interval min (
  const cxsc::Interval& A_iv1,
  const cxsc::Interval& A_iv2)
{
   return cxsc::Interval(
                  min(A_iv1.inf(),A_iv2.inf()),
                  min(A_iv1.inf(),A_iv2.inf())
          );
}






//////////////////////////////////////////////////////////////////////////
//        CONSTANTS
//////////////////////////////////////////////////////////////////////////
namespace cxsc{


inline Interval Interval::pi() {
  return Interval(::cxsc::Pi_interval);
}

inline Interval Interval::euler() {
   return Interval(::cxsc::E_interval);
}

}} // end of namespace capd::cxsc

#endif // _CAPD_FILIB_INTERVAL_H_
