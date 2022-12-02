//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file MpReal.h
///
///  Defines C++ wrapper for multiple precision floating point numbers
///  from mpfr library
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

#ifndef _CAPD_MULTIPREC_MPREAL_H_
#define _CAPD_MULTIPREC_MPREAL_H_


#include <cstdio>
#include <iostream>
#include <cmath>
#include <cstddef>
#include <stdexcept>

#include <gmp.h>
#include <mpfr.h>
#include "capd/basicalg/TypeTraits.h"
#include "capd/basicalg/minmax.h"
#include "capd/auxil/ignoreUnused.h"

namespace capd{
namespace multiPrec{

// Rounding modes
typedef mp_rnd_t MpRoundingMode;

/// Wrapper for mpfr precision type
class MpPrecision{
public:
  typedef mp_prec_t MpfrPrecisionType;
  MpPrecision(const MpfrPrecisionType precision) : m_precision(precision){
  }
  operator MpfrPrecisionType() const{
    return m_precision;
  }
private:
  MpfrPrecisionType m_precision;
};


/**
 * MpReal represents multiple precision real number with controlled rounding
 *
 * It is a C++ wrapper for a C library mpfr
 */
class MpReal {

  protected:
    mpfr_t mpfr_rep;  // representation of the real in a mpfr format

  public:
    typedef MpPrecision PrecisionType;
    typedef MpRoundingMode RoundingMode;

    // Rounding modes constants
    static const RoundingMode RoundUp=mp_rnd_t(GMP_RNDU);
    static const RoundingMode RoundDown=mp_rnd_t(GMP_RNDD);
    static const RoundingMode RoundNearest=mp_rnd_t(GMP_RNDN);
    static const RoundingMode RoundToZero=mp_rnd_t(GMP_RNDZ);

    // ----- Base functions ------     definitions in MpReal_Base.hpp
    // constructors and destructors
    MpReal ();
    MpReal (double d,
        RoundingMode rnd =  getDefaultRndMode(),
        PrecisionType prec = getDefaultPrecision());
    MpReal (long double d,
        RoundingMode rnd = getDefaultRndMode(),
        PrecisionType prec = getDefaultPrecision());
    MpReal (int i,
        RoundingMode rnd = getDefaultRndMode(),
        PrecisionType prec = getDefaultPrecision());
    MpReal (unsigned int i,
        RoundingMode rnd = getDefaultRndMode(),
        PrecisionType prec = getDefaultPrecision());
    MpReal (long int i,
        RoundingMode rnd = getDefaultRndMode(),
        PrecisionType prec = getDefaultPrecision());
    MpReal (unsigned long int i,
        RoundingMode rnd = getDefaultRndMode(),
        PrecisionType prec = getDefaultPrecision());
    MpReal(std::string s,
        RoundingMode rnd = getDefaultRndMode(),
        PrecisionType prec = getDefaultPrecision());
    MpReal(const char * str,
        RoundingMode rnd = getDefaultRndMode(),
        PrecisionType prec = getDefaultPrecision());
    MpReal (mpz_srcptr z,
        RoundingMode rnd = getDefaultRndMode(),
        PrecisionType prec = getDefaultPrecision());
    MpReal (mpq_srcptr q,
        RoundingMode rnd = getDefaultRndMode(),
        PrecisionType prec = getDefaultPrecision());
    MpReal (mpfr_t r,
        RoundingMode rnd = getDefaultRndMode(),
        PrecisionType prec = getDefaultPrecision());
    MpReal (const MpReal& r);

    ~MpReal ();

    /// Assignment operator
    MpReal& operator = (const MpReal& r);

	friend inline bool isinf(const MpReal& r)      { return (mpfr_inf_p(r.mpfr_rep) != 0);    }
    friend inline bool isnan(const MpReal& r)      { return (mpfr_nan_p(r.mpfr_rep) != 0);   }
	/// Produce a value with the magnitude of x and sign bit of y: (-1)^signbit(y) * abs(x)
    friend inline MpReal copysign(const MpReal & x, const MpReal &y) {
	  MpReal res;
	  mpfr_copysign(res.mpfr_rep, x.mpfr_rep, y.mpfr_rep, getDefaultRndMode() );
	  return res;
	};
	friend inline int signbit(const MpReal & x)    { return mpfr_signbit(x.mpfr_rep); }

    friend inline bool isZero(const MpReal& r)     { return (mpfr_zero_p(r.mpfr_rep) != 0);  }
    friend inline bool isInf(const MpReal& r)      { return isinf(r);   }
    friend inline bool isNaN(const MpReal& r)      { return isnan(r);   }
    friend inline bool isNumber(const MpReal& r)   { return (mpfr_number_p(r.mpfr_rep) != 0);}
    friend inline bool isSingular(const MpReal & r){ return (mpfr_zero_p(r.mpfr_rep) != 0);  }
    friend inline int  sign(const MpReal & r)      { return (mpfr_sgn(r.mpfr_rep));          }


    // Precision and rounding mode
    /// Sets precision for given MpReal object
    void setPrecision (PrecisionType newprec);
    /// returns precision of a given object
    PrecisionType getPrecision () const;

    /// sets default precision for all operations
    static void setDefaultPrecision (PrecisionType newprec);
    /// returns default precision of all operations
    static PrecisionType getDefaultPrecision ();
    /// sets default rounding mode for all operations
    static void setDefaultRndMode (RoundingMode newrndmode);
    /// return default rounding mode
    static RoundingMode getDefaultRndMode();
    /// sets rounding up mode
    static void roundUp();
    /// sets rounding down mode
    static void roundDown();
    /// sets rounding to nearest mode
    static void roundNearest();
    /// sets rounding to zero mode
    static void roundToZero();
    static void roundCut();


    // "Constants" (but they depend on the precision and the rounding mode)
    /// pi constant (3.1415...) with given precision and rounding direction
    static MpReal pi(RoundingMode rnd = getDefaultRndMode(),
        PrecisionType prec = getDefaultPrecision()) ;
    /// ln(2) - natural logarithm of 2 constant with given precision and rounding direction
    static MpReal log2(RoundingMode rnd = getDefaultRndMode(),
        PrecisionType prec = getDefaultPrecision()) ;
    /// euler constant (2.7...) with given precision and rounding direction
    static MpReal euler(RoundingMode rnd = getDefaultRndMode(),
        PrecisionType prec = getDefaultPrecision()) ;
    /// Positive infinity constant
    static MpReal positiveInfinity(PrecisionType prec = getDefaultPrecision());
    /// Negative infinity constant
    static MpReal negativeInfinity(PrecisionType prec = getDefaultPrecision());

    // ------------------- operators ---------------   definitions in MpReal_Op.hpp

    // Comparison operators
    /// compare two MpReals
    /// returns 0 if r1==r2, negative value if r1<r2 and positive value if r1>r2
    friend int compare (const MpReal& r1, const MpReal& r2);
    friend int compare (const MpReal& r1, const double r2);
    friend int compare (const MpReal& r1, const int r2);
    friend int compare (const MpReal& r1, const unsigned int r2);
    friend int compare (const MpReal& r1, const long int r2);
    friend int compare (const MpReal& r1, const unsigned long int r2);


    // Arithmetic operators
    /// adds r1 and r2 using default precision and rounding mode
    friend MpReal operator+ (const MpReal& r1, const MpReal & r2);
    /// adds r using default rounding mode
    MpReal& operator+= (const MpReal& r) ;
    /// adds r1 and r2 using precision of the res and given rounding mode
    static void add (MpReal& res,
        const MpReal& r1, const MpReal& r2,
        RoundingMode rnd = getDefaultRndMode());


    /// substract r1 and r2 using default precision and rounding mode
    friend MpReal operator- (const MpReal& r1, const MpReal& r2);
    /// substract r using default rounding mode
    MpReal& operator-= (const MpReal& r) ;
    /// substract r1 and r2 using precision of the res and given rounding mode
    static void sub (MpReal& res,
        const MpReal& r1, const MpReal& r2,
        RoundingMode rnd = getDefaultRndMode());
    MpReal operator- () const;
    static void neg (MpReal& res,
        const MpReal& r,
        RoundingMode rnd = getDefaultRndMode());

    /// multiply r1 and r2 using default precision and rounding mode
    friend MpReal operator* (const MpReal& r1, const MpReal& r2);
    /// multiply r using default rounding mode
    MpReal& operator*= (const MpReal& r) ;
    /// multiply r1 and r2 using precision of the res and given rounding mode
    static void mul (MpReal& res,
        const MpReal& r1, const MpReal& r2,
        RoundingMode rnd = getDefaultRndMode());

    /// divides r1 by r2 using default precision and rounding mode
    friend MpReal operator/ (const MpReal& r1, const MpReal& r2) ;
    /// divides by r using default rounding mode
    MpReal& operator/= (const MpReal& r) ;
    /// divides r1 by r2 using precision of the res and given rounding mode
    static void div (MpReal& res,
        const MpReal& r1, const MpReal& r2,
        RoundingMode rnd = getDefaultRndMode());

    /// fast multiplication and addition res = r1*r2 + r3
    static void fma (MpReal& res,
        const MpReal& r1, const MpReal& r2,
        const MpReal& r3,
        RoundingMode rnd = getDefaultRndMode());


    // --------------- Input/Output --------------------------------------   definitions in MpReal_Op.hpp

    /// writes to stream given number in the \b base with given number of digits and rounding mode
    std::ostream& put(std::ostream& o,
                 RoundingMode rnd = getDefaultRndMode(),
                 PrecisionType prec = getDefaultPrecision(),
                 int base = 10,
                 int nb_digits = 0) const;

    void get(const std::string & str,
    		RoundingMode rnd = getDefaultRndMode(),
    		PrecisionType prec = getDefaultPrecision(),
    		int base = 10);

    /// reads number from a stream
    friend std::istream& operator >> (std::istream &i, MpReal& r);
    /// writes number to a stream
    friend std::ostream& operator << (std::ostream &o, const MpReal& r);

    MpReal pow (const unsigned long int e, RoundingMode rnd = MpReal::getDefaultRndMode()) const;
    MpReal pow (const long int e, RoundingMode rnd = MpReal::getDefaultRndMode()) const;
    MpReal pow (const MpReal& e, RoundingMode rnd = MpReal::getDefaultRndMode()) const;

    // ----------- Mathematical functions: exp, sin... -------------------  definitions in MpReal_Fun.hpp

    void random (PrecisionType prec = getDefaultPrecision());

#include "capd/multiPrec/MpReal_FriendInternal.hpp"

    //  friend MpReal factorial (const unsigned long int e, RoundingMode rnd = MpReal::getDefaultRndMode());

    //---------------- end mathematical functions -----------------------

    // The four rounding modes as in IEEE-754 arithmetic and the fractional part
    friend MpReal round (const MpReal& r);
    friend MpReal floor (const MpReal& r);
    friend MpReal trunc (const MpReal& r);
    friend MpReal ceil (const MpReal& r);
    friend MpReal frac (const MpReal& r);


    // Conversion to built-in types (we do not want implicit conversions)


    // functions for compatibility when MpReal is treated as point interval
    // TODO:  library should be changed so that this functions are not needed (at least changed to friend function)
    MpReal leftBound()                             { return *this; }
    MpReal rightBound()                            { return *this; }
    bool contains(const MpReal & X)                { return (mpfr_cmp(mpfr_rep,X.mpfr_rep)==0); }
    bool containsInInterior(const MpReal & /*X*/)      { return false; }
};

inline const MpReal& leftBound(const MpReal& x){
  return x;
}

inline const MpReal& rightBound(const MpReal& x){
  return x;
}

inline MpReal abs(const MpReal& r, MpReal::RoundingMode rnd= MpReal::getDefaultRndMode()) {
  capd::auxil::ignoreUnused(rnd);
  if(compare(r,0)>=0)
    return r; // changes the precision of r
  // if getDefaultPrecision() != prec of r
    return -r;
}

}  // end of namespace multiPrec
}  // end of namespace capd


#include "capd/multiPrec/MpReal_Base.hpp"
#include "capd/multiPrec/MpReal_Op.hpp"
#include "capd/multiPrec/MpReal_Fun.hpp"

namespace capd{

template <>
class TypeTraits<capd::multiPrec::MpReal>{
public:
  typedef capd::multiPrec::MpReal Real;
  /// returns object set to zero
  static inline capd::multiPrec::MpReal zero(){
    return capd::multiPrec::MpReal(0.0);
  }
  /// returns object set to one
  static inline capd::multiPrec::MpReal one(){
    return capd::multiPrec::MpReal(1.0);
  }
  /// number of decimal digits
  static inline int numberOfDigits(){
    return capd::multiPrec::MpReal::getDefaultPrecision()* 3.010299956639812e-1; // constant = log(2)/log(10);
  }
  /// Machine epsilon (the difference between 1 and the least value greater than 1 that is representable).
  static inline capd::multiPrec::MpReal epsilon() {
      capd::multiPrec::MpReal r(1.0);
      return nextAbove(r) - r;
  }
  template <typename S>
  static inline Real convert(const S & obj){
    return static_cast<Real>(obj);
  }
  /// this flag is true for all interval types
  static const bool isInterval = false;

  static inline Real abs(const Real & r){
    return capd::multiPrec::abs(r);
  }

  static inline Real min(const Real & a, const Real& b){
    return std::min(a,b);
  }

  static inline Real max(const Real& a, const Real& b){
    return std::max(a,b);
  }
  static inline bool isSingular(const Real& x) {
    return ( x==zero());
  }
};

}

#endif // _CAPD_MULTIPREC_MPREAL_H_

#endif   // __HAVE_MPFR__


