//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file Complex.h
///
/// @author Tomasz Kapela   @date 2022-01-14
//
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) Tomasz Kapela 2022
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_FIELDS_COMPLEX_H_
#define _CAPD_FIELDS_COMPLEX_H_
#include <stdexcept>
#include <complex>
#include "capd/basicalg/power.h"
#include "capd/basicalg/TypeTraits.h"
#include "capd/basicalg/doubleFun.h"
#include "capd/basicalg/minmax.h"
namespace capd{
  namespace fields{

/**
 * Class Complex represents complex number.
 *
 * This is replacement for std::complex for user defined types especially intervals.
 * For built-in floating point types it uses std::complex.
 * It also adds TypeTraits and functions needed by vectalg algorithms (like capd::abs, capd::min).
 *
 * @tparam T type of real and imaginary part
 */

template <typename T>
class Complex{
    T re, im;
public:
    typedef T value_type;
    Complex( const T& re = T(), const T& im = T() ): re(re), im(im) {}
    Complex( const Complex & other ) = default;
    template <typename R>
    explicit Complex(const Complex<R> & other) : re(other.real()), im(other.imag()) {}
	  T real() const       { return re; }
    void real( T value ) { re = value; }
    T imag() const       { return im; }
    void imag( T value ) { im = value; }
    Complex& operator+=(const T& other){
        re += other;
        return *this;
    }
     Complex& operator+=(const Complex & other){
        re += other.re;
        im += other.im;
        return *this;
    }
     Complex& operator-=(const T& other){
        re -= other;
        return *this;
    }
     Complex& operator-=(const Complex & other){
        re -= other.re;
        im -= other.im;
        return *this;
    }
     Complex& operator*=(const T& other){
        re *= other;
        im *= other;
        return *this;
    }
     Complex& operator*=(const Complex & other){
        *this = *this * other;
        return *this;
    }
     Complex& operator/=(const T& other){
        re /= other;
        im /= other;
        return *this;
    }
     Complex& operator/=(const Complex & other){
        *this = *this / other;
        return *this;
    }


    friend  Complex operator+( const Complex & val) {
        return val;
    }
    friend  Complex operator-( const Complex & val){
        return Complex(-val.re, -val.im);
    }
    /// +
    friend   Complex operator+( const Complex & a, const Complex & b){
        return Complex( a.re+b.re, a.im+b.im );
    }
    friend   Complex operator+( const Complex & a, const T & b){
        return Complex( a.re+b, a.im);
    }
    friend   Complex operator+( const T & a, const Complex & b){
        return Complex( a+b.re, b.im );
    }
    /// -
    friend   Complex operator-( const Complex & a, const Complex & b){
        return Complex( a.re-b.re, a.im-b.im );
    }
    friend   Complex operator-( const Complex & a, const T & b){
        return Complex( a.re-b, a.im);
    }
    friend   Complex operator-( const T & a, const Complex & b){
        return Complex( a-b.re, -b.im );
    }
    /// *
    friend   Complex operator*( const Complex & a, const Complex & b){
        return Complex( a.re * b.re - a.im * b.im, a.im * b.re + b.im * a.re);
    }
    friend   Complex operator*( const Complex & a, const T & b){
        return Complex( a.re*b, a.im*b);
    }
    friend   Complex operator*( const T & a, const Complex & b){
        return Complex( a*b.re, a*b.im );
    }
    /// /
    friend   Complex operator/( const Complex & x, const Complex & y){
        Complex z = x * conj(y);// /
        Complex::value_type r = (y.re*y.re + y.im*y.im);
        z.re /= r;
        z.im /= r;
        return z;
    }
    friend   Complex operator/( const Complex & a, const T & b){
        return Complex( a.re/b, a.im/b);
    }
    friend   Complex operator/( const T & a, const Complex & b){
        return Complex(a)/b;
    }
    friend  bool operator==(const Complex & a, const Complex & b){
        return (a.re == b.re) and (a.im == b.im);
    }
    friend  bool operator!=(const Complex & a, const Complex & b){
        return !(a==b);
    }

    friend std::ostream & operator<< (std::ostream & out, const Complex & z){
//        if(z.re != 0.0 or z.im == 0.0){
//            out << z.re;
//        }
//        if(z.im != 0.0){
//            if(z.im > 0 and z.re != 0.0){
//                out << "+";
//            }
//            out << z.im << "i";
//        }
       out <<"(" << z.re <<", " << z.im << ")";
        return out;
    }

};

template<>
class Complex<double> : public std::complex<double> {
  typedef std::complex<double> Base;
 public:
  using Base::Base;
  Complex(const Base & x) : Base(x){}
};

template<>
class Complex<float> : public std::complex<float> {
  typedef std::complex<float> Base;
 public:
  using Base::Base;
  Complex(const Base & x) : Base(x){}
};

template<>
class Complex<long double> : public std::complex<long double> {
  typedef std::complex<long double> Base;
 public:
  using Base::Base;
  Complex(const Base & x) : Base(x){}
};

    /// returns the real component
    template< class T >
    T real( const Complex<T>& z ){
        return z.real();
    }
    /// returns the imaginary component
    template< class T >
    T imag( const Complex<T>& z ){
        return z.imag();
    }
    ///  returns the magnitude of a complex number
    template< class T >
    T abs( const Complex<T>& z ){
        return sqrt(z.real()*z.real() + z.imag()*z.imag());
    }
    /// returns the phase angle
    template< class T >
    T arg( const Complex<T>& z ){
        //TODO
        throw std::runtime_error("arg(Complex) not implemented!");
    }
    /// returns the squared magnitude
    template< class T >
    T norm( const Complex<T>& z ){
        return (z.real()*z.real() + z.imag()*z.imag());
    }
    /// returns the complex conjugate
    template< class T >
    Complex<T> conj( const Complex<T>& z ){
        return Complex<T>(z.real(), -z.imag());
    }
    /// constructs a complex number from magnitude r and phase angle theta
    template< class T >
    Complex<T> polar( const T& r, const T& theta = T()){
        return Complex<T>(r*cos(theta), r*sin(theta));
    }

    template< typename T>
    bool intersection(const Complex<T> & a, const Complex<T> & b, Complex<T> & result ){
        T re, im;
        bool exist = intersection(real(a),real(b),re) &&  intersection(imag(a),imag(b), im);
        if(exist)
            result = Complex<T>(re, im);
        return exist;
    }

    template< typename T>
    const Complex<T> intervalHull(const Complex<T> & a, const Complex<T> & b ){
      return Complex<T>(intervalHull(a.real(),b.real()),intervalHull(a.imag(),b.imag()));
    }

/// TODO: this is a naive implementation. Should be improved.
    template <typename T>
    inline Complex<T> power(const Complex<T> & x, int n){
        if(n<0)
            throw std::runtime_error("power implemented only for n>=0");
        Complex<T> res = capd::TypeTraits<Complex<T> >::one();
        for(int i=0; i<n; i++){
            res *= x;
        }
        return res;
    }

    template <typename T>
    Complex<T> mid(const Complex<T> & x){
        return  Complex<T>(mid(x.real()),mid(x.imag()));
    }
}  // end of fields namespace

template <typename T>
class TypeTraits< fields::Complex<T> > {

public:
  typedef typename fields::Complex<T> Type;
  typedef typename capd::TypeTraits<T>::Real Real;
  typedef typename capd::TypeTraits<T> RealTraits;

  static inline Type zero() {
    return Type(RealTraits::zero(), RealTraits::zero());
  }
  static inline Type one() {
    return Type(RealTraits::one(), RealTraits::zero());
  }
  static inline int numberOfDigits(){
    return RealTraits::numberOfDigits();
  }
  /// Machine epsilon (the difference between 1 and the least value greater than 1 that is representable).
  static inline T epsilon() throw(){
    return RealTraits::epsilon();
  }
  template <typename S>
  static inline Type convert(const S & obj){
    return static_cast<Type>(RealTraits::convert(obj));
  }
  static inline Type convert(const T & obj){
    return static_cast<Type>(obj);
  }
  /// this flag is true for all interval types
  static const bool isInterval = RealTraits::isInterval;


  static inline T abs(const fields::Complex<T> & x){
	return  sqrt(sqr(capd::abs(x.real()))+sqr(capd::abs(x.imag())));
  }

  static inline bool isSingular(const fields::Complex<T> & x) {
	return ( TypeTraits<T>::isSingular(x.real())) && (TypeTraits<T>::isSingular(x.imag()));
  }
};


namespace vectalg {
// Computes euclidean norm of a complex number
    template<typename T>
    T euclNorm(::capd::fields::Complex<T> &x) {
      return abs(x);
    }
  } // end of namespace vectalg
} // end of namespace capd



#endif // _CAPD_FIELDS_COMPLEX_H_
