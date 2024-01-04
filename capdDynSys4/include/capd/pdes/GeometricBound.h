/// @addtogroup pdes
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file GeometricBound.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2008-2020 by the CAPD Group.
//
// Distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_PDES_GeometricBound_H_
#define _CAPD_PDES_GeometricBound_H_

#include <stdexcept>
#include <iostream>

#include "capd/basicalg/minmax.h"
#include "capd/basicalg/power.h"
#include "capd/vectalg/lib.h"
#include "capd/vectalg/Vector_Interval.hpp"

namespace capd {
namespace pdes {
/**
 * The class is class represents a subset of a countable infinite dimensional space.
 * The representation splits into two parts.
 * 1. A finite dimensional interval vector holds finite number of coefficients 1,1,2,...,M called the main variables.
 * 2. The remaining coefficients x_i, x>M are bounded uniformly by a geometric series
 *    	|x_i| <= C*q^{-i}
 *    where
 *    - C is a non-negative constant
 *    - q is a geometric decay satisfying q>1
 *
 * The condition q>1 guarantees that the function \f$ x(t) = \sum_{i=0}^\infty x_i*t^i \f$ is analytic near t=0.
 *
 * The implementation provides operators and methods for acting on such functional objects, such as
 * - algebraic operations (addition, subtraction, multiplication)
 * - automatic differentiation with respect to time variable
*/
template<class ScalarT = capd::interval> // ScalarT can be complex interval
class GeometricBound {
public:
  typedef ScalarT ScalarType;
  typedef capd::IVector::size_type size_type;
  typedef capd::vectalg::Matrix<ScalarT,0,0> MatrixType;
  typedef GeometricBound VectorType;
  typedef capd::vectalg::Vector<ScalarT,0> FiniteVectorType;
  typedef typename MatrixType::RefRowVectorType RefVectorType;

  GeometricBound(){}

  /// constructs a Geometric series with C=0, q=1 and s=0 and given number of explicitly stored coefficients.
  explicit GeometricBound(size_type dim);

  /// Constructs a Geometric series with given bound on tail C/(q^i).
  /// The coefficients 1,...,dim will be stored explicitly.
  GeometricBound(size_type dim, capd::interval C, capd::interval q);

  /// Constructs a Geometric series with given bound on tail C/(q^i).
  /// The coefficients 1,...,dim will be stored explicitly and initialized from an array coeff.
  GeometricBound(size_type dim, capd::interval C, capd::interval q, const ScalarType* coeff);

  /// Constructs a Geometric series with given bound on tail C/(q^i).
  /// The coefficients 1,...,dim will be stored explicitly and initialized from an array coeff.
  GeometricBound(capd::interval C, capd::interval q, const FiniteVectorType& x);

  GeometricBound& operator+=(const GeometricBound& x);
  GeometricBound& operator-=(const GeometricBound& x);
  GeometricBound& operator*=(const ScalarType& x);

  GeometricBound partialDerivative() const;

  /// Returns i-th coordinate of a GeometricBound. The argument is an arbitrary natural number
  ScalarType getCoefficient(size_type i) const;

  /// Assigns a value to i-th coordinate of a GeometricBound.
  /// The argument i must be less or equal to the number of explicitly stored coefficients.
  void setCoefficient(size_type i, const ScalarType& s);

  /// Returns constant used in the bound of the infinite dimensional tail.
  capd::interval getConstant() const { return m_C; }

  /// Sets new value of constant used in the bound of the infinite dimensional tail.
  /// It must be positive number. Otherwise an exception is thrown.

  void setConstant(capd::interval C){
    if( C>=0 )
      m_C = C;
    else
      throw std::runtime_error("SetConstant error - negative argument");
  }

  /// Returns the constant q used in the bound of the infinite dimensional tail: C/(decay^i).
  capd::interval getGeometricDecay() const { return m_decay; };

  /// Sets new value of q used in the bound of the infinite dimensional tail: C/(q^i).
  /// If q is out of (0,1) an exception is thrown.
  void setGeometricDecay(capd::interval decay){
    if( decay>1 )
      m_decay = decay.rightBound();
    else
      throw std::runtime_error("SetGeometricDecay error - q is not bigger than 1.");
  }

  size_type dimension() const { return m_x.dimension(); }
  void setExplicitCoefficients(const FiniteVectorType& x) {
    m_x = x;
  }

  FiniteVectorType& getExplicitCoefficients() { return m_x; }
  const FiniteVectorType& getExplicitCoefficients() const { return m_x; }

  operator FiniteVectorType() const {return m_x;}

  RefVectorType projection(size_type dim) const {
    return RefVectorType(m_x.begin(),dim);
  }

  ScalarType operator[] (size_type i) const { return m_x[i]; }
  ScalarType& operator[] (size_type i) { return m_x[i]; }

  void clear(){
    m_x.clear();
    m_C = 0.;
  }

  const static size_type csDim = capd::IVector::csDim;

private:

  void check();
  FiniteVectorType m_x;
  capd::interval m_C;
  capd::interval m_decay;
}; // end of class GeometricBound

//*****************************************************************************/
/** operator+
 *
 * this operator realizes addition of two GeometricBound
 * using the following formula
 *     result_i = x_i + y_i
 *  Since object is represented as a finite dimensional vector and a tail
 *  we do explicit summation on main variables only.
 *  Exponent of result is computed as
 *     t = min(x.exponent(),y.exponent())
 *
 *  The constant used in representation of tail in result is computed as
 *  M = M_result = min(M_x,M_y) - number of exactly represented coefficients
 *  C_result = C_x/(M+1)^{x.exponent()-t} + C_y/(M+1)^{y.exponent()-t}
 *
 * @param[in]  x      object of class GeometricBound
 * @param[in]  y      object of class GeometricBound
 * @returns    sum of x and y
*/
template<class ScalarT>
GeometricBound<ScalarT> operator+(const GeometricBound<ScalarT>& x, const GeometricBound<ScalarT>& y);


/*****************************************************************************/
/** operator-
 *
 * this operator realizes subtraction of two GeometricBound
 * using the following formula
 *     result_i = x_i - y_i
 *  Since object is represented as a finite dimensional vector and a tail
 *  we do explicit summation on main variables only.
 *  Exponent of result is computed as
 *     t = min(x.exponent(),y.exponent())
 *
 *  The constant used in representation of tail in result is computed as
 *  M = M_result = min(M_x,M_y) - number of exactly represented coefficients
 *  C_result = C_x/(M+1)^{x.exponent()-t} + C_y/(M+1)^{y.exponent()-t}
 *
 * @param[in]  x      object of class GeometricBound
 * @param[in]  y      object of class GeometricBound
 * @returns    subtraction of x and y
*/
template<class ScalarT>
GeometricBound<ScalarT> operator-(const GeometricBound<ScalarT>& x, const GeometricBound<ScalarT>& y);


//*****************************************************************************/
/** operator*
  *
  *  this operator realizes multiplication of any coefficient
  *  in a GeometricBound by some scalar
  *
  *  using the following formula
  *     result_i = s * x
  *
  * @param[in]  s      object of type Scalar
  * @param[in]  x      object of class GeometricBound
  *
  * @returns    GeometricBound p multiplied by s
*/
template<class ScalarT>
GeometricBound<ScalarT> operator*(const ScalarT& s, const GeometricBound<ScalarT>& x);



//*****************************************************************************/
/** operator*
  *
  *  this operator realizes multiplication of any coefficient
  *  in a GeometricBound by some scalar
  *
  *  using the following formula
  *     result_i = s * x_i
  *
  * @param[in]  s      object of type Scalar
  * @param[in]  x      object of class GeometricBound
  *
  * @returns    GeometricBound x multiplied by s
*/
template<class ScalarT>
GeometricBound<ScalarT> operator*(const GeometricBound<ScalarT>& x, const ScalarT& s);



//*****************************************************************************/
/** operator*
  *
  *  this operator realizes linear transformation of a GeometricBound
  *
  *  if A is square matrix m\times m then the first m coordinates of result
  *  are given by
  *  \f$ result_i =  \sum_{k=1}^m A_{i,k} x_k \f$
  * The remaining (infinite number) of coordinates remain unchanged.
  *
  * @param[in]  A      square matrix of dimension smaller than p.dimension()
  * @param[in]  x      object of class GeometricBound
  *
  * @returns    GeometricBound x in transformed by (A,Id)
  *
  * @note The dimension m of the matrix A must be less or equal than number of exactly represented coefficients in x.
  *       Exception is thrown if this requirement is violated.
*/

template<class ScalarT>
GeometricBound<ScalarT> operator*(const capd::vectalg::Matrix<ScalarT,0,0>& A, const GeometricBound<ScalarT>& x);


//*****************************************************************************/
/** operator<<
  *
  * This operator writes a GeometricBound object to a given stream
  * in the following form
  * {{x_1,x_2,...,x_M},C,exponent}
  * where M is the number of exactly represented coefficients of x.
  *
  * @param[in]  out - a stream to which the object x is to be written
  * @param[in]  x   - an instance of class GeometricBound
  *
  * @returns  a reference to stream out
*/
template<class ScalarT>
std::ostream& operator << (std::ostream& s, const GeometricBound<ScalarT>& x);

template<class ScalarT>
GeometricBound<ScalarT> intersection(const GeometricBound<ScalarT>& x, const GeometricBound<ScalarT>& y);

template<class ScalarT>
void split(const GeometricBound<ScalarT>& X, GeometricBound<ScalarT>& x, GeometricBound<ScalarT>& dx);

// --------------------- inline definitions -----------------------------------

template<class ScalarT>
inline ScalarT GeometricBound<ScalarT>::getCoefficient(size_type i) const {
  if(i <= this->m_x.dimension())
    return m_x[i-1];
  return m_C*ScalarType(-1,1) / power(ScalarType(m_decay),i);
}

template<class ScalarT>
inline void GeometricBound<ScalarT>::setCoefficient(size_type i, const ScalarType& s) {
  if(i <= this->m_x.dimension())
    m_x[i-1] = s;
  else
    throw std::runtime_error("GeometricBound::setCoefficient - index out of bound");
}

template<class ScalarT>
inline GeometricBound<ScalarT> midVector(const GeometricBound<ScalarT>& x){
  const auto& r = x.getExplicitCoefficients();
  return GeometricBound<ScalarT>(0., x.getGeometricDecay(), capd::vectalg::midVector(r));
}

template<class ScalarT>
inline void swap(GeometricBound<ScalarT>& a, GeometricBound<ScalarT>& b){
  swap(a.getExplicitCoefficients(),b.getExplicitCoefficients());
  auto q = a.getGeometricDecay();
  a.setGeometricDecay(b.getGeometricDecay());
  b.setGeometricDecay(q);

  auto C = a.getConstant();
  a.setConstant(b.getConstant());
  b.setConstant(C);
}
}} // namespace capd::pdes


#endif // _CAPD_PDES_GeometricBound_H_


/// @}
