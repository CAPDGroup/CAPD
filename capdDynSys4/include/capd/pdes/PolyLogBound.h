/// @addtogroup pdes
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file PolyLogBound.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2008-2016 by the CAPD Group.
//
// Distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_PDES_POLYLOGBOUND_H_
#define _CAPD_PDES_POLYLOGBOUND_H_

#include <stdexcept>
#include <iostream>

#include "capd/basicalg/minmax.h"
#include "capd/basicalg/power.h"
#include "capd/vectalg/lib.h"

namespace capd {
namespace pdes {
/**
 * The class is class represents a subset of a countable infinite dimensional space.
 * The representation splits into two parts.
 * 1. A finite dimensional interval vector stores finite number of coefficients 0,1,2,...,M called the main variables.
 * 2. The remaining coefficients x_i, x>M are bounded uniformly by a polylog series
 *    	|x_i| < = C*q^{-i}i^{-d}
 *    where
 *    - C is a positive constant
 *    - q is a constant satisfying q>1
 *    - d is a natural number
 *
 * The condition q>1 guarantees that the function \f$ x(t) = \sum_{i=0}^\infty x_i*t^i \f$ is analytic near t=0.
 *
 * The implementation provides operators and methods for acting on such functional objects, such as
 * - algebraic operations (addition, subtraction, multiplication)
 * - automatic differentiation with respect to time variable
*/
class PolyLogBound {
public:
  typedef capd::interval ScalarType;
  typedef capd::IVector::size_type size_type;
  typedef capd::IMatrix MatrixType;
  typedef PolyLogBound VectorType;
  typedef capd::IVector FiniteVectorType;
  typedef MatrixType::RefRowVectorType RefVectorType;

  PolyLogBound(){}

  /// constructs a Geometric series with C=0, q=2 and d=0 and given number of explicitly stored coefficients.
  PolyLogBound(size_type dim);

  /// Constructs a Geometric series with given bound on tail C/(q^i).
  /// The coefficients 1,...,dim will be stored explicitly.
  PolyLogBound(size_type dim, ScalarType C, ScalarType q, size_type d);

  /// Constructs a Geometric series with given bound on tail C/(q^i).
  /// The coefficients 1,...,dim will be stored explicitly and initialized from an array coeff.
  PolyLogBound(size_type dim, ScalarType C, ScalarType q, size_type d, const ScalarType* coeff);

  /// Constructs a Geometric series with given bound on tail C/(q^i).
  /// The coefficients 1,...,dim will be stored explicitly and initialized from an array coeff.
  PolyLogBound(ScalarType C, ScalarType q, size_type d, const FiniteVectorType& x);

  PolyLogBound& operator+=(const PolyLogBound& x);
  PolyLogBound& operator-=(const PolyLogBound& x);
  PolyLogBound& operator*=(const ScalarType& x);

  PolyLogBound partialDerivative() const;
  friend std::ostream& operator << (std::ostream& s, const PolyLogBound& x);

  /// Returns i-th coordinate of a PolyLogBound. The argument is an arbitrary natural number
  ScalarType getCoefficient(size_type i) const;

  /// Assigns a value to i-th coordinate of a PolyLogBound.
  /// The argument i must be less or equal to the number of explicitly stored coefficients.
  void setCoefficient(size_type i, const ScalarType& s);

  /// Returns constant used in the bound of the infinite dimensional tail.
  ScalarType getConstant() const { return m_C; }

  /// Sets new value of constant used in the bound of the infinite dimensional tail.
  /// It must be positive number. Otherwise an exception is thrown.
  void setConstant(ScalarType C){
    if( C>=0 )
      m_C = C;
    else{
      std::ostringstream out;
      out << "PolyLogBound: setConstant error - negative argument" << C << "\n";
      throw std::runtime_error(out.str());
    }
  }
  void setDegree(size_type d) { m_degree = d; }
  size_type getDegree() const { return m_degree; }

  /// Returns the constant q used in the bound of the infinite dimensional tail: C/(decay^i).
  ScalarType getGeometricDecay() const { return m_decay; };

  /// Sets new value of q used in the bound of the infinite dimensional tail: C/(q^i).
  /// If q is out of (0,1) an exception is thrown.
  void setGeometricDecay(ScalarType decay){
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
  ScalarType m_decay;
  size_type m_degree;
  ScalarType m_C;
}; // end of class PolyLogBound

//*****************************************************************************/
/** operator+
 *
 * this operator realizes addition of two PolyLogBound
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
 * @param[in]  x      object of class PolyLogBound
 * @param[in]  y      object of class PolyLogBound
 * @returns    sum of x and y
*/
PolyLogBound operator+(const PolyLogBound& x, const PolyLogBound& y);


/*****************************************************************************/
/** operator-
 *
 * this operator realizes subtraction of two PolyLogBound
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
 * @param[in]  x      object of class PolyLogBound
 * @param[in]  y      object of class PolyLogBound
 * @returns    subtraction of x and y
*/
PolyLogBound operator-(const PolyLogBound& x, const PolyLogBound& y);


//*****************************************************************************/
/** operator*
  *
  *  this operator realizes multiplication of any coefficient
  *  in a PolyLogBound by some scalar
  *
  *  using the following formula
  *     result_i = s * x
  *
  * @param[in]  s      object of type Scalar
  * @param[in]  x      object of class PolyLogBound
  *
  * @returns    PolyLogBound p multiplied by s
*/
PolyLogBound operator*(const interval& s, const PolyLogBound& x);



//*****************************************************************************/
/** operator*
  *
  *  this operator realizes multiplication of any coefficient
  *  in a PolyLogBound by some scalar
  *
  *  using the following formula
  *     result_i = s * x_i
  *
  * @param[in]  s      object of type Scalar
  * @param[in]  x      object of class PolyLogBound
  *
  * @returns    PolyLogBound x multiplied by s
*/
PolyLogBound operator*(const PolyLogBound& x, const interval& s);



//*****************************************************************************/
/** operator*
  *
  *  this operator realizes linear transformation of a PolyLogBound
  *
  *  if A is square matrix m\times m then the first m coordinates of result
  *  are given by
  *  \f$ result_i =  \sum_{k=1}^m A_{i,k} x_k \f$
  * The remaining (infinite number) of coordinates remain unchanged.
  *
  * @param[in]  A      square matrix of dimension smaller than p.dimension()
  * @param[in]  x      object of class PolyLogBound
  *
  * @returns    PolyLogBound x in transformed by (A,Id)
  *
  * @note The dimension m of the matrix A must be less or equal than number of exactly represented coefficients in x.
  *       Exception is thrown if this requirement is violated.
*/
PolyLogBound operator*(const IMatrix& A, const PolyLogBound& x);


//*****************************************************************************/
/** operator<<
  *
  * This operator writes a PolyLogBound object to a given stream
  * in the following form
  * {{x_1,x_2,...,x_M},C,exponent}
  * where M is the number of exactly represented coefficients of x.
  *
  * @param[in]  out - a stream to which the object x is to be written
  * @param[in]  x   - an instance of class PolyLogBound
  *
  * @returns  a reference to stream out
*/
std::ostream& operator << (std::ostream& s, const PolyLogBound& x);
PolyLogBound intersection(const PolyLogBound& x, const PolyLogBound& y);
void split(const PolyLogBound& X, PolyLogBound& x, PolyLogBound& dx);

// --------------------- inline definitions -----------------------------------

inline
PolyLogBound::ScalarType PolyLogBound::getCoefficient(size_type i) const {
  if(i <= this->m_x.dimension())
    return m_x[i-1];
  return m_C*ScalarType(-1,1) / (std::pow(i,m_degree)*power(ScalarType(m_decay),i));
}

inline void PolyLogBound::setCoefficient(size_type i, const ScalarType& s) {
  if(i <= this->m_x.dimension())
    m_x[i-1] = s;
  else
    throw std::runtime_error("PolyLogBound::setCoefficient - index out of bound");
}

inline
PolyLogBound midVector(const PolyLogBound& x){
  return PolyLogBound(0., x.getGeometricDecay(), 0, midVector(x.getExplicitCoefficients()));
}

}} // namespace capd::pdes


#endif // _CAPD_PDES_POLYLOGBOUND_H_


/// @}
